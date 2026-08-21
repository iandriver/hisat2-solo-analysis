# `ht2wg` — a whole graph index off the external path

`ht2path` and `ht2emit5` reproduce every section of a HISAT2 graph index byte
for byte, but both hold the whole path graph in RAM. That is the 625 GiB the
C++ build needs, reproduced rather than avoided. `ht2wg` runs the same
construction over the fragmented graph builder and the external doubling loop,
and writes the files instead of comparing against an existing index.

```
ht2wg <reference.fa> <snp> <haplotype> <workdir> <out_prefix> [budget_records] [verify_prefix]
  HT2_LARGE=1   build a 64-bit (.ht2l) index
  HT2_CHUNK=n   graph fragment size in bases (default 262144)
  HT2_KEEP=1    leave the workdir in place
```

## What had to stop being an array

The in-memory version is the specification; the differences are all about what
is allowed to be resident.

| in `ht2path` | at whole-genome scale | in `ht2wg` |
|---|---|---|
| `start[]`, a CSR over reference-graph node ids | 3.05e9 + 1 entries, 12 GB | a sort-merge join of the edges by their `to` against the path nodes by their `from` |
| `outdeg[]`, `indeg[]`, `floc[]` | 5.9e9 entries each | three streams, consumed in node order |
| `occ`/`rank_M`/`select_F` precomputed over every row | 5.9e9 × 24 B | side-local queries against the block just written |

Only one thing stays in RAM: an F-bit rank per side. `select_F` is the single
primitive with no per-side tally to stand on, and one `u64` per side is 228 MB
for a human graph.

## The ftab is where the transcription stops being literal

The graph `ftab` is built by *querying the index that was just written* — for
each of the `4^ftabChars` prefixes, walk the GFM backward and record the row
range (gfm.h:4993). Transcribed directly that is 10,485,760 LF steps for
`ftabChars` 10, and every one of them a random read into an 11 GB file.

But `mapGLF` consumes a prefix from its **last** character first (`c = q & 3;
q >>= 2`), so prefixes sharing low-order bits share a walk. Done as a depth-first
traversal of that trie it costs at most `4 + 4^2 + ... + 4^10` ≈ 1.4M steps, and
prunes a whole subtree the moment a range goes empty:

| fixture | trie steps | flat walk | ratio |
|---|---|---|---|
| `tinyidx` (200 bp) | 4,916 | 10,485,760 | 2,133× |
| `multiidx` (15.7 kb) | 181,732 | 10,485,760 | 58× |
| `cleanidx` (509 kb) | 954,660 | 10,485,760 | 11× |
| `exidx` (900 kb) | 1,103,756 | 10,485,760 | 9.5× |

The ceiling is 1.4M however large the reference gets, because the trie is over
the prefix space and not over the index. A whole-genome ftab is therefore ~1.4M
random side reads, not 10.5M.

## Two bugs the in-memory oracle could not have caught

**`generateEdges` re-sorts its nodes by rank** between building the label buckets
and the merge-join that assigns out-degrees (`nodes.sort_by_key(|n| n.k0)`).
The buckets need them by `from`; the merge-join needs them by rank, because the
`ni` it assigns is a rank index. Walking them in `from` order makes the
positional scan stall on its first mismatch — 1 of 242 edges on the 200 bp
graph. In RAM the two orderings are the same `Vec` sorted twice and the second
sort is easy to read past; as separate files it has to be deliberate.

**A node with no genomic position carries `(index_t)INDEX_MAX`.** The graph is
built with 32-bit ids, so the sentinel arrives at the writer as `u32::MAX` and a
large index needs it **widened**, not zero-extended. At 32 bits the two are the
same number, so only a `.ht2l` fixture can show it: 14 of 33,180 SA samples on
`cleanidx`, and the `.1.ht2l` around them byte-identical.

That second bug is why the fixture set now has a 64-bit half. A whole-genome
index is necessarily large, and `index_t = uint64_t` changes every `writeIndex`
field, moves the default line rate from 7 to 8 (six 8-byte tallies per side
rather than six 4-byte ones), and renames the files — none of which the 32-bit
fixtures exercise. The local indexes are `local_index_t` = `u16` either way, so
`.5.ht2l` differs from `.5.ht2` by exactly four fields: the index count and each
window's `tidx`, `localOffset` and `joinedOffset`.

## Two things that only break at scale

**Merge fan-in.** 5.9e9 records at an 18 MB sort budget is ~5,900 runs, and a
single k-way merge holds one open file and one buffer per run. That hits the
file-descriptor limit long before it hits the memory budget, so the sort merges
in passes of at most 64 runs.

**Slurping the inputs.** `graph::parse` read its FASTA with `read_to_string`,
which for a 3.2 GB reference doubles peak memory for the length of the parse and
buys nothing. All three input files are streamed now.

## `.7` / `.8` are not part of the graph

They are a straight serialisation of the alt and haplotype lists (gfm.h:1912) —
the variant database `Zs:Z` gets its rsIDs from — and they need none of the
above. Two details decided them: the repeat block that follows the haplotypes is
empty unless `--repeat-ref` was given, which is why the fixtures' `.7` ends
exactly at the last haplotype; and `_alts` is **sorted** before anything is built
with it, with the names and every haplotype's alt indices permuted to match.
That sort is a no-op on every fixture here because variant files arrive
position-sorted, but a file that did not would otherwise build a different index
than `hisat2-build` does.

## Scratch: the files have to shrink while they are read

Almost every file in the pipeline is written once and read once. Held as whole
files, each is a full copy of the node list sitting on disk for the length of the
pass that consumes it. Tracing what was live at the peak on a 20 Mb reference
found five of them at once:

```
joined.bin  sorted.bin  by_to.bin  by_from.bin  cur.bin      2.82 GB
```

— 133 bytes per path node, or 7.4 copies of an 18-byte node. Three changes, in
increasing order of how much they had to change:

1. **Delete at the point a file stops being readable, not at the end of the
   phase.** `cur` is dead the moment both orderings of it exist; `by_to` and
   `by_from` are dead the moment the join returns. → 1.34 GB.

2. **Free each sort run as it empties**, and let a sort consume its source once
   the run phase has read it. → the peak moved into `generateEdges`, where
   deleting `nbf` after the join and building the rank-ordered copy *after* it
   rather than beside it took the same treatment.

3. **Segment every record file**, so a one-pass reader hands the space back as it
   goes (`SegReader` with `consume`). A source now shrinks at the rate its output
   grows instead of both being resident. → **1.03 GB**.

Two things went wrong on the way there and are worth keeping:

- **`seg_remove` probing 0, 1, 2, ... stops at the first gap.** A consuming
  reader deletes a *prefix* of the segments, so a base read part-way has a hole
  at the front — and the join does leave `by_from` part-way, because it stops
  when the `by_to` side runs out. Probing from zero found nothing and leaked
  everything behind the hole. It lists the directory now; `seg_truncate` is the
  cheap forward walk, for files that are either untouched or fully consumed.
- **A run file is only `budget` records long**, which was shorter than the
  default segment size — so segmentation was silently off for exactly the files
  a merge needs to shrink, and eleven half-drained runs still occupied their full
  size. Run files get their own segment cap.

What is left is the floor for this structure: the join needs the nodes ordered by
`to` and by `from` at the same time, and both are full copies. 36 bytes per path
node, two copies of an 18-byte node, plus the reference graph on disk.

## Threading, and what it was worth

Timing the stages first turned out to matter more than the threading did. The
profile said 50% of the build was external sorts inside the doubling and 4% was
`mergeUpdateRank` -- but the generation curve said something louder: unsorted
nodes collapse after generation 5, so five of eleven generations were sorting 21
million records to move 99.9% of them unchanged. Not sorting them is worth more
than sorting them on more cores.

| | 20 Mb build | what it was |
|---|---|---|
| start | 62.7s | |
| skip the pass-through nodes | ~51s | no threads involved |
| one local index per core | ~50s | 7.5s -> 1.2s on that stage |
| ftab: a side cache, one thread | ~47s | threading it made it *slower* |
| sort runs on every core | ~36s | keep the chunk size, not the total |
| merge split by key range | ~19s | the big one |
| join split by its join value | **18.7s** | modest; `join_late` is untouched |

**3.35x, and two of the six steps were not threading at all.** Peak RSS is now
177 MB on one thread and ~930 MB on eight, because `budget` is the chunk size
per worker; whoever runs it divides. Peak scratch went 1.03 -> 1.20 GB, which is
the partitioned merge's second copy of the runs.

Three findings worth more than the speedup:

- **The ftab got slower with threads** -- 2.9s on one against 4.3s on eighteen.
  Six `pread` calls per step on one file, 8.3M of them, and at this scale the
  block is page-cached so the syscall contention is all there is. Caching four
  sides per worker removed most of the calls and the stage went 4.9s -> 2.9s
  without any concurrency at all.
- **Dividing the memory budget among workers made the sort slower at every
  thread count.** Smaller chunks mean more runs, and the merge fan-in that buys
  costs more than the parallel sort saves.
- **`mergeUpdateRank` stays serial**, as it is in `gbwt_graph.h`, and is now 20%
  of the doubling. It is the one stage with a genuine prefix dependency.

## Resume

A build that runs for hours has to survive a failure at hour five. It resumes by
being re-run with the same arguments; `HT2_FRESH=1` starts over, and
`HT2_NO_RESUME=1` turns the whole thing off.

The doubling was already checkpointed by its own shape -- at the end of a
generation the node file is a complete, self-describing state -- so most of the
work was making that true at finer grain and then proving it.

**The node file alternates between two names.** A generation reads one and writes
the other, so the state a checkpoint names is never the state the next step is
overwriting. That cannot be had from a single `cur.bin` each generation consumes
in place: once generation g's sorts begin, generation g-1 is gone and there is
nothing to fall back to.

**Each phase's input outlives the step that replaces it.** A generation is
`sorted` -> `joined` -> `keyed` -> done, and every artifact stays until the mark
for the phase after it lands. That is the cost: peak scratch 1.20 -> 1.73 GB at
20 Mb, 82 bytes per path node against 57, which is why turning it off is an
option rather than an argument.

**The fingerprint includes the builder's own mtime and size.** Finishing half a
build with a different version of the construction is the stale-binary failure
again -- everything looks consistent and nothing is.

### Proving it, rather than hoping

Racing a `kill -9` finds bugs but cannot say which boundary it found them at, or
prove a boundary was ever reached. `HT2_CRASH_AT` names a checkpoint boundary and
`HT2_CRASH_GEN` narrows it to a generation, so `verify_resume.sh` visits every
place a restart can land, on purpose, and finishes each one. It found four bugs
that the random kills had only hinted at:

- **`join_par` deleted its own inputs.** Harmless in a straight run, where the
  caller deletes them anyway; fatal on a restart, which re-runs the join from
  files that are already gone. The join produces nothing, the generation reports
  zero nodes, and zero equals zero so the loop calls it converged.
- **Convergence was recorded by a second save.** In the window between them the
  state said "generation 10 done, keep going" -- and going on runs an eleventh
  generation over a converged file, which converges again and puts a generation
  in the curve that `hisat2-build` never had.
- **`cur` was deleted before the checkpoint saying it was no longer needed.**
  The other order says the node file is expendable before anything has recorded
  what replaced it.
- **Renaming `joined` onto the generation's output is not idempotent.** A restart
  landing after it found `joined` already moved, cleared the destination, and
  moved nothing onto it. Generations 1-3 have no second phase at all, so the join
  now writes straight into the output and there is nothing to hand over.

The guard that turned all four from silent into loud: **a generation that
produces no nodes panics.** Without it they all look like early convergence.

One thing is deliberately not resumable: the reference-graph build. It is a
single stage with no internal boundary, so a crash inside it restarts it -- 30-60
minutes at whole-genome scale against the hours the doubling would have cost.

## Measured

All eight files, both widths, byte-identical against `hisat2-build` on every
fixture plus a 20 Mb reference:

| reference | path nodes | C++ wall / RSS | Rust wall / RSS | scratch |
|---|---|---|---|---|
| 8 fixtures, 200 bp - 900 kb | 241 - 1.0M | — | 32-82 MB | — |
| 20 Mb, 77,843 variants | 21,171,995 | 14.9 s / 3,277 MB | 59 s / **182 MB** | **1.03 GB** |

16.6x less memory, single-threaded against a build that had 96 vCPU available,
at 2.7x the disk it used before.

**On the wall time.** It measured 45 s before the scratch work and 59 s
after, on an otherwise idle machine, but the attribution is unsettled and the honest reading is "somewhere
under 1.5x, direction not fully established". Segmentation itself is free --
forcing one segment per file (`HT2_SEG=1000000000`) gives 61.8 s against 61.2 s
with segments, and the difference is inside the noise. Making the record reader
and writer block-based instead of 18-byte `read_exact`/`write_all` calls through
`BufReader` moved user time by 0.03 s, so the per-record path was not it either.
What remains is the delete-as-you-go work itself and the machine: this is an
Apple M5 Pro with heterogeneous cores, the 45 s baseline was measured on an idle
machine, and every run since has shared it. Worth re-measuring cleanly before
anyone treats the slowdown as real. Peak RSS is not a function of the path-node count: it is
`graph::parse`'s joined text plus one F-bit rank per side.

`verify_wg.sh <fixture_dir>` runs the sweep; `HT2_LARGE=1` does the 64-bit half.
