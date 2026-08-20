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

## Measured

All eight files, both widths, byte-identical against `hisat2-build` on every
fixture plus a 20 Mb reference:

| reference | path nodes | C++ wall / RSS | Rust wall / RSS | scratch |
|---|---|---|---|---|
| 8 fixtures, 200 bp - 900 kb | 241 - 1.0M | — | 32-82 MB | — |
| 20 Mb, 77,843 variants | 21,171,995 | 14.9 s / 3,277 MB | 44.9 s / **197 MB** | 2.82 GB |

3.0x the wall for 16.6x less memory, single-threaded against a build that had 96
vCPU available. Peak RSS is not a function of the path-node count: it is
`graph::parse`'s joined text plus one F-bit rank per side.

`verify_wg.sh <fixture_dir>` runs the sweep; `HT2_LARGE=1` does the 64-bit half.
