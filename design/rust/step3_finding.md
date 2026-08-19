# Step 3 result — the doubling is solved, and the bottleneck moved somewhere worse

## What worked

`ht2ext` reproduces HISAT2's generation curve exactly at three scales, with the
sort buffer set small enough to force real spill runs and a real k-way merge:

| graph | generations | sort buffer | result |
|---|---|---|---|
| 200 bp, 3 variants | 6 | 4 KB | exact |
| 509,431 bp, 1,881 variants | 11 | 64 KB | exact |
| 900,000 bp, 3,502 variants | 11 | 64 KB | exact |

**The doubling now costs nothing.** Peak RSS of a full external run is identical
to graph construction alone, to within 50 KB:

| case | graph construction alone | full external run |
|---|---|---|
| 509,431 bp | 148,897,792 | 148,946,944 |
| 900,000 bp | 208,633,856 | 208,617,472 |

That was the whole point of step 3, and it holds: `createNewNodesMaker`'s
random-access probe is gone, replaced by a sort-merge join over sequential
streams, and the three resident `PathNode` arrays with it.

Two rules had to be recovered that a careful port still misses, both caught by
the curve rather than by reading:

1. **Generation 4 does not use the block walk.** `mergeUpdateRank` has a
   separate body for that generation built on `nextMaximalSet`. Using the block
   walk keeps 248 nodes where HISAT2 keeps 241.
2. **The block walk's lookahead discards a record.** After a multi-node block the
   next record is dropped outright when it is a block of one, the previous
   written node is sorted, and they share a `from`.

## What that exposed

With the doubling free, **`RefGraph` construction and `reverse_determinize` are
now the entire memory cost**, and they scale worse than what was removed:

| reference | peak RSS | bytes/bp | nodes removed by determinizing |
|---|---|---|---|
| 509,431 bp | 142 MB | 279 | 124 |
| 900,000 bp | 199 MB | 221 | 268 |
| 20,000,040 bp | **5.0 GB** | **262** | **0** |

Flat at roughly 260 bytes/bp. At 3,099,750,718 bp that projects to **~812 GB —
worse than the 625 GiB the C++ build needed.** Step 3 did not reduce the
whole-genome requirement; it relocated it.

The 20 Mb row is the diagnostic: **zero nodes were removed**, and it still cost
5 GB. The overhead is structural, not a consequence of merging. It is
`reverse_determinize`'s `cnode_map: HashMap<Vec<u32>, u32>` — one heap-allocated
`Vec` and one hash entry per composite node, ~100 bytes to carry what is usually
a single `u32`. Against ~15 MB of actual nodes, edges and text at 900 kb, better
than 90% of the peak is bookkeeping.

## The correction this forces

Earlier, on finding that the single-chain construction reproduces the node and
edge counts of HISAT2's fragmented-automaton branch exactly, this was recorded as
*"the determinisation is what makes the branch invisible."*

**That was right about the counts and wrong about the memory, and the memory is
why the branch exists.** `RefGraph` switches to per-fragment automata at
`jlen >= 1 << 16` (`gbwt_graph.h:383`) and determinises each fragment
separately. It is not a parallelisation detail. It is what bounds the
determinisation working set to one chunk instead of one genome — and reproducing
only its *output* while ignoring its *structure* is what put an 812 GB
projection into a plan whose entire purpose was to get under 625 GiB.

## Next step, which is now step 3.5 rather than step 1

Fragmented determinisation, matching `RefGraph`'s own chunking: split at ~1 Mb on
variant-free boundaries, determinise each fragment independently, and stream the
resulting nodes and edges to disk instead of accumulating them.

Two separate things have to be bounded, and only the first is the determiniser:

- **the determinisation working set** — one fragment at ~260 bytes/bp is ~260 MB
  per 1 Mb chunk, which is fine;
- **the graph itself** — 3.05e9 nodes and edges at 8 bytes each is ~49 GB even
  packed, so it cannot stay in a `Vec`. It has to be a record file, the same
  treatment the path nodes just received.

Only once both are external does a whole-genome run become a workstation job,
and only then is the step-4 byte-diff against `s3://rustar-bench/hisat2-wg64/`
reachable.

## Honest status against the four steps

| step | state |
|---|---|
| 1 — byte-identical graph emission | blocked: the graph `ftab` needs GFM search primitives (`mapGLF`), not layout |
| 2 — packed 18-byte node | done |
| 3 — external doubling | **done and verified at three scales; doubling memory is now ~0** |
| 3.5 — fragmented determinisation + external graph | **new, and now the blocker** |
| 4 — whole-genome run and byte diff | needs 1, 3.5 |

The plan said step 3 was the only research risk. It was the only risk *in the
doubling*; measuring it honestly surfaced a second one that the earlier
count-level agreement had hidden.


## Item 3 done: the determiniser's bookkeeping, and where the peak really sits

`reverse_determinize` now keeps member lists in one arena instead of a `Vec` per
composite node, and resolves the 99.97% singleton case through a flat table
indexed by node id instead of hashing a `Vec`. All correctness held — generation
curves, GFM rows and F bits still exact on all three graphs, external curve still
exact.

| reference | before | after |
|---|---|---|
| 509,431 bp | 279 B/bp | **106 B/bp** |
| 900,000 bp | 221 B/bp | **92 B/bp** |
| 20,000,040 bp | 262 B/bp | **88 B/bp** |

Whole-genome projection: **812 GB -> 272 GB**. Real, and not enough.

### The attribution that reorders the remaining work

Splitting the 20 Mb peak with `HT2_NO_DET=1`:

| | peak | bytes/bp |
|---|---|---|
| graph construction only | 827 MB | **43.4** |
| plus determinisation | 1,671 MB | **87.6** |

Almost exactly half each: ~135 GB for the graph arrays at whole-genome scale and
~137 GB for the determiniser. **Neither item alone reaches a workstation, and
they are not independent.**

The plan had item 1 (stream the graph to disk) before item 3.5 (fragment the
determiniser). That ordering does not work: `reverse_determinize` needs random
access to the node labels and to the predecessor index, so the graph cannot be
streamed out from under it. Fragmentation is the precondition for streaming, not
a later optimisation — which is precisely why `RefGraph` fragments in the first
place.

**So items 1 and 3.5 are one job**: split at ~1 Mb on variant-free boundaries,
determinise each fragment with everything resident (~90 MB at 88 B/bp), append
its nodes and edges to record files with a global id offset, and drop it. Peak
becomes one fragment plus the doubling's sort buffer.

## Item 1: the graph streams to disk, and the memory curve finally bends

`build_fragmented_to_disk` appends each fragment's nodes and edges to record
files and drops the fragment, so nothing about the graph accumulates.
`ht2ext`'s `HT2_DISK=1` builds generation 0 by streaming those files — edges are
in `sortEdgesFrom` order and nodes in id order, so the label each edge needs is a
single forward pass, not a random probe.

### One trap on the way, worth recording

With the graph on disk, 20 Mb still cost 400 MB and the per-bp figure barely
moved when the fragment chunk shrank — so it was not the fragment.

**The external sort's merge phase was using memory proportional to `n`.** Every
run was opened with a 1 MB `BufReader`, and the run count is `n / budget`: at
20 Mb that is ~305 runs, so **305 MB of merge buffers behind a 1.2 MB sort
budget.** The sort phase respected its budget; the merge phase silently did not.
Splitting one buffer budget across the runs fixes it.

| reference | before | after |
|---|---|---|
| 509,431 bp | 18.4 MB | 19.8 MB |
| 900,000 bp | 25.5 MB | 20.7 MB |
| 20,000,040 bp | **399.6 MB** | **84.6 MB** |

### Where that leaves the whole-genome projection

The curve is now a fixed ~18 MB plus roughly **3.3 bytes/bp**, fitted between the
900 kb and 20 Mb points. At 3,099,750,718 bp:

| | |
|---|---|
| **projected peak RSS** | **~10 GB** |
| C++ measured peak | 625 GiB |
| disk, scaled from 2.0 GB at 20 Mb | ~300 GB |

**That is a workstation.** Which was the point.

The remaining ~3.3 bytes/bp is mostly the reference text at 1 byte/bp — 2-bit
packing would take it to 0.25 — plus the variant tables at ~0.2 bytes/bp at
whole-genome density. Neither is a wall.

### Everything still exact

Generation curves, GFM rows, F bits, fragmented-equals-global, and the external
curve, on all three graphs, through every change above.

## Status against the four steps

| step | state |
|---|---|
| 1 — graph streamed to disk | **done** |
| 2 — packed 18-byte node | **done** |
| 3 — external doubling | **done, verified at three scales** |
| 3.5 — fragmented determinisation | **done, verified identical to global** |
| 1' — byte-identical graph emission | still blocked on GFM search primitives for the graph `ftab` |
| 4 — whole-genome run and byte diff | needs 1' |

The memory work is finished and measured. What stands between here and a
whole-genome byte diff is `mapGLF` and the side packing, not scale.
