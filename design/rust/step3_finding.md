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
