# Plan: a whole-genome graph index builder in Rust

## The target, measured rather than assumed

The C++ build of exactly this index (E3, 2026-08-18):

| | |
|---|---|
| generations | 23 |
| peak path nodes | 5,920,332,832 (generation 11) |
| final path nodes | 5,917,131,871 |
| **peak RSS** | **625.2 GiB** |
| wall | 4 h 25 m on 96 vCPU |
| index | 11 GB |
| cost | ~$28 on an r7i.24xlarge |

**Goal: produce that index, byte for byte, on a workstation.** Not faster — the
C++ is already fast. The win is that a 768 GB instance stops being a
precondition, which is what currently makes a real-phasing panel a thing only
someone with cloud budget can rebuild.

We have the C++ output in S3, so the acceptance test is a byte comparison
against a known-correct 11 GB artifact. Most builders never get that.

## Where the 625 GiB actually is

`lateGeneration` holds three live `EList<PathNode>`: `past_nodes`, `from_table`,
`nodes`. `PathNode` is `{from, to, key.first, key.second}` — four `index_t`, so
**32 bytes at 64-bit**. At the generation-11 peak:

```
3 x 5.92e9 x 32 B = 568 GB payload
```

which the measured 625 GiB brackets with allocator and thread overhead. The
suffix array is **not** on this path — with variants present, `hisat2-build`
takes `buildToDisk(PathGraph&)`, not the blockwise-SA branch. The naive SA
measured in the rung-2 benchmark (15.3 bytes/bp, 10.7x slower at 50 Mb) is a
problem for the *linear* builder and for variant-free local indexes; it is not
what makes the whole-genome graph build need 768 GB.

## Step 1 — finish byte-identical graph emission (small scale)

Everything semantic is already reproduced and verified against real indexes:
`RefGraph` + reverse-determinisation, the full generation curve, `generateEdges`,
and the GFM rows character for character with every F bit correct, on graphs of
200 bp / 509 kb / 900 kb.

What remains is layout, which rung 2 showed is mechanical once the ordering is
right: side packing (2 rows/byte with the F and M bitvectors interleaved, 6
`index_t` of tallies per side), `zOffs`, the graph-mode `ftab`/`eftab`, the
`.2.ht2` offsets emitted on `M == 1` boundaries, and then `.5`/`.6`.

Gate: `.1`–`.6` byte-identical on the chr22-scale graph.

## Step 2 — widen the ids, and shrink the node

`u32` cannot hold 5.92e9 path nodes; that is not negotiable. But going to a flat
`u64` PathNode is the wrong reflex, because the four fields have different ranges:

| field | meaning | max at whole-genome scale | bits |
|---|---|---|---|
| `from`, `to` | reference-graph node ids | ~3.05e9 | **32** |
| `key.first` | rank among path nodes | 5.92e9 | **40** |
| `key.second` | rank / CSR offset | 5.92e9 | **40** |

A packed 18-byte node instead of 32 cuts the peak from 568 GB to **~320 GB** on
its own. That alone does not fit a workstation, but it halves the I/O volume for
step 3, which is where it pays.

Risk to watch: `from`/`to` at 3.05e9 leave only 28% headroom under `2^32`. A
larger panel (more haplotypes, a pangenome) overflows it, so the width should be
a compile-time parameter, not a hardcoded `u32`.

## Step 3 — external-memory doubling, which is the actual project

`createNewNodesMaker` is an equi-join:

```
past_nodes  ⋈  from_table   ON  past_nodes.to = from_table.from
```

In C++ it is a random-access probe into `from_table`, which is why all three
arrays must be resident. As an **external sort-merge join** both sides become
sequential streams:

1. sort `past_nodes` by `to` (external merge sort, bounded RAM)
2. stream-merge against `from_table` sorted by `from`, emitting new nodes
3. sort the output by key
4. stream `mergeUpdateRank` over it — it is already a single ordered pass

Per generation that is ~2 external sorts and ~2 streaming merges over ~6e9
records. At 18 bytes/node the working set is ~108 GB per pass; 23 generations at
roughly 4 passes each is on the order of **10 TB of sequential I/O**, which is
~3 hours on a single NVMe at 1 GB/s — comparable to the 4 h 25 m the C++ takes
in RAM, on a machine with 32–64 GB instead of 768 GB.

`mergeUpdateRank` is the one part that must be re-derived rather than ported: its
block walk assumes random access within a block, and the "can this node fold into
the previous one" lookahead needs a bounded window in a streaming form. The rule
itself is already pinned down — *a run sharing one `from` is one path node seen
several ways and collapses; a run spanning several `from` values is several nodes
at the same rank and all survive* — and it is verified generation-for-generation
on three graphs, so the streaming rewrite has an exact oracle to test against.

Gate: generation curve identical to the C++ log on chr22, then on chr1, with
peak RSS held under a declared budget.

## Step 4 — run it, and diff against the artifact

Build the whole genome from the same inputs (`genome.snp`, `genome.haplotype` in
S3) and byte-compare `.1`–`.8` against `s3://rustar-bench/hisat2-wg64/index/`.

Anything short of byte-identity is a bug with a known address, because every
intermediate — reference graph size, per-generation node and rank counts, GFM row
characters, F bits — already has its own oracle.

## What is deliberately not in this plan

**Order bounding.** Stopping the doubling at generation 8–9 would cut nodes from
5.9e9 to ~3.6e9 and the index by ~47%, and E1 showed no query HISAT2 issues
exceeds 131 bp, so the precision above order 256 is never used. But E2 tested it:
a truncated loop produces an index the aligner cannot query at all — 100 reads
unfinished in 120 s — because the GFM and the search assume ranks correspond
one-to-one with nodes. Capturing it means implementing merged-node semantics in
both construction and search, which is GCSA2's design and a larger project than
this one. It stays a later, bigger prize.

**Beating the C++ on speed.** Not a goal. If the external version lands within
2x of 4 h 25 m while running on a workstation, that is the win.

## Order of work

1. graph `.1`–`.6` byte-identical at chr22 scale — finishes the rung already in flight
2. packed node + parameterised id width — small, and halves step 3's I/O
3. external sort-merge doubling, with `mergeUpdateRank` re-derived for streaming
4. whole-genome run, byte-diff against the S3 index

Steps 1 and 2 are bounded and well-understood. Step 3 is the research risk, and
it is the only one.
