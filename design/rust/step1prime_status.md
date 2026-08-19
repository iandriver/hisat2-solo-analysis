# Step 1′ — graph `.1.ht2` emission: everything but the M bitvector

## What is exact

Emitting the graph-mode gbwt block and comparing byte for byte against the index
`hisat2-build` wrote:

| graph | BWT rows | F bits | `fchr` | `zOffs` | `.2.ht2` samples compared |
|---|---|---|---|---|---|
| 200 bp | 242/242 | 242/242 | exact | exact | 15/15 |
| 509,431 bp | 532,723/532,723 | all | exact | exact | 31,998/31,998 |
| 900,000 bp | 957,345/957,345 | all | exact | exact | 56,621/56,621 |

The side layout is confirmed: within each side's `sideGbwtSz` bytes, BWT at 4
rows/byte low-pair-first in `[0, sz/2)`, the F bitvector at 8 rows/byte in
`[sz/2, 3sz/4)` with `F_bpi = bpi + ((sideCur & 1) << 2)`, M the same way after
it, and six `index_t` of tallies — `F_locSave`, `M_occSave`, `occSave[0..3]` —
holding the values as of the side's *start*.

`fchr`, `zOffs` and the `.2.ht2` offsets all fall out of the same walk and all
match, which is a strong check on the row semantics: `zOffs` is where the `'Z'`
row lands, and the SA samples are taken on `M == 1` boundaries at
`2^offRate` intervals.

## What is not

The **M bitvector**. The gbwt block matches everywhere except the M region:

| graph | gbwt bytes matching | first differing byte |
|---|---|---|
| 200 bp | 240/256 | 78 — the M region starts at 78 |
| 509,431 bp | 284,824/327,936 | 89 |
| 900,000 bp | 501,222/589,184 | 80 |

The **count** of `M == 1` rows is right — one per path node, 240 of them on the
200 bp graph. The **run lengths** are not.

Two rules were recovered correctly along the way:

- A node with out-degree **0 still emits one row**. `nextRow` sets `M` before
  testing whether to advance (`gbwt_graph.h:1632`), so the run length is
  `max(1, outdegree)`. That is what makes 241 out-degrees over 240 nodes produce
  242 rows.
- `F_loc` is the running sum of in-degrees, sampled once per `M == 1` row
  (`nextFLocation`, `:1640`).

What is still wrong is **which path node each edge's out-degree is credited to**.
Extracting HISAT2's actual run lengths from the stored bitvector: **238 nodes
have a run of 1 and exactly two have a run of 2, with no zeros.** Grouping edges
by `from` — assigning every edge with a given `from` to the first path node
carrying it — instead produces zeros and twos in different places, because
several path nodes legitimately share a `from` after full sorting (two distinct
paths of the same length starting at the same reference node, one per branch of
a variant).

The C++ writes this as a sequential merge-join (`gbwt_graph.h:2561`), and that
loop cannot be read as a group-by: replaying it literally against edges in
label-bucket order **stalls after exactly one bucket, 81 of 242 edges**, because
it only advances the edge cursor while `edge->from == node->from`. So either the
edges are in a different order at that point than the surrounding code suggests,
or `edge->from` means something other than a reference-graph node id there. That
is the single open question, and it is worth answering by instrumenting a debug
build of `hisat2-build` on the 200 bp graph rather than by more reading — the
graph is small enough to dump both lists in full.

## Consequence

`.1.ht2` cannot be emitted byte-identically until M is right, and the graph
`ftab` — which needs `mapGLF` walking a *correct* index — sits behind it. So
step 4's whole-genome byte diff is blocked on this one bitvector, not on scale.

Everything else in the file is now reproduced.
