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


---

# Resolved — the whole gbwt block, byte for byte

Instrumenting a throwaway `hisat2-build` to dump `nodes` and `edges` at the
merge-join answered it in one run, after source reading had gone in circles:

```
@@N 0  from=196     @@E 0  from=196  ranking=4  A
@@N 1  from=48      @@E 1  from=48   ranking=5  A
@@N 2  from=17      @@E 2  from=17   ranking=7  A
@@N 3  from=128     @@E 3  from=128  ranking=9  A
```

**Edge `i`'s `from` equals node `i`'s `from`.** The two lists arrive positionally
aligned, which is why a sequential scan is the right algorithm and why reading it
as a group-by was wrong. The scan pairs them off and advances the node only on a
mismatch, so a node collects two rows exactly when two edges land on it.

That also exposed the actual defect. `PathEdge(edge->from, nodes[j].key.first,
label)` stores the **reference edge's** `from`, and we were storing the target
path node's own `from`. Only the merge-join reads that field, so the BWT rows and
F bits stayed perfect while every out-degree was attributed to the wrong node —
an error no aggregate check could see.

A second, independent bug sat behind it. The C++ writes the six per-side tallies
when a side *fills*, using `occSave`: the counts as they stood when that side
*opened*. Writing them at the side's start instead — which is where they belong —
the live counters already hold exactly those values, so carrying `*_save`
variables across was not merely redundant but **wrong by a full side**. Side 0
matched (all zero) and every later side did not.

| graph | gbwt block | `fchr` | `zOffs` | `.2.ht2` SA sample |
|---|---|---|---|---|
| 200 bp | **256/256** | exact | exact | 15/15 |
| 509,431 bp | **327,936/327,936** | exact | exact | 33,180/33,180 |
| 900,000 bp | **589,184/589,184** | exact | exact | 59,615/59,615 |

Every byte of the graph-mode gbwt block, plus `fchr`, `zOffs`, and the complete
`.2.ht2` offset array — the sample counts now match too, where before ours fell
short because the M runs were wrong.

## What is left in `.1.ht2`

Only the **graph `ftab`/`eftab`**, which is built by querying the finished index:
`mapGLF`/`mapGLF1` walking the GFM backward for all 1,048,576 prefixes
(`gfm.h:4993`). That needs `SideLocus` and occ navigation over the block — and
the block is now known to be correct, so those primitives can be developed
against a verified target rather than a hypothesis.

Then the header, `nPat`/`plen`/`nFrag`/`rstarts` and `refnames` — all of which
rung 2 already emits byte-identically for the linear case and which are
unchanged here.


---

# The graph `ftab`/`eftab` — exact

| graph | gbwt block | `ftab` | `eftab` |
|---|---|---|---|
| 200 bp | 256/256 | **1,048,577/1,048,577** | **40/40** |
| 509,431 bp | 327,936/327,936 | **1,048,577/1,048,577** | **1,110/1,110** |
| 900,000 bp | 589,184/589,184 | **1,048,577/1,048,577** | **1,238/1,238** |

Unlike the linear `ftab`, which falls out of the suffix-array walk, the graph one
is built by **querying the finished index**: for each of the 1,048,576 prefixes,
walk the GFM backward and record the row range (`gfm.h:4993`). Which is why it
could only be attempted once the block itself was known correct.

`mapGLF` does not need `SideLocus`'s bit machinery to be reproduced, only its
semantics. It is an LF step over the BWT followed by a hop through the node
structure:

```
top = fchr[c] + occ_c(top)            bot = fchr[c] + occ_c(bot)
node_top = rank_M(top + 1) - 1        node_bot = rank_M(bot)
top = select_F(node_top + 1)          bot = select_F(node_bot + 1)
```

`F` marks where each node's incoming edges begin, so `select_F` converts a node
index back into a row; `M` marks where each node's rows begin, so `rank_M` does
the reverse. The C++ reaches the same values through per-side `F_locSave` and
`M_occSave` tallies and a backward walk over sides, which is the streaming form
of exactly this.

Two details decided it:

- **`occ` must exclude the `'Z'` row.** It is stored as `'A'` but was never
  counted when the block was written, so counting the stored bytes naively
  shifts every LF step after it.
- **`rank_M` is exclusive.** `rank_M(initFromRow_bit(x))` counts M bits in
  `[0, x)`, so `node_top = rank_M(top + 1) - 1` needs the exclusive prefix sum.
  Using an inclusive one left the search finding empty ranges almost everywhere
  — 8,242 of 1,048,577 entries matching, and the ones that did match were the
  degenerate propagated entries rather than real hits.

The diagnostic that made this quick was decoding *their* `ftab` back into ranges
and printing them beside ours. Theirs were width-1 and consecutive — 1-2, 2-3,
3-4 — which is what an F-column partition looks like; ours were empty. That said
"the search is failing", not "the search is subtly off", and pointed straight at
the rank rather than at the LF step.

## `.1.ht2` is now fully reproduced

Header, front end and `refnames` are unchanged from the linear case that rung 2
already emits byte-identically. Every graph-specific section — the gbwt block
with its F and M bitvectors and per-side tallies, `zOffs`, `fchr`, `ftab`,
`eftab` — now matches exactly, as does the whole `.2.ht2` offset array.

What remains for a complete graph index is `.5`/`.6`, the hierarchical local
indexes, which are themselves graph FM indexes over ~57 kb windows and so reuse
everything above.
