# `.5.ht2` / `.6.ht2` — the hierarchical local indexes, layout accounted for

Every byte of both files parsed, on four indexes including the distributed
`grch38_snp`:

| index | local indexes | `.5` bytes | `.6` bytes |
|---|---|---|---|
| 200 bp | 1 graph | 8,557 | 64 |
| 509,431 bp | 10 graph | 372,925 | 131,142 |
| 900,000 bp | 18 graph | 662,997 | 232,192 |
| **`grch38_snp`** | **50,782 graph + 1,852 linear** | **2,154,980,035** | **786,599,058** |

All exact. The test has the same teeth as rung 1: every section size is derived,
so a wrong formula does not sum.

## Geometry

From `hier_idx_common.h`: windows of **57,344 bp** (`(1<<16) - (1<<13)`) at an
interval of **56,320**, so consecutive windows overlap by **1,024** bp. Local
indexes use `lineRate` 7, `offRate` 3 and `ftabChars` 6 — so a local `ftab` is
4,097 entries rather than 1,048,577.

`local_index_t` is **2 bytes**. That is the reason for the window size:
`local_max_gbwt` is `(1<<16) - (1<<11)`, so a window is chosen small enough that
its whole GBWT indexes in 16 bits. Every offset, `plen`, `rstart`, `fchr`, `ftab`
and `eftab` entry inside a local index is a `u16`, while `tidx`, `localOffset`
and `joinedOffset` stay full width.

```text
.5:  u32 sentinel
     index_t nLocalGFMs
     i32 lineRate, i32 unused, i32 offRate, i32 ftabChars, i32 flags
     per local index:
       full_index_t  tidx, localOffset, joinedOffset
       local_index_t len, gbwtLen, numNodes, eftabLen
       local_index_t nPat, plen[nPat]
       local_index_t nFrag, rstarts[nFrag*3]
       bytes         gbwt[gbwtTotLen]
       local_index_t numZOffs, zOffs[]
       local_index_t fchr[5], ftab[4097], eftab[eftabLen]
     u8 '\0'
.6:  u32 sentinel, then each local index's offs[] at local width
```

A local index with `len == 0` contributes only its header — `readIntoMemory`
returns early — which matters for windows that fall entirely inside an N run.

## The mixed structure, which was not obvious

`grch38_snp` has **1,852 linear local indexes among 50,782 graph ones**. A window
containing no variants is built through the *linear* FM path, not the graph one,
and the two differ in reserved bytes per side (4 `index_t` against 6) and in
packing (4 rows per byte against 2). So `.5` is not a homogeneous array of graph
indexes; the mode is decided per window and inferred on read from
`len + 1 == gbwtLen`.

Any emitter has to make that same per-window choice, and both branches are
already implemented — rung 2 emits the linear layout byte-identically and rung 3
now emits the graph one.

## What remains for emission

The pieces are all built: a local index is a GFM over a 57,344 bp window, so it
reuses the graph construction, doubling, `generateEdges`, row order, side
packing and `ftab` machinery already verified — at `u16` width with
`ftabChars` 6.

What is genuinely new is the **windowing**: splitting the reference's
`RefRecord`s per window (`hgfm.h:2128`), which is fiddly around N runs and
sequence boundaries, and deciding graph vs linear per window. That is the next
piece of work, and it is bookkeeping rather than algorithm.

---

# Emission: it runs, and one window is byte-identical

`ht2emit5` windows the reference, builds a GFM per window at `u16` width with
`offRate` 3 and `ftabChars` 6, and byte-compares against the real files.

```
200 bp -> 1 local index
  .5.ht2: BYTE-IDENTICAL (8,557 bytes)
  .6.ht2: BYTE-IDENTICAL (64 bytes)
```

First try, which is the payoff for the format work already being done: the local
index for a 200 bp reference is the *same graph* as the global one, so nothing
but the parameters and the width changed.

## Multi-window does not match yet — two distinct causes

| reference | our windows | theirs |
|---|---|---|
| 509,431 bp, no N | 10 | 10 |
| 900,000 bp, one 100 kb N run | **16** | **18** |

**1. Windows are laid out in REFERENCE coordinates, not joined-text ones.**
The 900,000 bp case is a 1,000,000 bp sequence with a 100,000 bp N run, and
`ceil(1000000 / 56320) = 18` while `ceil(900000 / 56320) = 16`. The C++ accrues
`rec.off + rec.len` per `RefRecord` — ambiguous characters included — and asserts
exactly this count (`hgfm.h:2172`). So a window can span or even fall entirely
inside an N run, which is what the `len == 0` early return in
`LocalGFM::readIntoMemory` exists for.

**2. Our per-window graph is slightly too large.** On the no-N reference, where
the window count is right, the first difference is at byte 42 — `gbwtLen` of the
first local index. Totals: ours ~527,300 GBWT rows against 525,146, roughly 200
extra per window against ~190 variants per window, so on the order of one extra
edge per variant. The likely cause is the haplotype admission rule at a window
edge: `build_range` takes haplotypes with `left >= a && right < b`, and the C++
clips its `RefRecord`s to `local_index_interval` rather than `local_index_size`
in places (`hgfm.h:2185`), so a variant in the 1,024 bp overlap may belong to
only one of the two windows rather than both.

Both are bookkeeping in the windowing, not the index format: every byte of a
correctly-scoped window already emits exactly.
