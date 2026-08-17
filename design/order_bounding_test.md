# E2 — bounding the doubling order: the build wins are real, the index is unusable

**Result: truncating the doubling loop is not sufficient. It produces an index
the aligner cannot query — 100 reads unfinished after 120 s, against 1 s for the
stock index — even though every build-side number improves substantially.**

This is the risk flagged at the end of `query_lengths.md` ("a naive `break` out
of the doubling loop would produce a wrong index, not a bounded one"), tested
rather than assumed. It is confirmed.

## The patch

`PathGraph`'s constructor runs `while(!isSorted()) lateGeneration();` — full
disambiguation. Added an env-var cap, gated on graph size so that local indexes
(which also build `PathGraph`s, `hgfm.h:1921`, and are <= 63,488 edges) still
fully sort:

```cpp
index_t maxGen_ = 0;
{ const char* e_ = getenv("HISAT2_MAX_GENERATION");
  if(e_ != NULL && temp_nodes > 1000000) maxGen_ = (index_t)atoi(e_); }
while(!isSorted()) {
    if(maxGen_ > 0 && generation >= maxGen_) { /* report and */ break; }
    lateGeneration();
}
```

Kept out of the fork; the patch lives in `e2/`.

## Build side: everything improves

chr22, 50,818,468 bp, 224,119 variants, 254,497 haplotypes.

| | stock (full) | order 256 (gen 8) |
|---|---|---|
| generations | 14 | 9 |
| final nodes | 57,776,238 | ~44,386,691 (**-23%**) |
| peak RSS | 11.55 GB | 9.09 GB (**-21%**) |
| **index on disk** | **111.1 MB** | **59.1 MB (-47%)** |
| nodes left unsorted | 0 | 265,616 (0.6%) |
| build exit | 0 | 0 |

The build completes cleanly and the index shrinks by nearly half. The
*opportunity* identified in E1 is real: chr22 needs order 8192 to fully sort,
and stopping at order 256 discards a great deal of structure that no query uses.

## Query side: it does not work

40,074 simulated 150 bp reads (exact substrings of chr22, true position encoded
in the read name) and the same reads with 2 mismatches.

| index | 100 reads | 40,074 reads |
|---|---|---|
| stock | **1 s**, 104 records | completes; 100.00% / 98.88% aligned |
| order 256 | **>120 s, 0 records** | >570 s at 99.4% CPU, **0 records, 0-byte SAM** |

Not even the SAM header was written, so it is stuck inside the first read's
search rather than producing output slowly. This is degenerate behaviour, not a
throughput regression.

**Why**: after an early stop, `ranks < nodes.size()` — nodes sharing a 2^g
prefix collapse to the same rank. The GFM is built over ranks, and the search
resolves a BWT range to node positions assuming that correspondence is sound.
With merged nodes it is not, and the range walk does not terminate usefully.
Handling that is exactly what GCSA2 implements and what the stock pipeline,
which runs `generateEdges` and the GFM build only after `while(!isSorted())`,
has never had to do.

## Baseline sanity

The stock index on the simulated set: 97.123% of exact reads and 94.470% of
2-mismatch reads placed within 1 bp of their true position. The residual is
expected — chr22 contains repeats where a 150 bp exact substring genuinely has
more than one equally good placement, so "wrong" here includes reads that are
legitimately ambiguous. It is recorded as the reference point any future
order-bounded index has to match, not as a defect.

## What this changes

E1 concluded that order-bounding was "the primary lever" and demoted 64-bit ids
and external memory to optimisations. **That conclusion was premature.** The
correct statement now:

- **The over-precision is real and worth attacking** — 47% of the chr22 index is
  structure no query reaches, and on the whole genome the same effect is what
  pushes generation 10 past the 2^32 ceiling.
- **Capturing it requires work in the GFM construction and the search**, not a
  bounded loop. Concretely: represent merged nodes explicitly, and make `locate`
  and the extension path handle a rank that maps to several graph positions.
- **So the builder plan's ordering stands as originally written.** 64-bit ids
  plus external-memory doubling remain the route to a whole-genome index that
  works with today's aligner, because they change only where data lives, not
  what the index means. Order-bounding is a larger, later change with a bigger
  payoff.

The cheap experiment that would have been misleading here is the one that stops
at "the build got 47% smaller and exited 0". It did. The index is still unusable.

## Reproducing

```
HISAT2_MAX_GENERATION=8 hisat2-build -p 8 --snp 22.snp --haplotype 22.haplotype chr22.fa idx_o256/c22
hisat2 -x idx_o256/c22 -U tiny.fq -p 1 --seed 0 --reorder -S out.sam    # will not finish
```
