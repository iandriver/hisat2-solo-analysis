# E1 — how long are the queries against the global index?

The question decides whether an order-bounded (GCSA2-style) index is a free
optimisation or a semantic change. HISAT2 searches the global graph FM index
with `ftabLoHi` followed by repeated `mapLF`/`mapGLF` — a maximal exact match
whose length is bounded only by the read, since `maxHitLen` defaults to
`INDEX_MAX`. An order-k index is exact for queries of length <= k and may return
false positives beyond it, so everything turns on the distribution.

Measured by instrumenting the two global entry points — `partialSearch`
(`hi_aligner.h:6361`, the main seed search) and `globalGFMSearch` (`:6606`,
spliced extension) — and, separately, `localGFMSearch` (`:6751`) as a contrast.
Instrumentation lives in a throwaway copy of the tree, never in the fork.
200,000 reads each against `grch38_snp`, `-p 1`.

## Result

| reads | scope | n | max | mean | p50 | p99 | p99.9 | p99.99 |
|---|---|---|---|---|---|---|---|---|
| 150 bp | **global** | 411,021 | **131** | 24.1 | 18 | 91 | 121 | 130 |
| 150 bp | local | 1,295,774 | 130 | 9.2 | 9 | 27 | 56 | 88 |
| 91 bp | **global** | 266,879 | **86** | 22.3 | 17 | 62 | 75 | 83 |
| 91 bp | local | 959,483 | 70 | 9.1 | 8 | 25 | 39 | 55 |

Global queries reaching 32 bp: 20.0% (150 bp reads) / 19.2% (91 bp).
Reaching 64 bp: 3.6% / 0.85%. Reaching 128 bp: 0.027% / 0.000%.

**The maximum is bounded by read length** — 131 on 150 bp reads, 86 on 91 bp
reads — which is what a maximal exact match must do. Nothing queries deeper than
the read.

## What that means, against the construction curve

`generation` in `PathGraph` means "sorted by paths of length 2^generation"
(`gbwt_graph.h:1694`), and construction runs `while(!isSorted())` — to **full**
disambiguation. Cross-referencing the whole-genome curve from `E0`:

| generation | order = 2^g | nodes | % of 2^32 |
|---|---|---|---|
| 8 | **256** | 3,390,374,563 | **78.9%** — fits |
| 9 | **512** | 3,595,956,954 | **83.7%** — fits |
| 10 | 1024 | overflow | — |

**Stopping the doubling at generation 8 or 9 keeps the whole human graph index
inside the 32-bit ceiling, and is exactly equivalent for every query HISAT2
actually issues.** Order 256 exceeds the observed maximum by 2x; order 512 by
nearly 4x.

For scale on how much precision is being bought and never used: chr22 required
14 generations to fully sort — order 8192 — and the whole genome needs more than
10. The index is disambiguated to a resolution roughly two orders of magnitude
beyond the longest query.

## Consequence: this reorders the whole plan

The builder design assumed the fix was 64-bit ids plus external memory, because
the node count exceeds 2^32. That assumption came from measuring the *converged*
count. It is now clear the convergence itself is the thing that is unnecessary.

Ranked by leverage:

1. **Bound the order.** Stop at generation 9 (order 512, 83.7% of ceiling) or 8
   (order 256, 78.9%). A human whole-genome graph index becomes buildable at
   **32 bits, in ~176 GB** — the measured peak RSS of the 32-bit run before it
   overflowed. No new integer width, no external memory, no format change to the
   node ids.
2. **External memory** becomes an optimisation for build-machine cost rather
   than a precondition for the index existing at all.
3. **64-bit** becomes unnecessary for human, though still required for larger
   collections (pangenomes, many haplotypes).

## What still has to be proved

This is a measurement plus an inference, and the inference has an implementation
cost that has not been paid:

- **Downstream construction assumes a fully sorted `PathGraph`.** `generateEdges`
  and the GFM build run after `while(!isSorted())` completes. Stopping early
  leaves nodes that share a 2^g prefix merged, and the rest of the pipeline has
  to handle that. This is exactly what GCSA2 implements, so it is known-possible,
  but it is real work and it is where a naive "just break out of the loop" would
  produce a wrong index rather than a bounded one.
- **The equivalence claim is for queries <= 2^g.** It holds for these read
  lengths, with margin. It would need re-checking before anyone ran 300 bp reads
  against an order-256 index — order 512 covers those, and the curve says order
  512 fits.
- **Only two read lengths and one sample type were measured.** The shape is
  strongly constrained by "a MEM cannot exceed its read", so this is unlikely to
  move much, but longer-read RNA-seq would want its own check.

## Reproducing

`qlen_probe.h` plus hooks at the seven global sites and three local ones. The
patch is deliberately not applied to the fork; it lives here and is applied to a
scratch copy.

```
HISAT2_QLEN_OUT=qlen.tsv hisat2 -x grch38_snp/genome_snp -U reads.fq -p 1 --no-unal -S /dev/null
```

`-p 1` matters: the histogram counters are deliberately unsynchronised, since
adding atomics to the hot path would distort exactly what is being measured.
