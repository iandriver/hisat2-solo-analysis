# Does variant retention affect alignment-stage reference bias?

Two experiments. The first (A) tests the premise cheaply and fails to detect an
effect; the second (C) is powered to resolve the expected effect size and finds
it real but negligible.

Both exist because the binary-search backoff in `hgfm.h` was written on the
assumption that variants dropped from local indexes cost alignment sensitivity.
That assumption was never tested before the patch was committed.

## Design

Reference bias is measured with matched read pairs. For each sampled SNV, two
100 bp reads identical except for the base at the variant — one reference
allele, one alternate — carrying the **same** three extra mismatches at the
**same** offsets. The pair differs in exactly one position, so any difference in
whether the two align is attributable to the variant and nothing else.

ALT fraction = ALT-aligned / (ALT-aligned + REF-aligned). **0.5 is unbiased by
construction.** A read counts only if it aligns at its true locus.

The three extra mismatches are load-bearing. A lone SNV sits inside HISAT2's
default mismatch budget, so a bias test without them measures almost nothing.
With the budget partly spent, knowing the variant is what decides whether the
alt-carrying read aligns at all. At five extra mismatches almost nothing aligns
(16 reads of 10,000) and at seven nothing does, so the usable window is narrow:
default scoring allows roughly 3-4 mismatches per 100 bp.

## A — 3 Mb chr1 slice, 4,000 sites

| extra mismatches | linear | graph 53.1% | graph 64.2% |
|---|---|---|---|
| 3 | 0.0285 | 0.5084 | 0.5084 |
| 5 | 0.0000 | 0.4848 | 0.4857 |
| 7 | 0.0000 | 0.0000 | 0.0000 |

No detectable difference between retention levels (net +3 sites of 4,000).
Stratifying by local SNV density did not surface one either — at 11+ SNVs per
200 bp, 0.5116 vs 0.5127.

Underpowered by design: standard error ~0.8 points against an expected effect of
~0.25 points, extrapolated from the mouse chr19 datapoint where a 78% instance
loss cost 1.9 points of alignment rate.

## C — whole chr1, 200,000 sites, three retention levels

Same reference, same `.snp`/`.haplotype` (dbSNP b157 common, 1,985,316 variant
instances). The indexes differ only in how the builder backs off when a local
graph exceeds the edge budget.

| index | retention | REF aln | ALT aln | ALT fraction | 95% CI |
|---|---|---|---|---|---|
| linear | — | 171,207 | 5,043 | 0.02861 | ±0.00078 |
| quarter-decimation | 75.1% | 170,026 | 170,574 | 0.50080 | ±0.00168 |
| halving (stock) | 80.9% | 170,021 | 170,590 | 0.50084 | ±0.00168 |
| binary search | 87.5% | 170,018 | 170,608 | 0.50087 | ±0.00168 |

Paired on the same sites (McNemar):

| contrast | lost | gained | net | chi2 | |
|---|---|---|---|---|---|
| linear -> 75.1% | 77 | 165,608 | +165,531 | 165,375 | p<0.001 |
| 75.1% -> 80.9% | 23 | 39 | +16 | 3.63 | n.s. |
| 80.9% -> 87.5% | 54 | 72 | +18 | 2.29 | n.s. |
| 75.1% -> 87.5% | 51 | 85 | +34 | 8.01 | p<0.05 |

## Conclusions

**The variant-aware claim is confirmed emphatically.** Linear alignment recovers
2.9% of alt-carrying reads; the graph recovers 50.8%. On the paired test that is
165,531 sites out of 200,000, chi2 = 165,375. Since the REF and ALT reads differ
in exactly one base, nothing else can explain it.

**Retention barely matters over the range the patch operates on.** Across the
full 12.4-point span the ALT fraction moves 0.50080 to 0.50087 — 34 sites out of
200,000. The direction is positive in all three contrasts so it is not noise,
but the graph-vs-linear effect is roughly 6,700x larger.

**Mechanism, now measured rather than inferred.** At 75.1% retention, with a
quarter of variant instances missing from the local indexes, the ALT fraction is
already 0.5008 — unbiased. The global ALT list keeps every variant, and for
reads inside the mismatch budget the global graph does essentially all the work.
Local-index retention is close to irrelevant for this class of read. This is the
same fact that makes the loss undetectable from outside: `hisat2-inspect --snp`
reports all 29,435 variants of a slice where only 53.1% reached the local
indexes.

**Consequence for the patch.** The binary search moves retention 80.9% to 87.5%,
which is squarely inside the flat part of the curve. Its justification is that
users get the variants they supplied, and that silent dropping should be
reported — not that alignment improves. That is what was reported upstream in
DaehwanKimLab/hisat2#473.

**Untested regime.** The mouse chr19 measurement that motivated the original
instrumentation was at roughly 22% retention and did show a real effect
(alignment rate 97.6% -> 95.7%). The curve may be flat across 75-87.5% and turn
sharply below it. No threshold search was run.

## Reproducing

```
mkbias.py <genome.fa> <variants.snp> <reads.fa> <n_sites> <n_mismatch> <seed>
hisat2 -x IDX -f -U reads.fa --no-spliced-alignment --no-unal -S out.sam
analyze_bias.py LINEAR=a.sam QUARTER=b.sam BASE=c.sam BSEARCH=d.sam
```

The three backoff variants are built from `hgfm.h`: stock halving, the same with
`>> 1` replaced by `>> 2` (quarter), and the binary search on
`upstream/modernize` at ba686ac. SAM files and indexes are not kept here — they
regenerate from the above in a few minutes.
