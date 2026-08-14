Ran the powered version I said was needed in the previous comment: a dose-response across retention levels on whole chr1, 200,000 sampled sites (400,000 matched reads) instead of 4,000.

Three indexes over the same reference and the same `.snp`/`.haplotype` (dbSNP b157 common, 1,985,316 variant instances), differing only in how the builder backs off when a local graph exceeds the edge budget — quarter-decimation, the current halving, and a binary search for the largest subset that fits.

| index | retention | REF aligned | ALT aligned | ALT fraction | 95% CI |
|---|---|---|---|---|---|
| linear | — | 171,207 | 5,043 | 0.02861 | ±0.00078 |
| quarter | 75.1% | 170,026 | 170,574 | 0.50080 | ±0.00168 |
| halving (current) | 80.9% | 170,021 | 170,590 | 0.50084 | ±0.00168 |
| binary search | 87.5% | 170,018 | 170,608 | 0.50087 | ±0.00168 |

Paired on the same sites (McNemar, ALT aligned at its true locus or not):

| contrast | lost | gained | net | chi2 | |
|---|---|---|---|---|---|
| 75.1% -> 80.9% | 23 | 39 | +16 | 3.63 | n.s. |
| 80.9% -> 87.5% | 54 | 72 | +18 | 2.29 | n.s. |
| 75.1% -> 87.5% | 51 | 85 | +34 | 8.01 | p<0.05 |

## Conclusion

The effect exists and is positive in every contrast, but it is negligible in size. Across the full 12.4-point retention span the ALT fraction moves 0.50080 to 0.50087 — seven ten-thousandths of a point, or 34 sites out of 200,000. The extreme contrast reaches significance only because n is large.

For scale, the linear-to-graph difference on the same reads is +0.472. The retention difference is +0.00007, roughly 6,700x smaller.

At 75.1% retention — a quarter of variant instances missing from the local indexes — the ALT fraction is already 0.5008, i.e. unbiased. This supports the mechanism I guessed at earlier: the global ALT list keeps every variant, and for reads inside the mismatch budget the global graph is doing essentially all of the work. Local-index retention is close to irrelevant for this kind of read.

## What I am now claiming, and what I am not

- **The reporting problem is the substantive part of this issue.** Variants are dropped and `hisat2-inspect --snp` reports the full supplied set regardless. Whatever the practical cost, a user cannot currently find out what their index contains, and I would still like to see the count reported.
- **The binary-search patch improves retention and nothing else that I can measure.** It operates over 80.9% to 87.5%, which this experiment shows is a flat part of the curve. If you would rather have only the one-line summary and not the algorithm change, that is a reasonable call and I would not argue with it.
- **One regime remains untested.** The measurement that originally motivated my instrumentation was mouse chr19 at roughly 22% retention, where alignment rate moved 97.6% to 95.7% — a real effect far below anything tested here. So the curve may be flat across 75-87.5% and turn sharply somewhere below it. I have not looked for that threshold.

Method detail, in case it matters: each site contributes two 100 bp reads identical except for the base at the variant, with the same three extra mismatches at the same offsets, so the pair differs in exactly one position. The extra mismatches are deliberate — a lone SNV sits inside the default mismatch budget, so without them the test measures nothing. A read counts only if it aligns at its true locus.
