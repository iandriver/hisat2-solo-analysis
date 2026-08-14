Following up with a result that argues **against** part of my own report, since it seems more useful to have it here than not.

I tested whether the retention difference actually changes alignment outcomes, and on this test it does not.

## Method

Matched read pairs at 4,000 sampled SNVs on the 3 Mb chr1 slice above. For each site, two 100 bp reads identical except for the base at the variant — one reference allele, one alternate — with the same extra mismatches at the same offsets, so the pair differs in exactly one position and any difference in alignment is attributable to the variant alone. Verified: 4,000 of 4,000 pairs differ at exactly one base. A read counts only if it aligns at its true locus.

Three indexes over the same reference and the same `.snp`/`.haplotype`: one linear, one built with the current halving backoff (53.1% retention), one with the binary-search backoff (64.2%).

## Result

ALT fraction — the share of aligned reads at a site carrying the alternate allele. 0.5 is unbiased.

| extra mismatches | linear | graph, 53.1% retained | graph, 64.2% retained |
|---|---|---|---|
| 3 | 0.0285 | 0.5084 | 0.5084 |
| 5 | 0.0000 | 0.4848 | 0.4857 |
| 7 | 0.0000 | 0.0000 | 0.0000 |

**Graph vs linear is large**: 2.9% of alt-carrying reads align against a linear index, 50.8% against the graph.

**The retention difference is not detectable**: 0.5084 vs 0.5084, a net gain of 3 sites out of 4,000. Stratifying by local variant density does not surface it either — even at 11+ SNVs per 200 bp the figures are 0.5116 vs 0.5127. At 5 extra mismatches almost nothing aligns (16 vs 17 reads) and at 7 nothing does, so there is no harsher regime available to separate them: default scoring allows roughly 3-4 mismatches per 100 bp.

## Caveat on power, since this is a null

With 4,000 sites the standard error on a proportion is about 0.8 points. The best estimate I have of the true effect size is from mouse chr19, where losing 78% of window-level variant instances cost 1.9 points of overall alignment rate — scaling that, an 11-point retention gain is worth roughly 0.25 points, about 3x below this test's noise floor. So this is "no detectable effect, and the test could not have detected the expected one", not "no effect". Resolving it properly needs ~10^5-10^6 reads across a range of retention levels.

## What I think this does and does not change

- **The reporting problem stands unchanged.** Variants are still silently dropped and `hisat2-inspect --snp` still reports the full set that was supplied. That is true regardless of how much the loss costs in practice, and it is the part I would most like to see fixed.
- **The binary-search change is justified on retention grounds only.** It puts more of the user's variants into the index; it does not have a demonstrated downstream benefit, and I am no longer claiming one. If that makes it less interesting as a patch, that is a fair reading.
- A plausible mechanism for the null, though I have not proven it: the global ALT list retains every variant and only the local indexes lose them. Since local indexes handle extension, reads that fit within the mismatch budget can be served by the global graph alone.
