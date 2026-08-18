# P1 — reference bias: linear vs invented haplotypes vs real phasing

**Result for a European donor: both graph indexes remove about half the
reference bias, and real 1000 Genomes phasing adds essentially nothing over the
distributed index's invented haplotypes — +0.0007 against the +0.021 both
achieve over linear.** But the mechanism check is unambiguous, and it says where
the value actually is: **an index helps at a site if and only if it carries that
site's variant.**

## Design

Alternate-allele fraction at sites where the donor is known heterozygous.
Unbiased is 0.5; a linear reference sits below because alt-carrying reads pay a
mismatch penalty. Three arms, same 50,000,000 cDNA reads (SRR8551677), same
sites, fully paired — only sites at depth >= 20 in **all three** arms count, so
a difference between arms cannot come from which sites were covered.

| arm | index | variants | haplotypes |
|---|---|---|---|
| `linear` | `grch38` | none | — |
| `snp_invented` | `grch38_snp` (JHU) | 14,460,407 | greedy graph colouring |
| `snp_real` | `wg64` (ours) | 14,947,745 | real 1000 Genomes phasing |

Truth is GIAB HG001 v4.2.1 for NA12878: 1,950,933 het SNVs inside the
high-confidence regions, of which **8,689** are covered at depth >= 20 in all
three arms. Pileups use `-q 60`, so multi-mappers cannot let mapping ambiguity
masquerade as allelic imbalance.

The two indexes are near-matched on coverage — 93.83% and 94.76% of this donor's
alt alleles — so this is a test of **haplotype quality, not variant count**.

Alignment rates: 88.37% / 88.30% / 87.77%. The first reproduces H2's recorded
88.37% exactly, on the same reads and index, which is what validates the
pipeline before any new arm runs.

## Result

| arm | pooled alt frac | 95% CI | per-site mean | vs linear (paired) |
|---|---|---|---|---|
| linear | 0.46245 | 0.46172–0.46318 | 0.46921 | — |
| `snp_invented` | 0.46442 | 0.46372–0.46513 | 0.49030 | **+0.02108**, 5797 up / 363 down, z = 69.2 |
| `snp_real` | 0.46556 | 0.46485–0.46626 | 0.49098 | **+0.02176**, 5924 up / 374 down, z = 69.9 |

Both graph arms beat linear overwhelmingly. **They do not meaningfully differ
from each other**: 0.49098 against 0.49030, a thirtieth of the gap either one
opens over linear.

### The two estimators disagree tenfold, and both are right

Per-site mean improves by **+0.021**; the pooled fraction moves only **+0.002**.
Pooled weights by read mass and is dominated by a handful of very deep sites
where bias persists; per-site weights every site once and sees the graph
rescuing alt reads across many shallow-to-moderate sites. Reporting either alone
misrepresents the effect. H2's median-based 0.475 -> 0.500 was closer to the
per-site view, and its apparent cleanliness came from quantisation.

## The mechanism: carrying the variant is the whole effect

| site group | n | `snp_invented` vs linear | `snp_real` vs linear |
|---|---|---|---|
| variant **in** `snp_invented` | 8042 | +0.02262, z = 70.2 | +0.02285, z = 70.0 |
| variant **not in** `snp_invented` | 647 | +0.00194, z = −0.2, **p = 0.82** | **+0.00824, z = 7.4, p = 1.4e-13** |
| variant **in** `snp_real` | 8194 | +0.02202, z = 69.9 | +0.02304, z = 70.8 |
| variant **not in** `snp_real` | 495 | +0.00553, z = 1.5, p = 0.14 | **+0.00069, z = 0.1, p = 0.94** |

Read the diagonal. **Where an index does not carry the variant, that index gives
no benefit** — `snp_invented` p = 0.82 on its own missing sites, `snp_real`
p = 0.94 on its own. Where it does carry the variant, the benefit is enormous.
The effect is not a general property of graph alignment; it is variant-specific,
exactly as the mechanism predicts.

And the off-diagonal is the reason to build a better panel: at the 647 sites
`grch38_snp` misses, **our index still delivers a real gain (+0.008, z = 7.4)**,
because it carries variants the distributed index lacks.

## What this means, stated against the prediction made before the run

The prediction registered before the data existed was: *"Both graph arms beat
linear, and match each other... the distributed index already captures the
recoverable bias here, and real phasing earns its cost elsewhere — which makes
the ancestry arm the decisive experiment rather than a bonus."*

That is what happened. For a **European** donor, whose variation both panels
cover well (93.8% / 94.8%), real phasing is not worth 4 h 25 m and $28. The
honest headline is not "our index is better"; it is **"an index helps exactly
where it knows the variant, so the value of a real-phasing panel is entirely a
question of whose variants it contains."**

That relocates the claim rather than weakening it, and it makes the NA18502 arm
— a Yoruba donor, whose variation dbSNP common-variant sets cover least well —
the experiment that decides whether the whole-genome build was worth doing.

## Caveats

- **Both donors are IN the 1000 Genomes panel our index was built from**, so
  `snp_real`'s coverage of them is not a held-out measurement. The
  variant-not-carried rows are the uncircular control, and they behave correctly.
- `snp_real` aligns 0.6 points fewer reads than linear (87.77% vs 88.37%). More
  variants means more graph ambiguity; that cost is real and not yet quantified
  against its benefit.
- One donor, one cell type, one library. The ancestry question needs the cohort.

## Prediction for the Yoruba donor, registered before the alignment ran

Index coverage of each donor's heterozygous sites, on the same four chromosomes
(chr1, 6, 17, 19) so the contrast is not confounded by which chromosomes each
donor's truth list covers:

| index | NA12878 (CEU) | NA18502 (YRI) | disparity |
|---|---|---|---|
| `grch38_snp` (invented) | 93.69% | 86.83% | **6.87 pts** |
| `wg64` (real phasing) | 94.73% | 91.13% | **3.60 pts** |

NA18502 also carries **50.6% more heterozygous sites** than NA12878 over the
same chromosomes (552,456 against 366,890) — more variation, and a smaller
fraction of it known to either panel.

Real phasing roughly **halves the ancestry disparity**, and its advantage over
the distributed index is **four times larger for the Yoruba donor**: +4.30
points (91.13 vs 86.83) against +1.04 points (94.73 vs 93.69).

Given the mechanism established above — an index helps at a site if and only if
it carries that site's variant — this predicts:

> For NA18502, `snp_real` should beat `snp_invented` by roughly **four times**
> the margin seen for NA12878, with the gain concentrated in the ~24,000 sites
> where `grch38_snp` lacks the variant but ours carries it. If the two arms match
> again, the mechanism story is wrong: the coverage gap would not be translating
> into recovered reads.

Recorded before the Yoruba alignment was run, so the interpretation cannot be
chosen after seeing the number.
