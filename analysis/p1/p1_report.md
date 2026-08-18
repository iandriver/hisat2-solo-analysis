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

## Result for the Yoruba donor

50,000,000 reads (SRR8551676), 564,631 known het sites on chr1/6/17/19, **3,874**
covered at depth >= 20 in all three arms. Alignment rates 87.87% / 87.76% /
87.38%, again reproducing H2's recorded 87.87%.

| donor | arm | pooled | dist. from 0.5 | per-site | dist. from 0.5 | bias removed |
|---|---|---|---|---|---|---|
| CEU | linear | 0.46245 | 0.0376 | 0.46921 | 0.0308 | — |
| | `snp_invented` | 0.46442 | 0.0356 | 0.49030 | 0.0097 | 68.5% |
| | `snp_real` | 0.46556 | 0.0344 | 0.49098 | 0.0090 | **70.7%** |
| YRI | linear | 0.40127 | 0.0987 | 0.45236 | 0.0476 | — |
| | `snp_invented` | 0.50254 | 0.0025 | 0.47705 | 0.0229 | 51.8% |
| | `snp_real` | 0.52915 | 0.0292 | 0.47909 | 0.0209 | **56.1%** |

### The prediction held, at 3x rather than 4x

Real phasing's per-site advantage over invented haplotypes is **+0.00204** for
the Yoruba donor against **+0.00068** for the European — **3.0x**, against the
~4x registered before the run. Direction and order of magnitude confirmed.

The mechanism prediction held exactly. At the 507 sites where `grch38_snp` lacks
the variant, `snp_invented` does nothing (+0.002, p = 0.29) while `snp_real`
delivers **+0.010, z = 3.7, p = 1.8e-4** — the gain lands precisely where the
coverage gap said it would.

### Two findings that were not predicted, and matter more

**1. The graph narrows the ancestry gap but does not close it.** Linear bias is
1.5x worse for the Yoruba donor (0.0476 against 0.0308). After the best graph
arm the residual is **2.3x worse** (0.0209 against 0.0090), because the graph
removes only 56% of the bias for the Yoruba donor against 71% for the European.
**Variant-aware alignment reduces reference bias for everyone and reduces the
ancestry disparity in it, but a donor from an underrepresented population is
still left with more than twice the residual bias.** That is the honest headline
and it is not the flattering one.

**2. Pooled and per-site disagree about which index is better, and pooled says
`snp_invented`.** For the Yoruba donor `snp_invented` lands at 0.5025 — almost
exactly unbiased — while `snp_real` overshoots to 0.5292. By "closest to 0.5" on
the pooled estimator, the distributed index wins.

This is not a small discrepancy and it cannot be waved away. The reason to
prefer per-site here is that **the top 1% of sites carry 42% of all reads**
(57% for the European donor), and those are the most highly expressed genes,
where true allele-specific expression is a real biological signal rather than a
mapping artifact. A read-mass-weighted statistic therefore measures ASE and
mapping bias together, which is not the quantity under test. Excluding the top
1% of sites moves `snp_real` from 0.5292 only to 0.5159, so the overshoot is not
a handful of loci — it is broad, and it is a real property of the pooled view.

**The conservative reading is that `snp_real` beats `snp_invented` on the
estimator that isolates mapping bias, and loses on the estimator that does not,
and that a single donor cannot settle which matters more.**

### A cost, measured

At sites where `snp_real` does *not* carry the variant, it is slightly **worse**
than linear for the Yoruba donor: −0.0029 per site, 52 up / 80 down, z = −2.4,
p = 0.015. Carrying 14.9M variants buys recall where the variant is known and
costs a little precision where it is not. `snp_real` also aligns the fewest
reads of the three arms in both donors (87.77% and 87.38% against linear's
88.37% and 87.87%).

## Standing caveats

- **Truth quality differs between the donors.** NA12878 uses GIAB HG001 v4.2.1;
  NA18502 uses 1000 Genomes panel calls, which carry more genotype error. The
  diagnostic argues against this driving the overshoot — sites with alt fraction
  above 0.9 are 2.17% for the Yoruba donor against 1.50% for the European, not
  the spike a hom-alt miscall would produce — but it is not eliminated.
- **Both donors are IN the panel our index was built from.** The
  variant-not-carried rows are the uncircular control and behave correctly.
- **n = 2 donors.** The ancestry claim rests on one comparison. The 3x ratio and
  the 2.3x residual gap are point estimates from a single pair, not a cohort.
