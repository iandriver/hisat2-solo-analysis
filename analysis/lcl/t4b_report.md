# T4b: does variant-aware alignment cut the scRNA-seq genotype error rate?

**Answer: no, not meaningfully — the effect is ~0.1 percentage point and its
sign depends on the genotype caller.** But the experiment was worth running,
because it closes T2's standing caveat and it exposes a cost of the graph index
we had previously called negligible.

## Design

Yao & Gazal (doi:10.1101/2025.03.25.645175) call genotypes from scRNA-seq reads
aligned with **STARsolo — a linear reference** — and measure an 8% genotype
error rate against SNP-array truth, noting ADMIXTURE's ancestry error rises
exponentially above ~6%. Our measured 8-10% relative undercount of the alternate
allele points the right way to explain some of that.

**GSE126321** supplies the test: two lymphoblastoid lines, deeply sequenced,
raw reads open in SRA, and — unlike the multiplexed alternatives — one donor per
library, so no genetic demultiplexing.

| line | 1000G ID | population | truth |
|---|---|---|---|
| GM12878 | NA12878 | CEU / European | **GIAB HG001 v4.2.1** |
| GM18502 | NA18502 | YRI / African | 1000 Genomes (not used here) |

50,000,000 cDNA reads (150 bp) per line, aligned to `grch38_snp` and `grch38`.
Genotyping is pseudobulk, so barcodes are irrelevant and only the cDNA mate was
downloaded. Alignment rates 88.30% / 88.37% (GM12878) and 87.76% / 87.87%
(GM18502).

Genotypes are called **by allele fraction rather than with GATK**, deliberately:
the question is what the aligner contributes, and a threshold caller adds no
model of its own. Sites are common dbSNP SNVs inside the GIAB high-confidence
regions; positions the benchmark VCF touches with an indel or non-PASS record
are dropped rather than scored against a wrong truth.

## Result: a wash, and threshold-dependent

GM12878, 54,754 sites callable at depth >=20 by both (42,063 hom-ref, 7,991
het, 4,700 hom-alt):

| truth | n | graph wrong | linear wrong |
|---|---|---|---|
| hom-ref | 42,063 | 126 (**0.30%**) | 55 (**0.13%**) |
| het | 7,991 | 331 (**4.14%**) | 359 (**4.49%**) |
| hom-alt | 4,700 | 43 (0.91%) | 53 (1.13%) |
| **all** | 54,754 | 500 (**0.91%**) | 467 (**0.85%**) |

The graph is better where the donor carries an alternate allele and worse where
it does not. Hom-ref sites outnumber het sites 5:1, so the two roughly cancel.

**And the sign flips with the calling threshold:**

| het calling rule | graph | linear | difference |
|---|---|---|---|
| AF >= 0.10 | 0.91% | 0.85% | **+0.06** (graph worse) |
| AF >= 0.15 | 1.00% | 1.06% | **-0.06** (graph better) |
| AF >= 0.20 | 1.22% | 1.34% | **-0.12** (graph better) |

Requiring a minimum alternate-read count changes nothing, because at depth >=20
an AF of 0.10 already implies two alternate reads.

**At truth-het sites the graph is better under every rule**, and increasingly so
as the threshold tightens:

| het calling rule | graph | linear | difference |
|---|---|---|---|
| AF >= 0.10 | 4.14% | 4.49% | -0.35 |
| AF >= 0.15 | 5.47% | 6.21% | -0.74 |
| AF >= 0.20 | 7.31% | 8.25% | -0.94 |

**Our ~1% error rate is not comparable to the paper's 8%** and this is not a
reduction of it. Different caller (threshold versus GATK de novo), different
site set (known dbSNP inside high-confidence regions), different depth. The
honest statement is that *at the aligner level* there is no lever here worth
0.1 of a percentage point either way.

## What the experiment did establish

### 1. T2's central measurement, validated against real genotypes

This was the main reason to run it. At sites GIAB confirms are heterozygous:

| | median ALT fraction |
|---|---|
| graph | **0.5000** |
| linear | **0.4762** (0.4735 at depth >=10) |

Identical to what T2 measured from data-derived het sites in two other donors
(0.4500 and 0.4615). **The circularity caveat is closed** — the reference-allele
undercount is real when the heterozygous sites are supplied by an external
benchmark rather than called from the same reads.

### 2. Allele dropout is real, and the graph reduces it

Truth-het sites miscalled homozygous reference:

| depth | graph | linear |
|---|---|---|
| >=10 | 2.80% | 3.40% |
| >=20 | **2.47%** | **2.93%** |

About one in six dropouts is recovered. This is the mechanism working as
designed and it is the same effect T2 measures, expressed as a genotype error.

### 3. A cost we had wrongly called negligible

Truth hom-ref sites miscalled heterozygous:

| depth | graph | linear |
|---|---|---|
| >=10 | 0.359% | 0.165% |
| >=20 | **0.247%** | **0.081%** |

**The graph makes three times as many false heterozygous calls.** At the 73
graph-only false-het sites (depth >=20) the median allele fraction is 0.143
under the graph and 0.000 under linear, with 7% more depth — the graph is
recruiting a handful of alternate-carrying reads at positions where the donor
has no alternate allele.

Earlier reports called the graph's spurious ALT signal "negligible for ASE at
these depths" on the basis of a mean ALT fraction of 0.00042 versus 0.00008.
That was true as an average and **wrong as a conclusion**: concentrated at a
small number of sites, it crosses genotype-calling thresholds. This is the
mirror image of the benefit — alternate paths let reads align that should not.

## Ancestry arm: bias versus distance from the reference

The same two lines, same protocol, same 50M-read depth, same indexes. Het sites
called from the data identically for both, so this is a sample-to-sample
comparison needing no external genotypes.

| | GM12878 (CEU) | GM18502 (YRI) |
|---|---|---|
| sites callable | 60,239 | 62,774 |
| heterozygous sites | 8,457 | **10,856** |
| het as % of callable | 14.04% | **17.29%** |
| non-ref alleles per 1k sites | 306.3 | **359.6** |
| ALT fraction, graph | 0.5000 | 0.5000 |
| ALT fraction, linear | 0.4750 | 0.4706 |
| **linear's deficit from 0.5** | **0.0250** | **0.0294** |
| depth ratio, 1 alt copy | 1.0339 | 1.0335 |
| depth ratio, 2 alt copies | 1.0673 | 1.0678 |

The African-ancestry line is **1.23x more heterozygous** — further from a
reference that is predominantly of European ancestry — and linear alignment's
per-site deficit is **1.18x larger**. Combined, the total burden of reference
bias is roughly **1.45x greater**. The graph sits at exactly 0.5000 for both.

Note the per-site depth ratios are essentially identical (1.034 vs 1.034;
1.067 vs 1.068). The extra burden comes mostly from *more affected sites*, not
from each site being worse.

**This is n=2 individuals, one per population.** It is directionally consistent
with the hypothesis and the within-sample estimates are precise, but two cell
lines cannot establish a trend across ancestries. It is an illustration, not a
demonstration.

## Verdict

- **T4b as posed: fails.** Variant-aware alignment does not meaningfully change
  the scRNA-seq genotype error rate, and does not offer a route to the 8% the
  paper reports.
- **T2 is stronger than before.** Its core measurement now holds against GIAB
  truth in a third donor, and the reference-allele undercount is confirmed at
  0.4762 versus an exact 0.5000.
- **The graph has a real cost**, previously understated: 3x the false
  heterozygous calls, from spurious alternate-allele recruitment.
- **The ancestry direction is right but modest**, and n=2.

## Caveats

- One donor with benchmark truth. GIAB high-confidence regions exclude the
  hardest parts of the genome, including much of the MHC — so the sites where
  the graph helps most for expression are largely *absent* from this test.
- LCLs are not primary tissue, and the library is not 10x-standard: 150 bp
  single-end cDNA reads at very high depth.
- The threshold caller is not GATK. It isolates the aligner, but it also means
  these error rates should not be compared with any published pipeline's.
- The 1000 Genomes truth for GM18502 was not used; the ancestry arm calls its
  het sites from data, identically for both samples.
