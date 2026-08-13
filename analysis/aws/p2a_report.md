# Stage 2a — chr1 probe on GENCODE v50 + dbSNP b157

`PROBE_PASS`. r7i.8xlarge (32 vCPU / 247 GB), us-east-1, 109 min on-demand ≈ **$3.85**.
Instance `i-0480589ef97d15d0a` stopped itself; 500 GB volume preserved.
Raw trace: `s3://rustar-bench/hisat2-p2a/pass2/trace.txt`.

## Reference

| | source | chr1 result |
|---|---|---|
| genome | GENCODE v50 `GRCh38.primary_assembly.genome.fa` | 248,956,422 bp, **not soft-masked** |
| annotation | GENCODE v50 `primary_assembly.annotation.gtf` | 1,119,185 chr1 lines |
| variants | dbSNP **b157** (`GCF_000001405.40`) | 93,124,941 chr1 records |
| phased panel | 1000G 30x, 20220422 SNV+INDEL+SV | 1,273,092 @ MAF≥1%, 3,202 samples |

`incompatible=0` on all three extractions: GENCODE v50 and dbSNP b157 agree base-for-base
across 1.24M sites, so the assembly guard is real rather than permanently-firing.

## "Common" is no longer a defined set

dbSNP dropped the `COMMON` flag after b151, so it must be derived from `FREQ=`.
The choice of frequency source swings density **3.6×**:

| rule | chr1 | per bp | → genome | vs shipped 12.3M |
|---|---|---|---|---|
| any source | 4,477,685 | 1/55.6 | ~56M | 4.5× |
| **big-4** (1000G_30X, gnomAD gen+ex, TOPMED) | 1,241,644 | 1/200 | ~15.5M | **1.26×** |
| 1000G_30X only | 1,073,447 | 1/232 | ~13.4M | 1.09× |

"Any source" is not usable: on a test window, 36% of its survivors rested solely on
SGDP_PRJ (~280 samples, where one heterozygote reads as 50%), plus 1,206 on KOREAN
alone and 240 on Siberian.

## The headline: a human graph index silently discards ~19% of its variants

| build | instances | retained | over-budget graphs | per 1k inst. |
|---|---|---|---|---|
| 1000G phased (real haplotypes) | 1,204,387 | **92.4%** | 472 | 0.392 |
| dbSNP 1000G_30X | 1,731,945 | **86.7%** | 692 | 0.400 |
| dbSNP big-4 | 1,985,316 | **80.9%** | 1,135 | 0.572 |

Previously measured: 98.9% on a 19.2 Mb targeted reference. The loss is therefore a
property of genome-scale variant density, not of pathological input.

**Upstream consequence.** That retention line only exists because this fork added it to
`hgfm.h:2428`. Stock `hisat2-build` drops these variants **silently**. So the shipped
`grch38_snp` index very likely contains materially fewer variants than the `.snp` file it
was built from, and no user has ever been told. Strong candidate for upstream issue #2.

**Project consequence.** The T2 reference-bias result was obtained with roughly four
fifths of the variant set actually present in the index. The effect was real and
replicated; this says it was measured with a handicap, not that it was overstated.

## Density costs memory and retention, not time or disk

| | 1000G_30X | big-4 | Δ |
|---|---|---|---|
| variant instances | 1,731,945 | 1,985,316 | +14.6% |
| retained instances | 1,501,926 | 1,605,978 | **+6.9%** |
| peak RSS | 12.65 GB | 16.81 GB | +33% |
| wall | 2:40.4 | 2:47.7 | +4.5% |
| index size | 579 MB | 594 MB | +2.6% |

**The marginal variant survives at 41%.** Adding 14.6% more variants buys 6.9% more
retained variants, for a third more memory and 64% more blown edge budgets.

## Haplotype realism is roughly free; density is what costs

Invented haplotypes (greedy colouring, `hisat2_extract_snps_haplotypes_VCF.py:298-330`,
which is what a sample-less VCF gets) are a *partition* — 0.78 haplotype rows per variant.
Real phasing gives **1.09**, because a site genuinely appears in several observed
combinations.

I predicted before the build that this would cost extra edge budget. It does not:
per 1,000 variant instances the over-budget rate is **0.392 (real) vs 0.400 (invented)**.
Only variant count drives it — big-4 sits at 0.572.

Caveat, stated rather than buried: the three arms differ in density *and* haplotype
structure, so this is not a controlled comparison. The density-normalised column is the
part that survives the confound.

## Build-parameter titration

Baseline = big-4, `-p 32`: 2:47.7, 16.81 GB, 607,788 KB index.

| config | wall | peak RSS | index | verdict |
|---|---|---|---|---|
| `-p 8` | 3:39.4 (+31%) | 17,604,432 KB (−0.1%) | −8 KB | **4× threads buys 1.31×** |
| `--offrate 5` | 2:48.6 (+0.5%) | +0.004% | −5.5% | weak lever in a graph index |
| `--ftabchars 12` | 3:14.7 (+16%) | ~0 | +10.1% | **misleading at chr1 scale — see below** |
| `--bmaxdivn 8` | 2:48.4 (+0.4%) | +0.002% | identical | **no effect at all** |

Retention was bit-identical (80.9%, 1,135 over-budget) across every config including
`-p 8`, so variant dropping is thread-invariant and parameter-invariant.

Three findings:

1. **`hisat2-build` scales terribly with threads** — 33% efficiency from 8→32. It is
   dominated by serial blockwise suffix-array construction. Paying for vCPUs is waste;
   the only reason to take a 32-vCPU box is that it is the only way to get 256 GB in the
   r7i family.

2. **The documented memory knobs do not control a graph build's peak RSS.**
   `--bmaxdivn 8` moved peak RSS by 384 KB out of 16.8 GB. This corrects a guess I made
   earlier in the run — that the high RSS was auto-tuning being generous because the box
   was large. It is not: the peak is in variant-graph construction, which `--bmax`
   does not govern. Linear chr1 peaks at 1.41 GB, so the entire suffix-array term is
   too small to matter. **There is no flag that trades memory for time here.** The only
   lever on memory is fewer variants.

3. **`--ftabchars 12` inverts at genome scale.** The ftab is `4^k` entries — an absolute
   size, not proportional to the genome. The observed +61.5 MB matches 4^12 vs 4^10
   exactly. On chr1's 594 MB index that is +10%; on an ~8 GB whole-genome index the same
   63 MB is **+0.8%**. Worth taking for the real index, which is the opposite of what the
   raw row says. General warning for reading this table: chr1 ratios only transfer for
   quantities that scale with sequence length.

## What the probe did *not* determine

**Whole-genome peak RSS is still not bounded.** Memory does not scale smoothly with
variant count — 1.20M and 1.73M instances both land at ~12.5 GB, then 1.99M steps to
16.8 GB. With a single chromosome length there is no way to separate the
length-dependent term from the variant-dependent one. JHU's ~160 GB figure is consistent
with these numbers but is not confirmed by them.

Practical consequence: **do not downsize to 128 GB on the strength of the thread-scaling
result.** The two findings point opposite ways and memory wins.

## Costs measured, for stage 2b

| step | chr1 | genome-wide, serial | parallel per chromosome |
|---|---|---|---|
| dbSNP remote tabix slice | 45 s / 2.13 GB | ~20 min | network-bound |
| `FREQ=` filter (gawk) | 13 m 24 s | ~5.4 h | ~20 min on 32 cores |
| dbSNP haplotype extraction | 23 s | ~9 min | trivial |
| 1000G phased extraction | 17 m 48 s | ~7 h | ~40 min on 32 cores |
| graph build | 2 m 48 s | not extrapolable (single global index) | n/a |

Remote tabix at 45 s/chromosome is the notable one: whole-genome variant prep is a
~20 minute network pull, not a 29.5 GB download plus a 1.2-billion-line scan.

## Recommendation for 2b

- **Instance: r7i.8xlarge again.** Not for the cores — for the 256 GB, which is
  unbounded-but-necessary and unavailable at 16 vCPU in this family.
- **Variant set: 1000G 30x phased.** Best retention (92.4%), lowest memory, real
  haplotypes, same build time and index size. dbSNP big-4 is the wrong operating point:
  more than half its marginal variants never reach the index.
- **Flags: `--ftabchars 12`.** Skip `--offrate 5`; skip the memory knobs, which do nothing.
- **Estimate: 5–8 h ≈ $11–17**, dominated by the single global graph build, which is the
  one term chr1 cannot price.
