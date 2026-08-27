# XCI — X-inactivation escape under linear vs graph alignment

**Result: a graph index that carries chrX variants recovers escape genes a linear
reference misses — ZFX in the European donor, EIF1AX in the Yoruba donor, with
none lost either way. The per-site effect is +0.013 (CEU) and +0.006 (YRI) on the
alternate-allele fraction, and the mechanism control is exact: an index with no
non-PAR chrX variants moves nothing (z = -1.3 to +0.9), while the same index with
them moves everything (z = +5.9 to +8.6).**

The measurement is thin — 208-349 usable sites, 1-2 genes gained per donor — and
the chrX genotype truth is demonstrably unreliable. This is a pilot with a
working mechanism, not a corrected escape catalogue.

## Why chrX, and why now

Escape is detected as biallelic expression at heterozygous sites: an
allele-fraction measurement, which is the quantity P1 showed a linear reference
distorts by 0.031 (CEU) to 0.048 (YRI) per site. The failure mode is directional
and does not average out. If the **inactive** X carries the ALT allele, reference
bias suppresses the very reads that prove escape, so the gene reads as silenced.
Which X is inactivated is random per individual, so this injects
individual-specific false negatives rather than a constant offset that
calibration could remove. The published escape catalogues are built on linear
alignment.

This became runnable only after repairing chrX in our panel: the whole-genome
index carried 184 variants/Mb on chrX against an autosome mean of 5,158, all of
them in PAR1, because `hisat2_extract_snps_haplotypes_VCF.py` cannot read the
haploid male genotypes that non-PAR chrX carries in a phased panel. See
`extractor_haploid_fix.patch` and `../../design/rust/README.md`.

## Design

Four arms, same 50,000,000 cDNA reads per donor, same sites, fully paired. Only
sites covered in all four arms count, so a difference between arms cannot come
from which sites happened to be covered. `-q 60` keeps uniquely-mapped reads
only: a multi-mapper at a het site would let mapping ambiguity masquerade as
allelic imbalance, which is the thing being measured.

| arm | index | chrX variants |
|---|---|---|
| `linear` | `grch38` | 0 |
| `snp_invented` | `grch38_snp` (JHU) | 465,201 |
| `snp_real_Xbroken` | `wg64_idx` | 28,711 — **PAR1 only** |
| `snp_real_Xfixed` | rebuilt with the extractor fix | 558,026 |

The last two are the same index family differing **only** in chrX content. That
is the control: arm 4 should move non-PAR chrX and arm 3 should not.

Both donors are female — verified from diploid `0|0` calls on non-PAR chrX in the
panel — so both have an inactive X and escape is measurable in principle. chrX
het sites come from the 1000 Genomes phased panel; GIAB HG001 stops at chr22, so
there is no gold-standard chrX truth for either donor. That limitation drives the
caveats below.

**PAR is excluded everywhere.** PAR genes escape X-inactivation by definition, so
they are a positive control rather than a discovery set — and PAR1 is precisely
where the "broken" index does carry variants. Including PAR sites makes the
control arm appear to gain +0.074 and destroys the experiment; see The
correction.

## Result — per-site allele fraction, non-PAR chrX

| donor | arm | depth>=10 (n=208 / 349) | depth>=20 (n=86 / 176) |
|---|---|---|---|
| CEU | `snp_invented` vs linear | +0.01308, z=+8.7 | +0.01321, z=+6.4 |
| | `snp_real_Xbroken` vs linear | **-0.00104, z=-1.3** | **+0.00039, z=+0.2** |
| | `snp_real_Xfixed` vs linear | **+0.01290, z=+8.6** | **+0.01220, z=+6.3** |
| YRI | `snp_invented` vs linear | +0.00620, z=+6.9 | +0.00716, z=+5.6 |
| | `snp_real_Xbroken` vs linear | **+0.00013, z=+0.8** | **+0.00040, z=+0.9** |
| | `snp_real_Xfixed` vs linear | **+0.00600, z=+7.1** | **+0.00734, z=+5.9** |

**The control is exact.** The index that carries no non-PAR chrX variants moves
the allele fraction by 0.0001 to 0.001 at |z| <= 1.3 — nothing. The same index
with those variants added moves it by 0.006 to 0.013 at z = +5.9 to +8.6. This is
P1's mechanism reproduced on a different chromosome with a within-family control
P1 did not have: an index helps at a site if and only if it carries that site's
variant.

Isolating the chrX contribution directly, `X-fixed` minus `X-broken`:

| donor | depth>=10 | depth>=20 |
|---|---|---|
| CEU | +0.01369, 87 up / 1 down, z=+9.2 | +0.01154, 41 up / 0 down, z=+6.4 |
| YRI | +0.00519, 61 up / 5 down, z=+6.9 | +0.00755, 45 up / 4 down, z=+5.9 |

Sites move up and essentially never down — 1 and 5 down against 87 and 61 up —
which is the signature of recovering alt-carrying reads a variant-less reference
was losing, not of noise.

Real phasing and invented haplotypes are indistinguishable here (+0.0129 vs
+0.0131 for CEU; +0.0060 vs +0.0062 for YRI), matching P1's autosomal finding for
the European donor. The value of our index on chrX is that it *has* chrX at all,
not that its haplotypes are real.

## X-inactivation skewing: neither line is clonally skewed

| donor | mean hap2 fraction | monoallelic sites | split |
|---|---|---|---|
| GM12878 | 0.4732 | 57% | 33% hap1 / 24% hap2 |
| GM18502 | 0.4741 | 75% | 39% hap1 / 35% hap2 |

No consistent haplotype preference, so neither LCL is clonal. But most sites read
monoallelic, **splitting both ways**, which fits no X-inactivation model: clonal
skewing would favour one haplotype throughout, and random inactivation would
average to 0.5 in bulk.

The monoallelic fraction is **depth-independent** — 26% to 35% (CEU) and 72% to
76% (YRI) as depth rises from 20 to 100 — so sampling noise is excluded. The
remaining explanation is **genotype error**: sites called heterozygous in the
panel that are not heterozygous in these lines. P1 flagged this risk for the
Yoruba donor; on chrX it applies to both, and three times harder for YRI.

**Site-level chrX conclusions are therefore not trustworthy.** Gene-level
aggregation over multiple sites is, because a single bad site cannot dominate.

## Gene-level escape calls

A gene is called escape-like with >= 2 usable sites, >= 20 reads, and a
haplotype fraction in 0.30-0.70.

| donor | linear | invented | X-broken | X-fixed |
|---|---|---|---|---|
| GM12878 | 5 | 6 | 5 | **6** |
| GM18502 | 8 | 8 | 7 | **9** |

```
GM12878  X-fixed vs linear   gained ZFX            lost none
GM18502  X-fixed vs linear   gained EIF1AX         lost none
         X-fixed vs X-broken gained CHM, EIF1AX    lost none
```

The gene sets validate themselves. **KDM6A, ZFX, EIF2S3, EIF1AX, SMC1A, PHF8,
PRKX and TRAPPC2 are all documented escape genes**, recovered without being told
about them, and ZFX and PRKX appear in both donors. The two genes that only the
repaired index finds — ZFX and EIF1AX — are established escapees, not novel
calls, which is the right kind of result for a method check: it recovers known
biology that linear alignment missed, in the direction the bias mechanism
predicts, with nothing lost.

## The correction

The first pass at this analysis reported that the control **failed** — that
`snp_real_Xbroken` gained +0.074 (CEU) and +0.058 (YRI) on chrX despite carrying
no non-PAR chrX variants, and that the design therefore did not isolate what it
claimed.

That was an analysis error, not a result. PAR sites were included in the
comparison, and PAR1 is exactly where `wg64_idx` does carry its 28,711 variants.
The apparent gain was real, legitimate, and entirely from PAR. Excluding PAR —
which the design always required, since PAR escapes X-inactivation by definition
— the control arm moves nothing, as the mechanism predicts.

Recorded because the failure mode is instructive: a control region chosen for one
reason (PAR escapes XCI, so it is not a discovery set) turned out to matter for a
second, unrelated reason (PAR is where the broken index has its variants), and
missing that inverted the headline.

## Caveats

- **No gold-standard chrX truth.** GIAB HG001 is chr1-22. Both donors' chrX het
  calls come from the 1000 Genomes panel, and the skewing analysis shows they
  carry substantial error — 26-35% (CEU) and 72-76% (YRI) of "het" sites are
  monoallelic at high depth. Gene-level aggregation tolerates this; site-level
  work does not.
- **Coverage is the binding constraint.** 84 and 149 chrX genes have any usable
  site out of 2,484 annotated; 33 and 63 have two or more. One and two genes
  gained is not a rate.
- **No ancestry claim.** The per-site effect is larger for the European donor
  (+0.013) than the Yoruba one (+0.006), the opposite of the autosomal
  expectation, but with this n and this much genotype error the comparison
  carries no weight.
- **Escape calls here are not escape.** An allele fraction in 0.30-0.70 is
  consistent with escape, with mixed inactivation, and with a mis-called
  genotype. Distinguishing them needs reliable phased truth and deeper coverage.
- **Both donors are in the 1000 Genomes panel our index was built from**, so the
  index's coverage of their variants is not a held-out measurement. The
  `X-broken` arm is the uncircular control and behaves correctly.

## Files

```
run_escape.sh          the four-arm run, extending analysis/p1/run_p1.sh
mk_xtruth.sh           chrX het sites for both donors from the phased panel
analyse_xci.py         every number in this report, from counts/ and the phased TSVs
results.txt            its output
counts/                per-arm base counts at chrX het sites
*.chrX.phased.tsv.gz   phased chrX heterozygous genotypes per donor
extractor_haploid_fix.patch  the one-line fix that made a populated chrX possible
```
