# T1 + T2: does variant-aware alignment change single-cell conclusions?

Two donors, two chemistries. Sections 1-4 use 10x PBMC 1k **3' v3**, 66.6M
reads; section 5 adds a second donor on **5' v2**, 79.6M reads. In both, the
same gene model, the same counting code, the same machine, and **the aligner
index as the only variable** — prebuilt `grch38_snp` (14.46M SNPs) against
prebuilt linear `grch38`.

| | 3' v3 | 5' v2 |
|---|---|---|
| alignment, graph vs linear | 87.20% / 87.07% | 88.33% / 88.23% |
| UMIs in cells | +0.39% | +0.75% |
| median log2(graph/linear), expressed genes | **0.0000** (n=4,042) | **0.0000** (n=2,672) |
| sequencing saturation | 68.8% | 86.3% |

The graph does not add signal everywhere. In both libraries the median gene is
unchanged, so the interesting question is entirely about *where* the differences
land — and, it turns out, about how deeply the library was sequenced.

---

# T1 — HLA expression under linear alignment

## Result: confirmed for HLA-DQ and HLA-C, not for HLA-A/B/DP

| gene | graph | linear | ratio |
|---|---|---|---|
| HLA-DQA1 | 2,565 | 80 | **32.1x** |
| HLA-DRB5 | 77 | 4 | 19.3x |
| HLA-DQB1 | 2,459 | 820 | **3.00x** |
| HLA-DRB1 | 6,435 | 4,171 | 1.54x |
| HLA-C | 16,913 | 11,198 | **1.51x** |
| HLA-DRA | 18,619 | 18,091 | 1.03x |
| HLA-A | 15,419 | 15,216 | 1.01x |
| HLA-B | 23,473 | 23,626 | 0.99x |
| HLA-DPA1 / DPB1 | 10,162 / 7,513 | 10,162 / 7,517 | 1.00x |
| B2M (control) | 73,097 | 73,097 | 1.00x |

Which genes are affected is a property of **this donor's haplotypes**, not of
HLA in general: DP is untouched, DQ is transformed.

## It is recovery, not reallocation

Two checks, because a gene can gain UMIs either by recovering reads that were
lost, or by stealing them from a paralog. Only the first is an advantage.

**Family sums** rule out reallocation within the locus:

| family | graph | linear | delta |
|---|---|---|---|
| HLA-DQA (DQA1+DQA2) | 2,982 | 475 | +2,507 (6.28x) |
| HLA-DQB (DQB1+DQB2+DQB3) | 2,463 | 827 | +1,636 (2.98x) |
| HLA-DRB (DRB1,5,6,9) | 7,996 | 4,446 | +3,550 (1.80x) |
| HLA class I (A,B,C,E,F,G) | 64,287 | 58,514 | +5,773 (1.10x) |
| HLA-DP (DPA1-3, DPB1-2) | 17,676 | 17,680 | -4 (1.00x) |

**Read tracing** shows what the linear index actually did with them. Taking the
7,048 reads the graph places uniquely (MAPQ 60) in HLA-DQA1's divergent 3'
window (chr6:32,642,179-32,644,179) and looking each one up in the linear
alignment:

| where linear put them | n | % |
|---|---|---|
| **unaligned** | 6,777 | **96.2%** |
| aligned, MAPQ 60 | 252 | 3.6% (113 of them at the paralog HLA-DQA2) |
| aligned, ambiguous | 19 | 0.3% |

Under the graph, 4,401 of the same 7,048 reads (62.4%) align with **zero
mismatches**. Linear GRCh38 does not misplace these reads — it cannot align
them at all. The donor's DQA1 haplotype is simply too far from the reference.

Repeating that for four genes, over the whole gene body rather than the
divergent window — reads the graph places uniquely, looked up in the linear
alignment:

| gene | graph-unique reads | **unaligned in linear** | unique in linear | ambiguous |
|---|---|---|---|---|
| HLA-DQA1 | 11,406 | **9,959 (87.3%)** | 1,414 | 33 |
| HLA-DQB1 | 8,597 | **7,177 (83.5%)** | 1,417 | 3 |
| HLA-DRB1 | 20,404 | **11,541 (56.6%)** | 8,853 | 10 |
| HLA-C | 66,993 | **32,124 (48.0%)** | 34,377 | 492 |

In every case the "unique in linear" column is within ~1.5% of linear's *own*
MAPQ-60 read count at that gene (1,258 / 1,429 / 8,969 / 34,932). That is the
clean statement: **the graph keeps everything linear already had and adds reads
linear could not align at all.** Nothing is taken from anywhere else.

The size of the gain tracks how far the donor's haplotype sits from GRCh38 —
87% new reads at DQA1, 48% at HLA-C.

## Where in the gene — untestable on 3', since answered on 5'

The plan predicted the gain would concentrate in the polymorphic exons (2 and 3
for class I, exon 2 for class II beta). **That cannot be tested on 3' 10x data**:
essentially all coverage is in the 3'-terminal exon, so exons 2-3 carry too few
reads either way. Coverage of HLA-DQA1's exon 2 is 916 (graph) / 276 (linear);
its terminal exon carries 603,967 / 12,430.

**A 5' library has since been run and does test it — see section 5. Short
version: confirmed for class II, one class I gene out of three.**

What the data does show is arguably more relevant to single-cell work — **the
divergence is in the 3' UTR, which is exactly the window 3' chemistry
sequences**:

| gene | terminal-exon coverage, graph | linear | ratio |
|---|---|---|---|
| HLA-DQA1 | 603,967 | 12,430 | **48.6x** |
| HLA-DQB1 | 634,431 | 102,539 | 6.19x |
| HLA-DRB1 | 1,518,837 | 666,347 | 2.28x |
| HLA-C | 5,439,694 | 2,496,054 | 2.18x |
| HLA-B | 7,705,194 | 7,163,682 | 1.08x |
| **HLA-DRA (invariant alpha chain)** | 6,695,196 | 5,813,744 | **1.15x** |

HLA-DRA is the internal control: it is the non-polymorphic partner of DRB1, in
the same locus, in the same cells. It moves 1.15x while its polymorphic partner
moves 2.28x and DQB1 moves 6.19x.

## Does it change a conclusion?

Partly. Per cell type, expression per 10k UMIs:

| | class II total, graph/linear | class I total | **HLA-DQ only** |
|---|---|---|---|
| B | 132.4 / 109.4 = 1.21x | 1.11x | 19.2 / 3.4 = **5.65x** |
| CD14 Mono | 70.7 / 61.5 = 1.15x | 1.12x | 5.8 / 1.1 = **5.39x** |
| CD16 Mono | 84.9 / 75.0 = 1.13x | 1.13x | **6.60x** |
| T | 1.2 / 1.0 = 1.16x | 1.10x | **3.59x** |
| NK | 3.7 / 3.1 = 1.19x | 1.12x | **8.89x** |

- **Aggregate class I and class II gains are near-uniform across cell types**, so
  ratios between cell types barely move: monocyte:T class II ratio 60.06 vs
  60.53 (**-0.8%**), class I 0.93 vs 0.91 (+2.0%). A paper comparing class II
  between cell types would reach the same answer either way.
- **HLA-DQ specifically does move the comparison**: monocyte:T DQ ratio 56.19
  (graph) vs 37.41 (linear), **+50.2%**.
- Clustering derived from each matrix: 14 clusters both, **ARI 0.873, NMI 0.910**
  — same cell types, some boundary movement (see the caveat below for the likely
  cause).

**Verdict: partial success.** The absolute quantification of specific
polymorphic genes is badly wrong under linear alignment — HLA-DQA1 is
effectively invisible (80 UMIs) when it should be comparable to DQB1, and the
mechanism is proven to be unalignable reads rather than misplaced ones. But the
gain is broad enough across cell types that most *comparative* conclusions
survive; only DQ-specific comparisons shift materially.

---

# T2 — does the bias produce false allele-specific-expression calls?

## Result: confirmed, and one-sided

Read-level REF/ALT counts at dbSNP SNV positions from both alignments (MAPQ 60,
BQ >= 20), 1,020,853 positions covered in both.

### Dose-response control first

Sites classified by *pooled* genotype, so neither aligner picks its own sites:

| genotype | n | ALT frac graph | ALT frac linear | depth ratio graph/linear |
|---|---|---|---|---|
| hom-ref | 71,529 | 0.0000 | 0.0000 | **1.0000** |
| het | 12,923 | 0.5000 | 0.4500 | **1.0603** |
| hom-alt | 7,459 | 1.0000 | 1.0000 | **1.1290** |

The graph's depth advantage scales with alternate-allele dosage: none at 0
copies, +6.0% at 1, +12.9% at 2. That is the signature of reference bias and
nothing else. At homozygous-reference sites the graph's mean ALT fraction is
0.00042 against linear's 0.00008 — it does not manufacture alternate reads
where the donor carries none.

### The ASE test

12,580 heterozygous sites, **selected by the linear aligner** (the conservative
choice — only sites linear itself sees both alleles at). Binomial test vs
p=0.5 per site, BH-corrected at FDR 0.05, **depth-matched** by subsampling both
aligners to the same per-site depth so power is identical (mean of 5 replicates):

|  | graph balanced | graph imbalanced |
|---|---|---|
| **linear balanced** | 9,865 | 348 |
| **linear imbalanced** | **986** | 1,379 |

| | n | ref-skewed | alt-skewed |
|---|---|---|---|
| **linear-only calls** | 986 | **950 (96.3%)** | 36 |
| **graph-only calls** | 340 | 23 | **317 (93.2%)** |

- **986 sites (7.8%) get a false allelic-imbalance call from linear alignment**,
  and 96.3% of them are skewed toward the reference allele — the exact
  signature of alignment-stage bias rather than biology.
- **340 sites (2.7%) have real alternate-biased imbalance masked** by linear,
  93.2% of them alt-skewed.
- Together, **~10.5% of heterozygous sites get the wrong ASE answer**.
- At the linear-only sites, median ALT fraction is 0.268 under linear and 0.388
  under the graph. Linear does not just cross a threshold; it reports a
  substantially different allelic ratio.

Without depth matching the effect is larger still (linear-only 842 vs graph-only
509 at raw depth on the same sites, and 1,191 vs 436 on graph-selected sites),
so matching is the conservative presentation.

**Verdict: success.** This is the claim in the form a reviewer wants: not "8-10%
bias" but "one call in ten is wrong, and the errors point one way."

---

# The countervailing finding: ribosomal-protein pseudogenes

This is the most important thing found today and it is **not** in the plan's
favour. It should be fixed before the feature is promoted.

The +0.39% global UMI change is not a small uniform gain. It is:

| biotype | genes | net delta | graph | linear |
|---|---|---|---|---|
| processed_pseudogene | 10,148 | **+159,541** | 220,590 | 61,049 |
| transcribed_processed_pseudogene | 512 | +24,364 | 33,365 | 9,001 |
| protein_coding | 20,070 | **-154,943** | 7,051,598 | 7,206,541 |

Almost all of it is one family:

| | graph | linear | delta |
|---|---|---|---|
| ribosomal protein genes (n=83) | 1,422,115 | 1,525,147 | **-103,032** |
| their pseudogenes (n=1,498) | 177,632 | 42,959 | **+134,673** |

Worst individual cases: RPSA 8,411 -> 341 UMIs (-96%), RPS10 3,500 -> 219,
RPS27 18,465 -> 3,825, RPL10 29,419 -> 13,830.

## Mechanism, traced

Taking the 125,480 reads that linear places uniquely (MAPQ 60) at RPS27 and
looking them up in the graph alignment:

| under the graph | n | % |
|---|---|---|
| MAPQ 60 (still unique) | 19,856 | 15.8% |
| **MAPQ 1 (multi-mapping)** | 104,216 | **83.1%** |
| MAPQ 0 | 1,408 | 1.1% |

Mismatch counts are **unchanged** (111,188 reads at NM=0 under the graph vs
111,118 under linear). The reads did not find a better home — the graph created
*equally good* alternative placements at processed-pseudogene copies
(chr3:40.76M, chr1:202.47M, chr12:3.21M, chr11:117.03M ...), destroying
uniqueness.

### CORRECTION — the cause proposed here was wrong

This section originally hypothesised that dbSNP entries inside processed
pseudogenes are paralogous sequence variants, and that removing them would fix
the problem. **A rebuild has since falsified that** (see `rebuild_report.md`).
The paragraph is kept because the measurements above are still correct; only
the explanation was wrong.

The filter was written, and a three-way index comparison built on a 19.2 Mb
targeted reference reproduces the defect faithfully (RPS27 0.128 of its linear
unique reads, against 0.16 on the real index) with controls at 1.000. Removing
every pseudogene variant recovers only **9.7%** of the deficit.

Why: of the competing alignments under the graph, **88% are exact matches that
use no ALT path at all** (no `Zs:Z` tag). They cannot be caused by variants.

Ground truth, counting exact occurrences of read sequences in the reference
independently of any aligner — 2,000 reads the linear index calls uniquely
mapped at RPS27:

| exact copies in the reference | reads |
|---|---|
| 1 (genuinely unique) | 145 (7.2%) |
| **2 or more** | **1,574 (78.7%)** |
| 0 (read carries a sequencing error) | 281 (14.1%) |

Verified against the whole GRCh38 assembly for a subset: 18 of 20 reads occur
exactly 3 times genome-wide. On the real prebuilt indexes those same 20 reads
come back **20/20 at MAPQ 60 with `NH:i:1` from linear `grch38`**, and 19/20 at
MAPQ 1 with `NH:i:4` from `grch38_snp`. `NH:i:1` is the aligner reporting one
alignment for a read that matches three places exactly.

**So the linear index is assigning MAPQ 60 to reads with several perfect genomic
matches, and the graph index is right to call them ambiguous.** The linear
numbers that looked correct — RPSA 8,411 UMIs, RPS27 18,465 — are inflated by
falsely-unique assignments.

This is still a practical problem for quantification, because scattering a real
gene's reads across identical pseudogene copies collapses its count under
`--solo-multi-mappers Unique`. But the remedy is a multi-mapper policy (EM or
Uniform), not a variant filter, and the underlying issue belongs to the
reference and the annotation rather than to variant-aware alignment. HISAT2's
linear index reporting MAPQ 60 for reads with 2-4 exact matches is worth an
upstream report on its own.

**Superseded.** The proposed filter (`hisat2_filter_snps.py`) was built and
tested and drops 78,935 of 14,460,407 variants (0.55%) on the full GRCh38
variant set — *not* the 181,160 / 1.40% quoted in an earlier draft, which came
from a regex that also matched **un**processed pseudogenes. It does not fix
this problem; see the correction above.

This also plausibly explains part of the ARI 0.873 clustering movement, since
ribosomal protein genes are a large fraction of every cell's UMIs.

---

# 5. Second donor, 5' chemistry — the exon hypothesis, tested

10x PBMC 1k **5' v2**, 79.6M reads, different donor, same indexes and counting
code. Three parameters differ from 3' v3 and each was confirmed empirically:
the 737K barcode list (91.0% valid vs 6.5% for 3M), a **10 bp** UMI, and
`--gene-strand Reverse` (48.2% unique-gene vs 4.1% for Forward).

Alignment 88.33% vs 88.23%; 944 cells from both, Jaccard 1.0000; UMIs in cells
+0.75%; median log2(graph/linear) across 2,672 expressed genes again **0.0000**.

## The chemistry moves coverage onto the polymorphism

Share of CDS coverage inside the peptide-binding domain (class I alpha1+alpha2,
class II beta1/alpha1), defined **by codon on the MANE Select CDS**:

| gene | 3' | 5' |
|---|---|---|
| HLA-DQA1 | **0.3%** | **79.9%** |
| HLA-DRB1 | 13.2% | 75.9% |
| HLA-DQB1 | 23.9% | 81.9% |
| HLA-A | 41.7% | 92.3% |
| HLA-B | 44.4% | 88.5% |
| HLA-C | 48.8% | 92.9% |

*(Exon indices could not be used. These HLA transcripts do not follow the
textbook class I exon structure, and picking the longest transcript gives HLA-B
a 1270 bp first exon where the real leader exon is ~73 bp.)*

## Result: class II confirmed, class I one gene in three

Graph/linear read coverage, domain versus the rest of the same CDS. Same donor,
same gene, same library — this contrast is not confounded by anything:

| gene | | domain | rest | enrichment |
|---|---|---|---|---|
| **HLA-DRB1** | 5' | **16.49** | 1.92 | **8.6x** |
| **HLA-DQA1** | 5' | **14.44** | 2.67 | **5.4x** |
| **HLA-B** | 5' | **2.28** | 1.05 | **2.2x** |
| HLA-DQB1 | 5' | 6.15 | 4.68 | 1.3x |
| HLA-A | 5' | 1.17 | 1.43 | 0.8x |
| HLA-C | 5' | 1.02 | 0.99 | 1.0x |
| **HLA-DRA** *(invariant control)* | 5' | **1.00** | 0.99 | — |
| **B2M** *(control)* | 5' | — | **1.00** | — |

HLA-DRA is the control that gives this meaning: the invariant alpha chain, same
locus, same cells, same library, **1.00 in its own domain**. The effect tracks
polymorphism, not the MHC region and not expression level.

HLA-B shows a genuine domain-localised class I effect that is invisible on 3'
data. HLA-A and HLA-C show nothing — this donor is near-reference there.

## Reads are lost far more than molecules

| gene | reads (MAPQ 60) | UMIs |
|---|---|---|
| HLA-DQA1 | **7.19x** | 2.89x |
| HLA-DRB1 | **4.48x** | 2.45x |
| HLA-B | **1.62x** | 1.07x |
| HLA-A | 1.17x | 1.11x |
| HLA-DRA *(control)* | 1.00x | 1.00x |

Sequencing saturation is **86.3%** here against 68.8% on the 3' library, so most
recovered reads are duplicates of molecules the linear index already captured.

**This reconciles the 3' result too.** Reference bias costs *reads* wherever it
acts, but costs *molecules* only where the loss is near-total. HLA-DQA1 on 3'
lost 87-96% of its reads — the molecules went with them, hence 32x. HLA-B on 5'
loses about a third of its domain reads, the molecules survive on the rest,
hence 1.07x.

**Scoping consequence: the UMI-level benefit is largest at low sequencing
saturation and for severely divergent genes.** A deeply saturated library of a
near-reference gene will show almost nothing.

## T2 replicates on the second donor

The whole T2 pipeline re-run on this library, 45,234 sites with >=20 reads in
both alignments.

**Dose-response control**, sites classified by pooled genotype:

| alt copies | n | AF graph | AF linear | depth ratio |
|---|---|---|---|---|
| 0 | 34,740 | 0.0000 | 0.0000 | **1.0000** |
| 1 (het) | 5,955 | **0.5000** | 0.4615 | **1.0458** |
| 2 | 3,532 | 1.0000 | 1.0000 | **1.1000** |

The graph sits exactly on 0.5000 again, and the depth advantage again scales
with dosage: nothing at 0 copies, +4.6% at 1, +10.0% at 2. Mean ALT fraction at
hom-ref sites 0.00024 (graph) vs 0.00009 (linear) — it still does not
manufacture alternate reads.

**The ASE test**, 5,781 linear-selected het sites, depth-matched, mean of 5
replicates:

| | 3' donor | 5' donor |
|---|---|---|
| het sites tested | 12,580 | 5,781 |
| **linear-only calls (false imbalance)** | 986 (**7.8%**) | 319 (**5.5%**) |
| ... skewed to the reference allele | **96.3%** | **97.2%** |
| graph-only calls (imbalance linear missed) | 340 (2.7%) | 258 (4.5%) |
| ... skewed to the alternate allele | 93.2% | 79.3% |
| **total discordant** | **10.5%** | **10.0%** |

**The headline result holds on an independent donor and chemistry**: about one
heterozygous site in ten gets a different ASE answer, and the linear aligner's
false calls are ~97% reference-skewed in both.

## The ribosomal-protein defect reproduces

Top gainers are again processed pseudogenes (RPS28P7 +8.61, RPS7P14 +8.08,
RPS3AP6 12,452 vs 247), and graph multi-mapping is 18.71% against linear's
15.98% — a wider gap than 3' (8.47% vs 6.46%). Chemistry-independent, as the
splice-versus-retrocopy mechanism predicts.

---

# Caveats

- **Two donors, two chemistries — and they are confounded with each other.**
  The 3' and 5' libraries are different individuals, so no single gene's ratio
  can be compared across them: HLA-C at 1.51x on 3' and 0.88x on 5' is a donor
  difference, not a chemistry effect. Only the within-library domain-vs-rest
  contrast is clean, which is why that is what section 5 reports.
- **Which HLA genes are affected depends on the donor's haplotypes.** DP was
  untouched in the 3' donor; HLA-A and HLA-C were untouched in the 5' donor.
  Neither says the gene is never affected.
- ~~**Het sites were called from the data, not from an independent genotype.**~~
  **Closed.** A third donor with GIAB HG001 benchmark genotypes (GM12878, see
  `lcl/t4b_report.md`) gives median ALT fraction **0.4762 linear against an
  exact 0.5000 for the graph** at externally confirmed heterozygous sites —
  matching the 0.4500 and 0.4615 measured here from data-derived sites.
- **The graph's spurious ALT signal is not negligible after all.** Against GIAB
  truth it produces **3x the false heterozygous calls** at homozygous-reference
  sites (0.247% vs 0.081% at depth >=20). Averaged over all sites it looks tiny,
  but it concentrates at a few positions and crosses genotype-calling
  thresholds. The bullet below understated this.
- **The exon-level hypothesis is untestable on 3' data.** It was tested on 5'
  instead (section 5); the 3' sections stand as written.
- **The graph accepts a small amount of spurious ALT signal** — mean ALT
  fraction 0.00042 at hom-ref sites vs linear's 0.00008. Negligible for ASE at
  these depths, but real.
- **MAPQ 60 filtering** was applied throughout the read-level work. It is the
  right filter for ASE and is applied identically to both, but it interacts with
  the pseudogene finding above, where the graph's reads lose MAPQ.

# What this supports claiming, and what it does not

**Claim:** variant-aware alignment recovers polymorphic-gene expression that
linear alignment cannot align at all (HLA-DQA1 32x, 96% of the recovered reads
previously unalignable), the recovery is localised to the peptide-binding domain
where the polymorphism is (HLA-DRB1 16.5x inside against 1.9x outside, with the
invariant HLA-DRA at 1.00 as a control), and it removes a reference-allele bias
that produces a false allele-specific-expression call at roughly one
heterozygous site in thirteen.

**Scope it honestly:** the effect is on class II in both donors and on one class
I gene in one donor. The UMI-level benefit is much smaller than the read-level
one and shrinks as sequencing saturation rises, so it is largest for shallow
libraries and severely divergent genes.

**Do not claim:** that it improves single-cell analysis generally. 99% of genes
are unchanged in both donors, cell-type comparisons of aggregate HLA barely
move, clustering agreement is 0.873 rather than 1.0, and the current prebuilt
SNP index degrades ribosomal-protein quantification on both chemistries.
