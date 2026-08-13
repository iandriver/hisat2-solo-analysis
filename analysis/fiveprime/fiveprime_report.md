# 5' chemistry: testing the hypothesis 3' data could not

10x PBMC 1k **5' v2** (79.6M reads, one healthy 27-year-old male donor), graph
`grch38_snp` against linear `grch38`, aligner the only variable. Same gene
model, same counting code as the 3' work.

Three 5'-specific parameters, each confirmed empirically before the run rather
than assumed:

| | 3' v3 | 5' v2 |
|---|---|---|
| barcode list | 3M-february-2018 | **737K-august-2016** (91.0% valid vs 6.5%) |
| UMI length | 12 bp | **10 bp** |
| `--gene-strand` | Forward | **Reverse** (48.2% unique-gene vs 4.1%) |

Alignment: graph 88.33%, linear 88.23%. 944 cells from both, Jaccard 1.0000.
Total UMIs in cells 4,999,060 vs 4,961,669 (**+0.75%**). Median
log2(graph/linear) across 2,672 expressed genes is again **0.0000**.

---

## 1. The chemistry shift did what it was supposed to

Percentage of CDS coverage falling inside the peptide-binding domain (class I
alpha1+alpha2, class II beta1/alpha1), defined by codon on the MANE Select CDS:

| gene | 3' | 5' |
|---|---|---|
| HLA-A | 41.7% | **92.3%** |
| HLA-B | 44.4% | **88.5%** |
| HLA-C | 48.8% | **92.9%** |
| HLA-DRB1 | 13.2% | **75.9%** |
| HLA-DQB1 | 23.9% | **81.9%** |
| HLA-DQA1 | **0.3%** | **79.9%** |

The original hypothesis — that the gain concentrates in the polymorphic exons —
is now testable. On 3' data HLA-DQA1's domain carried 0.3% of coverage; on 5'
it carries 80%.

*(Exon indices could not be used for this. The HLA transcripts in this
annotation do not follow the textbook class I exon structure, and choosing the
longest transcript gives HLA-B a 1270 bp first exon where the real leader exon
is ~73 bp. Defining the window by codon on the MANE Select CDS avoids it.)*

## 2. The result: confirmed for class II, mostly not for class I

Graph/linear **read coverage**, peptide-binding domain versus the rest of the
same CDS. The domain-vs-rest contrast is the clean comparison — same donor,
same gene, same library:

| gene | | domain | rest | enrichment |
|---|---|---|---|---|
| **HLA-DRB1** | 5' | **16.49** | 1.92 | **8.6x** |
| | 3' | 7.02 | 1.18 | 5.9x |
| **HLA-DQA1** | 5' | **14.44** | 2.67 | **5.4x** |
| | 3' | 3.44 | 41.27 | 0.08x |
| **HLA-DQB1** | 5' | **6.15** | 4.68 | 1.3x |
| | 3' | 3.57 | 2.64 | 1.4x |
| **HLA-B** | 5' | **2.28** | 1.05 | **2.2x** |
| | 3' | 2.01 | 1.44 | 1.4x |
| HLA-A | 5' | 1.17 | 1.43 | 0.8x |
| | 3' | 2.25 | 1.62 | 1.4x |
| HLA-C | 5' | 1.02 | 0.99 | 1.0x |
| | 3' | 1.01 | 0.87 | 1.2x |
| **HLA-DRA** *(invariant control)* | 5' | **1.00** | 0.99 | — |
| | 3' | 1.01 | 1.02 | — |
| **B2M** *(control)* | 5' | — | **1.00** | — |

- **Class II is confirmed and sharpened.** DRB1 gains 16.5x inside the
  peptide-binding groove against 1.9x outside it. DQA1 14.4x against 2.7x. The
  gain is where the polymorphism is.
- **HLA-DRA is the control that makes it mean something** — the invariant alpha
  chain, same locus, same cells, same library: **1.00 in its own domain**. The
  effect tracks polymorphism, not the MHC region or expression level.
- **HLA-B does show domain-localised gain on 5'**: 2.28x in the groove, 1.05x
  outside. That is a real class I signal and it is only visible on 5' data.
- **HLA-A and HLA-C do not.** HLA-C is 1.02 in the domain on both chemistries.

**So the prediction was half right.** Class I is not uniformly rescued by 5'
chemistry; one of three genes shows a clean domain-localised effect in this
donor. The claim should stay centred on class II, with HLA-B as a demonstration
that class I *can* be affected when the donor's haplotype is divergent.

## 3. The thing worth knowing: reads are lost far more than molecules

Read-level and UMI-level gains for the same genes, 5' data:

| gene | reads (MAPQ 60) | UMIs |
|---|---|---|
| HLA-DQA1 | **7.19x** | 2.89x |
| HLA-DRB1 | **4.48x** | 2.45x |
| HLA-B | **1.62x** | 1.07x |
| HLA-A | 1.17x | 1.11x |
| HLA-DRA *(control)* | 1.00x | 1.00x |

Sequencing saturation is **86.3%** here (68.8% on the 3' library). At that
depth most of the reads the graph recovers are PCR duplicates of molecules the
linear index already captured, so a large read-level recovery collapses into a
small UMI-level one.

**This reconciles the whole picture, including the 3' result.** Reference bias
costs *reads* everywhere it acts, but it only costs *molecules* where the loss
is close to total. HLA-DQA1 on 3' data lost 87-96% of its reads — enough to
lose the molecules — and the UMI gain was 32x. HLA-B on 5' loses about a third
of its domain reads, the molecules survive on what remains, and the UMI gain is
1.07x.

Practical consequence: **the UMI-level benefit is largest at low sequencing
saturation and for severely divergent genes.** A shallow library, or a gene the
donor's haplotype puts far from the reference, is where this matters. A deeply
saturated library of a near-reference gene will show almost nothing.

## 4. The ribosomal-protein defect reproduces on 5'

Top gainers are again processed pseudogenes: RPS28P7 +8.61, RPS7P14 +8.08,
RPL39P3 +7.38, RPS3AP6 (12,452 vs 247). Graph multi-mapping is 18.71% against
linear's 15.98%, a wider gap than 3' (8.47% vs 6.46%).

Consistent with the mechanism traced earlier — junction-spanning reads losing to
contiguous retrocopies — and chemistry-independent, as that explanation predicts.

## 5. Caveats

- **Different donor.** The 3' and 5' libraries are different individuals, so
  cross-chemistry comparison of any single gene's ratio confounds chemistry with
  haplotype. HLA-C at 1.51x on 3' and 0.88x on 5' is a donor difference, not a
  chemistry effect. **The domain-vs-rest contrast within one library is the
  comparison that is not confounded**, which is why it is the one reported.
- One donor per chemistry. HLA-A and HLA-C being flat here says this donor is
  near-reference at those genes, not that class I is never affected.
- Coverage is read-level at MAPQ 60; UMI counts are from the solo matrices.
  The two answer different questions, as section 3 shows.
