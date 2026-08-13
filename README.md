# hisat2-solo — analysis

Measurement scripts, reports and result tables behind the HISAT2 work in
[iandriver/hisat2](https://github.com/iandriver/hisat2). Kept out of that repo so
the source tree stays clean for a possible pull request.

The code being measured lives on two branches there:

| branch | what |
|---|---|
| `upstream/modernize` | fixes that apply to HISAT2 generally — Python 3, C++17, `std::thread`, zlib, BAM output, build portability, variant retention |
| `solo/gene-model` | the above plus single-cell: gene model, barcode/UMI, counting, allele-specific counting |

## Findings

### Variant-aware alignment reduces reference bias (`analysis/human`, `analysis/fiveprime`)

At heterozygous sites an unbiased aligner reports ALT at 0.5. Measured on three
donors, GIAB-verified:

| | linear | graph |
|---|---|---|
| ALT fraction at het sites | 0.4500 / 0.4615 / 0.4762 | **0.5000** |

About one het site in ten gets a different allele-specific-expression answer,
and linear's false calls are ~97% reference-skewed. Replicated independently on
5′ chemistry, where the enrichment concentrates in the peptide-binding domain
(HLA-DRB1 16.49 vs 1.92) while the invariant control HLA-DRA sits at 1.00.

HLA recovery is real recovery, not reallocation: 87.3% of graph-unique HLA-DQA1
reads are unaligned under linear, not aligned elsewhere.

### A human graph index silently discards ~19% of its variants (`analysis/hap`)

`hisat2-build` halves a window's variant set whenever the local graph exceeds
the edge budget, and reports nothing. Measured on GRCh38 chr1:

| variant set | instances | retained |
|---|---|---|
| dbSNP b157 common (gnomAD/TOPMED/1000G) | 1,985,316 | 80.9% |
| dbSNP b157, 1000G 30x only | 1,731,945 | 86.7% |
| 1000G 30x phased | 1,204,387 | 92.4% |

`hisat2-inspect --snp` reports every variant that was supplied regardless, so
there is no way to detect this from a finished index. Filed upstream as
[hisat2#473](https://github.com/DaehwanKimLab/hisat2/issues/473).

Windows explode 1.15–1.38 times each, so the first halving overshoots. Replacing
it with a binary search for the largest subset that fits recovers most of the
loss with no format change:

| chr1 variant set | halving | binary search |
|---|---|---|
| dbSNP b157 common | 80.9% | **87.5%** |
| 1000G 30x phased | 92.4% | **97.8%** |

### The 2^32 node bound is governed by haplotypes, not variants (`analysis/aws`)

A whole-genome graph index over the 1000 Genomes phased panel does not build.
Real phasing carries ~1.11 haplotypes per variant against ~0.78 for the greedy
colouring partition every distributed index is built from — about 46% more
distinct paths. Cutting variants 15% moved haplotypes only 14% and did not help,
because the reference alone accounts for ~72% of the 32-bit path-node space
before any variant is added.

| attempt | variants | haplotypes | outcome |
|---|---|---|---|
| 32-bit | 14.95M | 16.34M | exceeded integer bounds |
| 64-bit | 14.95M | 16.34M | OOM at 368.8 GB |
| 32-bit | 12.64M SNVs | 14.01M | exceeded integer bounds |
| 64-bit | 12.64M SNVs | 14.01M | OOM at 368.8 GB |

Both 64-bit runs were killed at the machine ceiling, so the requirement is
unbounded above 368.8 GB rather than measured. This is a limitation of HISAT2's
unbounded prefix-doubling construction; GCSA2 solves the same problem with
order-bounded, external-memory construction.

### Where the advantage does not appear

Recorded because they bound the claim:

- **Genotyping accuracy (`analysis/lcl`)** — graph 0.91% error vs linear 0.85%.
  The graph makes 3x more false het calls (0.247% vs 0.081%). Variant-aware
  alignment helps allele *fractions*, not genotype calling.
- **Pseudogene variant filtering (`analysis/filt`)** — the hypothesis that
  paralogous-sequence variants caused ribosomal-protein gene loss was wrong;
  removing them recovered 9.7%. The reads are spliced, and processed
  pseudogenes have no introns.
- **Memory (`analysis/t3`)** — the floor is 8 GB at `-p 2`, 10 GB at `-p 8`;
  wall time degrades 13.3x at 10 GB.

## Layout

```
analysis/human       T1 (HLA recovery) and T2 (allele-specific expression), 3' chemistry
analysis/fiveprime   5' replication, peptide-binding-domain enrichment
analysis/lcl         T4b, genotype error against GIAB HG001
analysis/t3          memory-cap titration
analysis/filt        pseudogene variant filter and its falsification
analysis/cmp         three-way comparison, HISAT2 / STAR / rustar
analysis/hap         variant-retention A/B: halving vs binary search
analysis/aws         cloud build stages and reports (chr1 probe, whole-genome attempts)
upstream/            the two upstream issue writeups
```

## Caveats

- Several scripts were recovered from the session transcript after the working
  directory was cleared; they are the last version written, and their primary
  outputs (large matrices, FASTQs, indexes) are not preserved here. The reports
  carry the numbers.
- Absolute timings come from a mix of an Apple M5 Pro laptop and AWS r7i
  instances and are not comparable across those; within-table comparisons were
  run on one machine.
