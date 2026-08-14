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

### Variant-aware alignment reduces reference bias (`analysis/human`, `analysis/fiveprime`, `analysis/bias`)

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

**But retention barely affects alignment** (`analysis/bias`). A dose-response
over 200,000 matched read pairs on chr1 moves the ALT fraction 0.50080 (75.1%
retention) to 0.50087 (87.5%) — 34 sites in 200,000. The same test puts linear
at 0.02861 against the graph's 0.50080, an effect ~6,700x larger. At 75.1%
retention the index is already unbiased, because the global ALT list keeps every
variant and only the local indexes lose them. So the patch is justified by
"users get the variants they supplied" and by the reporting bug, **not** by
better alignment — which is what was reported upstream.

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

### Parity and throughput (`analysis/cmp`)

Mouse 10x, 10M reads, same machine and input: HISAT2-solo **51.2 s / 5.18 GB**
against rustar's 52.1–313 s / 25.6–28.2 GB and STARsolo's ~28.3 GB. rustar's
6x runtime spread is page-cache thrash from a 25 GB index on a 51.5 GB machine;
HISAT2 varied 0.6 s and took 78 page faults on its best run. Concordance with
STARsolo is Jaccard 0.944 / per-gene r 0.977, against a same-family control
ceiling of 0.983 / 0.9999.

Gene and GeneFull now come from a single alignment pass (1.90x over two runs on
a 111,600-read fixture). The three-way report predates that change and is
annotated where superseded. Its inputs are in `s3://rustar-bench/` — FASTQs,
gene model, GTF, whitelist, STAR and rustar indexes — so it can be re-run after
rebuilding the mouse HISAT2 index from a public GRCm39 FASTA.

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
analysis/bias        reference-bias dose-response across retention levels
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
