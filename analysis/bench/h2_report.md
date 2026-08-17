# H2 — does variant-aware recovery scale with distance from the reference?

**Answer: directionally yes, but this pair does not establish it (p = 0.092),
and one gene carries most of the difference. What the experiment does establish
is that the HLA benefit is large and reproducible in *both* donors.**

Two lymphoblastoid lines from GSE126321, 50,000,000 reads each, aligned against
the `grch38_snp` graph index and the linear `grch38` index. Same reads, same
gene model, same machine; the index is the only variable.

| line | 1000G | population | why |
|---|---|---|---|
| GM12878 | NA12878 | CEU / European | close to GRCh38, which is largely European-derived |
| GM18502 | NA18502 | YRI / Yoruba | the ancestry contrast |

Reads are 150 bp paired RNA-seq with no barcode read, so this is pseudobulk —
the same treatment T4b used, which keeps the two directly comparable. Counts are
uniquely-mapped reads (`NH:i:1`) assigned to exactly one gene.

## The pipeline reproduces the earlier work exactly

All four alignment rates match T4b to the second decimal:

| | graph | linear |
|---|---|---|
| GM12878 | 88.30% | 88.37% |
| GM18502 | 87.76% | 87.87% |

## The baseline is not 1.00

The graph index **loses** uniquely-assigned reads overall:

| donor | graph | linear | ratio |
|---|---|---|---|
| GM12878 | 30,991,258 | 33,297,327 | **0.9307** |
| GM18502 | 30,890,657 | 33,429,802 | **0.9240** |

Alternate alleles give more loci a chance to match, so some reads that were
unique under the linear index become multimappers and are excluded. Any HLA gain
must be read against ~0.93, not against 1.00 — which makes the gains larger, not
smaller. Housekeeping controls agree: B2M 0.995/0.981, ACTB 0.998/1.001,
GAPDH 0.965/0.991, TMSB4X 0.971/0.971 (mean 0.982 / 0.986).

## The ribosomal-pseudogene effect, for the third time

Two genes were tried as controls and rejected:

| gene | GM12878 | GM18502 |
|---|---|---|
| RPL13A | **0.371** | **0.391** |
| EEF1A1 | **0.567** | **0.494** |

Processed pseudogenes of ribosomal-protein and translation-factor genes turn
these reads into multimappers under the graph index. This is the same mechanism
that produced H1's `OLFM3`/`RPSAP19` disagreement with CellRanger, and the same
one behind the mouse three-way comparison's ribosomal-protein disagreements.
**Three independent datasets, one cause.** It is a real cost of the graph index
under a unique-reads-only counting rule, and it is not small.

## HLA, per donor

Ratios below are graph/linear normalised by each donor's transcriptome-wide
ratio above.

| gene | CEU | YRI |
|---|---|---|
| HLA-DRB1 | 1.782 | **2.601** |
| HLA-C | 1.706 | **1.918** |
| HLA-DQB1 | **1.909** | 1.863 |
| HLA-DQA1 | **1.556** | 1.340 |
| HLA-DPB1 | 1.378 | **1.493** |
| HLA-DPA1 | 1.167 | 1.169 |
| HLA-A | 1.069 | 1.157 |
| HLA-B | 1.135 | 1.145 |
| HLA-E | 1.060 | 1.086 |
| HLA-F | 1.006 | 1.089 |
| HLA-DMA | 1.072 | 1.081 |
| HLA-DMB | 1.074 | 1.081 |
| HLA-DRA | 1.074 | 1.070 |

Below the >=100-read cutoff but worth recording, because they show how
donor-specific this is: **HLA-DRB5 gains 8.5x in both donors** (324 vs 38, and
425 vs 50), while **HLA-DQA2 gains 8.7x in CEU and *loses* (0.68x) in YRI**.

The pattern replicates T1 on a third and fourth donor and a different library
type: class II DR/DQ and class I C transformed, while HLA-DRA, HLA-DMA and
HLA-DMB — the invariant chains — sit flat at 1.07. That those three land within
0.01 of each other in both donors is the strongest internal check here.

## The comparison H2 exists to make

| | CEU | YRI |
|---|---|---|
| mean HLA graph/linear, raw | 1.2162 | 1.2862 |
| housekeeping-normalised | 1.2383 | 1.3042 |
| **transcriptome-normalised** | **1.3067** | **1.3919** |
| genes where YRI gains more | \- | **10 / 13** |

**Two-sided exact sign test: p = 0.092.** Mean paired difference +0.085.

That is a consistent direction and a conventional miss. Worse for the
hypothesis, the mean is not evenly spread: **excluding HLA-DRB1 alone, the mean
difference falls from +0.085 to +0.024** across the remaining 12 genes. HLA-DQA1
runs the other way entirely (CEU 1.556, YRI 1.340).

So the honest statement is:

- **Established:** variant-aware alignment recovers HLA reads substantially in
  both donors — 1.31x and 1.39x on average across 13 HLA genes, up to 2.6x for
  a single gene, against flat invariant controls. This is the T1 result
  replicated on two more donors.
- **Not established:** that the benefit scales with ancestry distance. The
  direction is right in 10 of 13 genes but p = 0.092, one gene dominates, and
  n = 1 donor pair. A cohort of ~10 donors per population would settle it; two
  lines cannot.

The competing explanation remains the one T1 already raised and this experiment
strengthens: **which HLA genes benefit is a property of the individual's
haplotypes, not of their population label.** HLA-DQA2 gaining 8.7x in one donor
and losing in the other, from the same reference and the same code, is hard to
explain any other way. Ancestry may shift the *distribution* of how divergent a
donor's haplotypes are, which would produce exactly the weak-but-present trend
seen here, but the per-donor variance swamps it at n = 2.

## Reproducing

```
h2_grab.sh              # stream 50M reads per donor from ENA (SRP184746)
h2_run.sh               # 4 alignments, counts straight from GX/GN tags, no SAM stored
h2_analyze.py lcl/      # controls, HLA table, sign test
```

## Caveats

- Pseudobulk, not single-cell: these runs carry no barcodes. The single-cell
  version of this result is T1 (`analysis/human`), on a different donor.
- `--gene-strand Unstranded`, since the library's strandedness is not declared.
  It applies identically to both indexes, so it cannot create a graph/linear
  difference.
- Multi-gene reads are excluded rather than distributed. That is what makes the
  ribosomal-pseudogene effect visible; a `--solo-multi-mappers EM` run would
  redistribute those reads and is a separate experiment.
- The 13-gene cutoff (>=100 linear reads in both donors) excludes HLA-DRB5 and
  HLA-DQA2, the two most donor-specific genes in the panel. Including them would
  raise both means and widen the spread, not narrow it.
