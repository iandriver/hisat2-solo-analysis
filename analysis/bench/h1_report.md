# H1 — HISAT2-solo on human 10x, against CellRanger and against STAR's footprint

10x PBMC 1k v3, **66,601,887 read pairs**, GRCh38. The prebuilt JHU `grch38_snp`
graph index (12.3M SNPs) against GENCODE v50, on an 18-core Apple M5 Pro with
48 GB.

The intent was a head-to-head against STARsolo. Half of that survived contact
with the machine; the half that did not is documented below rather than
patched over.

## What is measured, and what is not

| | status |
|---|---|
| memory and index footprint | **measured, stable, trustworthy** |
| count concordance vs CellRanger | **measured** |
| wall-clock runtime | **not trustworthy on this machine** — see "Timing" |
| STARsolo count matrix | **not obtained** — STAR would not run reliably here |

## Memory and index footprint

This is the claim the project rests on, and it holds.

| | HISAT2-solo (graph) | STAR / STARsolo (linear) |
|---|---|---|
| index on disk | **6.5 GB** | **29 GB** |
| peak RSS, mapping | **9.04–9.43 GB** | **31.3 GB** |
| peak RSS, index build | n/a (see below) | 30.6 GB |
| index build wall | n/a | 8,776 s (2 h 26 m) |

**≈3.4x less memory to run, ≈4.5x smaller on disk.** The HISAT2 figure was
stable at 9.04, 9.06 and 9.43 GB across three runs whose wall times spanned 10x,
so it is not an artifact of machine load.

The STAR mapping RSS is legitimate even though no STAR run produced counts:
STAR's footprint is dominated by loading the index, which every attempt did
complete.

**The build row is blank for HISAT2 for a reason, and it cuts against us.**
A human *graph* index cannot be built at this scale at all — stage 2b established
that over four attempts and ~$24 (32-bit `exceeded integer bounds`; 64-bit
OOM-killed above 368.8 GB). Users take the prebuilt JHU index or go without.
STAR's index is 4.5x larger and takes 2.5 hours, but you can rebuild it yourself
for any assembly or annotation. HISAT2 wins the footprint comparison decisively
and loses the flexibility one outright.

## Concordance with CellRanger

CellRanger's own published output for this exact sample is the comparator
(`pbmc_1k_v3_filtered_feature_bc_matrix`, 1,222 cells). CellRanger uses STAR
internally and is what the field benchmarks against, so this substitutes for the
STARsolo run that could not be completed.

Annotation restricted to CellRanger's own 33,538 gene ids (32,364 survive in
v50) so that the annotation stops being a variable:

| | HISAT2-solo | CellRanger |
|---|---|---|
| cells called | 1,132 | 1,222 |
| UMIs on shared cells | 8,618,514 | 9,130,347 (-5.6%) |

| agreement over the 1,127 shared cells | |
|---|---|
| **cell Jaccard** | **0.9185** |
| per-cell UMI totals | Pearson(log1p) **0.9950**, Spearman 0.9941 |
| per-gene totals | Pearson(log1p) **0.9550**, Spearman 0.9552 |
| genes within 2x | **94.4%** (10,505 / 11,123 with >=50 UMIs) |

### The annotation is worth more than it looks

Running the *full* GENCODE v50 gene set (78,941 genes) instead of CellRanger's
filtered one costs real signal:

| | full v50 (78,941 genes) | CellRanger set (32,364) |
|---|---|---|
| reads multi-gene, discarded | **10.04%** | **1.58%** |
| reads mapped to a unique gene | 44.10% | 48.24% |
| total UMIs | 9,130,654 | 9,558,062 |
| UMI gap vs CellRanger | -14.8% | **-5.6%** |
| per-gene correlation | 0.9496 | **0.9550** |

Overlapping readthrough and lncRNA annotations make reads ambiguous, and
`--solo-multi-mappers Unique` then throws them away. **This is why CellRanger
ships a filtered reference**, and anyone pointing HISAT2-solo at a raw GENCODE
GTF is silently discarding about a tenth of their reads.

### Where the two still disagree, and why

Neither cause is a defect in HISAT2.

**1. Gene ids drift between annotation versions.** CellRanger 3.0.0 is built on
GENCODE v28; ours is v50. `CAST` (ENSG00000153113) is the clearest case: v50 has
two genes named CAST, an lncRNA at `5:95,962,001-96,631,085` and the
protein-coding gene at `96,247,756-96,779,595`. Our model carries the former
coordinates for that id, CellRanger the latter — hence 1 UMI against 1,938.
Matching gene *ids* does not match their *coordinates*. Same story for PYURF and
FAM89B.

**2. Unannotated pseudogenes sitting inside annotated exons.** `OLFM3` scores
12,992 UMIs here and 0 in CellRanger — impossible for a brain-specific gene in
PBMCs. The locus explains it: **`RPSAP19`, a processed pseudogene of the highly
expressed ribosomal protein gene RPSA, lies at `1:101,786,340-101,787,219`,
entirely inside OLFM3's first exon** (`1:101,786,269-101,787,278`). RPSAP19 is
in neither reference, so reads from RPSA that land on its pseudogene have only
OLFM3 to go to. STAR discards them as multimappers; HISAT2 places them.

That is the **same failure mode the mouse three-way comparison hit**, where the
only genuine count disagreements were ribosomal protein genes. Two species, two
comparators, one mechanism: processed pseudogenes of ribosomal proteins.

## Cross-validation against the earlier human work

Run independently, months apart, and it reproduces:

| | this run | T1/T2 report |
|---|---|---|
| overall alignment rate | 87.20% | 87.20% |
| sequencing saturation | 68.91% | 68.8% |
| HLA-DQA1 UMIs | 2,565 | 2,565 |
| HLA-DRA UMIs | 18,625 | 18,619 |

HLA-DQA1 matching to the UMI is the strongest single check that the pipeline,
gene model and counting path are all behaving as they did before.

## Timing — why there is no number here

Three identical runs, same annotation, same machine, nothing else changed:

| run | wall | peak RSS |
|---|---|---|
| 1 | 431.9 s | 9.06 GB |
| 2 | 2,175.7 s | 9.04 GB |
| 3 | 4,398.7 s | 9.43 GB |

**A 10x spread, monotonically increasing.** The cause is environmental: writing
and deleting ~80 GB on this volume set macOS Spotlight (`mdbulkimport`,
`mds_stores`) and `StorageManagementService` indexing hard, with load average
reaching 24 against 18 cores. No thermal warnings were recorded, so it is
contention rather than throttling. Disabling Spotlight needs system settings
changes that were out of scope.

Earlier Gene-only runs on the same input gave 560 s and 961 s and 1,653 s under
varying load. **No runtime claim should be made from any of these**, and none is
made. The memory numbers are unaffected — they held to within 4% across the same
10x wall spread.

Two side observations that survive the noise, both provisional:

- Adding `GeneFull` to `Gene` cost 2.4x wall (560 s -> 1,359 s) on the full
  78,941-gene set, against the +4.3% measured on mouse with 33,696 genes.
  GeneFull queries gene *bodies* rather than exon unions, so its cost should
  scale with how much annotated gene bodies overlap — much worse on full
  GENCODE. Worth confirming on a quiet machine.
- HISAT2 has no `--no-sam`. Counting runs still format SAM records and discard
  them to `/dev/null`, where STAR has `--outSAMtype None`. The plan specified
  `suppressAlignments` as the lever for exactly this and it was never wired to a
  flag. This is real overhead on the tool's primary use case.

## STARsolo could not be run here

STAR 2.7.11b, Homebrew ARM64 build. `--soloFeatures GeneFull` and eventually
every configuration failed at

```
Transcriptome.cpp:18: could not open input file /geneInfo.tab
```

an empty directory prefix, while `geneInfo.tab` was present in `genomeDir` the
whole time. The behaviour is **nondeterministic**: two early runs succeeded with
byte-identical arguments, then 6 of 6 consecutive attempts failed. Editing the
stored `sjdbGTFfile` in `genomeParameters.txt` changed nothing.

Worse, one run that got past initialisation **mapped 0 reads and printed
`ALL DONE!`**, producing an empty matrix and a plausible-looking 114 s runtime.
Had that landed in a rotation unexamined it would have become a benchmark
number.

Identical inputs producing different resolved paths points at uninitialised
memory rather than configuration, and STAR is not well supported on Apple
Silicon. **No STAR-derived count in this report comes from a completed mapping
run** — only the index build and the genome-load footprint, both of which
completed normally.

## Reproducing

```
prep.sh        # fastqs, references, gene model, index
star_index.sh  # STAR genome (2.5 h, 38 GB headroom)
run3.sh        # three HISAT2 runs against the CellRanger-matched gene set
compare.py <hisat2_filtered_dir> <cellranger_dir> HISAT2 CellRanger
```

`compare.py` needs no numpy/scipy — rank correlation is computed in stdlib, so
it runs against a bare system python.

## Caveats

- The reference given to STAR was dumped out of HISAT2's own index
  (`hisat2-inspect`) so both tools saw byte-identical sequence. Verified: all
  194 sequences match in name and length, chr6 is byte-identical. The one
  divergence found was chrY, where the Ensembl-derived index hard-masks the
  pseudoautosomal regions (2,778,688 more Ns) and GENCODE does not.
- The gene-model `#ref` guard fired on the first human run, catching a GENCODE
  (`chr1`) versus Ensembl (`1`) naming mismatch against the prebuilt index. It
  did its job; the annotation was renamed rather than the guard relaxed.
- Cell counts differ partly because the filters differ in implementation
  (`CellRanger2.2` knee here vs CellRanger 3.0.0's EmptyDrops-based calling),
  not only because of alignment. 95 of CellRanger's cells are absent here.
