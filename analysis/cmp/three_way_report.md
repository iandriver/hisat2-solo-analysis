# HISAT2-solo vs STARsolo vs rustar — mouse 10x, 10M reads

> **Snapshot, partially superseded.** Measured against `solo/gene-model` @
> cc0e998. Finding 1 (no single-pass Gene+GeneFull) has since been fixed and is
> struck through below; every other number here still stands as measured.
>
> **Re-running it is possible.** The input FASTQs, the `.ht2gm` gene model, the
> GTF, the whitelist and the STAR/rustar indexes are all in
> `s3://rustar-bench/` (`fastq/`, `hisat2-bench/refs/`, `idx/`). Only the GRCm39
> genome FASTA — a public download — and the HISAT2 mouse index built from it
> are absent, so reproducing this needs one index rebuild and nothing else.

## Setup

Identical input for all three: `5k_Mouse_PBMCs_5p_gem-x_GEX_S1_L001` R1+R2, 10,000,000
reads, the same 3.7M-barcode 5'GEX whitelist, and the same 33,696-gene GRCm39
annotation. CB 1–16, UMI 17–28, `--soloStrand Reverse`, `1MM_CR` UMI dedup,
16 threads.

Machine: Apple M5 Pro, 51.5 GB RAM, 18 cores (6 "Super" + 12 "Performance"),
native arm64.

| Tool | Build | Provenance |
|---|---|---|
| HISAT2-solo | `solo/gene-model` @ cc0e998 (fork + U1–U9) | run natively, this machine |
| rustar | `origin/main` 6fb99a1 | run natively, this machine |
| STARsolo | 2.7.11b | **pre-existing output**, x86 container; not runnable natively here |

## Runtime and memory

Native, same machine, same input, 16 threads.

| Tool | Wall | Peak RSS | Page faults | Features per pass |
|---|---|---|---|---|
| HISAT2-solo | **51.2 / 51.5 / 51.8 s** | **5.18 GB** | **78–576** | 1 at the time; 2 now, see below |
| rustar | 52.1 s … 95.0 s | 25.6–28.2 GB | 1.5–3.9 M | 2 |
| STARsolo | not runnable natively | (28.3 GB reported previously) | — | 2 |

Two things matter more than the headline numbers.

**~~HISAT2 needs two passes for Gene + GeneFull.~~ Superseded — this gap is
closed.** When this was measured, `--gene-feature` took one of `Gene` or
`GeneFull`, so the common exonic-plus-intronic use case cost ~103 s against
rustar's ~52 s. It now accepts a comma-separated list and counts both features
from a single alignment pass, as STARsolo and rustar do with
`--soloFeatures Gene GeneFull`: per-feature state lives in its own arena and
only the interval query repeats, so the alignment, CIGAR and reference-block
work is shared.

**Measured on this dataset: 1.92x.** The ~103 s two-pass figure is retired. Same
10M reads, same annotation and whitelist, rebuilt index, 16 threads:

| config | 3 runs (s) | median |
|---|---|---|
| Gene | 47.4 / 46.8 / 46.5 | 46.8 |
| GeneFull | 46.7 / 46.8 / 46.8 | 46.8 |
| Gene,GeneFull | 49.2 / 48.8 / 48.6 | 48.8 |

Two passes 93.6 s, one pass 48.8 s, **1.92x**. The second feature costs **+2.0 s,
4.3% over Gene alone** — an interval query and a second record arena, not another
alignment.

**These absolute seconds are not comparable to the 51.2 s above.** They were taken
while another workload held the machine at load 13-20; the original was measured
on a quiet one. The ratio is what survives contention, and three checks say it
did: the runs were rotated Gene -> GeneFull -> Gene,GeneFull so each config
sampled the same conditions, per-rotation speedups were 1.91 / 1.92 / 1.92 while
load climbed from 13.5 to 19.7, and the median and minimum agree exactly. Gene
and GeneFull do near-identical work and came out at 46.8 s each, which is the
built-in check that the sampling was sound.

An earlier attempt with the wrong barcode whitelist put the overhead at 0.2%.
That run had 3.9% valid barcodes, so 96% of reads never reached the counter and
the second feature had nothing to do. With the correct 5' GEM-X whitelist and
3.9M UMIs actually counted, the cost is 4.3%. Correctness at that scale held:
Gene 3,907,687 UMIs < GeneFull 4,624,940, and both combined-pass matrices were
byte-identical to their single-feature runs. Alignment rate 87.40% and index size
3.8 GB both reproduce the original run exactly.

Correctness of the combined pass is covered in the test suite: each matrix is
byte-identical to the one a single-feature run produces, output does not depend
on the order features are listed, and `sum(Gene) < sum(GeneFull)` is asserted
strictly so the identity checks cannot pass vacuously on a fixture where the two
features coincide.

**rustar's runtime is unstable on this machine, HISAT2's is not.** rustar was
measured at 52 s, 95 s, 96 s and 313 s on identical inputs; the 313 s run showed
16.1M page faults and 30.4M involuntary context switches. Its 25 GB index sits at
the edge of 51.5 GB of RAM, so wall time swings ~6x with page-cache state. HISAT2
varied by 0.6 s across runs and took **78 page faults** on its best run. Run
back-to-back under matched conditions the two are the same speed (96.6 s vs
95.0 s). The variance is the finding, not the mean: this is what a 5 GB index buys.

## Concordance — Gene

| Pair | Jaccard | per-cell UMI r | per-gene UMI r | depth ratio |
|---|---|---|---|---|
| STARsolo ↔ rustar | 0.98326 | 0.99999 | 0.99992 | 0.9991 |
| HISAT2 ↔ STARsolo | 0.94425 | 0.99993 | 0.97707 | 0.9606 |
| HISAT2 ↔ rustar | 0.94250 | 0.99993 | 0.97759 | 0.9597 |

## Concordance — GeneFull

| Pair | Jaccard | per-cell UMI r | per-gene UMI r | depth ratio |
|---|---|---|---|---|
| STARsolo ↔ rustar | 0.98904 | 1.00000 | 0.99988 | 0.9969 |
| HISAT2 ↔ STARsolo | 0.93042 | 0.99981 | 0.97907 | 0.9315 |

The STARsolo↔rustar row is the control, and it is the most useful number in the
table: two tools of the same algorithmic family agree at Jaccard 0.983 and
per-gene r 0.9999. That is the ceiling. HISAT2 sits below it because it aligns to
a different index with different scoring, which is expected and not a defect —
but it does bound how much of the gap is attributable to anything fixable.

## Called cells

| Pair | Jaccard | shared |
|---|---|---|
| STARsolo ↔ rustar | 0.99864 | 3,674 |
| HISAT2 ↔ rustar | 0.99647 | 3,674 |
| HISAT2 ↔ STARsolo | 0.99620 | 3,673 |

All three call the same ~3,674 cells. Downstream, the three matrices are
interchangeable at the level of which cells exist.

## Where HISAT2's missing 4% comes from

| Tool | unique | multi | total mapped |
|---|---|---|---|
| STARsolo | 75.26% | 16.69% | 91.95% |
| rustar | 75.23% | 16.81% | 92.04% |
| HISAT2-solo | 72.41% | 14.99% | **87.40%** |

HISAT2 maps 4.6 points fewer reads, which accounts for essentially all of the
0.96 Gene depth ratio. Per-gene medians are 1.0000 (vs STARsolo) with IQR
0.99–1.01, and 95.6% of genes with ≥100 UMIs fall within 10% of the overall depth
ratio — so the shortfall is a uniform depth effect, not a systematic bias against
particular genes.

GeneFull is worse (0.9315 vs 0.9606). Intronic and intergenic-adjacent reads are
where HISAT2's graph index and scoring differ most from a linear-reference
suffix-array search, and GeneFull is precisely the feature that depends on them.
Worth a look if GeneFull is a target use case.

## Findings

**1. ~~HISAT2-solo cannot emit Gene and GeneFull in one pass.~~ Fixed.**
`--gene-feature` now takes a comma-separated list and produces both features
from one alignment pass; per-feature state is separate and only the interval
query repeats. Measured at 1.90x over two separate runs on a 111,600-read
fixture, and covered by the test suite (byte-identical to single-feature runs,
order-independent, `sum(Gene) < sum(GeneFull)` asserted strictly). Not
re-measured on this dataset — see the runtime section.

**2. STAR is not runnable natively on this machine.** Two separate faults: the
prebuilt index resolves its transcript info to `/geneInfo.tab` (an empty
directory prefix, from `sjdbInsertSave Basic` in the index metadata — patching
that metadata does not help), and `genomeGenerate` was previously found to write
a 0-byte SA on this arm64 build. So a native three-way runtime comparison is not
possible here; STARsolo appears above on outputs only.

**3. rustar: the checked-out worktree is on an older branch than `origin/main`,
and the difference is a real counting bug that `origin/main` has already fixed.**

The local worktree is at 5c020ff (Jul 11); `origin/main` is 6fb99a1 (Jul 26), and
they have diverged — 44 commits on the local line are not in main, and main's
solo fixes are not in the local line. Built from the worktree, rustar produces
2,650,591 non-zero Gene entries; built from `origin/main`, 2,543,293 — a 4.2%
difference. GeneFull is unaffected (3,201,700 vs 3,202,076, 0.01%).

The cause is in `src/solo/gene.rs`. The worktree branch assigns a read to a gene
on **exon overlap**; `origin/main` requires **containment**:

```rust
// STAR's Gene feature requires the alignment to be *concordant* with the gene,
// not merely to overlap an exon: every aligned block must lie wholly within the
// gene's exons. [...] Pure exon-overlap over-assigns such boundary reads — the
// measured rustar-vs-STAR divergence.
let concordant = tr.exons.iter()
    .all(|b| gene_ann.block_is_exonic(g, b.genome_start, b.genome_end));
```

This is already diagnosed and fixed upstream, so it is not a new bug — but
anything built from the current worktree will over-count Gene by ~4% and detect
~600 spurious genes. Measured against STARsolo: the worktree build gives a depth
ratio of 0.9641 and Jaccard 0.94575, the `origin/main` build 0.9991 and 0.98326.

This is the same containment-vs-overlap distinction HISAT2-solo was validated
against in S1, where `GX:Z` matched `htseq-count -m intersection-strict`
(containment) on 100.000% of uniquely-mapped reads but only 96.7% against `union`
(overlap). Two independent implementations, same trap.

**4. The June container benchmark used the pre-fix rustar.** Its recorded Gene
matrix (2,650,612 non-zeros) matches the worktree build, not `origin/main`. Any
conclusion drawn from those numbers about rustar-vs-STAR Gene agreement is
measuring the bug.

## Caveats

- STARsolo numbers are outputs only; it was not re-run.
- STARsolo used `--clipAdapterType CellRanger4 --outFilterScoreMin 30`; HISAT2 has
  no equivalent, and rustar was run without them. This is part of the residual
  difference and is not separated out here.
- STARsolo's `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts` is more
  permissive than HISAT2's `--solo-cb-match 1MM`.
- Human was not attempted: no human FASTQs, reference, or index of any kind are
  present locally, and STAR could not participate regardless.
