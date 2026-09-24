# hisat2-solo — analysis

Measurement scripts, reports and result tables behind the HISAT2 work in
[iandriver/hisat2](https://github.com/iandriver/hisat2). Kept out of that repo so
the source tree stays clean for a possible pull request.

The code being measured lives on two branches there:

| branch | what |
|---|---|
| `upstream/modernize` | fixes that apply to HISAT2 generally — Python 3, C++17, `std::thread`, zlib, BAM output, build portability, variant retention |
| `solo/gene-model` | the above plus single-cell: gene model, barcode/UMI, counting, allele-specific counting |

## The index builder (`design/rust/ht2fmt`)

`hisat2-build` needs the whole path graph resident. For a human graph index
that measured 671 GB of RSS, which is the real reason variant-aware references
stay unused: you cannot build one on hardware you have.

`ht2wg` runs the same construction against disk. It is a Rust reimplementation
of HISAT2's index format, written against `gfm.h` and `hgfm.h` rather than
copied from them, and it is checked by byte-equality rather than by inspection.

```sh
cd design/rust/ht2fmt && cargo build --release      # no dependencies
cd .. && ./mkfixtures.sh /tmp/fix ../../../hisat2   # derive the fixtures
     ./verify.sh /tmp/fix ../../../hisat2           # ALL FIXTURES BYTE-IDENTICAL
```

A whole-genome build then looks like:

```sh
HT2_LARGE=1 HT2_THREADS=8 \
  ht2wg genome.fa genome.snp genome.haplotype scratchdir out/genome 134217728
```

Peak RAM is bounded by `threads * budget * 18 bytes` plus about 4 GB for the
reference stage, so roughly 19 GB at the settings above. Scratch is the cost
instead: about 325 GB at the `generateEdges` peak, not at the generation-11
node peak, so size the disk from the emit stage. It resumes by being re-run.

`HT2_SS` and `HT2_EXON` add `--ss`/`--exon` annotation. Splice sites and exons
enter the ALT table, splice sites add one backbone edge each, and exons never
reach the graph. On GENCODE v32 that is 382,106 junctions and 324,668 exons,
and it takes the longest uncuttable fragment from 1.06 Mb to 3.06 Mb without
changing how many doublings the build needs.

Verified against `hisat2-build` on the 1 Mb example reference and on chr22 with
real GENCODE v32: all eight index files byte-identical, and the doubling curve
matches generation for generation.

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

### Variant-aware alignment recovers X-inactivation escape genes linear misses (`analysis/xci`)

Escape is detected as biallelic expression at heterozygous sites — an
allele-fraction measurement, and the published catalogues are built on linear
alignment. A graph index carrying chrX variants recovers **ZFX** (European donor)
and **EIF1AX** (Yoruba donor), both established escapees, and loses none.

The mechanism control is exact. An index with **no non-PAR chrX variants** moves
the allele fraction by 0.0001–0.001 (|z| ≤ 1.3); the same index with them added
moves it by 0.006–0.013 (z = +5.9 to +8.6). Sites move up and essentially never
down — 1 and 5 down against 87 and 61 up.

This needed chrX repaired first: the panel carried 184 variants/Mb on chrX
against an autosome mean of 5,158, all in PAR1, because
`hisat2_extract_snps_haplotypes_VCF.py` hardcodes diploid genotypes and cannot
read the haploid male calls that non-PAR chrX carries in a phased panel. It
crashes *after* writing what it processed, leaving a plausible file holding PAR1
only.

Thin, and honestly so: 1–2 genes per donor, and the 1000 Genomes chrX genotypes
carry enough error that 26–35% (CEU) and 72–76% (YRI) of "het" sites read
monoallelic at high depth. Gene-level aggregation tolerates that; site-level work
does not. A pilot with a working mechanism, not a corrected catalogue.

### The whole-genome index builds in ~11 GB, byte for byte (`design/rust`)

The external builder reproduces that index — all eight files, 8.9 GB, 64-bit —
**byte-identically, in ~11 GB against `hisat2-build`'s measured 671 GB**, in
17.7 h on one external SSD. No cloud instance, no vCPU quota increase. (RSS was
sampled during the run, not maximised over it; the sort budget caps it at 19 GB.)

That answers the construction limit noted above inside HISAT2's own format
rather than by moving to GCSA2: the prefix doubling runs as an external
sort-merge join, so peak memory follows the sort budget instead of the node
count.

Six construction bugs had to be fixed first, none of them reachable below
chromosome scale — every fixture is a single contig with no N runs and no
co-located variants. Two were in the global graph, three in which variants a
local-index window keeps, and one was not a defect at all: `.7`'s haplotype
order comes from an unstable `std::sort` under a comparator that admits ties,
so `hisat2-build`'s own `.7` is not reproducible across standard library
implementations.

### Parity and throughput (`analysis/cmp`)

Mouse 10x, 10M reads, same machine and input: HISAT2-solo **51.2 s / 5.18 GB**
against rustar's 52.1–313 s / 25.6–28.2 GB and STARsolo's ~28.3 GB. rustar's
6x runtime spread is page-cache thrash from a 25 GB index on a 51.5 GB machine;
HISAT2 varied 0.6 s and took 78 page faults on its best run. Concordance with
STARsolo is Jaccard 0.944 / per-gene r 0.977, against a same-family control
ceiling of 0.983 / 0.9999.

Gene and GeneFull now come from a single alignment pass: **1.92x** over two
separate runs on the same 10M mouse reads (93.6 s -> 48.8 s), the second feature
costing 4.3%. The three-way report predates that change and is annotated where
superseded. Its inputs are in `s3://rustar-bench/` — FASTQs, gene model, GTF,
whitelist, STAR and rustar indexes — so it can be re-run after rebuilding the
mouse HISAT2 index from a public GRCm39 FASTA.

### Variant loss concentrates in the MHC (`analysis/probe`)

`hisat2_index_probe.py` (in the hisat2 repo) predicts both of hisat2-build's
failure modes from the `.snp`/`.haplotype` files in seconds. On the whole-genome
phased set it returns **WILL FAIL at 72.5% of the 2^32 ceiling** — the outcome
stage 2b took four builds and ~$24 to establish.

It also shows where the dropped variants would be:

| | over-budget windows | variants at risk |
|---|---|---|
| genome-wide | 19.1% | 28.9% |
| extended MHC | 55.1% | **83.9%** |

The MHC carries 2.4x the genome-wide variant density, 27 of the 100 densest
over-budget windows sit in it (0.18% of windows), and **every HLA gene measured
in the T1 work is in an over-budget window** — the densest window in the genome
covers HLA-DQA1, the gene with the 32x gain.

The constraint and the payoff share an address. But over-budget status does not
predict benefit — those same genes range 1.00x to 32.1x — and the dose-response
above says retention barely moves alignment anyway. The defensible reading is
that the T1 numbers were obtained with a partial variant set in those windows,
making them a floor rather than a ceiling.

### Human 10x benchmark: memory holds, runtime unmeasurable here (`analysis/bench`)

PBMC 1k v3, 66.6M read pairs, GRCh38 + GENCODE v50, against CellRanger's own
published output for the same sample.

| | HISAT2-solo (graph) | STAR (linear) |
|---|---|---|
| index on disk | **6.5 GB** | 29 GB |
| peak RSS, mapping | **9.04-9.43 GB** | 31.3 GB |
| index build | impossible at this scale | 2 h 26 m, 30.6 GB |

Concordance with CellRanger: **cell Jaccard 0.9185**, per-cell UMI totals
r = 0.9950, per-gene r = 0.9550, 94.4% of genes within 2x.

Two things the benchmark turned up that matter more than the ratios:

- **Pointing HISAT2-solo at a raw GENCODE GTF silently discards ~10% of reads.**
  Full v50 (78,941 genes) leaves 10.04% of reads multi-gene and dropped; the
  CellRanger-filtered set (32,364) leaves 1.58%. This is why CellRanger ships a
  filtered reference.
- **Processed pseudogenes of ribosomal proteins are the recurring disagreement.**
  `OLFM3` scores 12,992 UMIs against CellRanger's 0 because `RPSAP19` — an RPSA
  pseudogene in neither reference — sits entirely inside OLFM3's first exon. The
  mouse three-way comparison hit the same mechanism from the other direction.

No runtime claim is made: three identical runs spanned 431.9 s to 4,398.7 s
under macOS Spotlight load, while peak RSS moved less than 4%. STARsolo itself
could not be run — the Homebrew ARM64 build fails nondeterministically resolving
`geneInfo.tab`, and one run mapped 0 reads while printing `ALL DONE!`.

### HLA recovery replicates on two more donors; ancestry scaling does not resolve (`analysis/bench`)

GM12878 (CEU) and GM18502 (YRI), 50M reads each, graph vs linear index, all four
alignment rates reproducing T4b exactly.

Normalised by each donor's transcriptome-wide graph/linear ratio (0.931 / 0.924
— the graph *loses* ~7% of unique reads overall, so 1.00 is the wrong baseline):

| gene | CEU | YRI |
|---|---|---|
| HLA-DRB1 | 1.78 | **2.60** |
| HLA-C | 1.71 | 1.92 |
| HLA-DQB1 | 1.91 | 1.86 |
| HLA-DRA / DMA / DMB (invariant) | 1.07 | 1.07-1.08 |

**Established:** the T1 HLA result replicates on two further donors — mean 1.31x
(CEU) and 1.39x across 13 HLA genes, against invariant chains flat at 1.07.

**Not established:** that the benefit tracks ancestry. YRI gains more in 10 of
13 genes, but the sign test gives **p = 0.092**, and dropping HLA-DRB1 alone
collapses the mean difference from +0.085 to +0.024. HLA-DQA2 gains 8.7x in CEU
and *loses* in YRI. Which genes benefit still looks like a property of the
individual's haplotypes rather than their population label — as T1 suspected.
Settling it needs a cohort, not a pair.

Third appearance of one mechanism: RPL13A (0.37/0.39) and EEF1A1 (0.57/0.49)
lose reads under the graph index, because processed pseudogenes turn them into
multimappers. Same cause as H1's OLFM3/RPSAP19 and the mouse comparison's
ribosomal-protein disagreements.

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
analysis/bench       human benchmarks: footprint vs STAR, concordance vs CellRanger,
                     HLA recovery across two ancestries (H1, H2)
analysis/probe       build-outcome prediction; where variant loss lands
analysis/xci         X-inactivation escape, linear vs graph; the chrX panel repair
analysis/aws         cloud build stages and reports (chr1 probe, whole-genome attempts)
design/rust          external index builder: byte-identical whole-genome construction
design/e3            the 64-bit whole-genome C++ build this is verified against
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
