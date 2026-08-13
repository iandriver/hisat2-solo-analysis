# Downstream analysis, allelic counting, and StringTie — HISAT2-solo

All on the same 10M reads (5k mouse PBMC 5' GEX, GRCm39, 33,696-gene model),
same whitelist, same parameters.

---

## 1. Side-by-side downstream analysis

Identical scanpy workflow for all three — no per-tool tuning — restricted to the
3,671 barcodes all three called, so clustering is compared on the same cells
rather than on different cell sets.

| | HISAT2 | STARsolo | rustar |
|---|---|---|---|
| barcodes called | 3,684 | 3,676 | 3,677 |
| clusters (leiden 0.5) | 8 | 9 | 8 |

**Clustering agreement (ARI / NMI on the shared cells)**

| Pair | ARI | NMI |
|---|---|---|
| HISAT2 ↔ STARsolo | **0.9315** | 0.9240 |
| HISAT2 ↔ rustar | 0.8925 | 0.8849 |
| STARsolo ↔ rustar | 0.8604 | 0.8901 |

HISAT2 agrees with STARsolo *more* than STARsolo agrees with rustar. Given
rustar is a STAR reimplementation and HISAT2 uses a completely different index,
that is worth stating plainly: at the level a analyst actually works at, the
choice among these three does not change the answer.

**Differential expression.** Wilcoxon markers per cluster, matched between tools
by marker overlap. Mean best-match Jaccard over the top 25 markers: **0.825**
(vs STARsolo), **0.752** (vs rustar). Three clusters match at Jaccard 1.000.
The canonical populations come out the same everywhere — B (Cd79a, Cd19, Bank1),
T (Cd3d/e/g), monocyte/macrophage (Csf1r, Cst3, Ctss), neutrophil (Csf3r, Cxcr2),
cytotoxic (Ccl5, Nkg7).

**Pseudobulk** (sum over shared cells, keyed by gene ID):

| Pair | raw r | log r |
|---|---|---|
| HISAT2 ↔ STARsolo | 0.9786 | 0.9931 |
| HISAT2 ↔ rustar | 0.9791 | 0.9940 |
| STARsolo ↔ rustar | 0.9999 | 0.9992 |

Figures: `comparison.png` (UMAPs, Cd79a on the same cells, pseudobulk scatters),
`markers.png` (per-cluster marker heatmaps), `umap_side_by_side.png`.

### A bug this turned up

The first pseudobulk correlation came out at **-0.29**, which is impossible next
to a per-gene r of 0.977. The cause was a positional join in the analysis
script: HISAT2 wrote `features.tsv` in the gene model's coordinate-sorted order
while STARsolo and rustar write GTF order, so column *i* was not the same gene.

Nothing was incorrect — `features.tsv` ships with the matrix and any reader that
resolves by name was always right — but the layout is meant to be
STARsolo-compatible and anyone comparing matrices reaches for the column index
first. Fixed in `hisat2_extract_genes.py`: gene order is now GTF order, which is
no less deterministic. The 33,696 IDs now match STARsolo's `features.tsv` row
for row, and all 2,496,209 (gene, barcode) counts are unchanged.

---

## 2. `--solo-allelic`

**Does it work: yes, exactly.** Validated against reads whose true REF/ALT
allele is known by construction (exact reference substrings, with the alternate
base substituted at one SNP).

- **160 of 160** (cell, variant, allele) entries match the hand-computed
  expectation. Zero mismatches.
- Totals: 175 REF / 120 ALT observations across 21 variants.

The 21st variant is not an error. `rs147349046` sits 39 bp from a tested SNP, so
100 bp reads span both; carrying no ALT edit there, they are correctly recorded
as REF. That is the documented rule — any tracked variant inside the alignment
without a matching edit was seen as reference.

**What it adds:** `Solo.out/Allelic/raw/` with `ref.mtx`, `alt.mtx`,
`barcodes.tsv` and `features.tsv`, where features are variant IDs (real rsIDs
from the index) rather than genes — a variants × cells pair of matrices. Plus
three `Summary.csv` lines.

**Does it break integer counts: no.**

| matrix | header | non-integer values |
|---|---|---|
| `Gene/raw/matrix.mtx` | `coordinate integer` | 0 |
| `Gene/filtered/matrix.mtx` | `coordinate integer` | 0 |
| `Allelic/raw/ref.mtx` | `coordinate integer` | 0 |
| `Allelic/raw/alt.mtx` | `coordinate integer` | 0 |

The gene matrix is **byte-identical** with and without `--solo-allelic`; the
allelic stream is additive and does not perturb counting. The only fractional
output HISAT2-solo ever writes is `UniqueAndMult-{Uniform,EM}.mtx` under
`--solo-multi-mappers`, which is `coordinate real` by design because split
multimappers are genuinely fractional. The unique-only `matrix.mtx` beside it
stays integer.

### A real bug, found and fixed

`Summary.csv` reported `Variants Observed,0`, `Allelic UMIs: Reference,0` and
`Allelic UMIs: Alternate,0` while the matrices held 21 / 175 / 120. Same for the
velocyto lines.

Cause: the single-pass Gene+GeneFull change moved `writeSummary()` inside the
per-feature loop — correct, each feature needs its own summary — but that put it
ahead of `writeVelocyto()`/`writeAllelic()`, which are what fill in those
tallies. A regression I introduced, not pre-existing. Both writes now run before
the loop. Verified the new test catches it by reintroducing the bug.

### Caveat

`--solo-allelic` needs a SNP-aware index. The GRCm39 index here has none
(`.7.ht2`/`.8.ht2` are empty), so this was validated on the bundled chr22 SNP
index. **Its cost and behaviour on a real dataset are untested** — that needs a
SNP-aware mouse or human index, which is a separate build.

---

## 3. HISAT2 → StringTie

The classic bulk workflow, on the same reads and the same aligner.

| Step | Wall | Peak RSS | Output |
|---|---|---|---|
| `hisat2 --bam` | 47.2 s | 4.71 GB | 708 MB BAM |
| `samtools sort` + index | 4.4 s | — | 388 MB BAM |
| `stringtie -e -G` | 21.9 s | 0.40 GB | gene abundances |
| **total** | **73.5 s** | | ~1.1 GB intermediates |
| *vs* `hisat2 --solo-out-dir` | **51.2 s** | 5.18 GB | cell × gene matrix, no intermediates |

**Agreement at gene level**, all from the same BAM and annotation:

| Pair | log r | Spearman ρ |
|---|---|---|
| solo UMIs ↔ htseq reads | **0.9655** | **0.9680** |
| solo UMIs ↔ StringTie TPM | 0.8797 | 0.8873 |
| htseq reads ↔ StringTie TPM | 0.8653 | 0.8612 |

The third row is the control and it settles the interpretation: **StringTie
disagrees with plain read counting slightly more than it disagrees with solo.**
The ~0.88 is StringTie's model, not anything solo is doing.

That is expected rather than a defect. StringTie is built for full-length
RNA-seq: it estimates transcript-level abundance from coverage and normalises by
length. 10x 5' GEX is tag-based — reads pile at one end of the molecule — so
coverage-based length normalisation is out of domain. The signature is visible
directly: the log ratio of solo UMIs to StringTie-implied reads correlates with
gene length at **r = -0.65**. UMI counts are per-molecule and length-independent;
StringTie's are not.

Solo agreeing with htseq at ρ = 0.968 is the number that matters, and it is
consistent with the earlier S1 result where `GX:Z` matched
`htseq-count -m intersection-strict` on **100.000%** of uniquely-mapping reads.
The residual is UMI deduplication, which is a real biological difference — the
median solo/StringTie ratio of 0.193 is PCR duplicate collapse, not disagreement.

**When to use which.** StringTie remains the right tool for novel transcript
discovery and isoform-level quantification on full-length data; it is not a
single-cell quantifier and has no concept of barcodes or UMIs. For 10x-style
data the solo path is both faster end-to-end and the only one of the two that
produces per-cell output.
