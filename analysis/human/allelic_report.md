# SNP-aware index build cost, and `--solo-allelic` on real data

## 1. Why the mouse dataset could not answer the question

Before building anything, I checked whether the 5k mouse PBMC sample actually
carries non-reference alleles. It does not.

| Evidence | Result |
|---|---|
| Reads aligning with zero mismatches | **90.3%** (NM=0), 6.9% NM=1 |
| chr19 sites ≥20x with non-ref fraction >0.8 | 2 of 118,464 (**0.002%**) |
| chr11 sites ≥20x with non-ref fraction >0.8 | 20 of 246,970 (**0.008%**) |

A divergent strain shows on the order of 0.5%. This is C57BL/6 — the reference
strain itself — so allelic counting on it would correctly find nothing. That is
why the demonstration moved to human, where individuals are naturally
heterozygous.

## 2. SNP-aware index build cost, measured

Built on mouse chr19 (61,420,004 bp, 2.26% of the genome), `-p 8`, on this
machine (M5 Pro, 51.5 GB RAM).

| SNP set | SNVs supplied | SNVs in index | Wall | Peak RSS | "graph exploded" warnings |
|---|---|---|---|---|---|
| none (linear) | 0 | 0 | 7.2 s | **0.43 GB** | 0 |
| 1 per 214 bp | 287,009 | 287,009 | 25.1 s | **12.43 GB** | 5 |
| 1 per 30 bp (all Ensembl) | 2,039,608 | **2** | 41.7 s | 19.62 GB | 2,285 |

The middle row matches the SNP density of the published human index
(14.46M SNPs over 3.1 Gb ≈ 1 per 214 bp), so it is the realistic case.
**Variants cost about 29x the memory of a linear build** of the same sequence.

### A silent failure worth knowing about

Feeding the *full* Ensembl mouse variant set — which pools all strains, giving
one variant every 30 bp — produces an index containing **2 SNPs**. The build
exits 0 and writes a normal-looking index. It emits 2,285
`Warning: a local graph exploded` lines, but nothing fails, and
`hisat2-inspect --snp` is the only way to discover the index is effectively
linear.

Anyone pointing `hisat2-build --snp` at a whole dbSNP or Ensembl VCF without
filtering to common variants will get a "SNP-aware" index with no SNPs in it,
and no error to tell them.

### Whole-genome extrapolation

I would not extrapolate the 12.43 GB linearly. Scaling by genome size and SNP
count (~45x each) gives ~550 GB, against the ~160 GB HISAT2 documents for
human — so the builder's memory is clearly sublinear, partly bounded by the
blockwise suffix-array construction. What the measurement does establish is
that a whole-genome SNP-aware build is out of reach on 51.5 GB of RAM, which is
why JHU ships prebuilt indexes (`grch38_snp`, 5.0 GB; `grcm38_snp`, 3.8 GB).
Those are what a user should actually use, and what the rest of this used.

## 3. `--solo-allelic` on real human data

10x PBMC 1k v3, both lanes, 66.6M reads, against the prebuilt GRCh38 SNP-aware
index (14,460,407 SNPs) and an Ensembl GRCh38 gene model (62,754 genes).
Strand confirmed empirically as `Forward` (42.3% unique-gene vs 5.0% for
`Reverse`).

| | Wall | Peak RSS |
|---|---|---|
| with `--solo-allelic` | 475.5 s | 12.11 GB |
| without | 382.3 s | 10.00 GB |
| **overhead** | **+93 s (+24.4%)** | **+2.1 GB (+21%)** |

The Gene matrix is **byte-identical** with and without the flag. This is the
real-data cost that was previously untested.

**Output:** 1,134 cells, 23,896 genes, 8.9M UMIs, and

```
Variants Observed,1241667
Allelic UMIs: Reference,9480076
Allelic UMIs: Alternate,1665132
```

### The counts are real genotype, not noise

Aggregated per variant across cells, the allele fraction is trimodal — the
signature of a diploid genome — and the proportions are stable as depth
increases, which noise would not be:

| variants with ≥N UMIs | N=5 | N=10 | N=20 |
|---|---|---|---|
| n | 252,821 | 125,374 | 57,665 |
| AF<0.1 (hom-ref) | 72.5% | 73.3% | 73.6% |
| 0.2–0.8 (het) | 15.7% | 15.8% | 15.7% |
| AF>0.9 (hom-alt) | 8.8% | 8.7% | 8.7% |
| median AF of het sites | 0.474 | 0.478 | 0.481 |

If these were sequencing errors the "het" class would shrink toward zero as
depth rose, and the median would sit near 0 rather than near 0.5.

### Reference bias

Measured at read level from a pileup of the graph alignment, at the 6,971
heterozygous sites with ≥20 reads:

- **median ALT fraction = 0.5000**
- mean = 0.4969

That is an unbiased allelic ratio. Reads carrying the alternate allele are
recovered at the same rate as reference-carrying reads, which is what aligning
to the graph is supposed to buy — the alternate path costs no mismatch penalty.

## 4. What is missing

The direct comparison against a **linear** GRCh38 index — same reads, same
sites, aligner as the only variable — did not complete. The prebuilt linear
index downloaded and extracted, but partway through this work my shell lost
read access to the repository directory:

```
$ ls /Users/iandriver/Downloads/hisat2
ls: /Users/iandriver/Downloads/hisat2: Operation not permitted
```

so `hisat2` can no longer be invoked. Everything above was already on disk and
is unaffected. To finish that comparison the directory permission needs
restoring; the linear index, the het-site list with coordinates, and the graph
pileup are all staged and ready, so it is one alignment plus one pileup.

Without it, the 0.5000 figure stands on its own as "the graph alignment is
unbiased" but cannot be quoted as "better than linear by X" — linear aligners
are commonly reported at 0.45–0.47 on this metric, but I have not measured it
here and will not claim it.
