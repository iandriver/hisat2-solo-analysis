# Phase 2: build everything from a current reference, then measure all three

## The workflow, and what is HISAT2-specific

Two steps, both HISAT2's own tooling. Nothing shared with the wider field — no
bcftools, no vg, no external phasing:

```
hisat2_extract_snps_haplotypes_UCSC.py genome.fa snp151Common.txt genome
   -> genome.snp        <id> <type> <chrom> <pos> <allele>
   -> genome.haplotype  <id> <chrom> <left> <right> <comma-separated snp ids>

hisat2-build genome.fa --snp genome.snp --haplotype genome.haplotype genome_snp
```

**The haplotype file is load-bearing**, which this project learned the hard way:
without it the local graphs enumerate the power set of nearby variants instead
of the combinations that actually occur. A 19.2 Mb targeted reference OOMed at
~200 GB with singleton haplotypes. And `hisat2-inspect` can dump `--snp` but
**not** haplotypes, so an existing index cannot be reverse-engineered into a
rebuildable form — the reason HLA had to be dropped from the pseudogene rebuild.

Genotypes/haplotypes are derived by that script from the UCSC table's own
genotype columns, with a fallback assigning genotypes where missing. It is a
HISAT2 heuristic, not real phased data.

## Sources — verified live, and mutually consistent

| | source | size | naming |
|---|---|---|---|
| genome | GENCODE v39 `GRCh38.primary_assembly.genome.fa.gz` | 0.79 GB | `chr1` |
| annotation | GENCODE v39 `primary_assembly.annotation.gtf.gz` | 0.04 GB | `chr1` |
| variants | UCSC `snp151Common.txt.gz` | 0.70 GB | `chr1` |

All three are `chr`-prefixed and agree, so **no name munging** — unlike
`make_grch38_snp.sh`, which had to strip `chr` to match Ensembl.

Two things to decide, flagged rather than assumed:

- **`snp144Common` no longer exists.** UCSC retired it; the URL in the shipped
  build script is dead. `snp151` is denser, so **this will not reproduce the
  shipped `grch38_snp`** — it is a newer, denser index. That is a confound
  against every prior number, and also the most interesting thing here (below).
- **GENCODE v39 is from 2021-12-09; the current release is v50.** v39 is a fine
  pin if it matches other work, but it is not "current". Say the word and it is
  a one-line change.

## What gets built

| artifact | tool | notes |
|---|---|---|
| `genome` (linear) | `hisat2-build` | baseline |
| `genome_snp` | `hisat2-build --snp --haplotype` | **the expensive one** |
| STAR index | `STAR --runMode genomeGenerate --sjdbGTFfile` | GTF baked in |
| rustar index | rustar's own | |
| `genes.ht2gm` | `hisat2_extract_genes.py` | regenerated from the v39 GTF |
| splice sites | `hisat2_extract_splice_sites.py` | fed at runtime, see below |

**Deliberate asymmetry, stated up front.** STAR bakes the annotation into its
index; HISAT2's equivalent is a `_tran` index, which JHU documents at ~200 GB
and would be a second enormous build. Instead HISAT2 gets the same splice sites
at runtime via `--known-splicesite-infile`, which is functionally equivalent for
alignment. The index-build comparison (D5) is therefore
**STAR-with-annotation vs HISAT2-without**, and must be reported that way.

## Instance and cost

JHU documents ~160 GB for a human SNP-aware build, so the r7i.4xlarge from the
original plan (128 GB) **cannot do it**.

| type | vCPU | RAM | us-east-1 on-demand |
|---|---|---|---|
| r7i.4xlarge | 16 | 128 GiB | $1.058 |
| **r7i.8xlarge** | 32 | **256 GiB** | **$2.117** |
| r7i.12xlarge | 48 | 384 GiB | $3.175 |

**Duration is the real unknown.** The only datapoint is mouse chr19: 61 Mb /
287k SNVs -> 25 s / 12.4 GB. Human is 3.1 Gb / ~15M variants, and the builder is
sublinear in memory but not in time. Several hours is likely; it cannot be
bounded from here.

### Staged, so the unknown is priced before it is bought

| stage | instance | est. | cost |
|---|---|---|---|
| **2a — chr1 probe** | r7i.8xlarge | 30–60 min | **~$1.50** |
| 2b — whole-genome SNP build | r7i.8xlarge | 3–8 h | **$6.50–17** |
| 2c — other indexes (linear, STAR, rustar) | same box | ~1.5 h | $3.20 |
| 2d — the measurements (D1–D4) | same box | ~1.5 h | $3.20 |
| EBS 300 GB gp3 + S3 | | | ~$0.60 |
| | | **total** | **$15–26** |

**Stage 2a is the decision point.** chr1 is 249 Mb, ~8% of the genome, with real
haplotypes — enough to measure peak RSS, wall time and variant retention, and to
say whether the whole-genome build is a 4-hour or a 12-hour proposition before
committing. If it extrapolates past ~10 h, better options are a bigger box
(r7i.12xlarge finishes sooner and may cost less overall) or thinning the variant
set.

Disk: STAR index ~30 GB, rustar ~25 GB, HISAT2 SNP ~8 GB + build scratch,
genome+GTF+variants ~6 GB, FASTQs 5 GB, outputs ~10 GB -> **300 GB gp3**.

## What this buys beyond the head-to-head

The head-to-head (D1–D4) is the point, but the build itself yields two numbers
nobody has:

1. **Variant retention on a real human build.** The guard added to `hgfm.h` this
   session prints `Local indexes: N of M variant instances retained (X%)`. That
   has never been measured on human. We have 98.9% on a 19.2 Mb targeted
   reference and catastrophic loss on over-dense mouse chr19 — the human figure
   is unknown and goes directly to how variant-aware the index really is.
2. **Whether a denser variant set loses more.** `snp151` > `snp144`. If
   retention drops materially against the shipped index's density, that says the
   shipped index sits near a density ceiling — which bears on the pseudogene
   filter question too, since that filter also changes variant density.

## Risks

- **An 8-hour build that hangs is expensive.** Dead-man timer is mandatory, and
  it must be re-armed after any stop/start — user-data only runs at first boot,
  which this session already tripped over.
- **The build may OOM at 256 GB.** If ~160 GB is optimistic for `snp151`'s
  density, stage 2a will show it in the extrapolation rather than 6 hours in.
- **Nothing about this reproduces the shipped index**, so prior Mac numbers stay
  non-comparable. They were already going to be retired in favour of AWS numbers
  for all three tools; this makes that mandatory rather than tidy.
