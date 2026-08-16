# Variant loss concentrates in the MHC — the region the project is about

`hisat2_index_probe.py` run on the whole-genome 1000 Genomes phased set
(14,947,745 variants, 16,338,502 haplotypes, GENCODE v50 primary assembly).
Seven seconds, no build.

## It calls the failure

```
--- global node bound (32-bit) ---
  reference alone      3,099,750,718 nodes  (72.2% of the ceiling)
  plus variants        3,115,407,639 nodes  (72.5% of the ceiling)
  verdict              WILL FAIL -- GRCh38 failed twice at this fraction (72%)
```

That is what stage 2b established over four builds and roughly $24: 32-bit died
with `exceeded integer bounds`, 64-bit was OOM-killed above 368.8 GB. The probe
reaches the same verdict before an instance is launched.

## Where the dropped variants would be

| | windows | variants | over budget | variants at risk |
|---|---|---|---|---|
| genome-wide | 48,705 | 14,947,745 | 9,309 (19.1%) | 28.9% |
| extended MHC (chr6:28.5-33.5 Mb) | 89 | 66,177 | 49 (55.1%) | **83.9%** |

The MHC carries **2.4x** the genome-wide variant density, and **84% of its
variants sit in windows the builder cannot hold** against 29% genome-wide. Of
the 100 densest over-budget windows in the genome, **27 are in the MHC** — a
region that is 89 of 48,705 windows, or 0.18% of them.

Every HLA gene measured in the reference-bias work sits in an over-budget
window:

| gene | variants in its window(s) | over budget | measured T1 gain |
|---|---|---|---|
| HLA-DQA1 | 3,145 | 1/1 | 32.10x |
| HLA-DQB1 | 5,922 | 2/2 | 3.00x |
| HLA-DRB1 | 2,717 | 1/1 | 1.54x |
| HLA-C | 2,192 | 1/1 | 1.51x |
| HLA-A | 2,314 | 1/1 | 1.00x |
| HLA-B | 2,576 | 1/1 | 1.00x |

The densest window in the genome, `chr6:32,609,280-32,666,624` at 3,145 variants
(8.2x the effective capacity of 384), covers HLA-DQA1 — the gene with the
largest measured advantage of the graph index over a linear one.

## What this does and does not mean

**It does mean the constraint and the payoff share an address.** The windows the
builder is least able to represent are the windows where variant-aware alignment
was worth the most. Whatever variant set is chosen, the MHC is where the losses
land, so a build that reports 87% retention genome-wide is losing far more than
13% of exactly the variants that motivate the feature.

**It does not mean the measured gains were suppressed by variant loss.** Being
in an over-budget window does not predict the benefit: all six genes above are
over budget, and their gains range from 1.00x to 32.1x. The dose-response in
`../bias` is the direct test, and it found retention between 14.6% and 87.5%
moves the ALT fraction by 0.023 points — so recovering the dropped MHC variants
would probably not move these numbers much either.

The right reading is narrower and still worth having: **the T1 results were
obtained with a partial variant set in precisely those windows**, and they were
large anyway. That makes them a floor rather than a ceiling, without any claim
about how much headroom is left.

**It sharpens what to do about density.** Thinning variants uniformly to fit the
node bound would take them disproportionately from the MHC, because that is
where the density is. If a future build has to lose variants, losing them
*outside* chr6:28.5-33.5 Mb costs less of what the index is for. Nothing in
HISAT2 offers that control today — the backoff is per-window and blind to what
the window contains.

## Caveats

- Retention percentages are calibrated on chr1 with the stock halving backoff.
  The binary-search backoff now in `upstream/modernize` retains more, so the
  87.6% genome-wide figure is a floor for that build.
- Moot for this variant set regardless: it cannot complete a 32-bit build at all.
- Window assignment uses `pos // local_index_interval`, which matches how the
  builder partitions, but the 1,024 bp overlap between consecutive local indexes
  means a variant near a boundary is present in two windows. That inflates
  per-window counts slightly at boundaries and is not corrected here.

## Reproducing

```
hisat2_index_probe.py --fai genome.fa.fai --snp genome.snp \
    --haplotype genome.haplotype -v
```

`mhc_concentration.py` recomputes the tables above from the same `.snp` file.
