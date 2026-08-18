# E3 — a whole-genome HISAT2 graph index over the real 1000 Genomes panel

**Result: it exists. 23 generations, 5,917,131,871 path nodes, 671 GB peak RSS,
4 h 24 m on 96 vCPU, an 11 GB index.** The 32-bit format cannot represent it —
the node count is 138% of the 2^32 ceiling — and no such index has been
distributed.

The distributed `grch38_snp` is not a smaller version of this. Its haplotypes are
*invented*: `hisat2_extract_snps_haplotypes_VCF.py` partitions variants by greedy
graph colouring, producing ~0.78 haplotype rows per variant, where real phasing
carries ~1.11. That is the whole difference between an index that fits in 32 bits
at 77% of the ceiling and one that does not fit at all.

Inputs: GRCh38 primary assembly (3,099,750,718 bp, 194 contigs), 14,947,745 SNPs
and 16,338,502 haplotypes from the phased 1000 Genomes panel.

## The curve

`generation` g means "sorted by paths of length 2^g"; construction runs
`while(!isSorted())`, i.e. until `ranks == nodes`. Full curve in
`generations.txt`.

| gen | temp nodes | nodes | ranks | unsorted | closed |
|---|---|---|---|---|---|
| 9 | 3,635,922,500 | 3,595,956,954 | 3,457,621,509 | 138,335,445 | |
| **10** | **5,001,382,296** | 4,806,521,075 | 4,430,810,128 | 375,710,947 | |
| 11 | 6,181,237,301 | 5,920,332,832 | 5,902,565,765 | 17,767,067 | |
| 12 | 5,937,472,450 | 5,921,881,630 | 5,905,447,680 | 16,434,150 | 1,332,917 |
| 13 | 5,940,151,418 | 5,922,573,534 | 5,906,798,477 | 15,775,057 | 659,093 |
| 14 | 5,941,358,942 | 5,922,875,186 | 5,907,573,627 | 15,301,559 | 473,498 |
| 15 | 5,941,714,274 | 5,922,955,913 | 5,908,169,023 | 14,786,890 | 514,669 |
| 16 | 5,941,630,236 | 5,922,967,653 | 5,908,749,491 | 14,218,162 | 568,728 |
| 17 | 5,941,323,780 | 5,922,935,862 | 5,909,352,780 | 13,583,082 | 635,080 |
| 18 | 5,940,568,779 | 5,922,854,509 | 5,910,523,521 | 12,330,988 | 1,252,094 |
| 19 | 5,938,460,934 | 5,922,056,975 | 5,911,639,470 | 10,417,505 | 1,913,483 |
| 20 | 5,934,574,065 | 5,920,772,743 | 5,913,736,622 | 7,036,121 | 3,381,384 |
| 21 | 5,927,763,068 | 5,917,713,295 | 5,916,946,835 | 766,460 | 6,269,661 |
| **22** | 5,918,466,931 | **5,917,131,871** | **5,917,131,871** | **0** | 766,460 |

**Generation 10 is where the 32-bit build died**, and now there is a number for
it rather than an inference: 5,001,382,296 temporary nodes, **116.4% of
4,294,967,295**. The peak is generation 11 at 6,181,237,301 — **143.9% of the
ceiling**.

### The tail decelerates, then collapses

Between generations 12 and 14 the closure rate fell by half each time — 1.33M,
659k, 473k — which extrapolates to never finishing. That extrapolation was
wrong, and predictably so: prefix doubling resolves repeats in *bursts*, as 2^g
crosses the length of a repeat class. Generation 14 is 16 kb, below LINEs and
segmental duplications, so nothing was being resolved yet.

From generation 15 the rate rose every single generation — 515k, 569k, 635k,
1.25M, 1.91M, 3.38M, 6.27M — and generation 22 finished the remaining 766,460 in
one step. The same shape appears on a 900 kb graph, at 1/6000 the scale
(20.3k → 9.6k → 2.0k → 0), which is what made the prediction of ~20–22
generations safe while the rate was still falling.

**A decelerating tail in prefix doubling is not evidence of non-convergence.**
Reading it that way would have argued for killing a run that was 6 generations
from done.

## Memory, and why the instance size was the whole decision

**Peak RSS 655,593,232 kB = 625.2 GiB = 671.3 GB**, on a 743 GiB machine — 84%.

`lateGeneration` holds three live arrays at once (`past_nodes`, `from_table`,
`nodes`), so peak tracks the peak node count rather than the final one. At
generation 11's 5.92e9 nodes that is 3 × 5.92e9 × 32 B ≈ 568 GB of payload,
which the measured 671 GB brackets with allocator and thread overhead.

The earlier plan was to fall back to an r7i.16xlarge (512 GiB) because the
X-instance quota was 0. **That would have died**: 625 GiB does not fit in 512
GiB, and it would have died late, after hours of paid compute. Requesting the
Standard vCPU increase to 192 and waiting for an r7i.24xlarge was the correct
call, and the instruction that produced it — "don't waste money on something
that will probably OOM" — is now confirmed by measurement rather than by
judgement.

## The artifact

`s3://rustar-bench/hisat2-wg64/index/`, 11 GB, 64-bit (`.ht2l`):

| file | size |
|---|---|
| `genome.1.ht2l` | 3.4 GiB |
| `genome.2.ht2l` | 2.8 GiB |
| `genome.3.ht2l` | 21.0 KiB |
| `genome.4.ht2l` | 703.0 MiB |
| `genome.5.ht2l` | 2.0 GiB |
| `genome.6.ht2l` | 756.3 MiB |
| `genome.7.ht2l` | 936.1 MiB |
| `genome.8.ht2l` | 228.3 MiB |

## Cost and hygiene

4 h 27 m of r7i.24xlarge at $6.35/h ≈ **$28**. Every exit path wrote to the
watched S3 key and then powered off; the instance was launched with
`--instance-initiated-shutdown-behavior terminate` and an 11 h `shutdown -h +660`
dead-man switch, which was never needed — the build finished with 6.5 h to
spare. The instance terminated itself and no instance is running in the region.

## What this is for

The three live arrays are the entire argument for the Rust builder in
`../rust/README.md`: this is a property of the doubling loop, not of the
alignment or the format, and it is what an external-memory implementation would
change. The Rust prototype now reproduces that loop generation for generation on
graphs up to 900 kb, so the loop whose memory profile this measures is the loop
it matches.
