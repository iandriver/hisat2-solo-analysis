# Rung 2 full-build benchmark — Rust emitter against `hisat2-build`

**Byte-identical at every size, up to 50 Mb of real human sequence.** All four
files (`.1`–`.4.ht2`) match `hisat2-build`'s output exactly, so this compares two
builders producing *the same bytes* rather than two different indexes — which is
the only way a build benchmark means anything.

References are prefixes of GRCh38 chr20 taken from the linear index's own dump,
with ambiguous characters removed so each size label is exact. Single-threaded
both sides, M5 Pro.

| size | builder | wall (s) | peak RSS (MB) | bytes/bp | identical |
|---|---|---|---|---|---|
| 1 Mb | `hisat2-build` | 0.27 | 96 | 96 | |
| | **`ht2emit`** | 0.50 | **38** | 38 | **4/4** |
| 5 Mb | `hisat2-build` | 1.12 | 110 | 22 | |
| | **`ht2emit`** | 5.58 | 94 | 19 | **4/4** |
| 20 Mb | `hisat2-build` | 4.63 | **165** | 8.3 | |
| | `ht2emit` | 39.19 | 330 | 16.5 | **4/4** |
| 50 Mb | `hisat2-build` | 15.91 | **229** | 4.6 | |
| | `ht2emit` | 169.89 | 766 | 15.3 | **4/4** |

## The memory curves cross, and that is the finding

At 1 Mb the Rust builder uses **2.5x less** memory than `hisat2-build`. At 50 Mb
it uses **3.3x more**. The crossover is around 5 Mb.

The reason is visible in the bytes-per-bp column. `hisat2-build` falls from 96 to
**4.6 bytes/bp** as the reference grows — its blockwise Kärkkäinen construction
with a difference cover holds a bounded working set and amortises its fixed
overhead. `ht2emit` sits flat at **15–16 bytes/bp**, because prefix doubling
keeps `sa`, `rank` and `tmp` as full `u32` arrays over the whole text, plus the
text itself: 12 bytes/bp of unavoidable state, and no way to bound it.

Time tells the same story: the ratio grows 1.9x -> 5.0x -> 8.5x -> **10.7x**.
`hisat2-build` scales close to n log n; the naive doubling is O(n log^2 n) with a
full sort every round.

## What that means for the plan

**The suffix array has to be replaced before the Rust builder is useful for
anything but verification.** Extrapolating 15.3 bytes/bp to a 3.1 Gb genome is
~47 GB for the *linear* index alone, against `hisat2-build`'s measured 4.6
bytes/bp, and `u32` runs out at 4.3 Gb regardless. SA-IS or divsufsort is linear
time in ~5n bytes and is a drop-in for this stage; the blockwise approach
`hisat2-build` already uses is the other option and has the better memory
profile.

This does not touch the reason the Rust builder exists. The 625 GiB measured on
the whole-genome graph build comes from `PathGraph`'s three live arrays, not from
the suffix array — different code, different rung. But it does say the naive SA
cannot be carried into that work, and it puts a number on it rather than a
suspicion.

The honest summary: **`hisat2-build` is better on both axes at scale today.**
What the Rust side has is byte-identity — a verification harness that can prove
any future construction, external-memory or otherwise, produces exactly the
index HISAT2 produces. That was the goal of the rung, and the benchmark's job was
to say what it cost.

## Reproducing

```
analysis/rung2bench/bench.sh 1M 5M 20M 50M
```
