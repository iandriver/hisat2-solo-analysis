# Rust builder prototype

Staged against the validation ladder in `../rust_builder_plan.md`. Each rung must
pass before the next is attempted.

## S1 — format round-trip (rung 1): **passing**

`ht2fmt` reads a `.ht2` header and re-emits it byte-identically. No construction
code; the only claim under test is that the on-disk layout is understood well
enough to reproduce it. Field order is taken from `GFM::readIntoMemory`
(`gfm.h:5905`), the authoritative reader.

```
cargo run --release -- /path/to/genome_snp.1.ht2
```

On the distributed `grch38_snp`:

```
  index version      2.0.2-beta      <- hisat2-inspect -s agrees
  ftabChars          10              <- agrees
  offRate            4               <- agrees (SA-Sample 1 in 16)
  len                2,945,849,067
  gbwtLen            3,315,031,465
  numNodes           3,300,850,145
  ROUND-TRIP OK: 44 header bytes reproduced exactly
```

### An observation that closes a loop

`numNodes = 3,300,850,145` is **76.9% of the 2^32 ceiling**. The E0 run against
the real 1000 Genomes phased panel is heading to ~4.5e9, about 105%.

The difference is the haplotype effect measured in stage 2b: real phasing
carries ~1.11 haplotype rows per variant, the greedy graph-colouring partition
in `hisat2_extract_snps_haplotypes_VCF.py` produces ~0.78. **The distributed
index exists because its haplotypes are invented and simpler than real ones** —
and even so it sits at 77% of the ceiling.

Anyone rebuilding a human graph index from real phased data is therefore not
doing a slightly harder version of what JHU did; they are doing something the
32-bit format cannot represent.

## Next rungs

2. `22_20-21M.fa`, no variants — full index byte-identical to `hisat2-build`.
3. Same with SNPs and haplotypes.
4. chr1 against the three `analysis/hap` variant sets, matching measured retention.
5. Alignment equivalence: SAM byte-identical at `-p 1 --seed 0 --reorder`.
6. Whole genome, where no C++ reference output exists.
