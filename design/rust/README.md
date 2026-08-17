# Rust builder prototype

Staged against the validation ladder in `../rust_builder_plan.md`. Each rung must
pass before the next is attempted.

## S1 — format layout (rung 1): **passing**

`ht2fmt` parses a `.1.ht2` in full and accounts for every byte. No construction
code; the claim under test is only that the on-disk layout is understood
exactly, which is the prerequisite for emitting anything.

**The test has teeth because every section size is derived, not stored.** The
bulk `gbwt` block, `ftab`, `offsLen` and the side geometry all come from
`GFMParams::init` (`gfm.h:138`). Reproduce any formula wrongly and the sections
will not sum to the file length.

```
cargo run --release -- /path/to/genome_snp.1.ht2
```

Verified against three real indexes spanning both modes and three orders of
magnitude:

| index | mode | file bytes | sections sum |
|---|---|---|---|
| `22_20-21M_snp` | graph FM | 4,789,374 | exact |
| `grch38_snp` | graph FM | 2,044,241,775 | exact |
| `grch38` | linear FM | 986,172,031 | exact |

Both branches of `GFMParams::init` are exercised: graph indexes reserve 6
`index_t` per side and pack 2 nucleotides per byte, linear ones reserve 4 and
pack 4. The mode is inferred from `len + 1 == gbwtLen`, not stored.

Header fields also cross-check against an independent tool — `hisat2-inspect -s`
agrees on index version (2.0.2-beta), `ftabChars` (10) and `offRate` (4, shown
there as "SA-Sample 1 in 16").

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
