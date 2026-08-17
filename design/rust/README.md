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

## S2 — construction (rung 2): **front end passing, BWT layout not yet solved**

`ht2build` builds the front end from a FASTA and checks each piece against a
reference index produced by `hisat2-build` on `example/reference/22_20-21M.fa`
(1,000,000 bp containing a 100,000 bp N run, so `len` = 900,000 while
`plen[0]` = 1,000,000).

**Verified exact:**

| quantity | status |
|---|---|
| `len`, `nPat`, `plen[0]`, `nFrag` | match |
| `rstarts[]` — both fragments, all three fields | match |
| `fchr[5]` (A/C/G/T cumulative histogram) | match |
| suffix array | independently brute-force verified |

The suffix array was checked against a direct count of suffixes lexicographically
smaller than the whole string (815,267), which agrees with the constructed SA.
The indexed text was also confirmed byte-identical to what `hisat2-inspect`
reconstructs from the index, so the input to the BWT is not in question.

**Not yet solved:**

- `zOffs` — we compute 815,268 under sentinel-first numbering, 815,267 under
  sentinel-last. The stored value is **815,266**. `gfm.h:2722`
  (`if(elt == _zOffs[i]) return 0;`) confirms `zOffs` is the row whose text
  offset is 0, which is what we compute, so the gap is not a definitional one.
- **The row-to-BWT mapping.** Unpacking the stored `gbwt` with the layout from
  `GFM::postReadInit` (`gfm.h:2783`) — 4 characters per byte, low bits first,
  in the first `sideGbwtSz` bytes of each `sideSz` side — matches our BWT for
  the first ~20,000 rows but **disagrees on 20.5% of all 899,999 rows**.

  The disagreements are not at the fragment boundary and are not a constant
  offset; they appear as short runs of consecutive rows whose values look
  permuted. That pattern suggests the row mapping drifts rather than being
  uniformly shifted.

  **Recorded as a caution:** an earlier pass sampled only the first 20,000 rows,
  found 99.99% agreement, and concluded the layout was solved. It was not. The
  full-range check is the one that counts.

Since the SA, the text and the character histogram are all independently
confirmed, the remaining discrepancy is in how rows are laid out in the packed
sides, not in the index content. That is the next thing to work out.

## Next rungs

2. `22_20-21M.fa`, no variants — full index byte-identical to `hisat2-build`.
3. Same with SNPs and haplotypes.
4. chr1 against the three `analysis/hap` variant sets, matching measured retention.
5. Alignment equivalence: SAM byte-identical at `-p 1 --seed 0 --reorder`.
6. Whole genome, where no C++ reference output exists.
