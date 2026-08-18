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

## S2 — construction (rung 2): front end passing, SA semantics identified

`ht2build` builds the front end from a FASTA and checks each piece against a
reference index from `hisat2-build` on `example/reference/22_20-21M.fa`
(1,000,000 bp with a 100,000 bp N run, so `len` = 900,000, `plen[0]` = 1,000,000).

### Verified exact

| quantity | how it was checked |
|---|---|
| `len`, `nPat`, `plen[0]`, `nFrag` | compared to the stored header |
| `rstarts[]`, both fragments, all three fields | compared to the stored records |
| `fchr[5]` | compared to the stored histogram |
| **BWT unpacking** | **4688/4688 sides** — recounted our unpacked rows against every side's stored occ tally |
| forward orientation | reversed text scores 25.7%, i.e. chance for a 4-letter alphabet |

### Layout, settled from `GFM::buildToDisk` (`gfm.h:5148`)

Reading the writer rather than inferring from the reader resolved the packing:

- rows are laid out sequentially, **4 per byte, low bit-pair first** —
  `pack_2b_in_8b(c, byte, bpi)` is `byte |= c << (bpi*2)`, and `fw` is `true`
  throughout the linear path, so the "backward bucket" branch never runs;
- each side holds `sideGbwtSz*4` rows, with the four occ tallies in the final
  `4*sizeof(index_t)` bytes;
- **those tallies are `occSave`, the counts as of the *start* of that side**,
  not the running total;
- padding rows past the end of the SA are written as `A` and *are* counted;
  the `zOffs` row is written as `A` and is *not* counted;
- `zOffs` is the row where `saElt == 0`, confirmed twice: `zOffs.push_back(si)`
  in the writer and `if(elt == _zOffs[i]) return 0;` at `gfm.h:2722`.

### The open question, now precise

`.2.ht2` stores HISAT2's own SA values (one per `2^offRate` rows), which is
direct ground truth. Against ours:

```
their row k  vs  our row k+0:      0/56251 =  0.00%
their row k  vs  our row k+1:  36938/56250 = 65.67%
their row k  vs  our row k+2:  19192/56250 = 34.12%
```

The two offsets sum to ~99.8% and the changeover is sharp: **offset +1 holds up
to their row 194,000, offset +2 after it.** So their array is ours minus the
leading empty suffix, and then ours contains **one extra element around our row
194,002** (suffix 877,124) that theirs does not.

Our suffix array is itself sound — an independent brute-force count found
815,267 suffixes smaller than the whole string, matching it — so this is not a
sorting bug.

### What `nextSuffix()` emits, from `blockwise_sa.h`

`KarkkainenBlockwiseSA::nextBlock` (`blockwise_sa.h:924`) fills each bucket with
positions strictly less than `len` — `bucket.resize(len); for(i < len)` in the
all-inclusive branch, and `for(... i < len ...) bucket.push_back(i)` with
`assert_lt(i, len)` in the bucketed one. The empty suffix is never a bucket
member. Then, after sorting:

```cpp
if(hi != OFF_MASK) bucket.push_back(hi);   // non-final: append the RHS pivot
else               bucket.push_back(len);  // final: append the $ suffix
```

So the emitted sequence is: sorted bucket, its pivot, sorted bucket, its
pivot, ..., final sorted bucket, then **`len` — the empty suffix — last**.
Appending the pivot after sorting is order-preserving because the pivot bounds
its bucket above. Total is `len` real positions plus one, matching
`size() == text.length()+1`.

**HISAT2's suffix array is therefore the ordinary one with `$` at the END, not
the beginning.** That is consistent with `.2.ht2` row 0 holding a real suffix
(176,766) and with the BWT's best shift being −1.

### SOLVED — HISAT2 pads past the end with a character *larger* than all others

Recovering HISAT2's **full** suffix array settled it. `.2.ht2` samples only one
row in 16, too coarse to localise the difference, so `ht2sa` walks the stored
BWT instead: `$` is the last row (`len`), and LF-stepping from there visits the
row of every suffix in turn, labelling all 900,001 rows. Two checks confirm the
walk: the reconstructed text matches the reference exactly, and after 900,000
steps it lands on the stored `zOffs` (815,266).

Against that full array, ours differed in 309,050 of 900,000 rows — but the two
held the **same multiset** of positions, and HISAT2's array had only **11
adjacent pairs out of lexicographic order**. A handful of far-displaced elements
shifts a long run of rows by one, which is exactly the "+1 then +2" pattern seen
earlier through the sampled view.

All 11 violations have the same shape:

| their row | longer suffix | shorter suffix | len | shorter is a prefix of longer |
|---|---|---|---|---|
| 17,670 | 669,678 | 899,989 | 11 | yes |
| 48,277 | 154,949 | 899,990 | 10 | yes |
| 189,243 | 154,950 | 899,991 | 9 | yes |
| … | … | … | … | yes |
| 899,998 | 64,679 | 899,999 | 1 | yes |

The suffixes involved are precisely the 11 shortest (lengths 11 down to 1), and
in every case HISAT2 places the **longer** suffix first where true
lexicographic order puts the shorter (prefix) first.

**The rule: HISAT2 sorts suffixes as though the text were padded past its end
with a character greater than any real one, rather than terminated by a
sentinel smaller than all.** That single rule accounts for both observations —
the empty suffix (length 0) sorting last, and every proper prefix sorting after
its extension.

Rebuilding our suffix array under that rule reproduces HISAT2's array **exactly:
0 mismatches of 900,001 rows** (`sa_rule_check.rs`).

Note this is a real deviation from a textbook suffix array, affecting the 11
shortest suffixes here. It is almost certainly harmless for alignment — those
are the shortest possible matches — but **byte-identical construction requires
reproducing it**, which is the sort of thing only a byte-level target would
have surfaced.

### Rung 2 status

Solved: section geometry, the front end (`nPat`, `plen`, `nFrag`, `rstarts`,
`fchr`), the BWT row layout and packing, `zOffs`, and now the suffix-array
order. Still to do before a byte-identical build: `ftab`, `eftab`, the `.2.ht2`
SA sample, `refnames`, and emitting the file.

## Next rungs

2. `22_20-21M.fa`, no variants — full index byte-identical to `hisat2-build`.
3. Same with SNPs and haplotypes.
4. chr1 against the three `analysis/hap` variant sets, matching measured retention.
5. Alignment equivalence: SAM byte-identical at `-p 1 --seed 0 --reorder`.
6. Whole genome, where no C++ reference output exists.
