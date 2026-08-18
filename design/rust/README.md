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

### Where the two arrays actually diverge

With `$`-last understood, the mapping is `their row r == our row r+1`, and it
holds **exactly** up to their row 194,000:

```
k=    0  theirRow=      0  theirs=176766  our row      1 =176766  ok
k=    1  theirRow=     16  theirs=759847  our row     17 =759847  ok
k=12125  theirRow= 194000  theirs= 13331  our row 194001 = 13331  ok
k=12126  theirRow= 194016  theirs=304543  our row 194017 =723517  DIFF
```

Two sortedness tests were run, and they are **not** equivalent — worth stating
plainly because it changes which array is suspect:

| array | test | result |
|---|---|---|
| ours | **all 900,000** adjacent pairs compared as suffixes | 0 out of order |
| theirs | only **sampled** pairs, 16 rows apart (`.2.ht2` stores 1 in 16) | 0 out of order |

Ours is therefore a verified suffix array. Theirs passes a test that **cannot
see a permutation inside a 16-row window**, so it is much weaker evidence. The
earlier reading — that ours must be wrong because a suffix array is unique —
does not follow from it.

### The part that still does not add up

That structure predicts a *uniform* offset of +1: their row k should equal our
row k+1 for every k below `len`. The `.2.ht2` sample says otherwise — +1 holds
for 65.7% of sampled rows and +2 for 34.1%, with a sharp changeover at their row
194,000. A uniform relabelling cannot produce that.

So the emission order is understood and is not sufficient on its own. The
divergence starts abruptly at one row rather than drifting, which is the
signature of a **bucket boundary**: buckets are sorted independently and their
pivots are appended after sorting. If a pivot is not in fact the maximum of its
bucket range — for instance because bucket membership is decided by a
bounded-depth comparison — the emitted order would be locally wrong there and
everything after it would shift.

Concrete next step: recover the `_sampleSuffs` pivot positions for this build
and check whether 194,001 is one of them. `KarkkainenBlockwiseSA::qsort`
(`blockwise_sa.h:436`) and the difference-cover tie-breaking are where a
bounded-depth comparison would live.

Until that is settled, no BWT we generate can match: the best full-range
agreement is 79.5% at shift −1, which is exactly what a single insertion
partway through would produce.

### A method note worth keeping

An earlier pass gated on the first 20,000 rows, saw 99.99% agreement, and
recorded the layout as solved. The full-range check gave 79.5%. Prefix samples
are worthless here precisely because the divergence is a single insertion two
thirds of the way in — the first 20,000 rows agree under *any* hypothesis that
gets the early offset right.

## Next rungs

2. `22_20-21M.fa`, no variants — full index byte-identical to `hisat2-build`.
3. Same with SNPs and haplotypes.
4. chr1 against the three `analysis/hap` variant sets, matching measured retention.
5. Alignment equivalence: SAM byte-identical at `-p 1 --seed 0 --reorder`.
6. Whole genome, where no C++ reference output exists.
