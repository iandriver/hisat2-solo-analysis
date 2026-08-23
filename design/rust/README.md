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

### ftab and eftab — exact

`ht2ftab` reproduces both tables: **1,048,577/1,048,577 ftab entries and 20/20
eftab entries**.

Walking the SA in row order, a suffix with `len - saElt >= ftabChars`
contributes `ftab[sufInt+1]++`, packing the first `ftabChars` characters two
bits each, leftmost most significant. A shorter suffix bumps `absorbCnt`, which
lands in `absorbFtab[sufInt]` at the next long suffix — or in
`absorbFtab[ftabLen-1]` if the walk ends with `absorbCnt > 0`. Finalisation then
turns `ftab` into a running prefix sum, with absorbing entries storing `[lo,hi]`
in `eftab` and keeping the marker `eftabCur ^ INDEX_MAX` in `ftab`, resolved by
`ftabLo`/`ftabHi` (`gfm.h:2617`).

One subtlety cost a wrong `eftab[17]` (900,000 against 900,001) before being
found: **the empty suffix is not skipped.** `buildToDisk` sees
`len - saElt == 0`, which is below `ftabChars`, so the `$` row counts as a short
suffix and bumps `absorbCnt` like any other. That single row is what lifts the
top of the table from `len` to `len + 1`, matching the writer's own
`assert_eq(ftabHi(..., ftabLen-1), len+1)`.

### Rung 2 — **PASSING: byte-identical**

`ht2emit` writes the whole index from a FASTA and byte-compares it against
`hisat2-build`'s own output. This is the only test on the rung that cannot be
passed by accident: a misunderstanding anywhere shows up as a differing offset
rather than as a number that happens to agree.

```
ht2emit <reference.fa> <out_prefix> <reference.1.ht2>
```

| reference | `.1.ht2` | `.2.ht2` | `.3.ht2` | `.4.ht2` |
|---|---|---|---|---|
| `example/reference/22_20-21M.fa` (1 seq, 900,000 bp) | 4,494,550 B exact | 225,008 B exact | 26 B exact | 225,000 B exact |
| `testdata/multi.fa` (5 seqs, 15,700 bp) | 4,199,874 B exact | 3,932 B exact | 89 B exact | 3,925 B exact |

`.2.ht2` is trivial once the suffix array is right — a `1` sentinel followed by
`saElt` for every `2^offRate`-th row — and it agreeing is the independent
confirmation that the SA rule from the previous section is exactly right, at
every sampled row rather than in aggregate.

Two header fields are written twice and only the second value survives:
`gbwtLen`/`numNodes` (offsets 12/16, `gfm.h:5166`) and `eftabLen` (offset 36,
`gfm.h:5471`). All three go out as `0` in `writeFromMemory(justHeader=true)` and
are seeked back over by `buildToDisk`. `eftabLen` is also always `ftabChars*2`
regardless of how many entries actually absorb — the code computes the true
count and then discards it (`gfm.h:5439`).

### The second reference is what made the test worth running

The single-sequence example passed on the first attempt. `multi.fa` — five
sequences with a header description, leading Ns, a trailing N run, lowercase
`n`, and a wholly lowercase sequence — failed on two counts, both of which the
example is structurally unable to detect:

1. **`refnames` stores the whole header line.** `_refnames_nospace`
   (`gfm.h:1400`) is a *separate* whitespace-truncated copy, used only for
   matching chromosome names against SNP and splice-site files. A header with no
   description cannot tell the two apart.
2. **`.3.ht2` holds more records than `.1.ht2` has `rstarts` entries.** A
   sequence that ends in ambiguity contributes a record with `len == 0` carrying
   the trailing count — `seqC` and `seqD` each add one, 9 records against 7
   fragments. `joinToDisk` keeps only the non-empty records for `rstarts`, so
   the two counts are genuinely different quantities and an index whose
   reference happens not to end in N never reveals it.

Neither would have shown up as a wrong number; both showed up as a wrong byte at
a known offset.

### Scope

`.5`–`.8` are out of scope at this rung by design. `.5`/`.6` hold the
hierarchical local indexes (`hgfm.h:2068`), which only exist once there is a
graph to localise, and belong with the SNP/haplotype rung that introduces
`PathGraph`.

## Rung 3, step 1 — the reference graph, and where it stops matching

`ht2graph` builds `RefGraph` from a FASTA plus `.snp`/`.haplotype` files and
checks its size against HISAT2's own first log line. `PathGraph::makeFromRef`
(`gbwt_graph.h:1821`) sets `temp_nodes = base.edges.size() + 1`, and `printInfo`
reports `temp_nodes` as the left-hand number, so

```
Generation 0 (236 -> 236 nodes, 0 ranks)
```

is a free exact oracle for the edge count — available before any of the
prefix-doubling machinery exists.

```
ht2graph <reference.fa> <snp> <haplotype> [<build.log>]
```

The construction itself is straightforward, and reproduced exactly: node 0 is a
head labelled `Y`, reference position `p` is node `p + 1`, the tail `Z` is node
`len + 1`, and each haplotype **duplicates its whole `[left, right]` span**
rather than just the variant. Duplicating the span is what preserves phase — a
read can only walk an allele combination some haplotype actually carries.

| case | ours | HISAT2 | |
|---|---|---|---|
| 200 bp, 3 singles, 2 haplotypes | 236 | 236 | exact |
| 509,431 bp, 86 deletions | 509,519 | 509,519 | exact |
| 509,431 bp, 1,721 singles | 512,875 | 512,811 | **+64** |
| 509,431 bp, 74 insertions | 509,884 | 509,715 | **+169** |
| 900,000 bp, all 3,502 variants | 907,398 | 906,917 | **+481** |

*(the +N rows are the naive walk alone; with reverse-determinisation added, every
row below is exact — see the end of this section)*

### The gap is reverse-determinisation, and it is not an edge case

`RefGraph` finishes by testing `isReverseDeterministic` — *no node may have two
incoming edges from nodes carrying the same label* (`gbwt_graph.h:191`) — and
running a subset construction if it fails. Both the fragmented and the simple
branch do this, so it is not a scale artifact.

Reducing to one variant on a 200 bp reference isolates it exactly:

| single variant | ours | HISAT2 | delta |
|---|---|---|---|
| one SNP | 204 | 204 | 0 |
| one 3 bp deletion | 203 | 203 | 0 |
| one 1 bp insertion | 206 | 204 | **−2** |

**Every insertion costs exactly two edges, whatever base is inserted.** The
mechanism is in the haplotype loop: an insertion consumes no reference position,
so `for(j...) { if(prev_ALT_type == ALT_SNP_INS) j--; ... }` (`gbwt_graph.h:679`)
revisits the same `j` and emits a *duplicate* of the reference node the
insertion lands on. That duplicate carries the same label as the original and
points at the same successor, which is precisely the violation
`isReverseDeterministic` looks for. Collapsing it removes one node and two
edges: the insertion's exit edge and the duplicate's in-edge become the same
edge, as do the two exit edges.

So the redundancy is generated deliberately by the simple walk and cleaned up
afterwards by the determinisation, rather than being avoided. Reproducing
`RefGraph` means reproducing both halves — a builder that emits the tidier graph
directly would not match, and on the full example set the correction is **13.7%
of all variant edges**, not a rounding error.

Deletions needing no correction is the useful control: they add an edge but no
node, so they cannot create a duplicate-label predecessor.

Singles contribute 64 of the 481 (3.7% of 1,721), which the one-SNP case shows
is not intrinsic to a SNP — it comes from variants close enough to interact.

### With the subset construction added: exact everywhere

`reverse_determinize` implements `gbwt_graph.h:1015` — walk backward from the
tail `Z`, and at each composite node collect the predecessors of every member,
group them by label, and share any group whose member set already exists. That
sharing is what collapses the duplicates.

| case | ours | HISAT2 |
|---|---|---|
| 200 bp, 3 singles | 236 | 236 |
| 509,431 bp, 1,721 singles | 512,811 | 512,811 |
| 509,431 bp, 74 insertions | 509,715 | 509,715 |
| 509,431 bp, 86 deletions | 509,519 | 509,519 |
| 509,431 bp, all 1,881 | 513,179 | 513,179 |
| **900,000 bp, all 3,502, two fragments** | **906,917** | **906,917** |

The last row is the one that carries weight: at 900,000 bp HISAT2 takes its
*fragmented automaton* branch (`jlen >= 1 << 16`), building per-fragment
automata and merging them, which is different code from the single-chain walk
reproduced here. The two agree exactly. The determinisation is evidently what
makes the branch invisible — whatever redundancy the fragmenting introduces is
removed by the same subset construction.

The tool also asserts the post-condition the construction exists to establish:
no node with two same-label predecessors, checked on the output rather than
assumed.

## Rung 3, step 2 — the prefix doubling, matched generation for generation

`ht2path` runs the doubling over the graph `ht2graph` builds and compares the
**whole** curve, not just its first line. `generation` means "sorted by paths of
length 2^g", construction runs `while(!isSorted())`, and each generation prints
`temp_nodes -> nodes, ranks` — three numbers per generation, all of which depend
on the exact pruning rule. A wrong rule diverges within a generation or two and
never recovers, so this is not a test that can be passed approximately.

```
ht2path <reference.fa> <snp> <haplotype> <build.log>
```

| case | generations | result |
|---|---|---|
| 200 bp, 3 singles | 6 | all three numbers exact, every generation |
| 509,431 bp, 1,881 variants | 11 | exact |
| 900,000 bp, 3,502 variants, two fragments | 11 | exact |

99 numbers across 33 generation lines, none off by one.

```
  gen |          ours          |         HISAT2         |
    4 |    958371    939876 802762 |    958371    939876 802762 |
    5 |    948781    945155 882419 |    948781    945155 882419 |
   10 |    954009    953833 953833 |    954009    953833 953833 |
```

### The four-way split is arithmetic, not just memory management

HISAT2 divides the loop into `generationOne`, `earlyGeneration`,
`firstPruneGeneration` and `lateGeneration`, which reads like a memory
optimisation. It is not only that — the stages compute different things, and
reproducing the curve means reproducing each:

- **Generations 1–3 do no pruning at all**, so `nodes == temp_nodes` and
  `ranks == 0`. Keys are packed into a single integer, the left key shifted up
  by `3 * 2^(g-1)` bits. That is why the stage ends at 3: at generation 4 the
  packed key would need 24 bits per side and no longer fits, so the
  representation switches to a pair of ranks.
- **Generation 4** does the same join, then a full sort by key and the first
  pruning pass — and `mergeUpdateRank` has a *completely different body* for
  this generation, built on `nextMaximalSet`, collapsing each maximal run that
  shares a `from`.
- **Generations 5+** skip already-sorted nodes entirely and re-join only the
  unsorted ones, against a `from_table` built by sorting a copy by `from` while
  the main list stays in rank order. The join output therefore arrives already
  grouped by `key.first`, which is why there is no full sort in this stage — only
  the within-block sort by `key.second` that `mergeUpdateRank` does itself.

The pruning rule that matters is the one in the block walk: **a run of nodes
sharing a single `from` is one path node seen several ways and collapses to one;
a run spanning several `from` values is genuinely several nodes at the same
rank** and all of them survive. Getting that backwards changes `nodes` in the
first pruned generation and the curve never rejoins.

## Rung 3, step 3 — `generateEdges`, checked against the index header

With the nodes sorted, `generateEdges` (`gbwt_graph.h:2367`) turns them into GFM
rows: each reference edge contributes one row per path node whose `from` is that
edge's `to`, labelled by the character of the edge's `from` node, bucketed by
label and sorted by the ranking it points at. That ordering is the BWT.

This has its own oracle, and a better one than the log — the header of the index
HISAT2 actually wrote, which comes from a different part of the pipeline
entirely:

| case | `numNodes` | `gbwtLen` |
|---|---|---|
| 200 bp, 3 singles | 240 = 240 | 242 = 242 |
| 509,431 bp, 1,881 variants | 530,869 = 530,869 | 532,723 = 532,723 |
| 900,000 bp, 3,502 variants | 953,832 = 953,832 | 957,345 = 957,345 |

Better still, the per-label row counts match the stored `fchr`, which is a
**content** check rather than a total — it says how many rows carry each
character, so it fails if the rows are right in number but wrong in kind:

| case | A | C | G | T |
|---|---|---|---|---|
| 200 bp | 81 | 58 | 42 | 60 |
| 509,431 bp | 117,761 | 146,867 | 147,258 | 120,836 |
| 900,000 bp | 224,717 | 251,539 | 252,066 | 229,022 |

Exact on all three.

`numNodes` is **path nodes − 1** and `gbwtLen` is the path-edge count exactly.
The off-by-one is worth stating rather than absorbing: `makeFromRef` adds a final
self-looping node for the tail `Z`, and it survives the doubling but is not a
GFM row, which is also why the `Z` label bucket is always empty (`Y` holds
exactly one row, the head).

## Rung 3, step 4 — the GFM rows themselves, character for character

The row *order* is not the `(label, ranking)` order `generateEdges` leaves
behind. `nextRow` (`gbwt_graph.h:1609`) walks **path nodes in rank order and,
within each node, its incoming edges**; `F` marks each node's first row, and a
second independent cursor emits `M` by walking the same nodes and advancing by
out-degree.

Reaching that order needs the tail of `generateEdges` too: rewriting `edge.from`
from a reference-graph node id to a path-node index while accumulating each
node's out-degree (`:2561`), relabelling `Y` to `Z` and dropping the
second-to-last node (`:2576`), then re-sorting and building the CSR (`:2604`).

Checked by unpacking the BWT the index actually stores — graph sides pack 2 rows
per byte with the F and M bitvectors interleaved after the characters — and
comparing row by row:

| case | rows | characters | F bits |
|---|---|---|---|
| 200 bp, 3 singles | 242 | **242/242** | **242/242** |
| 509,431 bp, 1,881 variants | 532,723 | **532,723/532,723** | **532,723/532,723** |
| 900,000 bp, 3,502 variants | 957,345 | **957,345/957,345** | **957,345/957,345** |

### The tie-break the source misreads

`PathEdgeToCmp` compares `(to, from)`, so the obvious implementation sorts the
edges by `(ranking, from)`. That is wrong, and wrong in a way only a row-level
comparison finds: on the 200 bp graph it swapped exactly **two** rows out of 242
— 117 and 118, both inside a single node's edge group (F=1 then F=0).

The stored order there is A then C, i.e. **label order**. The final
`radix_sort_copy` is fed the edges in label-bucket order and is stable within a
ranking, so ties resolve by label and the `from` field never participates.
Sorting by `(ranking, from)` reorders pairs inside a node's group.

Two rows in 242 is the kind of error that survives every aggregate check: the
node count, the edge count, `fchr`, and every F bit were already exact while
those two characters were transposed.

### What this leaves

Everything the graph `.1.ht2` needs is now reproduced except the byte-level
emission: the side packing (2 rows/byte, F and M bitvectors, 6 `index_t` of
tallies per side) and the graph-mode `ftab`. Then the `.5`/`.6` local indexes.

The whole-genome question the AWS run is answering (three live arrays, ~570 GB
at 5.9e9 path nodes) is a property of exactly this loop, so an implementation
that matches it generation for generation is the precondition for changing where
those arrays live.

## Next rungs

2. ~~`22_20-21M.fa`, no variants — full index byte-identical to `hisat2-build`.~~
   **Done**, on two references, for `.1`–`.4`.
3. Same with SNPs and haplotypes — the first rung that needs `PathGraph`, and
   the one that brings `.5`/`.6` into scope.
4. chr1 against the three `analysis/hap` variant sets, matching measured retention.
5. Alignment equivalence: SAM byte-identical at `-p 1 --seed 0 --reorder`.
6. ~~Whole genome, where no C++ reference output exists.~~ **Done** — all
   eight files byte-identical against the E3 index (`design/e3`).

## Fixtures

`mkfixtures.sh <outdir>` regenerates the whole rung-3 fixture set, and
`verify.sh <outdir>` runs every emitter against it.

Each fixture input is *derived* in the script from files tracked in a git repo —
`example/reference/22_20-21M.fa` and its variant list, plus a 200 bp
hand-written case and `testdata/mkmulti.py` — so nothing is carried over from a
previous fixture directory. The derivations were checked against the original
set: all fifteen inputs reproduce byte for byte, including an `ex.snp` that had
been lost and turned out to be the repo's own `22_20-21M.snp`.

| prefix | reference | variants | exercises |
|---|---|---|---|
| `tinyidx` | 200 bp | 3 SNPs, 2 haplotypes | one window, one multi-variant haplotype |
| `xidx` | 200 bp | 1 deletion | the deletion path at minimum size |
| `multiidx` | 5 sequences, 15,700 bp | 55 SNPs | full-header refnames, leading/trailing/lowercase Ns |
| `cleanidx` | 509,431 bp, no ambiguity | 1,881 mixed | 10 windows, no N handling |
| `t_single` | same | 1,721 SNPs | substitutions alone |
| `t_insertion` | same | 74 insertions | **two variant-free windows → linear local indexes** |
| `t_deletion` | same | 86 deletions | one variant-free window |
| `exidx` | 1,000,000 bp with a 100,000 bp N run | 3,502 mixed | 18 windows, chromosome-vs-joined coordinates |

`mkfixtures.sh` rebuilds `hisat2-build-s` before it does anything else, and
refuses to run if the binary is older than any source file beside it. That guard
exists because the first fixture set was built by a binary months older than its
sources: `git status` was clean, repeated builds agreed with each other, and a
correct emitter looked wrong for two rounds of debugging. The manifest records
the source commit and the builder's hash next to the fixtures for the same
reason.

Regenerating the set immediately paid for itself — the two variant-free-window
cases had never been built before, and both crashed the emitter. See
`local_index_layout.md`.

## Whole-genome path

`ht2wg` builds a complete graph index — `.1` through `.8` — off the fragmented
graph builder and the external doubling loop, and writes it rather than
comparing against one. It is verified byte for byte on every fixture at both
index widths, and on a 20 Mb reference where it uses 197 MB against
`hisat2-build`'s 3.28 GB.

```bash
ht2wg <reference.fa> <snp> <haplotype> <workdir> <out_prefix> [budget_records] [verify_prefix]
```

It resumes: re-run the same command after a failure and it picks up from the
last checkpoint. `HT2_FRESH=1` starts over, `HT2_NO_RESUME=1` turns it off and
gives back the disk it costs.

`verify_wg.sh <fixture_dir>` checks the output; `verify_resume.sh <fixture_dir>`
crashes the build at every checkpoint boundary in every generation and checks
that finishing it lands in the same place. `HT2_LARGE=1` on any of these switches
to the 64-bit (`.ht2l`) half, which is what a whole-genome index is.

See `external_emitter.md` for the design and `whole_genome_plan.md` for the
measurements behind it.

### Whole genome, done

GRCh38 primary assembly with 14.9M SNPs and 16.3M phased haplotypes, 64-bit:
**all eight files byte-identical** to the index `hisat2-build` produced, in
17.7 h on one external SSD in ~11 GB against that build's measured 671 GB —
RSS sampled during the run rather than maximised over it, with the sort budget
capping it at 19 GB. All
23 generations matched its curve, ending at 5,917,131,871 path nodes. Scratch
peaked at ~325 GB, and in `generateEdges` rather than at the generation-11 node
peak (~139 GB) — size a disk from the emit stage, not the doubling.

Six construction bugs had to be fixed to get there, none reachable below
chromosome scale, since every fixture is a single contig with no N runs and no
co-located variants. Two in the global graph, three in local-index window
selection, and one that is not a defect: `.7`'s haplotype order comes from an
unstable `std::sort` under a comparator that admits `(left, right)` ties, so
`hisat2-build`'s own `.7` is not reproducible across standard library
implementations. `gnusort` reproduces libstdc++'s introsort to match it;
`HT2_HAPSORT=stable` opts out.

Two tools made that tractable, both exploiting stages that depend only on the
parse and not on the graph, so they answer in seconds against an index that
already exists instead of after a 17.7 h build:

```bash
ht2predict <reference.fa> <snp> <haplotype> <index.7.ht2l>       # rebuild and diff .7
ht2win     <reference.fa> <snp> <haplotype> <index.5.ht2l> <win>[,<win>...]
```

`HT2_REF_ONLY=1` on `ht2wg` stops after `.3`/`.4`, which is enough to prove the
reference going in is the one an existing index was built from before spending
the hours on the rest.
