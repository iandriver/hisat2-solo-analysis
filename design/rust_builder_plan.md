# A Rust graph-index builder for HISAT2 — design

## Why the builder and not the aligner

The aligner works. It is fast, its index is 4.5x smaller than STAR's, its memory
held to within 4% across runs whose wall time spanned 10x, and its alignment
rate reproduces to four decimal places across months (87.20%, and all four LCL
rates in H2). Rewriting working code buys risk.

The builder does not work. It cannot build a human graph index at all, it
silently discards ~19% of the variants it is given, and the bound it hits is
spelled `(index_t)-1`.

## What must not change

**The `.ht2` format.** Byte-compatible output means the existing aligner runs
unmodified and every prebuilt JHU index stays valid. It also gives the strongest
possible correctness test: build the same input with both builders and `cmp`
the files. Format is a fixed binary layout with an endianness sentinel
(`gfm.h:1397`), `index_t`-width fields, files `.1.ht2`/`.2.ht2` for the global
index and `.7.ht2`/`.8.ht2` for variants and haplotypes.

Anything that breaks the format needs an aligner change and loses that test, so
it belongs in a later stage, behind evidence.

## The three failure modes, with their code sites

| # | failure | site | measured |
|---|---|---|---|
| 1 | global node-id overflow | `gbwt_graph.h:2038`, `:2120` | fails at 72.5% of 2^32 |
| 2 | 64-bit build exhausts RAM | `PathGraph` holds `nodes` + `past_nodes` | OOM at 368.8 GB, twice |
| 3 | local windows silently drop variants | `hgfm.h:1941` vs `local_max_gbwt` | 80.9-92.4% retention; 83.9% of MHC variants at risk |

## Experiment 0 — partly done, and it changed the question

**Update: the convergence half is measured. See `convergence_trace.md`.**

`printInfo` already logs the per-generation curve under `--verbose`, so the
convergence question was answerable on a laptop rather than a rented 1 TB box.
Six chromosomes plus a combined build establish:

- doubling **converges cleanly**, 14-18 generations per chromosome, 22 combined;
- the construction **peak never exceeds the finished size by more than 4.6%** —
  there is no transient explosion to absorb, which is what makes the
  external-memory design viable;
- combining chromosomes costs **8 more generations but only 0.84% more nodes**,
  so per-chromosome numbers are additive;
- chr6, the densest chromosome, is **not** an outlier.

It also produced a contradiction: the measured curve extrapolates to ~3.43e9
nodes, **80% of the 2^32 ceiling**, when the real 32-bit build overflowed. So
the main graph should have fit and did not. The leading hypothesis is that the
overflow is in the **repeat index** (`rfm.h`, its own `PathGraph`, sharing the
error string in `gbwt_graph.h`) rather than in the main graph.

**E0 therefore still needs to run, with a changed purpose:** not "does it
converge" but "where does it overflow". The decisive cheap form is a single
whole-genome 32-bit build with `--no-repeat-index`; if that completes, the
target is the repeat index, not the main doubling loop. It needs ~144 GB, so it
is one instance-hour, not an exploratory rental.

Note also that measured peak RSS (10-21 GB for single chromosomes) sits far
above the node arrays alone, so `PathGraph` does not dominate `hisat2-build`
peak at these scales. **Bounding `PathGraph` will not by itself bound the
build** — the GFM and local index construction have to be accounted for too.

## Experiment 0 (original framing, superseded above)

**Rent a 1-2 TB machine and run the existing 64-bit C++ builder to completion.**

Both prior 64-bit attempts were killed at the *machine* ceiling (386,797,904 and
386,805,112 KB — the same number twice), so we never learned the actual
requirement; we only learned it exceeds 368.8 GB. Everything below assumes the
count is bounded by graph structure rather than exploding, and that assumption is
worth $50 to test rather than designing around.

- `x2iedn.8xlarge` (1 TB) or `x2iedn.16xlarge` (2 TB), a few hours, ~$25-100.
- Instrument peak RSS per generation, and log `temp_nodes` at each doubling step.

Three possible outcomes, each changing the plan:

| outcome | what it means | consequence |
|---|---|---|
| completes under ~600 GB | bounded, engineering problem | build the Rust builder for external memory; **also we get a usable index today** |
| completes but needs >1 TB | bounded but brutal | external memory is mandatory, order-bounding attractive |
| does not converge | node count explodes after all | order-bounding (GCSA2) becomes necessary, not optional |

The per-generation `temp_nodes` curve is the real deliverable here: it tells us
whether doubling is converging and what the peak actually is. That single trace
determines which of the two architectures below we build.

## Architecture

Three phases, matching the existing pipeline so byte-identity stays testable.

### Phase A — reference graph

Stream FASTA + `.snp` + `.haplotype` into the node/edge graph that
`RefGraph::RefGraph` builds today (`gbwt_graph.h:636-767` for the haplotype
path construction). Linear in genome size, ~3.1e9 nodes for human, trivially
parallel per chromosome. Serialise to disk in a compact form (2-bit labels,
delta-coded edge lists) so Phase B can stream it.

Nothing clever is needed here; it is the cheap part and it is where the input
validation belongs.

### Phase B — global GFM

This is the whole problem. Two designs, chosen by Experiment 0.

**B1 — external-memory prefix doubling, exact (preferred).**

Keep HISAT2's semantics exactly: full disambiguation, same output. Change only
where the data lives.

- `index_t` = `u64` throughout; no 32-bit path at all for the global index.
- `PathNode` is `{from: u64, to: u64, key: (u64, u64)}` = 32 B. Sortable by a
  plain u64 key, which means **LSD radix sort over memory-mapped chunks**, not a
  comparison sort.
- Each doubling generation is a sort + a join. Both are external: sort runs of
  `chunk_bytes` in RAM, spill, then k-way merge. Peak RAM becomes a *parameter*
  (`--build-memory`), not a consequence.
- Disk: ~2 x 137 GB for the node arrays at 4.3e9 nodes, so budget 400-600 GB.
  Cheap and, unlike RAM, elastic.
- `nodes`/`past_nodes` never coexist fully in RAM, which is exactly what kills
  the C++ builder today.

Risk: I/O bound. Mitigated by radix sort (sequential passes), memory-mapped
chunk files, and the fact that a 2.5 hour build is already acceptable — STAR
takes 2 h 26 m for a *linear* index.

**B2 — order-bounded, GCSA2-style (fallback / optimisation).**

Sirén's GCSA2 caps the order at k, indexing all paths of length <= k rather than
disambiguating fully. Construction becomes bounded by construction, and the index
shrinks.

The cost is semantic: queries longer than k can return false positives, which
must be verified against the graph. **Whether that matters here is an empirical
question with a specific answer.** HISAT2 searches the global index by
`ftabLoHi` followed by repeated `mapLF` (`hi_aligner.h:6438-6469`) — a maximal
exact match, so query length is bounded by read length, not by a seed constant.

So before adopting B2:

**Experiment 1 — DONE, and it reorders this whole plan. See
`query_lengths.md`.**

The longest query ever issued against the global index was **131 bp** (150 bp
reads; 86 bp on 91 bp reads), mean 24, median 18. A maximal exact match cannot
exceed its read, so the distribution is structurally bounded.

Against the E0 construction curve, where `generation` means "sorted by paths of
length 2^generation" and construction runs to *full* disambiguation:

| generation | order | nodes | % of 2^32 |
|---|---|---|---|
| 8 | 256 | 3,390,374,563 | **78.9% — fits** |
| 9 | 512 | 3,595,956,954 | **83.7% — fits** |
| 10 | 1024 | overflow | — |

**Stopping the doubling at generation 8 or 9 puts a whole human graph index
inside the 32-bit ceiling and is exactly equivalent for every query HISAT2
issues.** Order 256 clears the observed maximum by 2x. For contrast, chr22 needs
order 8192 to fully sort — the index is built to roughly two orders of magnitude
more resolution than anything ever asks for.

So B2 is not a fallback. It is the primary lever, and B1 demotes to an
optimisation:

1. **Bound the order** — human whole-genome graph index becomes buildable at
   32 bits in ~176 GB, no new integer width, no external memory, no id format
   change.
2. **External memory** — now about build-machine cost, not about whether the
   index can exist.
3. **64-bit** — unnecessary for human; still needed for pangenomes.

The unpaid cost: `generateEdges` and the GFM build run *after*
`while(!isSorted())` and assume a fully sorted `PathGraph`. Stopping early leaves
nodes sharing a 2^g prefix merged, and the rest of the pipeline must handle it —
which is precisely what GCSA2 implements, so known-possible, but real work. A
naive `break` out of the doubling loop would produce a wrong index, not a
bounded one.

### Phase C — local indexes, and the MHC problem

55,000-odd windows of 57,344 bp, offsets in `u16`, budget `local_max_gbwt` =
63,488 edges. Embarrassingly parallel; each window is tiny. This is where the
silent dropping lives, and where the MHC finding lands: 55.1% of MHC windows
over budget, 83.9% of MHC variants at risk, against 19.1%/28.9% genome-wide.

Four options, in increasing order of intrusiveness:

| option | fix | aligner change | keeps format |
|---|---|---|---|
| C0 | binary-search backoff + **report every dropped variant** | none | yes |
| C1 | adaptive window size — shrink dense windows until they fit | window lookup table | no |
| C2 | per-window `u16`/`u32` offset type | dispatch on a flag | no |
| C3 | haplotype-restricted local graphs | none | yes |

C0 is already implemented in the C++ fork (87.5% and 97.8% retention, up from
80.9% and 92.4%) and is the zero-risk baseline the Rust builder must at least
match. **The reporting half is the part upstream #473 is actually about** — a
finished index today gives no way to discover what it lost.

C2 is the surgical fix: `LocalGFM` is already templated on `local_index_t`
(`hgfm.h:1557`), so a `u32` instantiation costs a type parameter, not an
algorithm. Only ~19% of windows need it and 0.18% of windows are MHC, so the
size cost is small and lands where the value is. This is my preferred target
once B1 works.

C1 is more elegant — the MHC gets *more* local indexes rather than fewer
variants — but `getLocalGFM` computes `offset / local_index_interval`
(`hgfm.h:1719`), a fixed stride, so variable windows need a lookup structure in
the index and a hot-path change in the aligner. Later, if C2 proves insufficient.

## What I am deliberately not proposing

- **A full HISAT2 rewrite.** rustar showed a faithful Rust port matches on
  accuracy (r = 0.99887) and gains modestly on speed, but used the *same* memory
  as STAR — because memory is set by the data structure. There is no memory win
  available from the language.
- **r-index / run-length BWT.** A human graph over one reference is extremely
  repetitive, so an RLBWT would likely be dramatically smaller. It also replaces
  the index wholesale, which means a new aligner. Interesting, out of scope, and
  worth revisiting only if B1 fails.
- **Prefix-free parsing.** Same objection, less maturity for graphs.
- **GBWT as a replacement.** GBWT indexes haplotype *paths*, not all substrings;
  in vg it sits alongside GCSA2 rather than replacing it. HISAT2 already ingests
  haplotypes into the graph, so the idea is partly present.

## Validation — the byte-identity ladder

Each rung must pass before the next is attempted.

1. **Round-trip.** Read an existing `.ht2` in Rust, re-emit it, `cmp` clean. No
   construction. Proves format mastery and gives the writer a test harness.
2. **`example/reference/22_20-21M.fa`, no variants.** Rust builder output
   byte-identical to `hisat2-build`. Proves Phase A + B on the linear case.
3. **Same, with SNPs and haplotypes.** Byte-identical. Proves the graph path.
4. **chr1 with the three variant sets from `analysis/hap`.** Byte-identical,
   *and* retention matching the measured 80.9 / 86.7 / 92.4% (stock) or
   87.5 / 97.8% (binary search). Those numbers are already banked, so this rung
   is a regression test against real measurements.
5. **Alignment equivalence.** Same reads through Rust-built and C++-built chr1
   indexes; SAM byte-identical at `-p 1 --seed 0 --reorder`. Catches anything
   `cmp` on the index would miss.
6. **Whole genome.** No C++ output to compare against — that is the point. Falls
   back to: the probe's predicted retention, alignment rate sanity against the
   87.20% PBMC figure, and the H2 HLA ratios.

Rungs 1-5 are all reachable on a laptop in minutes. That matters: the current
C++ builder's feedback loop for a real test is hours.

## Tooling

- `memmap2` for chunk files; `rayon` for the parallel phases.
- **No off-the-shelf succinct library.** We must emit HISAT2's exact byte
  layout, so `sucds`/`vers` serialization is unusable. The structures (packed
  2-bit BWT, occ tables, `ftab` of 4^k entries, `offs`) are simple enough to
  write directly, and writing them is what rung 1 tests.
- Custom LSD radix sort on u64 keys — the standard `sort_unstable_by_key` is the
  wrong tool at 4e9 elements.
- `.snp` / `.haplotype` / `.fa` parsing is trivial TSV/FASTA; no `noodles`
  dependency needed for the builder proper.
- `proptest` for the graph invariants, which is where Rust earns its keep over
  the C++ original's zero unit tests.

## Staging

| stage | content | gate |
|---|---|---|
| **E0** | Experiment 0 on a 1-2 TB box; per-generation `temp_nodes` trace | picks B1 vs B2 |
| **E1** | Instrument global-index query lengths | tells us if B2 is free |
| **S1** | Format round-trip (rung 1) | `cmp` clean |
| **S2** | Phase A + linear global index (rungs 2) | `cmp` clean |
| **S3** | Graph global index, in-memory, 32-bit (rungs 3-5) | `cmp` clean + SAM identical |
| **S4** | Local indexes with full drop reporting (C0) | retention matches banked chr1 numbers |
| **S5** | External-memory, u64 (B1) | whole-genome build completes |
| **S6** | C2 `u32` local windows | MHC retention measured, aligner patched |

S1-S4 are laptop work and deliver value even if S5 fails: a builder that
*reports* what it drops closes #473 regardless of whether it can build a human
graph.

## Risks

1. **Experiment 0 shows non-convergence.** Then B1 is dead and the project
   becomes "port GCSA2's order-bounded construction", which is a research task,
   not an engineering one. This is why E0 comes first and costs $50.
2. **Byte-identity may be unreachable** if the C++ builder's output depends on
   thread count or iteration order anywhere. Rungs 2-3 will expose it early. If
   it is unreachable, fall back to semantic equivalence (rung 5) and accept a
   weaker test.
3. **`-fno-strict-aliasing` is load-bearing in the C++ original** (`gfm.h:3180`
   type-puns). When reimplementing, the punned layouts must be reproduced
   explicitly rather than inferred from reading the C++ — a place to be
   suspicious of one's own translation.
4. **Scope creep into the aligner.** C1/C2 and B2 all eventually touch it. The
   discipline is that S1-S5 must not, so the builder stays independently
   shippable.
5. **Two builders to maintain** if this lands upstream. Mitigated by the format
   being frozen and the C++ builder remaining the reference for small inputs.
