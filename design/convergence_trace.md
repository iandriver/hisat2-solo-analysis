# E0 — does prefix doubling converge, and how big does it get?

Measured instead of rented. `PathGraph::printInfo` (`gbwt_graph.h:2301`) already
logs `Generation N (temp -> nodes, ranks)` under `--verbose`, and node counts are
a property of the graph rather than of `index_t` width — so 32-bit builds on
subsets give the same curve a 64-bit whole-genome build would, at a fraction of
the memory. Input is the whole-genome 1000 Genomes phased set (14,947,745
variants, 16,338,502 haplotypes) filtered per chromosome.

## It converges, and there is no transient explosion

| chrom | bp | variants | initial | peak temp | final | **peak/final** | gens | nodes/bp |
|---|---|---|---|---|---|---|---|---|
| 1 | 248,956,422 | 1,182,916 | 238,338,543 | 265,072,635 | 261,676,691 | **1.0130** | 16 | 1.051 |
| 2 | 242,193,529 | 1,236,990 | 248,283,321 | 272,740,968 | 268,626,757 | **1.0153** | 14 | 1.109 |
| 6 | 170,805,979 | 964,230 | 178,980,281 | 200,513,642 | 195,440,675 | **1.0260** | 17 | 1.144 |
| 20 | 64,444,167 | 343,906 | 66,385,047 | 72,940,441 | 71,490,614 | **1.0203** | 15 | 1.109 |
| 21 | 46,709,983 | 216,745 | 41,785,256 | 47,411,040 | 45,320,912 | **1.0461** | 18 | 0.970 |
| 22 | 50,818,468 | 224,119 | 40,925,912 | 58,659,139 | 57,776,238 | **1.0153** | 14 | 1.137 |

**The construction peak never exceeds the finished size by more than 4.6%**,
across a 5x range of chromosome length and a 28% range of variant density.
Generation counts sit in a tight 14-18 band.

chr6 — the MHC chromosome, highest density at 5.65 variants/kb — is not an
outlier (1.0260). Density does not produce a transient blowup.

**This is the result that matters for an external-memory builder.** There is no
intermediate spike to absorb: the algorithm holds three copies of a set that
grows ~8-18% over the whole doubling process. Streaming two of those copies to
disk costs I/O and essentially nothing in space.

## Combining chromosomes costs generations, not nodes

Building chr20+21+22 as one graph versus three:

| | separate | combined | ratio |
|---|---|---|---|
| initial nodes | 149,096,215 | 149,096,211 | 1.0000 |
| final nodes | 174,587,764 | 176,045,586 | **1.0084** |
| peak temp | 179,010,620 | 177,005,097 | 0.9888 |
| generations | 14 / 15 / 18 | **22** | — |
| peak/final | 1.0253 | **1.0055** | — |

Sequence repeated across chromosomes has to be disambiguated against the whole
collection, so the combined graph needs **eight more generations** — but only
**0.84% more nodes**, and its peak/final is *tighter*. Per-chromosome numbers
are additive to within 1%.

## The extrapolation contradicts the observed failure

Summing the six chromosomes gives 1.0927 nodes/bp. Applied to the 3,099,750,718
bp primary assembly and scaled by the cross-chromosome inflation factor:

```
naive projection    3,387,097,609 nodes   78.9% of 2^32
x inflation         3,415,380,149 nodes   79.5% of 2^32
x inflation, peak   3,433,995,184 nodes   80.0% of 2^32
```

**Under the ceiling.** Yet the real 32-bit whole-genome build failed with
`exceeded integer bounds` twice, at 14.95M and at 12.64M variants.

The measurements say the main graph should have fit. Rather than retune the
growth constant until the prediction matches the known outcome, record the
contradiction. Three candidate explanations:

1. **The overflow is in the repeat index, not the main graph.** `hisat2-build`
   builds both; `rfm.h` runs its own `PathGraph`; the error string lives in
   `gbwt_graph.h`, which serves both paths. Repeat families are precisely where
   path counts would explode, and a per-chromosome trace of the main graph would
   never show it.
2. **Growth turns superlinear past a few sequences** in a way three chromosomes
   cannot reveal.
3. Those particular stage-2b runs differed in some parameter now lost — the S3
   scratch holding their logs was deleted during teardown.

Discriminating between these needs more memory than this machine has: chr2 alone
peaked at 19.65 GB, and a whole-genome 32-bit build needs ~144 GB for the node
arrays alone.

## The whole-genome curve — measured, and it resolves the contradiction

Run on AWS (r7i.12xlarge, 48 cpu / 371 GB, ~$1.20 of instance time), 32-bit,
GRCh38 primary assembly + the full 14,947,745-variant phased set, exactly the
stage-2b configuration.

| gen | temp nodes | nodes | % of 2^32 | ranks/nodes |
|---|---|---|---|---|
| 0 | 3,051,980,501 | 3,051,980,501 | 71.1% | — |
| 3 | 3,166,890,352 | 3,166,890,352 | 73.7% | — |
| 5 | 3,498,007,929 | 3,315,396,267 | 81.4% | 0.848 |
| 8 | 3,436,077,182 | 3,390,374,563 | 80.0% | 0.980 |
| 9 | 3,635,922,500 | 3,595,956,954 | 84.7% | 0.962 |
| **10** | **> 4,294,967,295** | — | **overflow** | — |

**The extrapolation above was wrong, and the error is quantified.** Whole-genome
growth is **1.453 nodes/bp against 1.0927 per-chromosome — 33% higher**.
Superlinear growth in collection size is confirmed; the three-chromosome test
(+0.84%) sampled far too little of the genome's repeat structure to see it.

Projecting chr22's generation-10 jump (1.2524x) onto the measured generation-9
count gives a final of **~4.50e9 nodes, ~105% of the 2^32 ceiling**. The 32-bit
build misses by roughly 5%.

### What this settles

- **`index_t` must be 64-bit.** Measured, not assumed, and the margin is small
  enough that measuring was worth it.
- **No variant filtering rescues the 32-bit build.** Variants and haplotypes
  contribute only ~3.3% of the initial node count; the reference dominates. That
  is why stage 2b failed identically at 14.95M and at 12.64M variants.
- **The 64-bit requirement is now sized.** At ~4.5e9 nodes, 32 B per `PathNode`,
  three live arrays in `lateGeneration`: **~400 GB for `PathGraph` alone**. That
  is a precise explanation of the observed 368.8 GB OOM — it died inside that
  allocation. With GFM and local-index construction on top, an in-memory 64-bit
  build wants 500-700 GB.
- **The property the external-memory design depends on holds at full scale**:
  within-generation temp/nodes stayed between 1.01 and 1.06. No spike to absorb.

Measured 32-bit peak RSS before the overflow: **176 GB** (184,548,364 KB).

### Method note

Three launches were needed. The first invoked `hisat2-build-s` (the binary) with
`--verbose`, which only the Python wrapper accepts; the second failed because
the output directory did not exist. Both died in seconds and self-terminated, so
the waste was ~$0.40. The third added a pre-flight that builds a 1.2 Mb chr1
slice with real variants and requires >=3 `Generation` lines before committing to
the long run — a gate that would have caught both.

Worth recording separately: `--verbose` is a *wrapper* option and is not passed
to the binary at all. The binary is verbose by default (`-q` disables it), which
is why the local traces produced curves despite the wrapper stripping the flag.

## Consequence for the plan

**Answered.** The overflow is in the main graph, at generation 10, needing
~4.5e9 path nodes against a 4.295e9 ceiling. Hypothesis 1 (the repeat index) was
wrong at the premise: `hisat2-build` does not build a repeat index at all unless
one is supplied via `--repeat-ref`, and stage 2b supplied none. Hypothesis 2,
superlinear growth with collection size, is confirmed and measured at +33%
nodes/bp over the per-chromosome figure.

## Peak RSS observed (32-bit, `-p 8`, includes GFM and local index construction)

| input | peak RSS |
|---|---|
| chr21 | 10.05 GB |
| chr22 | 11.53 GB |
| chr20 | 13.58 GB |
| chr6 | 17.40 GB |
| chr1 | 18.66 GB |
| chr2 | 19.65 GB |
| chr20+21+22 | 20.90 GB |

Note these are far above the node arrays alone (3 x N x 16 B), so `hisat2-build`
peak is not dominated by `PathGraph` at these scales — the GFM and local index
construction contribute substantially. That is worth remembering when sizing the
external-memory design: bounding `PathGraph` does not by itself bound the build.

## Reproducing

```
trace.sh 22 21 20 6 1 2          # per-chromosome curves
trace_combo.sh combo202122 20 21 22   # additivity test
```

Raw per-generation curves are in `trace/gen_*.txt`.
