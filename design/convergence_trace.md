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

## Consequence for the plan

**E0 on a large-memory machine is now more valuable, not less.** Its purpose has
changed: not "does doubling converge" (answered here — it does, cleanly) but
"where does the overflow actually occur". Building the external-memory sort-merge
join around an assumption about which structure blows up would be building
around a guess.

The cheap version of that question, worth trying first: run `hisat2-build` with
`--no-repeat-index` on the whole genome at 32-bit and see whether it still
overflows. If it completes, hypothesis 1 is confirmed for free and the repeat
index becomes the target. That still needs ~144 GB, so it is a one-instance
experiment rather than a laptop one — but it is a single decisive run rather than
an exploratory rental.

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
