# The pseudogene filter: built, rebuilt, and the hypothesis falsified

The T1/T2 report proposed that the SNP-aware index degrades ribosomal-protein
quantification because dbSNP entries inside processed pseudogenes are
paralogous sequence variants, and that removing them would fix it. The filter
is written and tested. **The rebuild shows the hypothesis is wrong**, and the
underlying finding was misdiagnosed.

## 1. What was built

`hisat2_filter_snps.py` — drops variants inside gene bodies of chosen biotypes
from a `.snp`/`.haplotype` pair before `hisat2-build --snp`. Four checks added
to `tests/run_tests.sh` (now 23 checks, all passing).

On the full prebuilt GRCh38 variant set (14,460,407 variants):

| | |
|---|---|
| masked | 10,660 processed / transcribed-processed pseudogenes, 16.9 Mb |
| dropped | **78,935 variants (0.55%)** |

*(An earlier figure of 181,160 / 1.40% in the T1 report was wrong — the regex
used there, `[a-z_]*processed_pseudogene`, also matched **un**processed
pseudogenes, which are 2.6k intron-retaining genes rather than retrocopies.)*

### A real hazard the tool has to defend against

Pseudogenes are not disjoint from real genes. **DHFRP2**, a transcribed
processed pseudogene of DHFR at chr6:31,360,865-31,367,354, lies across the 3'
end of **HLA-B** (31,353,872-31,367,067) and covers HLA-B's terminal exon — the
exact window 3' single-cell chemistry sequences. The naive filter stripped
**223 of HLA-B's 747 variants (29.9%)** to fix a ribosomal-protein problem.

So the tool protects exons of protein-coding genes by default (`--protect`).
After that fix:

| gene | variants | dropped | |
|---|---|---|---|
| HLA-A, C, DQA1, DQB1, DRA, DRB1, DRB5 | 82-2,011 each | **0** | |
| HLA-B | 747 | 193 | **all intronic — 0 inside any HLA-B exon** |
| HLA-DPA1 / DPB1 | 473 / 490 | 4 / 4 | |

Verified directly: 0 dropped variants fall inside HLA-B's 3,706 exonic bases.

## 2. What was rebuilt

A whole-genome SNP-aware build needs ~160 GB. A chromosome-scale build cannot
substitute, because processed pseudogenes are almost never on their parent's
chromosome (89 of 1,501 parent/pseudogene pairs). So the rebuild used a
**targeted 19.2 Mb reference** holding every locus involved — 83 ribosomal
protein genes, their 1,502 pseudogenes, and 250 unaffected expressed genes as
controls — with variant and GTF coordinates translated into it.

- Translation verified: reference base agrees with GRCh38 at **200/200** sampled
  variants, none equal to the ALT allele.
- HLA had to be left out: at 22-49 variants per 100 bp its local graphs are
  combinatorial, and `hisat2-inspect` cannot dump the haplotype file a prebuilt
  index was made with. Without haplotypes the build OOMs at ~200 GB. HLA is
  instead checked against the dropped-variant list, above.
- Three indexes on that reference — linear, all variants (92,975), filtered
  (86,523) — each retaining **98.9%** of variants in local indexes, so the arms
  are not confounded by different graph explosion.
- 16,715,717 reads, being every read either full-genome alignment placed in
  these loci.

**The mini reference reproduces the defect faithfully**: RPS27 keeps 0.128 of
its linear unique reads under the all-variants index, against 0.16 measured on
the real prebuilt index. Controls sit at 1.000.

### Result: the filter does not work

Reads whose unique (MAPQ 60) alignment starts inside each gene:

| group | genes | linear | snp-all | snp-filtered | all/lin | filt/lin |
|---|---|---|---|---|---|---|
| ribosomal protein genes | 83 | 5,066,190 | 4,173,844 | 4,260,252 | 0.824 | **0.841** |
| their pseudogenes | 1,498 | 191,107 | 842,591 | 823,202 | 4.409 | **4.308** |
| controls | 250 | 3,935,279 | 3,936,935 | 3,936,989 | 1.000 | 1.000 |

Across the 64 well-expressed RP genes the deficit is 891,783 reads and the
filter restores **86,416 of them — 9.7%**.

## 3. Why: the finding was misdiagnosed

Taking the 125,523 reads the linear index places uniquely at RPS27 and looking
at their competing alignments under the all-variants index:

| competing alignment | reads |
|---|---|
| **exact match, no ALT path used (no `Zs:Z`)** | **234,169 (88%)** |
| exact match via an ALT path | 30,737 (12%) |

Almost all the competition needs no variant at all. So it cannot be caused by
variants, and removing them cannot fix it — which is exactly the 9.7%.

**Ground truth settles it.** Counting exact occurrences of the read sequences in
the reference, independent of any aligner — 2,000 reads that the linear index
calls uniquely mapped at RPS27:

| exact copies in the reference | reads | |
|---|---|---|
| 1 (genuinely unique) | 145 | 7.2% |
| **2** | 1,210 | 60.5% |
| **3** | 172 | 8.6% |
| **4** | 192 | 9.6% |
| 0 (read carries a sequencing error) | 281 | 14.1% |

**78.7% of them match two or more places in the genome exactly.** Confirmed
against the whole GRCh38 assembly rather than the mini reference for a subset:
18 of 20 reads occur exactly 3 times genome-wide, 20/20 agreeing with the mini
reference count.

And on the **real prebuilt indexes**, whole genome, for those same 20 reads:

| index | result |
|---|---|
| `grch38` (linear) | **20/20 at MAPQ 60 with `NH:i:1`** |
| `grch38_snp` (graph) | 19/20 at MAPQ 1 with `NH:i:4` |

This is not a MAPQ-calibration nuance. `NH:i:1` is the aligner reporting that it
found **one** alignment for a read that matches three places in the assembly
exactly. **The graph index is right and the linear index is over-confident.**

### What this means for the T1 report

The "ribosomal protein pseudogene defect" is not a defect of the SNP-aware
index. It is the linear index failing to find multi-mappings that exist, and
the graph index finding them. The linear numbers that looked correct — RPSA
8,411 UMIs, RPS27 18,465 — are inflated by falsely-unique assignments; the
graph's lower numbers reflect genuine ambiguity.

That does not make the graph's output *useful* as-is: scattering a real gene's
reads across pseudogene copies as ties still collapses its count under
`--solo-multi-mappers Unique`. But the fix is a multi-mapper policy (EM, or
Uniform), not a variant filter, and the problem belongs to the annotation and
the reference, not to variant-aware alignment.

**This should be reported upstream**: HISAT2's linear index reporting MAPQ 60
for reads with 2-4 exact genomic matches is a correctness issue independent of
anything single-cell.

## 4. Status of the tool

`hisat2_filter_snps.py` is correct, tested, and does what it says. It is
**not** justified by the ribosomal-protein finding, which turned out to be
something else. Keep it for what it genuinely provides — the ability to exclude
variants from chosen regions, with protein-coding exons protected — and do not
apply it to the prebuilt index expecting a quantification improvement.

## 4b. It reproduces on a second donor and a different chemistry

A 10x PBMC 1k **5' v2** library (79.6M reads, different donor) run through the
same two prebuilt indexes shows the same thing: the top gainers are again
processed pseudogenes (RPS28P7 +8.61, RPS7P14 +8.08, RPL39P3 +7.38, RPS3AP6
12,452 vs 247 UMIs), and graph multi-mapping is **18.71% against linear's
15.98%** — a wider gap than the 3' library's 8.47% vs 6.46%.

That is what the splice-versus-retrocopy mechanism predicts. The competition is
between a spliced alignment at the parent gene and contiguous matches at its
retrocopies; nothing about it depends on which end of the transcript was
sequenced.

## 5. Caveats

- The mini reference is 19.2 Mb of extracted loci, not a genome. It reproduces
  the parent/pseudogene competition and the control genes behave, but absolute
  alignment rates (98.6%) are not comparable to a whole-genome run.
- The two SNP indexes were built without haplotypes, which the prebuilt index
  has. Both arms are affected identically, and both retained 98.9% of variants.
- The 14.1% of reads with 0 exact copies cannot be classified by exact matching;
  they may or may not be multi-mapping with mismatches, so 78.7% is a lower
  bound on the false-uniqueness rate only among reads that match exactly
  somewhere.
