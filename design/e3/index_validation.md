# P0 — is the whole-genome index actually usable?

**Result: yes, on every axis tested.** 1,000,000 reads in 10.5 s at `-p 8`,
99.93% aligned, 96.4% placed within 1 bp of truth, 13.65 GB peak RSS, and
alt-allele reads aligning at `NM:i:0` with their 1000 Genomes IDs attached.

The question is not rhetorical. `order_bounding_test.md` records a chr22 build
that exited 0, shrank the index 47%, and produced something the aligner could
not query — 100 reads unfinished after 120 s. "rc=0 and 11 GB on disk" is not
evidence of a working index, so this is the check before anything depends on it.

## Results

Reads are exact 150 bp substrings of the genome **dumped from the index itself**
via `hisat2-inspect`, with the true position in the read name. Drawing from the
index's own sequence removes every assembly-version question (see the trap
below). Machine: M5 Pro, 52 GB.

| test | reads | aligned | placed within 1 bp | wall | peak RSS |
|---|---|---|---|---|---|
| smoke, `-p 1` | 100 | 100.00% | — | 5.6 s | 6.67 GB |
| genome-wide, exact | 9,200 | 99.95% | **96.15%** | 4.0 s | — |
| genome-wide, 2 mismatches | 9,200 | 98.63% | **93.21%** | 3.6 s | — |
| genome-wide, exact | **1,000,000** | 99.93% | **96.36%** | **10.5 s** | **13.65 GB** |

~95,000 reads/s at `-p 8`. Against E2's failure mode — 100 reads unfinished in
120 s — this is 1,000,000 reads in 10.5 s.

**The residual ~4% is repeat ambiguity, not error.** A 150 bp exact substring
sometimes has more than one equally good placement. The reference point is the
known-good 900 kb chr22 graph index on the same style of test: 95.92% exact and
93.80% at 2 mismatches. The whole-genome index scores 96.15% and 93.21% — **no
degradation, despite 3,400x more sequence to be confused by.** Misplacements
cluster where they should: 1q21 segmental duplications, and chr1 pericentromere
to chr19.

## The variants work, which is the entire point

2,000 chr1 SNPs isolated by 200 bp, each read built twice — once carrying the
reference base, once the alternate — with the variant at read offset 75.

| | aligned | placed within 1 bp | `NM:i:0` | carrying `Zs:Z` |
|---|---|---|---|---|
| reference allele | 99.90% | 95.40% | **100.00%** | 0.15% |
| alternate allele | 99.95% | **99.25%** | **100.00%** | **99.65%** |

```
v0_chr1_10508  0  chr1  10508  60  150M  ...  NM:i:0  MD:Z:75G74  Zs:Z:75|S|1:10583:G:A
```

`MD:Z:75G74` says the reference base at offset 75 differs; `NM:i:0` says nothing
was charged for it. The read carries a non-reference allele and pays no
mismatch penalty, and `Zs:Z` names the variant with its real 1000 Genomes ID.
This is the mechanism a linear reference cannot provide and that STAR's
`--varVCFfile`/WASP can only detect after the fact.

An unforced bonus: **alt-allele reads place *better* than reference-allele ones**
(99.25% against 95.40%). Carrying a variant makes a read more unique, so it
distinguishes its true locus from paralogous copies.

## Memory

13.65 GB peak at `-p 8`, against 9.0–9.4 GB measured for the distributed
`grch38_snp` in the head-to-head benchmark. The extra ~4.5 GB buys 64-bit ids
and real phasing. Index load alone is 6.65 GB (`.1` + `.2`); the rest is local
indexes and per-thread buffers.

## Two traps, both worth recording

**1. `example/reference/22_20-21M.fa` is not current GRCh38.** The first attempt
scored 0.0% placed-within-1-bp against the whole-genome index, which looks like
catastrophic failure. It was not: the offsets were dominated by a handful of
*constants* (+12,477 at 36.8%, −354,287 at 25.1%, −1,849,867 at 11.2%), which is
the signature of a coordinate error, not a broken index. Locating the slice in
this build:

```
slice starts at chr22:20,012,478   (not 20,000,001 -- off by 12,477)
matches the genome over its full length:  False
```

So it is from a different assembly revision and only partially matches. It also
comes from 22q11.2, which is full of low-copy repeats — hence the several
distinct constant offsets. **Do not use the bundled example slice as positional
ground truth against a GRCh38 index.** Dump the reference from the index.

**2. The realism check is still open.** The local read set (GSE236579) aligned at
0.02% and 0.93% for its two mates. Those reads are not human genomic — the names
end `/2` and `/3`, a three-read layout — so this says nothing about the index and
the file was set aside. A known-good human RNA-seq sample is still needed for
throughput and alignment rate on real data; the numbers above are from perfect
simulated reads and are optimistic (no splicing, adapters, or quality issues).

## Reproducing

```
aws s3 sync s3://rustar-bench/hisat2-wg64/index/ wg64_idx/
hisat2-inspect wg64_idx/genome > genome_dump.fa      # 2.9 GB, 194 sequences
# sample 150 bp windows, encode chrom+pos in the read name, align, score
```
