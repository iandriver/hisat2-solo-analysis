# A SNP-aware index never reports the spliced alignment when an equal-scoring contiguous match exists elsewhere

**HISAT2 2.2.3** (also reproduces with the prebuilt `grch38` / `grch38_snp` indexes from the HISAT2 download page), macOS arm64 and the same behaviour on the in-repo build.

## Summary

Build two indexes from the *same* FASTA, one plain and one with `--snp`. Take a
read that spans an exon-exon junction, where the gene also has a processed
pseudogene — a retrotransposed copy in which the two exons are contiguous. The
read then has two exact, equal-scoring placements: a spliced one at the parent
gene and a contiguous one at the retrocopy. Both alignments score `AS:i:0`.

- the **plain** index reports the spliced alignment at the parent gene, and only that
- the **`--snp`** index reports the contiguous alignment at the retrocopy, and only that — **the spliced alignment is never reported at all**, not even with `--known-splicesite-infile`

The `--snp` index is perfectly capable of the spliced alignment: delete the
retrocopy from the reference and it produces it (see control 3 below). It just
never surfaces it when a contiguous tie exists.

This is not a scoring difference — the alignments tie — but the consequence is
one-sided. Genes with processed pseudogenes lose their junction-spanning reads
to the pseudogene copies, which for RNA quantification means the real gene is
undercounted and a pseudogene that is generally not transcribed is credited
instead.

## Impact

Measured on 10x PBMC 1k v3 (66.6M reads, public), same reads through both
prebuilt human indexes, counting reads whose unique (MAPQ 60) alignment starts
inside each gene:

| gene set | plain index | `--snp` index | ratio |
|---|---|---|---|
| 83 ribosomal protein genes | 5,066,190 | 4,173,844 | **0.824** |
| their 1,498 processed pseudogenes | 191,107 | 842,591 | **4.41** |
| 250 expressed control genes | 3,935,279 | 3,936,935 | 1.000 |

Ribosomal protein genes are the worst affected because almost all of them have
many retrocopies, but nothing about the mechanism is specific to them.

Worst individual cases (fraction of the plain index's unique reads retained):

| gene | plain | `--snp` | retained |
|---|---|---|---|
| RPL18A | 53,019 | 4,989 | 0.094 |
| RPL17 | 17,337 | 1,724 | 0.099 |
| RPL9 | 88,412 | 10,886 | 0.123 |
| RPS27 | 125,523 | 16,084 | 0.128 |
| RPL3 | 100,165 | 13,680 | 0.137 |

At UMI level in a single-cell run this is severe, because a multi-mapping read
is discarded outright under a unique-only counting policy: RPSA drops from 8,411
UMIs to 341 (-96%), RPS10 3,500 to 219, RPL10 29,419 to 13,830.

## Reproducing

Everything below uses the prebuilt `grch38` and `grch38_snp` indexes and the
public 10x `pbmc_1k_v3` FASTQs, so no local files are needed.

Twenty reads that show it (91 bp, single-end, from `pbmc_1k_v3` R2). One of them:

```
@A00228:279:HFWFVDMXX:1:1214:10438:37043
CCAGAAGAGGAGAAGAGGAAACACAAGAAGAAACGCCTGGTGCAGAGCCCCAATTCCTACTTCATGGATGTGAAATGCCCAGGATGCTATA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

That sequence occurs **exactly three times** in GRCh38 (verified by direct
string search of the assembly, both strands), at chr1:37,556,279,
chr1:202,471,993 and chr12:3,211,793 — all processed pseudogenes. It does *not*
occur contiguously at RPS27, because at RPS27 it spans a 342 bp intron.

```bash
hisat2 -x grch38/genome     -U probe20.fq --no-unal -S plain.sam
hisat2 -x grch38_snp/genome_snp -U probe20.fq --no-unal -S snp.sam
```

| | result |
|---|---|
| `grch38` (plain) | 20/20 at `chr1:153,991,142`, CIGAR **`82M342N9M`**, `NM:i:0`, MAPQ 60, `NH:i:1` — spliced at RPS27 |
| `grch38_snp` | 19/20 at MAPQ 1 with `NH:i:4`, all CIGAR `91M` at the pseudogene loci; **no alignment at RPS27** |

### Controls

**1. It is not `-k` truncation.** `-k 5`, `-k 10`, `-k 50` all give the plain
index the same `NH:i:1`, and the `--snp` index the same set without RPS27.

**2. It is not the annotation.** Supplying the junction with
`--known-splicesite-infile` (from `hisat2_extract_splice_sites.py`) makes the
plain index 20/20 spliced, and changes the `--snp` result not at all —
still zero alignments at RPS27.

**3. The `--snp` index can splice here.** Rebuild both indexes from a reference
containing only the RPS27 locus, so the pseudogene copies are absent:

| index | default | with `--known-splicesite-infile` |
|---|---|---|
| plain | 19 spliced, 1 contiguous | **20 spliced** |
| `--snp` | 19 spliced, 1 contiguous | **20 spliced** |

Identical. The spliced alignment is found and reported as soon as nothing ties
with it.

**4. Same sequence, only the index differs.** All of the above was also run on a
19.2 Mb reference holding the RP genes, their pseudogenes and 250 controls, with
a plain index and a `--snp` index built from that one FASTA. Same split: plain
19/20 spliced at RPS27, `--snp` 0/20 at RPS27, 77 contiguous records elsewhere,
every alignment `AS:i:0`. Both indexes retained 98.9% of variants in their local
indexes, so this is not graph explosion.

## What I could not isolate

A minimal synthetic case did not reproduce it. Building a two-contig reference
with an exon1-intron-exon2 gene and a contiguous retrocopy, and a read spanning
the junction, **both** indexes prefer the contiguous copy — with and without
`--known-splicesite-infile`. So the trigger involves something about the real
loci (three or four competing copies rather than one, or the surrounding
sequence) that I have not pinned down. The real-data reproduction above is
reliable and fully controlled, but a maintainer will see the mechanism faster
than I can guess at it.

A plain duplicate is handled correctly, for what it is worth: a 1 kb block
placed on two contigs, with a read from inside it, gives `NH:i:2` and MAPQ 1
from both index types.

## Why it matters beyond this case

Whichever alignment is preferred, the reported set is incomplete on both sides:
the plain index reports one of several exact matches with `NH:i:1`, and the
`--snp` index omits an equal-scoring spliced alignment. Reporting all
equal-scoring alignments — or at least not systematically excluding the spliced
one — would make MAPQ meaningful here and let downstream tools apply their own
multi-mapper policy. As it stands, users of a SNP-aware index silently lose a
large fraction of the expression of every gene that has a processed pseudogene,
with no warning and nothing in the output to indicate it.
