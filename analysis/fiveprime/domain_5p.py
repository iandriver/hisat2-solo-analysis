#!/usr/bin/env python3
"""Does the graph's HLA gain sit in the peptide-binding domain? 5' vs 3'.

Exon indices are unreliable here -- the HLA transcripts in this annotation do
not all follow the textbook class I exon structure, and picking the "longest"
transcript gives a 1270 bp first exon for HLA-B where the real leader exon is
~73 bp. So the polymorphic region is defined on the MANE Select CDS by codon,
which is what the domain boundaries actually mean:

  class I  (A/B/C)      leader 24 aa, alpha1 1-90, alpha2 91-182
                        -> precursor codons 25-206
  class II beta (DRB1, DQB1)  leader ~29 aa, beta1 1-95  -> codons 30-124
  class II alpha (DQA1, DRA)  leader ~23 aa, alpha1 1-84 -> codons 24-107

Everything outside that window, within the same CDS, is the internal control.
"""
import subprocess, re
import numpy as np

H5 = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/fiveprime"
H3 = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
GTF = f"{H3}/genes.gtf"
BAMS = {"3'": (f"{H3}/graph.sorted.bam", f"{H3}/linear.sorted.bam"),
        "5'": (f"{H5}/5p_graph.bam", f"{H5}/5p_linear.bam")}

GENES = [("HLA-A", "6", "classI"), ("HLA-B", "6", "classI"), ("HLA-C", "6", "classI"),
         ("HLA-DRB1", "6", "classII_b"), ("HLA-DQB1", "6", "classII_b"),
         ("HLA-DQA1", "6", "classII_a"), ("HLA-DRA", "6", "classII_a"),
         ("B2M", "15", "control")]
WINDOW = {"classI": (25, 206), "classII_b": (30, 124), "classII_a": (24, 107),
          "control": None}


def mane_cds(gene, chrom):
    """CDS intervals of the MANE Select transcript, ordered 5'->3'."""
    out = subprocess.run(
        ["awk", "-F", "\t", '$1=="%s" && $3=="CDS" && /gene_name "%s"/' % (chrom, gene), GTF],
        capture_output=True, text=True).stdout
    tx, strand = {}, {}
    for line in out.splitlines():
        f = line.split("\t")
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        rank = 2 if 'tag "MANE_Select"' in f[8] else (1 if 'tag "Ensembl_canonical"' in f[8] else 0)
        tx.setdefault(m.group(1), {"rank": 0, "cds": []})
        tx[m.group(1)]["rank"] = max(tx[m.group(1)]["rank"], rank)
        tx[m.group(1)]["cds"].append((int(f[3]), int(f[4])))
        strand[m.group(1)] = f[6]
    if not tx:
        return None, None, None
    best = max(tx, key=lambda t: (tx[t]["rank"], sum(b - a + 1 for a, b in tx[t]["cds"])))
    cds = sorted(tx[best]["cds"])
    if strand[best] == "-":
        cds = cds[::-1]
    return best, cds, strand[best]


def codon_window_to_genomic(cds, strand, c1, c2):
    """Genomic intervals covering precursor codons c1..c2 inclusive."""
    lo, hi = (c1 - 1) * 3, c2 * 3          # CDS base offsets, half-open
    out, off = [], 0
    for a, b in cds:
        ln = b - a + 1
        s, e = max(lo, off), min(hi, off + ln)
        if s < e:
            if strand == "+":
                out.append((a + (s - off), a + (e - off) - 1))
            else:
                out.append((b - (e - off) + 1, b - (s - off)))
        off += ln
    return out


def cov(bam, chrom, ivs):
    if not ivs:
        return 0
    tot = 0
    for a, b in ivs:
        out = subprocess.run(
            ["samtools", "depth", "-a", "-Q", "60", "-d", "0",
             "-r", f"{chrom}:{a}-{b}", bam], capture_output=True, text=True).stdout
        tot += sum(int(l.split("\t")[2]) for l in out.splitlines())
    return tot


print(f"{'gene':<10}{'chem':<5}"
      f"{'--- peptide-binding domain ---':^32}{'--- rest of the CDS ---':^30}")
print(f"{'':<15}{'graph':>12}{'linear':>11}{'ratio':>8}{'graph':>12}{'linear':>11}{'ratio':>8}")
res = {}
for gene, chrom, cls in GENES:
    tx, cds, strand = mane_cds(gene, chrom)
    if cds is None:
        print(f"{gene:<10} no CDS found")
        continue
    L = sum(b - a + 1 for a, b in cds)
    if WINDOW[cls] is None:
        poly, rest = [], cds
    else:
        c1, c2 = WINDOW[cls]
        poly = codon_window_to_genomic(cds, strand, c1, c2)
        allb = codon_window_to_genomic(cds, strand, 1, L // 3)
        pol = set()
        for a, b in poly:
            pol.update(range(a, b + 1))
        rest = []
        for a, b in allb:
            run = None
            for p in range(a, b + 1):
                if p in pol:
                    if run:
                        rest.append(run); run = None
                else:
                    run = (run[0], p) if run else (p, p)
            if run:
                rest.append(run)
    for chem, (bg, bl) in BAMS.items():
        pg, pl = cov(bg, chrom, poly), cov(bl, chrom, poly)
        og, ol = cov(bg, chrom, rest), cov(bl, chrom, rest)
        res[(gene, chem)] = (pg, pl, og, ol)
        pr = pg / pl if pl else float("nan")
        orr = og / ol if ol else float("nan")
        print(f"{gene if chem=='3'+chr(39) else '':<10}{chem:<5}"
              f"{pg:>12,}{pl:>11,}{pr:>8.2f}{og:>12,}{ol:>11,}{orr:>8.2f}")
    print(f"{'':<10}     CDS {L} bp ({L//3} codons), transcript {tx}, strand {strand}")

print("\n=== how much of the CDS coverage each chemistry puts in the domain ===")
print(f"{'gene':<10}{'3-prime':>10}{'5-prime':>10}")
for gene, chrom, cls in GENES:
    if WINDOW[cls] is None or (gene, "3'") not in res:
        continue
    v = []
    for chem in ("3'", "5'"):
        pg, pl, og, ol = res[(gene, chem)]
        v.append(100 * pg / (pg + og) if (pg + og) else float("nan"))
    print(f"{gene:<10}{v[0]:>9.1f}%{v[1]:>9.1f}%")
