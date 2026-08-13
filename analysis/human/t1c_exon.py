#!/usr/bin/env python3
"""T1c: within an HLA gene, where do the graph's extra reads land?

For class I the polymorphic residues are encoded by exons 2 and 3; for class II
beta chains by exon 2. If reference bias is the mechanism, the coverage gain
should concentrate there rather than spreading evenly over the gene.
"""
import subprocess, re
import numpy as np

D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f549cf1f37/scratchpad/human"
D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"

GENES = [
    ("HLA-A",    "6", 29941260, 29949572, "+", "classI"),
    ("HLA-B",    "6", 31353872, 31367067, "-", "classI"),
    ("HLA-C",    "6", 31268749, 31272130, "-", "classI"),
    ("HLA-DQA1", "6", 32628179, 32647062, "+", "classII_a"),
    ("HLA-DQB1", "6", 32659467, 32668383, "-", "classII_b"),
    ("HLA-DRB1", "6", 32577902, 32589848, "-", "classII_b"),
    ("HLA-DRA",  "6", 32439878, 32445046, "+", "invariant"),
]
POLY = {"classI": {2, 3}, "classII_a": {2}, "classII_b": {2}, "invariant": set()}


def canonical_exons(gene, chrom, start, end):
    awk = ('$1=="%s" && $3=="exon" && $4>=%d && $5<=%d && /gene_name "%s"/'
           % (chrom, start, end, gene))
    out = subprocess.run(["awk", "-F", "\t", awk, f"{D}/genes.gtf"],
                         capture_output=True, text=True).stdout
    tx = {}
    for line in out.splitlines():
        f = line.split("\t")
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if m:
            tx.setdefault(m.group(1), []).append((int(f[3]), int(f[4])))
    if not tx:
        return []
    return sorted(max(tx.values(), key=lambda ex: sum(b - a + 1 for a, b in ex)))


def depth(bam, chrom, start, end):
    out = subprocess.run(
        ["samtools", "depth", "-a", "-Q", "60", "-d", "0",
         "-r", f"{chrom}:{start}-{end}", f"{D}/{bam}.sorted.bam"],
        capture_output=True, text=True).stdout
    d = np.zeros(end - start + 1, dtype=np.int64)
    for line in out.splitlines():
        f = line.split("\t")
        d[int(f[1]) - start] = int(f[2])
    return d


for gene, chrom, s, e, strand, cls in GENES:
    exons = canonical_exons(gene, chrom, s, e)
    if not exons:
        print(f"{gene}: no exons found\n")
        continue
    if strand == "-":
        exons = exons[::-1]
    dg, dl = depth("graph", chrom, s, e), depth("linear", chrom, s, e)
    rows = []
    for k, (a, b) in enumerate(exons, 1):
        i, j = a - s, b - s + 1
        rows.append((k, b - a + 1, int(dg[i:j].sum()), int(dl[i:j].sum())))
    tot_gain = sum(g - l for _, _, g, l in rows)
    print(f"== {gene}  ({cls}, {len(exons)} exons, strand {strand})   "
          f"total exonic coverage gain {tot_gain:+,}")
    print(f"{'exon':>6}{'len':>7}{'graph':>13}{'linear':>13}{'gain':>12}"
          f"{'ratio':>8}{'%gain':>8}")
    for k, ln, g, l in rows:
        star = "  <- polymorphic" if k in POLY[cls] else ""
        r = g / l if l else float("nan")
        pct = 100 * (g - l) / tot_gain if tot_gain else float("nan")
        print(f"{k:>6}{ln:>7}{g:>13,}{l:>13,}{g-l:>+12,}{r:>8.2f}{pct:>7.1f}%{star}")
    pol = POLY[cls]
    if pol:
        pg = sum(g - l for k, _, g, l in rows if k in pol)
        pl = sum(ln for k, ln, _, _ in rows if k in pol)
        al = sum(ln for _, ln, _, _ in rows)
        print(f"       polymorphic exon(s) = {100*pl/al:.1f}% of exonic bases, "
              f"{100*pg/tot_gain:.1f}% of the gain")
    print()
