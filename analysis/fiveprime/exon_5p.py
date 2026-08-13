#!/usr/bin/env python3
"""5' vs 3': does the graph's HLA gain move to the polymorphic exons?

The 3' result could not test the plan's original hypothesis, because 3' chemistry
puts essentially all coverage in the terminal exon. Class I polymorphism sits in
exons 2 and 3 (the peptide-binding groove), near the 5' end. This runs the same
per-exon comparison on both chemistries so the difference is visible directly.
"""
import subprocess, re, sys
import numpy as np

H5 = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/fiveprime"
H3 = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
GTF = f"{H3}/genes.gtf"

BAMS = {"5'": (f"{H5}/5p_graph.bam", f"{H5}/5p_linear.bam"),
        "3'": (f"{H3}/graph.sorted.bam", f"{H3}/linear.sorted.bam")}

GENES = [
    ("HLA-A",    "6", 29941260, 29949572, "+", "classI"),
    ("HLA-B",    "6", 31353872, 31367067, "-", "classI"),
    ("HLA-C",    "6", 31268749, 31272130, "-", "classI"),
    ("HLA-DQA1", "6", 32628179, 32647062, "+", "classII_a"),
    ("HLA-DQB1", "6", 32659467, 32668383, "-", "classII_b"),
    ("HLA-DRB1", "6", 32577902, 32589848, "-", "classII_b"),
    ("HLA-DRA",  "6", 32439878, 32445046, "+", "invariant"),
    ("B2M",      "15", 44711487, 44718877, "+", "control"),
]
POLY = {"classI": {2, 3}, "classII_a": {2}, "classII_b": {2},
        "invariant": set(), "control": set()}


def canonical_exons(gene, chrom, start, end):
    awk = ('$1=="%s" && $3=="exon" && $4>=%d && $5<=%d && /gene_name "%s"/'
           % (chrom, start, end, gene))
    out = subprocess.run(["awk", "-F", "\t", awk, GTF],
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
         "-r", f"{chrom}:{start}-{end}", bam],
        capture_output=True, text=True).stdout
    d = np.zeros(end - start + 1, dtype=np.int64)
    for line in out.splitlines():
        f = line.split("\t")
        d[int(f[1]) - start] = int(f[2])
    return d


rows = {}
for chem, (bg, bl) in BAMS.items():
    for gene, chrom, s, e, strand, cls in GENES:
        exons = canonical_exons(gene, chrom, s, e)
        if not exons:
            continue
        if strand == "-":
            exons = exons[::-1]                 # number 5'->3'
        dg, dl = depth(bg, chrom, s, e), depth(bl, chrom, s, e)
        per = []
        for k, (a, b) in enumerate(exons, 1):
            i, j = a - s, b - s + 1
            per.append((k, b - a + 1, int(dg[i:j].sum()), int(dl[i:j].sum())))
        rows[(chem, gene)] = (cls, per)

# ---- 1. where does each chemistry put its reads? -------------------------
print("=== fraction of exonic coverage falling in exons 2-3 (class I) or "
      "exon 2 (class II) ===")
print(f"{'gene':<10}{'class':<11}{'3-prime':>10}{'5-prime':>10}")
for gene, chrom, s, e, strand, cls in GENES:
    if not POLY[cls]:
        continue
    out = []
    for chem in ("3'", "5'"):
        if (chem, gene) not in rows:
            out.append(float("nan")); continue
        _, per = rows[(chem, gene)]
        tot = sum(g for _, _, g, _ in per)
        pol = sum(g for k, _, g, _ in per if k in POLY[cls])
        out.append(100 * pol / tot if tot else float("nan"))
    print(f"{gene:<10}{cls:<11}{out[0]:>9.1f}%{out[1]:>9.1f}%")

# ---- 2. the test: graph/linear inside the polymorphic exons --------------
print("\n=== graph / linear coverage ratio, polymorphic exons vs the rest ===")
print(f"{'gene':<10}{'chem':<6}{'poly exons':>26}{'other exons':>26}")
print(f"{'':<16}{'graph':>11}{'linear':>8}{'ratio':>7}{'graph':>11}{'linear':>8}{'ratio':>7}")
for gene, chrom, s, e, strand, cls in GENES:
    for chem in ("3'", "5'"):
        if (chem, gene) not in rows:
            continue
        _, per = rows[(chem, gene)]
        pg = sum(g for k, _, g, _ in per if k in POLY[cls])
        pl = sum(l for k, _, _, l in per if k in POLY[cls])
        og = sum(g for k, _, g, _ in per if k not in POLY[cls])
        ol = sum(l for k, _, _, l in per if k not in POLY[cls])
        pr = pg / pl if pl else float("nan")
        orr = og / ol if ol else float("nan")
        print(f"{gene if chem=='3' + chr(39) else '':<10}{chem:<6}"
              f"{pg:>11,}{pl:>8,}{pr:>7.2f}{og:>11,}{ol:>8,}{orr:>7.2f}"
              + ("" if POLY[cls] else "   (no polymorphic exon defined)"))

# ---- 3. full per-exon profile, 5' only -----------------------------------
print("\n=== 5' per-exon detail ===")
for gene, chrom, s, e, strand, cls in GENES:
    if ("5'", gene) not in rows:
        continue
    _, per = rows[("5'", gene)]
    print(f"-- {gene} ({cls}, {strand})")
    for k, ln, g, l in per:
        star = " *polymorphic" if k in POLY[cls] else ""
        r = g / l if l else float("nan")
        print(f"     exon {k:<3}{ln:>6} bp{g:>12,}{l:>10,}{r:>8.2f}{star}")

import json
json.dump({f"{c}|{g}": v[1] for (c, g), v in rows.items()},
          open(f"{H5}/exon_profiles.json", "w"))
