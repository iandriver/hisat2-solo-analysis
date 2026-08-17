#!/usr/bin/env python3
"""Concordance between two Solo.out-style matrix directories.

Equality is not the target and never will be: the two tools align to different
indexes with different scoring, so a difference is only interesting if it is
larger than what two runs of the *same* family disagree by. What this reports is
the shape of the disagreement -- which cells, which genes, and by how much.

Usage: compare.py <dirA> <dirB> [labelA] [labelB]
where each dir contains barcodes.tsv, features.tsv, matrix.mtx
"""

import sys
import os
from math import log, log1p, sqrt


def read_mtx(d):
    """-> (genes, barcodes, {(gi,bi): count})  with gene ids version-stripped."""
    def lines(p):
        with open(p) as fh:
            for ln in fh:
                ln = ln.rstrip('\n')
                if ln:
                    yield ln

    feats = []
    for ln in lines(os.path.join(d, 'features.tsv')):
        f = ln.split('\t')
        gid = f[0].split('.')[0]          # ENSG00000123456.7 -> ENSG00000123456
        name = f[1] if len(f) > 1 else f[0]
        feats.append((gid, name))
    bcs = [ln.split('\t')[0].split('-')[0] for ln in lines(os.path.join(d, 'barcodes.tsv'))]

    m = {}
    with open(os.path.join(d, 'matrix.mtx')) as fh:
        hdr_seen = False
        for ln in fh:
            if ln.startswith('%'):
                continue
            f = ln.split()
            if not hdr_seen:
                hdr_seen = True                     # rows cols nnz
                continue
            r, c, v = int(f[0]) - 1, int(f[1]) - 1, float(f[2])
            if v:
                m[(r, c)] = v
    return feats, bcs, m


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float('nan')
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sxx = syy = 0.0
    for x, y in zip(xs, ys):
        dx, dy = x - mx, y - my
        sxy += dx * dy
        sxx += dx * dx
        syy += dy * dy
    return sxy / sqrt(sxx * syy) if sxx > 0 and syy > 0 else float('nan')


def ranks(vs):
    order = sorted(range(len(vs)), key=lambda i: vs[i])
    r = [0.0] * len(vs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vs[order[j + 1]] == vs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(xs, ys):
    return pearson(ranks(xs), ranks(ys))


def main(da, db, la, lb):
    fa, ba, ma = read_mtx(da)
    fb, bb, mb = read_mtx(db)

    print("=" * 66)
    print("%-24s %-18s %-18s" % ("", la, lb))
    print("=" * 66)
    print("%-24s %18s %18s" % ("cells called", f"{len(ba):,}", f"{len(bb):,}"))
    print("%-24s %18s %18s" % ("features", f"{len(fa):,}", f"{len(fb):,}"))
    print("%-24s %18s %18s" % ("total UMIs", f"{int(sum(ma.values())):,}",
                               f"{int(sum(mb.values())):,}"))

    # ---- cell overlap ------------------------------------------------------
    sa, sb = set(ba), set(bb)
    inter, union = sa & sb, sa | sb
    print("\n--- cells ---")
    print("  shared            %s" % f"{len(inter):,}")
    print("  %s only%s%s" % (la, " " * max(1, 12 - len(la)), f"{len(sa - sb):,}"))
    print("  %s only%s%s" % (lb, " " * max(1, 12 - len(lb)), f"{len(sb - sa):,}"))
    print("  Jaccard           %.4f" % (len(inter) / len(union) if union else 0))
    if not inter:
        print("\nNo shared cells -- nothing further to compare.")
        return

    # ---- align gene and cell axes -----------------------------------------
    gia = {g: i for i, (g, _) in enumerate(fa)}
    gib = {g: i for i, (g, _) in enumerate(fb)}
    shared_genes = [g for g in gia if g in gib]
    name_of = {g: n for g, n in fa}
    bia = {b: i for i, b in enumerate(ba)}
    bib = {b: i for i, b in enumerate(bb)}
    cells = sorted(inter)
    print("\n--- features ---")
    print("  shared gene ids   %s" % f"{len(shared_genes):,}")

    # per-gene totals over shared cells
    ta = dict.fromkeys(shared_genes, 0.0)
    tb = dict.fromkeys(shared_genes, 0.0)
    cell_a = dict.fromkeys(cells, 0.0)
    cell_b = dict.fromkeys(cells, 0.0)
    colsa = {bia[b] for b in cells}
    colsb = {bib[b] for b in cells}
    rowa_gene = {gia[g]: g for g in shared_genes}
    rowb_gene = {gib[g]: g for g in shared_genes}
    inv_a = {bia[b]: b for b in cells}
    inv_b = {bib[b]: b for b in cells}

    for (r, c), v in ma.items():
        if c in colsa:
            g = rowa_gene.get(r)
            if g is not None:
                ta[g] += v
                cell_a[inv_a[c]] += v
    for (r, c), v in mb.items():
        if c in colsb:
            g = rowb_gene.get(r)
            if g is not None:
                tb[g] += v
                cell_b[inv_b[c]] += v

    expressed = [g for g in shared_genes if ta[g] + tb[g] > 0]
    xa = [ta[g] for g in expressed]
    xb = [tb[g] for g in expressed]
    print("  expressed in either %s" % f"{len(expressed):,}")

    print("\n--- agreement on shared cells ---")
    print("  per-gene totals   Pearson(log1p) %.4f   Spearman %.4f"
          % (pearson([log1p(v) for v in xa], [log1p(v) for v in xb]),
             spearman(xa, xb)))
    ca = [cell_a[b] for b in cells]
    cb = [cell_b[b] for b in cells]
    print("  per-cell totals   Pearson(log1p) %.4f   Spearman %.4f"
          % (pearson([log1p(v) for v in ca], [log1p(v) for v in cb]),
             spearman(ca, cb)))
    sa_umi, sb_umi = sum(xa), sum(xb)
    print("  UMIs on shared cells  %s vs %s  (%+.2f%%)"
          % (f"{int(sa_umi):,}", f"{int(sb_umi):,}",
             100.0 * (sa_umi - sb_umi) / sb_umi if sb_umi else 0))

    # ---- where they disagree ----------------------------------------------
    # Restricted to genes with real signal, so the ranking is not dominated by
    # 1-vs-0 counts that carry no information.
    MIN = 50
    cand = [g for g in expressed if ta[g] + tb[g] >= MIN]
    lr = sorted(((log((ta[g] + 1) / (tb[g] + 1), 2), g) for g in cand),
                reverse=True)
    print("\n--- largest disagreements (genes with >=%d UMIs total) ---" % MIN)
    print("  %d genes qualify" % len(cand))
    print("\n  %-14s %-12s %10s %10s %8s" % ("gene", "id", la, lb, "log2"))
    for v, g in lr[:15]:
        print("  %-14s %-12s %10d %10d %+8.2f"
              % (name_of.get(g, '?')[:14], g[:12], ta[g], tb[g], v))
    print("  %s" % ("." * 60))
    for v, g in lr[-15:]:
        print("  %-14s %-12s %10d %10d %+8.2f"
              % (name_of.get(g, '?')[:14], g[:12], ta[g], tb[g], v))

    within2x = sum(1 for v, _ in lr if abs(v) < 1.0)
    print("\n  within 2x         %d / %d  (%.1f%%)"
          % (within2x, len(lr), 100.0 * within2x / len(lr) if lr else 0))

    # HLA specifically -- the locus the variant-aware work is about
    hla = [(v, g) for v, g in lr if name_of.get(g, '').startswith('HLA-')]
    if hla:
        print("\n--- HLA genes ---")
        print("  %-14s %10s %10s %8s" % ("gene", la, lb, "log2"))
        for v, g in sorted(hla, reverse=True):
            print("  %-14s %10d %10d %+8.2f" % (name_of.get(g, '?'), ta[g], tb[g], v))


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2],
         sys.argv[3] if len(sys.argv) > 3 else 'A',
         sys.argv[4] if len(sys.argv) > 4 else 'B')
