#!/usr/bin/env python3
"""Alignment-stage reference bias across indexes, plus a paired site-level test.

ALT fraction is the share of aligned reads at a site that carry the alternate
allele. Because every site contributes one REF and one ALT read differing in
exactly one base, 0.5 is unbiased by construction and any departure is the
aligner's.

The paired McNemar test is the sensitive one: comparing two indexes on the same
sites, only the sites where they disagree carry information, so it detects
differences far below what the pooled fraction's confidence interval suggests.

usage: analyze_bias.py <label>=<sam> [<label>=<sam> ...]
"""
import sys, collections, math

def load(path):
    """site -> {'REF','ALT'} for reads that aligned at their true locus."""
    hits = collections.defaultdict(set)
    for line in open(path):
        if line[0] == '@':
            continue
        f = line.split('\t')
        flag = int(f[1])
        if flag & 4 or flag & 256 or flag & 2048:      # unmapped / secondary / supplementary
            continue
        vid, site, typ, start = f[0].rsplit('_', 3)
        if int(f[3]) != int(start) + 1:                # wrong locus does not count
            continue
        hits[(vid, site)].add(typ)
    return hits

args = [a.split('=', 1) for a in sys.argv[1:]]
H = {lab: load(p) for lab, p in args}

print("%-12s %10s %10s %11s %11s" % ("index", "REF aln", "ALT aln", "ALT frac", "95% CI"))
for lab, _ in args:
    h = H[lab]
    a = sum(1 for s in h.values() if 'ALT' in s)
    r = sum(1 for s in h.values() if 'REF' in s)
    n = a + r
    frac = a / n if n else float('nan')
    ci = 1.96 * math.sqrt(frac * (1 - frac) / n) if n else float('nan')
    print("%-12s %10d %10d %11.5f %11.5f" % (lab, r, a, frac, ci))

print("\npaired (McNemar) on ALT alignment, same sites:")
labs = [l for l, _ in args]
for i in range(len(labs)):
    for j in range(i + 1, len(labs)):
        x, y = labs[i], labs[j]
        hx, hy = H[x], H[y]
        keys = set(hx) | set(hy)
        lost = sum(1 for k in keys if 'ALT' in hx.get(k, ()) and 'ALT' not in hy.get(k, ()))
        gain = sum(1 for k in keys if 'ALT' not in hx.get(k, ()) and 'ALT' in hy.get(k, ()))
        chi = (abs(lost - gain) - 1) ** 2 / (lost + gain) if lost + gain else 0.0
        sig = "p<0.001" if chi > 10.83 else ("p<0.05" if chi > 3.84 else "n.s.")
        print("  %-10s -> %-10s lost=%-6d gained=%-6d net=%+6d chi2=%7.2f  %s"
              % (x, y, lost, gain, gain - lost, chi, sig))
