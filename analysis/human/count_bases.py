#!/usr/bin/env python3
"""Read `samtools mpileup` (no -f) on stdin, emit chrom pos A C G T depth."""
import sys

out = sys.stdout
for line in sys.stdin:
    f = line.split('\t')
    if len(f) < 5:
        continue
    bases = f[4]
    n = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    i = 0
    L = len(bases)
    while i < L:
        c = bases[i]
        if c == '^':          # read start, next char is mapping quality
            i += 2
            continue
        if c == '$':
            i += 1
            continue
        if c == '+' or c == '-':   # indel: +<len><seq>
            j = i + 1
            num = ''
            while j < L and bases[j].isdigit():
                num += bases[j]
                j += 1
            i = j + int(num) if num else j
            continue
        u = c.upper()
        if u in n:
            n[u] += 1
        i += 1
    tot = n['A'] + n['C'] + n['G'] + n['T']
    if tot == 0:
        continue
    out.write("%s\t%s\t%d\t%d\t%d\t%d\n" % (f[0], f[1], n['A'], n['C'], n['G'], n['T']))
