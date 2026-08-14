#!/usr/bin/env python3
"""Matched REF/ALT reads for measuring alignment-stage reference bias.

For each sampled SNV, emit two reads that are identical except for the base at
the variant: one carrying the reference allele, one the alternate. Both get the
SAME extra mismatches at the SAME offsets, so the pair differs in exactly one
position and any difference in whether they align is attributable to the variant
alone.

The extra mismatches are the point. A lone SNV is usually tolerated by the
aligner's mismatch budget, so a bias test without them measures almost nothing.
With the budget partly spent, knowing the variant is what decides whether the
alt-carrying read aligns at all -- which is the regime the graph index exists
for, and the one where dropped variants should hurt.
"""
import sys, random

fa, snp, out, n_sites, n_mm, seed = sys.argv[1:7]
n_sites, n_mm = int(n_sites), int(n_mm)
random.seed(int(seed))

seq = []
with open(fa) as f:
    f.readline()
    for line in f:
        seq.append(line.strip())
seq = ''.join(seq)

sites = []
with open(snp) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 5 or p[1] != 'single':
            continue
        pos, alt = int(p[3]), p[4]
        if len(alt) != 1 or alt not in 'ACGT':
            continue
        sites.append((p[0], pos, alt))

READ = 100
HALF = READ // 2
usable = []
for vid, pos, alt in sites:
    s = pos - HALF
    if s < 0 or s + READ > len(seq):
        continue
    w = seq[s:s + READ]
    if 'N' in w or seq[pos] == alt or seq[pos] not in 'ACGT':
        continue
    usable.append((vid, pos, alt, s, w))

random.shuffle(usable)
usable = usable[:n_sites]
usable.sort(key=lambda x: x[1])

nt = 'ACGT'
with open(out, 'w') as o:
    for vid, pos, alt, s, w in usable:
        vi = pos - s                                   # variant offset in the read
        offs = [i for i in range(READ) if abs(i - vi) > 3]
        mm = random.sample(offs, n_mm)                 # identical for both mates of the pair
        base = list(w)
        for i in mm:
            base[i] = random.choice([c for c in nt if c != base[i]])
        ref_read = base[:]                             # reference allele at vi
        alt_read = base[:]
        alt_read[vi] = alt
        o.write('>%s_%d_REF_%d\n%s\n' % (vid, pos, s, ''.join(ref_read)))
        o.write('>%s_%d_ALT_%d\n%s\n' % (vid, pos, s, ''.join(alt_read)))

print('sites usable=%d sampled=%d  reads=%d  mismatches/read=%d'
      % (len(usable), len(usable), 2 * len(usable), n_mm))
