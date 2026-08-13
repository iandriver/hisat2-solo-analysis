#!/usr/bin/env python3
"""Reduce a .snp/.haplotype pair to SNVs only.

The 64-bit whole-genome build OOM-killed at 368.8 GB during PathGraph
construction, and the 32-bit build had already refused 14.95M variants with
"exceeded integer bounds". The 2.3M indels are the expensive part: insertions
alone contribute 3,013,044 extra graph nodes on top of one node per variant.
Dropping them leaves 12.64M SNVs -- within 2.8% of the shipped grch38_snp's
12.3M, which is proof by existence that a set this size fits a 32-bit index.

Haplotype pruning is delegated to hisat2_filter_snps.filter_haplotypes rather
than reimplemented: it already drops haplotypes left with no surviving variant
and retightens left/right, and it is covered by the repo's test suite. Those
spans drive local graph construction, so a stale one is not cosmetic.
"""
import sys, gzip
sys.path.insert(0, '/Users/iandriver/Downloads/hisat2')
from hisat2_filter_snps import filter_haplotypes, variant_span

snp_in, snp_out, ht_in, ht_out = sys.argv[1:5]

opener = gzip.open if snp_in.endswith('.gz') else open
kept, snp_pos = set(), {}
n_in = n_out = 0
by_type = {}
with opener(snp_in, 'rt') as fi, open(snp_out, 'w') as fo:
    for line in fi:
        f = line.rstrip('\n').split('\t')
        if len(f) < 5:
            continue
        n_in += 1
        by_type[f[1]] = by_type.get(f[1], 0) + 1
        if f[1] != 'single':
            continue
        kept.add(f[0])
        snp_pos[f[0]] = variant_span(f[1], int(f[3]), f[4])
        fo.write(line)
        n_out += 1

opener = gzip.open if ht_in.endswith('.gz') else open
with opener(ht_in, 'rt') as hi, open(ht_out, 'w') as ho:
    h_in, h_out, h_trim = filter_haplotypes(hi, ho, kept, snp_pos)

print('variants  %d in -> %d out' % (n_in, n_out))
for k in sorted(by_type):
    print('   %-10s %d' % (k, by_type[k]))
print('haplotypes %d in -> %d out (%d trimmed)' % (h_in, h_out, h_trim))
