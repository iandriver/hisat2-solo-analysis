#!/usr/bin/env python3
"""H2: does variant-aware recovery scale with distance from the reference?

Two donors, graph vs linear index, same reads and gene model each time.

  GM12878 = NA12878, CEU -- close to GRCh38
  GM18502 = NA18502, YRI -- the ancestry contrast

The hypothesis is that the YRI donor gains more, because GRCh38 is largely
European-derived and so represents that donor's haplotypes less well.

The controls matter as much as the HLA genes. B2M and ACTB should not move for
either donor; if they do, something other than reference divergence is going on
and the HLA numbers cannot be read as recovery.
"""

import sys
import os
from math import log

HLA = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G',
       'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DQA1', 'HLA-DQA2',
       'HLA-DQB1', 'HLA-DQB2', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DMA', 'HLA-DMB']
# Housekeeping genes that are NOT paralog-prone. RPL13A and EEF1A1 were tried
# and rejected as controls: both drop to 0.37 and 0.57 under the graph index,
# because processed pseudogenes of ribosomal-protein and translation-factor
# genes turn their reads into multimappers, which the NH:i:1 filter then
# excludes. That is the same mechanism H1 found behind OLFM3/RPSAP19 and the
# mouse three-way comparison found behind its ribosomal-protein disagreements --
# a real effect, but the opposite of a control.
CONTROLS = ['B2M', 'ACTB', 'GAPDH', 'TMSB4X']
PARALOG_PRONE = ['RPL13A', 'EEF1A1']


def load(path):
    d = {}
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) == 2:
                d[f[0]] = int(f[1])
    return d


def rate(path):
    """overall alignment rate from a hisat2 log"""
    if not os.path.exists(path):
        return None
    for ln in open(path):
        if 'overall alignment rate' in ln:
            return float(ln.split('%')[0].strip())
    return None


def main(L):
    donors = [('GM12878', 'CEU / European'), ('GM18502', 'YRI / Yoruba')]
    data = {}
    for s, _ in donors:
        g, l = load(f'{L}/gx_{s}_graph.tsv'), load(f'{L}/gx_{s}_linear.tsv')
        if g is None or l is None:
            print(f'missing counts for {s}; run h2_run.sh first')
            return 1
        data[s] = (g, l)

    print('=' * 72)
    print('H2 -- variant-aware recovery by donor ancestry')
    print('=' * 72)

    print('\n--- alignment rate (all reads) ---')
    print('%-10s %-16s %8s %8s %8s' % ('donor', 'population', 'graph', 'linear', 'delta'))
    for s, pop in donors:
        rg, rl = rate(f'{L}/aln_{s}_graph.log'), rate(f'{L}/aln_{s}_linear.log')
        if rg is not None and rl is not None:
            print('%-10s %-16s %7.2f%% %7.2f%% %+7.2f' % (s, pop, rg, rl, rg - rl))

    # The graph index loses uniquely-assigned reads overall, because alternate
    # alleles give more loci a chance to match and some reads that were unique
    # under the linear index become multimappers. So the right baseline for "no
    # effect" is this global ratio, not 1.00 -- otherwise every HLA gain is
    # understated.
    print('\n--- uniquely-assigned reads in genes (the normalising baseline) ---')
    print('%-10s %14s %14s %8s' % ('donor', 'graph', 'linear', 'ratio'))
    globalr = {}
    for s, pop in donors:
        g, l = data[s]
        tg, tl = sum(g.values()), sum(l.values())
        globalr[s] = tg / tl if tl else 1.0
        print('%-10s %14s %14s %7.4f' % (s, f'{tg:,}', f'{tl:,}', globalr[s]))

    print('\n--- paralog-prone genes, rejected as controls ---')
    print('%-10s %-10s %10s %10s %8s' % ('donor', 'gene', 'graph', 'linear', 'ratio'))
    for s, _ in donors:
        g, l = data[s]
        for gene in PARALOG_PRONE:
            a, b = g.get(gene, 0), l.get(gene, 0)
            if b >= 100:
                print('%-10s %-10s %10d %10d %8.3f' % (s, gene, a, b, a / b))

    print('\n--- controls (must be ~1.00 for the HLA numbers to mean anything) ---')
    print('%-10s %-10s %10s %10s %8s' % ('donor', 'gene', 'graph', 'linear', 'ratio'))
    ctrl_ratios = {s: [] for s, _ in donors}
    for s, _ in donors:
        g, l = data[s]
        for gene in CONTROLS:
            a, b = g.get(gene, 0), l.get(gene, 0)
            if b >= 100:
                r = a / b
                ctrl_ratios[s].append(r)
                print('%-10s %-10s %10d %10d %8.3f' % (s, gene, a, b, r))
    print()
    for s, _ in donors:
        rs = ctrl_ratios[s]
        if rs:
            print('  %s control mean ratio %.4f (n=%d)' % (s, sum(rs) / len(rs), len(rs)))

    print('\n--- HLA, graph vs linear, per donor ---')
    print('%-10s %10s %10s %8s | %10s %10s %8s' %
          ('gene', 'CEU grf', 'CEU lin', 'ratio', 'YRI grf', 'YRI lin', 'ratio'))
    g1, l1 = data['GM12878']
    g2, l2 = data['GM18502']
    rows = []
    for gene in HLA:
        a1, b1 = g1.get(gene, 0), l1.get(gene, 0)
        a2, b2 = g2.get(gene, 0), l2.get(gene, 0)
        r1 = a1 / b1 if b1 else float('nan')
        r2 = a2 / b2 if b2 else float('nan')
        rows.append((gene, a1, b1, r1, a2, b2, r2))
        print('%-10s %10d %10d %8.3f | %10d %10d %8.3f'
              % (gene, a1, b1, r1, a2, b2, r2))

    # The headline: does the YRI donor gain more than the CEU donor?
    print('\n--- the comparison H2 exists to make ---')
    both = [(gene, r1, r2) for gene, _, b1, r1, _, b2, r2 in rows
            if b1 >= 100 and b2 >= 100 and r1 == r1 and r2 == r2]
    if both:
        m1 = sum(r for _, r, _ in both) / len(both)
        m2 = sum(r for _, _, r in both) / len(both)
        print('  HLA genes compared (>=100 linear reads in both): %d' % len(both))
        print('  mean graph/linear ratio   CEU %.4f   YRI %.4f' % (m1, m2))
        c1 = sum(ctrl_ratios['GM12878']) / len(ctrl_ratios['GM12878']) if ctrl_ratios['GM12878'] else 1.0
        c2 = sum(ctrl_ratios['GM18502']) / len(ctrl_ratios['GM18502']) if ctrl_ratios['GM18502'] else 1.0
        print('  housekeeping-normalised   CEU %.4f   YRI %.4f' % (m1 / c1, m2 / c2))
        print('  transcriptome-normalised  CEU %.4f   YRI %.4f'
              % (m1 / globalr['GM12878'], m2 / globalr['GM18502']))
        # compare AFTER normalising: the two donors have different baselines
        # (0.9307 vs 0.9240), so a raw r2 > r1 test miscounts.
        n_yri = sum(1 for _, r1, r2 in both
                    if r2 / globalr['GM18502'] > r1 / globalr['GM12878'])
        print('  HLA genes where YRI gains more than CEU: %d / %d' % (n_yri, len(both)))
        print('\n  per-gene, both normalisers applied:')
        print('  %-10s %10s %10s' % ('gene', 'CEU norm', 'YRI norm'))
        for gene, r1, r2 in sorted(both, key=lambda x: -(x[2] / globalr['GM18502'])):
            print('  %-10s %10.3f %10.3f'
                  % (gene, r1 / globalr['GM12878'], r2 / globalr['GM18502']))
        print()
        if (m2 / globalr['GM18502']) > (m1 / globalr['GM12878']):
            print('  -> YRI gains more, consistent with recovery scaling with')
            print('     distance from a largely European-derived reference.')
        else:
            print('  -> YRI does NOT gain more. The ancestry-distance hypothesis is')
            print('     not supported by this pair.')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else '.'))
