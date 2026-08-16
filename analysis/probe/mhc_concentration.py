#!/usr/bin/env python3
"""Where in the genome would hisat2-build drop variants?

The probe reports a genome-wide retention figure. That average hides the thing
that matters here: the loss is not spread evenly, it concentrates wherever
variant density is highest, and for a human variant set that is the MHC -- the
same region where variant-aware alignment showed its largest advantage.

usage: mhc_concentration.py <genome.snp>
"""
import sys

INTERVAL = (1 << 16) - (1 << 13) - 1024        # hier_idx_common.h
CAPACITY = 384.0                                # 420 / 1.09 haplotypes-per-variant
MHC = (28510120, 33480577)                      # extended MHC, GRCh38

# Coordinates and the gains measured in the T1 reference-bias work.
GENES = {
    'HLA-DQA1': (32628179, 32647062, 32.10),
    'HLA-DQB1': (32659467, 32668383, 3.00),
    'HLA-DRB1': (32578775, 32589848, 1.54),
    'HLA-C':    (31268749, 31272130, 1.51),
    'HLA-A':    (29942532, 29945870, 1.00),
    'HLA-B':    (31353875, 31357188, 1.00),
}


def load(path):
    w = {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 5:
            continue
        k = (f[2], int(f[3]) // INTERVAL)
        w[k] = w.get(k, 0) + 1
    return w


def region(w, chrom, lo, hi):
    ks = [(c, i) for (c, i) in w
          if c == chrom and i * INTERVAL < hi and (i + 1) * INTERVAL > lo]
    n = sum(w[k] for k in ks)
    over = [k for k in ks if w[k] > CAPACITY]
    return len(ks), n, len(over), sum(w[k] for k in over)


def main(snp):
    w = load(snp)
    tot_w, tot_v = len(w), sum(w.values())
    over = {k: v for k, v in w.items() if v > CAPACITY}
    print("genome-wide   %6d windows  %10d variants  %5d over (%.1f%%)  %.1f%% at risk"
          % (tot_w, tot_v, len(over), 100 * len(over) / tot_w,
             100 * sum(over.values()) / tot_v))
    nk, n, no, at = region(w, 'chr6', *MHC)
    print("extended MHC  %6d windows  %10d variants  %5d over (%.1f%%)  %.1f%% at risk"
          % (nk, n, no, 100 * no / nk, 100 * at / n))
    print("\nMHC density is %.1fx the genome-wide average"
          % ((n / nk) / (tot_v / tot_w)))

    print("\n%-10s %8s %9s %12s %9s" % ("gene", "windows", "variants", "over budget", "T1 gain"))
    for g, (lo, hi, gain) in sorted(GENES.items(), key=lambda x: -x[1][2]):
        nk, n, no, _ = region(w, 'chr6', lo, hi)
        print("%-10s %8d %9d %12s %8.2fx" % (g, nk, n, "%d/%d" % (no, nk), gain))

    top = sorted(over.items(), key=lambda x: -x[1])[:100]
    inm = sum(1 for (c, i), _ in top
              if c == 'chr6' and i * INTERVAL < MHC[1] and (i + 1) * INTERVAL > MHC[0])
    print("\n%d of the 100 densest over-budget windows are in the MHC,"
          " which is %.2f%% of all windows" % (inm, 100 * 89 / tot_w))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
