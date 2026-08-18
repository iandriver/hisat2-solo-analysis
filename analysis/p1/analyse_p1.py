#!/usr/bin/env python3
"""P1 -- alternate-allele fraction at known heterozygous sites, three arms.

Unbiased alignment puts the alt fraction at 0.5. A linear reference pushes it
below, because a read carrying the non-reference base pays a mismatch penalty
and sometimes loses to a paralogue or falls below threshold.

Everything is PAIRED: only sites covered at sufficient depth in all three arms
are used, so a difference between arms cannot come from which sites happened to
be covered. Reported three ways, because they answer different questions:

  pooled      sum(alt)/sum(alt+ref) -- the population-level ratio, dominated by
              deep sites; this is the number a downstream ASE analysis feels
  per-site    mean over sites of alt/(alt+ref) -- each site counts once
  paired      per-site difference between arms, with a sign test; this is the
              one with power, since it removes site-to-site variation

The presence split is the honest core. Both donors are IN the 1000 Genomes panel
this index was built from, so a bare "graph beats linear" number is partly
circular. Splitting sites on whether the donor's alt allele is actually carried
by each index separates "the variant is in the index, so bias is removed" from
"the variant is absent, so bias remains" -- and the second group is the
uncircular control that says what happens to a donor the panel does not cover.
"""
import sys, math, os, gzip

MINDEP = int(os.environ.get('MINDEP', '20'))
L = '/Users/iandriver/Downloads/p1_lcl/out'
T = '/Users/iandriver/Downloads/p1_truth'
ARMS = ['linear', 'snp_invented', 'snp_real']

def load_sites(path):
    d = {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) != 4:
            continue          # a stream still being written leaves a partial last line
        c, p, r, a = f
        d[(c.replace('chr', ''), int(p))] = (r, a)
    return d

def load_counts(path):
    d = {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        c = f[0].replace('chr', '')
        d[(c, int(f[1]))] = dict(zip('ACGT', map(int, f[2:6])))
    return d

def load_index_variants(path):
    """chrom,pos,alt for single SNVs an index carries (hisat2-inspect --snp)."""
    s = set()
    op = gzip.open if path.endswith('.gz') else open
    for line in op(path, 'rt'):
        f = line.rstrip('\n').split('\t')
        if len(f) < 5 or f[1] != 'single':
            continue
        c = f[2].split()[0].replace('chr', '')
        s.add((c, int(f[3]) + 1, f[4]))   # stored 0-based -> 1-based
    return s

def wilson(k, n):
    if n == 0: return (0.0, 0.0)
    z = 1.96; p = k / n
    d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return (c - h, c + h)

def sign_test(diffs):
    """Returns (up, down, z, p). z is reported alongside p because at these
    counts p underflows to 0.0, and printing "p = 0" would be a false precision."""
    pos = sum(1 for d in diffs if d > 0); neg = sum(1 for d in diffs if d < 0)
    n = pos + neg
    if n == 0: return pos, neg, 0.0, 1.0
    z = (pos - n/2) / math.sqrt(n/4)
    p = math.erfc(abs(z)/math.sqrt(2))
    return pos, neg, z, p

def report(donor, sites_path, present):
    sites = load_sites(sites_path)
    cnt = {}
    for a in ARMS:
        p = f'{L}/cnt_{donor}_{a}.tsv'
        if not os.path.exists(p):
            print(f'  {donor}: missing {p}'); return
        cnt[a] = load_counts(p)
    common = [k for k in sites
              if all(k in cnt[a] and sum(cnt[a][k].values()) >= MINDEP for a in ARMS)]
    print(f'\n=== {donor} ===')
    print(f'  {len(sites)} known het sites; {len(common)} covered at depth >= {MINDEP} in all three arms')
    if not common: return

    def frac(a, keys):
        pooled_alt = pooled_tot = 0; per = []
        for k in keys:
            r, al = sites[k]
            n = cnt[a][k]; ra, aa = n[r], n[al]
            if ra + aa == 0: continue
            pooled_alt += aa; pooled_tot += ra + aa
            per.append(aa / (ra + aa))
        return pooled_alt, pooled_tot, per

    groups = [('all sites', common)]
    for arm in ('snp_invented', 'snp_real'):
        if arm in present:
            inx = [k for k in common if (k[0], k[1], sites[k][1]) in present[arm]]
            out = [k for k in common if (k[0], k[1], sites[k][1]) not in present[arm]]
            groups.append((f'variant IN {arm} (n={len(inx)})', inx))
            groups.append((f'variant NOT in {arm} (n={len(out)})', out))

    for label, keys in groups:
        if len(keys) < 50:
            print(f'\n  -- {label}: too few sites ({len(keys)})'); continue
        print(f'\n  -- {label}')
        print(f'     {"arm":<14} {"pooled alt frac":>18} {"95% CI":>18} {"per-site mean":>15}')
        base = None
        for a in ARMS:
            pa, pt, per = frac(a, keys)
            lo, hi = wilson(pa, pt)
            m = sum(per)/len(per) if per else float('nan')
            print(f'     {a:<14} {pa/pt:>18.5f} {f"{lo:.5f}-{hi:.5f}":>18} {m:>15.5f}')
        for a in ('snp_invented', 'snp_real'):
            diffs = []
            for k in keys:
                r, al = sites[k]
                x, y = cnt['linear'][k], cnt[a][k]
                if x[r]+x[al] == 0 or y[r]+y[al] == 0: continue
                diffs.append(y[al]/(y[r]+y[al]) - x[al]/(x[r]+x[al]))
            pos, neg, z, p = sign_test(diffs)
            md = sum(diffs)/len(diffs) if diffs else 0.0
            ps = f'p = {p:.3g}' if p > 0 else 'p < 1e-300'
            print(f'     {a} vs linear: mean per-site change {md:+.5f}, '
                  f'{pos} up / {neg} down, sign test z = {z:.1f}, {ps}')

if __name__ == '__main__':
    present = {}
    for arm, path in (('snp_invented', f'{T}/vars_grch38_snp.tsv.gz'),
                      ('snp_real',     f'{T}/vars_wg64.tsv.gz')):
        if os.path.exists(path):
            present[arm] = load_index_variants(path)
            print(f'{arm}: {len(present[arm])} single-SNV variants carried')
    report('GM12878', f'{T}/na12878.het.tsv', present)
    if os.path.exists(f'{T}/na18502.het.tsv'):
        report('GM18502', f'{T}/na18502.het.tsv', present)
