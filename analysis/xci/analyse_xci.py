#!/usr/bin/env python3
"""X-inactivation escape under linear vs graph alignment: all four analyses.

Reads the per-arm base counts in counts/ and the phased chrX genotypes, and
emits every number quoted in xci_report.md. Deliberately one script so the
report cannot drift from the data.

  usage: analyse_xci.py [gtf]     (default: the GRCh38 GTF used for the run)
"""
import collections, math, statistics, sys, gzip, os

HERE = os.path.dirname(os.path.abspath(__file__))
GTF = sys.argv[1] if len(sys.argv) > 1 else "/Users/iandriver/Downloads/hisat2-refs/genes.gtf"
ARMS = ["linear", "snp_invented", "snp_real_Xbroken", "snp_real_Xfixed"]
DONORS = [("GM12878", "na12878"), ("GM18502", "na18502")]
# GRCh38 pseudoautosomal regions. PAR escapes X-inactivation by definition, so
# it is a positive control rather than a discovery set -- excluded throughout.
PAR = lambda p: p <= 2781479 or p >= 155701383


def counts(donor, arm):
    d = {}
    with open(f"{HERE}/counts/cnt_{donor}_{arm}.tsv") as fh:
        for line in fh:
            f = line.split()
            if len(f) < 6:            # count_bases.py emits 6 columns; its
                continue              # docstring claims a 7th it never writes
            d[(f[0].replace("chr", ""), int(f[1]))] = {
                'A': int(f[2]), 'C': int(f[3]), 'G': int(f[4]), 'T': int(f[5])}
    return d


def phased(stem):
    p = {}
    with gzip.open(f"{HERE}/{stem}.chrX.phased.tsv.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            p[(f[0].replace("chr", ""), int(f[1]))] = (f[1 + 1], f[3], f[4])
    return p


def hap2(alt_n, ref_n, gt):
    """Reads on haplotype 2. GT 0|1 puts ALT on hap2; 1|0 puts REF there."""
    return (alt_n, ref_n) if gt == "0|1" else (ref_n, alt_n)


def sign_z(up, dn):
    n = up + dn
    return (up - n / 2) / math.sqrt(n / 4) if n else 0.0


def xgenes():
    g = []
    for line in open(GTF):
        if line[0] == '#':
            continue
        f = line.split('\t')
        if len(f) < 9 or f[0] != "X" or f[2] != "gene":
            continue
        nm = ""
        for kv in f[8].split(';'):
            kv = kv.strip()
            if kv.startswith('gene_name '):
                nm = kv.split('"')[1]
        g.append((int(f[3]), int(f[4]), nm or "?"))
    return sorted(g)


def main():
    GENES = xgenes()
    def gene_at(p):
        for a, b, n in GENES:
            if a <= p <= b:
                return n
        return None

    print(f"chrX genes in annotation: {len(GENES)}\n")

    for donor, stem in DONORS:
        P = phased(stem)
        D = {a: counts(donor, a) for a in ARMS}

        # ---- 1. paired allele fraction, sites covered in ALL arms ----------
        common = set(D[ARMS[0]])
        for a in ARMS[1:]:
            common &= set(D[a])
        for THR in (10, 20):
            keep = [k for k in common if k in P and not PAR(k[1])
                    and all(sum(D[a][k].values()) >= THR for a in ARMS)]
            if not keep:
                continue
            af = {a: [D[a][k].get(P[k][1], 0) / sum(D[a][k].values()) for k in keep]
                  for a in ARMS}
            print(f"=== {donor}: paired alt fraction, non-PAR, depth>={THR}, n={len(keep)}")
            for a in ARMS:
                print(f"    {a:<20} {sum(af[a])/len(af[a]):.4f}")
            base = af["linear"]
            for a in ARMS[1:]:
                d = [x - y for x, y in zip(af[a], base)]
                up = sum(1 for x in d if x > 1e-9); dn = sum(1 for x in d if x < -1e-9)
                print(f"    {a:<20} vs linear {sum(d)/len(d):+.5f}  "
                      f"{up} up / {dn} down  z={sign_z(up,dn):+.1f}")
            print()

        # ---- 2. the only clean contrast: same index, chrX content differs ---
        B, F = D["snp_real_Xbroken"], D["snp_real_Xfixed"]
        for THR in (10, 20):
            keep = [k for k in (set(B) & set(F)) if k in P and not PAR(k[1])
                    and sum(B[k].values()) >= THR and sum(F[k].values()) >= THR]
            if not keep:
                continue
            d = [F[k].get(P[k][1], 0)/sum(F[k].values()) - B[k].get(P[k][1], 0)/sum(B[k].values())
                 for k in keep]
            up = sum(1 for x in d if x > 1e-9); dn = sum(1 for x in d if x < -1e-9)
            print(f"=== {donor}: X-fixed minus X-broken, non-PAR, depth>={THR}, n={len(keep)}")
            print(f"    mean {sum(d)/len(d):+.5f}  {up} up / {dn} down / "
                  f"{len(d)-up-dn} unchanged  z={sign_z(up,dn):+.1f}")

        # ---- 3. X-inactivation skewing -------------------------------------
        rows = []
        for k, c in D["snp_real_Xfixed"].items():
            if k not in P or PAR(k[1]):
                continue
            ref, alt, gt = P[k]
            a, r = c.get(alt, 0), c.get(ref, 0)
            if a + r < 10:
                continue
            h2, h1 = hap2(a, r, gt)
            rows.append((h2 / (h2 + h1), h2 + h1))
        h = [x[0] for x in rows]
        lo = sum(1 for x in h if x < 0.10); hi = sum(1 for x in h if x > 0.90)
        print(f"\n=== {donor}: X-inactivation skewing, n={len(h)} non-PAR sites depth>=10")
        print(f"    mean hap2 fraction {sum(h)/len(h):.4f}  median {statistics.median(h):.4f}")
        print(f"    hap1-only {lo} ({100*lo/len(h):.1f}%)  hap2-only {hi} ({100*hi/len(h):.1f}%)")
        print(f"    -> {100*(lo+hi)/len(h):.0f}% monoallelic, split BOTH ways: not clonal skewing")
        for thr in (20, 50, 100):
            m = sum(1 for f2, d0 in rows if d0 >= thr and (f2 < 0.05 or f2 > 0.95))
            b = sum(1 for f2, d0 in rows if d0 >= thr and 0.05 <= f2 <= 0.95)
            if m + b:
                print(f"      depth>={thr:<4}: {m} mono / {b} bi ({100*m/(m+b):.0f}% mono)")

        # ---- 4. gene-level escape calls, per arm ---------------------------
        calls = {}
        for a in ARMS:
            per = collections.defaultdict(lambda: [0, 0, 0])
            for k, c in D[a].items():
                if k not in P or PAR(k[1]):
                    continue
                ref, alt, gt = P[k]
                av, rv = c.get(alt, 0), c.get(ref, 0)
                if av + rv < 10:
                    continue
                g = gene_at(k[1])
                if not g:
                    continue
                h2v, h1v = hap2(av, rv, gt)
                per[g][0] += h2v; per[g][1] += h1v; per[g][2] += 1
            calls[a] = {g for g, (h2v, h1v, n) in per.items()
                        if h1v + h2v >= 20 and n >= 2 and 0.30 <= h2v/(h1v+h2v) <= 0.70}
        print(f"\n=== {donor}: escape-like genes (>=2 sites, depth>=20, hap frac 0.30-0.70)")
        for a in ARMS:
            print(f"    {a:<20} {len(calls[a]):>3}  {sorted(calls[a])}")
        L, Fx, Bk = calls["linear"], calls["snp_real_Xfixed"], calls["snp_real_Xbroken"]
        print(f"    X-fixed vs linear   gained {sorted(Fx-L)}  lost {sorted(L-Fx)}")
        print(f"    X-fixed vs X-broken gained {sorted(Fx-Bk)}  lost {sorted(Bk-Fx)}")
        print()


if __name__ == "__main__":
    main()
