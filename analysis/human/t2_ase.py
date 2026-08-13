#!/usr/bin/env python3
"""T2: does linear-reference bias produce false allele-specific-expression calls?

Read-level REF/ALT counts at dbSNP SNV positions from the graph and linear
alignments of the same reads, then a per-site binomial test for allelic
imbalance under each aligner, BH-corrected.

The depth-matched arm is the important one: linear recovers fewer reads at
variant sites, so it has less power; downsampling both aligners to the same
per-site depth removes that confound.
"""
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from scipy.stats import false_discovery_control as bh

D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
MINDEPTH = 20
AFLO, AFHI = 0.15, 0.85
FDR = 0.05
BASES = list("ACGT")

print("loading snp table...", flush=True)
snp = pd.read_csv(f"{D}/snp_pos.tsv", sep="\t", header=None,
                  names=["rsid", "chrom", "pos0", "alt"], dtype={"chrom": str})
snp["pos"] = snp.pos0 + 1
snp = snp[snp.alt.isin(BASES)]
snp = snp.drop_duplicates(subset=["chrom", "pos"], keep=False)   # unambiguous only

def load_counts(tag):
    c = pd.read_csv(f"{D}/cnt_{tag}.tsv", sep="\t", header=None,
                    names=["chrom", "pos", "A", "C", "G", "T"], dtype={"chrom": str})
    return c

print("loading pileups...", flush=True)
cg = load_counts("graph")
cl = load_counts("linear")
print(f"  graph {len(cg):,} covered sites, linear {len(cl):,}")

m = snp.merge(cg, on=["chrom", "pos"], suffixes=("", "")) \
       .merge(cl, on=["chrom", "pos"], suffixes=("_g", "_l"))
print(f"  {len(m):,} SNV positions covered in both\n")

Ag = m[["A_g", "C_g", "G_g", "T_g"]].values.astype(np.int64)
Al = m[["A_l", "C_l", "G_l", "T_l"]].values.astype(np.int64)
alt_i = m.alt.map({b: i for i, b in enumerate(BASES)}).values

pool = Ag + Al
rows = np.arange(len(m))
alt_pool = pool[rows, alt_i]
tmp = pool.copy()
tmp[rows, alt_i] = -1
ref_i = tmp.argmax(axis=1)                 # REF = commonest non-ALT base, pooled
ref_pool = pool[rows, ref_i]

# biallelic: REF+ALT must dominate
biallelic = (ref_pool + alt_pool) >= 0.95 * pool.sum(axis=1)

ref_g, alt_g = Ag[rows, ref_i], Ag[rows, alt_i]
ref_l, alt_l = Al[rows, ref_i], Al[rows, alt_i]
dep_g, dep_l = ref_g + alt_g, ref_l + alt_l

with np.errstate(invalid="ignore", divide="ignore"):
    af_g = alt_g / dep_g
    af_l = alt_l / dep_l

df = pd.DataFrame(dict(rsid=m.rsid, chrom=m.chrom, pos=m.pos,
                       ref_g=ref_g, alt_g=alt_g, dep_g=dep_g, af_g=af_g,
                       ref_l=ref_l, alt_l=alt_l, dep_l=dep_l, af_l=af_l,
                       biallelic=biallelic))


def binom_p(k, n):
    """Two-sided binomial p vs 0.5, vectorised over unique (k,n) pairs."""
    key = pd.MultiIndex.from_arrays([k, n])
    uniq = pd.MultiIndex.from_frame(pd.DataFrame({"k": k, "n": n}).drop_duplicates())
    lut = {(int(a), int(b)): binomtest(int(a), int(b), 0.5).pvalue
           for a, b in zip(uniq.get_level_values(0), uniq.get_level_values(1))}
    return np.array([lut[(int(a), int(b))] for a, b in zip(k, n)])


def analyse(sel, label, kg, ng, kl, nl):
    n = sel.sum()
    print(f"\n{'='*72}\n{label}\n  {n:,} sites\n{'='*72}")
    if n == 0:
        return None
    pg = binom_p(kg[sel], ng[sel])
    pl = binom_p(kl[sel], nl[sel])
    qg, ql = bh(pg), bh(pl)
    sg, sl_ = qg < FDR, ql < FDR
    afg = kg[sel] / ng[sel]
    afl = kl[sel] / nl[sel]

    print(f"  median depth   graph {np.median(ng[sel]):.0f}   linear {np.median(nl[sel]):.0f}")
    print(f"  median ALT frac graph {np.median(afg):.4f}   linear {np.median(afl):.4f}")
    print(f"\n  imbalanced at FDR<{FDR}:  graph {sg.sum():,} ({100*sg.mean():.2f}%)"
          f"   linear {sl_.sum():,} ({100*sl_.mean():.2f}%)")
    print(f"\n  {'':<22}{'graph bal':>12}{'graph imbal':>13}")
    print(f"  {'linear balanced':<22}{(~sl_&~sg).sum():>12,}{(~sl_&sg).sum():>13,}")
    print(f"  {'linear imbalanced':<22}{(sl_&~sg).sum():>12,}{(sl_&sg).sum():>13,}")

    lonly = sl_ & ~sg
    gonly = sg & ~sl_
    print(f"\n  linear-only calls: {lonly.sum():,}"
          f"   ref-skewed {int((afl[lonly]<0.5).sum()):,}"
          f"   alt-skewed {int((afl[lonly]>0.5).sum()):,}")
    print(f"  graph-only  calls: {gonly.sum():,}"
          f"   ref-skewed {int((afg[gonly]<0.5).sum()):,}"
          f"   alt-skewed {int((afg[gonly]>0.5).sum()):,}")
    if lonly.sum():
        print(f"\n  at linear-only sites: median AF linear {np.median(afl[lonly]):.4f}"
              f"   graph {np.median(afg[lonly]):.4f}")
    if gonly.sum():
        print(f"  at graph-only  sites: median AF linear {np.median(afl[gonly]):.4f}"
              f"   graph {np.median(afg[gonly]):.4f}")
    return dict(n=n, g=int(sg.sum()), l=int(sl_.sum()),
                lonly=int(lonly.sum()), gonly=int(gonly.sum()),
                lonly_ref=int((afl[lonly] < 0.5).sum()),
                lonly_alt=int((afl[lonly] > 0.5).sum()),
                gonly_ref=int((afg[gonly] < 0.5).sum()),
                gonly_alt=int((afg[gonly] > 0.5).sum()))


# ---- site selection ------------------------------------------------------
het_l = df.biallelic & (df.dep_l >= MINDEPTH) & (df.af_l >= AFLO) & (df.af_l <= AFHI)
het_g = df.biallelic & (df.dep_g >= MINDEPTH) & (df.af_g >= AFLO) & (df.af_g <= AFHI)
both = het_l & het_g & (df.dep_g >= MINDEPTH) & (df.dep_l >= MINDEPTH)
het_l = het_l & (df.dep_g >= MINDEPTH)
het_g = het_g & (df.dep_l >= MINDEPTH)

res = {}
res["linear-selected"] = analyse(
    het_l.values, "A. Sites called heterozygous BY THE LINEAR ALIGNER "
                  "(conservative: favours linear)\n   raw depth, unmatched power",
    alt_g.astype(np.int64), dep_g.astype(np.int64),
    alt_l.astype(np.int64), dep_l.astype(np.int64))

res["graph-selected"] = analyse(
    het_g.values, "B. Sites called heterozygous BY THE GRAPH\n   raw depth, unmatched power",
    alt_g.astype(np.int64), dep_g.astype(np.int64),
    alt_l.astype(np.int64), dep_l.astype(np.int64))

# ---- depth-matched: equalise power --------------------------------------
print(f"\n\n{'#'*72}\n# DEPTH-MATCHED ARM\n"
      "# Linear recovers fewer reads at variant sites, so it has less power and\n"
      "# would make fewer calls of any kind. Subsample both aligners to the same\n"
      "# per-site depth, preserving each one's own allele proportion.\n"
      f"{'#'*72}")
sel = het_l.values
nmin = np.minimum(dep_g, dep_l)
reps = []
for seed in range(5):
    rng = np.random.default_rng(seed)
    kg = rng.hypergeometric(alt_g, ref_g, nmin)
    kl = rng.hypergeometric(alt_l, ref_l, nmin)
    r = analyse(sel, f"C. Linear-selected sites, depth-matched (replicate seed={seed})",
                kg, nmin, kl, nmin)
    reps.append(r)

print(f"\n{'='*72}\nDEPTH-MATCHED SUMMARY over 5 replicates\n{'='*72}")
for k in ["g", "l", "lonly", "gonly", "lonly_ref", "lonly_alt", "gonly_ref", "gonly_alt"]:
    v = [r[k] for r in reps]
    print(f"  {k:<12} {np.mean(v):>10.1f}  (range {min(v)}-{max(v)})")

df.to_csv(f"{D}/t2_sites.tsv.gz", sep="\t", index=False, compression="gzip")
print("\nwrote t2_sites.tsv.gz")
