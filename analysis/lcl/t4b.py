#!/usr/bin/env python3
"""T4b: does variant-aware alignment cut the scRNA-seq genotype error rate?

Yao & Gazal (doi:10.1101/2025.03.25.645175) call genotypes from scRNA-seq reads
aligned with STARsolo -- a linear reference -- and measure an 8% error rate
against SNP-array truth, noting that ADMIXTURE's ancestry error rises
exponentially above ~6%. Reference bias undercounts the alternate allele, which
is the right direction to cause some of that. This measures it directly against
the GIAB HG001 v4.2.1 benchmark for GM12878.

Genotypes are called by allele fraction rather than with GATK, deliberately:
the question is what the aligner contributes, and a threshold caller adds no
model of its own. Sites are restricted to the GIAB high-confidence regions.
"""
import sys
import numpy as np
import pandas as pd

L = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/lcl"
D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "GM12878"
MINDEPTH = int(sys.argv[2]) if len(sys.argv) > 2 else 10
LO, HI = 0.10, 0.90
BASES = list("ACGT")

print(f"=== {SAMPLE}, min depth {MINDEPTH} ===\n", flush=True)

print("loading dbSNP positions...", flush=True)
snp = pd.read_csv(f"{D}/snp_pos.tsv", sep="\t", header=None,
                  names=["rsid", "chrom", "pos0", "alt"], dtype={"chrom": str})
snp["pos"] = snp.pos0 + 1
snp = snp[snp.alt.isin(BASES)].drop_duplicates(subset=["chrom", "pos"], keep=False)

print("loading GIAB truth...", flush=True)
truth = pd.read_csv(f"{L}/hg001_truth.tsv", sep="\t", header=None,
                    names=["chrom", "pos", "ref", "alt", "tgt"], dtype={"chrom": str})

print("loading high-confidence regions...", flush=True)
bed = pd.read_csv(f"{L}/HG001.bed", sep="\t", header=None,
                  names=["chrom", "s", "e"], dtype={"chrom": str})
bed["chrom"] = bed.chrom.str.replace("^chr", "", regex=True)
hc = {c: (g.s.values, g.e.values) for c, g in bed.sort_values(["chrom", "s"]).groupby("chrom")}


def in_hc(chrom, pos):
    out = np.zeros(len(pos), dtype=bool)
    for c in np.unique(chrom):
        if c not in hc:
            continue
        s, e = hc[c]
        m = chrom == c
        i = np.searchsorted(s, pos[m], side="right") - 1
        ok = i >= 0
        v = np.zeros(m.sum(), dtype=bool)
        v[ok] = pos[m][ok] < e[np.clip(i[ok], 0, len(e) - 1)]
        out[m] = v
    return out


def counts(tag):
    return pd.read_csv(f"{L}/cnt_{SAMPLE}_{tag}.tsv", sep="\t", header=None,
                       names=["chrom", "pos", "A", "C", "G", "T"], dtype={"chrom": str})


print("loading pileups...", flush=True)
m = snp.merge(counts("graph"), on=["chrom", "pos"]) \
       .merge(counts("linear"), on=["chrom", "pos"], suffixes=("_g", "_l"))
print(f"  {len(m):,} dbSNP positions covered by both alignments", flush=True)

# truth: het / hom-alt from the VCF, otherwise hom-ref if inside high-confidence
m = m.merge(truth[["chrom", "pos", "ref", "alt", "tgt"]].rename(columns={"alt": "talt"}),
            on=["chrom", "pos"], how="left")
# A position inside the high-confidence BED that the benchmark VCF does not
# touch is genuinely hom-ref. But the truth table above keeps only PASS
# biallelic SNVs, so indel and multiallelic positions would otherwise be
# mislabelled hom-ref -- drop them rather than score against a wrong truth.
allvar = pd.read_csv(f"{L}/hg001_allvar_pos.txt", sep="\t", header=None,
                     names=["chrom", "pos"], dtype={"chrom": str})
allvar["v"] = True
m = m.merge(allvar, on=["chrom", "pos"], how="left")
drop = m["tgt"].isna() & m.v.fillna(False)
print(f"  dropping {int(drop.sum()):,} positions the VCF touches with a "
      f"non-SNV or non-PASS record", flush=True)
m = m[~drop].reset_index(drop=True)
m["tgt"] = m["tgt"].fillna("hom-ref")
keep = in_hc(m.chrom.values, m.pos.values)
m = m[keep].reset_index(drop=True)
print(f"  {len(m):,} inside GIAB high-confidence regions", flush=True)

# where the VCF gives an ALT it wins; elsewhere use the dbSNP ALT
m["ALT"] = np.where(m.talt.notna(), m.talt, m.alt)
m = m[m.ALT.isin(BASES)].reset_index(drop=True)

Ag = m[["A_g", "C_g", "G_g", "T_g"]].values.astype(np.int64)
Al = m[["A_l", "C_l", "G_l", "T_l"]].values.astype(np.int64)
rows = np.arange(len(m))
ai = m.ALT.map({b: i for i, b in enumerate(BASES)}).values
pool = Ag + Al
tmp = pool.copy(); tmp[rows, ai] = -1
ri = tmp.argmax(axis=1)

res = {}
for tag, A in (("graph", Ag), ("linear", Al)):
    alt, ref = A[rows, ai], A[rows, ri]
    dep = alt + ref
    with np.errstate(invalid="ignore", divide="ignore"):
        af = alt / dep
    call = np.full(len(m), "no-call", dtype=object)
    ok = dep >= MINDEPTH
    call[ok & (af < LO)] = "hom-ref"
    call[ok & (af >= LO) & (af <= HI)] = "het"
    call[ok & (af > HI)] = "hom-alt"
    res[tag] = (call, dep, af)

both = (res["graph"][1] >= MINDEPTH) & (res["linear"][1] >= MINDEPTH)
print(f"  {both.sum():,} sites callable at depth >={MINDEPTH} by both\n")

t = m["tgt"].values
print(f"{'':16}{'graph':>22}{'linear':>22}")
print(f"{'truth':<16}{'n':>8}{'wrong':>7}{'rate':>7}{'n':>8}{'wrong':>7}{'rate':>7}")
tot_err = {}
for cls in ["hom-ref", "het", "hom-alt"]:
    sel = both & (t == cls)
    row = f"{cls:<16}"
    for tag in ("graph", "linear"):
        c = res[tag][0][sel]
        wrong = int((c != cls).sum())
        tot_err.setdefault(tag, [0, 0])
        tot_err[tag][0] += wrong
        tot_err[tag][1] += int(sel.sum())
        row += f"{sel.sum():>8,}{wrong:>7,}{100*wrong/max(sel.sum(),1):>6.2f}%"
    print(row)
print(f"{'-'*60}")
row = f"{'ALL':<16}"
for tag in ("graph", "linear"):
    w, n = tot_err[tag]
    row += f"{n:>8,}{w:>7,}{100*w/n:>6.2f}%"
print(row)

# the mechanism: heterozygous sites miscalled as homozygous reference
print("\n=== allele dropout: truth het, called hom-ref ===")
sel = both & (t == "het")
for tag in ("graph", "linear"):
    c = res[tag][0][sel]
    drop = int((c == "hom-ref").sum())
    print(f"  {tag:<8}{drop:>8,} of {sel.sum():,} het sites ({100*drop/max(sel.sum(),1):.2f}%)")
print("\n=== median ALT fraction at truth-het sites ===")
for tag in ("graph", "linear"):
    af = res[tag][2][sel]
    print(f"  {tag:<8}{np.nanmedian(af):.4f}")

# and false heterozygosity at truth hom-ref
print("\n=== false het: truth hom-ref, called het ===")
sel = both & (t == "hom-ref")
for tag in ("graph", "linear"):
    c = res[tag][0][sel]
    fp = int((c == "het").sum())
    print(f"  {tag:<8}{fp:>8,} of {sel.sum():,} ({100*fp/max(sel.sum(),1):.3f}%)")

out = m[["chrom", "pos", "tgt"]].copy()
out["call_g"], out["dep_g"], out["af_g"] = res["graph"]
out["call_l"], out["dep_l"], out["af_l"] = res["linear"]
out[both].to_csv(f"{L}/t4b_{SAMPLE}_d{MINDEPTH}.tsv.gz", sep="\t",
                 index=False, compression="gzip")
