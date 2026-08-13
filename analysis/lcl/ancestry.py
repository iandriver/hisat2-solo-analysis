#!/usr/bin/env python3
"""Does reference bias scale with genetic distance from GRCh38?

Two lymphoblastoid lines, same protocol, same depth (50M cDNA reads each),
same indexes: GM12878 (NA12878, CEU, European) and GM18502 (NA18502, YRI,
African). GRCh38 is predominantly of European ancestry, so if bias tracks
distance from the reference the African-ancestry line should show more of it.

No external genotypes needed: het sites are called from the data, identically
for both samples, and the comparison is between samples rather than against a
truth set.
"""
import numpy as np
import pandas as pd

L = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/lcl"
D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
MINDEPTH = 20
BASES = list("ACGT")

snp = pd.read_csv(f"{D}/snp_pos.tsv", sep="\t", header=None,
                  names=["rsid", "chrom", "pos0", "alt"], dtype={"chrom": str})
snp["pos"] = snp.pos0 + 1
snp = snp[snp.alt.isin(BASES)].drop_duplicates(subset=["chrom", "pos"], keep=False)

rows = []
per_sample = {}
for s, pop in (("GM12878", "CEU / European"), ("GM18502", "YRI / African")):
    def c(tag):
        return pd.read_csv(f"{L}/cnt_{s}_{tag}.tsv", sep="\t", header=None,
                           names=["chrom", "pos", "A", "C", "G", "T"], dtype={"chrom": str})
    m = snp.merge(c("graph"), on=["chrom", "pos"]).merge(
        c("linear"), on=["chrom", "pos"], suffixes=("_g", "_l"))
    Ag = m[["A_g", "C_g", "G_g", "T_g"]].values.astype(np.int64)
    Al = m[["A_l", "C_l", "G_l", "T_l"]].values.astype(np.int64)
    r = np.arange(len(m))
    ai = m.alt.map({b: i for i, b in enumerate(BASES)}).values
    pool = Ag + Al
    tmp = pool.copy(); tmp[r, ai] = -1
    ri = tmp.argmax(axis=1)
    bi = (pool[r, ri] + pool[r, ai]) >= 0.95 * pool.sum(axis=1)

    ag, rg = Ag[r, ai], Ag[r, ri]
    al, rl = Al[r, ai], Al[r, ri]
    dg, dl = ag + rg, al + rl
    with np.errstate(invalid="ignore", divide="ignore"):
        afg, afl = ag / dg, al / dl

    cov = bi & (dg >= MINDEPTH) & (dl >= MINDEPTH)
    afp = (ag + al) / (dg + dl)
    cls = pd.cut(afp, [-.01, .02, .15, .85, .98, 1.01],
                 labels=["hom-ref", "low", "het", "high", "hom-alt"])
    het = cov & np.asarray(cls == "het")
    homalt = cov & np.asarray(cls == "hom-alt")
    nonref = het.sum() + 2 * homalt.sum()

    per_sample[s] = dict(
        pop=pop, covered=int(cov.sum()), het=int(het.sum()),
        homalt=int(homalt.sum()),
        het_frac=100 * het.sum() / cov.sum(),
        nonref_per_1k=1000 * nonref / cov.sum(),
        af_g=float(np.nanmedian(afg[het])), af_l=float(np.nanmedian(afl[het])),
        deficit=float(0.5 - np.nanmedian(afl[het])),
        depth_ratio_het=float(np.nanmedian(dg[het] / dl[het])),
        depth_ratio_homref=float(np.nanmedian(
            dg[cov & np.asarray(cls == "hom-ref")] / dl[cov & np.asarray(cls == "hom-ref")])),
        depth_ratio_homalt=float(np.nanmedian(dg[homalt] / dl[homalt])),
        graph_higher=100 * float(np.nanmean(afg[het] > afl[het])),
    )

print(f"{'':<26}{'GM12878':>16}{'GM18502':>16}")
print(f"{'population':<26}{'CEU / EUR':>16}{'YRI / AFR':>16}")
lab = [("sites callable (both)", "covered", "{:,}"),
       ("heterozygous sites", "het", "{:,}"),
       ("het as % of callable", "het_frac", "{:.2f}%"),
       ("non-ref alleles per 1k sites", "nonref_per_1k", "{:.1f}"),
       ("", None, None),
       ("ALT fraction, graph", "af_g", "{:.4f}"),
       ("ALT fraction, linear", "af_l", "{:.4f}"),
       ("linear's deficit from 0.5", "deficit", "{:.4f}"),
       ("graph higher at", "graph_higher", "{:.1f}%"),
       ("", None, None),
       ("depth ratio, 0 alt copies", "depth_ratio_homref", "{:.4f}"),
       ("depth ratio, 1 (het)", "depth_ratio_het", "{:.4f}"),
       ("depth ratio, 2 (hom-alt)", "depth_ratio_homalt", "{:.4f}")]
for name, k, fmt in lab:
    if k is None:
        print()
        continue
    print(f"{name:<26}" + "".join(
        f"{fmt.format(per_sample[s][k]):>16}" for s in ("GM12878", "GM18502")))

a, b = per_sample["GM12878"], per_sample["GM18502"]
print(f"\nheterozygosity ratio YRI/CEU        {b['het_frac']/a['het_frac']:.3f}x")
print(f"linear's ALT deficit ratio YRI/CEU  {b['deficit']/a['deficit']:.3f}x")
pd.DataFrame(per_sample).to_csv(f"{L}/ancestry_summary.tsv", sep="\t")
