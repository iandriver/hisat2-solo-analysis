#!/usr/bin/env python3
"""Does dropping pseudogene variants restore unique placement at the parent gene?

Three indexes over the same 19.2 Mb mini reference -- linear, all variants,
variants outside processed pseudogenes -- and the same reads. Counts reads
whose unique (MAPQ 60) alignment starts inside each gene.
"""
import subprocess, sys
import numpy as np, pandas as pd

O = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/filt"

# translated gene coordinates, grouped
genes = []
for line in open(f"{O}/mini.gtf"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "gene":
        continue
    a = f[8]
    def at(k):
        i = a.find(k + ' "')
        return a[i + len(k) + 2:a.find('"', i + len(k) + 2)] if i >= 0 else None
    genes.append((f[0], int(f[3]) - 1, int(f[4]), at("gene_name") or at("gene_id")))
gdf = pd.DataFrame(genes, columns=["c", "s", "e", "gname"]).drop_duplicates()

sel = pd.read_csv(f"{O}/selected_genes.tsv", sep="\t")
gdf = gdf.merge(sel[["gname", "grp"]].drop_duplicates("gname"), on="gname", how="left")
gdf["grp"] = gdf.grp.fillna("other")
gdf = gdf.sort_values(["c", "s"]).reset_index(drop=True)

by_c = {c: (g.s.values, g.e.values, g.index.values)
        for c, g in gdf.groupby("c")}


def count(bam):
    n = np.zeros(len(gdf), dtype=np.int64)
    p = subprocess.Popen(["samtools", "view", "-F", "0x100", "-q", "60", bam],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t", 5)
        v = by_c.get(f[2])
        if v is None:
            continue
        s, e, idx = v
        pos = int(f[3]) - 1
        i = np.searchsorted(s, pos, side="right") - 1
        # gene spans here do not nest, so at most one candidate
        if i >= 0 and pos < e[i]:
            n[idx[i]] += 1
    p.wait()
    return n


res = {}
for tag in ("linear", "snpall", "snpfilt"):
    res[tag] = count(f"{O}/aln_{tag}.bam")
    print(f"counted {tag}: {res[tag].sum():,} uniquely-placed reads in genes",
          flush=True)

for t in res:
    gdf[t] = res[t]
gdf.to_csv(f"{O}/compare_genes.tsv", sep="\t", index=False)

print(f"\n{'group':<12}{'genes':>7}{'linear':>12}{'snp-all':>12}{'snp-filt':>12}"
      f"{'all/lin':>9}{'filt/lin':>10}")
for grp in ["rp_parent", "rp_pseudo", "control"]:
    g = gdf[gdf.grp == grp]
    l, a, f = g.linear.sum(), g.snpall.sum(), g.snpfilt.sum()
    print(f"{grp:<12}{len(g):>7}{l:>12,}{a:>12,}{f:>12,}"
          f"{a/l if l else 0:>9.3f}{f/l if l else 0:>10.3f}")

print(f"\nworst-hit ribosomal protein genes (by snp-all / linear):")
rp = gdf[(gdf.grp == "rp_parent") & (gdf.linear >= 2000)].copy()
rp["all_lin"] = rp.snpall / rp.linear
rp["filt_lin"] = rp.snpfilt / rp.linear
rp = rp.sort_values("all_lin")
print(f"  {'gene':<12}{'linear':>10}{'snp-all':>10}{'snp-filt':>10}"
      f"{'all/lin':>9}{'filt/lin':>10}{'recovered':>11}")
for _, r in rp.head(15).iterrows():
    rec = (r.snpfilt - r.snpall) / (r.linear - r.snpall) if r.linear > r.snpall else float("nan")
    print(f"  {r.gname:<12}{r.linear:>10,}{r.snpall:>10,}{r.snpfilt:>10,}"
          f"{r.all_lin:>9.3f}{r.filt_lin:>10.3f}{100*rec:>10.1f}%")

lost = rp.linear.sum() - rp.snpall.sum()
back = rp.snpfilt.sum() - rp.snpall.sum()
print(f"\nacross those {len(rp)} genes: linear-vs-all-variants deficit {lost:,} reads; "
      f"the filter restores {back:,} ({100*back/lost if lost else 0:.1f}%)")
