#!/usr/bin/env python3
"""T1: is HLA expression systematically underestimated by linear alignment?

Compares per-cell UMIs from the graph (SNP-aware) and linear GRCh38 indexes,
everything else held constant.
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import wilcoxon

D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
GRAPH = f"{D}/out_noallelic/Solo.out/Gene"
LINEAR = f"{D}/out_linear/Solo.out/Gene"


def load_mtx(d):
    feats = pd.read_csv(f"{d}/raw/features.tsv", sep="\t", header=None,
                        names=["gid", "gname", "kind"])
    bcs = pd.read_csv(f"{d}/raw/barcodes.tsv", sep="\t", header=None)[0].values
    hdr = []
    with open(f"{d}/raw/matrix.mtx") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            hdr = [int(x) for x in line.split()]
            break
    nr, nc, nnz = hdr
    coo = pd.read_csv(f"{d}/raw/matrix.mtx", sep=" ", skiprows=3, header=None,
                      names=["r", "c", "v"], dtype=np.int64)
    M = sparse.csr_matrix((coo.v.values, (coo.r.values - 1, coo.c.values - 1)),
                          shape=(nr, nc))
    return feats, bcs, M


print("loading graph...", flush=True)
fg, bg, Mg = load_mtx(GRAPH)
print("loading linear...", flush=True)
fl, bl, Ml = load_mtx(LINEAR)

assert (fg.gid.values == fl.gid.values).all(), "feature order differs"
assert (bg == bl).all(), "barcode order differs"
genes = fg.gname.values
gids = fg.gid.values

# --- cell set: graph's filtered cells, present in both --------------------
cells_g = set(pd.read_csv(f"{GRAPH}/filtered/barcodes.tsv", header=None)[0])
cells_l = set(pd.read_csv(f"{LINEAR}/filtered/barcodes.tsv", header=None)[0])
print(f"filtered cells: graph {len(cells_g)}, linear {len(cells_l)}, "
      f"intersection {len(cells_g & cells_l)}, jaccard "
      f"{len(cells_g & cells_l)/len(cells_g | cells_l):.4f}")

keep = cells_g & cells_l
idx = np.where(pd.Series(bg).isin(keep).values)[0]
Mg = Mg[:, idx].tocsc()
Ml = Ml[:, idx].tocsc()
ncell = len(idx)
print(f"using {ncell} cells common to both\n")

tot_g = np.asarray(Mg.sum(axis=0)).ravel()
tot_l = np.asarray(Ml.sum(axis=0)).ravel()
print(f"total UMIs   graph {tot_g.sum():,}   linear {tot_l.sum():,}   "
      f"ratio {tot_g.sum()/tot_l.sum():.4f}")
print(f"median UMIs/cell  graph {np.median(tot_g):,.0f}  linear {np.median(tot_l):,.0f}\n")

# --- per-gene totals ------------------------------------------------------
sg = np.asarray(Mg.sum(axis=1)).ravel()
sl = np.asarray(Ml.sum(axis=1)).ravel()

HLA1 = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"]
HLA2 = ["HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1",
        "HLA-DPA1", "HLA-DPB1"]
B2M = ["B2M"]

print("=== per-gene totals across %d cells ===" % ncell)
print(f"{'gene':<12}{'graph':>10}{'linear':>10}{'delta':>9}{'ratio':>8}")
for g in HLA1 + HLA2 + B2M + ["ACTB", "LYZ", "CD3D", "MALAT1"]:
    i = np.where(genes == g)[0]
    if len(i) == 0:
        continue
    i = i[0]
    r = sg[i] / sl[i] if sl[i] else float("nan")
    print(f"{g:<12}{sg[i]:>10,}{sl[i]:>10,}{sg[i]-sl[i]:>9,}{r:>8.3f}")

# --- is HLA an outlier, or does everything gain? -------------------------
expressed = (sg + sl) >= 500
lr = np.full(len(sg), np.nan)
lr[expressed] = np.log2((sg[expressed] + 1) / (sl[expressed] + 1))
n_expr = expressed.sum()
print(f"\n=== background: {n_expr:,} genes with >=500 UMIs pooled ===")
print(f"median log2(graph/linear) = {np.nanmedian(lr):+.4f}  "
      f"(ratio {2**np.nanmedian(lr):.4f})")
print(f"IQR = {np.nanpercentile(lr[expressed],25):+.4f} .. "
      f"{np.nanpercentile(lr[expressed],75):+.4f}")

order = np.argsort(-np.where(np.isnan(lr), -np.inf, lr))
rank_of = {genes[j]: k for k, j in enumerate(order)}
print(f"\n{'gene':<12}{'log2FC':>9}{'pctile':>9}   rank of {n_expr:,}")
for g in HLA1 + HLA2 + B2M:
    i = np.where(genes == g)[0]
    if len(i) == 0 or not expressed[i[0]]:
        continue
    i = i[0]
    pct = 100.0 * (1 - rank_of[g] / n_expr)
    print(f"{g:<12}{lr[i]:>+9.4f}{pct:>8.2f}%   {rank_of[g]+1}")

print("\n=== top 25 genes by log2(graph/linear) ===")
shown = 0
for j in order:
    if not expressed[j]:
        continue
    print(f"  {genes[j]:<16}{lr[j]:>+8.4f}  graph {sg[j]:>8,}  linear {sl[j]:>8,}")
    shown += 1
    if shown >= 25:
        break

# --- paired per-cell test -------------------------------------------------
print("\n=== paired per-cell test (Wilcoxon signed-rank, %d cells) ===" % ncell)
print(f"{'gene':<12}{'med graph':>11}{'med lin':>10}{'mean d':>9}{'up':>7}{'dn':>6}{'p':>12}")
for g in HLA1 + HLA2 + B2M:
    i = np.where(genes == g)[0]
    if len(i) == 0:
        continue
    i = i[0]
    a = np.asarray(Mg[i, :].todense()).ravel().astype(float)
    b = np.asarray(Ml[i, :].todense()).ravel().astype(float)
    d = a - b
    if (d != 0).sum() < 10:
        continue
    p = wilcoxon(a, b, zero_method="wilcox").pvalue
    print(f"{g:<12}{np.median(a):>11.0f}{np.median(b):>10.0f}{d.mean():>9.2f}"
          f"{(d>0).sum():>7}{(d<0).sum():>6}{p:>12.2e}")

np.save(f"{D}/t1_sg.npy", sg)
np.save(f"{D}/t1_sl.npy", sl)
sparse.save_npz(f"{D}/t1_Mg.npz", Mg.tocsr())
sparse.save_npz(f"{D}/t1_Ml.npz", Ml.tocsr())
pd.DataFrame({"gid": gids, "gname": genes}).to_csv(f"{D}/t1_genes.tsv", sep="\t", index=False)
np.save(f"{D}/t1_cells.npy", bg[idx])
print("\nsaved matrices for the cell-type breakdown")
