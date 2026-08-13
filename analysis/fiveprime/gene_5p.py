#!/usr/bin/env python3
"""Gene-level graph-vs-linear UMIs on 5' data, alongside the 3' result."""
import numpy as np
import pandas as pd
from scipy import sparse

H5 = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/fiveprime"
H3 = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"


def load(d):
    feats = pd.read_csv(f"{d}/raw/features.tsv", sep="\t", header=None,
                        names=["gid", "gname", "kind"])
    bcs = pd.read_csv(f"{d}/raw/barcodes.tsv", sep="\t", header=None)[0].values
    with open(f"{d}/raw/matrix.mtx") as fh:
        for line in fh:
            if not line.startswith("%"):
                nr, nc, _ = [int(x) for x in line.split()]
                break
    coo = pd.read_csv(f"{d}/raw/matrix.mtx", sep=" ", skiprows=3, header=None,
                      names=["r", "c", "v"], dtype=np.int64)
    M = sparse.csr_matrix((coo.v.values, (coo.r.values - 1, coo.c.values - 1)),
                          shape=(nr, nc))
    cells = set(pd.read_csv(f"{d}/filtered/barcodes.tsv", header=None)[0])
    return feats, bcs, M, cells


fg, bg, Mg, cg = load(f"{H5}/out_5p_graph/Gene")
fl, bl, Ml, cl = load(f"{H5}/out_5p_linear/Gene")
assert (fg.gid.values == fl.gid.values).all() and (bg == bl).all()
keep = cg & cl
print(f"cells: graph {len(cg)}, linear {len(cl)}, common {len(keep)}, "
      f"jaccard {len(keep)/len(cg | cl):.4f}")
idx = np.where(pd.Series(bg).isin(keep).values)[0]
Mg, Ml = Mg[:, idx], Ml[:, idx]
sg = np.asarray(Mg.sum(axis=1)).ravel()
sl = np.asarray(Ml.sum(axis=1)).ravel()
genes = fg.gname.values
print(f"UMIs in cells: graph {sg.sum():,}  linear {sl.sum():,}  "
      f"({100*(sg.sum()-sl.sum())/sl.sum():+.3f}%)\n")

expressed = (sg + sl) >= 500
lr = np.log2((sg[expressed] + 1) / (sl[expressed] + 1))
print(f"background: {expressed.sum():,} genes >=500 UMIs, "
      f"median log2(graph/linear) {np.median(lr):+.4f}\n")

# 3' numbers for the same genes
t3 = pd.read_csv(f"{H3}/t1_gene_delta.tsv", sep="\t").set_index("gname")

SEL = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
       "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1",
       "HLA-DPA1", "HLA-DPB1", "B2M", "ACTB", "CD3D", "LYZ"]
print(f"{'gene':<11}{'--------- 5-prime ---------':^30}{'--- 3-prime ---':^20}")
print(f"{'':<11}{'graph':>10}{'linear':>9}{'ratio':>8}{'ratio':>12}")
for g in SEL:
    i = np.where(genes == g)[0]
    if len(i) == 0:
        continue
    i = i[0]
    r5 = sg[i] / sl[i] if sl[i] else float("nan")
    if g in t3.index:
        row = t3.loc[g]
        row = row.iloc[0] if isinstance(row, pd.DataFrame) else row
        r3 = row.graph / row.linear if row.linear else float("nan")
    else:
        r3 = float("nan")
    print(f"{g:<11}{sg[i]:>10,}{sl[i]:>9,}{r5:>8.3f}{r3:>12.3f}")

order = np.argsort(-np.where(expressed, np.log2((sg + 1) / (sl + 1)), -np.inf))
print("\n=== top 15 gainers on 5' data ===")
n = 0
for j in order:
    if not expressed[j]:
        continue
    print(f"  {genes[j]:<16}{np.log2((sg[j]+1)/(sl[j]+1)):>+8.3f}"
          f"  graph {sg[j]:>8,}  linear {sl[j]:>8,}")
    n += 1
    if n >= 15:
        break

pd.DataFrame({"gname": genes, "graph": sg, "linear": sl}).to_csv(
    f"{H5}/gene_delta_5p.tsv", sep="\t", index=False)
