#!/usr/bin/env python3
"""T1d: does the HLA gain change a cell-type-level conclusion?

Cluster on the graph matrix, label by markers, then compare HLA class II
expression per cell type between the two alignments. Also check whether the
clustering itself moves.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
sc.settings.verbosity = 0

g = pd.read_csv(f"{D}/t1_genes.tsv", sep="\t")
cells = np.load(f"{D}/t1_cells.npy", allow_pickle=True)
Mg = sparse.load_npz(f"{D}/t1_Mg.npz").T.tocsr()      # cells x genes
Ml = sparse.load_npz(f"{D}/t1_Ml.npz").T.tocsr()


def build(M, tag):
    a = ad.AnnData(M.astype(np.float32))
    a.var_names = pd.Index(g.gname.values).astype(str)
    a.var_names_make_unique()
    a.obs_names = pd.Index(cells).astype(str)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


A, B = build(Mg, "graph"), build(Ml, "linear")


def cluster(a, seed=0):
    x = a.copy()
    sc.pp.highly_variable_genes(x, n_top_genes=2000)
    x = x[:, x.var.highly_variable].copy()
    sc.pp.scale(x, max_value=10)
    sc.tl.pca(x, n_comps=30, svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(x, n_neighbors=15, random_state=seed)
    sc.tl.leiden(x, resolution=1.0, key_added="cl", flavor="igraph",
                 n_iterations=2, directed=False, random_state=seed)
    return x.obs.cl.values.astype(str)


cl_g, cl_l = cluster(A), cluster(B)
print("=== does the clustering move? (graph-derived vs linear-derived) ===")
print(f"  clusters: graph {len(set(cl_g))}, linear {len(set(cl_l))}")
print(f"  ARI {adjusted_rand_score(cl_g, cl_l):.4f}   "
      f"NMI {normalized_mutual_info_score(cl_g, cl_l):.4f}")

MARKERS = {
    "CD14 Mono": ["LYZ", "CD14", "S100A9", "VCAN"],
    "CD16 Mono": ["FCGR3A", "MS4A7", "CDKN1C"],
    "T": ["CD3D", "CD3E", "IL7R", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "DC": ["FCER1A", "CD1C", "CLEC10A"],
    "Platelet": ["PPBP", "PF4"],
}
A.obs["cl"] = pd.Categorical(cl_g)
score = {}
for name, mk in MARKERS.items():
    mk = [m for m in mk if m in A.var_names]
    sc.tl.score_genes(A, mk, score_name=f"s_{name}")
    score[name] = A.obs[f"s_{name}"].values
S = pd.DataFrame(score, index=A.obs_names)
lab = S.groupby(pd.Series(cl_g, index=S.index)).mean().idxmax(axis=1)
ct = pd.Series(cl_g, index=S.index).map(lab)
print("\n=== cell types (graph clustering, marker-scored) ===")
print(ct.value_counts().to_string())

tot_g = np.asarray(Mg.sum(axis=1)).ravel().astype(float)
tot_l = np.asarray(Ml.sum(axis=1)).ravel().astype(float)
gi = {n: i for i, n in enumerate(g.gname.values)}

CLASSII = ["HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1",
           "HLA-DPA1", "HLA-DPB1"]
CLASSI = ["HLA-A", "HLA-B", "HLA-C"]


def cpm(M, tot, names):
    idx = [gi[n] for n in names if n in gi]
    s = np.asarray(M[:, idx].sum(axis=1)).ravel().astype(float)
    return 1e4 * s / tot


for label, names in [("HLA class II (DR/DQ/DP)", CLASSII),
                     ("HLA class I (A/B/C)", CLASSI),
                     ("HLA-DQ only (DQA1+DQB1)", ["HLA-DQA1", "HLA-DQB1"])]:
    xg, xl = cpm(Mg, tot_g, names), cpm(Ml, tot_l, names)
    print(f"\n=== {label}: mean expression per 10k UMIs, by cell type ===")
    print(f"{'cell type':<14}{'n':>6}{'graph':>10}{'linear':>10}{'ratio':>8}")
    for t in ct.value_counts().index:
        m = (ct == t).values
        a, b = xg[m].mean(), xl[m].mean()
        print(f"{t:<14}{m.sum():>6}{a:>10.1f}{b:>10.1f}{(a/b if b else np.nan):>8.3f}")
    # the comparison a paper would actually make
    if "CD14 Mono" in set(ct) and "T" in set(ct):
        mm, tt = (ct == "CD14 Mono").values, (ct == "T").values
        rg = xg[mm].mean() / xg[tt].mean()
        rl = xl[mm].mean() / xl[tt].mean()
        print(f"  monocyte : T-cell ratio    graph {rg:.2f}   linear {rl:.2f}"
              f"   ({100*(rg-rl)/rl:+.1f}%)")

pd.DataFrame({"cell": A.obs_names, "celltype": ct.values,
              "cl_graph": cl_g, "cl_linear": cl_l}).to_csv(
    f"{D}/t1_celltypes.tsv", sep="\t", index=False)
