#!/usr/bin/env python3
"""Side-by-side scanpy analysis of the three Gene matrices.

The question is not whether the matrices are identical -- they cannot be, the
tools align to different indexes -- but whether an analyst following a standard
workflow would reach the same conclusions from any of them.
"""
import sys, os, warnings
warnings.filterwarnings("ignore")
import numpy as np, scipy.io, scipy.sparse as sp
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 0
OUT = sys.argv[1]
RUNS = [("HISAT2", sys.argv[2]), ("STARsolo", sys.argv[3]), ("rustar", sys.argv[4])]


def load(raw_dir):
    """Load a STARsolo-layout raw/ directory as cells x genes."""
    m = scipy.io.mmread(os.path.join(raw_dir, "matrix.mtx")).tocsc()
    genes, names = [], []
    for line in open(os.path.join(raw_dir, "features.tsv")):
        f = line.rstrip("\n").split("\t")
        genes.append(f[0]); names.append(f[1] if len(f) > 1 else f[0])
    bcs = [l.strip() for l in open(os.path.join(raw_dir, "barcodes.tsv"))]
    a = ad.AnnData(sp.csr_matrix(m.T))          # rows = cells
    a.var_names = genes
    a.var["gene_name"] = names
    a.obs_names = bcs
    a.var_names_make_unique()
    return a


def pipeline(a, seed=0):
    """A deliberately ordinary workflow -- no per-tool tuning."""
    sc.pp.filter_cells(a, min_genes=200)
    sc.pp.filter_genes(a, min_cells=3)
    a.var["mt"] = [n.startswith(("mt-", "MT-")) for n in a.var["gene_name"]]
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    a = a[a.obs.pct_counts_mt < 20].copy()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    a.raw = a
    a = a[:, a.var.highly_variable].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=30, svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    sc.tl.leiden(a, resolution=0.5, random_state=seed, flavor="igraph", n_iterations=2)
    return a


ads = {}
for name, path in RUNS:
    a = load(path)
    # Cell calling differs slightly between tools; restrict to barcodes all
    # three called so clustering is compared on the same cells, not on
    # different cell sets.
    ads[name] = a
    print(f"{name:9} raw: {a.n_obs} barcodes x {a.n_vars} genes", flush=True)

common = set(ads["HISAT2"].obs_names)
for n in ads:
    common &= set(ads[n].obs_names)
common = sorted(common)
print(f"\nbarcodes called by all three: {len(common)}")

res = {}
for name in ads:
    a = ads[name][common].copy()
    res[name] = pipeline(a)
    print(f"{name:9} after QC: {res[name].n_obs} cells, "
          f"{len(res[name].obs.leiden.cat.categories)} clusters", flush=True)

# ---- side-by-side UMAP ----
fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))
for ax, name in zip(axes, [r[0] for r in RUNS]):
    a = res[name]
    sc.pl.umap(a, color="leiden", ax=ax, show=False, legend_loc="on data",
               legend_fontsize=8, frameon=False, title=f"{name}  ({a.n_obs} cells)")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "umap_side_by_side.png"), dpi=140, bbox_inches="tight")
print("\nwrote umap_side_by_side.png")

# ---- cluster agreement ----
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
print("\n-- clustering agreement on the shared cells --")
names = [r[0] for r in RUNS]
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = res[names[i]], res[names[j]]
        shared = sorted(set(a.obs_names) & set(b.obs_names))
        la = a[shared].obs.leiden.astype(str).values
        lb = b[shared].obs.leiden.astype(str).values
        print(f"  {names[i]:9} vs {names[j]:9}  ARI {adjusted_rand_score(la, lb):.4f}"
              f"   NMI {normalized_mutual_info_score(la, lb):.4f}   (n={len(shared)})")

# ---- differential expression ----
print("\n-- top marker genes per cluster (Wilcoxon) --")
tops = {}
for name in names:
    a = res[name]
    sc.tl.rank_genes_groups(a, "leiden", method="wilcoxon")
    d = {}
    for cl in a.obs.leiden.cat.categories:
        ids = [x[cl] for x in a.uns["rank_genes_groups"]["names"][:25]]
        d[cl] = [a.raw.var["gene_name"][i] if i in a.raw.var_names else i for i in ids]
    tops[name] = d

# Match clusters between tools by marker overlap, then report agreement.
ref = names[0]
for other in names[1:]:
    print(f"\n  {ref} cluster -> best-matching {other} cluster (top-25 marker Jaccard)")
    jac_all = []
    for cl, genes in tops[ref].items():
        best, bestj = None, -1.0
        for cl2, g2 in tops[other].items():
            j = len(set(genes) & set(g2)) / len(set(genes) | set(g2))
            if j > bestj: best, bestj = cl2, j
        jac_all.append(bestj)
        print(f"    {cl:>3} -> {best:>3}   Jaccard {bestj:.3f}   shared: "
              + ",".join(sorted(set(genes) & set(tops[other][best]))[:6]))
    print(f"    mean best-match Jaccard: {np.mean(jac_all):.3f}")

# ---- pseudobulk correlation ----
print("\n-- pseudobulk (sum over shared cells) --")
pb = {}
for name in names:
    a = ads[name][common]
    pb[name] = np.asarray(a.X.sum(axis=0)).ravel()
gene_index = {g: i for i, g in enumerate(ads[names[0]].var_names)}
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        x, y = pb[names[i]], pb[names[j]]
        n = min(len(x), len(y))
        keep = (x[:n] + y[:n]) > 0
        r = np.corrcoef(np.log1p(x[:n][keep]), np.log1p(y[:n][keep]))[0, 1]
        print(f"  {names[i]:9} vs {names[j]:9}  log pseudobulk r = {r:.5f}")
np.save(os.path.join(OUT, "pseudobulk.npy"),
        {n: pb[n] for n in names}, allow_pickle=True)
for name in names:
    res[name].write(os.path.join(OUT, f"{name}.h5ad"))
print("\nwrote h5ad per tool")
