#!/usr/bin/env python3
"""Multi-panel comparison figure from the per-tool h5ad files."""
import sys, os, warnings
warnings.filterwarnings("ignore")
import numpy as np, scanpy as sc, anndata as ad
import scipy.io, scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 0
OUT = sys.argv[1]
NAMES = ["HISAT2", "STARsolo", "rustar"]
res = {n: sc.read_h5ad(os.path.join(OUT, f"{n}.h5ad")) for n in NAMES}

# Canonical PBMC markers; these should light up the same cells whichever
# tool produced the matrix.
MARKERS = {"B": "Cd79a", "T": "Cd3e", "Mono/Mac": "Csf1r", "Neutrophil": "Csf3r"}


def gene_vec(a, sym):
    ids = [i for i, s in zip(a.raw.var_names, a.raw.var["gene_name"]) if s == sym]
    if not ids:
        return None
    x = a.raw[:, ids[0]].X
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()


fig = plt.figure(figsize=(17, 13.5))
gs = fig.add_gridspec(3, 3, hspace=0.28, wspace=0.12)

for j, name in enumerate(NAMES):
    a = res[name]
    ax = fig.add_subplot(gs[0, j])
    sc.pl.umap(a, color="leiden", ax=ax, show=False, legend_loc="on data",
               legend_fontsize=9, frameon=False,
               title=f"{name}\n{a.n_obs} cells, {len(a.obs.leiden.cat.categories)} clusters")

for j, name in enumerate(NAMES):
    a = res[name]
    ax = fig.add_subplot(gs[1, j])
    v = gene_vec(a, "Cd79a")
    u = a.obsm["X_umap"]
    s = ax.scatter(u[:, 0], u[:, 1], c=v, s=3, cmap="viridis")
    ax.set_title(f"{name} — Cd79a (B cell)", fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values(): sp_.set_visible(False)
    plt.colorbar(s, ax=ax, fraction=0.03)

# pseudobulk scatters, keyed by gene ID
def pb(name):
    a = res[name]
    x = a.raw.X
    tot = np.asarray(x.sum(axis=0)).ravel()
    return dict(zip(a.raw.var_names, tot))

P = {n: pb(n) for n in NAMES}
pairs = [("HISAT2", "STARsolo"), ("HISAT2", "rustar"), ("STARsolo", "rustar")]
for j, (A, B) in enumerate(pairs):
    ax = fig.add_subplot(gs[2, j])
    shared = sorted(set(P[A]) & set(P[B]))
    x = np.array([P[A][g] for g in shared]); y = np.array([P[B][g] for g in shared])
    keep = (x + y) > 0
    r = np.corrcoef(x[keep], y[keep])[0, 1]
    ax.scatter(x[keep], y[keep], s=2, alpha=0.25, color="#2b6cb0")
    lim = [min(x[keep].min(), y[keep].min()) or 1e-3, max(x[keep].max(), y[keep].max())]
    ax.plot(lim, lim, "--", color="grey", lw=1)
    ax.set_xlabel(f"{A} (normalised, log)"); ax.set_ylabel(f"{B} (normalised, log)")
    ax.set_title(f"pseudobulk  r = {r:.4f}  (n={keep.sum()})", fontsize=11)

fig.suptitle("HISAT2-solo vs STARsolo vs rustar — 10M reads, 5k mouse PBMC, GRCm39",
             fontsize=14, y=0.995)
plt.savefig(os.path.join(OUT, "comparison.png"), dpi=130, bbox_inches="tight")
print("wrote comparison.png")

# ---- marker dotplot, one per tool, same genes and same cluster order ----
fig2, axes = plt.subplots(1, 3, figsize=(19, 4.6))
syms = ["Cd79a", "Cd19", "Cd3e", "Cd3d", "Csf1r", "Cst3", "Csf3r", "Cxcr2", "Ccl5", "Nkg7"]
for ax, name in zip(axes, NAMES):
    a = res[name]
    present = [s for s in syms if gene_vec(a, s) is not None]
    mat = np.array([gene_vec(a, s) for s in present])
    cl = a.obs.leiden.astype(str).values
    order = sorted(set(cl), key=int)
    M = np.array([[mat[i][cl == c].mean() for c in order] for i in range(len(present))])
    im = ax.imshow(M, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, fontsize=9)
    ax.set_yticks(range(len(present))); ax.set_yticklabels(present, fontsize=9)
    ax.set_xlabel("leiden cluster"); ax.set_title(name)
    plt.colorbar(im, ax=ax, fraction=0.035)
fig2.suptitle("Mean marker expression per cluster — same genes, each tool's own clustering", y=1.04)
plt.savefig(os.path.join(OUT, "markers.png"), dpi=130, bbox_inches="tight")
print("wrote markers.png")
