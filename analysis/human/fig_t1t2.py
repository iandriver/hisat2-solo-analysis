#!/usr/bin/env python3
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
G, L = "#2166ac", "#b2182b"

df = pd.read_csv(f"{D}/t2_sites.tsv.gz", sep="\t", dtype={"chrom": str})
df = df[df.biallelic & (df.dep_g >= 20) & (df.dep_l >= 20)]
afp = (df.alt_g + df.alt_l) / (df.dep_g + df.dep_l)
cls = pd.cut(afp, [-.01, .02, .15, .85, .98, 1.01],
             labels=["hom-ref", "low", "het", "high", "hom-alt"])

fig, ax = plt.subplots(2, 2, figsize=(12, 8.6))

# (a) dose response
a = ax[0, 0]
ks = ["hom-ref", "het", "hom-alt"]
vals = [np.median(df.dep_g[cls == k] / df.dep_l[cls == k]) for k in ks]
ns = [(cls == k).sum() for k in ks]
b = a.bar(range(3), vals, color=["#999999", G, "#08306b"], width=.55)
a.axhline(1.0, color="k", lw=.9, ls="--")
for i, (v, n) in enumerate(zip(vals, ns)):
    a.text(i, v + .004, f"{v:.4f}\nn={n:,}", ha="center", va="bottom", fontsize=9)
a.set_xticks(range(3)); a.set_xticklabels(["0 alt copies", "1 (het)", "2 (hom-alt)"])
a.set_ylim(.97, 1.16); a.set_ylabel("median read depth, graph / linear")
a.set_title("(a) The depth advantage scales with alt-allele dosage", fontsize=11, loc="left")

# (b) ALT fraction at het sites
a = ax[0, 1]
h = cls == "het"
bins = np.linspace(0, 1, 61)
a.hist(df.af_g[h], bins=bins, histtype="step", lw=2, color=G, label="graph (SNP-aware)")
a.hist(df.af_l[h], bins=bins, histtype="step", lw=2, color=L, label="linear")
a.axvline(.5, color="k", lw=.9, ls="--")
a.set_xlabel("alternate-allele fraction"); a.set_ylabel(f"heterozygous sites (n={h.sum():,})")
a.legend(frameon=False, fontsize=9)
a.set_title(f"(b) median {np.median(df.af_g[h]):.4f} vs {np.median(df.af_l[h]):.4f}",
            fontsize=11, loc="left")

# (c) discordant ASE calls, depth-matched replicate means
a = ax[1, 0]
lab = ["linear-only calls\n(false imbalance)", "graph-only calls\n(imbalance linear missed)"]
ref = [949.6, 23.4]; alt = [36.2, 316.8]
x = np.arange(2)
a.bar(x, ref, .5, label="skewed to REFERENCE allele", color=L)
a.bar(x, alt, .5, bottom=ref, label="skewed to ALTERNATE allele", color=G)
for i in range(2):
    a.text(i, ref[i] + alt[i] + 15, f"{ref[i]+alt[i]:.0f}", ha="center", fontsize=10)
a.set_xticks(x); a.set_xticklabels(lab, fontsize=9)
a.set_ylabel("sites (mean of 5 depth-matched replicates)")
a.legend(frameon=False, fontsize=9)
a.set_title("(c) Discordant ASE calls at 12,580 het sites, FDR<0.05", fontsize=11, loc="left")

# (d) per-gene effect
a = ax[1, 1]
gd = pd.read_csv(f"{D}/t1_gene_delta.tsv", sep="\t")
sel = ["HLA-DQA1", "HLA-DQB1", "HLA-DRB1", "HLA-C", "HLA-DRA", "HLA-A", "HLA-B",
       "B2M", "RPS28", "RPS27", "RPL10", "RPSA"]
sub = gd.set_index("gname").loc[sel]
lr = np.log2((sub.graph + 1) / (sub.linear + 1))
cols = [G if v > 0 else L for v in lr]
a.barh(range(len(sel))[::-1], lr.values, color=cols, height=.62)
a.axvline(0, color="k", lw=.9)
for i, (g, v) in enumerate(zip(sel, lr.values)):
    a.text(v + (.12 if v > 0 else -.12), len(sel) - 1 - i, f"{2**v:.2f}x",
           va="center", ha="left" if v > 0 else "right", fontsize=8.5)
a.set_yticks(range(len(sel))[::-1]); a.set_yticklabels(sel, fontsize=9)
a.set_xlabel("log2 (graph / linear) UMIs"); a.set_xlim(-2.6, 6.6)
a.set_title("(d) Polymorphic HLA genes gain; ribosomal proteins lose", fontsize=11, loc="left")

for r in ax:
    for a in r:
        a.spines[["top", "right"]].set_visible(False)
fig.suptitle("Variant-aware vs linear alignment, same 66.6M reads (10x PBMC 1k v3)",
             fontsize=12.5, y=.985)
fig.tight_layout(rect=[0, 0, 1, .96])
fig.savefig(f"{D}/t1t2.png", dpi=160)
print("wrote t1t2.png")
