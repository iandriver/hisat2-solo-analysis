#!/usr/bin/env python3
"""T1b: where does the graph's extra signal actually come from?

Two competing explanations for a gene gaining UMIs under the graph index:
  (i)  recovery  - reads that the linear index failed to align or assign
  (ii) reallocation - reads the linear index gave to a paralog, or called
       multi-gene, that the graph resolves elsewhere
Only (i) is an advantage. This separates them.
"""
import numpy as np
import pandas as pd
from scipy import sparse

D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"

g = pd.read_csv(f"{D}/t1_genes.tsv", sep="\t")
sg = np.load(f"{D}/t1_sg.npy")
sl = np.load(f"{D}/t1_sl.npy")
bt = pd.read_csv(f"{D}/gid_biotype.tsv", sep="\t", header=None,
                 names=["gid", "biotype"])
g = g.merge(bt, on="gid", how="left")
g["graph"], g["linear"] = sg, sl
g["delta"] = sg - sl

tot_g, tot_l = sg.sum(), sl.sum()
print(f"total UMIs in cells:  graph {tot_g:,}  linear {tot_l:,}  "
      f"delta {tot_g-tot_l:+,} ({100*(tot_g-tot_l)/tot_l:+.3f}%)\n")

# ---- 1. how concentrated is the delta? ----------------------------------
up = g[g.delta > 0].sort_values("delta", ascending=False)
dn = g[g.delta < 0].sort_values("delta")
print(f"genes gaining: {len(up):,}  (+{up.delta.sum():,} UMIs)")
print(f"genes losing : {len(dn):,}  ({dn.delta.sum():,} UMIs)")
print(f"genes equal  : {(g.delta==0).sum():,}\n")
print(f"top 10 gainers account for {100*up.delta.head(10).sum()/up.delta.sum():.1f}% "
      f"of all gains")
print(f"top 50 gainers account for {100*up.delta.head(50).sum()/up.delta.sum():.1f}%\n")

# ---- 2. gains by biotype -------------------------------------------------
print("=== net delta by gene biotype (top by |delta|) ===")
bb = g.groupby("biotype").agg(n=("delta", "size"), ngain=("delta", lambda x: (x > 0).sum()),
                              delta=("delta", "sum"), graph=("graph", "sum"),
                              linear=("linear", "sum"))
bb = bb.reindex(bb.delta.abs().sort_values(ascending=False).index)
print(f"{'biotype':<34}{'genes':>7}{'gaining':>9}{'delta':>10}{'graph':>12}{'linear':>12}")
for k, r in bb.head(10).iterrows():
    print(f"{str(k):<34}{r.n:>7,}{r.ngain:>9,}{r.delta:>+10,}{r.graph:>12,}{r.linear:>12,}")
ps = g.biotype.fillna("").str.contains("pseudogene")
print(f"\nall pseudogene classes:  delta {g.delta[ps].sum():+,}  "
      f"({100*g.delta[ps].sum()/(tot_g-tot_l):.1f}% of the net gain)")
print(f"protein_coding:          delta {g.delta[g.biotype=='protein_coding'].sum():+,}")

# ---- 3. paralog families: recovery or reallocation? ---------------------
FAMS = {
    "HLA class I (A,B,C,E,F,G)": ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"],
    "HLA-DQA (DQA1,DQA2)": ["HLA-DQA1", "HLA-DQA2"],
    "HLA-DQB (DQB1,DQB2,DQB3)": ["HLA-DQB1", "HLA-DQB2", "HLA-DQB3"],
    "HLA-DRB (DRB1,5,6,9)": ["HLA-DRB1", "HLA-DRB5", "HLA-DRB6", "HLA-DRB9"],
    "HLA-DR alpha (DRA)": ["HLA-DRA"],
    "HLA-DP (DPA1-3,DPB1-2)": ["HLA-DPA1", "HLA-DPA2", "HLA-DPA3", "HLA-DPB1", "HLA-DPB2"],
}
print("\n=== HLA families: parent + paralogs summed ===")
print(f"{'family':<30}{'graph':>10}{'linear':>10}{'delta':>9}{'ratio':>8}")
for name, members in FAMS.items():
    sub = g[g.gname.isin(members)]
    a, b = sub.graph.sum(), sub.linear.sum()
    print(f"{name:<30}{a:>10,}{b:>10,}{a-b:>+9,}{(a/b if b else np.nan):>8.3f}")

# ---- 4. ribosomal-protein pseudogene families ---------------------------
rp = g[g.gname.str.match(r"^RP[LS]\d+[A-Z]*$", na=False)]
rpp = g[g.gname.str.match(r"^RP[LS]\d+[A-Z]*P\d+$", na=False)]
print(f"\n=== ribosomal protein genes vs their pseudogenes ===")
print(f"{'RP parent genes':<26}n={len(rp):<5} graph {rp.graph.sum():>10,}  "
      f"linear {rp.linear.sum():>10,}  delta {rp.delta.sum():>+9,}")
print(f"{'RP pseudogenes':<26}n={len(rpp):<5} graph {rpp.graph.sum():>10,}  "
      f"linear {rpp.linear.sum():>10,}  delta {rpp.delta.sum():>+9,}")
print(f"{'both':<26}      graph {rp.graph.sum()+rpp.graph.sum():>10,}  "
      f"linear {rp.linear.sum()+rpp.linear.sum():>10,}  "
      f"delta {rp.delta.sum()+rpp.delta.sum():>+9,}")

print("\n=== top 20 losers ===")
for _, r in dn.head(20).iterrows():
    print(f"  {r.gname:<18}{r.delta:>+8,}  graph {r.graph:>8,}  linear {r.linear:>8,}  "
          f"{str(r.biotype)}")

g.to_csv(f"{D}/t1_gene_delta.tsv", sep="\t", index=False)
