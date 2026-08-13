#!/usr/bin/env python3
"""Build a targeted mini-reference for the pseudogene-filter experiment.

A whole-genome SNP-aware rebuild needs ~160 GB, and processed pseudogenes are
almost never on their parent's chromosome (89 of 1,501 pairs), so a
chromosome-scale build cannot reproduce the parent/pseudogene competition
either. Instead extract every locus involved -- ribosomal protein genes, their
pseudogenes, the HLA locus as a positive control, and ordinary expressed genes
as neutral controls -- into one small reference, translating variant and gene
coordinates along with it.
"""
import re, random
import numpy as np, pandas as pd

D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
O = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/filt"
FLANK_BIG, FLANK_SMALL = 5000, 1500

gd = pd.read_csv(f"{D}/t1_gene_delta.tsv", sep="\t")
co = pd.read_csv(f"{D}/gene_coords.tsv", sep="\t", header=None,
                 names=["c", "s", "e", "strand", "gname"], dtype={"c": str})
fai = pd.read_csv(f"{D}/genome.fa.fai", sep="\t", header=None,
                  usecols=[0, 1], names=["c", "len"], dtype={"c": str})
chrlen = dict(zip(fai.c, fai.len))

m = gd.merge(co, on="gname", how="inner")
m = m[m.c.isin(chrlen)]

rp_parent = m.gname.str.match(r"^RP[LS]\d+[A-Z]*$", na=False)
rp_pseudo = m.gname.str.match(r"^RP[LS]\d+[A-Z]*P\d+$", na=False)
hla = m.gname.str.startswith("HLA-", na=False)

# neutral controls: expressed protein-coding genes that moved by <1% either way
expr = (m.graph + m.linear) >= 2000
neutral = m[(m.biotype == "protein_coding") & expr & ~rp_parent & ~hla &
            ((m.graph - m.linear).abs() <= 0.01 * (m.linear + 1))]
random.seed(0)
neutral = neutral.iloc[sorted(random.sample(range(len(neutral)),
                                            min(250, len(neutral))))]

groups = {"rp_parent": m[rp_parent], "rp_pseudo": m[rp_pseudo],
          "hla": m[hla], "control": neutral}
sel = pd.concat([g.assign(grp=k) for k, g in groups.items()])
sel = sel.drop_duplicates(subset=["gname", "c", "s", "e"])
for k, g in groups.items():
    print(f"  {k:<12} {len(g):>5} genes")

# 0-based half-open, flanked, merged per chromosome
rows = []
for _, r in sel.iterrows():
    fl = FLANK_SMALL if r.grp == "rp_pseudo" else FLANK_BIG
    rows.append((r.c, max(0, r.s - 1 - fl), min(chrlen[r.c], r.e + fl)))
byc = {}
for c, s, e in rows:
    byc.setdefault(c, []).append((s, e))

regions = []
for c in sorted(byc):
    ivs = sorted(byc[c])
    cur = list(ivs[0])
    for s, e in ivs[1:]:
        if s <= cur[1]:
            cur[1] = max(cur[1], e)
        else:
            regions.append((c, cur[0], cur[1])); cur = [s, e]
    regions.append((c, cur[0], cur[1]))
total = sum(e - s for _, s, e in regions)
print(f"\n{len(regions):,} merged regions, {total/1e6:.2f} Mb")

with open(f"{O}/regions.bed", "w") as fh:
    for c, s, e in regions:
        fh.write(f"{c}\t{s}\t{e}\tR{c}_{s}\n")

# lookup for coordinate translation
idx = {}
for c, s, e in regions:
    idx.setdefault(c, []).append((s, e))
starts = {c: np.array([s for s, _ in v]) for c, v in idx.items()}
ends = {c: np.array([e for _, e in v]) for c, v in idx.items()}

def locate(c, a, b):
    """Region fully containing [a,b), as (name, start), else None."""
    if c not in starts:
        return None
    i = np.searchsorted(starts[c], a, side="right") - 1
    if i < 0 or b > ends[c][i]:
        return None
    return f"R{c}_{starts[c][i]}", int(starts[c][i])

# ---- translate the .snp file ------------------------------------------------
kept = dropped_edge = 0
with open(f"{D}/genome_snp.snp") as fin, open(f"{O}/mini.snp", "w") as fout:
    for line in fin:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        c = f[2].split()[0]
        pos = int(f[3])
        if f[1] == "deletion":
            span = int(f[4])
        elif f[1] == "insertion":
            span = 1
        else:
            span = len(f[4])
        hit = locate(c, pos, pos + span)
        if hit is None:
            dropped_edge += 1
            continue
        name, off = hit
        fout.write(f"{f[0]}\t{f[1]}\t{name}\t{pos-off}\t{f[4]}\n")
        kept += 1
print(f"variants translated into the mini reference: {kept:,}")

# ---- translate the GTF ------------------------------------------------------
kept_g = 0
with open(f"{D}/genes.gtf") as fin, open(f"{O}/mini.gtf", "w") as fout:
    for line in fin:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        hit = locate(f[0], int(f[3]) - 1, int(f[4]))
        if hit is None:
            continue
        name, off = hit
        f[0] = name
        f[3] = str(int(f[3]) - off)
        f[4] = str(int(f[4]) - off)
        fout.write("\t".join(f) + "\n")
        kept_g += 1
print(f"GTF records translated: {kept_g:,}")

sel[["gname", "grp", "c", "s", "e", "strand", "graph", "linear", "delta"]] \
    .to_csv(f"{O}/selected_genes.tsv", sep="\t", index=False)
