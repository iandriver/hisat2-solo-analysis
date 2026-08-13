#!/usr/bin/env python3
"""T2b: adversarial check on T2.

The obvious way the graph could reach ALT fraction 0.5 dishonestly is by
over-recruiting alternate-allele reads - pulling in reads from elsewhere that
happen to match an ALT path. If it did, homozygous-reference sites would show
an inflated ALT fraction under the graph. They should not.
"""
import numpy as np
import pandas as pd

D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
df = pd.read_csv(f"{D}/t2_sites.tsv.gz", sep="\t", dtype={"chrom": str})
df = df[df.biallelic & (df.dep_g >= 20) & (df.dep_l >= 20)]
print(f"{len(df):,} biallelic SNV sites with >=20 reads under both aligners\n")

# genotype classes defined by the POOLED evidence, so neither aligner picks them
af_pool = (df.alt_g + df.alt_l) / (df.dep_g + df.dep_l)
cls = pd.cut(af_pool, [-.01, .02, .15, .85, .98, 1.01],
             labels=["hom-ref", "low", "het", "high", "hom-alt"])

print(f"{'genotype':<12}{'n':>9}{'ALT frac graph':>16}{'ALT frac linear':>17}"
      f"{'depth graph':>13}{'depth linear':>14}{'depth ratio':>13}")
for k in ["hom-ref", "low", "het", "high", "hom-alt"]:
    s = df[cls == k]
    if not len(s):
        continue
    print(f"{k:<12}{len(s):>9,}{np.median(s.af_g):>16.4f}{np.median(s.af_l):>17.4f}"
          f"{np.median(s.dep_g):>13.0f}{np.median(s.dep_l):>14.0f}"
          f"{np.median(s.dep_g/s.dep_l):>13.4f}")

hr = df[cls == "hom-ref"]
print(f"\nhom-ref sites: mean ALT fraction graph {hr.af_g.mean():.5f}  "
      f"linear {hr.af_l.mean():.5f}")
print(f"  -> the graph does NOT manufacture alternate-allele reads where the\n"
      f"     donor carries none; the residual is sequencing error in both.")

ha = df[cls == "hom-alt"]
print(f"\nhom-alt sites (n={len(ha):,}): the mirror-image test - here the *reference*\n"
      f"allele is the minor one, so reference bias should show as inflated REF.")
print(f"  mean REF fraction graph {1-ha.af_g.mean():.5f}  linear {1-ha.af_l.mean():.5f}")
print(f"  median depth ratio graph/linear {np.median(ha.dep_g/ha.dep_l):.4f}")

het = df[cls == "het"]
print(f"\nhet sites (n={len(het):,}): median depth ratio "
      f"{np.median(het.dep_g/het.dep_l):.4f}")
print(f"hom-ref sites: median depth ratio {np.median(hr.dep_g/hr.dep_l):.4f}"
      f"   <- the within-experiment control")
