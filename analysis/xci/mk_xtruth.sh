#!/bin/bash
# chrX heterozygous sites for the two donors, from the 1000G phased panel.
# Both are female (verified: diploid 0|0 calls on non-PAR chrX), so a het call
# means one allele on the active X and one on the inactive X -- which is what
# makes escape measurable at all.
set -euo pipefail
cd /Volumes/IanSSD/xval
for pair in NA12878:na12878 NA18502:na18502; do
  s=${pair%%:*}; o=${pair##*:}
  bcftools view -s "$s" -v snps -Ou chrX.vcf.gz \
   | bcftools view -i 'GT="het"' -Ov \
   | awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && length($4)==1 && length($5)==1 {print $1,$2,$4,$5}' \
   > ${o}.chrX.het.tsv
  n=$(wc -l < ${o}.chrX.het.tsv)
  # PAR genes escape X-inactivation by definition, so they are a positive
  # control rather than a discovery set; count them separately.
  par=$(awk '$2<=2781479 || $2>=155701383' ${o}.chrX.het.tsv | wc -l)
  echo "  $s: $n chrX het SNVs  (PAR $par, non-PAR $((n-par)))"
done
