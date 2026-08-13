#!/bin/bash
# A/B on the local-graph explosion retry path. Same fasta, same .snp/.haplotype,
# same machine -- the only variable is whether the retry preserves haplotypes.
set -u
FA=/tmp/chr1.up.fa
cd "$(dirname "$0")"
: > ab_results.txt
for set in kg_phased dbsnp_big4; do
  for bin in BASE PATCHED; do
    rm -rf idx; mkdir -p idx
    S=$(date +%s)
    /usr/bin/time -l /tmp/hisat2-build-s-$bin -p 18 \
        --snp $set.snp --haplotype $set.haplotype "$FA" idx/chr1 \
        > log_${set}_${bin}.txt 2> time_${set}_${bin}.txt
    rc=$?
    E=$(( $(date +%s) - S ))
    ret=$(grep -o 'Local indexes: .*' log_${set}_${bin}.txt | tail -1)
    rss=$(awk '/maximum resident set size/{printf "%.2f", $1/1073741824}' time_${set}_${bin}.txt)
    sz=$(du -sk idx | cut -f1)
    echo -e "$set\t$bin\trc=$rc\tsec=$E\trss_gb=$rss\tindex_kb=$sz\t$ret" >> ab_results.txt
    cat ab_results.txt
  done
done
rm -rf idx
echo AB_DONE >> ab_results.txt
