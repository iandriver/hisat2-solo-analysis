#!/bin/bash
# Does building chromosomes TOGETHER cost more than building them separately?
# Sequence repeated across chromosomes must be disambiguated against the whole
# collection, so a combined graph may need more generations and more temp nodes
# than the sum of its parts. That difference is the term missing from any
# per-chromosome extrapolation to the whole genome.
set -uo pipefail
B="$(cd "$(dirname "$0")" && pwd)"
H=/Users/iandriver/Downloads/hisat2
S=/private/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad
T=$B/trace
tag=$1; shift
chrs="$*"
[ -s "$T/gen_$tag.txt" ] && { echo "$tag already done"; exit 0; }
rm -f "$T/$tag.fa" "$T/$tag.snp" "$T/$tag.haplotype"
for c in $chrs; do
  awk -v want=">$c" '$1==want {f=1; print; next} /^>/{if(f)exit} f' "$B/ref/ensembl_from_index.fa" >> "$T/$tag.fa"
  awk -F'\t' -v c="chr$c" -v n="$c" 'BEGIN{OFS="\t"} $3==c {$3=n; print}' "$S/hap/genome.snp"       >> "$T/$tag.snp"
  awk -F'\t' -v c="chr$c" -v n="$c" 'BEGIN{OFS="\t"} $2==c {$2=n; print}' "$S/hap/genome.haplotype" >> "$T/$tag.haplotype"
done
echo "[$(date -u +%H:%M:%S)] $tag ($chrs): $(grep -v '^>' "$T/$tag.fa" | tr -d '\n' | wc -c | tr -d ' ') bp, $(wc -l < "$T/$tag.snp" | tr -d ' ') variants"
/usr/bin/time -l "$H/hisat2-build" -p 8 --verbose \
    --snp "$T/$tag.snp" --haplotype "$T/$tag.haplotype" \
    "$T/$tag.fa" "$T/idx_$tag" > "$T/build_$tag.out" 2> "$T/build_$tag.err"
echo "rc=$?"
grep -E '^Generation ' "$T/build_$tag.err" > "$T/gen_$tag.txt"
echo "[$(date -u +%H:%M:%S)] $tag: $(wc -l < "$T/gen_$tag.txt" | tr -d ' ') generations, peak RSS $(grep 'maximum resident' "$T/build_$tag.err" | awk '{printf "%.2f GB",$1/1073741824}')"
rm -f "$T/idx_$tag".*.ht2 "$T/$tag.fa"
