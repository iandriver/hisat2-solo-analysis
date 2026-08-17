#!/bin/bash
# E0 (revised): measure the prefix-doubling convergence curve without renting a
# 1 TB machine.
#
# PathGraph::printInfo already logs "Generation N (temp -> nodes, ranks)" under
# --verbose (gbwt_graph.h:2301). Node counts are a property of the graph, not of
# index_t width, so a 32-bit build on a subset gives the same curve shape a
# 64-bit whole-genome build would -- for a fraction of the memory.
#
# The questions this answers:
#   1. Does doubling converge, and in how many generations?
#   2. What is peak_nodes / final_nodes?  That ratio times the whole-genome
#      final size is the real memory requirement.
#   3. Is the ratio stable across chromosomes of different size and variant
#      density?  If it is, extrapolation to the whole genome is defensible.
#
# chr6 is deliberately included: it carries the MHC, 2.4x the genome-wide
# variant density, and is where the local budget blows out.
set -uo pipefail
B="$(cd "$(dirname "$0")" && pwd)"
H=/Users/iandriver/Downloads/hisat2
S=/private/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad
T=$B/trace; mkdir -p $T
P=${P:-8}
say(){ echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$T/trace.log"; }

prep(){ # chrom (ensembl name, e.g. 22)
  local c=$1
  [ -s "$T/$c.fa" ] || awk -v want=">$c" '
      $1==want {f=1; print; next} /^>/{if(f)exit} f' "$B/ref/ensembl_from_index.fa" > "$T/$c.fa"
  # variant files use chr-prefixed names; the fasta does not
  [ -s "$T/$c.snp" ] || awk -F'\t' -v c="chr$c" -v n="$c" \
      'BEGIN{OFS="\t"} $3==c {$3=n; print}' "$S/hap/genome.snp" > "$T/$c.snp"
  [ -s "$T/$c.haplotype" ] || awk -F'\t' -v c="chr$c" -v n="$c" \
      'BEGIN{OFS="\t"} $2==c {$2=n; print}' "$S/hap/genome.haplotype" > "$T/$c.haplotype"
}

run(){ # chrom
  local c=$1
  [ -s "$T/gen_$c.txt" ] && { say "chr$c already traced"; return 0; }
  prep $c
  local bp=$(grep -v '^>' "$T/$c.fa" | tr -d '\n' | wc -c | tr -d ' ')
  local nv=$(wc -l < "$T/$c.snp" | tr -d ' ')
  local nh=$(wc -l < "$T/$c.haplotype" | tr -d ' ')
  say "chr$c: $bp bp, $nv variants, $nh haplotypes -- building"
  /usr/bin/time -l "$H/hisat2-build" -p $P --verbose \
      --snp "$T/$c.snp" --haplotype "$T/$c.haplotype" \
      "$T/$c.fa" "$T/idx_$c" > "$T/build_$c.out" 2> "$T/build_$c.err"
  local rc=$?
  grep -E '^Generation ' "$T/build_$c.err" > "$T/gen_$c.txt" 2>/dev/null
  say "chr$c: rc=$rc, $(wc -l < "$T/gen_$c.txt" | tr -d ' ') generations, peak RSS $(grep 'maximum resident' "$T/build_$c.err" | awk '{printf "%.2f GB", $1/1073741824}')"
  # keep the logs, drop the index and the fasta copy -- disk is tight
  rm -f "$T/idx_$c".*.ht2 "$T/$c.fa"
  echo "$c	$bp	$nv	$nh" >> "$T/sizes.tsv"
}

for c in "$@"; do run $c; done
say "TRACE DONE for: $*"
