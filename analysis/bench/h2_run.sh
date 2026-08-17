#!/bin/bash
# H2: align each donor against the graph and the linear index, count reads per
# gene from the GX:Z tag, and never store a SAM file (they would be ~30 GB each).
#
# Rotation is donor-major and aligner-alternating so that if the machine drifts
# it does not land entirely on one aligner. Timing is NOT the point here and no
# timing claim is made -- see the H1 report for why this machine cannot be timed.
set -uo pipefail
B="$(cd "$(dirname "$0")" && pwd)"
H=/Users/iandriver/Downloads/hisat2
L=$B/lcl
GM=$B/ref/genes_crset.ht2gm
say(){ echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$B/logs/h2.log"; }

run(){ # donor aligner
  local s=$1 a=$2
  local idx
  [ "$a" = graph ] && idx=$B/idx/grch38_snp/genome_snp || idx=$B/idx/grch38/genome
  [ -s "$L/gx_${s}_${a}.tsv" ] && { say "$s $a already done"; return 0; }
  say "$s $a start"
  # Count only uniquely-mapped reads (NH:i:1) so that the comparison is about
  # reads placed confidently, not about multimapper reporting differences.
  "$H/hisat2" -x "$idx" -U "$L/${s}_R2.fq.gz" \
      --gene-annotation "$GM" --gene-strand Unstranded \
      -p 14 --no-unal -S /dev/stdout 2> "$L/aln_${s}_${a}.log" \
    | awk -F'\t' '
        /^@/ { next }
        {
          nh = 0; gx = ""
          for (i = 12; i <= NF; i++) {
            if ($i ~ /^NH:i:/) nh = substr($i, 6) + 0
            else if ($i ~ /^GN:Z:/) gx = substr($i, 6)
          }
          if (nh == 1 && gx != "" && gx != "-") cnt[gx]++
        }
        END { for (g in cnt) print g "\t" cnt[g] }' \
    | sort -k2,2nr > "$L/gx_${s}_${a}.tsv"
  say "$s $a done: $(wc -l < "$L/gx_${s}_${a}.tsv") genes, $(grep 'overall alignment rate' "$L/aln_${s}_${a}.log" | awk '{print $1}')"
}

for s in GM12878 GM18502; do
  for a in graph linear; do run $s $a; done
done
say "H2 RUNS COMPLETE"
