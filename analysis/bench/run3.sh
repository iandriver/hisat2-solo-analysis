#!/bin/bash
# Three HISAT2 runs against the CellRanger-matched gene set, on a quiet machine
# (STAR's 29 GB index has been removed, so nothing competes for page cache).
set -uo pipefail
B="$(cd "$(dirname "$0")" && pwd)"
H=/Users/iandriver/Downloads/hisat2
say(){ echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$B/logs/run3.log"; }
# don't start timing while the index download is still writing to this disk
while pgrep -f 'grch38_linear.tar.gz' >/dev/null 2>&1; do sleep 20; done
say "download settled, starting"
R1=$(ls "$B"/fq/*_R1_001.fastq.gz | sort | paste -sd, -)
R2=$(ls "$B"/fq/*_R2_001.fastq.gz | sort | paste -sd, -)
for i in 1 2 3; do
  out=$B/cr_h_$i; rm -rf "$out"; mkdir -p "$out"
  say "run $i start"
  /usr/bin/time -l "$H/hisat2" \
    -x "$(cat "$B/idx/ht2_prefix")" -1 "$R1" -2 "$R2" \
    --solo-barcode-mate 1 --solo-cb-whitelist "$B/ref/whitelist_v3.txt" \
    --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 12 \
    --solo-cb-match 1MM --solo-umi-dedup 1MM_CR \
    --solo-cell-filter CellRanger2.2 --solo-expected-cells 3000 \
    --gene-annotation "$B/ref/genes_crset.ht2gm" --gene-strand Forward \
    --gene-feature Gene \
    --solo-out-dir "$out" \
    -p 16 --no-unal -S /dev/null \
    > "$B/logs/cr_h_$i.out" 2> "$B/logs/cr_h_$i.time"
  say "run $i done: $(grep -E '^ *[0-9.]+ real' "$B/logs/cr_h_$i.time" | awk '{print $1"s"}') rss=$(grep 'maximum resident' "$B/logs/cr_h_$i.time" | awk '{printf "%.2fGB", $1/1073741824}')"
  [ $i -lt 3 ] && rm -rf "$B/cr_h_$((i-1))" 2>/dev/null
done
say "RUN3 COMPLETE"
