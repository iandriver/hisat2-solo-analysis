#!/bin/bash
# Rotated timing: H,S,H,S,H,S rather than HHH,SSS.
#
# This is not pedantry. The mouse benchmark ran blocked and a step change in
# another process's load landed entirely on one config, producing a 7x spread
# for identical work. Rotating spreads any drift across both tools, and three
# rounds show whether the machine was quiet.
#
# Intermediate outputs are deleted as they are superseded; only round 3 is kept
# for the concordance comparison, because disk is the binding constraint here.
set -uo pipefail
B="$(cd "$(dirname "$0")" && pwd)"
H=/Users/iandriver/Downloads/hisat2
P=16
WL=$B/ref/whitelist_v3.txt
GM=$B/ref/genes.ht2gm
R1=$(ls "$B"/fq/*_R1_001.fastq.gz | sort | paste -sd, -)
R2=$(ls "$B"/fq/*_R2_001.fastq.gz | sort | paste -sd, -)
say(){ echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$B/logs/rot.log"; }
free_gb(){ df -g /System/Volumes/Data | tail -1 | awk '{print $4}'; }

hisat2_run(){ # round
  local i=$1 out=$B/rot_h_$1
  rm -rf "$out"; mkdir -p "$out"
  say "H$i start ($(free_gb)GB free)"
  /usr/bin/time -l "$H/hisat2" \
     -x "$(cat "$B/idx/ht2_prefix")" -1 "$R1" -2 "$R2" \
     --solo-barcode-mate 1 --solo-cb-whitelist "$WL" \
     --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 12 \
     --solo-cb-match 1MM --solo-umi-dedup 1MM_CR \
     --solo-cell-filter CellRanger2.2 --solo-expected-cells 3000 \
     --gene-annotation "$GM" --gene-strand Forward \
     --gene-feature Gene \
     --solo-out-dir "$out" \
     -p $P --no-unal -S /dev/null \
     > "$B/logs/h_$i.out" 2> "$B/logs/h_$i.time"
  say "H$i done rc=$?"
}

star_run(){ # round
  local i=$1 out=$B/rot_s_$1
  rm -rf "$out"; mkdir -p "$out"
  say "S$i start ($(free_gb)GB free)"
  # STAR 2.7.11b resolves the transcript-info directory nondeterministically on
  # this genome: identical commands fail ~half the time with
  #   "could not open input file /geneInfo.tab"
  # (an empty directory prefix, while geneInfo.tab is present in genomeDir).
  # The failure happens during init, ~30s in and always before mapping starts,
  # and each attempt is its own process -- so retrying costs nothing and does
  # not contaminate the timing of the attempt that succeeds.
  local att=0
  while [ $att -lt 8 ]; do
    att=$((att+1))
    rm -rf "$out"; mkdir -p "$out"
  /usr/bin/time -l STAR \
     --genomeDir "$B/idx/star" \
     --readFilesIn "$R2" "$R1" --readFilesCommand gzcat \
     --soloType CB_UMI_Simple --soloCBwhitelist "$WL" \
     --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
     --soloCBmatchWLtype 1MM --soloUMIdedup 1MM_CR \
     --soloFeatures Gene --soloStrand Forward \
     --soloCellFilter CellRanger2.2 3000 0.99 10 \
     --runThreadN $P --outSAMtype None \
     --outFileNamePrefix "$out/" \
       > "$B/logs/s_$i.out" 2> "$B/logs/s_$i.time"
    if grep -q 'Transcriptome: exiting' "$B/logs/s_$i.time" 2>/dev/null; then
      say "S$i attempt $att hit the STAR geneInfo.tab bug, retrying"
      continue
    fi
    break
  done
  say "S$i done after $att attempt(s)"
}

for i in 1 2 3; do
  hisat2_run $i
  star_run $i
  # keep only the newest pair
  [ $i -gt 1 ] && rm -rf "$B/rot_h_$((i-1))" "$B/rot_s_$((i-1))"
done
say "ROTATION COMPLETE"
