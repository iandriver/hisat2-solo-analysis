#!/bin/bash
# H2: two LCLs of different ancestry, graph vs linear.
#
#   GM12878 = NA12878, CEU (European)  -- close to GRCh38, which is itself
#                                         largely European-derived
#   GM18502 = NA18502, YRI (Yoruba)    -- the ancestry contrast
#
# The full runs are 71-77 GB each, so only a prefix is taken. head closes the
# pipe and curl stops transferring, so this costs the prefix and not the file.
# Reads are 150 bp paired RNA-seq with no barcode read, so this is pseudobulk --
# the same treatment T4b used, which keeps the two comparable.
set -uo pipefail
B="$(cd "$(dirname "$0")" && pwd)"
L=$B/lcl; mkdir -p $L
N=${N:-50000000}          # reads per donor

grab(){ # name srr subdir mate
  local name=$1 srr=$2 sub=$3 mate=$4
  local out=$L/${name}_R${mate}.fq.gz
  if [ -f "$out.done" ]; then echo "skip $name"; return 0; fi
  echo "[$(date -u +%H:%M:%S)] $name mate$mate -> $N reads"
  curl -sS --http1.1 --max-time 10800 \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR855/$sub/$srr/${srr}_${mate}.fastq.gz" 2>/dev/null \
    | gzip -dc 2>/dev/null | head -n $((N*4)) | gzip -1 > "$out"
  local got=$(( $(gzip -dc "$out" | wc -l) / 4 ))
  echo "[$(date -u +%H:%M:%S)] $name got $got reads, $(ls -l "$out" | awk '{printf "%.2f GB", $5/1073741824}')"
  # a short prefix means the transfer died early; leave it unstamped so a rerun retries
  if [ "$got" -ge $((N * 9 / 10)) ]; then touch "$out.done"; else echo "  SHORT, not stamping"; fi
}

grab GM12878 SRR8551677 007 2
grab GM18502 SRR8551676 006 2
echo "GRAB COMPLETE"
