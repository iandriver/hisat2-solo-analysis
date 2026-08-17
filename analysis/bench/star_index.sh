#!/bin/bash
# Human STAR index at stock settings -- deliberately not tuned down with
# --genomeSAsparseD, because the memory figure being compared is the one a user
# actually gets when they follow the STARsolo docs.
set -euo pipefail
B="$(cd "$(dirname "$0")" && pwd)"
P=${P:-16}
free_gb(){ df -g /System/Volumes/Data | tail -1 | awk '{print $4}'; }

[ -f "$B/idx/star/SA" ] && { echo "star index already present"; exit 0; }
f=$(free_gb)
[ "$f" -ge 38 ] || { echo "ABORT: STAR index needs ~38GB headroom, only ${f}GB free"; exit 1; }

mkdir -p "$B/idx/star"
echo "[$(date -u +%H:%M:%S)] STAR index build start, ${f}GB free"
/usr/bin/time -l STAR \
  --runMode genomeGenerate \
  --genomeDir "$B/idx/star" \
  --genomeFastaFiles "$B/ref/ensembl_from_index.fa" \
  --sjdbGTFfile "$B/ref/gencode.v50.ensnames.gtf" \
  --sjdbOverhang 100 \
  --runThreadN "$P" \
  --outFileNamePrefix "$B/logs/starbuild_" \
  > "$B/logs/star_index.out" 2> "$B/logs/star_index.time"
echo "[$(date -u +%H:%M:%S)] STAR index done, $(free_gb)GB free"
du -sh "$B/idx/star"
