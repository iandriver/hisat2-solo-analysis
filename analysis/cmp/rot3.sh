#!/bin/bash
# Three rotations of G,F,B back-to-back. Rotating rather than blocking means
# each config samples the same stretch of machine conditions, so a burst of
# someone else's load hits all three roughly equally instead of landing on one.
set -u
H=/Users/iandriver/Downloads/hisat2
cd "$(dirname "$0")"
: > rot3.txt
run () {
  local lab=$1 spec=$2
  rm -rf "rot_$lab"
  local s=$(python3 -c 'import time;print(time.time())')
  "$H/hisat2" -x idx/mouse -1 sub_R1.fq -2 sub_R2.fq \
      --solo-barcode-mate 1 --solo-cb-whitelist whitelist.txt \
      --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 12 \
      --gene-annotation genes.ht2gm --gene-strand Reverse \
      --gene-feature "$spec" --solo-out-dir "rot_$lab" \
      -p 16 --no-unal -S /dev/null > /dev/null 2>&1
  local e=$(python3 -c 'import time;print(time.time())')
  local t=$(python3 -c "print(f'{$e-$s:.1f}')")
  local l=$(uptime | sed 's/.*load averages: //' | awk '{print $1}')
  printf '%s\t%s\t%s\n' "$lab" "$t" "$l" >> rot3.txt
  echo "$(date -u +%H:%M:%S)  $lab ${t}s (load $l)"
}
for i in 1 2 3; do
  run "gene$i" Gene
  run "full$i" GeneFull
  run "both$i" Gene,GeneFull
done
echo "ROT3_DONE"
