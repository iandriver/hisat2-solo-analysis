#!/bin/bash
# T3: what is the least memory in which HISAT2-solo completes a human
# single-cell run?
#
# ulimit -v is not enforced on macOS, so the cap has to be a real one: a Linux
# VM (colima, 20 GiB) with a per-run cgroup limit. --memory-swap equal to
# --memory disables swap, so the limit is hard -- the kernel OOM-kills rather
# than swapping, and docker reports exit 137.
set -u
S=/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad
F=/data/t3/pbmc_1k_v3_fastqs
IDX=${IDX:-/data/human/grch38_snp/genome_snp}
THREADS=${THREADS:-8}


for GB in 18 12 10 8 6 5 4; do
  out=$S/t3/run_${GB}g
  rm -rf "$out"; mkdir -p "$out"
  start=$(date +%s)
  docker run --rm --memory=${GB}g --memory-swap=${GB}g \
    -v "$S":/data hisat2-t3 \
    /usr/bin/time -v hisat2 -x $IDX \
      -1 $F/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,$F/pbmc_1k_v3_S1_L002_R1_001.fastq.gz \
      -2 $F/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,$F/pbmc_1k_v3_S1_L002_R2_001.fastq.gz \
      --solo-barcode-mate 1 \
      --solo-cb-whitelist /data/human/whitelist_v3.txt \
      --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 12 \
      --gene-annotation /data/human/genes.ht2gm --gene-strand Forward \
      --solo-out-dir /data/t3/run_${GB}g/Solo.out \
      -p $THREADS --no-unal -S /dev/null \
      > $S/t3/log_${GB}g.txt 2>&1
  rc=$?
  el=$(( $(date +%s) - start ))
  peak=$(awk '/Maximum resident set size/{printf "%.2f", $NF/1048576}' $S/t3/log_${GB}g.txt)
  cells=$(awk -F, '/Estimated Number of Cells/{print $2}' $S/t3/run_${GB}g/Solo.out/Gene/Summary.csv 2>/dev/null)
  rate=$(sed -n 's/^\([0-9.]*\)% overall alignment rate/\1/p' $S/t3/log_${GB}g.txt)
  printf "%-6s exit=%-4s %5ss  peakRSS=%-7s aligned=%-8s cells=%s\n" \
    "${GB}g" "$rc" "$el" "${peak:-–}" "${rate:-–}" "${cells:-–}" | tee -a $S/t3/titrate.log
done
echo done > $S/t3/titrate.done
