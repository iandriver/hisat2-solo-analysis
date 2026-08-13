#!/bin/bash
# T1: solo counting against the LINEAR GRCh38 index, everything else identical
# to the graph run that produced out_noallelic/.
set -euo pipefail
D=/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human
H=/Users/iandriver/Downloads/hisat2
F=$D/pbmc_1k_v3_fastqs

/usr/bin/time -l $H/hisat2 \
  -x $D/grch38/genome \
  -1 $F/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,$F/pbmc_1k_v3_S1_L002_R1_001.fastq.gz \
  -2 $F/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,$F/pbmc_1k_v3_S1_L002_R2_001.fastq.gz \
  --solo-barcode-mate 1 \
  --solo-cb-whitelist $D/whitelist_v3.txt \
  --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 12 \
  --gene-annotation $D/genes.ht2gm \
  --gene-strand Forward \
  --solo-out-dir $D/out_linear \
  -p 16 --no-unal -S /dev/null 2> $D/run_linear.log
echo DONE >> $D/run_linear.log
