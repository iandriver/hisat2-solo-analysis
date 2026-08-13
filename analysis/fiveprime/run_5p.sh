#!/bin/bash
# 5' v2 PBMC 1k: graph (SNP-aware) vs linear GRCh38, everything else identical.
# 5' v2 differs from 3' v3 in two ways that matter: a 10 bp UMI (not 12), the
# 737K barcode list (not 3M), and reverse-strand cDNA. All three confirmed
# empirically before this run.
set -u
F=/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/fiveprime
D=/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human
H=/Users/iandriver/Downloads/hisat2
P=$F/sc5p_v2_hs_PBMC_1k_5gex_S1

run () {  # run <tag> <index>
  /usr/bin/time -l $H/hisat2 -x "$2" \
    -1 ${P}_L001_R1_001.fastq.gz,${P}_L002_R1_001.fastq.gz \
    -2 ${P}_L001_R2_001.fastq.gz,${P}_L002_R2_001.fastq.gz \
    --solo-barcode-mate 1 \
    --solo-cb-whitelist $F/whitelist_737K.txt \
    --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 10 \
    --gene-annotation $D/genes.ht2gm --gene-strand Reverse \
    --solo-out-dir $F/out_5p_$1 \
    -p 14 --no-unal -S /dev/stdout 2> $F/run_5p_$1.log \
    | samtools sort -@ 4 -m 1G -T $F/srt_$1 -o $F/5p_$1.bam -
  samtools index -@ 4 $F/5p_$1.bam
  echo "$1 done" >> $F/run5p.status
}

run graph  $D/grch38_snp/genome_snp
run linear $D/grch38/genome
echo done > $F/run5p.done
