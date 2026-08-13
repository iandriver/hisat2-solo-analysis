#!/bin/bash
# T4b: genotype calling from scRNA-seq reads, graph vs linear alignment.
#
# GM12878 (NA12878, CEU) has a GIAB v4.2.1 truth set; GM18502 (NA18502, YRI)
# gives the ancestry contrast. Only the cDNA mate is used -- genotyping here is
# pseudobulk, so barcodes are irrelevant.
#
# BAMs are deleted as soon as their allele counts are extracted; disk is tight.
set -u
L=/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/lcl
D=/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human
H=/Users/iandriver/Downloads/hisat2

for s in GM12878 GM18502; do
  for a in graph linear; do
    [ "$a" = graph ] && IDX=$D/grch38_snp/genome_snp || IDX=$D/grch38/genome
    [ -s $L/cnt_${s}_${a}.tsv ] && continue
    /usr/bin/time -l $H/hisat2 -x $IDX -U $L/${s}_R2.fq.gz -p 14 --no-unal \
      -S /dev/stdout 2> $L/aln_${s}_${a}.log \
      | samtools sort -@ 4 -m 1G -T $L/srt_${s}_${a} -o $L/${s}_${a}.bam -
    samtools index -@ 4 $L/${s}_${a}.bam
    # Pile up at all common dbSNP SNVs, not just the truth-VCF variants, so that
    # truth hom-ref sites are represented too and the error rate is comparable
    # to a published whole-genotype figure rather than variant-sites-only.
    samtools mpileup -B -q 60 -Q 20 -d 100000 -l $D/snv_sites.txt $L/${s}_${a}.bam \
      2> $L/mp_${s}_${a}.err | python3 $D/count_bases.py > $L/cnt_${s}_${a}.tsv
    rm -f $L/${s}_${a}.bam $L/${s}_${a}.bam.bai
    echo "$s $a done $(wc -l < $L/cnt_${s}_${a}.tsv) sites" >> $L/lcl.status
  done
done
echo done > $L/lcl.done
