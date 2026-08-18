#!/bin/bash
# P1 -- reference bias at heterozygous sites: linear vs invented-haplotype graph
#       vs real-phasing graph.
#
# The measurement is the alternate-allele fraction at sites where the donor is
# known to be heterozygous. With no bias it is 0.5; a linear reference pushes it
# below 0.5 because alt-carrying reads pay a mismatch penalty. H2 measured that
# deficit at 0.025 (CEU) and 0.029 (YRI) against `grch38_snp`, whose haplotypes
# are INVENTED by greedy graph colouring. This adds the third arm the whole
# 4 h 25 m whole-genome build existed to make possible: the same measurement
# against real 1000 Genomes phasing.
#
# Deliberately paired: the same reads and the same sites go through all three
# indexes, so every site contributes a within-site difference and the comparison
# does not depend on which sites happen to be covered.
#
# BAMs are deleted as soon as their counts are extracted; disk is the binding
# constraint, not compute.
set -u
L=/Users/iandriver/Downloads/p1_lcl
T=/Users/iandriver/Downloads/p1_truth
R=/Users/iandriver/Downloads/ref_idx
H=/Users/iandriver/Downloads/hisat2
NP=${NP:-10}

mkdir -p $L/out
for s in GM12878 GM18502; do
  case $s in
    GM12878) SITES=$T/na12878.het.tsv ;;
    GM18502) SITES=$T/na18502.het.tsv ;;
  esac
  for a in linear snp_invented snp_real; do
    case $a in
      linear)       IDX=$R/grch38/genome            ; PFX=none ;;
      snp_invented) IDX=$R/grch38_snp/genome_snp    ; PFX=none ;;
      snp_real)     IDX=/Users/iandriver/Downloads/wg64_idx/genome ; PFX=chr ;;
    esac
    OUT=$L/out/cnt_${s}_${a}.tsv
    [ -s $OUT ] && { echo "skip $s $a"; continue; }

    # The JHU indexes name contigs the Ensembl way (1, 2, ...) and ours the
    # GENCODE way (chr1, chr2, ...). Same coordinates, different RNAME, so the
    # position list has to be written in the naming the BAM will use.
    POS=$L/pos_${s}_${a}.txt
    if [ "$PFX" = chr ]; then awk -v OFS='\t' '{print $1,$2}' $SITES > $POS
    else                      awk -v OFS='\t' '{sub(/^chr/,"",$1); print $1,$2}' $SITES > $POS
    fi

    /usr/bin/time -l $H/hisat2 -x $IDX -U $L/${s}_R2.fq.gz -p $NP --no-unal \
        -S /dev/stdout 2> $L/out/aln_${s}_${a}.log \
      | samtools sort -@ 3 -m 1G -T $L/srt_${s}_${a} -o $L/${s}_${a}.bam -
    samtools index -@ 3 $L/${s}_${a}.bam

    # -q 60: uniquely-mapped only. Multi-mappers at a het site would let mapping
    # ambiguity masquerade as allelic imbalance, which is the thing being measured.
    samtools mpileup -B -q 60 -Q 20 -d 100000 -l $POS $L/${s}_${a}.bam \
        2> $L/out/mp_${s}_${a}.err \
      | python3 $(dirname $0)/count_bases.py > $OUT

    rm -f $L/${s}_${a}.bam $L/${s}_${a}.bam.bai $POS
    echo "$s $a done: $(wc -l < $OUT) covered sites" >> $L/out/status
  done
done
echo ALL_DONE >> $L/out/status
