#!/bin/bash
# X-inactivation escape under variant-aware alignment.
#
# Escape is detected as biallelic expression at heterozygous sites on chrX, i.e.
# an allele-fraction measurement -- exactly the quantity P1 showed a linear
# reference distorts by 0.031 (CEU) to 0.048 (YRI) per site. If the INACTIVE X
# carries the ALT allele, reference bias suppresses the very signal that proves
# escape, so the gene reads as silenced. Which X is inactivated is random per
# individual, so this is not a constant offset that calibration could remove.
#
# Four arms. The last two are the same index family differing ONLY in chrX
# content, which is the internal control: arm 4 should gain on chrX and be
# identical to arm 3 on the autosomes. A gain on the autosomes would mean the
# experiment is wrong, not the biology.
set -u
L=/Users/iandriver/Downloads/p1_lcl
X=/Volumes/IanSSD/xval
R=/Users/iandriver/Downloads/ref_idx
H=/Users/iandriver/Downloads/hisat2
NP=${NP:-10}
mkdir -p $X/out

for s in ${@:-GM12878 GM18502}; do
  case $s in
    GM12878) SITES=$X/na12878.chrX.het.tsv ;;
    GM18502) SITES=$X/na18502.chrX.het.tsv ;;
  esac
  [ -s "$SITES" ] || { echo "$s: no truth file"; continue; }
  for a in linear snp_invented snp_real_Xbroken snp_real_Xfixed; do
    case $a in
      linear)            IDX=$R/grch38/genome                         ; PFX=none ;;
      snp_invented)      IDX=$R/grch38_snp/genome_snp                 ; PFX=none ;;
      snp_real_Xbroken)  IDX=/Users/iandriver/Downloads/wg64_idx/genome ; PFX=chr ;;
      snp_real_Xfixed)   IDX=/Volumes/IanSSD/ht2wg/out_x/genome       ; PFX=chr ;;
    esac
    OUT=$X/out/cnt_${s}_${a}.tsv
    [ -s $OUT ] && { echo "skip $s $a"; continue; }

    POS=$X/pos_${s}_${a}.txt
    if [ "$PFX" = chr ]; then awk -v OFS='\t' '{print $1,$2}' $SITES > $POS
    else                      awk -v OFS='\t' '{sub(/^chr/,"",$1); print $1,$2}' $SITES > $POS
    fi

    N=$(gzip -dc $L/${s}_R2.fq.gz | wc -l)
    if [ $((N/4)) -ne 50000000 ]; then
      echo "$s: only $((N/4)) reads, expected 50000000 -- skipping" >> $X/out/status; break
    fi

    /usr/bin/time -l $H/hisat2 -x $IDX -U $L/${s}_R2.fq.gz -p $NP --no-unal \
        -S /dev/stdout 2> $X/out/aln_${s}_${a}.log \
      | samtools sort -@ 3 -m 1G -T $X/srt_${s}_${a} -o $X/${s}_${a}.bam -
    samtools index -@ 3 $X/${s}_${a}.bam

    # -q 60 keeps uniquely-mapped reads only: a multi-mapper at a het site lets
    # mapping ambiguity masquerade as allelic imbalance, which is the very thing
    # being measured.
    samtools mpileup -B -q 60 -Q 20 -d 100000 -l $POS $X/${s}_${a}.bam \
        2> $X/out/mp_${s}_${a}.err \
      | python3 /Users/iandriver/Downloads/hisat2-solo-analysis/analysis/p1/count_bases.py > $OUT

    rm -f $X/${s}_${a}.bam $X/${s}_${a}.bam.bai $POS
    echo "$s $a done: $(wc -l < $OUT) covered chrX sites" >> $X/out/status
  done
done
echo ESCAPE_ALL_DONE >> $X/out/status
