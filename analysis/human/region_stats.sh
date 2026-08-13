#!/bin/bash
# Per-region alignment statistics from both BAMs: depth, uniqueness, mismatches.
D=/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human
cd $D
printf "%-12s %-24s %8s %8s %8s %7s %7s\n" gene region reads mapq60 mapq1 meanNM pctNM0
while read -r gene chrom start end; do
  for b in graph linear; do
    stats=$(samtools view -F 0x100 $b.sorted.bam "$chrom:$start-$end" | awk '
      { n++; if($5==60) q60++; else if($5==1) q1++;
        nm=-1; for(i=12;i<=NF;i++) if($i ~ /^NM:i:/){ split($i,a,":"); nm=a[3] }
        if(nm>=0){ sn+=nm; sc++; if(nm==0) z++ } }
      END{ printf "%d %d %d %.3f %.1f", n, q60, q1, (sc?sn/sc:0), (sc?100*z/sc:0) }')
    printf "%-12s %-24s %8s %8s %8s %7s %7s   (%s)\n" "$gene" "$chrom:$start-$end" $stats "$b"
  done
done
