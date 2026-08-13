#!/bin/bash
# trace.sh <from> <to> <chrom> <start> <end> <label>
# Take reads that <from> places uniquely (MAPQ 60) in the region, and report
# where <to> puts the same reads. Recovery looks like "unaligned in <to>";
# reallocation looks like "confidently placed elsewhere in <to>".
D=/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human
from=$1; to=$2; chrom=$3; start=$4; end=$5; label=$6
tmp=$(mktemp)
samtools view -F 0x100 -q 60 $D/$from.sorted.bam "$chrom:$start-$end" | awk '{print $1}' | sort -u > $tmp
n=$(wc -l < $tmp)
samtools view -F 0x100 $D/$to.sorted.bam | awk -v T="$tmp" '
  BEGIN{ while((getline l < T) > 0) w[l] }
  ($1 in w){ if($3=="*") un++; else if($5==60) q60++; else amb++ }
  END{ printf "%d %d %d\n", un, q60, amb }' > $tmp.res
read un q60 amb < $tmp.res
printf "%-34s %8d reads unique in %-6s ->  in %-6s: unaligned %6d (%5.1f%%)  unique %6d (%5.1f%%)  multi %6d (%5.1f%%)\n" \
  "$label" $n "$from" "$to" $un $(echo "100*$un/$n"|bc -l) $q60 $(echo "100*$q60/$n"|bc -l) $amb $(echo "100*$amb/$n"|bc -l)
rm -f $tmp $tmp.res
