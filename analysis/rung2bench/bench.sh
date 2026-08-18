#!/bin/bash
# Rung 2 full-build benchmark: hisat2-build against the Rust emitter, on the
# same reference, at four sizes.
#
# The point is not that the Rust side is faster -- it uses a deliberately simple
# O(n log^2 n) prefix-doubling suffix array, chosen because it is easy to argue
# is correct, which mattered more than speed while establishing byte-identity.
# The point is to put real numbers on where a Rust builder stands today and where
# its cost actually goes, with byte-identity re-verified at every size so the
# comparison is between two builders producing THE SAME BYTES.
set -u
W=/Users/iandriver/Downloads/rung2_bench
H=/Users/iandriver/Downloads/hisat2
E=/Users/iandriver/Downloads/hisat2-solo-analysis/design/rust/ht2fmt/target/release/ht2emit
cd $W
printf "%-6s %10s  %-14s %9s %9s  %s\n" size bp builder wall_s peak_MB identical
for tag in "$@"; do
  fa=ref_$tag.fa
  [ -s "$fa" ] || { echo "missing $fa"; continue; }
  bp=$(grep -v '^>' $fa | tr -d '\n' | wc -c | tr -d ' ')

  rm -f cpp_$tag.*.ht2
  /usr/bin/time -l $H/hisat2-build -p 1 $fa cpp_$tag > /dev/null 2> t_cpp_$tag.log
  cw=$(awk '/ real /{print $1}' t_cpp_$tag.log)
  cm=$(awk '/maximum resident/{printf "%.0f", $1/1048576}' t_cpp_$tag.log)

  rm -f rs_$tag.*.ht2
  /usr/bin/time -l $E $fa rs_$tag cpp_$tag.1.ht2 > o_rs_$tag.log 2> t_rs_$tag.log
  rw=$(awk '/ real /{print $1}' t_rs_$tag.log)
  rm_=$(awk '/maximum resident/{printf "%.0f", $1/1048576}' t_rs_$tag.log)
  id=$(grep -c 'BYTE-IDENTICAL' o_rs_$tag.log)

  printf "%-6s %10s  %-14s %9s %9s  %s\n" $tag $bp hisat2-build $cw $cm "-"
  printf "%-6s %10s  %-14s %9s %9s  %s\n" ""   ""   ht2emit      $rw $rm_ "$id/4 files"
  rm -f rs_$tag.*.ht2 cpp_$tag.5.ht2 cpp_$tag.6.ht2 cpp_$tag.7.ht2 cpp_$tag.8.ht2
done
