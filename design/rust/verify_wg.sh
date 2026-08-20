#!/usr/bin/env bash
#
# Build every rung-3 fixture through the external/fragmented path and byte-diff
# it against the index `hisat2-build` wrote.
#
# HT2_LARGE=1 does the same against the 64-bit fixture set.
#
# HT2_SMALL=1 runs only the three sub-16 kb fixtures, which is what makes a
# stressed configuration usable -- `HT2_SEG=7 verify_wg.sh <dir> "" 13` forces
# dozens of segments and real spilling on a 200 bp graph, and would take hours
# on the 900 kb ones.
#
# usage: verify_wg.sh <fixture_dir> [hisat2_source_dir] [budget_records]

set -uo pipefail

FIX=$(cd "${1:?usage: verify_wg.sh <fixture_dir> [hisat2_source_dir] [budget]}" && pwd)
SRC=${2:-}; [ -z "$SRC" ] && SRC=/Users/iandriver/Downloads/hisat2
BUDGET=${3:-200000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BIN="$HERE/ht2fmt/target/release"
REF="$SRC/example/reference/22_20-21M.fa"
WORK=${HT2_WORK:-/tmp/ht2wg}

( cd "$HERE/ht2fmt" && cargo build --release 2>&1 | grep -E '^error' ) && exit 1
EXT=ht2; [ -n "${HT2_LARGE:-}" ] && EXT=ht2l

CASES=(
  "tinyidx:tiny.fa:tiny.snp:tiny.haplotype"
  "xidx:tiny.fa:x.snp:x.haplotype"
  "multiidx:multi.fa:m.snp:m.haplotype"
)
[ -z "${HT2_SMALL:-}" ] && CASES+=(
  "cleanidx:clean.fa:clean.snp:clean.haplotype"
  "t_single:clean.fa:t_single.snp:t_single.haplotype"
  "t_insertion:clean.fa:t_insertion.snp:t_insertion.haplotype"
  "t_deletion:clean.fa:t_deletion.snp:t_deletion.haplotype"
  "exidx:$REF:ex.snp:ex.haplotype"
)

fail=0
cd "$FIX"
for c in "${CASES[@]}"; do
  IFS=: read -r p fa snp hap <<< "$c"
  printf '\n=== %s ===\n' "$p"
  rm -rf "$WORK" "$WORK.out"*
  out=$(/usr/bin/time -l "$BIN/ht2wg" "$fa" "$snp" "$hap" "$WORK" "$WORK.out" "$BUDGET" "$p" 2>&1)
  echo "$out" | grep -qE 'EXTERNAL BUILD MATCHES' || fail=1
  echo "$out" | grep -E 'ftab:|BYTE-IDENTICAL|differing|CURVE|WARNING|NOTE|maximum resident' \
    | sed 's/^ *//;s/^/  /'
  rm -rf "$WORK" "$WORK.out"*
done

printf '\n%s\n' "$([ $fail -eq 0 ] && echo "EXTERNAL PATH REPRODUCES EVERY FIXTURE ($EXT)" || echo 'SOME FIXTURES DIFFER')"
exit $fail
