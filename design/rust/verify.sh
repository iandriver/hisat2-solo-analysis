#!/usr/bin/env bash
#
# Run every emitter against every rung-3 fixture and report byte-equality.
#
# usage: verify.sh <fixture_dir> [hisat2_source_dir]

set -uo pipefail

FIX=${1:?usage: verify.sh <fixture_dir> [hisat2_source_dir]}
SRC=${2:-/Users/iandriver/Downloads/hisat2}
FIX=$(cd "$FIX" && pwd)
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BIN="$HERE/ht2fmt/target/release"
REF="$SRC/example/reference/22_20-21M.fa"

( cd "$HERE/ht2fmt" && cargo build --release 2>&1 | grep -E '^error' ) && exit 1

# prefix : fasta : snp : haplotype
CASES=(
  "tinyidx:tiny.fa:tiny.snp:tiny.haplotype"
  "xidx:tiny.fa:x.snp:x.haplotype"
  "multiidx:multi.fa:m.snp:m.haplotype"
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

  # .5 / .6 -- the hierarchical local indexes
  out=$("$BIN/ht2emit5" "$fa" "$snp" "$hap" "$p" 2>&1)
  echo "$out" | grep -E '^\.[56]\.ht2:' | sed 's/^/  emit5 /'
  echo "$out" | grep -qE '^\.5\.ht2: BYTE-IDENTICAL' || fail=1
  echo "$out" | grep -qE '^\.6\.ht2: BYTE-IDENTICAL' || fail=1

  # .1 / .2 -- the global graph index
  out=$("$BIN/ht2path" "$fa" "$snp" "$hap" "$p.log" "$p.1.ht2" 2>&1)
  echo "$out" | grep -E 'GRAPH GFM|GENERATION CURVE|ftab:|eftab:|SA sample:' | sed 's/^ *//;s/^/  path  /'
  echo "$out" | grep -q 'GRAPH GFM MATCHES'      || fail=1
  echo "$out" | grep -q 'GENERATION CURVE MATCHES' || fail=1

  # layout parse -- every byte of .5/.6 accounted for
  "$BIN/ht2local" "$p" 2>&1 | tail -1 | sed 's/^/  local /'
done

printf '\n%s\n' "$([ $fail -eq 0 ] && echo 'ALL FIXTURES BYTE-IDENTICAL' || echo 'SOME FIXTURES DIFFER')"
exit $fail
