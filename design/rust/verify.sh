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

  # A 64-bit fixture set is `.ht2l`. Hardcoding `.ht2` here made every check
  # below silently miss on rung3l: the greps never matched, so `fail` was set
  # unconditionally and the sweep reported DIFFER while every emitter passed.
  ext=ht2; [ -f "$p.1.ht2l" ] && ext=ht2l

  # .5 / .6 -- the hierarchical local indexes
  # ht2emit5 picks its width from HT2_LARGE, not from the files on disk, so a
  # 64-bit fixture set needs it set or the emitter builds a 32-bit index and
  # then looks for a `.5.ht2` that is not there.
  [ "$ext" = ht2l ] && export HT2_LARGE=1 || unset HT2_LARGE
  out=$("$BIN/ht2emit5" "$fa" "$snp" "$hap" "$p" 2>&1)
  echo "$out" | grep -E "^\.[56]\.$ext:" | sed 's/^/  emit5 /'
  echo "$out" | grep -qE "^\.5\.$ext: BYTE-IDENTICAL" || fail=1
  echo "$out" | grep -qE "^\.6\.$ext: BYTE-IDENTICAL" || fail=1

  # .1 / .2 -- the global graph index.
  # ht2path reads index_t as 4 bytes throughout and has no HT2_LARGE concept, so
  # it cannot check a `.ht2l`. Say so rather than counting it as a failure: the
  # 64-bit cover for .1/.2 is `HT2_LARGE=1 verify_wg.sh`, which builds all eight
  # files end to end through the same emitters this sweep exercises piecewise.
  if [ "$ext" = ht2l ]; then
    echo "  path  SKIPPED (ht2path is 32-bit only; see HT2_LARGE=1 verify_wg.sh)"
  else
    out=$("$BIN/ht2path" "$fa" "$snp" "$hap" "$p.log" "$p.1.$ext" 2>&1)
    echo "$out" | grep -E 'GRAPH GFM|GENERATION CURVE|ftab:|eftab:|SA sample:' | sed 's/^ *//;s/^/  path  /'
    echo "$out" | grep -q 'GRAPH GFM MATCHES'      || fail=1
    echo "$out" | grep -q 'GENERATION CURVE MATCHES' || fail=1
  fi

  # layout parse -- every byte of .5/.6 accounted for. This one never fed into
  # `fail`, so a layout regression would have printed and been ignored.
  out=$("$BIN/ht2local" "$p" 2>&1)
  echo "$out" | tail -1 | sed 's/^/  local /'
  echo "$out" | grep -q 'LOCAL INDEX LAYOUT OK' || fail=1
done

printf '\n%s\n' "$([ $fail -eq 0 ] && echo 'ALL FIXTURES BYTE-IDENTICAL' || echo 'SOME FIXTURES DIFFER')"
exit $fail
