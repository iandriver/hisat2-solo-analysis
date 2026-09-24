#!/usr/bin/env bash
#
# Crash the build at every checkpoint boundary, in every generation, and check
# that finishing it produces the same index as never having been interrupted.
#
# Racing a `kill -9` finds bugs but cannot say which boundary it found them at,
# and cannot prove a boundary was ever reached. `HT2_CRASH_AT` names one and
# `HT2_CRASH_GEN` narrows it to a generation, so every place a restart can land
# is visited on purpose. A crash point that is never reached is skipped and
# counted, so the total is a claim about coverage rather than a hope.
#
# usage: verify_resume.sh <fixture_dir> [hisat2_source_dir]
#        HT2_LARGE=1 does the same against the 64-bit fixtures

set -uo pipefail

FIX=$(cd "${1:?usage: verify_resume.sh <fixture_dir> [hisat2_source_dir]}" && pwd)
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC=${2:-}
if [ -z "${SRC:-}" ]; then
    for c in "$HERE/../../../hisat2" "$HERE/../../hisat2" "$PWD/../hisat2"; do
        [ -d "$c/example/reference" ] && SRC=$(cd "$c" && pwd) && break
    done
fi
if [ -z "${SRC:-}" ] || [ ! -d "$SRC/example/reference" ]; then
    echo "need a HISAT2 checkout: pass it as an argument, or put one beside this repo" >&2
    exit 2
fi
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
B="$HERE/ht2fmt/target/release/ht2wg"
REF_FA="$SRC/example/reference/22_20-21M.fa"
EXT=ht2; [ -n "${HT2_LARGE:-}" ] && EXT=ht2l
W=${HT2_WORK:-/tmp/ht2resume}; OUT=$W.out

( cd "$HERE/ht2fmt" && cargo build --release 2>&1 | grep -E '^error' ) && exit 1

CASES=(
  "cleanidx:clean.fa:clean.snp:clean.haplotype"
  "exidx:$REF_FA:ex.snp:ex.haplotype"
)
[ -n "${HT2_QUICK:-}" ] && CASES=("cleanidx:clean.fa:clean.snp:clean.haplotype")

cd "$FIX"
total=0; fails=0; skipped=0
for c in "${CASES[@]}"; do
  IFS=: read -r p fa snp hap <<< "$c"
  printf '\n=== %s ===\n' "$p"
  for at in graph sorted joined keyed gen edges; do
    for g in 1 2 3 4 5 6 7 8 9 10 11 12; do
      case "$at" in graph|edges) [ "$g" != 1 ] && continue;; esac
      rm -rf "$W"; for n in 1 2 3 4 5 6 7 8; do rm -f "$OUT.$n.$EXT"; done
      HT2_CRASH_AT=$at HT2_CRASH_GEN=$g "$B" "$fa" "$snp" "$hap" "$W" "$OUT" 5000 "$p" \
        >/tmp/resume1.log 2>&1
      [ $? -ne 70 ] && { skipped=$((skipped+1)); continue; }   # never reached
      total=$((total+1))
      "$B" "$fa" "$snp" "$hap" "$W" "$OUT" 5000 "$p" >/tmp/resume2.log 2>&1
      if ! grep -q 'EXTERNAL BUILD MATCHES' /tmp/resume2.log; then
        fails=$((fails+1))
        echo "  FAIL at=$at gen=$g"
        grep -E 'differing|panicked|Error|resuming' /tmp/resume2.log | head -3 | sed 's/^/      /'
      fi
    done
  done
done
rm -rf "$W"
printf '\n%d crash points exercised, %d unreachable, %d failures\n' "$total" "$skipped" "$fails"
[ "$fails" -eq 0 ] && echo "EVERY INTERRUPTED BUILD FINISHED BYTE-IDENTICAL ($EXT)"
exit "$fails"
