#!/usr/bin/env bash
#
# Byte-verify the linearFM (no-variant) path against hisat2-build.
#
# This case had no coverage at all: every entry in mkfixtures.sh, all eight
# rung-3 fixtures and wg64_idx are built WITH a .snp/.haplotype, so the whole
# validation history exercised only the graph encoding. HISAT2 silently picks a
# different one when there are no variants (gfm.h:149,
# `_linearFM = (len + 1 == gbwtLen || gbwtLen == 0)`), and ht2wg emitted the
# graph form regardless -- producing an index hisat2-inspect could not read.
#
# NOTE the comparator: `hisat2-build` with an EMPTY --snp is not the same as
# `hisat2-build` with no --snp at all (the two differ in .1 and .5). The empty
# form is the one ht2wg is reproducing, so compare against that or the result
# looks wrong when it is right.
#
# usage: verify_linear.sh [hisat2_source_dir]
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC=${1:-}
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
BIN="$HERE/ht2fmt/target/release"
REF="$SRC/example/reference/22_20-21M.fa"
[ -f "$REF" ] || { echo "no reference at $REF"; exit 1; }

( cd "$HERE/ht2fmt" && cargo build --release --bin ht2wg 2>&1 | grep -E '^error' ) && exit 1

W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
: > "$W/empty.snp"; : > "$W/empty.haplotype"
mkdir -p "$W/tmp"

"$SRC/hisat2-build" -p 4 --snp "$W/empty.snp" --haplotype "$W/empty.haplotype" \
  "$REF" "$W/cpp" > "$W/cpp.log" 2>&1 || { echo "hisat2-build failed"; tail -5 "$W/cpp.log"; exit 1; }
HT2_THREADS=4 "$BIN/ht2wg" "$REF" "$W/empty.snp" "$W/empty.haplotype" \
  "$W/tmp" "$W/rust" 134217728 > "$W/rust.log" 2>&1 || { echo "ht2wg failed"; tail -8 "$W/rust.log"; exit 1; }

fail=0
for i in 1 2 3 4 5 6 7 8; do
  if cmp -s "$W/cpp.$i.ht2" "$W/rust.$i.ht2"; then
    printf '  .%s BYTE-IDENTICAL (%s bytes)\n' "$i" "$(wc -c < "$W/cpp.$i.ht2" | tr -d ' ')"
  else
    printf '  .%s DIFFERS\n' "$i"; fail=1
  fi
done

# The index must also be readable: the original symptom was an assertion in
# hisat2-inspect, not a byte diff.
if "$SRC/hisat2-inspect" -s "$W/rust" > "$W/insp.txt" 2>&1; then
  echo "  hisat2-inspect OK ($(grep -c '^Sequence' "$W/insp.txt") sequence(s))"
else
  echo "  hisat2-inspect FAILED to read the index"; fail=1
fi

[ $fail -eq 0 ] && echo "LINEAR PATH BYTE-IDENTICAL" || echo "LINEAR PATH DIFFERS"
exit $fail
