#!/usr/bin/env bash
#
# Regenerate the rung-3 test fixtures from current source.
#
# Every fixture input is *derived* here from files tracked in a git repo -- the
# HISAT2 example reference and its variant list, plus two small hand-written
# cases -- so the whole set is reproducible from scratch. Nothing is copied
# from a previous fixture directory.
#
# Why this script exists
# ----------------------
# The first fixture set was built ad hoc by a `hisat2-build-s` binary that was
# months older than the sources beside it. `git status` was clean, both fixture
# builds agreed with each other, and nothing looked wrong -- but the binary did
# not match the source anyone was reading, so a correct emitter appeared to
# have a bug for two rounds of debugging. This script therefore *rebuilds the
# binary first* and records its hash next to the fixtures.
#
# usage: mkfixtures.sh <outdir> [hisat2_source_dir]

set -euo pipefail

OUT=${1:?usage: mkfixtures.sh <outdir> [hisat2_source_dir]}
SRC=${2:-/Users/iandriver/Downloads/hisat2}
SRC=$(cd "$SRC" && pwd)

REF="$SRC/example/reference/22_20-21M.fa"
SNP="$SRC/example/reference/22_20-21M.snp"
TESTDATA="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/testdata"

for f in "$REF" "$SNP" "$TESTDATA/mkmulti.py"; do
  [ -r "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done

# --- 1. rebuild the builder from current source ------------------------------
echo "== rebuilding hisat2-build-s from $SRC"
make -C "$SRC" hisat2-build-s >/dev/null
BUILDER="$SRC/hisat2-build-s"

# Guard against the failure this script was written for: a binary older than
# the sources it was supposedly built from.
newer=$(find "$SRC" -maxdepth 1 \( -name '*.cpp' -o -name '*.h' \) -newer "$BUILDER" -print -quit)
if [ -n "$newer" ]; then
  echo "hisat2-build-s is older than $newer -- refusing to build stale fixtures" >&2
  exit 1
fi

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)

# --- 2. derive every fixture input -------------------------------------------
echo "== deriving fixture inputs into $OUT"
REF="$REF" SNP="$SNP" OUT="$OUT" TESTDATA="$TESTDATA" python3 - <<'PY'
import os

ref_path, snp_path, out = os.environ['REF'], os.environ['SNP'], os.environ['OUT']
J = lambda *p: os.path.join(out, *p)

def read_fa(path):
    recs, name, buf = [], None, []
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name is not None: recs.append((name, ''.join(buf)))
            name, buf = line[1:], []
        else:
            buf.append(line)
    if name is not None: recs.append((name, ''.join(buf)))
    return recs

def write_fa(path, recs, width=60):
    with open(path, 'w') as f:
        for name, s in recs:
            f.write('>%s\n' % name)
            for i in range(0, len(s), width):
                f.write(s[i:i+width] + '\n')

def haplotypes(snps, numbering=None):
    """One single-variant haplotype per variant, as `hisat2-build` reads them.

    A deletion spans [pos, pos+len-1]; everything else is a point. `numbering`
    is a parallel list of ht indices, letting a filtered subset keep the numbers
    it had in the unfiltered file -- which is what the type-restricted fixtures
    do. It must be positional, not keyed by rsID: 32 of the 1,881 clean.snp
    entries are multi-allelic and share an rsID with a sibling record."""
    out = []
    for i, r in enumerate(snps):
        pos = int(r[3])
        right = pos + int(r[4]) - 1 if r[1] == 'deletion' else pos
        idx = numbering[i] if numbering is not None else i
        out.append(['ht%d' % idx, r[2], str(pos), str(right), r[0]])
    return out

def write_tsv(path, rows):
    with open(path, 'w') as f:
        for r in rows:
            f.write('\t'.join(r) + '\n')

# -- ex: the example reference unchanged, with its own variant list -----------
ex_snps = [l.rstrip('\n').split('\t') for l in open(snp_path)]
write_tsv(J('ex.snp'), ex_snps)
write_tsv(J('ex.haplotype'), haplotypes(ex_snps))

# -- clean: the example reference truncated at its first N --------------------
# 509,431 bp with no ambiguity at all, so window boundaries and rstarts stay
# trivial -- the control against which the N-run reference is read.
full = read_fa(ref_path)[0][1]
cut = full.index('N')
write_fa(J('clean.fa'), [('chrC', full[:cut])])
clean_snps = [[r[0], r[1], 'chrC', r[3], r[4]] for r in ex_snps if int(r[3]) < cut]
write_tsv(J('clean.snp'), clean_snps)
clean_haps = haplotypes(clean_snps)
write_tsv(J('clean.haplotype'), clean_haps)

# -- t_single / t_insertion / t_deletion: one variant type at a time ----------
# Each keeps its ht numbering from clean.haplotype, so a diff against the
# unfiltered fixture shows only removals.
for typ in ('single', 'insertion', 'deletion'):
    keep = [(i, r) for i, r in enumerate(clean_snps) if r[1] == typ]
    write_tsv(J('t_%s.snp' % typ), [r for _, r in keep])
    write_tsv(J('t_%s.haplotype' % typ),
              haplotypes([r for _, r in keep], numbering=[i for i, _ in keep]))

# -- tiny / x: hand-written minimal cases ------------------------------------
# 200 bp, checked in verbatim so the smallest fixture never drifts
tiny = ('TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACG'
        'GAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC'
        'CCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAA'
        'TGACGGATATATATTAAAAA')
assert len(tiny) == 200
write_fa(J('tiny.fa'), [('chrT', tiny)])
write_tsv(J('tiny.snp'), [['rs0', 'single', 'chrT', '50', 'A'],
                          ['rs1', 'single', 'chrT', '80', 'A'],
                          ['rs2', 'single', 'chrT', '120', 'A']])
write_tsv(J('tiny.haplotype'), [['ht0', 'chrT', '50', '80', 'rs0,rs1'],
                                ['ht1', 'chrT', '120', '120', 'rs2']])
# a lone deletion -- the one variant type the SNP-only cases never reach
write_tsv(J('x.snp'), [['v', 'deletion', 'chrT', '50', '3']])
write_tsv(J('x.haplotype'), [['h', 'chrT', '50', '52', 'v']])

# -- multi: five sequences, each breaking a different reference-reading rule --
import subprocess, sys
subprocess.run([sys.executable, os.path.join(os.environ['TESTDATA'], 'mkmulti.py')],
               cwd=out, check=True)
multi = read_fa(J('multi.fa'))
m_snps, m_haps = [], []
for name, s in multi:
    short = name.split()[0]
    for pos in range(300, len(s) - 300, 250):
        ref = s[pos].upper()
        if ref == 'N': continue
        m_snps.append(['rs%d' % len(m_snps), 'single', short, str(pos),
                       'C' if ref == 'A' else 'A'])
write_tsv(J('m.snp'), m_snps)
write_tsv(J('m.haplotype'), haplotypes(m_snps))
PY

# --- 3. build every index ----------------------------------------------------
build () {  # build <prefix> <fasta> <snp> <haplotype>
  echo "== $1"
  ( cd "$OUT" && "$BUILDER" --snp "$3" --haplotype "$4" "$2" "$1" > "$1.log" 2>&1 )
}

build cleanidx    clean.fa clean.snp        clean.haplotype
build t_single    clean.fa t_single.snp     t_single.haplotype
build t_insertion clean.fa t_insertion.snp  t_insertion.haplotype
build t_deletion  clean.fa t_deletion.snp   t_deletion.haplotype
build tinyidx     tiny.fa  tiny.snp         tiny.haplotype
build xidx        tiny.fa  x.snp            x.haplotype
build multiidx    multi.fa m.snp            m.haplotype
( cd "$OUT" && echo "== exidx" && "$BUILDER" --snp ex.snp --haplotype ex.haplotype "$REF" exidx > exidx.log 2>&1 )

# --- 4. manifest -------------------------------------------------------------
{
  echo "# rung-3 fixtures"
  echo "generated-by:  $(basename "${BASH_SOURCE[0]}")"
  echo "hisat2-source: $SRC"
  echo "hisat2-commit: $(git -C "$SRC" rev-parse HEAD)"
  echo "hisat2-dirty:  $(git -C "$SRC" status --porcelain | wc -l | tr -d ' ') modified files"
  echo "builder-sha:   $(shasum -a 256 "$BUILDER" | cut -d' ' -f1)"
  echo
  ( cd "$OUT" && shasum -a 256 $(ls | grep -v '^MANIFEST') )
} > "$OUT/MANIFEST.txt"

echo "== done; manifest at $OUT/MANIFEST.txt"
