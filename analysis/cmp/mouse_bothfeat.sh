#!/bin/bash
# Close the last stale number in the three-way report: the both-feature runtime
# on the actual mouse dataset, rather than inferred from a 111,600-read fixture.
#
# Same input as the original benchmark -- 5k_Mouse_PBMCs_5p_gem-x, first 10M
# reads, the same .ht2gm gene model and 3.7M whitelist, 5' chemistry
# (CB 1-16, UMI 17-28, Reverse strand), 16 threads.
#
# Only the HISAT2 index has to be rebuilt; everything else comes from S3. The
# .ht2gm carries #ref name/length lines and the loader fails loudly on a
# mismatch, so a wrong GENCODE release cannot silently produce plausible counts.
set -euo pipefail
W=${1:-/private/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/cmp/mouse}
H=/Users/iandriver/Downloads/hisat2
NREADS=10000000
THREADS=16
mkdir -p "$W" && cd "$W"

say() { printf '%s  %s\n' "$(date -u +%H:%M:%S)" "$*"; }

# ---- inputs -------------------------------------------------------------
[ -s R1.fq.gz ] || aws s3 cp s3://rustar-bench/fastq/5k_Mouse_PBMCs_5p_gem-x_GEX_S1_L001_R1_001.fastq.gz R1.fq.gz --quiet
[ -s R2.fq.gz ] || aws s3 cp s3://rustar-bench/fastq/5k_Mouse_PBMCs_5p_gem-x_GEX_S1_L001_R2_001.fastq.gz R2.fq.gz --quiet
# NOT hisat2-bench/refs/genes.ht2gm -- that one is human (ref 1 = 248,956,422 bp).
# The mouse model is regenerated from the GENCODE M33 GTF that the benchmark used.
[ -s mouse_genes.gtf.gz ] || aws s3 cp s3://rustar-bench/mouse_genes.gtf.gz . --quiet
# The 3' and 5' GEM-X whitelists are both 3,686,400 barcodes and 62,668,800
# bytes, so size proves nothing. hisat2-bench/refs/whitelist_v3.txt is the 3'
# list and matches 1 of 902 observed barcodes here; the benchmark's own copy
# matches 733. Take the one the original run used.
[ -s whitelist.txt ] || aws s3 cp s3://rustar-bench/results/3way-20260706T145112Z/whitelist.txt whitelist.txt --quiet
say "inputs ready ($(du -sh . | cut -f1))"

# 10M reads to match the original run, not the whole file.
if [ ! -s sub_R1.fq ] || [ ! -s sub_R2.fq ]; then
  gzip -dc R1.fq.gz | head -n $((NREADS * 4)) > sub_R1.fq
  gzip -dc R2.fq.gz | head -n $((NREADS * 4)) > sub_R2.fq
fi
say "subset: $(( $(wc -l < sub_R1.fq) / 4 )) reads"

# ---- reference ----------------------------------------------------------
# gzip -t, not -s: EBI truncated this twice and a partial .gz passes a size test.
if [ ! -s genome.fa ]; then
  gzip -t genome.fa.gz || { echo "genome.fa.gz incomplete"; exit 1; }
  gzip -dc genome.fa.gz > genome.fa
fi
[ -s genome.fa.fai ] || samtools faidx genome.fa
[ -s genes.ht2gm ] || python3 "$H/hisat2_extract_genes.py" --fai genome.fa.fai \
    -o genes.ht2gm <(gzip -dc mouse_genes.gtf.gz)
say "genome: $(grep -c '^>' genome.fa) contigs"

# The gene model's #ref lines must match the genome, or the run is meaningless.
python3 - <<'PY'
import re, sys
lens = {}
for line in open('genome.fa'):
    if line.startswith('>'):
        cur = line[1:].split()[0]; lens[cur] = 0
    else:
        lens[cur] += len(line.strip())
bad = 0; n = 0
for line in open('genes.ht2gm'):
    if not line.startswith('#ref'): continue
    _, name, ln = line.split('\t')[:3]
    n += 1
    if lens.get(name) != int(ln):
        bad += 1
        if bad < 4: print("  mismatch %s: model=%s genome=%s" % (name, ln.strip(), lens.get(name)))
print("ref check: %d refs, %d mismatched" % (n, bad))
sys.exit(1 if bad else 0)
PY
say "gene model matches the genome"

if [ ! -s idx/mouse.1.ht2 ]; then
  mkdir -p idx
  /usr/bin/time -l "$H/hisat2-build" -p "$THREADS" genome.fa idx/mouse > build.log 2> build.time
fi
say "index built ($(du -sh idx | cut -f1))"

# ---- the measurement ----------------------------------------------------
run () {  # run <label> <feature-spec>
  local lab=$1 spec=$2
  rm -rf "out_$lab"
  local s=$(python3 -c 'import time;print(time.time())')
  "$H/hisat2" -x idx/mouse -1 sub_R1.fq -2 sub_R2.fq \
      --solo-barcode-mate 1 --solo-cb-whitelist whitelist.txt \
      --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 12 \
      --gene-annotation genes.ht2gm --gene-strand Reverse \
      --gene-feature "$spec" --solo-out-dir "out_$lab" \
      -p "$THREADS" --no-unal -S /dev/null > "log_$lab.txt" 2>&1
  local e=$(python3 -c 'import time;print(time.time())')
  python3 -c "print(f'{$e-$s:.1f}')"
}

# The machine is shared, so run the configs in a palindrome and take the best of
# each: monotone load drift cancels, and the minimum is the least-contaminated
# sample rather than an average of interference.
g1=$(run gene Gene);       f1=$(run full GeneFull);  b1=$(run both Gene,GeneFull)
b2=$(run both2 Gene,GeneFull); f2=$(run full2 GeneFull); g2=$(run gene2 Gene)
t_gene=$(python3 -c "print(min($g1,$g2))")
t_full=$(python3 -c "print(min($f1,$f2))")
t_both=$(python3 -c "print(min($b1,$b2))")
say "Gene only:      ${g1}s / ${g2}s -> ${t_gene}s"
say "GeneFull only:  ${f1}s / ${f2}s -> ${t_full}s"
say "Gene,GeneFull:  ${b1}s / ${b2}s -> ${t_both}s"

python3 - "$t_gene" "$t_full" "$t_both" <<'PY'
import sys
g, f, b = (float(x) for x in sys.argv[1:4])
print()
print("two passes (Gene then GeneFull): %.1fs" % (g + f))
print("one pass  (Gene,GeneFull):       %.1fs" % b)
print("speedup:                          %.2fx" % ((g + f) / b))
print("overhead of the second feature:   %+.1fs (%.1f%% over Gene alone)"
      % (b - g, 100 * (b - g) / g))
PY

# Same invariant the test suite asserts, on real data this time.
for f in gene full both; do :; done
gsum=$(awk 'NR>3{s+=$3} END{print s+0}' out_both/Gene/raw/matrix.mtx)
fsum=$(awk 'NR>3{s+=$3} END{print s+0}' out_both/GeneFull/raw/matrix.mtx)
echo "Gene UMIs=$gsum  GeneFull UMIs=$fsum"
cmp -s out_gene/Gene/raw/matrix.mtx out_both/Gene/raw/matrix.mtx \
  && echo "Gene matrix identical to the single-feature run" \
  || echo "WARNING: Gene matrix differs from the single-feature run"
cmp -s out_full/GeneFull/raw/matrix.mtx out_both/GeneFull/raw/matrix.mtx \
  && echo "GeneFull matrix identical to the single-feature run" \
  || echo "WARNING: GeneFull matrix differs from the single-feature run"
say "DONE"
