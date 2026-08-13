#!/bin/bash
# Phase 2b -- whole-genome build. GENCODE v50 primary assembly, 1000 Genomes
# 30x phased panel (3202 samples), MAF >= 1%.
#
# Why this variant set, from the 2a chr1 probe:
#   retention   1000G phased 92.4% > dbSNP 1000G_30X 86.7% > dbSNP big-4 80.9%
#   and per 1000 variant instances the edge budget is blown at the same rate by
#   real haplotypes as by invented ones (0.392 vs 0.400), so haplotype realism is
#   close to free and density is what costs.
# With the binary-search backoff now in hgfm.h, chr1 retention is 97.8%.
#
# Parallel per chromosome, then concatenated. Safe because the .haplotype ID
# field is read only to skip comments (gfm.h:1563-1571) -- haplotypes are keyed
# by chrom/left/right/alt-list -- so duplicate htN ids across files do not
# collide. Verified in the source before relying on it.
#
# chrX is NOT the same URL shape as the autosomes: it carries a .v2. infix.
# Building the URL by pattern would silently drop it.
set -eux
set -o pipefail
export HOME=/root
export DEBIAN_FRONTEND=noninteractive
BUCKET=s3://rustar-bench
RUN=hisat2-p2b
REGION=us-east-1
MAIN=$$
OUT=/root/p2b
D=/root/data
V=/root/data/v2b
T0=$(date -u +%s)
mkdir -p "$OUT" "$V" /root/bin
: > $OUT/trace.txt
: > $OUT/progress
note() { echo "$(date -u +%H:%M:%S) $*" >> $OUT/progress; }
# The failure path has to publish where the watcher actually looks, and then
# end the box. Previously it killed the heartbeat and synced to out/ only, while
# the watcher read the top-level progress that only the heartbeat wrote -- so a
# 58-second failure was invisible and the instance idled for 2h20m at $2.12/hr.
fail() { rc=$?; note "STAGE_FAILED rc=$rc line=${BASH_LINENO[0]}"
         kill ${HEARTBEAT:-0} 2>/dev/null || true
         aws s3 sync $OUT "$BUCKET/$RUN/out/" --quiet 2>/dev/null || true
         aws s3 cp $OUT/progress   "$BUCKET/$RUN/progress"   --quiet 2>/dev/null || true
         aws s3 cp $OUT/trace.txt  "$BUCKET/$RUN/trace.txt"  --quiet 2>/dev/null || true
         printf '{"ts":%s,"count":0,"busy":0,"script":0,"done":1,"disk_free_gb":0,"phase":"STAGE_FAILED rc=%s"}\n' \
           "$(date -u +%s)" "$rc" > $OUT/status.json
         aws s3 cp $OUT/status.json "$BUCKET/$RUN/status.json" --quiet 2>/dev/null || true
         sleep 5
         shutdown -h now                       # stop, not terminate: the volume survives for debugging
         exit $rc; }
trap fail ERR
elapsed() { echo $(( $(date -u +%s) - T0 )); }
have_big() { [ -f "$1" ] && [ "$(stat -c %s "$1")" -ge "$2" ]; }

note "STAGE_START"
export PATH=/root/bin:/usr/local/bin:$PATH
CHRS="$(seq 1 22) X"

heartbeat() {
  set +e
  while true; do
    ts=$(date -u +%s)
    count=$(($(wc -l < $OUT/trace.txt 2>/dev/null || echo 0)))
    script=$(kill -0 "$MAIN" 2>/dev/null && echo 1 || echo 0)
    done_flag=$(($(grep -c RESULTS_COLLECTED $OUT/progress 2>/dev/null || true)))
    free=$(df -BG --output=avail /root | tail -1 | tr -dc '0-9')
    phase=$(tail -1 $OUT/progress 2>/dev/null | tr -d '"\\' | tr -cd '[:print:]')
    printf '{"ts":%s,"count":%s,"busy":0,"script":%s,"done":%s,"disk_free_gb":%s,"phase":"%s"}\n' \
      "$ts" "$count" "$script" "$done_flag" "${free:-0}" "$phase" > $OUT/status.json
    echo "$ts,$count,${free:-0}" >> $OUT/heartbeat.csv
    for f in trace.txt progress status.json heartbeat.csv; do
      aws s3 cp $OUT/$f "$BUCKET/$RUN/$f" --quiet 2>/dev/null || true
    done
    aws cloudwatch put-metric-data --region "$REGION" --namespace RunKit \
      --metric-data "MetricName=Heartbeat,Dimensions=[{Name=Run,Value=$RUN}],Value=1" \
                    "MetricName=TraceLines,Dimensions=[{Name=Run,Value=$RUN}],Value=$count" \
      2>/dev/null || true
    sleep 60
  done
}
heartbeat & HEARTBEAT=$!
note "HEARTBEAT_UP"

trace() {
  local label=$1; shift
  local tf=$OUT/.t.$$.$RANDOM rc=0
  /usr/bin/time -v -o "$tf" "$@" > $OUT/${label}.log 2>&1 || rc=$?
  local rt=$(awk -F': ' '/Elapsed \(wall clock\)/{print $NF}' "$tf" 2>/dev/null)
  local rss=$(awk '/Maximum resident set size/{print $NF}' "$tf" 2>/dev/null)
  printf '%s\trealtime=%s\tpeak_rss_kb=%s\texit=%s\n' \
    "$label" "${rt:-MISSING}" "${rss:-MISSING}" "$rc" >> $OUT/trace.txt
  note "TRACE $label rt=${rt:-MISSING} rss=${rss:-MISSING} rc=$rc elapsed=$(elapsed)"
  rm -f "$tf"
  return 0
}

# ------------------------------------------------- source, with the 2b fix
cd /root/src && rm -rf hisat2 && mkdir -p hisat2
aws s3 cp "$BUCKET/src/hisat2-src-2b.tar.gz" . --quiet
tar xzf hisat2-src-2b.tar.gz -C hisat2 --strip-components=1
cd hisat2
make -j"$(nproc)" hisat2-align-s hisat2-build-s hisat2-inspect-s >/dev/null 2>&1
cp hisat2 hisat2-build hisat2-inspect hisat2-align-s hisat2-build-s hisat2-inspect-s /root/bin/
chmod +x /root/bin/*
# The whole point of this run is the new backoff; refuse to burn hours without it.
grep -q "selectAlts" hgfm.h
note "HISAT2_READY $(hisat2-build-s --version | head -1)"

# ------------------------------------------------- genome
cd $D
have_big genome.fa 3000000000 || pigz -dc genome.fa.gz > genome.fa
[ -s genome.fa.fai ] || samtools faidx genome.fa
note "GENOME bp=$(awk '{s+=$2} END{print s}' genome.fa.fai) contigs=$(wc -l < genome.fa.fai)"
for c in $CHRS; do
  have_big $V/chr$c.fa 1000000 || samtools faidx genome.fa chr$c > $V/chr$c.fa
done
note "CHR_FASTAS_READY"

# ------------------------------------------------- per-chromosome variant prep
KGBASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
cat > /root/prep_chr.sh <<'PREP'
#!/bin/bash
set -euo pipefail
c=$1; V=/root/data/v2b; OUT=/root/p2b
KGBASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
# Resume: a chromosome that already produced variants is not redone. Each writes
# its own varset file rather than appending to a shared one, so the result is
# idempotent and free of concurrent-append games.
if [ -s $V/chr$c.varset ] && [ -s $V/chr$c.snp ] && [ -s $V/chr$c.haplotype ]; then
  exit 0
fi
# chrX is published with a .v2. infix; the autosomes are not.
if [ "$c" = "X" ]; then
  U=$KGBASE/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
else
  U=$KGBASE/1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
fi
# Eight parallel pulls from the 1000G FTP truncated three files, and "larger
# than 10 MB" happily accepted a half-finished 1.5 GB download -- bcftools then
# died with "Inflate operation failed". bgzip -t validates the whole BGZF
# stream including its EOF block, which is exactly what bcftools needs, so test
# the file rather than its size.
f=$V/chr$c.raw.vcf.gz
for attempt in 1 2 3; do
  if [ -s "$f" ] && bgzip -t "$f" 2>/dev/null; then break; fi
  rm -f "$f"
  curl --retry 5 --retry-delay 10 --connect-timeout 30 -sSfL -o "$f" "$U" || true
done
bgzip -t "$f"                        # still bad after three tries: fail loudly
bcftools view -v snps,indels -q 0.01:minor -Ov "$f" \
 | gawk -F'\t' 'BEGIN{OFS="\t"} /^#/{print; next}
     { if (length($4) > 50) next
       n = split($5, a, ","); for (i = 1; i <= n; i++) if (length(a[i]) > 50) next
       print }' > $V/chr$c.vcf
# --non-rs: the panel names variants 1:10583:G:A, not rsNNN. Without it the
# extractor's only_rs default silently discards every record.
python3 /root/src/hisat2/hisat2_extract_snps_haplotypes_VCF.py --non-rs \
    $V/chr$c.fa $V/chr$c.vcf $V/chr$c > $V/chr$c.extract.log 2>&1
nsnp=$(wc -l < $V/chr$c.snp); nhap=$(wc -l < $V/chr$c.haplotype)
bad=$(grep -c 'seems to be incompatible' $V/chr$c.extract.log || true)
[ "$bad" -eq 0 ]
[ "$nsnp" -gt 0 ]
printf 'chr%s\tvcf=%s\tsnp=%s\thap=%s\tincompatible=%s\n' \
  "$c" "$(grep -vc '^#' $V/chr$c.vcf)" "$nsnp" "$nhap" "$bad" > $V/chr$c.varset
rm -f "$f"                           # 36 GB of panel is not worth keeping
PREP
chmod +x /root/prep_chr.sh
# $CHRS holds newlines (seq), and interpolating it into a bash -c string turned
# them into statement separators: the inner script ran `printf 1`, then `2` as a
# command, and only the last line reached xargs -- with empty input, exiting 0
# having done nothing. Pass the list as a file instead of through the string.
# The failed attempt wrote its rows to one shared file, so the 20 chromosomes
# that succeeded have .snp/.haplotype but no .varset and would be redone (~70
# min). Backfill a varset for anything already complete, so the resume check
# sees them.
for c in $CHRS; do
  if [ -s $V/chr$c.snp ] && [ -s $V/chr$c.haplotype ] && [ ! -s $V/chr$c.varset ]; then
    printf 'chr%s\tvcf=%s\tsnp=%s\thap=%s\tincompatible=0\n' \
      "$c" "$(grep -vc '^#' $V/chr$c.vcf 2>/dev/null || echo 0)" \
      "$(wc -l < $V/chr$c.snp)" "$(wc -l < $V/chr$c.haplotype)" > $V/chr$c.varset
  fi
done
note "BACKFILLED existing=$(ls $V/*.varset 2>/dev/null | wc -l)"
printf '%s\n' $CHRS > /root/chrs.txt
[ "$(wc -l < /root/chrs.txt)" -eq 23 ]
# Canary: one small chromosome first, so a bug in prep_chr.sh costs a minute
# rather than a full fan-out.
trace variant_prep_canary /root/prep_chr.sh 21
[ -s $V/chr21.varset ]
note "CANARY_OK $(cat $V/chr21.varset)"
grep -v '^21$' /root/chrs.txt > /root/chrs_rest.txt
trace variant_prep bash -c \
  "xargs -P 6 -I{} /root/prep_chr.sh {} < /root/chrs_rest.txt"
cat $V/chr*.varset 2>/dev/null | sort -V > $OUT/varsets.txt
note "VARPREP_DONE n=$(wc -l < $OUT/varsets.txt) elapsed=$(elapsed)"
[ "$(wc -l < $OUT/varsets.txt)" -eq 23 ]
sort -V $OUT/varsets.txt >> $OUT/progress

# ------------------------------------------------- concatenate
cd $V
: > $D/genome.snp; : > $D/genome.haplotype
for c in $CHRS; do cat chr$c.snp >> $D/genome.snp; cat chr$c.haplotype >> $D/genome.haplotype; done
NSNP=$(wc -l < $D/genome.snp); NHAP=$(wc -l < $D/genome.haplotype)
note "VARSET_TOTAL snp=$NSNP hap=$NHAP"
[ "$NSNP" -gt 5000000 ]
# Variant ids must be unique genome-wide: the .haplotype rows reference them by
# name. Chromosome-qualified panel ids make that automatic, but check rather
# than trust, because a collision would silently mis-wire haplotypes.
DUP=$(cut -f1 $D/genome.snp | sort | uniq -d | head -5)
[ -z "$DUP" ] || { note "DUPLICATE_SNP_IDS $DUP"; exit 1; }
note "SNP_IDS_UNIQUE"

# ------------------------------------------------- annotation
[ -s $D/genome.gtf ] || pigz -dc $D/gencode.v50.gtf.gz > $D/genome.gtf
trace extract_ss python3 /root/src/hisat2/hisat2_extract_splice_sites.py $D/genome.gtf
cp $OUT/extract_ss.log $D/genome.ss
trace extract_genes python3 /root/src/hisat2/hisat2_extract_genes.py \
  --fai $D/genome.fa.fai -o $D/genome.ht2gm $D/genome.gtf
note "ANNOT ss=$(wc -l < $D/genome.ss) genes=$(grep -c '^G' $D/genome.ht2gm || true)"

# ------------------------------------------------- builds
NP=$(nproc)
build() {
  local label=$1 snpargs=$2
  rm -rf $D/idx_$label; mkdir -p $D/idx_$label
  trace build_$label hisat2-build-s -p $NP $snpargs $D/genome.fa $D/idx_$label/genome
  local sz=$(du -sk $D/idx_$label | cut -f1)
  local ret=$(grep -o 'Local indexes: .*' $OUT/build_$label.log | tail -1)
  printf 'SIZE\t%s\tindex_kb=%s\t%s\n' "$label" "$sz" "${ret:-no-variants}" >> $OUT/trace.txt
  note "SIZE $label index_kb=$sz ${ret:-no-variants}"
  note "DISK_FREE $(df -BG --output=avail /root | tail -1 | tr -dc '0-9')GB"
}

build linear ""
note "P_LINEAR_DONE elapsed=$(elapsed)"
build snp "--snp $D/genome.snp --haplotype $D/genome.haplotype"
note "P_SNP_DONE elapsed=$(elapsed)"

# ------------------------------------------------- gate + collect
rows=$(wc -l < $OUT/trace.txt)
missing=$(grep -c 'MISSING' $OUT/trace.txt || true)
core_ok=1
for k in build_linear build_snp variant_prep; do
  grep -q "^$k	" $OUT/trace.txt || core_ok=0
  gawk -F'exit=' -v k="$k" '$0 ~ "^"k"\t" && $2+0 != 0 {bad=1} END{exit bad+0}' $OUT/trace.txt || core_ok=0
done
# A graph build that carried no variants is the 2a failure mode; refuse to call it a pass.
grep -q 'variant instances retained' $OUT/build_snp.log || core_ok=0
note "GATE rows=$rows missing=$missing core_ok=$core_ok"

cp $D/genome.ss $D/genome.ht2gm $OUT/ 2>/dev/null || true
pigz -f $OUT/genome.ss $OUT/genome.ht2gm 2>/dev/null || true
pigz -kf $D/genome.snp $D/genome.haplotype
aws s3 cp $D/genome.snp.gz "$BUCKET/$RUN/refs/" --quiet
aws s3 cp $D/genome.haplotype.gz "$BUCKET/$RUN/refs/" --quiet
# The indexes are the deliverable: keep them off the box's fate.
aws s3 sync $D/idx_snp    "$BUCKET/$RUN/refs/grch38_v50_1kgp_snp/" --quiet
aws s3 sync $D/idx_linear "$BUCKET/$RUN/refs/grch38_v50/"          --quiet
note "INDEXES_UPLOADED"

kill $HEARTBEAT 2>/dev/null || true
aws s3 sync $OUT "$BUCKET/$RUN/out/" --quiet
aws s3 cp $OUT/trace.txt "$BUCKET/$RUN/trace.txt"
local_size=$(stat -c %s $OUT/trace.txt)
remote_size=$(aws s3api head-object --bucket "${BUCKET#s3://}" --key "$RUN/trace.txt" \
                --query ContentLength --output text 2>/dev/null || echo -1)
[ "$local_size" = "$remote_size" ] || { note "S3_UNVERIFIED local=$local_size remote=$remote_size"; exit 1; }
note "S3_VERIFIED size=$local_size"
if [ "$core_ok" -eq 1 ] && [ "$missing" -eq 0 ]; then
  note "BUILD_PASS"; note "RESULTS_COLLECTED"
else
  note "BUILD_FAIL missing=$missing core_ok=$core_ok"
fi
printf '{"ts":%s,"count":%s,"busy":0,"script":0,"done":1,"disk_free_gb":0,"phase":"%s"}\n' \
  "$(date -u +%s)" "$rows" "$(tail -1 $OUT/progress | tr -d '"\\')" > $OUT/status.json
aws s3 cp $OUT/status.json "$BUCKET/$RUN/status.json" --quiet
aws s3 cp $OUT/progress    "$BUCKET/$RUN/progress"    --quiet
sleep 10
shutdown -h now
