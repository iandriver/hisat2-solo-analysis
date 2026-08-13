#!/bin/bash
# Phase 2a, second pass. Idempotent: reuses anything pass 1 already put on disk
# (genome.fa.gz, chr1.fa, chr1.gtf, kg_chr1.vcf.gz) rather than re-downloading.
#
# Four bugs from pass 1, all fixed here:
#
#   1. Ubuntu's tabix is built WITHOUT libcurl, so the remote dbSNP range query
#      died with "Protocol not supported" -- the whole dbSNP arm produced zero
#      variants. htslib is now compiled with --enable-libcurl, and a 2 kb probe
#      query proves remote access works before the 2 GB slice is attempted.
#   2. `tabix ... | bgzip > out` reported rc=0 because the pipeline's status is
#      bgzip's. Every pipeline now runs under `set -o pipefail`.
#   3. The 1000G 30x panel uses POSITIONAL ids (1:10583:G:A), not rs ids, so the
#      extractor's only_rs default silently discarded all 1.27M variants. The
#      1000G extraction now passes --non-rs. dbSNP keeps the rs-only default,
#      where it is meaningful.
#   4. hisat2-build accepts an empty --snp file, builds a linear index and exits
#      0. So pass 1 "succeeded" at measuring nothing, and the gate would have
#      called it PROBE_PASS. Variant counts are now asserted before each build,
#      recorded in the SIZE row, and required by the gate.
#
# Also corrected: GENCODE's primary assembly is NOT soft-masked (pass 1 measured
# softmasked=0). The uppercase step was a no-op and is kept only as cheap
# insurance for a differently-sourced FASTA; the soft-mask concern was real for
# UCSC's per-chromosome files, not for GENCODE.
set -eux
set -o pipefail
export HOME=/root
export DEBIAN_FRONTEND=noninteractive
BUCKET=s3://rustar-bench
RUN=hisat2-p2a
REGION=us-east-1
MAIN=$$
OUT=/root/p2a2
D=/root/data
BUDGET_SEC=${BUDGET_SEC:-6000}
T0=$(date -u +%s)
mkdir -p "$OUT" /root/bin "$D"
: > $OUT/trace.txt
note() { echo "$(date -u +%H:%M:%S) $*" >> $OUT/progress; }
fail() { rc=$?; note "STAGE_FAILED rc=$rc line=${BASH_LINENO[0]}"
         kill ${HEARTBEAT:-0} 2>/dev/null || true
         aws s3 sync $OUT "$BUCKET/$RUN/out2/" --quiet 2>/dev/null || true
         exit $rc; }
trap fail ERR
elapsed() { echo $(( $(date -u +%s) - T0 )); }
budget_left() { echo $(( BUDGET_SEC - $(elapsed) )); }
have() { [ -s "$1" ]; }          # reuse guard: present AND non-empty
# `-s` is too weak to reuse across a failed pass: a dead `tabix | bgzip` leaves a
# 28-byte EOF block, which is non-empty and would silently skip a 2 GB download.
# Reuse of anything expensive must clear a plausible size floor.
have_big() { [ -f "$1" ] && [ "$(stat -c %s "$1")" -ge "$2" ]; }

note "STAGE_START pass=2"
export PATH=/root/bin:/usr/local/bin:$PATH

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
      aws s3 cp $OUT/$f "$BUCKET/$RUN/pass2/$f" --quiet 2>/dev/null || true
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
  local tf=$OUT/.t.$$ rc=0
  /usr/bin/time -v -o "$tf" "$@" > $OUT/${label}.log 2>&1 || rc=$?
  local rt=$(awk -F': ' '/Elapsed \(wall clock\)/{print $NF}' "$tf" 2>/dev/null)
  local rss=$(awk '/Maximum resident set size/{print $NF}' "$tf" 2>/dev/null)
  printf '%s\trealtime=%s\tpeak_rss_kb=%s\texit=%s\n' \
    "$label" "${rt:-MISSING}" "${rss:-MISSING}" "$rc" >> $OUT/trace.txt
  note "TRACE $label rt=${rt:-MISSING} rss=${rss:-MISSING} rc=$rc elapsed=$(elapsed)"
  rm -f "$tf"
  return 0
}

# ------------------------------------------------- fix 1: htslib with libcurl
apt-get update -qq
apt-get install -y -qq libcurl4-openssl-dev libssl-dev liblzma-dev libbz2-dev \
                       libdeflate-dev autoconf automake gawk >/dev/null
if ! /usr/local/bin/tabix --version 2>/dev/null | grep -q htslib; then
  cd /root && curl -sSL -o htslib.tar.bz2 \
    https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
  tar xjf htslib.tar.bz2 && cd htslib-1.21
  ./configure --enable-libcurl --prefix=/usr/local >/dev/null 2>&1
  make -j"$(nproc)" >/dev/null 2>&1 && make install >/dev/null 2>&1
  ldconfig
fi
hash -r
/usr/local/bin/tabix --version | head -1 >> $OUT/progress
# Fail fast and cheap: a 2 kb probe, not a 2 GB slice.
DBSNP=https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
PROBE_N=$(/usr/local/bin/tabix "$DBSNP" NC_000001.11:1000000-1002000 2>&1 | grep -c '^NC_' || true)
note "TABIX_REMOTE_PROBE rows=$PROBE_N"
[ "$PROBE_N" -gt 100 ]
note "HTSLIB_READY"

# ------------------------------------------------- inputs (reuse pass 1)
cd $D
have chr1.fa   || { echo "chr1.fa missing from pass 1"; exit 1; }
have chr1.gtf  || { echo "chr1.gtf missing from pass 1"; exit 1; }
[ -s chr1.fa.fai ] || samtools faidx chr1.fa
note "INPUTS_REUSED chr1_bp=$(cut -f2 chr1.fa.fai)"

# --- dbSNP b157 chr1 slice --------------------------------------------------
if ! have_big dbsnp157_chr1.raw.vcf.gz 100000000; then
  trace dbsnp_slice bash -o pipefail -c \
    "/usr/local/bin/tabix -h '$DBSNP' NC_000001.11 | /usr/local/bin/bgzip -@ 8 -c > $D/dbsnp157_chr1.raw.vcf.gz"
fi
RAW_N=$(/usr/local/bin/bgzip -dc dbsnp157_chr1.raw.vcf.gz | grep -vc '^#' || true)
note "DBSNP_RAW n=$RAW_N"
[ "$RAW_N" -gt 1000000 ]          # chr1 carries 93.5M RS; anything less is a broken slice

cat > /root/common.awk <<'AWK'
BEGIN{ FS="\t"; OFS="\t"
       split("1000Genomes_30X GnomAD_genomes GnomAD_exomes TOPMED", B, " ")
       for (i in B) BIG[B[i]] = 1 }
/^#/ { print $0 > OUT_BIG; print $0 > OUT_KG; next }
{
  if ($1 != "NC_000001.11") next
  $1 = "chr1"
  if (length($4) > 50) next
  if ($5 ~ /[<\[\]]/) next
  n = split($5, alts, ","); bad = 0
  for (a = 1; a <= n; a++) if (length(alts[a]) > 50) bad = 1
  if (bad) next
  if (!match($0, /FREQ=[^;\t]*/)) next
  freq = substr($0, RSTART+5, RLENGTH-5)
  ns = split(freq, srcs, "|"); okBig = 0; okAny = 0; okKg = 0
  for (s = 1; s <= ns; s++) {
    p = index(srcs[s], ":"); if (!p) continue
    src = substr(srcs[s], 1, p-1)
    nf = split(substr(srcs[s], p+1), f, ",")
    for (j = 2; j <= nf; j++) {
      if (f[j] == ".") continue
      if (f[j]+0 >= 0.01 && f[j]+0 <= 0.99) {
        okAny = 1
        if (src in BIG) okBig = 1
        if (src == "1000Genomes_30X") okKg = 1
      }
    }
  }
  if (okAny) nAny++
  if (okKg)  { nKg++;  print $0 > OUT_KG }
  if (okBig) { nBig++; print $0 > OUT_BIG }
}
END{ printf "COMMON_BIG4=%d COMMON_1KG30X=%d COMMON_ANY=%d\n", nBig, nKg, nAny > "/root/common.counts" }
AWK
trace dbsnp_filter bash -o pipefail -c \
  "/usr/local/bin/bgzip -dc $D/dbsnp157_chr1.raw.vcf.gz | \
   gawk -v OUT_BIG=$D/dbsnp_big4.vcf -v OUT_KG=$D/dbsnp_kg30x.vcf -f /root/common.awk"
note "DBSNP_COMMON $(cat /root/common.counts)"
grep -vc '^#' $D/dbsnp_big4.vcf | grep -qvx 0

# --- 1000G 30x phased chr1 (already filtered in pass 1) --------------------
if ! have_big kg_phased.vcf 10000000; then
  KG=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
  have_big kg_chr1.vcf.gz 100000000 || trace kg_fetch curl -sSL -o $D/kg_chr1.vcf.gz "$KG"
  trace kg_filter bash -o pipefail -c \
    "bcftools view -v snps,indels -q 0.01:minor -Ov $D/kg_chr1.vcf.gz | \
     gawk -F'\t' 'BEGIN{OFS=\"\t\"} /^#/{print; next}
       { if (length(\$4) > 50) next
         n = split(\$5, a, \",\"); for (i = 1; i <= n; i++) if (length(a[i]) > 50) next
         print }' > $D/kg_phased.vcf"
fi
note "KG_PHASED n=$(grep -vc '^#' kg_phased.vcf) samples=$(grep -m1 '^#CHROM' kg_phased.vcf | gawk '{print NF-9}')"

# ------------------------------------------------- fix 3+4: extract, then ASSERT
extract() {
  local label=$1 vcf=$2; shift 2
  trace extract_$label timeout 5400 python3 \
    /root/src/hisat2/hisat2_extract_snps_haplotypes_VCF.py "$@" $D/chr1.fa "$vcf" $D/$label
  local incompat=$(grep -c 'seems to be incompatible' $OUT/extract_$label.log || true)
  local nsnp=$(wc -l < $D/$label.snp 2>/dev/null || echo 0)
  local nhap=$(wc -l < $D/$label.haplotype 2>/dev/null || echo 0)
  local nvcf=$(grep -vc '^#' "$vcf" || true)
  printf 'VARSET\t%s\tvcf=%s\tsnp=%s\thap=%s\tincompatible=%s\n' \
    "$label" "$nvcf" "$nsnp" "$nhap" "$incompat" >> $OUT/trace.txt
  note "VARSET $label vcf=$nvcf snp=$nsnp hap=$nhap incompatible=$incompat"
  [ "$incompat" -eq 0 ]
  # An extraction that yields no variants is the pass-1 failure mode. Never
  # let it reach a build, where it would masquerade as a successful linear one.
  [ "$nsnp" -gt 0 ]
}
extract dbsnp_big4  $D/dbsnp_big4.vcf
extract dbsnp_kg30x $D/dbsnp_kg30x.vcf
# --non-rs: the 1000G panel names variants 1:10583:G:A, not rsNNN.
extract kg_phased   $D/kg_phased.vcf --non-rs || note "KG_EXTRACT_FAILED"

[ -s $D/chr1.ht2gm ] || trace extract_genes python3 /root/src/hisat2/hisat2_extract_genes.py \
  --fai $D/chr1.fa.fai -o $D/chr1.ht2gm $D/chr1.gtf
note "EXTRACT_DONE elapsed=$(elapsed)"

# ------------------------------------------------- builds
NP=$(nproc)
build() {
  local label=$1 snp=$2; shift 2
  local nsnp=0 nhap=0 args=""
  if [ -n "$snp" ]; then
    nsnp=$(wc -l < $D/$snp.snp); nhap=$(wc -l < $D/$snp.haplotype)
    [ "$nsnp" -gt 0 ]                       # fix 4: never build a fake graph index
    args="--snp $D/$snp.snp --haplotype $D/$snp.haplotype"
  fi
  rm -rf idx_$label; mkdir -p idx_$label
  trace build_$label hisat2-build-s -p $NP $args "$@" $D/chr1.fa idx_$label/chr1
  local sz=$(du -sk idx_$label | cut -f1)
  local ret=$(grep -o 'Local indexes: .*' $OUT/build_$label.log | tail -1)
  printf 'SIZE\t%s\tindex_kb=%s\tthreads=%s\tsnp=%s\thap=%s\t%s\n' \
    "$label" "$sz" "$NP" "$nsnp" "$nhap" "${ret:-no-variants}" >> $OUT/trace.txt
  note "SIZE $label index_kb=$sz threads=$NP snp=$nsnp hap=$nhap ${ret:-no-variants}"
  rm -rf idx_$label
}

build linear ""
note "P4_LINEAR elapsed=$(elapsed) left=$(budget_left)"
build dbsnp_big4 dbsnp_big4
note "P4_BIG4 elapsed=$(elapsed) left=$(budget_left)"

U=$(gawk -F'realtime=' '/^build_dbsnp_big4\t/{split($2,a,"\t"); n=split(a[1],t,":");
  if(n==3) print int(t[1]*3600+t[2]*60+t[3]); else print int(t[1]*60+t[2])}' $OUT/trace.txt | tail -1)
U=${U:-300}
if [ "$U" -lt 60 ]; then U=60; fi   # a bare `[ ] && x` at top level trips set -e
note "UNIT big4_sec=$U left=$(budget_left)"

# density arm, then the real-haplotype arm, then the knob titration
if [ "$(budget_left)" -gt $(( U * 2 )) ] && have $D/dbsnp_kg30x.snp; then
  build dbsnp_kg30x dbsnp_kg30x
else note "SKIP dbsnp_kg30x left=$(budget_left)"; fi

if [ "$(budget_left)" -gt $(( U * 3 )) ] && have $D/kg_phased.snp; then
  build kg_phased kg_phased
else note "SKIP kg_phased left=$(budget_left)"; fi

for cfg in "p8:" "offrate5:--offrate 5" "ftab12:--ftabchars 12" "bmaxdivn8:--bmaxdivn 8"; do
  lab=${cfg%%:*}; args=${cfg#*:}
  if [ "$(budget_left)" -lt $(( U * 2 )) ]; then note "SKIP tit_$lab left=$(budget_left)"; continue; fi
  if [ "$lab" = "p8" ]; then NP=8; else NP=$(nproc); fi
  build tit_$lab dbsnp_big4 $args
  NP=$(nproc)
done
note "WORK_DONE elapsed=$(elapsed)"

# ------------------------------------------------- gate
rows=$(wc -l < $OUT/trace.txt)
missing=$(grep -c 'MISSING' $OUT/trace.txt || true)
core_ok=1
for k in extract_dbsnp_big4 build_linear build_dbsnp_big4; do
  grep -q "^$k	" $OUT/trace.txt || core_ok=0
  gawk -F'exit=' -v k="$k" '$0 ~ "^"k"\t" && $2+0 != 0 {bad=1} END{exit bad+0}' $OUT/trace.txt || core_ok=0
done
# fix 4, at the gate: the graph build must actually have carried variants.
gawk -F'\t' '$1=="SIZE" && $2=="dbsnp_big4" {for(i=1;i<=NF;i++) if($i ~ /^snp=/){split($i,a,"="); if(a[2]+0>0) ok=1}} END{exit !ok}' \
  $OUT/trace.txt || core_ok=0
failed=$(gawk -F'exit=' '/^(build_|extract_|dbsnp_|kg_)/ && $2 != "" && $2+0 != 0 {n++} END{print n+0}' $OUT/trace.txt)
note "GATE rows=$rows missing=$missing failed=$failed core_ok=$core_ok"

cp /root/common.counts $OUT/ 2>/dev/null || true
cp $D/*.snp $D/*.haplotype $D/chr1.ht2gm $OUT/ 2>/dev/null || true
pigz -f $OUT/*.snp $OUT/*.haplotype $OUT/chr1.ht2gm 2>/dev/null || true
kill $HEARTBEAT 2>/dev/null || true
aws s3 sync $OUT "$BUCKET/$RUN/out2/" --quiet
aws s3 cp $OUT/trace.txt "$BUCKET/$RUN/pass2/trace.txt"
local_size=$(stat -c %s $OUT/trace.txt)
remote_size=$(aws s3api head-object --bucket "${BUCKET#s3://}" --key "$RUN/pass2/trace.txt" \
                --query ContentLength --output text 2>/dev/null || echo -1)
[ "$local_size" = "$remote_size" ] || { note "S3_UNVERIFIED local=$local_size remote=$remote_size"; exit 1; }
note "S3_VERIFIED size=$local_size"
if [ "$core_ok" -eq 1 ] && [ "$missing" -eq 0 ]; then
  note "PROBE_PASS"; note "RESULTS_COLLECTED"
else
  note "PROBE_FAIL missing=$missing core_ok=$core_ok"
fi
printf '{"ts":%s,"count":%s,"busy":0,"script":0,"done":1,"disk_free_gb":0,"phase":"%s"}\n' \
  "$(date -u +%s)" "$rows" "$(tail -1 $OUT/progress | tr -d '"\\')" > $OUT/status.json
aws s3 cp $OUT/status.json "$BUCKET/$RUN/pass2/status.json" --quiet
aws s3 cp $OUT/progress    "$BUCKET/$RUN/pass2/progress"    --quiet
sleep 10
shutdown -h now
