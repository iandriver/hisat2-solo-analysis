#!/bin/bash
# Phase 2a -- chr1 probe on a current reference (GENCODE v50 + dbSNP b157).
#
# Purpose, in priority order:
#   1. Price the whole-genome SNP build before committing: peak RSS, wall time,
#      index size and variant retention on ~8% of the genome. chr1 extrapolates
#      well because the local-index cost is driven by variant DENSITY per 56 kb
#      window, not by chromosome length.
#   2. Measure retention against variant density. dbSNP b157 has no COMMON flag
#      (dropped after b151), so "common" must be derived from FREQ=. The choice
#      of frequency source swings density 5x, so build two densities:
#        big4  = MAF>=1% in 1000Genomes_30X | GnomAD_genomes | GnomAD_exomes |
#                TOPMED  (~2x the shipped grch38_snp density)
#        kg30x = MAF>=1% in 1000Genomes_30X only  (~1.5x shipped)
#      Rejected: "any source", where 36% of survivors rest solely on SGDP_PRJ
#      (~280 samples, one het reads as 50%) -- ~10x shipped density, unbuildable.
#   3. Compare invented haplotypes against real ones. A sample-less VCF sends
#      hisat2_extract_snps_haplotypes_VCF.py into its greedy graph-colouring
#      branch (:298-330), so dbSNP haplotypes are fiction. 1000G 30x phased
#      carries 3202 real phased samples on the same population basis.
#   4. Titrate hisat2-build's size/time knobs, but only with leftover budget --
#      the baseline must land even if the titration does not.
#
# Verified locally before this ran (runkit rule 1: never debug on the expensive
# box): the FREQ= filter, the extractor on a real dbSNP slice, and the soft-mask
# behaviour below.
set -eux
export HOME=/root
export DEBIAN_FRONTEND=noninteractive
BUCKET=s3://rustar-bench
RUN=hisat2-p2a
REGION=us-east-1
MAIN=$$
OUT=/root/p2a
BUDGET_SEC=${BUDGET_SEC:-6600}          # hard cap on the build section
T0=$(date -u +%s)
mkdir -p "$OUT" /root/bin
: > $OUT/trace.txt
note() { echo "$(date -u +%H:%M:%S) $*" >> $OUT/progress; }
fail() { rc=$?; note "STAGE_FAILED rc=$rc line=${BASH_LINENO[0]}"
         kill ${HEARTBEAT:-0} 2>/dev/null || true
         aws s3 sync $OUT "$BUCKET/$RUN/out/" --quiet 2>/dev/null || true
         exit $rc; }
trap fail ERR
elapsed() { echo $(( $(date -u +%s) - T0 )); }
budget_left() { echo $(( BUDGET_SEC - $(elapsed) )); }

# ---------------------------------------------------------------- phase 0
note "STAGE_START"
apt-get update -qq
apt-get install -y -qq unzip curl build-essential zlib1g-dev python3 \
                       tabix bcftools pigz time samtools gawk >/dev/null
curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscli.zip
unzip -q /tmp/awscli.zip -d /tmp && /tmp/aws/install >/dev/null
aws s3 ls "$BUCKET/" >/dev/null       # prove credentials, not merely presence
note "AWSCLI_READY"

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

# trace <label> <cmd...>  -- always returns 0; the gate reads exit= from the row.
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

# ---------------------------------------------------------------- phase 1
cd /root && mkdir -p src && cd src
aws s3 cp "$BUCKET/src/hisat2-src.tar.gz" . --quiet
mkdir -p hisat2 && tar xzf hisat2-src.tar.gz -C hisat2 --strip-components=1
cd hisat2
make -j$(nproc) hisat2-align-s hisat2-build-s hisat2-inspect-s >/dev/null 2>&1
cp hisat2 hisat2-build hisat2-inspect hisat2-align-s hisat2-build-s hisat2-inspect-s /root/bin/
chmod +x /root/bin/*
export PATH=/root/bin:$PATH
hisat2-build-s --version | head -1 >> $OUT/progress
note "HISAT2_READY"

# ---------------------------------------------------------------- phase 2
D=/root/data; mkdir -p $D; cd $D
G=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50
trace fetch_gencode bash -c "
  curl -sSL -o $D/genome.fa.gz '$G/GRCh38.primary_assembly.genome.fa.gz' &&
  curl -sSL -o $D/gencode.v50.gtf.gz '$G/gencode.v50.primary_assembly.annotation.gtf.gz'"

# chr1 only, uppercased. GENCODE ships soft-masked: on UCSC chr1, 119.2M of
# 249M bases are lowercase, and the extractor's reference check is
# case-sensitive (:103), so it fires ~823 errors per 200 kb of variants.
# Verified locally that it prints without `continue`, so output is byte-
# identical either way -- but clearing the noise is what lets a REMAINING
# incompatibility be treated as a real assembly-mismatch failure, below.
pigz -dc genome.fa.gz | awk '/^>/{k=($1==">chr1")} k' > chr1.raw.fa
grep -c '^>' chr1.raw.fa | grep -qx 1
head -1 chr1.raw.fa | grep -q '^>chr1'
LOWER=$(grep -v '^>' chr1.raw.fa | tr -cd 'acgt' | wc -c)
CHR1_BP=$(grep -v '^>' chr1.raw.fa | tr -d '\n' | wc -c)
awk '/^>/{print ">chr1"; next}{print toupper($0)}' chr1.raw.fa > chr1.fa
rm -f chr1.raw.fa
samtools faidx chr1.fa
note "CHR1_FA bp=$CHR1_BP softmasked=$LOWER"

pigz -dc gencode.v50.gtf.gz | awk -F'\t' '$1=="chr1" || /^#/' > chr1.gtf
note "CHR1_GTF lines=$(grep -vc '^#' chr1.gtf)"

# --- dbSNP b157 chr1 slice: remote tabix range, not a 25 GB download --------
DBSNP=https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
trace dbsnp_slice bash -c \
  "tabix -h '$DBSNP' NC_000001.11 | bgzip -@ 8 -c > $D/dbsnp157_chr1.raw.vcf.gz"
RAW_N=$(bgzip -dc dbsnp157_chr1.raw.vcf.gz | grep -vc '^#' || true)
note "DBSNP_RAW n=$RAW_N"

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
trace dbsnp_filter bash -c \
  "bgzip -dc $D/dbsnp157_chr1.raw.vcf.gz | \
   gawk -v OUT_BIG=$D/dbsnp_big4.vcf -v OUT_KG=$D/dbsnp_kg30x.vcf -f /root/common.awk"
note "DBSNP_COMMON $(cat /root/common.counts)"

# --- 1000G 30x phased chr1: 3202 real phased samples -----------------------
KG=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
trace kg_fetch curl -sSL -o $D/kg_chr1.vcf.gz "$KG"
# Length filtering is done in awk rather than a bcftools expression: -v snps,indels
# already drops symbolic SV alleles, and an untested filter expression is not worth
# the risk of discovering its syntax on a $2/hr box.
trace kg_filter bash -c \
  "bcftools view -v snps,indels -q 0.01:minor -Ov $D/kg_chr1.vcf.gz | \
   gawk -F'\t' 'BEGIN{OFS=\"\t\"} /^#/{print; next}
     { if (length(\$4) > 50) next
       n = split(\$5, a, \",\"); for (i = 1; i <= n; i++) if (length(a[i]) > 50) next
       print }' > $D/kg_phased.vcf"
KG_N=$(grep -vc '^#' $D/kg_phased.vcf || true)
KG_SAMPLES=$(grep -m1 '^#CHROM' $D/kg_phased.vcf | awk '{print NF-9}')
note "KG_PHASED n=$KG_N samples=$KG_SAMPLES"
[ "${KG_SAMPLES:-0}" -gt 100 ]   # if this is 0 the whole phasing comparison is vacuous

# ---------------------------------------------------------------- phase 3
# extract <label> <vcf>  -- and treat any surviving reference incompatibility as
# fatal. After uppercasing there is no legitimate source of these, so a nonzero
# count means the VCF and the FASTA are different assemblies.
extract() {
  local label=$1 vcf=$2 lim=${3:-3600}
  trace extract_$label timeout "$lim" python3 \
    /root/src/hisat2/hisat2_extract_snps_haplotypes_VCF.py $D/chr1.fa "$vcf" $D/$label
  local incompat=$(grep -c 'seems to be incompatible' $OUT/extract_$label.log || true)
  local nsnp=$(wc -l < $D/$label.snp 2>/dev/null || echo 0)
  local nhap=$(wc -l < $D/$label.haplotype 2>/dev/null || echo 0)
  printf 'VARSET\t%s\tsnp=%s\thap=%s\tincompatible=%s\n' \
    "$label" "$nsnp" "$nhap" "$incompat" >> $OUT/trace.txt
  note "VARSET $label snp=$nsnp hap=$nhap incompatible=$incompat"
  [ "$incompat" -eq 0 ]
}
extract dbsnp_big4  $D/dbsnp_big4.vcf
extract dbsnp_kg30x $D/dbsnp_kg30x.vcf
# 1000G is the slow one: 3202 samples means the haplotype builder does
# O(6404 x variants) string work in pure Python. Capped so it cannot eat the
# budget; a timeout here is a finding, not a failure.
extract kg_phased   $D/kg_phased.vcf 5400 || note "KG_EXTRACT_INCOMPLETE"

trace extract_ss python3 /root/src/hisat2/hisat2_extract_splice_sites.py $D/chr1.gtf
trace extract_genes python3 /root/src/hisat2/hisat2_extract_genes.py \
  --fai $D/chr1.fa.fai -o $D/chr1.ht2gm $D/chr1.gtf
note "EXTRACT_DONE elapsed=$(elapsed)"

# ---------------------------------------------------------------- phase 4
NP=$(nproc)
build() {
  local label=$1; shift
  rm -rf idx_$label; mkdir -p idx_$label
  trace build_$label hisat2-build-s -p $NP "$@" $D/chr1.fa idx_$label/chr1
  local sz=$(du -sk idx_$label | cut -f1)
  local ret=$(grep -o 'Local indexes: .*' $OUT/build_$label.log | tail -1)
  printf 'SIZE\t%s\tindex_kb=%s\tthreads=%s\t%s\n' \
    "$label" "$sz" "$NP" "${ret:-no-variants}" >> $OUT/trace.txt
  note "SIZE $label index_kb=$sz threads=$NP ${ret:-no-variants}"
  rm -rf idx_$label          # sizes are recorded; the bytes are not the deliverable
}

build linear
note "P4_LINEAR elapsed=$(elapsed) left=$(budget_left)"
build dbsnp_big4 --snp $D/dbsnp_big4.snp --haplotype $D/dbsnp_big4.haplotype
note "P4_BIG4 elapsed=$(elapsed) left=$(budget_left)"

# The big4 build's own wall time is the unit of account for everything optional.
unit_sec() {
  awk -F'realtime=' -v k="$1" '$0 ~ "^"k"\t" {split($2,a,"\t"); n=split(a[1],t,":");
    if(n==3) print int(t[1]*3600+t[2]*60+t[3]); else print int(t[1]*60+t[2])}' $OUT/trace.txt | tail -1
}
U=$(unit_sec build_dbsnp_big4); U=${U:-900}
note "UNIT big4_sec=$U left=$(budget_left)"

# The density arm. Higher value than the knob titration, so it goes first.
if [ "$(budget_left)" -gt $(( U * 2 )) ] && [ -s $D/dbsnp_kg30x.snp ]; then
  build dbsnp_kg30x --snp $D/dbsnp_kg30x.snp --haplotype $D/dbsnp_kg30x.haplotype
else
  note "SKIP dbsnp_kg30x left=$(budget_left)"
fi

# The real-haplotype arm. Genuinely may explode; it runs where a failure costs
# only itself.
if [ "$(budget_left)" -gt $(( U * 2 )) ] && [ -s $D/kg_phased.snp ]; then
  build kg_phased --snp $D/kg_phased.snp --haplotype $D/kg_phased.haplotype
else
  note "SKIP kg_phased left=$(budget_left)"
fi

# --- build-parameter titration, strictly on leftovers ----------------------
for cfg in "p8:-p 8" "offrate5:--offrate 5" "ftab12:--ftabchars 12" "bmaxdivn8:--bmaxdivn 8"; do
  lab=${cfg%%:*}; args=${cfg#*:}
  if [ "$(budget_left)" -lt $(( U * 2 )) ]; then
    note "SKIP tit_$lab left=$(budget_left)"; continue
  fi
  if [ "$lab" = "p8" ]; then NP=8; args=""; else NP=$(nproc); fi
  build tit_$lab --snp $D/dbsnp_big4.snp --haplotype $D/dbsnp_big4.haplotype $args
  NP=$(nproc)
done
note "WORK_DONE elapsed=$(elapsed)"

# ---------------------------------------------------------------- phase 5
# The probe's job is the big4 measurement plus a linear baseline. Everything
# else is upside, so the gate asks only for those.
rows=$(wc -l < $OUT/trace.txt)
missing=$(grep -c 'MISSING' $OUT/trace.txt || true)
core_ok=1
for k in extract_dbsnp_big4 build_linear build_dbsnp_big4; do
  grep -q "^$k	" $OUT/trace.txt || core_ok=0
  awk -F'exit=' -v k="$k" '$0 ~ "^"k"\t" && $2+0 != 0 {bad=1} END{exit bad+0}' $OUT/trace.txt || core_ok=0
done
failed=$(awk -F'exit=' '/^(build_|extract_|dbsnp_|kg_|fetch_)/ && $2 != "" && $2+0 != 0 {n++} END{print n+0}' $OUT/trace.txt)
note "GATE rows=$rows missing=$missing failed=$failed core_ok=$core_ok"

cp /root/common.counts $OUT/ 2>/dev/null || true
cp $D/*.snp $D/*.haplotype $D/chr1.ht2gm $OUT/ 2>/dev/null || true
pigz -f $OUT/*.snp $OUT/*.haplotype $OUT/chr1.ht2gm 2>/dev/null || true
kill $HEARTBEAT 2>/dev/null || true
aws s3 sync $OUT "$BUCKET/$RUN/out/" --quiet
aws s3 cp $OUT/trace.txt "$BUCKET/$RUN/trace.txt"
local_size=$(stat -c %s $OUT/trace.txt)
remote_size=$(aws s3api head-object --bucket "${BUCKET#s3://}" --key "$RUN/trace.txt" \
                --query ContentLength --output text 2>/dev/null || echo -1)
if [ "$local_size" != "$remote_size" ]; then
  note "S3_UNVERIFIED local=$local_size remote=$remote_size -- staying up"; exit 1
fi
note "S3_VERIFIED size=$local_size"
if [ "$core_ok" -eq 1 ] && [ "$missing" -eq 0 ]; then
  note "PROBE_PASS"; note "RESULTS_COLLECTED"
else
  note "PROBE_FAIL missing=$missing core_ok=$core_ok"
fi
printf '{"ts":%s,"count":%s,"busy":0,"script":0,"done":1,"disk_free_gb":0,"phase":"%s"}\n' \
  "$(date -u +%s)" "$rows" "$(tail -1 $OUT/progress | tr -d '"\\')" > $OUT/status.json
aws s3 cp $OUT/status.json "$BUCKET/$RUN/status.json" --quiet
aws s3 cp $OUT/progress    "$BUCKET/$RUN/progress"    --quiet
sleep 10
shutdown -h now
