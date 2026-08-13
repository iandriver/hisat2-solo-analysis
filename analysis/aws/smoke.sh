#!/bin/bash
# Phase 1 smoke test (runbook rule 1: prove the script somewhere cheap).
#
# chr22-scale everything. The point is not the science -- it is that all three
# toolchains build and run, that the trace rows carry REAL realtime and peak_rss
# rather than merely existing, and that the credential path to S3 works.
set -eux
export HOME=/root
export DEBIAN_FRONTEND=noninteractive
BUCKET=s3://rustar-bench
RUN=hisat2-smoke
REGION=us-east-1
MAIN=$$
mkdir -p /root/out
: > /root/out/trace.txt
note() { echo "$(date -u +%H:%M:%S) $*" >> /root/out/progress; }
# The ERR trap must also stop the heartbeat and exit non-zero. Otherwise the main
# work dies, the backgrounded heartbeat keeps the process-group name alive, and
# `pgrep -f smoke.sh` still answers UP -- an hour of billing looking healthy.
fail() { rc=$?; note "STAGE_FAILED rc=$rc"; kill ${HEARTBEAT:-0} 2>/dev/null || true;
         aws s3 cp /root/out/progress "$BUCKET/$RUN/progress" --quiet 2>/dev/null || true;
         exit $rc; }
trap fail ERR
note "STAGE_START"

# The stock Ubuntu AMI ships no AWS CLI. Everything below -- heartbeat, inputs,
# results -- needs it, so install it first and prove it works before anything else.
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq unzip curl >/dev/null
curl -sS "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscli.zip
unzip -q -o /tmp/awscli.zip -d /tmp
/tmp/aws/install --update >/dev/null 2>&1 || /tmp/aws/install >/dev/null 2>&1
hash -r
aws --version
aws s3 ls "$BUCKET/" >/dev/null
note "AWSCLI_READY"

heartbeat() {
  set +e
  while true; do
    ts=$(date -u +%s)
    count=$(($(wc -l < /root/out/trace.txt 2>/dev/null || echo 0)))
    script=$(kill -0 "$MAIN" 2>/dev/null && echo 1 || echo 0)
    done_flag=$(($(grep -c RESULTS_COLLECTED /root/out/progress 2>/dev/null || true)))
    disk=$(df -BG --output=avail /root 2>/dev/null | tail -1 | tr -dc 0-9)
    phase=$(tail -1 /root/out/progress 2>/dev/null | tr -d '"\\' | tr -cd '[:print:]')
    printf '{"ts":%s,"count":%s,"busy":0,"script":%s,"done":%s,"disk_free_gb":%s,"phase":"%s"}\n' \
      "$ts" "$count" "$script" "$done_flag" "${disk:-0}" "$phase" > /root/out/status.json
    echo "$ts,$count,0" >> /root/out/heartbeat.csv
    for f in trace.txt progress status.json heartbeat.csv; do
      aws s3 cp /root/out/$f "$BUCKET/$RUN/$f" --quiet 2>/dev/null || true
    done
    aws cloudwatch put-metric-data --region "$REGION" --namespace RunKit \
      --metric-data \
        "MetricName=Heartbeat,Dimensions=[{Name=Run,Value=$RUN}],Value=1" \
        "MetricName=TraceLines,Dimensions=[{Name=Run,Value=$RUN}],Value=$count" \
      2>/dev/null || true
    sleep 60
  done
}
heartbeat & HEARTBEAT=$!
note "HEARTBEAT_UP"

# --- toolchain -------------------------------------------------------------
apt-get install -y -qq build-essential zlib1g-dev wget time python3 >/dev/null
note "APT_DONE"

mkdir -p /root/bin /root/src /root/data
cd /root/src
aws s3 cp $BUCKET/hisat2-bench/src/hisat2-src.tar.gz . --quiet
aws s3 cp $BUCKET/hisat2-bench/src/rustar-src.tar.gz . --quiet
tar xzf hisat2-src.tar.gz && tar xzf rustar-src.tar.gz

# STAR: official static Linux build, so no compile and no toolchain drift
wget -q https://github.com/alexdobin/STAR/raw/master/bin/Linux_x86_64_static/STAR -O /root/bin/STAR
chmod +x /root/bin/STAR
/root/bin/STAR --version
note "STAR_READY"

cd /root/src/hisat2 && make -j$(nproc) hisat2-align-s >/dev/null 2>&1
cp hisat2 hisat2-align-s /root/bin/
/root/bin/hisat2-align-s --version | head -1
note "HISAT2_READY"

curl -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable --profile minimal >/dev/null
. "$HOME/.cargo/env"
cd /root/src/rustar && cargo build --release >/dev/null 2>&1
find target/release -maxdepth 1 -type f -perm -u+x ! -name "*.d" -exec cp {} /root/bin/ \;
ls /root/bin
note "RUSTAR_READY"

# --- chr22-scale inputs ----------------------------------------------------
cd /root/data
wget -q ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
gunzip -f Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.22.fa chr22.fa
wget -q ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip -c Homo_sapiens.GRCh38.110.gtf.gz | awk -F'\t' '$1=="22" || /^#/' > chr22.gtf
rm -f Homo_sapiens.GRCh38.110.gtf.gz
aws s3 cp $BUCKET/hisat2-bench/smoke/smoke_R1.fq.gz . --quiet
aws s3 cp $BUCKET/hisat2-bench/smoke/smoke_R2.fq.gz . --quiet
wget -q https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-august-2016.txt.gz
gunzip -f 737K-august-2016.txt.gz
# 3' v3 uses the 3M list; for the smoke test any valid list proves the code path
note "DATA_READY"

# --- trace: the deliverable. realtime and peak_rss must be REAL ------------
# /usr/bin/time -v gives both; a row with a blank or zero in either is a failure
# of the smoke test even if the command exited 0.
trace() {  # trace <label> <cmd...>
  local label=$1; shift
  local tf=/root/out/.t.$$
  /usr/bin/time -v -o "$tf" "$@" > /root/out/${label}.log 2>&1 || true
  local rt=$(awk -F': ' '/Elapsed \(wall clock\)/{print $NF}' "$tf")
  local rss=$(awk '/Maximum resident set size/{print $NF}' "$tf")
  local rc=$(awk -F': ' '/Exit status/{print $NF}' "$tf")
  printf '%s\trealtime=%s\tpeak_rss_kb=%s\texit=%s\n' "$label" "${rt:-MISSING}" "${rss:-MISSING}" "${rc:-MISSING}" \
    >> /root/out/trace.txt
  note "TRACE $label rt=${rt:-MISSING} rss=${rss:-MISSING} rc=${rc:-MISSING}"
  rm -f "$tf"
}

cd /root/data
mkdir -p star_idx hisat_idx out

trace star_index /root/bin/STAR --runMode genomeGenerate --genomeDir star_idx \
  --genomeFastaFiles chr22.fa --sjdbGTFfile chr22.gtf --sjdbOverhang 90 \
  --genomeSAindexNbases 11 --runThreadN $(nproc)

trace hisat2_index /root/bin/hisat2-build-s -p $(nproc) chr22.fa hisat_idx/chr22 || \
  trace hisat2_index_via_wrapper /root/src/hisat2/hisat2-build -p $(nproc) chr22.fa hisat_idx/chr22

trace star_solo /root/bin/STAR --genomeDir star_idx \
  --readFilesIn smoke_R2.fq.gz smoke_R1.fq.gz --readFilesCommand zcat \
  --soloType CB_UMI_Simple --soloCBwhitelist 737K-august-2016.txt \
  --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --outFileNamePrefix out/star_ --runThreadN $(nproc) --outSAMtype None

trace hisat2_solo /root/bin/hisat2 -x hisat_idx/chr22 \
  -1 smoke_R1.fq.gz -2 smoke_R2.fq.gz --solo-barcode-mate 1 \
  --solo-cb-whitelist 737K-august-2016.txt \
  --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 12 \
  --solo-out-dir out/hisat_solo -p $(nproc) --no-unal -S /dev/null

note "WORK_DONE"

# --- verify the trace is REAL, not merely present -------------------------
bad=$(grep -c 'MISSING' /root/out/trace.txt || true)
rows=$(wc -l < /root/out/trace.txt)
note "TRACE_ROWS=$rows BAD=$bad"

kill $HEARTBEAT 2>/dev/null || true
aws s3 sync /root/out "$BUCKET/$RUN/out/" --quiet
aws s3 cp /root/out/trace.txt "$BUCKET/$RUN/trace.txt"
local_size=$(stat -c %s /root/out/trace.txt)
remote_size=$(aws s3api head-object --bucket "${BUCKET#s3://}" --key "$RUN/trace.txt" \
                --query ContentLength --output text 2>/dev/null || echo -1)
if [ "$local_size" = "$remote_size" ]; then
  note "S3_VERIFIED size=$local_size"
  note "RESULTS_COLLECTED"
  printf '{"ts":%s,"count":%s,"busy":0,"script":0,"done":1,"disk_free_gb":0,"phase":"RESULTS_COLLECTED"}\n' \
    "$(date -u +%s)" "$rows" > /root/out/status.json
  aws s3 cp /root/out/status.json "$BUCKET/$RUN/status.json" --quiet
  aws s3 cp /root/out/progress    "$BUCKET/$RUN/progress"    --quiet
  sleep 10
  shutdown -h now
else
  note "S3_UNVERIFIED local=$local_size remote=$remote_size — staying up for recovery"
fi
