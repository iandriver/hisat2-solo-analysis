#!/bin/bash
# Phase 1, second pass. Reuses the preserved volume: toolchain, chr22 inputs and
# the STAR index are already on disk from the first run.
#
# Two bugs from pass 1, both fixed here:
#   1. only hisat2-align-s was built, so hisat2-build-s was command-not-found.
#   2. trace() swallowed failure -- `|| true` inside it meant a dead command still
#      produced a well-formed row, and the gate only checked that the timing
#      fields were non-empty. It reported BAD=0 with two rc!=0 rows in the trace.
#      A trace row is now only acceptable if exit=0.
set -eux
export HOME=/root
export DEBIAN_FRONTEND=noninteractive
BUCKET=s3://rustar-bench
RUN=hisat2-smoke4
REGION=us-east-1
MAIN=$$
mkdir -p /root/out4
: > /root/out4/trace.txt
note() { echo "$(date -u +%H:%M:%S) $*" >> /root/out4/progress; }
fail() { rc=$?; note "STAGE_FAILED rc=$rc"; kill ${HEARTBEAT:-0} 2>/dev/null || true;
         aws s3 cp /root/out4/progress "$BUCKET/$RUN/progress" --quiet 2>/dev/null || true;
         exit $rc; }
trap fail ERR
note "STAGE_START"
aws --version
note "AWSCLI_READY"

heartbeat() {
  set +e
  while true; do
    ts=$(date -u +%s)
    count=$(($(wc -l < /root/out4/trace.txt 2>/dev/null || echo 0)))
    script=$(kill -0 "$MAIN" 2>/dev/null && echo 1 || echo 0)
    done_flag=$(($(grep -c RESULTS_COLLECTED /root/out4/progress 2>/dev/null || true)))
    phase=$(tail -1 /root/out4/progress 2>/dev/null | tr -d '"\\' | tr -cd '[:print:]')
    printf '{"ts":%s,"count":%s,"busy":0,"script":%s,"done":%s,"disk_free_gb":0,"phase":"%s"}\n' \
      "$ts" "$count" "$script" "$done_flag" "$phase" > /root/out4/status.json
    echo "$ts,$count,0" >> /root/out4/heartbeat.csv
    for f in trace.txt progress status.json heartbeat.csv; do
      aws s3 cp /root/out4/$f "$BUCKET/$RUN/$f" --quiet 2>/dev/null || true
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

# --- fix 1: build the binaries actually needed -----------------------------
cd /root/src/hisat2
make -j$(nproc) hisat2-align-s hisat2-build-s hisat2-inspect-s >/dev/null 2>&1
cp hisat2 hisat2-align-s hisat2-build-s hisat2-inspect-s /root/bin/
cp hisat2-build hisat2-inspect /root/bin/ 2>/dev/null || true
chmod +x /root/bin/*
/root/bin/hisat2-build-s --version | head -1
note "HISAT2_BUILD_READY"

# --- fix 2: a trace row counts only if the command actually succeeded ------
trace() {  # trace <label> <cmd...>
  local label=$1; shift
  local tf=/root/out4/.t.$$
  local rc=0
  /usr/bin/time -v -o "$tf" "$@" > /root/out4/${label}.log 2>&1 || rc=$?
  local rt=$(awk -F': ' '/Elapsed \(wall clock\)/{print $NF}' "$tf" 2>/dev/null)
  local rss=$(awk '/Maximum resident set size/{print $NF}' "$tf" 2>/dev/null)
  printf '%s\trealtime=%s\tpeak_rss_kb=%s\texit=%s\n' \
    "$label" "${rt:-MISSING}" "${rss:-MISSING}" "$rc" >> /root/out4/trace.txt
  note "TRACE $label rt=${rt:-MISSING} rss=${rss:-MISSING} rc=$rc"
  rm -f "$tf"
  return 0                 # collect every row; the gate below decides pass/fail
}

cd /root/data
rm -rf out4; mkdir -p out4
aws s3 cp 's3://rustar-bench/hisat2-bench/refs/whitelist_v3.txt' /root/data/whitelist_v3.txt --quiet

# index already built and traced in pass 2 (0:34.21, 191964 kB)
trace gene_model python3 /root/src/hisat2/hisat2_extract_genes.py chr22.gtf -o chr22.ht2gm

trace hisat2_solo /root/bin/hisat2 -x hisat_idx/chr22 \
  -1 smoke_R1.fq.gz -2 smoke_R2.fq.gz --solo-barcode-mate 1 \
  --solo-cb-whitelist whitelist_v3.txt \
  --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 12 \
  --gene-annotation chr22.ht2gm --gene-strand Forward \
  --solo-out-dir out4/hisat_solo -p $(nproc) --no-unal -S /dev/null

note "WORK_DONE"

# --- the gate: real timings AND every command exited 0 --------------------
rows=$(wc -l < /root/out4/trace.txt)
missing=$(grep -c 'MISSING' /root/out4/trace.txt || true)
failed=$(awk -F'exit=' '$2 != 0 {n++} END{print n+0}' /root/out4/trace.txt)
matrix=$(ls out4/hisat_solo/Gene/raw/matrix.mtx 2>/dev/null | wc -l)
# A wrong barcode whitelist still exits 0 and still writes a matrix -- it just
# counts almost nothing. Pass 3 scored 2.2% valid barcodes against the 737K list
# on v3 data. Make the rate part of the gate, not something read afterwards.
cbrate=$(awk -F, '/Valid Barcodes/{print ($2>0.80)?1:0}' out4/hisat_solo/Gene/Summary.csv 2>/dev/null)
note "GATE rows=$rows missing=$missing failed=$failed matrix=$matrix cbrate_ok=${cbrate:-0}"
head -3 out4/hisat_solo/Gene/Summary.csv >> /root/out4/progress 2>/dev/null || true

kill $HEARTBEAT 2>/dev/null || true
aws s3 sync /root/out4 "$BUCKET/$RUN/out/" --quiet
aws s3 cp /root/out4/trace.txt "$BUCKET/$RUN/trace.txt"
local_size=$(stat -c %s /root/out4/trace.txt)
remote_size=$(aws s3api head-object --bucket "${BUCKET#s3://}" --key "$RUN/trace.txt" \
                --query ContentLength --output text 2>/dev/null || echo -1)
if [ "$local_size" != "$remote_size" ]; then
  note "S3_UNVERIFIED local=$local_size remote=$remote_size — staying up"
  exit 1
fi
note "S3_VERIFIED size=$local_size"
if [ "$missing" -eq 0 ] && [ "$failed" -eq 0 ] && [ "$matrix" -eq 1 ] && [ "${cbrate:-0}" -eq 1 ]; then
  note "SMOKE_PASS"
  note "RESULTS_COLLECTED"
else
  note "SMOKE_FAIL missing=$missing failed=$failed matrix=$matrix cbrate_ok=${cbrate:-0}"
fi
printf '{"ts":%s,"count":%s,"busy":0,"script":0,"done":1,"disk_free_gb":0,"phase":"%s"}\n' \
  "$(date -u +%s)" "$rows" "$(tail -1 /root/out4/progress | tr -d '"\\')" > /root/out4/status.json
aws s3 cp /root/out4/status.json "$BUCKET/$RUN/status.json" --quiet
aws s3 cp /root/out4/progress    "$BUCKET/$RUN/progress"    --quiet
sleep 10
shutdown -h now
