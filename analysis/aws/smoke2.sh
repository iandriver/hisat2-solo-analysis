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
RUN=hisat2-smoke2
REGION=us-east-1
MAIN=$$
mkdir -p /root/out2
: > /root/out2/trace.txt
note() { echo "$(date -u +%H:%M:%S) $*" >> /root/out2/progress; }
fail() { rc=$?; note "STAGE_FAILED rc=$rc"; kill ${HEARTBEAT:-0} 2>/dev/null || true;
         aws s3 cp /root/out2/progress "$BUCKET/$RUN/progress" --quiet 2>/dev/null || true;
         exit $rc; }
trap fail ERR
note "STAGE_START"
aws --version
note "AWSCLI_READY"

heartbeat() {
  set +e
  while true; do
    ts=$(date -u +%s)
    count=$(($(wc -l < /root/out2/trace.txt 2>/dev/null || echo 0)))
    script=$(kill -0 "$MAIN" 2>/dev/null && echo 1 || echo 0)
    done_flag=$(($(grep -c RESULTS_COLLECTED /root/out2/progress 2>/dev/null || true)))
    phase=$(tail -1 /root/out2/progress 2>/dev/null | tr -d '"\\' | tr -cd '[:print:]')
    printf '{"ts":%s,"count":%s,"busy":0,"script":%s,"done":%s,"disk_free_gb":0,"phase":"%s"}\n' \
      "$ts" "$count" "$script" "$done_flag" "$phase" > /root/out2/status.json
    echo "$ts,$count,0" >> /root/out2/heartbeat.csv
    for f in trace.txt progress status.json heartbeat.csv; do
      aws s3 cp /root/out2/$f "$BUCKET/$RUN/$f" --quiet 2>/dev/null || true
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
  local tf=/root/out2/.t.$$
  local rc=0
  /usr/bin/time -v -o "$tf" "$@" > /root/out2/${label}.log 2>&1 || rc=$?
  local rt=$(awk -F': ' '/Elapsed \(wall clock\)/{print $NF}' "$tf" 2>/dev/null)
  local rss=$(awk '/Maximum resident set size/{print $NF}' "$tf" 2>/dev/null)
  printf '%s\trealtime=%s\tpeak_rss_kb=%s\texit=%s\n' \
    "$label" "${rt:-MISSING}" "${rss:-MISSING}" "$rc" >> /root/out2/trace.txt
  note "TRACE $label rt=${rt:-MISSING} rss=${rss:-MISSING} rc=$rc"
  rm -f "$tf"
  return 0                 # collect every row; the gate below decides pass/fail
}

cd /root/data
rm -rf hisat_idx out2; mkdir -p hisat_idx out2

trace hisat2_index /root/bin/hisat2-build-s -p $(nproc) chr22.fa hisat_idx/chr22

trace hisat2_solo /root/bin/hisat2 -x hisat_idx/chr22 \
  -1 smoke_R1.fq.gz -2 smoke_R2.fq.gz --solo-barcode-mate 1 \
  --solo-cb-whitelist 737K-august-2016.txt \
  --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 12 \
  --solo-out-dir out2/hisat_solo -p $(nproc) --no-unal -S /dev/null

note "WORK_DONE"

# --- the gate: real timings AND every command exited 0 --------------------
rows=$(wc -l < /root/out2/trace.txt)
missing=$(grep -c 'MISSING' /root/out2/trace.txt || true)
failed=$(awk -F'exit=' '$2 != 0 {n++} END{print n+0}' /root/out2/trace.txt)
matrix=$(ls out2/hisat_solo/Gene/raw/matrix.mtx 2>/dev/null | wc -l)
note "GATE rows=$rows missing=$missing failed=$failed matrix=$matrix"
head -3 out2/hisat_solo/Gene/Summary.csv >> /root/out2/progress 2>/dev/null || true

kill $HEARTBEAT 2>/dev/null || true
aws s3 sync /root/out2 "$BUCKET/$RUN/out/" --quiet
aws s3 cp /root/out2/trace.txt "$BUCKET/$RUN/trace.txt"
local_size=$(stat -c %s /root/out2/trace.txt)
remote_size=$(aws s3api head-object --bucket "${BUCKET#s3://}" --key "$RUN/trace.txt" \
                --query ContentLength --output text 2>/dev/null || echo -1)
if [ "$local_size" != "$remote_size" ]; then
  note "S3_UNVERIFIED local=$local_size remote=$remote_size — staying up"
  exit 1
fi
note "S3_VERIFIED size=$local_size"
if [ "$missing" -eq 0 ] && [ "$failed" -eq 0 ] && [ "$matrix" -eq 1 ]; then
  note "SMOKE_PASS"
  note "RESULTS_COLLECTED"
else
  note "SMOKE_FAIL missing=$missing failed=$failed matrix=$matrix"
fi
printf '{"ts":%s,"count":%s,"busy":0,"script":0,"done":1,"disk_free_gb":0,"phase":"%s"}\n' \
  "$(date -u +%s)" "$rows" "$(tail -1 /root/out2/progress | tr -d '"\\')" > /root/out2/status.json
aws s3 cp /root/out2/status.json "$BUCKET/$RUN/status.json" --quiet
aws s3 cp /root/out2/progress    "$BUCKET/$RUN/progress"    --quiet
sleep 10
shutdown -h now
