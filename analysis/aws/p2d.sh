#!/bin/bash
# Phase 2d -- the graph index, 32-bit, SNVs only.
#
# Two failed attempts bound the problem from both sides:
#   32-bit, 14.95M variants -> "exceeded integer bounds" (over 2^32 graph nodes)
#   64-bit, 14.95M variants -> OOM-killed at 368.8 GB in PathGraph construction
# The kill happened after FINISHED RECURSIVE SORTS, i.e. in graph construction,
# which has no --bmax-style knob -- so the only levers were more RAM (unbounded
# above 368.8 GB) or fewer variants.
#
# Dropping the 2.3M indels leaves 12,643,877 SNVs, within 2.8% of the shipped
# grch38_snp's 12.3M, which is proof by existence that a set this size fits.
# Insertions alone contributed 3,013,044 extra nodes, so this removes more graph
# than it removes variants. Filtered and validated off-box (no dangling
# haplotype references, no empty haplotypes, spans retightened).
set -eux
set -o pipefail
export HOME=/root
export DEBIAN_FRONTEND=noninteractive
BUCKET=s3://rustar-bench
RUN=hisat2-p2d
REGION=us-east-1
MAIN=$$
OUT=/root/p2d
D=/root/data
T0=$(date -u +%s)
mkdir -p "$OUT" /root/bin
: > $OUT/trace.txt
: > $OUT/progress
note() { echo "$(date -u +%H:%M:%S) $*" >> $OUT/progress; }
fail() { rc=$?; note "STAGE_FAILED rc=$rc line=${BASH_LINENO[0]}"
         kill ${HEARTBEAT:-0} 2>/dev/null || true
         aws s3 sync $OUT "$BUCKET/$RUN/out/" --quiet 2>/dev/null || true
         aws s3 cp $OUT/progress  "$BUCKET/$RUN/progress"  --quiet 2>/dev/null || true
         aws s3 cp $OUT/trace.txt "$BUCKET/$RUN/trace.txt" --quiet 2>/dev/null || true
         printf '{"ts":%s,"count":0,"busy":0,"script":0,"done":1,"disk_free_gb":0,"phase":"STAGE_FAILED rc=%s"}\n' \
           "$(date -u +%s)" "$rc" > $OUT/status.json
         aws s3 cp $OUT/status.json "$BUCKET/$RUN/status.json" --quiet 2>/dev/null || true
         sleep 5; shutdown -h now; exit $rc; }
trap fail ERR
elapsed() { echo $(( $(date -u +%s) - T0 )); }
note "STAGE_START"
export PATH=/root/bin:/usr/local/bin:$PATH

heartbeat() {
  set +e
  while true; do
    ts=$(date -u +%s)
    count=$(($(wc -l < $OUT/trace.txt 2>/dev/null || echo 0)))
    script=$(kill -0 "$MAIN" 2>/dev/null && echo 1 || echo 0)
    done_flag=$(($(grep -c RESULTS_COLLECTED $OUT/progress 2>/dev/null || true)))
    free=$(df -BG --output=avail /root | tail -1 | tr -dc '0-9')
    # The build is one long opaque step, so report its live RSS: that is the
    # number this run exists to discover, and if it dies the last sample is the
    # only evidence of how close it got.
    rss=$(ps -eo rss,comm --sort=-rss | awk 'NR==2{printf "%.1f", $1/1048576}')
    phase=$(tail -1 $OUT/progress 2>/dev/null | tr -d '"\\' | tr -cd '[:print:]')
    printf '{"ts":%s,"count":%s,"busy":0,"script":%s,"done":%s,"disk_free_gb":%s,"rss_gb":%s,"phase":"%s"}\n' \
      "$ts" "$count" "$script" "$done_flag" "${free:-0}" "${rss:-0}" "$phase" > $OUT/status.json
    echo "$ts,$count,${free:-0},${rss:-0}" >> $OUT/heartbeat.csv
    for f in trace.txt progress status.json heartbeat.csv; do
      aws s3 cp $OUT/$f "$BUCKET/$RUN/$f" --quiet 2>/dev/null || true
    done
    aws cloudwatch put-metric-data --region "$REGION" --namespace RunKit \
      --metric-data "MetricName=Heartbeat,Dimensions=[{Name=Run,Value=$RUN}],Value=1" \
                    "MetricName=PeakRssGb,Dimensions=[{Name=Run,Value=$RUN}],Value=${rss:-0}" \
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

# Reclaim the aborted 32-bit attempt's scratch before starting.
rm -rf $D/idx_snp $D/idx_snp_l
note "DISK_FREE_START $(df -BG --output=avail /root | tail -1 | tr -dc '0-9')GB"

cd /root/src/hisat2
grep -q "selectAlts" hgfm.h             # the binary-search backoff must be present
note "BUILD_S_READY $(hisat2-build-s --version | head -1)"

# The SNV-only pair was filtered and validated off-box; fetch rather than redo.
cd $D
[ -s genome.snv.snp ] || { aws s3 cp "$BUCKET/hisat2-p2d/refs/genome.snv.snp.gz" - | pigz -dc > genome.snv.snp; }
[ -s genome.snv.haplotype ] || { aws s3 cp "$BUCKET/hisat2-p2d/refs/genome.snv.haplotype.gz" - | pigz -dc > genome.snv.haplotype; }
[ -s $D/genome.fa ]
NSNP=$(wc -l < genome.snv.snp); NHAP=$(wc -l < genome.snv.haplotype)
note "INPUTS snp=$NSNP hap=$NHAP"
[ "$NSNP" -eq 12643877 ]
[ "$NHAP" -eq 14005300 ]
# Indels here would mean the wrong file arrived.
[ "$(cut -f2 genome.snv.snp | sort -u | tr -d '\n')" = "single" ]

rm -rf $D/idx_snv; mkdir -p $D/idx_snv
trace build_snv hisat2-build-s -p "$(nproc)" \
  --snp $D/genome.snv.snp --haplotype $D/genome.snv.haplotype \
  $D/genome.fa $D/idx_snv/genome
RET=$(grep -o 'Local indexes: .*' $OUT/build_snv.log | tail -1)
SZ=$(du -sk $D/idx_snv | cut -f1)
printf 'SIZE\tsnv\tindex_kb=%s\t%s\n' "$SZ" "${RET:-no-variants}" >> $OUT/trace.txt
note "SIZE snv index_kb=$SZ ${RET:-no-variants}"
note "DISK_FREE $(df -BG --output=avail /root | tail -1 | tr -dc '0-9')GB"

# The gate: the build must have succeeded AND carried variants. A graph index
# with no retention line is a linear index wearing the wrong name.
core_ok=1
gawk -F'exit=' '/^build_snv\t/ && $2+0 != 0 {bad=1} END{exit bad+0}' $OUT/trace.txt || core_ok=0
grep -q 'variant instances retained' $OUT/build_snv.log || core_ok=0
ls $D/idx_snv/genome.1.ht2 >/dev/null 2>&1 || core_ok=0
note "GATE core_ok=$core_ok"

if [ "$core_ok" -eq 1 ]; then
  aws s3 sync $D/idx_snv "$BUCKET/$RUN/refs/grch38_v50_1kgp_snv/" --quiet
  note "INDEX_UPLOADED bytes=$(du -sb $D/idx_snv | cut -f1)"
fi
kill $HEARTBEAT 2>/dev/null || true
aws s3 sync $OUT "$BUCKET/$RUN/out/" --quiet
aws s3 cp $OUT/trace.txt "$BUCKET/$RUN/trace.txt"
lsz=$(stat -c %s $OUT/trace.txt)
rsz=$(aws s3api head-object --bucket "${BUCKET#s3://}" --key "$RUN/trace.txt" \
        --query ContentLength --output text 2>/dev/null || echo -1)
[ "$lsz" = "$rsz" ] || { note "S3_UNVERIFIED local=$lsz remote=$rsz"; exit 1; }
note "S3_VERIFIED size=$lsz"
if [ "$core_ok" -eq 1 ]; then note "BUILD_PASS"; note "RESULTS_COLLECTED"
else note "BUILD_FAIL core_ok=$core_ok"; fi
printf '{"ts":%s,"count":1,"busy":0,"script":0,"done":1,"disk_free_gb":0,"phase":"%s"}\n' \
  "$(date -u +%s)" "$(tail -1 $OUT/progress | tr -d '"\\')" > $OUT/status.json
aws s3 cp $OUT/status.json "$BUCKET/$RUN/status.json" --quiet
aws s3 cp $OUT/progress    "$BUCKET/$RUN/progress"    --quiet
sleep 10
shutdown -h now
