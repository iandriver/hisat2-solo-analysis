#!/bin/bash
# E0 -- whole-genome 32-bit graph build, instrumented for the growth curve.
#
# The question: per-chromosome traces extrapolate to ~3.43e9 path nodes (80% of
# 2^32), yet the real whole-genome build overflows. Either growth is superlinear
# in the number of sequences, or the extrapolation is wrong somewhere else.
# --verbose logs "Generation N (temp -> nodes, ranks)" per doubling step, so the
# curve right up to the overflow answers it.
#
# This build is EXPECTED to fail with "exceeded integer bounds". That is the
# measurement, not an error. What matters is the curve it prints first.
#
# Cost discipline, learned the hard way on stage 2b: a 58-second failure once sat
# idle for 2h20m at $2.12/hr because the failure path published to a location
# nothing was watching. Here every exit path writes to the watched key and then
# powers off, and the instance is launched with shutdown-behavior=terminate.
set -uo pipefail
exec > /var/log/e0.log 2>&1
set -x

BUCKET=s3://rustar-bench/hisat2-curve
D=/mnt/work
NP=$(nproc)

# Dead man's switch. If anything below wedges, the box still dies.
shutdown -h +330 &

note(){ echo "[$(date -u +%H:%M:%S)] $*" | tee -a $D/progress; aws s3 cp $D/progress $BUCKET/progress --only-show-errors || true; }

fail(){ rc=$?
  echo "STAGE_FAILED rc=$rc line=${BASH_LINENO[0]}" >> $D/progress
  aws s3 cp $D/progress $BUCKET/progress --only-show-errors || true
  aws s3 cp /var/log/e0.log $BUCKET/e0.log --only-show-errors || true
  [ -f $D/build.err ] && aws s3 cp $D/build.err $BUCKET/build.err --only-show-errors || true
  sleep 5
  shutdown -h now
  exit $rc; }
trap fail ERR

mkdir -p $D && cd $D
echo "starting" > $D/progress

dnf -y install gcc-c++ make git python3 pigz tar zlib-devel >/dev/null
note "DEPS_READY $(nproc) cpus, $(free -g | awk '/^Mem/{print $2}') GB RAM"

# ---- hisat2 (the fork, matching what stage 2b ran) ----
git clone -q --depth 1 -b upstream/modernize https://github.com/iandriver/hisat2 $D/hisat2
cd $D/hisat2
make -j$NP hisat2-build-s >/dev/null 2>&1
$D/hisat2/hisat2-build --version >/dev/null 2>&1 || { note 'WRAPPER_BROKEN'; false; }
grep -q selectAlts hgfm.h || { note "WRONG_BUILD selectAlts missing"; false; }
note "HISAT2_READY"

# ---- inputs ----
cd $D
curl -sSL -o genome.fa.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz
pigz -dc genome.fa.gz > genome.fa
BP=$(grep -v '^>' genome.fa | tr -d '\n' | wc -c)
CONTIGS=$(grep -c '^>' genome.fa)
note "GENOME bp=$BP contigs=$CONTIGS"

aws s3 cp $BUCKET/genome.snp.gz . --only-show-errors
aws s3 cp $BUCKET/genome.haplotype.gz . --only-show-errors
pigz -dc genome.snp.gz > genome.snp
pigz -dc genome.haplotype.gz > genome.haplotype
NV=$(wc -l < genome.snp); NH=$(wc -l < genome.haplotype)
[ "$NV" -gt 14000000 ] || { note "SNP_FILE_SHORT $NV"; false; }
note "VARIANTS snp=$NV haplotype=$NH"

mkdir -p $D/idx $D/smoke $D/smoke/idx
head -20000 $D/genome.fa > $D/smoke/smoke.fa
awk -F'\t' '$3=="chr1" && $4 < 1100000' $D/genome.snp       > $D/smoke/smoke.snp
awk -F'\t' '$2=="chr1" && $3 < 1100000' $D/genome.haplotype > $D/smoke/smoke.haplotype
note "SMOKE_INPUT snp=$(wc -l < $D/smoke/smoke.snp) hap=$(wc -l < $D/smoke/smoke.haplotype)"
trap - ERR
$D/hisat2/hisat2-build -p 4 --snp $D/smoke/smoke.snp --haplotype $D/smoke/smoke.haplotype \
    $D/smoke/smoke.fa $D/smoke/idx/smoke > $D/smoke/out 2> $D/smoke/err
SRC=$?
trap fail ERR
SGEN=$(grep -cE '^Generation ' $D/smoke/err || true)
note "SMOKE rc=$SRC generations=$SGEN"
if [ "$SRC" -ne 0 ] || [ "$SGEN" -lt 3 ]; then
  note "SMOKE_FAILED -- aborting before the long run"
  tail -20 $D/smoke/err >> $D/progress
  aws s3 cp $D/smoke/err $BUCKET/smoke.err --only-show-errors || true
  false
fi

# ---- the run ----
# Publish the generation curve as it appears, so the answer survives even if the
# instance dies unexpectedly.
( while sleep 60; do
    grep -E '^Generation ' $D/build.err > $D/generations.txt 2>/dev/null || true
    aws s3 cp $D/generations.txt $BUCKET/generations.txt --only-show-errors 2>/dev/null || true
    tail -c 200000 $D/build.err > $D/build.tail 2>/dev/null || true
    aws s3 cp $D/build.tail $BUCKET/build.tail --only-show-errors 2>/dev/null || true
    if [ ! -s $D/generations.txt ] && [ -s $D/build.err ]; then
      echo "[$(date -u +%H:%M:%S)] WARN no Generation lines yet" >> $D/progress
      aws s3 cp $D/progress $BUCKET/progress --only-show-errors 2>/dev/null || true
    fi
  done ) &
STREAM=$!

note "BUILD_START np=$NP"
# The overflow IS the measurement -- a nonzero rc must not trip the ERR trap.
trap - ERR
set +e
/usr/bin/time -v $D/hisat2/hisat2-build -p $NP --verbose \
    --snp $D/genome.snp --haplotype $D/genome.haplotype \
    $D/genome.fa $D/idx/genome > $D/build.out 2> $D/build.err
RC=$?
trap fail ERR
kill $STREAM 2>/dev/null || true

grep -E '^Generation ' $D/build.err > $D/generations.txt || true
note "BUILD_DONE rc=$RC generations=$(wc -l < $D/generations.txt)"
if grep -q 'exceeded integer bounds' $D/build.err; then note 'OVERFLOWED (the expected result)'
elif [ "$RC" -eq 0 ]; then note 'COMPLETED -- no overflow, index built'
else note "FAILED_OTHER rc=$RC: $(grep -iE 'Error|Could not|exception' $D/build.err | head -2 | tr '\n' ' ')"
fi
grep -E 'Maximum resident set size' $D/build.err | tail -1 >> $D/progress || true

for f in generations.txt build.err build.out progress; do
  [ -f $D/$f ] && aws s3 cp $D/$f $BUCKET/$f --only-show-errors || true
done
aws s3 cp /var/log/e0.log $BUCKET/e0.log --only-show-errors || true
note "ALL_DONE"
sleep 10
shutdown -h now
