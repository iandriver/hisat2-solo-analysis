#!/bin/bash
# Stage the benchmark inputs. Each step is guarded so a rerun after a failure
# resumes rather than redoing 30 minutes of work.
set -euo pipefail
B="$(cd "$(dirname "$0")" && pwd)"
H=/Users/iandriver/Downloads/hisat2
say(){ echo "[$(date -u +%H:%M:%S)] $*"; }
free_gb(){ df -g /System/Volumes/Data | tail -1 | awk '{print $4}'; }

need(){ # need <gb> <what>
  local f; f=$(free_gb)
  [ "$f" -ge "$1" ] || { echo "ABORT: $2 needs ~${1}GB, only ${f}GB free"; exit 1; }
}

# ---- fastqs ----------------------------------------------------------------
if [ ! -f "$B/fq/.extracted" ]; then
  need 12 "fastq extraction"
  say "extracting fastqs"
  tar xf "$B/fq/pbmc_1k_v3_fastqs.tar" -C "$B/fq"
  # The tar nests them one directory down; hoist so the paths are predictable.
  find "$B/fq" -name '*_R[12]_001.fastq.gz' -exec mv -n {} "$B/fq/" \; 2>/dev/null || true
  touch "$B/fq/.extracted"
  rm -f "$B/fq/pbmc_1k_v3_fastqs.tar" "$B/fq/pbmc_1k_v3_fastqs.tar.done"
  say "fastqs done, tar removed, $(free_gb)GB free"
fi
ls "$B"/fq/*_R1_001.fastq.gz >/dev/null || { echo "ABORT: no R1 fastqs"; exit 1; }

# ---- reference text --------------------------------------------------------
if [ ! -f "$B/ref/gencode.v50.gtf" ]; then
  say "unpacking GTF"
  gzip -dc "$B/ref/gencode.v50.gtf.gz" > "$B/ref/gencode.v50.gtf"
fi
if [ ! -f "$B/ref/GRCh38.fa" ]; then
  need 10 "genome fasta"
  say "unpacking genome fasta"
  gzip -dc "$B/ref/GRCh38.primary_assembly.fa.gz" > "$B/ref/GRCh38.fa"
fi

# ---- hisat2 gene model -----------------------------------------------------
if [ ! -f "$B/ref/genes.ht2gm" ]; then
  say "building .ht2gm gene model"
  python3 "$H/hisat2_extract_genes.py" "$B/ref/gencode.v50.gtf" > "$B/ref/genes.ht2gm"
fi
say "gene model: $(grep -c '^G' "$B/ref/genes.ht2gm") genes"

# ---- hisat2 graph index ----------------------------------------------------
if [ ! -f "$B/idx/.snp_extracted" ]; then
  need 20 "hisat2 index extraction"
  say "extracting grch38_snp"
  tar xzf "$B/idx/grch38_snp.tar.gz" -C "$B/idx"
  touch "$B/idx/.snp_extracted"
  rm -f "$B/idx/grch38_snp.tar.gz" "$B/idx/grch38_snp.tar.gz.done"
fi
HT2_PREFIX=$(ls "$B"/idx/grch38_snp/*.1.ht2 2>/dev/null | head -1 | sed 's/\.1\.ht2$//')
[ -n "$HT2_PREFIX" ] || { echo "ABORT: no .1.ht2 under $B/idx/grch38_snp"; exit 1; }
echo "$HT2_PREFIX" > "$B/idx/ht2_prefix"
say "hisat2 index prefix: $HT2_PREFIX"

say "PREP COMPLETE, $(free_gb)GB free"
