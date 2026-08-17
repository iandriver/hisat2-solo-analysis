#!/bin/bash
# Downloads for the human head-to-head. Resumable (-C -) so a dropped
# connection costs the remainder, not the whole file. Each writes a .done
# stamp only after an integrity check, so a truncated file is never reused.
set -uo pipefail
B="$(cd "$(dirname "$0")" && pwd)"
log(){ echo "[$(date -u +%H:%M:%S)] $*" >> "$B/logs/dl.log"; }

get(){ # url dest minbytes
  local url=$1 dest=$2 min=$3
  [ -f "$dest.done" ] && { log "skip $(basename $dest)"; return 0; }
  for try in 1 2 3 4 5; do
    curl -sSL -C - --max-time 7200 --retry 3 -o "$dest" "$url" 2>>"$B/logs/dl.log" && {
      local sz; sz=$(stat -f %z "$dest" 2>/dev/null || echo 0)
      if [ "$sz" -ge "$min" ]; then
        # integrity, not size: a truncated .gz passes a size floor
        case "$dest" in
          *.gz) gzip -t "$dest" 2>/dev/null || { log "corrupt gz $(basename $dest) try=$try"; rm -f "$dest"; continue; } ;;
          *.tar) tar tf "$dest" >/dev/null 2>&1 || { log "corrupt tar $(basename $dest) try=$try"; rm -f "$dest"; continue; } ;;
        esac
        touch "$dest.done"; log "ok $(basename $dest) $sz bytes"; return 0
      fi
      log "short $(basename $dest) $sz < $min try=$try"
    }
    sleep 10
  done
  log "FAILED $(basename $dest)"; return 1
}

get https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz \
    "$B/ref/gencode.v50.gtf.gz" 40000000 &
get https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz \
    "$B/ref/GRCh38.primary_assembly.fa.gz" 800000000 &
get https://genome-idx.s3.amazonaws.com/hisat/grch38_snp.tar.gz \
    "$B/idx/grch38_snp.tar.gz" 3000000000 &
get https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar \
    "$B/fq/pbmc_1k_v3_fastqs.tar" 5000000000 &
wait
log "ALL DOWNLOADS SETTLED"
