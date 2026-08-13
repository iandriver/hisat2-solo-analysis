# Minimal AWS run: close the human head-to-head gap

## What is actually missing

Everything about HISAT2-solo on human is measured. Nothing about its competitors on
human is. Specifically:

| | measured? |
|---|---|
| HISAT2-solo human: 382 s / 9.31 GB, floor 8–10 GB | yes, but on an **M5 Pro Mac** |
| STARsolo on human | **never run** — no index, needs ~64 GB to build |
| rustar on human | **never run** — same |
| rustar 25.6–28.2 GB | yes, but that is **mouse** |
| STAR ~32 GB for human | **documented, not measured** |

Two separate problems: the competitors have no human numbers at all, and our numbers
come from a machine none of theirs ever ran on. Both are fixed by one run: **all three
tools, same box, same data, matched concurrency.**

## Deliverables

| | what it settles |
|---|---|
| **D1** peak RSS, all three, human, same box | turns "2.5–3x less memory" from inference into measurement |
| **D2** wall time, all three, matched `-p`/`--runThreadN` | the speed claim, currently "a tie" on mouse only |
| **D3** hard memory cap per tool (cgroup) | which tools survive at 32/16/12/8 GB — T3 with competitors present |
| **D4** count concordance on human | we have ARI 0.9315 vs STARsolo on **mouse** only |
| **D5** index build cost | STAR/rustar RAM and wall to build; a secondary result, not part of D2 |

## Constraint that sets the instance

STAR `genomeGenerate` for human with a GTF is the binding requirement — commonly
~32 GB, comfortably 64 GB. Everything else is smaller. Disk: STAR index ~30 GB,
rustar ~25 GB, HISAT2 prebuilt 11 GB, genome+GTF 5 GB, FASTQs 5 GB, outputs ~10 GB
≈ **90 GB**, so 200 GB gp3.

## Two plans

Live on-demand/spot pricing, us-west-2, via `launch.hourly()`:

| type | vCPU | RAM | on-demand | spot |
|---|---|---|---|---|
| t3.large | 2 | 8 GB | $0.083 | $0.032 |
| r7i.2xlarge | 8 | 64 GB | $0.529 | $0.192 |
| **r7i.4xlarge** | 16 | 128 GB | **$1.058** | $0.431 |

### Plan A — single on-demand box (recommended)

Phase 1 on t3.large, chr22-scale, proves the script and the S3 path (runbook rule 1).
Phase 2 on r7i.4xlarge does everything.

| | hours | rate | cost |
|---|---|---|---|
| Phase 1 t3.large | 2 | $0.083 | $0.17 |
| Phase 2 r7i.4xlarge | 5 | $1.058 | $5.29 |
| EBS 200 GB gp3, 7 h | | | $0.15 |
| S3 + egress (~2 GB out) | | | $0.20 |
| | | **total** | **≈ $5.80** |

Phase 2 breakdown: install + build rustar 20 min; fetch 21 GB of inputs 20 min;
STAR index 45–60 min; rustar index 30–45 min; three aligner runs 20 min; cap matrix
90 min; upload 5 min. **~4.5 h, budget 5.**

### Plan B — split spot build / on-demand timing

Index building is not a timing result, so it can run on spot. Only the measurement
needs on-demand.

| | hours | rate | cost |
|---|---|---|---|
| Phase 1 t3.large | 2 | $0.083 | $0.17 |
| Build indexes, r7i.4xlarge **spot** | 2.5 | $0.431 | $1.08 |
| Time the runs, r7i.4xlarge on-demand | 1.5 | $1.058 | $1.59 |
| EBS + S3 + egress | | | $0.35 |
| | | **total** | **≈ $3.20** |

Saves ~$2.60 for one extra boot and an S3 round-trip of ~55 GB of indexes.

**Restoring indexes from S3 is safe for a timing run** — FINDINGS #5 is about EBS
snapshot lazy-loading, and an S3 download writes real blocks to a fresh volume. What
would *not* be safe is resuming the timing phase from a snapshot.

### Why not spot for the timing phase

FINDINGS #7 says spot is fine when work is resumable; FINDINGS #5 says a
snapshot-restored volume invalidates timing. Together they rule out spot for D2 —
recovery from an interruption is exactly the thing that would corrupt the number.
Spot would save $0.94 on a $5.80 run and put the primary deliverable at risk.

## Run outline

```
setup      apt: build-essential zlib1g-dev; STAR binary from GitHub releases
           cargo build --release rustar-aligner (github.com/iandriver, 5c020ff)
           make hisat2 (upstream/modernize @ d20f6b3 — carries -fsigned-char,
                        without which it does not compile on ARM; x86 here, but
                        build from the same tree we are claiming numbers for)
inputs     GRCh38 primary FASTA + Ensembl GTF; JHU prebuilt grch38 + grch38_snp;
           10x pbmc_1k_v3 FASTQs (the exact dataset every existing number uses)
D5         STAR genomeGenerate, rustar index — record /usr/bin/time -v each
D1/D2      each tool on 66.6M reads at --runThreadN/-p 16, 3 replicates,
           /usr/bin/time -v; HISAT2 twice (grch38_snp and grch38)
D3         re-run each under systemd-run --scope -p MemoryMax=N, N in
           {32,24,16,12,10,8} GB, swap off; record completion or OOM
D4         Solo.out matrices from all three -> S3; ARI computed locally
```

Everything streams to S3 as it completes (rule 3), heartbeat every 60 s, CloudWatch
alarm armed before the run starts, box stops itself on completion.

## What this does not buy

- It does not make HISAT2 faster or fix the pseudogene defect (#472).
- D3 will likely show STARsolo failing at 32 GB and HISAT2 surviving to 8–10 GB, but
  if STAR's shared-memory mode or a reduced `--genomeSAsparseD` gets it under, that is
  the honest answer and it weakens the memory claim.
- The Mac numbers stay non-comparable. After this run, **quote the AWS numbers for all
  three and retire the Mac figures from any comparative claim.**

## Housekeeping found while pricing

`teardown.sweep()` reports a pre-existing **stopped t4g.micro and an in-use 20 GB
volume in us-east-1** — not from this project. Stopped instances do not bill compute
but that volume does, ~$1.60/month. Left alone.
