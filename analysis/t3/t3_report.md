# T3: how little memory does HISAT2-solo need for a human single-cell run?

**Answer: 8 GB at low thread counts, 10 GB at 8 threads.** The plan's "completes
at 8 GB" is right only if you turn threads down; with `-p 8` it is OOM-killed
there. The 5.18 GB figure quoted elsewhere was **mouse**; human needs an
~8 GB working set either way.

## Getting a real cap

`ulimit -v` is **not enforced on macOS** — a process under a 2 GB `ulimit -v`
allocates and touches 4 GB without complaint. The test as written in the plan
could not have worked.

The cap therefore has to be a kernel one: a colima Linux VM (20 GiB, 8 CPU)
with a per-run cgroup limit, `--memory=Ng --memory-swap=Ng` so swap is disabled
and the kernel OOM-kills rather than swapping (docker then reports exit 137).

Two things had to be fixed to get there, both findings in their own right:

- **HISAT2 2.2.3 does not compile on aarch64 Linux.** `alphabet.cpp:404`,
  `char mask2iupac[16] = { -1, ... }` — plain `char` is unsigned on aarch64
  Linux, so `-1` will not narrow into it. It builds on this Mac because Apple's
  arm64 ABI keeps `char` signed. `EXTRA_FLAGS=-fsigned-char` fixes it and is a
  no-op on Darwin and x86 Linux. This is the latent ARM breakage the plan's
  Risk 7 anticipated and belongs with the U4 build-portability work.
- **`--solo-buffer-mb` does not exist.** It is in the plan's design notes for
  the external-merge-sort spill, but was never implemented — so the only lever
  for trading memory is `-p`, which turns out to be enough here (below).

## The curve

Human PBMC 1k v3, 66.6M reads, `grch38_snp`, 8 threads, full solo counting:

| cap | result | wall | peak RSS | aligned | cells |
|---|---|---|---|---|---|
| 18 GB | **completes** | 1,024 s | 8.12 GB | 87.20% | 1,134 |
| 12 GB | **completes** | 6,735 s (**6.6x**) | 8.12 GB | 87.20% | 1,134 |
| 10 GB | **completes** | 13,621 s (**13.3x**) | 8.12 GB | 87.20% | 1,134 |
| 8 GB | **OOM-killed** | died after 4,479 s | 7.96 GB | – | – |
| 6 GB | OOM-killed | 5 s | 5.97 GB | – | – |
| 5 GB | OOM-killed | 5 s | 4.98 GB | – | – |
| 4 GB | OOM-killed | 3 s | 3.98 GB | – | – |

Output is identical at every cap that completes — 87.20% aligned, 1,134 cells,
matching the native macOS run exactly, so the containerised Linux build is a
valid stand-in.

Below 8 GB the process dies in 3-5 seconds, during index load. At 8 GB it gets
to 7.96 GB and survives 75 minutes before being killed — it very nearly fits.

## Thread count is a lever, and it cuts both ways

Rerunning the failed 8 GB cap with `-p 2` instead of `-p 8`:

| cap | threads | result | wall | peak RSS |
|---|---|---|---|---|
| 10 GB | 8 | completes | 13,621 s | 8.12 GB |
| 8 GB | 8 | **OOM-killed** | died after 4,479 s | 7.96 GB |
| 8 GB | **2** | **completes** | **5,027 s** | **7.93 GB** |

Same output as everywhere else (87.20%, 1,134 cells). Two threads saves only
~0.19 GB of peak RSS, but that is the whole margin — it is the difference
between finishing and being killed.

**And it is 2.7x faster than the 8-thread run at 10 GB**, despite using a
quarter of the threads and less memory. Under memory pressure, cutting threads
helps twice: it lowers the working set *and* it leaves more room for page cache,
so the index stops being evicted. Adding threads to a memory-starved run makes
it slower, not faster.

So the honest floor is **8 GB at low thread counts, 10 GB at 8 threads** — and
if you are near the floor, turn threads *down*.

## The floor is a slope, not a cliff

The interesting result is not the floor but the approach to it. Peak RSS is
8.12 GB at every cap, yet wall time goes 1,024 → 6,735 → 13,621 s. The cgroup
counts page cache, so once the 8.12 GB working set and the 6.5 GB index can no
longer both be cached, the index is evicted and re-read continuously.

**The time penalty here is pessimistic.** Those re-reads go over virtiofs from
the macOS host, which is slow; on a native Linux box with local disk the same
eviction would cost far less. The 18 GB run is itself 2.7x slower than the same
job natively (1,024 s against 382 s), so the VM adds overhead before any memory
pressure. The memory floor is real; treat the multipliers as an upper bound.

Practically: **budget ~18 GB for full speed; at 8-10 GB it finishes, but only
if you also drop the thread count.**

## What this does and does not license claiming

**Does:** alignment-based human single-cell quantification in **8 GB** (`-p 2`)
or **10 GB** (`-p 8`), against
rustar's measured 25.6-28.2 GB peak on this machine and STAR's documented ~32 GB
for a human index. That is a 2.5-3x difference and it is the difference between
a 16 GB laptop and a dedicated node.

**Does not:** a head-to-head under the same cap. STAR and rustar were **not run
here** — neither has a human index on this machine, and building one needs more
memory and disk than the VM has. This is our measured floor against their
measured and documented requirements, which is weaker than the plan intended.

**And qualify the 8 GB claim wherever it appears.** The plan's "HISAT2-solo
completes at 8 GB; the others fail to load an index" holds only at low thread
counts — with `-p 8` it is OOM-killed at 8 GB. The claim needs the thread caveat
attached, and the 5.18 GB number replaced with ~8 GB for human.

## Caveats

- One dataset, one index. A larger annotation or `--solo-allelic` (measured at
  12.11 GB natively, against 10.0 GB without) pushes the floor higher still.
- The thread lever was tested only at 8 GB. Whether `-p 1` drops the floor
  below 8 GB is untested; the index alone is 6.5 GB, so there is not much room.
- The VM's virtiofs mount inflates every wall-clock number here, most severely
  the memory-constrained ones.
