#!/usr/bin/env python3
"""Pull the R1/R2 pairs whose reads touch the mini-reference loci."""
import gzip, hashlib, sys
import numpy as np

D = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/human"
O = "/tmp/claude-501/-Users-iandriver-Downloads-hisat2/2400b599-e2d4-4d43-95ad-d6f519cf1f37/scratchpad/filt"
F = f"{D}/pbmc_1k_v3_fastqs"


def h(name):
    return int.from_bytes(hashlib.blake2b(name, digest_size=8).digest(), "little")


print("hashing names...", flush=True)
keys = np.fromiter((h(l.rstrip(b"\n")) for l in open(f"{O}/names.txt", "rb")),
                   dtype=np.uint64)
keys.sort()
print(f"  {len(keys):,} names", flush=True)

out1 = gzip.open(f"{O}/sub_R1.fq.gz", "wb", compresslevel=1)
out2 = gzip.open(f"{O}/sub_R2.fq.gz", "wb", compresslevel=1)
kept = total = 0
for lane in ("L001", "L002"):
    f1 = gzip.open(f"{F}/pbmc_1k_v3_S1_{lane}_R1_001.fastq.gz", "rb")
    f2 = gzip.open(f"{F}/pbmc_1k_v3_S1_{lane}_R2_001.fastq.gz", "rb")
    while True:
        a = f1.readline()
        if not a:
            break
        rec1 = (a, f1.readline(), f1.readline(), f1.readline())
        rec2 = (f2.readline(), f2.readline(), f2.readline(), f2.readline())
        total += 1
        name = a[1:].split()[0]
        k = h(name)
        i = np.searchsorted(keys, k)
        if i < len(keys) and keys[i] == k:
            out1.write(b"".join(rec1))
            out2.write(b"".join(rec2))
            kept += 1
        if total % 10_000_000 == 0:
            print(f"  {total:,} read, {kept:,} kept", flush=True)
    f1.close(); f2.close()
out1.close(); out2.close()
print(f"kept {kept:,} of {total:,} pairs")
