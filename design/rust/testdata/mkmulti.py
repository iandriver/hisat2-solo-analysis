#!/usr/bin/env python3
"""Generate multi.fa -- a deliberately awkward reference for the emitter test.

The example reference (`22_20-21M.fa`) is a single sequence with one internal N
run, which leaves several rules untested. Every sequence here exists to exercise
one of them:

  seqA  header carries a description   -> refnames store the FULL line, not the
                                          whitespace-truncated `_refnames_nospace`
  seqB  leading Ns, then an internal N gap -> RefRecord.off on a first record
  seqC  trailing Ns                    -> a zero-length RefRecord at sequence end
  seqD  lowercase n, single trailing N -> ambiguity is case-insensitive
  seqE  wholly lowercase               -> bases are upper-cased before coding

Both defects the emitter had after passing on the example reference (refname
truncation, and the missing zero-length trailing records that make `.3.ht2` hold
more records than `.1.ht2` has rstarts) were caught by this file.
"""
import random

random.seed(7)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))

recs = [
    ("seqA some description here", seq(5000)),
    ("seqB", "N" * 137 + seq(3000) + "N" * 61 + seq(2500)),
    ("seqC", seq(1200) + "N" * 400),
    ("seqD", seq(900) + "n" * 3 + seq(1100) + "N"),
    ("seqE_lower", seq(2000).lower()),
]

with open('multi.fa', 'w') as f:
    for name, s in recs:
        f.write(">%s\n" % name)
        for i in range(0, len(s), 60):
            f.write(s[i:i + 60] + "\n")
