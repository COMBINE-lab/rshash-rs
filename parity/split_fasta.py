#!/usr/bin/env python3
"""Split every record of a FASTA file into chunks of random length (unitig-like).

usage: split_fasta.py <in.fa> <out.fa> <min_len> <max_len> <seed>
"""
import random
import sys

path, out, lo, hi, seed = sys.argv[1:6]
lo, hi = int(lo), int(hi)
rng = random.Random(int(seed))
seqs, cur = [], []
for line in open(path):
    if line.startswith(">"):
        if cur:
            seqs.append("".join(cur))
        cur = []
    else:
        cur.append(line.strip())
if cur:
    seqs.append("".join(cur))
i = 0
with open(out, "w") as f:
    for s in seqs:
        p = 0
        while p < len(s):
            n = rng.randint(lo, hi)
            f.write(f">c{i}\n{s[p:p + n]}\n")
            i += 1
            p += n
