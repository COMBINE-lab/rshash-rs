#!/usr/bin/env python3
"""Sample reads from a FASTA text: substrings (either strand) with substitution
errors, plus fully random reads. Deterministic for a given seed.

usage: sample_reads.py <text.fa> <out.fq> <n> <len> <error_rate> <random_fraction> <seed>
"""
import random
import sys

path, out, n, length, err, rnd_frac, seed = sys.argv[1:8]
n, length, err, rnd_frac = int(n), int(length), float(err), float(rnd_frac)
rng = random.Random(int(seed))
seqs, cur = [], []
for line in open(path):
    if line.startswith(">"):
        if cur:
            seqs.append("".join(cur))
        cur = []
    else:
        cur.append(line.strip().upper())
if cur:
    seqs.append("".join(cur))
seqs = [s for s in seqs if len(s) >= length]
weights = [len(s) for s in seqs]
comp = str.maketrans("ACGT", "TGCA")
with open(out, "w") as f:
    for i in range(n):
        if rng.random() < rnd_frac:
            r = "".join(rng.choice("ACGT") for _ in range(length))
        else:
            s = rng.choices(seqs, weights)[0]
            a = rng.randrange(len(s) - length + 1)
            r = s[a:a + length]
            if rng.random() < 0.5:
                r = r.translate(comp)[::-1]
            r = "".join(rng.choice("ACGT".replace(c, "")) if c in "ACGT" and rng.random() < err else c for c in r)
        f.write(f"@r{i}\n{r}\n+\n{'I' * len(r)}\n")
