#!/usr/bin/env bash
# Download / generate the parity and benchmark datasets into parity/work/data.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
d="$here/work/data"
mkdir -p "$d"
cd "$d"
# E. coli K-12 MG1655 (4.6 Mbp) from NCBI
[ -s ecoli.fa ] || curl -sS -o ecoli.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text"
# unitig-like chunks (40..3000 bp) to exercise sequence boundaries
[ -s ecoli_chunks.fa ] || python3 "$here/split_fasta.py" ecoli.fa ecoli_chunks.fa 40 3000 1
# 100k 150 bp reads: 1% substitutions, both strands, 20% random reads
[ -s reads.fq ] || python3 "$here/sample_reads.py" ecoli.fa reads.fq 100000 150 0.01 0.2 2
if [ "${1:-}" = "--large" ]; then
  # C. elegans WBcel235 (100 Mbp) from Ensembl + 500k reads
  [ -s celegans.fa ] || { curl -sS -o celegans.fa.gz "https://ftp.ensembl.org/pub/current_fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz" && gunzip -kf celegans.fa.gz; }
  [ -s celegans_reads.fq ] || python3 "$here/sample_reads.py" celegans.fa celegans_reads.fq 500000 150 0.01 0.2 3
fi
ls -la "$d"
