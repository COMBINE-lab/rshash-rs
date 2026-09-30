#!/usr/bin/env bash
# Parity harness: build the same index with the C++ reference and the Rust
# port, run streaming lookup of the same reads, and compare per-record hit
# counts (byte-identical TSV) plus the totals.
#
# usage: parity/parity.sh [text.fa] [reads.fq]
# Requires the C++ build (see parity/cpp/CMakeLists.txt) in parity/work/cpp-build
# and the data from parity/prepare_data.sh.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
work="$here/work"
text="${1:-$work/data/ecoli_chunks.fa}"
reads="${2:-$work/data/reads.fq}"
R="$here/../target/release/rshash"
C="$work/cpp-build/rshash"
D="$work/cpp-build/rshash_dump"
T="$here/../target/release/rshash-tools"
mkdir -p "$work/runs"
if command -v module >/dev/null 2>&1 || [ -f /usr/share/Modules/init/bash ]; then
    source /usr/share/Modules/init/bash && module load gcc/14.2.0 >/dev/null 2>&1 || true
fi

# Configurations where the C++ implementation is expected to be correct:
# streaming lookup with the default hash-table last level, contiguous k-mers
# or a single palindromic shape whose run equals the kernel. Levels > 1 need
# frequent minimisers (C++ crashes on an empty level), hence low thresholds.
configs=(
  "k31_l1|-k 31 --level 1"
  "k31_l1_t8|-k 31 --level 1 --t1 8"
  "k31_l2_t4|-k 31 --level 2 --m1 16 --m2 19 --t1 4 --t2 8"
  "k31_l3_t2|-k 31 --level 3 --m1 15 --m2 17 --m3 20 --t1 2 --t2 2 --t3 4"
  "k25_l2_t3|-k 25 --level 2 --m1 13 --m2 15 --t1 3 --t2 3"
  "k32_l1|-k 32 --level 1 --m1 20"
  "k21_l1_loc|-k 21 --level 1 --m1 11 --t1 16 --loc"
  "shape5FFFFFFD_l1|--shapes 1610612733 --level 1 --m1 16"
  "shape5FFFFFFD_l2_t4|--shapes 1610612733 --level 2 --m1 15 --m2 18 --t1 4 --t2 8"
  # expected to differ (DEVIATIONS.md D8/D9): palindromic shape whose kernel
  # run is not centred, and a non-palindromic shape
  "xfail_shape_offcentre|--shapes 15855 --level 1 --m1 3 --t1 64"
  "xfail_shape_nonpal|--shapes 268435445 --level 1 --m1 16"
)

status=0
printf "%-22s %14s %14s %14s %8s %12s %12s\n" config oracle cpp_hits rust_hits records cpp_ext rust_ext
for c in "${configs[@]}"; do
  name="${c%%|*}"; opts="${c#*|}"
  run="$work/runs/$name"
  # shellcheck disable=SC2086
  "$C" build -i "$text" -d "$run.cpp.idx" $opts > "$run.cpp.build.log" 2>&1 || { echo "$name: C++ build failed"; status=1; continue; }
  # shellcheck disable=SC2086
  "$R" build -i "$text" -d "$run.rs.idx" $opts > "$run.rs.build.log" 2>&1
  "$D" "$run.cpp.idx" "$reads" 2> "$run.cpp.err" | grep -v "^loaded" > "$run.cpp.tsv"
  "$R" lookup -d "$run.rs.idx" -q "$reads" --dump "$run.rs.tsv" > "$run.rs.log"
  ch=$(awk '{s+=$4} END {print s}' "$run.cpp.tsv")
  rh=$(awk '{s+=$4} END {print s}' "$run.rs.tsv")
  ce=$(sed -n 's/.*extensions \([0-9]*\)/\1/p' "$run.cpp.err")
  re=$(sed -n 's/^extensions = //p' "$run.rs.log")
  diffs=$(paste "$run.cpp.tsv" "$run.rs.tsv" | awk '$4 != $8' | wc -l)
  # brute-force oracle with the same k-mer / shape semantics
  if [[ "$opts" == *--shapes* ]]; then
    shp=$(sed -E 's/.*--shapes ([^ ]+).*/\1/' <<< "$opts")
    oracle=$("$T" hashtable -i "$text" -q "$reads" --shapes "$shp" --aligned | sed -n 's/^num_positive_kmers = \([0-9]*\).*/\1/p')
  else
    kk=$(sed -E 's/.*-k ([0-9]+).*/\1/' <<< "$opts")
    oracle=$("$T" hashtable -i "$text" -q "$reads" -k "$kk" | sed -n 's/^num_positive_kmers = \([0-9]*\).*/\1/p')
  fi
  [ "$oracle" = "$rh" ] || { echo "$name: Rust ($rh) != oracle ($oracle)"; status=1; }
  printf "%-22s %14s %14s %14s %8s %12s %12s" "$name" "$oracle" "$ch" "$rh" "$diffs" "$ce" "$re"
  if [[ "$name" == xfail_* ]]; then
    if [ "$diffs" -ne 0 ]; then echo "  differs (expected)"; else echo "  ok (expected to differ)"; fi
  elif [ "$diffs" -ne 0 ]; then echo "  MISMATCH"; status=1; else echo "  ok"; fi
done
exit $status
