#!/usr/bin/env bash
# Timing benchmark: single-thread vs OpenMP, adaptive vs fixed step, two resolutions.
# Every run appends one row to results/timing.csv (written by the renderer's --csv).
# Usage: scripts/bench.sh      (run from the project root, after building)
set -euo pipefail
cd "$(dirname "$0")/.."

BIN=./build/blackhole
CSV=results/timing.csv
REPEATS=${REPEATS:-3}
TMP=$(mktemp -d)
mkdir -p results renders
rm -f "$CSV"

# Checker sky: sky lookups are trivial either way, and it avoids rebuilding the
# procedural star texture each run (which isn't timed anyway). 1 spp = 1 ray per pixel.
COMMON="--sky checker --spp 1 --pos 0,-20,2"

for res in 640x360 1280x720; do
  for mode in "--step 0.02" "--fixed 0.03"; do
    for threads in 1 4 8; do
      for ((rep = 1; rep <= REPEATS; rep++)); do
        echo ">> $res  $mode  threads=$threads  rep=$rep"
        $BIN $COMMON --res "$res" $mode --threads "$threads" --csv "$CSV" -o "$TMP/out.png" | tail -1
      done
    done
  done
done

# Heatmaps for the report, with the same color range so they can be compared directly.
$BIN $COMMON --res 1280x720 --step 0.02  --heat-range 100,20000 -o "$TMP/a.png" --heatmap renders/heatmap_adaptive.png
$BIN $COMMON --res 1280x720 --fixed 0.03 --heat-range 100,20000 -o "$TMP/f.png" --heatmap renders/heatmap_fixed.png

# Accuracy vs. cost: compare each setting against a very fine adaptive reference
# (step 0.002) by counting pixels whose color differs. 640x360, all threads.
ACC=results/accuracy
mkdir -p "$ACC"
rm -f results/accuracy.csv
$BIN $COMMON --res 640x360 --step 0.002 -o "$ACC/ref.png" --csv results/accuracy.csv
for s in 0.01 0.02 0.05 0.1; do
  $BIN $COMMON --res 640x360 --step $s  -o "$ACC/adaptive_$s.png" --csv results/accuracy.csv
done
for s in 0.03 0.1 0.3 1; do
  $BIN $COMMON --res 640x360 --fixed $s -o "$ACC/fixed_$s.png"    --csv results/accuracy.csv
done

python3 scripts/summarize_bench.py
rm -rf "$TMP"
