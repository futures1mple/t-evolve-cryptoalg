#!/bin/sh
# Same run as scripts/run_windows.bat, for Linux/macOS.
set -e
make all
OUT=results/$(uname -s | tr A-Z a-z)/$(date +%Y%m%d_%H%M)
mkdir -p "$OUT"
{ uname -a; (lscpu 2>/dev/null || sysctl -n machdep.cpu.brand_string); cc --version; } > "$OUT/machine.txt" 2>&1
./build/test_core > "$OUT/test_core.txt"
./build/test_protocol -q > "$OUT/test_protocol.txt"
./build/exp_ballot_rej 1000 20260924 "$OUT" > "$OUT/ballot_rej_stdout.txt"
./build/bench_protocol -reps 50 -skip-agg -out "$OUT/bench.csv" > "$OUT/bench_stdout.txt"
./build/bench_protocol -only-agg -nv 1000 -out "$OUT/bench.csv" >> "$OUT/bench_stdout.txt"
./build/bench_protocol -only-agg -nv 10000 -out "$OUT/bench.csv" >> "$OUT/bench_stdout.txt"
echo "results in $OUT"
