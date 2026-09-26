#!/bin/sh
# Exactly the commands that produced results/linux_cloud/ and results/params_chain_signed.* for the
# paper, in this order, on one otherwise idle machine. The parameter chain needs SageMath and the
# lattice estimator (commit 53da598) on PYTHONPATH; it runs in parallel with the ballot experiment,
# after the timing benchmarks have finished.
set -e
make all
O=results/linux_cloud
mkdir -p $O/replications/seed2 $O/replications/seed3
rm -f $O/bench.csv
{ uname -a; lscpu | head -20; cc --version | head -1; free -g | head -2; } > $O/machine.txt 2>&1
./build/test_core > $O/test_core.txt
./build/test_protocol > $O/test_protocol.txt
./build/bench_protocol -reps 30 -skip-agg -out $O/bench.csv > $O/bench_stdout.txt
./build/bench_protocol -only-agg -nv 1000 -out $O/bench.csv >> $O/bench_stdout.txt
./build/bench_protocol -only-agg -nv 10000 -out $O/bench.csv >> $O/bench_stdout.txt
( python3 tools/chain.py results/params_chain_signed.jsonl signed 64 > results/params_chain_signed.log 2>&1 ) &
./build/exp_ballot_rej 1000 20260924 $O > $O/ballot_rej_stdout.txt
./build/exp_ballot_rej 500 2 $O/replications/seed2 > $O/replications/seed2/stdout.txt
./build/exp_ballot_rej 500 3 $O/replications/seed3 > $O/replications/seed3/stdout.txt
python3 scripts/summarize.py $O > /dev/null
# EVOLVE's ballot at d = 6, q = 2^36 - 303 (conclusion of the paper)
./build/bench_protocol -reps 10 -skip-agg -d 6 -q 68719476433 -out $O/bench_evolve_d6.csv > /dev/null
wait
# Comparison of two parameter chains (maximum change of every bound, in bits), e.g.
#   git show 31d0dff:results/params_chain_signed.jsonl > old.jsonl
#   python3 scripts/compare_chain.py old.jsonl results/params_chain_signed.jsonl
