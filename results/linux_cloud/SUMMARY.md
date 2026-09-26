# Measurements

Machine: Intel(R) Xeon(R) Processor @ 2.10GHz; compiler 13.3.0; flags `cc -O3 -march=native`; commit d03fd81.

All times are medians over the repetitions (single thread).

## Ballots

| Scheme | Vote (ms) | Verify (ms) | Check_k (ms) | Attempts | Size (KiB) |
|---|---|---|---|---|---|
| EVOLVE n=4 L=1 d=7 log2q=42.0 | 12.5 | 2.9 | nan | 2.47 | 58.6 |
| T-EVOLVE n=4 t=3 L=1 w=-1 d=7 log2q=42.0 | 28.7 | 9.3 | 1.19 | 2.07 | 108.9 |
| T-EVOLVE n=4 t=3 L=2 w=1 d=7 log2q=42.0 | 39.8 | 12.9 | 1.47 | 2.40 | 138.0 |
| T-EVOLVE n=5 t=3 L=1 w=-1 d=7 log2q=42.0 | 47.2 | 10.8 | 1.26 | 3.30 | 128.4 |

## Aggregation (per authority)

| N_V | E | Round 1 (s) | Round 2 (s) | VerAgg (s) | Published per ballot (KiB) | Status |
|---|---|---|---|---|---|---|
| 1000 | 6 | 25.3 | 7.8 | 4.9 | 47.4 | measured |
| 10000 | 24 | 136.3 | 36.4 | 29.9 | 28.2 | measured |
| 100000 | 217 | 1232 | 658 | 270 | — | extrapolated |
| 1000000 | 2142 | 12162 | 6491 | 2667 | — | extrapolated |

Rounds (kappa = 4, attempt success 1/M^2 = 0.333): an authority succeeds in a round with probability 0.80; expected rounds 1.22 with 4 honest authorities, at most 1.94 with 3. If fewer than t authorities announce an attempt, no responses are published. Misbehaviour cannot make a round fail that would otherwise succeed; an authority that announces and then posts no valid contribution only makes the honest authorities publish their responses in a round that fails anyway, and is excluded afterwards (at most n - t such rounds). The benchmark has no faulty authorities.

## Ballot rejection experiment

| n | L | log2 sigma_J | Ballots | Attempts (mean ± 95%) | max ‖s‖/T | Verified | GOF p | Independence p |
|---|---|---|---|---|---|---|---|---|
| 4 | 1 | 14.98 | 1000 | 2.87 ± 0.15 | 0.52 | 1000 | 0.577 | 0.892 |
| 4 | 2 | 15.21 | 1000 | 3.11 ± 0.16 | 0.50 | 1000 | 0.885 | 0.009 |
| 4 | 5 | 15.66 | 1000 | 2.99 ± 0.15 | 0.48 | 1000 | 0.558 | 0.228 |
| 4 | 10 | 16.11 | 1000 | 2.93 ± 0.15 | 0.48 | 1000 | 0.943 | 0.270 |
| 5 | 1 | 15.04 | 1000 | 2.92 ± 0.14 | 0.55 | 1000 | 0.866 | 0.337 |

Expected attempts: M = exp(43669/39762) = 3.00 for every L.

Over 15 tests of each kind (3 runs with different seeds), Fisher's combined p-value is 0.99 for goodness of fit and 0.09 for witness independence; the smallest single p-values are 0.200 and 0.009 (Bonferroni threshold 0.05/15 = 0.0033).

## Parameter chain

`results/params_chain_signed.jsonl` (every estimator result, all attacks) and `results/params_chain_signed.log` were produced by `tools/chain.py` with `build/param_report` of commit d03fd81 on the same machine (lattice estimator 53da598 in SageMath). All tests (`test_core.txt`, `test_protocol.txt`) passed on this build.
