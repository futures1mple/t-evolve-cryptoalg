# Measurements

Machine: Intel(R) Xeon(R) Processor @ 2.10GHz; compiler 13.3.0; flags `cc -O3 -march=native`; commit cbc1280.

All times are medians over the repetitions (single thread).

## Ballots

| Scheme | Vote (ms) | Verify (ms) | Check_k (ms) | Attempts | Size (KiB) |
|---|---|---|---|---|---|
| EVOLVE n=4 L=1 d=7 log2q=42.0 | 11.6 | 2.6 | nan | 2.73 | 58.1 |
| T-EVOLVE n=4 t=3 L=1 w=-1 d=7 log2q=42.0 | 30.8 | 9.1 | 1.12 | 2.43 | 107.8 |
| T-EVOLVE n=4 t=3 L=2 w=1 d=7 log2q=42.0 | 41.6 | 12.6 | 1.31 | 2.47 | 135.7 |
| T-EVOLVE n=5 t=3 L=1 w=-1 d=7 log2q=42.0 | 31.4 | 10.5 | 1.13 | 2.87 | 126.3 |

## Aggregation (per authority)

| N_V | E | Round 1 (s) | Round 2 (s) | VerAgg (s) | Published per ballot (KiB) | Status |
|---|---|---|---|---|---|---|
| 1000 | 6 | 24.6 | 6.2 | 4.6 | 46.7 | measured |
| 10000 | 24 | 130.7 | 81.7 | 29.2 | 27.9 | measured |
| 100000 | 217 | 1182 | 739 | 264 | — | extrapolated |
| 1000000 | 2142 | 11663 | 7292 | 2606 | — | extrapolated |

Rounds (kappa = 4, attempt success 1/M^2 = 0.333): an authority succeeds in a round with probability 0.80; expected rounds 1.22 with 4 honest authorities, at most 1.94 with 3. If fewer than t authorities announce an attempt, no responses are published. Misbehaviour cannot make a round fail that would otherwise succeed; an authority that announces and then posts no valid contribution only makes the honest authorities publish their responses in a round that fails anyway, and is excluded afterwards (at most n - t such rounds). The benchmark has no faulty authorities.

## Ballot rejection experiment

| n | L | log2 sigma_J | Ballots | Attempts (mean ± 95%) | max ‖s‖/T | Verified | GOF p | Independence p |
|---|---|---|---|---|---|---|---|---|
| 4 | 1 | 14.62 | 1000 | 2.83 ± 0.14 | 0.53 | 1000 | 0.394 | 0.107 |
| 4 | 2 | 14.85 | 1000 | 2.96 ± 0.15 | 0.50 | 1000 | 0.195 | 0.080 |
| 4 | 5 | 15.30 | 1000 | 3.01 ± 0.15 | 0.48 | 1000 | 0.903 | 0.889 |
| 4 | 10 | 15.74 | 1000 | 3.17 ± 0.15 | 0.48 | 1000 | 0.742 | 0.388 |
| 5 | 1 | 14.67 | 1000 | 2.96 ± 0.15 | 0.55 | 1000 | 0.522 | 0.208 |

Expected attempts: M = 2.99 for every L.

Over 15 tests of each kind (3 runs with different seeds), Fisher's combined p-value is 0.39 for goodness of fit and 0.52 for witness independence; the smallest single p-values are 0.061 and 0.080 (Bonferroni threshold 0.05/15 = 0.0033).

## Parameter chain

`results/params_chain_signed.jsonl` (every estimator result, all attacks) and `results/params_chain_signed.log` were produced by `tools/chain.py` with `build/param_report` of commit cbc1280 on the same machine (lattice estimator 53da598 in SageMath). All tests (`test_core.txt`, `test_protocol.txt`) passed on this build.
