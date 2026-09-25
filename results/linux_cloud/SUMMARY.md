# Measurements

Machine: Intel(R) Xeon(R) Processor @ 2.80GHz; compiler 13.3.0; flags `cc -O3 -march=native`; commit 3785b9d.

All times are medians over the repetitions (single thread).

## Ballots

| Scheme | Vote (ms) | Verify (ms) | Check_k (ms) | Attempts | Size (KiB) |
|---|---|---|---|---|---|
| EVOLVE n=4 L=1 d=7 log2q=42.0 | 14.5 | 2.9 | nan | 2.93 | 58.1 |
| T-EVOLVE n=4 t=3 L=1 w=-1 d=7 log2q=42.0 | 36.0 | 10.2 | 1.40 | 2.93 | 107.7 |
| T-EVOLVE n=4 t=3 L=2 w=1 d=7 log2q=42.0 | 51.3 | 14.4 | 1.63 | 3.03 | 135.7 |
| T-EVOLVE n=5 t=3 L=1 w=-1 d=7 log2q=42.0 | 58.6 | 12.0 | 1.40 | 3.07 | 126.3 |

## Aggregation (per authority)

| N_V | E | Round 1 (s) | Round 2 (s) | VerAgg (s) | Published per ballot (KiB) | Status |
|---|---|---|---|---|---|---|
| 1000 | 6 | 29.3 | 10.9 | 5.4 | 46.7 | measured |
| 10000 | 24 | 159.6 | 60.6 | 32.4 | 27.9 | measured |
| 100000 | 217 | 1443 | 730 | 293 | — | extrapolated |
| 1000000 | 2142 | 14245 | 7210 | 2895 | — | extrapolated |

Rounds (kappa = 4, attempt success 1/3): an authority succeeds in a round with probability 0.80; expected rounds 1.22 with 4 honest authorities, at most 1.94 with 3. Responses of rounds that fail are not published (announce-then-reveal): a failed round costs the computation of both rounds but no publication beyond round 1.

## Ballot rejection experiment

| n | L | log2 sigma_J | Ballots | Attempts (mean ± 95%) | max ‖s‖/T | Verified | GOF p | Independence p |
|---|---|---|---|---|---|---|---|---|
| 4 | 1 | 14.62 | 1000 | 3.09 ± 0.17 | 0.53 | 1000 | 0.653 | 0.105 |
| 4 | 2 | 14.85 | 1000 | 2.96 ± 0.15 | 0.50 | 1000 | 0.166 | 0.217 |
| 4 | 5 | 15.30 | 1000 | 3.07 ± 0.17 | 0.48 | 1000 | 0.585 | 0.028 |
| 4 | 10 | 15.74 | 1000 | 3.00 ± 0.15 | 0.48 | 1000 | 0.866 | 0.469 |
| 5 | 1 | 14.67 | 1000 | 3.02 ± 0.15 | 0.55 | 1000 | 0.952 | 0.616 |

Expected attempts: M = 2.99 for every L.

Over 15 tests of each kind (3 runs with different seeds), Fisher's combined p-value is 0.70 for goodness of fit and 0.15 for witness independence; the smallest single p-values are 0.020 and 0.020 (Bonferroni threshold 0.05/15 = 0.0033).
