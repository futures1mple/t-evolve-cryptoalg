# Measurements

Machine: Intel(R) Xeon(R) Processor @ 2.10GHz; compiler 13.3.0; flags `cc -O3 -march=native`; commit d8a88fb.

All times are medians over the repetitions (single thread).

## Ballots

| Scheme | Vote (ms) | Verify (ms) | Check_k (ms) | Attempts | Size (KiB) |
|---|---|---|---|---|---|
| EVOLVE n=4 L=1 d=7 log2q=42.0 | 13.0 | 2.7 | nan | 3.33 | 59.1 |
| T-EVOLVE n=4 t=3 L=1 w=-1 d=7 log2q=42.0 | 46.0 | 9.5 | 1.30 | 3.70 | 109.0 |
| T-EVOLVE n=4 t=3 L=2 w=1 d=7 log2q=42.0 | 62.5 | 14.3 | 1.61 | 3.27 | 137.1 |
| T-EVOLVE n=5 t=3 L=1 w=-1 d=7 log2q=42.0 | 52.7 | 11.1 | 1.34 | 2.80 | 127.8 |

## Aggregation (per authority)

| N_V | E | Round 1 (s) | Round 2 (s) | VerAgg (s) | Published per ballot (KiB) | Status |
|---|---|---|---|---|---|---|
| 1000 | 6 | 27.1 | 11.4 | 5.1 | 46.7 | measured |
| 10000 | 24 | 148.0 | 66.9 | 31.3 | 27.9 | measured |
| 100000 | 217 | 1338 | 806 | 283 | — | extrapolated |
| 1000000 | 2142 | 13205 | 7957 | 2791 | — | extrapolated |

Rounds (kappa = 4, attempt success 1/3): an authority succeeds in a round with probability 0.80; expected rounds 1.22 with 4 honest authorities, at most 1.94 with 3. Responses of rounds that fail are not published (announce-then-reveal): a failed round costs the computation of both rounds but no publication beyond round 1.

## Ballot rejection experiment

| n | L | log2 sigma_J | Ballots | Attempts (mean ± 95%) | max ‖s‖/T | Verified | GOF p | Independence p |
|---|---|---|---|---|---|---|---|---|
| 4 | 1 | 14.62 | 1000 | 2.96 ± 0.15 | 0.53 | 1000 | 0.455 | 0.834 |
| 4 | 2 | 14.85 | 1000 | 2.92 ± 0.15 | 0.50 | 1000 | 0.788 | 0.677 |
| 4 | 5 | 15.30 | 1000 | 2.95 ± 0.14 | 0.48 | 1000 | 0.591 | 0.634 |
| 4 | 10 | 15.74 | 1000 | 3.15 ± 0.16 | 0.48 | 1000 | 0.334 | 0.286 |
| 5 | 1 | 14.67 | 1000 | 3.03 ± 0.16 | 0.55 | 1000 | 0.978 | 0.211 |

Expected attempts: M = 2.99 for every L.

Over 15 tests of each kind (3 runs with different seeds), Fisher's combined p-value is 0.90 for goodness of fit and 0.37 for witness independence; the smallest single p-values are 0.126 and 0.062 (Bonferroni threshold 0.05/15 = 0.0033).
