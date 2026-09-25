# Measurements

Machine: Intel(R) Xeon(R) Processor @ 2.80GHz; compiler 13.3.0; flags `cc -O3 -march=native`; commit 30a7149 as reported by the binaries; they were built from the working tree that was then committed unchanged as 9f3a7bd (domain-separated seed expansion), so the code measured is that of 9f3a7bd.

All times are medians over the repetitions (single thread).

## Ballots

| Scheme | Vote (ms) | Verify (ms) | Check_k (ms) | Attempts | Size (KiB) |
|---|---|---|---|---|---|
| EVOLVE n=4 L=1 d=7 log2q=42.0 | 16.9 | 3.0 | nan | 3.13 | 59.1 |
| T-EVOLVE n=4 t=3 L=1 w=-1 d=7 log2q=42.0 | 50.8 | 10.5 | 1.46 | 3.43 | 109.0 |
| T-EVOLVE n=4 t=3 L=2 w=1 d=7 log2q=42.0 | 49.4 | 14.4 | 1.68 | 3.03 | 137.1 |
| T-EVOLVE n=5 t=3 L=1 w=-1 d=7 log2q=42.0 | 42.2 | 11.8 | 1.44 | 2.60 | 127.8 |

## Aggregation (per authority)

| N_V | E | Round 1 (s) | Round 2 (s) | VerAgg (s) | Published per ballot (KiB) | Status |
|---|---|---|---|---|---|---|
| 1000 | 6 | 28.8 | 10.9 | 5.6 | 46.7 | measured |
| 10000 | 24 | 157.2 | 62.6 | 33.5 | 27.9 | measured |
| 100000 | 217 | 1422 | 755 | 303 | — | extrapolated |
| 1000000 | 2142 | 14033 | 7453 | 2986 | — | extrapolated |

Rounds (kappa = 4, attempt success 1/3): an authority succeeds in a round with probability 0.80; expected rounds 1.22 with 4 honest authorities, at most 1.94 with 3. Responses of rounds that fail are not published (announce-then-reveal): a failed round costs the computation of both rounds but no publication beyond round 1.

## Ballot rejection experiment

| n | L | log2 sigma_J | Ballots | Attempts (mean ± 95%) | max ‖s‖/T | Verified | GOF p | Independence p |
|---|---|---|---|---|---|---|---|---|
| 4 | 1 | 14.62 | 1000 | 2.94 ± 0.15 | 0.53 | 1000 | 0.301 | 0.753 |
| 4 | 2 | 14.85 | 1000 | 2.94 ± 0.15 | 0.50 | 1000 | 0.757 | 0.562 |
| 4 | 5 | 15.30 | 1000 | 2.95 ± 0.15 | 0.48 | 1000 | 0.086 | 0.649 |
| 4 | 10 | 15.74 | 1000 | 3.10 ± 0.16 | 0.48 | 1000 | 0.096 | 0.270 |
| 5 | 1 | 14.67 | 1000 | 3.03 ± 0.16 | 0.55 | 1000 | 0.998 | 0.614 |

Expected attempts: M = 2.99 for every L.

Over 15 tests of each kind (3 runs with different seeds), Fisher's combined p-value is 0.48 for goodness of fit and 0.46 for witness independence; the smallest single p-values are 0.066 and 0.135 (Bonferroni threshold 0.05/15 = 0.0033).
