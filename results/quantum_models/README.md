# Hardness under four cost models

All values are log2 of the estimated cost of the best attack (lattice estimator 53da598, script `tools/quantum_models.py`), M-LWE (hiding) / M-SIS (binding). `MATZOV`: classical, the model of the parameter chain; `MATZOV_q`: the same with the quantum nearest-neighbour cost of AGPS20 (depth x width); `ADPS16`, `ADPS16_q`: core-SVP, 2^(0.292 beta) and 2^(0.265 beta). `kyber512` = parameters of Kyber-512 / ML-KEM-512 (n = 512, q = 3329, centered binomial eta = 3, m = 512), for comparison; its specification (version 3.02, Table 4) gives core-SVP 118 / 107 from the primal attack alone. `1e6_d8`: N_V = 10^6 with d = 8 at the same q (M-LWE only).

| set | MATZOV | MATZOV_q | ADPS16 | ADPS16_q |
|---|---|---|---|---|
| 1e4 | 138.4 / 144.5 | 133.3 / 138.9 | 110.7 / 116.8 | 100.7 / 106.0 |
| 1e5 | 135.0 / 138.2 | 130.3 / 133.2 | 107.2 / 110.1 | 97.5 / 99.9 |
| 1e6 | 128.9 / 134.6 | 124.7 / 129.9 | 101.0 / 106.3 | 91.7 / 96.5 |
| 1e6_d8 | 148.5 | 142.5 | - | - |
| kyber512 | 139.7 | 134.4 | 115.5 | 106.3 |

Produced 2026-10-16; each JSON file lists every attack with its block size beta.
