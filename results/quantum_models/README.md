# Hardness under four cost models

All values are log2 of the estimated cost of the best attack (lattice estimator 53da598, script `tools/quantum_models.py`), M-LWE (hiding) / M-SIS (binding) for the T-EVOLVE sets. `MATZOV`: classical, the model of the parameter chain; `MATZOV_q`: the same with the quantum nearest-neighbour cost of AGPS20 (depth x width); `ADPS16`, `ADPS16_q`: core-SVP, 2^(0.292 beta) and 2^(0.265 beta). `kyber512` = parameters of Kyber-512 / ML-KEM-512 (n = 512, q = 3329, centered binomial eta = 3, m = 512), for comparison; its specification (version 3.02, Table 4) gives core-SVP 118 / 107 from the primal attack alone. `1e6_d8`: N_V = 10^6 with d = 8 at the same q = 2^45 - 591 (beta_SIS = 2^44.6 from param_report). `eo_2048_rlwe`, `eo_1024_rlwe`: ciphertext RLWE instance of the NTRU-based key-release encryption of Appendix D (ternary secret and noise, N' + 256 samples, q' = 9123495937 resp. 2256363521); `eo_2048_ntru`: recovery of the NTRU master key (f, g of standard deviation 1.17 sqrt(q'/2N')), estimator NTRU module, circulant type: uSVP and the dense-sublattice attack return no finite cost for this instance (not overstretched), the value is the BDD attack.

| set | MATZOV | MATZOV_q | ADPS16 | ADPS16_q |
|---|---|---|---|---|
| 1e4 | 138.4 / 144.5 | 133.3 / 138.9 | 110.7 / 116.8 | 100.7 / 106.0 |
| 1e5 | 135.0 / 138.2 | 130.3 / 133.2 | 107.2 / 110.1 | 97.5 / 99.9 |
| 1e6 | 128.9 / 134.6 | 124.7 / 129.9 | 101.0 / 106.3 | 91.7 / 96.5 |
| 1e6_d8 | 148.5 / 155.0 | 142.5 / 148.4 | 121.1 / 127.6 | 110.2 / 115.8 |
| kyber512 | 139.7 | 134.4 | 115.5 | 106.3 |
| eo_2048_rlwe | 199.0 | 188.3 | 174.0 | 158.5 |
| eo_2048_ntru | 621.5 | 566.8 | 625.2 | 567.4 |
| eo_1024_rlwe | 103.0 | - | 75.0 | - |

Produced 2026-10-16 (d = 8 and EO rows 2026-10-18); each JSON file lists every attack with its block size beta.
