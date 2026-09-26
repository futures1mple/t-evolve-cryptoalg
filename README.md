# T-EVOLVE — reference implementation

Research implementation of **T-EVOLVE**, a threshold (t-out-of-n) variant of the lattice-based
e-voting protocol EVOLVE, and of EVOLVE's ballot on the same code base, for the paper

> J. Alizadeh, *Threshold Post-Quantum E-Voting from Lattice Commitments* (preprint, 2026).

It implements the protocol, runs the experiments and benchmarks reported in the paper, and
computes the parameters from the actual witnesses of the proofs (`docs/PARAMETERS.md`).

Repository: https://github.com/futures1mple/t-evolve-cryptoalg

    git clone https://github.com/futures1mple/t-evolve-cryptoalg.git

**This is research code.** It is not constant-time, it has not been audited, and the public-key
encryption of the seeds is a placeholder (see *Scope*). Do not use it to run elections.

## What is implemented

| Part | Files | Paper |
|---|---|---|
| SHAKE128/256 (FIPS 202) | `src/keccak.c` | — |
| `R_q = Z_q[X]/(X^256+1)`, `q = 17 mod 32` prime `< 2^50`, 3-level partial NTT, Montgomery arithmetic | `src/ring.c` | §7, choice of `q` |
| Discrete Gaussian samplers with bounded error (256-bit tables, convolution), exact rejection (up to a capped counter, deviation below 2^-(2^23) per call) with overflow-checked integer arithmetic, challenges, EVOLVE's permutation challenges | `src/sample.c` | §8 (samplers, rejection), §7 |
| BDLOP commitments (multi-message form), Shamir sharing, Reed–Solomon parity checks, Lagrange | `src/bdlop.c`, `src/shamir.c` | §3, §5 |
| Ballot with the EVOLVE-style ballot proof, one joint rejection step, OR-proofs, weight proof; seed-derived share randomness | `src/ballot.c` | §5 (Fig. 2), §7 |
| Authority aggregation: recommitment tree, `Pi_open` and `Pi_zero` after Baum–Lyubashevsky 2017 with one joint rejection per stage over all blocks, `kappa` attempts committed by `G(W, salt)`, announce-then-reveal | `src/agg.c` | §7 |
| Combine / Verify (interpolation at 0) | `src/tally.c` | §5 |
| EVOLVE ballot (additive shares, OR-proof on the sum) on the same core | `src/evolve.c` | §8, §9 (comparison) |
| Encoding (packed `Z_q`, Golomb–Rice for Gaussian vectors) — sizes are those of real encodings | `src/codec.c` | Tables 1, 3 |

## Build

Linux / macOS (gcc or clang):

    make            # builds everything into build/
    make test       # core and protocol tests (about 5 minutes; the aggregation test is the long part)

Windows: install [MSYS2](https://www.msys2.org), then in the *MSYS2 UCRT64* shell
`pacman -S mingw-w64-ucrt-x86_64-gcc git`, add `C:\msys64\ucrt64\bin` to `PATH`, and from the
repository root in `cmd`:

    scripts\build_windows.bat

The code needs a 64-bit compiler with `__int128` (gcc, clang, MinGW-w64). MSVC is not supported.

## Reproducing the paper's measurements

Windows (about one hour; plug in the power and select the best-performance power mode):

    scripts\run_windows.bat

Linux / macOS:

    scripts/run_linux.sh

Both write a dated folder under `results/` with the machine description, test logs, the ballot
rejection experiment and `bench.csv`. Every CSV row records CPU, compiler, flags, commit and seed.

The exact commands behind the numbers of the paper (`results/linux_cloud/` and
`results/params_chain_signed.*`) are in `scripts/reproduce_paper.sh`. `scripts/compare_chain.py`
compares two runs of the parameter chain and prints the largest change of every bound in bits,
e.g. between the versions of `results/params_chain_signed.jsonl` in two commits.

Individual programs:

| Program | What it does |
|---|---|
| `build/test_core` | SHAKE known-answer tests, NTT and products against a schoolbook reference, challenge permutations, chi-square tests of the Gaussian sampler |
| `build/test_protocol [-q]` | ballots (completeness, encoding round trip, 9 kinds of tampering rejected), share checks, interpolation, EVOLVE ballots, a complete aggregation and tally with tampering tests |
| `build/exp_ballot_rej [ballots] [seed] [outdir]` | experiment on the ballot proof: independent ballots, `sigma_J` fixed in advance, full proof and verification after encoding; attempts, `‖shift‖/T`, goodness-of-fit of the accepted responses to `D_{sigma_J}`, and a two-sample test of witness independence (real vs. simulated OR branch) |
| `build/bench_protocol [-nv N] [-reps R] [-d d -q q] [-out f.csv] [-skip-agg / -only-agg]` | micro benchmarks, ballots of T-EVOLVE (yes/no, 1-of-2, n=5) and EVOLVE, and a complete tally of `N` ballots by all authorities, with repeated rounds until `t` succeed |
| `build/param_report NV n t L w d q mode logQ` | every quantity of the parameter chain as JSON |
| `tools/chain.py` | search of `(d, q)` with the lattice estimator (needs SageMath and the estimator, commit 53da598) |

## Scope and limitations

* **PKE.** The paper needs labeled IND-CCA encryption with verifiable decryption. Here the seed
  of each share is encrypted with a placeholder (`seed xor SHAKE256(key || id)`), which fixes the
  data flow and the 32-byte payload but is not secure; the cost of a real KEM (e.g. ML-KEM, tens of
  microseconds) is not included, as in the paper's tables. Signatures and the proofs of correct
  decryption that authorities attach to complaints (`ProveDec`) are not implemented either. The
  measured sizes therefore count each ciphertext as its 32-byte payload and contain no signature;
  the paper gives an estimate with ML-KEM-768 and ML-DSA-44 separately.
* **Side channels.** Nothing is constant-time (in particular the table lookups of the samplers).
  The samplers use integer arithmetic only and are platform-independent; their statistical distance
  from the ideal distributions is bounded (`docs/PARAMETERS.md`, §5).
* **Scale.** Aggregation is measured for `N_V = 10^3` and `10^4`. For `10^5` and `10^6` the cost is
  linear in the number of blocks `E`, and the benchmark reports per-block costs from which the paper
  extrapolates; these numbers are marked as extrapolated.
* **Single thread.** All measurements are single-threaded; the blocks of the aggregation are
  independent and would parallelize.

## Layout

    include/, src/     library
    tests/             tests
    experiments/       ballot rejection experiment
    bench/             benchmarks
    tools/             parameter chain (C report + Python/Sage search)
    scripts/           build and run scripts (Windows, Linux)
    docs/              PARAMETERS.md: derivation of every parameter
    results/           measurements (one folder per machine and date)

## Parameters

Chosen by `tools/chain.py` (signed second challenge, seed-derived share randomness, grinding bound
`2^64`, lattice estimator 53da598; every estimate with all attacks is in
`results/params_chain_signed.jsonl`). `N = 256`, `sigma = 1`, `n = 4`, `t = 3`, fan-in 30, blocks of
`k = 500`, `ell = 517`, `kappa = 4`.

| N_V | d | q | beta_SIS | M-LWE (hiding) | M-SIS (binding) | Ballot, yes/no | Ballot, 1 of 2 | Authority, per ballot |
|---|---|---|---|---|---|---|---|---|
| 10^4 | 7 | 2^42 - 143 | 2^41.6 | 2^138.4 | 2^144.5 | 108.9 KiB | 138.0 KiB | 28.2 KiB |
| 10^5 | 7 | 2^43 - 175 | 2^42.9 | 2^135.0 | 2^138.2 | 110.1 KiB | 139.4 KiB | 28.3 KiB |
| 10^6 | 7 | 2^45 - 591 | 2^44.5 | 2^128.9 | 2^134.6 | 112.6 KiB | 142.3 KiB | 29.7 KiB |

`q` is the largest prime `q = 17 (mod 32)` below `2^b`, for the smallest `b` with `q > beta_SIS` and both
estimates at least `2^128`, so that elements of `Z_q` take exactly `b` bits. All entries are computed, not
measured; the estimates are rounded to one decimal and are not lower bounds (for one of two candidates the
binding estimates are `2^144.0`, `2^137.6` and `2^134.0`). They concern the primitives, not the security of
the protocol.

Measurements: `results/linux_cloud/` (cloud server, Linux; commit 996646f) and `results/windows/20260926_1300/`
(laptop, Windows 11; commit 1d97995), both built from commit d03fd81.

Sizes are those of the encoding in `src/codec.c` (for N_V = 10^4 they coincide with the measured encodings); the
seed ciphertexts are counted with their 32-byte payload only. The benchmarks use the set for
`N_V = 10^4` unless `-d` and `-q` are given.

## License and citation

MIT License (see `LICENSE`). If you use this code, please cite the paper (see `CITATION.cff`).
