# Parameter chain

This note derives every parameter of the implementation from the witnesses it actually
handles. The chain is: witnesses → bounds on their norms → rejection parameters → statistical
error and success probability → extraction bounds → `beta_SIS` → lattice estimator → `(d, q)` and
sizes. The C function that computes each step is named, and `build/param_report` prints all of
them for one configuration.

Throughout, `N = 256`, `mu = 2d + L` is the width of a commitment, `sigma = 1` is the parameter
of the commitment randomness. All statistical errors are budgeted so that their sum over at most
2^40 attempts (of ballot proofs and aggregation proofs together) stays below 2^-128 (Section 7): the
tail bounds use `tau = 171 ln 2` for the ballot proof and `tau = 187 ln 2` for the aggregation, and
each rejection step has error at most `e^{-r^2/2} < 2^-172` with `r = 309/20`.

## 1. Two facts used everywhere

* **Subgaussianity.** `D_{Z,sigma}` (with `rho(x) = exp(-x^2 / 2 sigma^2)`) is sigma-subgaussian:
  `E[exp(u X)] <= exp(sigma^2 u^2 / 2)` for every real `u` [Micciancio–Peikert 2012, Lemma 2.8,
  converted from their parameter `s = sqrt(2 pi) sigma`]. A sum of `v` independent copies is
  `sqrt(v) sigma`-subgaussian.
* **Quadratic forms.** If `x` has independent sigma-subgaussian coordinates and `Phi` is fixed,
  then with `Sigma = Phi^T Phi`
  `Pr[ ||Phi x||^2 > sigma^2 (tr Sigma + 2 sqrt(tr(Sigma^2) tau) + 2 ||Sigma|| tau) ] <= e^-tau`
  [Hsu, Kakade, Zhang 2012, Thm. 1]. We use `tr(Sigma^2) <= ||Sigma|| tr Sigma`. This is `hkz()` in
  `src/agg.c` and `shift_bound()` in `src/params.c`.

Rademacher vectors (entries ±1) are 1-subgaussian, so the same bound applies to them.

## 2. Ballot proof (`src/params.c`, `src/ballot.c`)

**Witnesses.** The voter's own randomness `r_0, ..., r_n`, each `D_sigma^{N mu}`; with the
seed variant, `r_k = SampleD_sigma(X(par, id, k, s_k))` for `k >= 1`, where `X` is SHAKE256
with domain byte 0x20 over (parameter digest || id || k || s_k). The inputs of `X` differ for
different `(id, k)`, so the leaves of one aggregation column are independent even if a voter
reuses a seed (Section 5 below).

**Shift.** `s = (c r_0, ..., c r_n, f_1 r_0, ..., f_L r_0)`, where `c` and the `f_a` have 60
coefficients in {-1, 1}. For fixed challenges `s = Phi r` with
`tr Sigma = 60 N mu (n + 1 + L)` exactly (each rotation of a challenge has squared Frobenius norm
`60 N`), and `||Sigma|| <= 3600 (1 + L)` (`||Rot(c)|| <= ||c||_1 = 60`; `r_0` appears in `1 + L`
blocks). Hence `||s|| <= T` except with probability `e^-tau`, with

    T^2 = sigma^2 (tr + 2 sqrt(3600 (1+L) tr tau) + 2 * 3600 (1+L) tau).

**Rejection.** `sigma_J` is the smallest value `>= alpha T` that the sampler produces (§5), with
`alpha = 141/10`, and `M = exp(r/alpha + 1/(2 alpha^2))` with `r = 309/20`, i.e.
`ln M = 43669/39762` exactly (`M = 2.999`). The acceptance test is exact: with the integers
`zs = <z, s>` and `ss = ||s||^2` (128-bit arithmetic with overflow checks) and the integer
`sigma_J^2`, the response is accepted with probability `min(1, exp((ss - 2 zs) / (2 sigma_J^2) - ln M))`,
drawn by the algorithm of Canonne–Kamath–Steinke (Alg. 1 and Prop. 33 of arXiv:2004.00010v6:
Bernoulli(exp(-g)) for rational `g`, integer arithmetic only; `reject_accept` in `src/sample.c`).
The only source of error is the PRG (and the cap of Section 7, below 2^-(2^23)). As in
[Lyubashevsky 2012, Lemma 4.5 and Thm. 4.6], with `<z, s>` one-sided subgaussian, the ratio of the
densities exceeds `M` with probability at most `e^{-r^2/2} = 2^-172.2`; conditioned on `||s|| <= T`,
the accepted responses are within statistical distance `e^{-r^2/2}/M` of `D_{sigma_J}` over the whole
concatenation, independently of the witness. Without the condition we add `e^-tau`. **Per attempt the
simulation error is at most `e^{-r^2/2}/M + 2^-171` plus the freshness term of the paper**, and it enters
`epsilon_zk` of the privacy theorem once per attempt of every honest ballot.

Compared with the heuristic choice `sigma_J = 11 sqrt(60 N mu (n+1+L))` used before, `T` is
larger by a factor 1.9–2.1 (about one bit), so the rigorous choice costs about 3% of the ballot
(6% of the proof). The measured ratio `||s|| / T` never exceeded 0.55 in 10 000 ballots
(`results/linux_cloud/ballot_rej_summary.csv`).

**Response bound.** `B_J = 2 sigma_J sqrt(N mu (n + 1 + 2L))`, checked exactly as
`||z||^2 <= B_J^2 = 4 sigma_J^2 N mu (n + 1 + 2L)` (an integer); the concatenation of all
responses exceeds it with probability below `2^{-N mu (n+1+2L)/4}` [BL17, (18)]. Measured:
`||responses|| / B_J <= 0.51`.

## 3. Aggregation (`src/agg.c`)

**Witnesses.** For authority `k`:

* `Pi_open`, one column per node: the leaf randomness `r_{x,k}` (seed-derived, so distributed as
  `D_sigma` up to the adversary's choice of seeds; the inputs `(par, id, k, s)` of `X` are pairwise
  distinct across leaves, so in the ROM different leaves are independent even for equal seeds) and
  the fresh randomness `rho_u` of inner nodes.
* `Pi_zero`, one column per inner node: `rho_u - sum_{children} r_child`, a sum of `1 + #children`
  independent vectors. A node's randomness appears in its own `Pi_open` column, its own `Pi_zero`
  column and its parent's `Pi_zero` column, so columns of different proofs and levels are
  dependent; the bounds below need only per-column bounds.
* Blocks: `Pi_open` takes all nodes (leaves first, then level by level) in consecutive blocks of
  `k`; `Pi_zero` takes the inner nodes of each level in consecutive blocks of `k`.

**Grinding.** A dishonest voter chooses its seeds and can try up to `Q = 2^64` of them. Every
bound below that involves a leaf takes a union bound over these choices by adding `g ln Q` to
`tau`, where `g` is the number of leaves in the column (1 for a leaf column of `Pi_open`, the
number of children for a level-1 column of `Pi_zero`, 0 otherwise).

**First shift** `B1 = (c_{1,e} s)_{columns}`: no mixing across columns. Per column, with
variance multiplicity `v` (1, or `1 + #children`):
`||c1 s_col||^2 <= sigma^2 v hkz(60 N mu, 3600, tau_col + ln Q)`; `T1^2` is the sum over all columns
(the union over columns is included in `tau_col`; the extra `ln Q` covers the adversary's ability to
steer `c1 = H_agg(par, BB_1)` through its own posts in round 1). This bounds the Frobenius norm directly and is
tighter than the route through `s_part` in [BL17, (22)].

**Second shift** `B2 = (S_e C_{2,e})_e`. With a binary `C2` (as in [BL17]) the columns of `C2`
have mean 1/2, so `||S C2||` contains `||S 1||/2`, which an adversary controlling many leaves in a
block can make up to `k` times larger than for independent columns by aligning its leaves; a
bound that holds for adversarial leaves has to go through `||S||_F s1(C2)` with
`s1(C2) ~ sqrt(k ell)/2`. With a **signed** `C2 in {-1,1}^{k x ell}` (`C2_SIGNED`), for every fixed
`S` the vector of all entries of all `C2` is Rademacher and independent of `S`, and
`||S C2||_F^2` is a quadratic form in it with trace `ell ||S||_F^2`:

    T2^2 = ell F + 2 sqrt(ell tau sum_e F_e^2) + 2 tau max_e F_e,   F_e = sum_{columns in e} b_col^2,

where `b_col^2 = sigma^2 v hkz(N mu, 1, tau_col)` bounds the squared norm of one column. This holds
for any alignment of the columns and needs only the per-column bounds.

**Class bounds and smaller trees.** The parameters are chosen for `N_V`, a hard upper bound on the
number of accepted ballots, while the adversary decides how many ballots are accepted. To make the
bounds valid for every tree with at most `N_V` leaves, whatever its block layout, `T1` and `T2` are
computed with one bound per class of columns: every column of `Pi_open` is counted as a leaf
(`v = 1`, `g = 1`), every column of `Pi_zero` as a node with the full fan-in (`v = 1 + l`, with
`g = l` at level 1 and `g = 0` above), and `tau_col` uses the number of columns of the tree for
`N_V`. Every block then weighs (number of its columns) x (class bound), all blocks of a sequence are
full except the last, and the number of nodes of every level is nondecreasing in the number of
leaves; hence `F`, `sum_e F_e^2`, `max_e F_e` and `T1^2` are nondecreasing in the number of leaves,
and the values for `N_V` bound those of every smaller tree. Compared with the exact per-column
values this costs at most 0.003 bit in `sigma1`, `sigma2` and `beta_SIS`.

**Rejection.** `sigma_i` is the smallest achievable value `>= alpha' T_i` with the rational
`alpha' = 281/10 = 28.1`, and `M' = exp(r/alpha' + 1/(2 alpha'^2))` with `r = 309/20`, so that
`ln M' = 86929/157922` exactly (`M' = 1.7340`, about `sqrt 3`). Each of the two steps is applied once
per attempt to the concatenation over all `E` blocks, with the exact test described for the ballot
proof. An attempt succeeds with probability about `1/M'^2 = 0.333`, independently of `E`; the
simulation error per attempt is at most `2 e^{-r^2/2}/M' + (l(E+1)+1) e^-tau + eps_fresh` with
`tau = 187 ln 2` (Lemma 13 of the paper; `(l(E+1)+1) e^-tau <= 2^-171` for `E <= 2^11`). The norm
checks of every block (`IsSmall`) are exact integer comparisons: rows of `Z1` with
`||.||^2 <= 2 k sigma1^2`; entries of `Z2` with `z^2 <= 81 sigma2^2`; each ring component of each
column of `Z2` with `||.||^2 <= 2 N sigma2^2`.

**Extraction.** From two accepting transcripts with different first challenges [BL17, Thm 3.5]:

* binary `C2`: openings with slack 1 and `||S_bar||_inf <= 2 ||Z2 - Z2'||_inf <= 36 sigma2`;
  M-SIS otherwise, with `l_inf` bound `2^{3/2} k sigma1 + 2^{5/2} sqrt(60 N) sigma2`;
* signed `C2`: `C2 - C2' = 2D` with `D` in {0, ±1}; Lemma 2.3 of [BL17] applied to `D` (whose entries
  are 0 or a fixed sign with probability 1/2 each) gives the same probability bound
  `q_H^2 2^{-ell-1} + 2^{-ell}`, and the extractor obtains `2 S_bar` with
  `||2 S_bar||_inf <= 2 ||Z2 - Z2'||_inf <= 36 sigma2`, i.e. openings with **slack 2**; the M-SIS
  bound becomes `2^{5/2} k sigma1 + 2^{5/2} sqrt(60 N) sigma2` because columns of `C2 - C2'` have norm
  at most `2 sqrt(k)`. Proved in the paper (Lemma 11 and Appendix B).

With `B_A = 36 sigma2 sqrt(N mu)` (`l_2` bound of an extracted opening), slack `s` (1 or 2),
`eta = 120`, `beta_h = 2 sigma sqrt(N mu)` and `B_bar = 2 B_J`:

    beta_SIS = max( s beta_h + B_A,  (l + 2) B_A,  s B_bar + eta B_A,  eta (B_bar + B_bar),  beta_BL sqrt(N mu) ).

**Rounds.** With `kappa = 4` attempts per round, an authority has no non-aborting attempt with
probability `(2/3)^4 = 0.198`. A round succeeds if at least `t` authorities have one: with `n = 4`,
`t = 3` and all authorities honest this happens with probability 0.82 (1.22 rounds on average); with
only `t` honest authorities, with probability at least `0.80^3 = 0.52` (at most 1.94 rounds). An
authority first announces the index of its non-aborting attempt by a deadline `Delta_1`; if at least
`t` authorities have announced one, they reveal their responses by `Delta_2`. If fewer than `t`
announce, no responses are published and the round costs the work of both rounds but no publication
beyond round 1. If at least `t` announce but fewer than `t` valid contributions arrive by `Delta_2`,
an announcing authority has misbehaved: it is recorded as faulty and excluded from later rounds,
and the honest authorities have already published their responses (one extra round-2 publication
each). Misbehaviour cannot make a round fail that would otherwise succeed (the contributions of `t`
honest authorities with a non-aborting attempt complete it), so it adds publications in at most
`n - t` rounds but no rounds. Contributions carry their round index, and only `t` contributions of
the same round are combined. The benchmark has no faulty authorities.

**Common leaves.** The leaves are the share commitments of the accepted ballots, which are the
same for every authority (a ballot with a valid complaint is removed for all). The hash of the leaf
commitments is bound into `G` and into the digest of `BB_1`.

## 4. Choice of (d, q) (`tools/chain.py`)

For every `N_V` and ballot format, the smallest `d` and, for it, the smallest bit length `b` such
that the largest prime `q = 17 (mod 32)` below `2^b` exceeds `beta_SIS` and both estimates are at
least `2^128` (elements of `Z_q` then take exactly `b` bits), using the lattice estimator
(commit 53da598), MATZOV cost model for M-LWE (`n = dN`, `m = (d+L)N`, secret and error `D_1`) and the
Euclidean norm for M-SIS (`n = dN`, `m = mu N`). Every estimate, with all attacks, is written to
`results/params_chain_*.jsonl`.

## 5. Sampling (`src/sample.c`, `tools/gen_cdt.py`)

Integer arithmetic only, so every platform produces the same output from the same stream; the
statistical distance of each sampler from the ideal distribution is bounded.

* `sigma = 1` (all commitment randomness, including the seed-derived randomness of the shares): a
  table of `floor(2^256 Pr[|X| <= x])` for `x < 20`, generated with 512-bit arithmetic by
  `tools/gen_cdt.py`; `r` is uniform in `[0, 2^256)`, drawn 64 bits at a time only as far as the
  comparisons need, and `|X| = #{x : r >= CDT[x]}`, with a uniform sign. Statistical distance at most
  `21 * 2^-256 + 2^-264 < 2^-251` per sample.
* Base sampler `D_{Z,256}`: the same construction with 4856 entries, distance `< 2^-243` per sample.
* Large `sigma` (masks): convolution as in Micciancio–Walter (CRYPTO 2017). Level 1 returns
  `a1 x + b1 x'` for two base samples, level 2 returns `a2 y + b2 y'` for two level-1 samples,
  with `gcd(a_l, b_l) = 1` and `max(a_l, b_l)^2 <= (pi / eta^2) sigma_in^2`, where
  `eta = sqrt(ln(2 + 2^221) / pi) >= eta_eps(Z)` for `eps = 2^-220` [MR04, Lemma 3.3]. By the
  convolution theorem [MP13, Thm. 3.3; MW17, Thm. 2.1 and Lemma 5.1] the output built from exact base
  samples has relative error at most `2^levels * 2 eps <= 2^-217` with respect to `D_{Z,sigma_out}`,
  `sigma_out = 256 sqrt((a1^2 + b1^2)(a2^2 + b2^2))`; replacing the (at most four) base samples by
  the table adds at most `4 * 2^-243`. Hence at most `2^-216` per sample. The largest achievable
  `sigma_out` is about `1.9 * 10^8 = 2^27.5`; the largest mask here is `sigma_1 = 2^26.8` (`N_V = 10^6`).
* `sigma_out` is the smallest achievable value `>= alpha T` (the search is in `plan()`); it exceeds
  the target by at most `0.02%` for all masks here. Every bound (`sigma_J`, `B_J`, `sigma_1`,
  `sigma_2`, `B_A`, `beta_BL`, rejection) is computed with `sigma_out`, via `gauss_round_sigma()`.
  A larger `sigma` keeps the rejection bound valid with the same `M`.
* Per attempt the ballot proof draws `(n + 1 + 2L) N (2d + L)` mask coefficients (26 880 for
  `n = 4, L = 1`; 153 600 for `n = 4, L = 10`), so the sampling error is at most
  `(n + 1 + 2L) N (2d + L) 2^-216`, i.e. `2^-201.3` and `2^-198.8`. The aggregation draws
  `E N (2d + L) (k + ell)` coefficients per attempt, at most `2^33` for `N_V = 10^6`, hence at most
  `2^-183`. The randomness of `sigma = 1` is replaced by ideal samples in one hybrid for all at most
  `Q = 2^64` outputs of `X` (Lemma 13 takes a union bound over them), at most
  `Q N (2d + L) 2^-251 <= 2^-174.4` once per experiment.

## 6. Status of the statements

| Statement | Status |
|---|---|
| Formulas above; the numbers printed by `param_report` and `chain.py` | established (given the cited lemmas) |
| Samplers | statistical distance bounded (see §5): `< 2^-251` (`sigma = 1`), `< 2^-216` (masks) per sample |
| Ballot proof: rejection error per attempt `e^{-r^2/2}/M + 2^-171` | established |
| Acceptance tests and norm checks | exact: integer arithmetic with overflow checks and exact Bernoulli(exp(-g)) sampling (cap of §7 below 2^-(2^23)); only the PRG remains |
| Total statistical error of the simulations over 2^40 attempts | below 2^-128 (§7) |
| Aggregation with signed `C2`: bounds on both shifts (Lemma 13), extraction with slack 2 (Lemma 11) | proved in the paper; new in this version, independent verification pending |
| Grinding: a union bound over `2^64` seeds per leaf | established in the ROM; `Q = 2^64` is an assumption on the adversary |
| Expected attempts (ballot ~3, aggregation ~3) and the resulting costs | measured; plausible for other machines |
| 128-bit security of the chosen `(d, q)` | estimate of the lattice estimator; the reduction of the paper is not tight |

## 7. Error budget and integer ranges

**Budget.** Per attempt, for the largest parameter sets (`N_V = 10^6`, `E = 2142`; ballot with
`n = 4`, `L = 10`):

| Source | Ballot attempt | Aggregation attempt |
|---|---|---|
| Tail bounds (Lemma 8, resp. Lemma 13: `(l(E+1)+1) e^-tau`) | `2^-171` | `2^-171.0` |
| Rejection (`e^{-r^2/2}/M` per step; one, resp. two steps) | `2^-173.8` | `2^-172.0` |
| Mask samplers | `2^-198.8` | `2^-183.0` |
| Freshness of the hashed first messages, programming | `< 2^-1000` | `2^-192` (salt, `q_G <= 2^64`, 256-bit salts) |
| Cap of the Bernoulli counter | `< 2^-(2^23)` | `< 2^-(2^23)` |
| **Sum** | **`< 2^-170.7`** | **`< 2^-170.3`** |

Over at most `2^40` attempts of both kinds together this is below `2^-130.3`; the hybrid for the
`sigma = 1` randomness adds `2^-174.4` once. The total statistical error of the simulations is
therefore below `2^-128`. The computational terms (hiding, IND-CCA, EUF-CMA, M-SIS) are separate.

**Integer ranges.** All acceptance decisions use integers only; every intermediate value is either
bounded below by the table or checked at run time.

| Quantity | Range | How it is ensured |
|---|---|---|
| Coefficients of masks, responses, witnesses | `|x| < 2^63` (int64) | Gaussian samples below `2^27.5 * 40`; decoded inputs are int64 |
| Products `x * y` of two coefficients | `< 2^126` | always, in `__int128` |
| Sums of products (`<z, s>`, `||s||^2`, norms) | `< 2^127` | prover side: `chk_madd`, aborts on overflow; verifier side: `sat_madd` saturates, so an oversized input fails the norm check |
| `sigma^2` | `< 2^56` | `sigma_out <= 2^27.5` (§5) |
| `cn/cd = ln M` | `cd < 2^18`, `cn < 2^17` | the two constants above |
| `2 sigma^2 cd` (denominator) | `< 2^106` | checked in `bern_exp` |
| `(2 zs - ss) cd + 2 sigma^2 cn` (numerator) | `< 2^127` | `|zs|, ss <= 2^100` and checked multiplications in `reject_accept` |
| Counter `K` of CKS Alg. 1 | `K <= 2^20`, so `den * K < 2^126` | cap; reached with probability `< 1/(2^20 - 1)! < 2^-(2^23)` per call, the only case in which the output may differ |
| Loop count of `bern_exp` (factors `exp(-1)`) | `<= num/den <= 2^101` | each factor continues only with probability `e^-1`, so fewer than 2 iterations are expected; the count affects running time only |

Values outside these ranges stop the program with an error message (`tv_range_abort`); they never
produce a silently wrong result. With the parameter sets of `tools/chain.py` the checks never fire.
