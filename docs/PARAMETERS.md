# Parameter chain

This note derives every parameter of the implementation from the witnesses it actually
handles. The chain is: witnesses → bounds on their norms → rejection parameters → statistical
error and success probability → extraction bounds → `beta_SIS` → lattice estimator → `(d, q)` and
sizes. The C function that computes each step is named, and `build/param_report` prints all of
them for one configuration.

Throughout, `N = 256`, `mu = 2d + L` is the width of a commitment, `sigma = 1` is the parameter
of the commitment randomness, and `tau = (128 + 40) ln 2` is the failure exponent used in every
tail bound: 2^-128 per event, with a union bound over at most 2^40 events.

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

**Rejection.** `sigma_J = 11 T`, `M = exp(12/11 + 1/242) = 2.989`. By [Lyubashevsky 2012,
Thm. 4.6/Lemma 4.7], conditioned on `||s|| <= T`, one attempt accepts with probability within
`2^-100` of `1/M`, and the accepted responses are within statistical distance `2^-100/M` of
`D_{sigma_J}` over the whole concatenation, independently of the witness. Without the condition
we add `e^-tau`. **Per attempt the simulation error is at most `2^-100 / M + 2^-168`**, and it
enters `epsilon_zk` of the privacy theorem once per attempt of every honest ballot.

Compared with the heuristic choice `sigma_J = 11 sqrt(60 N mu (n+1+L))` used before, `T` is
larger by a factor 1.9–2.1 (about one bit), so the rigorous choice costs about 8% of the ballot
size. The measured ratio `||s|| / T` never exceeded 0.57 in 5000 ballots
(`results/linux_cloud/ballot_rej_summary.csv`).

**Response bound.** `B_J = 2 sigma_J sqrt(N mu (n + 1 + 2L))`; the concatenation of all
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
  dependent; blocks of `Pi_zero` are formed within one level, so that columns inside a block of
  `Pi_zero` are independent.

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

**Rejection.** `sigma_i = alpha T_i` with `alpha = 12 / ln sqrt(3)`, `M = exp(12/alpha + 1/(2 alpha^2))`
for each of the two steps, applied once per attempt to the concatenation over all `E` blocks.
An attempt succeeds with probability about `1/M^2 = 1/3`, independently of `E`; the simulation error
per attempt is at most `2 * 2^-100 / M + 2 e^-tau`. The norm checks of every block (`IsSmall`)
use: rows of `Z1` at most `sqrt(2k) sigma1`; entries of `Z2` at most `9 sigma2`; each ring
component of each column of `Z2` at most `sqrt(2N) sigma2`.

**Extraction.** From two accepting transcripts with different first challenges [BL17, Thm 3.5]:

* binary `C2`: openings with slack 1 and `||S_bar||_inf <= 2 ||Z2 - Z2'||_inf <= 36 sigma2`;
  M-SIS otherwise, with `l_inf` bound `2^{3/2} k sigma1 + 2^{5/2} sqrt(60 N) sigma2`;
* signed `C2`: `C2 - C2' = 2D` with `D` in {0, ±1}; Lemma 2.3 of [BL17] applied to `D` (whose entries
  are 0 or a fixed sign with probability 1/2 each) gives the same probability bound
  `q_H^2 2^{-ell-1} + 2^{-ell}`, and the extractor obtains `2 S_bar` with
  `||2 S_bar||_inf <= 2 ||Z2 - Z2'||_inf <= 36 sigma2`, i.e. openings with **slack 2**; the M-SIS
  bound becomes `2^{5/2} k sigma1 + 2^{5/2} sqrt(60 N) sigma2` because rows of `C2 - C2'` have norm at
  most `2 sqrt(k)`. *Status: proof sketch; the paper must state and prove this variant.*

With `B_A = 36 sigma2 sqrt(N mu)` (`l_2` bound of an extracted opening), slack `s` (1 or 2),
`eta = 120`, `beta_h = 2 sigma sqrt(N mu)` and `B_bar = 2 B_J`:

    beta_SIS = max( s beta_h + B_A,  (l + 2) B_A,  s B_bar + eta B_A,  eta (B_bar + B_bar),  beta_BL sqrt(N mu) ).

**Rounds.** With `kappa = 4` attempts per round, an authority has no non-aborting attempt with
probability `(2/3)^4 = 0.198`. A round succeeds if at least `t` authorities have one: with `n = 4`,
`t = 3` and all authorities honest this happens with probability 0.82 (1.22 rounds on average); with
only `t` honest authorities, with probability at least `0.80^3 = 0.52` (at most 1.94 rounds). An
authority first announces the index of its non-aborting attempt and reveals the responses only if
the round succeeds, so the responses of failed rounds are never published; a failed round costs the
work of round 1 and round 2 but no publication beyond round 1.

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

## 5. Sampling

* `sigma = 1` (all commitment randomness, including the seed-derived randomness of the shares):
  integer cumulative-distribution table, statistical distance about `2^-60` per sample, identical on
  every platform, so the authority re-derives exactly the voter's randomness.
* Large `sigma` (masks): Box–Muller in double precision (with the radius tail refined to 106-bit
  resolution), rounding, and an acceptance step `I(0)/I(x)` that makes the rounded output exactly
  `D_{Z,sigma}` in exact arithmetic. The floating-point error of this sampler is not bounded formally
  and is **not** included in the simulation errors above; an integer sampler would remove this gap.

## 6. Status of the statements

| Statement | Status |
|---|---|
| Formulas above; the numbers printed by `param_report` and `chain.py` | established (given the cited lemmas) |
| Floating-point sampler for large `sigma` | error not bounded formally (see §5) |
| Ballot proof: rejection error per attempt `2^-100/M + 2^-168` | established |
| Aggregation with signed `C2`: bounds on both shifts, for any alignment of the leaves | established for the bounds; the extraction lemma for signed `C2` is a proof sketch |
| Grinding: a union bound over `2^64` seeds per leaf | established in the ROM; `Q = 2^64` is an assumption on the adversary |
| Expected attempts (ballot ~3, aggregation ~3) and the resulting costs | measured; plausible for other machines |
| 128-bit security of the chosen `(d, q)` | estimate of the lattice estimator; the reduction of the paper is not tight |
