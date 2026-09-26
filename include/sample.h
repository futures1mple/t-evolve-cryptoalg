/* Deterministic randomness (SHAKE256 streams) and samplers. */
#ifndef TV_SAMPLE_H
#define TV_SAMPLE_H

#include <stdint.h>
#include <stddef.h>
#include "keccak.h"
#include "ring.h"

#define TV_CH_W 60         /* nonzero coefficients of a challenge */

typedef struct {
    keccak_state ks;
    uint8_t buf[SHAKE256_RATE];
    unsigned pos;
} prg;

/* stream = SHAKE256(domain || seed) */
void prg_init(prg *p, uint8_t domain, const uint8_t *seed, size_t seedlen);
void prg_bytes(prg *p, uint8_t *out, size_t len);
uint64_t prg_u64(prg *p);
unsigned prg_below(prg *p, unsigned bound); /* uniform in [0, bound), bound <= 65536 */

void sample_uniform_poly(poly *r, prg *p);

/* Discrete Gaussian D_{Z,sigma}, rho(x) = exp(-x^2/(2 sigma^2)), with a proven bound on the
 * statistical distance and integer arithmetic only (identical output on every platform).
 *
 *  - sigma = 1: table of Pr[|X| <= x] with 256-bit entries (include/cdt_tables.h, 20 entries);
 *    distance at most 21 * 2^-256 + 2^-264 < 2^-251 per sample.
 *  - sigma > 1: convolution of samples of a base sampler D_{Z,256} (same kind of table, 4856 entries,
 *    distance < 4857 * 2^-256 + 2^-264 < 2^-243 per sample) as in Micciancio-Walter (CRYPTO 2017):
 *    level 1 returns a1 x + b1 x', level 2 returns a2 y + b2 y' for independent inputs, with
 *    gcd(a,b) = 1 and max(a,b) <= sqrt(pi) sigma_in / eta, eta = sqrt(ln(2 + 2^221)/pi) >= eta_eps(Z)
 *    for eps = 2^-220 [MR04, Lemma 3.3].  By [MP13, Thm. 3.3] (MW17, Thm. 2.1 and Lemma 5.1) the
 *    output with exact base samples has relative error at most 2^levels * 2 eps w.r.t.
 *    D_{Z,sigma_out}, sigma_out = 256 sqrt((a1^2+b1^2)(a2^2+b2^2)); with the tables, the
 *    statistical distance is at most 8 * 2^-220 + 4 * 2^-243 < 2^-216 per sample.
 *  sigma_out is the smallest achievable value >= the requested sigma; callers must use
 *  gauss_round_sigma() so that all bounds are computed with sigma_out. */
typedef struct {
    double sigma;              /* sigma_out */
    uint64_t sigma2;           /* sigma_out^2, an integer: 1, or 65536 (a1^2+b1^2)(a2^2+b2^2) */
    int levels;                /* -1: sigma = 1 table; 0: base; 1 or 2: convolution levels */
    int32_t a[2], b[2];
} gauss_sampler;

double gauss_round_sigma(double sigma);      /* smallest achievable sigma_out >= sigma; <0 if none */
uint64_t gauss_sigma2(double sigma);         /* sigma_out^2 as an integer for an achievable sigma_out */
int gauss_init(gauss_sampler *g, double sigma);   /* sigma must be 1 or >= 256 */
void gauss_free(gauss_sampler *g);
int64_t gauss_sample(gauss_sampler *g, prg *p);
void gauss_vec(gauss_sampler *g, prg *p, int64_t *out, size_t n);

/* Exact rejection sampling (integer arithmetic only, no rounding).
 * bern_exp returns 1 with probability exp(-num/den) for num >= 0, den > 0, by the algorithm of
 * Canonne, Kamath and Steinke (NeurIPS 2020, Alg. 1): Bernoulli(exp(-g)) for g in [0,1] from
 * Bernoulli(g/K) trials, and one Bernoulli(exp(-1)) trial per unit of g above 1. The only source of
 * error is the PRG. reject_accept(p, zs, ss, sigma2, cn, cd) accepts with probability
 * min(1, exp((ss - 2 zs)/(2 sigma2) - cn/cd)), the rule of Lyubashevsky with ln M = cn/cd, where
 * zs = <z, s> and ss = ||s||^2 are exact integers. */
typedef __int128 i128;
int bern_exp(prg *p, i128 num, i128 den);
int reject_accept(prg *p, i128 zs, i128 ss, uint64_t sigma2, int64_t cn, int64_t cd);

/* Integer ranges (docs/PARAMETERS.md, Section 7). The counter of CKS20 Alg. 1 is capped at
 * TV_BERN_KMAX; denominators must be below 2^TV_BERN_DEN_BITS, so that den * K < 2^126. Inner
 * products and squared norms passed to reject_accept must lie in [-TV_DOT_MAX, TV_DOT_MAX]. Every
 * violation aborts (tv_range_abort); it never yields a silently wrong result. */
#define TV_BERN_KMAX ((unsigned __int128)1 << 20)
#define TV_BERN_DEN_BITS 106
#define TV_DOT_MAX (((i128)1 << 100))
void tv_range_abort(const char *where);

/* Exact accumulation of products of 64-bit integers in 128 bits with overflow detection.
 * chk_madd (prover side, own data): aborts on overflow. sat_madd (verifier side, possibly
 * adversarial data, nonnegative sums of squares): saturates at I128_MAX, so that a norm check
 * with any bound below 2^127 rejects. */
#define I128_MAX ((i128)(((unsigned __int128)1 << 127) - 1))
static inline i128 chk_madd(i128 acc, int64_t a, int64_t b, const char *where) {
    i128 r;
    if (__builtin_add_overflow(acc, (i128)a * (i128)b, &r)) tv_range_abort(where);
    return r;
}
static inline i128 sat_add(i128 a, i128 b) {       /* a, b >= 0 */
    i128 r;
    return __builtin_add_overflow(a, b, &r) ? I128_MAX : r;
}
static inline i128 sat_madd(i128 acc, int64_t a, int64_t b) { return sat_add(acc, (i128)a * (i128)b); }

/* challenge c in C = { c : 60 coefficients in {-1,1}, the rest 0 } */
typedef struct {
    uint8_t pos[TV_CH_W];
    int8_t sgn[TV_CH_W];
} challenge;
void sample_challenge(challenge *c, prg *p);
void challenge_dense(int8_t out[TV_N], const challenge *c);
int challenge_from_dense(challenge *c, const int8_t in[TV_N]);   /* 0 if in C */

/* EVOLVE's OR challenge space Pi = Perm(N) x {0,1}^60 acting on C */
typedef struct {
    uint8_t s[TV_N];
    uint8_t b[TV_CH_W];
} perm_ch;
void sample_perm(perm_ch *pi, prg *p);
void perm_apply(challenge *out, const perm_ch *pi, const challenge *f);   /* out = pi(f) */
void perm_apply_inv(challenge *out, const perm_ch *pi, const challenge *g); /* out = pi^{-1}(g) */

/* c * v for an integer vector v of k polynomials: out[i] = c * v[i] over Z */
void ch_mul_int(int64_t *out, const challenge *c, const int64_t *v, size_t k);

#endif
