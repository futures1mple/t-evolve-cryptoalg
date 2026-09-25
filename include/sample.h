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
double prg_unif(prg *p);                  /* uniform in (0,1], 53-bit resolution */
unsigned prg_below(prg *p, unsigned bound); /* uniform in [0, bound), bound <= 65536 */

void sample_uniform_poly(poly *r, prg *p);

/* Discrete Gaussian D_{Z,sigma} (rho(x) = exp(-x^2/(2 sigma^2))).
 * Candidate: round(y), y ~ N(0, sigma^2) by Box-Muller; accepted with probability
 * I(0)/I(x), I(x) = int_{-1/2}^{1/2} exp(-(u^2+2xu)/(2 sigma^2)) du, which makes the
 * output D_{Z,sigma} up to floating-point precision. For sigma = 1, gauss_vec uses an integer
 * table instead (platform-independent). */
typedef struct {
    double sigma;
    int use_table;
    int xmax;
    double *acc;          /* acc[x] = I(0)/I(x) for 0 <= x <= xmax (small sigma) */
    double fast;          /* large sigma: accept at once if U < fast */
    double spare; int has_spare;
} gauss_sampler;

int gauss_init(gauss_sampler *g, double sigma);
void gauss_free(gauss_sampler *g);
int64_t gauss_sample(gauss_sampler *g, prg *p);
void gauss_vec(gauss_sampler *g, prg *p, int64_t *out, size_t n);
double gauss_accept_prob(double sigma, int64_t x);   /* exposed for tests */

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
