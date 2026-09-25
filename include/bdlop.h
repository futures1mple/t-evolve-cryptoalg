/* BDLOP commitments in the multi-message form of EVOLVE:
 *   Com(m; r) = C r + (0^d ; m),  C = (A ; B) uniform in R_q^{(d+L) x (2d+L)}. */
#ifndef TV_BDLOP_H
#define TV_BDLOP_H

#include "params.h"

typedef struct {
    tv_params prm;
    uint8_t seed[TV_SEEDBYTES];
    uint8_t digest[TV_HASHBYTES];   /* hash of all public parameters */
    poly *C;                        /* rows x mu, NTT + Montgomery form, row-major */
    uint64_t *H;                    /* (n+1-t) x (n+1) parity-check matrix, plain form */
    uint64_t *H_m;                  /* same, Montgomery form */
    poly *Bsum;                     /* sum of the L rows of B (mu polys, NTT+Montgomery), for the weight */
} tv_pub;

/* a commitment: rows polynomials in coefficient form; first d are c_1 = A r */
typedef struct { poly *p; } tv_com;

int tv_pub_init(tv_pub *pub, const tv_params *prm, const uint8_t seed[TV_SEEDBYTES]);
void tv_pub_free(tv_pub *pub);

static inline const poly *C_row(const tv_pub *pub, int i) { return pub->C + (size_t)i * pub->prm.mu; }

/* NTT of an integer vector of mu polynomials (out: mu polys in NTT form) */
void vec_ntt_from_int(poly *out, const int64_t *x, int k);
/* out[i] = sum_j C[row_i][j] * x[j] for rows row0..row0+nrows-1; xh in NTT form; out in coeff form */
void C_mul(poly *out, const tv_pub *pub, int row0, int nrows, const poly *xh);
/* commit to m (L polys) with integer randomness r (mu polys) */
void tv_commit(poly *c, const tv_pub *pub, const poly *m, const int64_t *r);

#endif
