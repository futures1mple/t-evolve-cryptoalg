/* Arithmetic in R_q = Z_q[X]/(X^256+1) for a prime q < 2^50 with q = 17 (mod 32).
 *
 * For such q, X^256+1 splits into 8 irreducible factors X^32 - zeta^e (zeta a primitive
 * 16th root of unity mod q), so a partial NTT with three levels maps R_q to 8 copies of
 * Z_q[X]/(X^32 - zeta^e). The same q makes all nonzero differences of challenges
 * invertible [Lyubashevsky-Seiler 2018, Cor. 1.2].
 *
 * Coefficients are stored as uint64_t in [0, q). Products use Montgomery reduction with
 * R = 2^64 and 128-bit intermediates (GCC/Clang/MinGW-w64 __int128). */
#ifndef TV_RING_H
#define TV_RING_H

#include <stdint.h>
#include <stddef.h>

#define TV_N 256
#define TV_NFAC 8                 /* number of NTT factors */
#define TV_FDEG (TV_N / TV_NFAC)  /* degree of each factor: 32 */

__extension__ typedef unsigned __int128 u128;

typedef struct { uint64_t c[TV_N]; } poly;   /* coefficient or NTT domain, values in [0,q) */

typedef struct {
    uint64_t q;
    uint64_t qinv;        /* -q^{-1} mod 2^64 */
    uint64_t r2;          /* 2^128 mod q */
    unsigned logq;        /* bit length of q */
    uint64_t zeta16;      /* primitive 16th root of unity */
    /* NTT roots (Montgomery form), level l (0..2), block b (0..2^l-1) */
    uint64_t root_m[3][4];
    uint64_t iroot_m[3][4];
    uint64_t fac_m[TV_NFAC];  /* zeta^e of each final factor, Montgomery form */
    uint64_t inv8_m;          /* 8^{-1}, Montgomery form */
} ring_ctx;

extern ring_ctx RING;

/* returns 0 on success, -1 if q is not a prime = 17 mod 32 below 2^50 */
int ring_init(uint64_t q);

static inline uint64_t mod_add(uint64_t a, uint64_t b) { uint64_t s = a + b; return s >= RING.q ? s - RING.q : s; }
static inline uint64_t mod_sub(uint64_t a, uint64_t b) { return a >= b ? a - b : a + RING.q - b; }
static inline uint64_t mod_neg(uint64_t a) { return a ? RING.q - a : 0; }
static inline uint64_t redc(u128 T) {           /* T < q * 2^64; returns T * 2^-64 mod q */
    uint64_t m = (uint64_t)T * RING.qinv;
    u128 t = (T + (u128)m * RING.q) >> 64;
    uint64_t r = (uint64_t)t;
    return r >= RING.q ? r - RING.q : r;
}
static inline uint64_t mont_mul(uint64_t a, uint64_t b) { return redc((u128)a * b); } /* a*b*2^-64 */
static inline uint64_t to_mont(uint64_t a) { return mont_mul(a, RING.r2); }           /* a*2^64 */
static inline uint64_t mod_mul(uint64_t a, uint64_t b) { return mont_mul(mont_mul(a, b), RING.r2); }
uint64_t mod_pow(uint64_t a, uint64_t e);
uint64_t mod_inv(uint64_t a);
static inline uint64_t mod_from_int(int64_t x) {
    int64_t r = x % (int64_t)RING.q;
    return (uint64_t)(r < 0 ? r + (int64_t)RING.q : r);
}
static inline int64_t mod_centered(uint64_t a) { return a > RING.q / 2 ? (int64_t)a - (int64_t)RING.q : (int64_t)a; }

/* polynomial operations (coefficient domain unless stated) */
void poly_zero(poly *r);
void poly_add(poly *r, const poly *a, const poly *b);
void poly_sub(poly *r, const poly *a, const poly *b);
void poly_neg(poly *r, const poly *a);
void poly_scale(poly *r, const poly *a, uint64_t s_mont);   /* r = s*a, s given in Montgomery form */
void poly_ntt(poly *a);
void poly_invntt(poly *a);
/* r = a o b in the NTT domain; b must be in Montgomery form (b*2^64) */
void poly_basemul_mont(poly *r, const poly *a, const poly *b);
void poly_to_mont(poly *a);
/* schoolbook negacyclic product, reference for tests */
void poly_mul_ref(poly *r, const poly *a, const poly *b);
/* r = c * a where c is sparse with coefficients in {-1,0,1}: pos[i], sgn[i] (+1/-1), w entries */
void poly_mul_sparse(poly *r, const poly *a, const uint8_t *pos, const int8_t *sgn, int w);
void poly_from_int(poly *r, const int64_t *x);

/* accumulate sum_j A[j] o y[j] for len terms, with A[j] in Montgomery NTT form and y[j] in
 * NTT form, using one reduction per coefficient; result in NTT domain */
void poly_dot_mont(poly *r, const poly *A, const poly *y, size_t len);

#endif
