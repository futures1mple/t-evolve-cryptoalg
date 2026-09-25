#include <string.h>
#include "ring.h"

ring_ctx RING;

static uint64_t mulmod_slow(uint64_t a, uint64_t b, uint64_t q) { return (uint64_t)(((u128)a * b) % q); }

static uint64_t powmod_slow(uint64_t a, uint64_t e, uint64_t q) {
    uint64_t r = 1 % q;
    a %= q;
    while (e) {
        if (e & 1) r = mulmod_slow(r, a, q);
        a = mulmod_slow(a, a, q);
        e >>= 1;
    }
    return r;
}

static int is_prime(uint64_t n) {                 /* deterministic Miller-Rabin for n < 2^64 */
    if (n < 2) return 0;
    static const uint64_t bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int i = 0; i < 12; i++) if (n % bases[i] == 0) return n == bases[i];
    uint64_t d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        uint64_t x = powmod_slow(bases[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int comp = 1;
        for (int r = 1; r < s; r++) {
            x = mulmod_slow(x, x, n);
            if (x == n - 1) { comp = 0; break; }
        }
        if (comp) return 0;
    }
    return 1;
}

uint64_t mod_pow(uint64_t a, uint64_t e) { return powmod_slow(a, e, RING.q); }
uint64_t mod_inv(uint64_t a) { return powmod_slow(a, RING.q - 2, RING.q); }

int ring_init(uint64_t q) {
    if (q >= (1ULL << 50) || q % 32 != 17 || !is_prime(q)) return -1;
    memset(&RING, 0, sizeof RING);
    RING.q = q;
    uint64_t inv = 1;                                  /* q^{-1} mod 2^64 by Newton iteration */
    for (int i = 0; i < 6; i++) inv *= 2 - q * inv;
    RING.qinv = (uint64_t)0 - inv;
    uint64_t r1 = (uint64_t)(((u128)1 << 64) % q);
    RING.r2 = mulmod_slow(r1, r1, q);
    RING.logq = 0;
    while ((q >> RING.logq) != 0) RING.logq++;
    /* primitive 16th root: x^((q-1)/16) with (.)^8 = -1 */
    for (uint64_t x = 2;; x++) {
        uint64_t z = powmod_slow(x, (q - 1) / 16, q);
        if (powmod_slow(z, 8, q) == q - 1) { RING.zeta16 = z; break; }
    }
    /* level l, block b: modulus X^{256/2^l} - zeta^{e}, split with root zeta^{e/2}.
       exponents: level 0: {8}; children of e: e/2 and e/2+8 */
    unsigned ex[4] = {8}, nex[8];
    unsigned cnt = 1;
    for (int l = 0; l < 3; l++) {
        for (unsigned b = 0; b < cnt; b++) {
            unsigned h = ex[b] / 2;
            uint64_t r = powmod_slow(RING.zeta16, h, q);
            RING.root_m[l][b] = mulmod_slow(r, r1, q);
            RING.iroot_m[l][b] = mulmod_slow(powmod_slow(r, q - 2, q), r1, q);
            nex[2 * b] = h;
            nex[2 * b + 1] = h + 8;
        }
        cnt *= 2;
        if (l < 2) for (unsigned b = 0; b < cnt; b++) ex[b] = nex[b];
    }
    for (unsigned f = 0; f < TV_NFAC; f++)
        RING.fac_m[f] = mulmod_slow(powmod_slow(RING.zeta16, nex[f], q), r1, q);
    RING.inv8_m = mulmod_slow(powmod_slow(8, q - 2, q), r1, q);
    return 0;
}

void poly_zero(poly *r) { memset(r, 0, sizeof *r); }
void poly_add(poly *r, const poly *a, const poly *b) { for (int i = 0; i < TV_N; i++) r->c[i] = mod_add(a->c[i], b->c[i]); }
void poly_sub(poly *r, const poly *a, const poly *b) { for (int i = 0; i < TV_N; i++) r->c[i] = mod_sub(a->c[i], b->c[i]); }
void poly_neg(poly *r, const poly *a) { for (int i = 0; i < TV_N; i++) r->c[i] = mod_neg(a->c[i]); }
void poly_scale(poly *r, const poly *a, uint64_t s_mont) { for (int i = 0; i < TV_N; i++) r->c[i] = mont_mul(a->c[i], s_mont); }
void poly_to_mont(poly *a) { for (int i = 0; i < TV_N; i++) a->c[i] = to_mont(a->c[i]); }
void poly_from_int(poly *r, const int64_t *x) { for (int i = 0; i < TV_N; i++) r->c[i] = mod_from_int(x[i]); }

/* Cooley-Tukey: (f_lo + X^m f_hi) mod (X^m -+ r) = f_lo +- r f_hi */
void poly_ntt(poly *a) {
    unsigned m = TV_N / 2, blocks = 1;
    for (int l = 0; l < 3; l++) {
        for (unsigned b = 0; b < blocks; b++) {
            uint64_t *f = a->c + 2 * m * b;
            uint64_t r = RING.root_m[l][b];
            for (unsigned i = 0; i < m; i++) {
                uint64_t t = mont_mul(f[i + m], r);
                uint64_t lo = f[i];
                f[i] = mod_add(lo, t);
                f[i + m] = mod_sub(lo, t);
            }
        }
        m >>= 1;
        blocks <<= 1;
    }
}

/* Gentleman-Sande inverse; the factor 1/8 is applied at the end */
void poly_invntt(poly *a) {
    unsigned m = TV_FDEG, blocks = 4;
    for (int l = 2; l >= 0; l--) {
        for (unsigned b = 0; b < blocks; b++) {
            uint64_t *f = a->c + 2 * m * b;
            uint64_t ir = RING.iroot_m[l][b];
            for (unsigned i = 0; i < m; i++) {
                uint64_t x = f[i], y = f[i + m];
                f[i] = mod_add(x, y);
                f[i + m] = mont_mul(mod_sub(x, y), ir);
            }
        }
        m <<= 1;
        blocks >>= 1;
    }
    for (int i = 0; i < TV_N; i++) a->c[i] = mont_mul(a->c[i], RING.inv8_m);
}

static inline uint64_t redc_wide(u128 T) {       /* T * 2^-64 mod q for any T < 2^128 */
    return mod_add(redc((uint64_t)T), (uint64_t)(T >> 64) % RING.q);
}

/* accumulate lo/hi convolution of one factor */
static inline void acc_factor(u128 *lo, u128 *hi, const uint64_t *a, const uint64_t *b) {
    for (int i = 0; i < TV_FDEG; i++) {
        uint64_t ai = a[i];
        u128 *l = lo + i;
        int j = 0;
        for (; j < TV_FDEG - i; j++) l[j] += (u128)ai * b[j];
        for (; j < TV_FDEG; j++) hi[i + j - TV_FDEG] += (u128)ai * b[j];
    }
}

static void finish_factor(uint64_t *out, const u128 *lo, const u128 *hi, uint64_t fac_m, int add) {
    for (int k = 0; k < TV_FDEG; k++) {
        uint64_t l = redc_wide(lo[k]);                 /* sum a*b_mont*2^-64 = sum a*b */
        uint64_t h = redc_wide(hi[k]);
        uint64_t v = mod_add(l, mont_mul(h, fac_m));
        out[k] = add ? mod_add(out[k], v) : v;
    }
}

void poly_basemul_mont(poly *r, const poly *a, const poly *b) {
    for (int f = 0; f < TV_NFAC; f++) {
        u128 lo[TV_FDEG] = {0}, hi[TV_FDEG] = {0};
        acc_factor(lo, hi, a->c + f * TV_FDEG, b->c + f * TV_FDEG);
        finish_factor(r->c + f * TV_FDEG, lo, hi, RING.fac_m[f], 0);
    }
}

void poly_dot_mont(poly *r, const poly *A, const poly *y, size_t len) {
    for (int f = 0; f < TV_NFAC; f++) {
        u128 lo[TV_FDEG] = {0}, hi[TV_FDEG] = {0};
        for (size_t j = 0; j < len; j++) acc_factor(lo, hi, y[j].c + f * TV_FDEG, A[j].c + f * TV_FDEG);
        finish_factor(r->c + f * TV_FDEG, lo, hi, RING.fac_m[f], 0);
    }
}

void poly_mul_ref(poly *r, const poly *a, const poly *b) {
    uint64_t out[TV_N] = {0};
    for (int i = 0; i < TV_N; i++)
        for (int j = 0; j < TV_N; j++) {
            uint64_t p = mulmod_slow(a->c[i], b->c[j], RING.q);
            if (i + j < TV_N) out[i + j] = mod_add(out[i + j], p);
            else out[i + j - TV_N] = mod_sub(out[i + j - TV_N], p);
        }
    memcpy(r->c, out, sizeof out);
}

void poly_mul_sparse(poly *r, const poly *a, const uint8_t *pos, const int8_t *sgn, int w) {
    uint64_t out[TV_N] = {0};
    for (int k = 0; k < w; k++) {
        unsigned p = pos[k];
        int s = sgn[k];
        for (unsigned i = 0; i < TV_N; i++) {
            unsigned idx = i + p;
            int neg = (s < 0) ^ (idx >= TV_N);
            if (idx >= TV_N) idx -= TV_N;
            out[idx] = neg ? mod_sub(out[idx], a->c[i]) : mod_add(out[idx], a->c[i]);
        }
    }
    memcpy(r->c, out, sizeof out);
}
