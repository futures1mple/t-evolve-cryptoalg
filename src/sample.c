#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sample.h"

void prg_init(prg *p, uint8_t domain, const uint8_t *seed, size_t seedlen) {
    shake256_init(&p->ks);
    keccak_absorb(&p->ks, &domain, 1);
    keccak_absorb(&p->ks, seed, seedlen);
    keccak_finalize(&p->ks);
    p->pos = SHAKE256_RATE;
}

void prg_bytes(prg *p, uint8_t *out, size_t len) {
    while (len) {
        if (p->pos == SHAKE256_RATE) {
            keccak_squeeze(&p->ks, p->buf, SHAKE256_RATE);
            p->pos = 0;
        }
        size_t k = SHAKE256_RATE - p->pos;
        if (k > len) k = len;
        memcpy(out, p->buf + p->pos, k);
        p->pos += (unsigned)k;
        out += k;
        len -= k;
    }
}

uint64_t prg_u64(prg *p) {
    uint8_t b[8];
    prg_bytes(p, b, 8);
    uint64_t x = 0;
    for (int i = 7; i >= 0; i--) x = (x << 8) | b[i];
    return x;
}

double prg_unif(prg *p) { return ((double)(prg_u64(p) >> 11) + 1.0) * (1.0 / 9007199254740992.0); }

unsigned prg_below(prg *p, unsigned bound) {
    unsigned lim = 65536u - 65536u % bound;
    for (;;) {
        uint8_t b[2];
        prg_bytes(p, b, 2);
        unsigned x = b[0] | ((unsigned)b[1] << 8);
        if (x < lim) return x % bound;
    }
}

void sample_uniform_poly(poly *r, prg *p) {
    uint64_t mask = (RING.logq >= 64) ? ~0ULL : ((1ULL << RING.logq) - 1);
    for (int i = 0; i < TV_N;) {
        uint64_t x = prg_u64(p) & mask;
        if (x < RING.q) r->c[i++] = x;
    }
}

/* ---------------- discrete Gaussian ---------------- */

/* I(x) = sum_k p_k (1/2)^{2k}/(2k+1), where sum_k p_k u^{2k} = exp(-a u^2) cosh(b u),
   a = 1/(2 sigma^2), b = x / sigma^2. Returns I(x) - 1 to keep precision. */
static long double I_minus_1(double sigma, int64_t x) {
    long double a = 1.0L / (2.0L * sigma * sigma), b = (long double)x / ((long double)sigma * sigma);
    enum { K = 48 };
    long double ea[K], cb[K];            /* (-a)^i/i!,  b^{2j}/(2j)! */
    ea[0] = 1; cb[0] = 1;
    for (int i = 1; i < K; i++) {
        ea[i] = ea[i - 1] * (-a) / i;
        cb[i] = cb[i - 1] * b * b / ((2.0L * i - 1) * (2.0L * i));
    }
    long double s = 0, pw = 1;           /* pw = (1/2)^{2k} */
    for (int k = 1; k < K; k++) {
        pw *= 0.25L;
        long double pk = 0;
        for (int i = 0; i <= k; i++) pk += ea[i] * cb[k - i];
        s += pk * pw / (2.0L * k + 1);
    }
    return s;
}

double gauss_accept_prob(double sigma, int64_t x) {
    long double s0 = I_minus_1(sigma, 0), sx = I_minus_1(sigma, x);
    return (double)((1.0L + s0) / (1.0L + sx));
}

int gauss_init(gauss_sampler *g, double sigma) {
    memset(g, 0, sizeof *g);
    if (!(sigma >= 1.0)) return -1;
    g->sigma = sigma;
    if (sigma < 64) {
        g->use_table = 1;
        g->xmax = (int)ceil(10 * sigma) + 2;
        g->acc = malloc(sizeof(double) * (size_t)(g->xmax + 1));
        if (!g->acc) return -1;
        for (int x = 0; x <= g->xmax; x++) g->acc[x] = gauss_accept_prob(sigma, x);
    } else {
        /* for |x| <= 9 sigma: I(0)/I(x) >= (1 - 1/(24 s^2)) / (1 + 3.4/s^2) >= 1 - 4/s^2 */
        g->fast = 1.0 - 4.0 / (sigma * sigma);
    }
    return 0;
}

void gauss_free(gauss_sampler *g) { free(g->acc); g->acc = NULL; }

static double normal01(gauss_sampler *g, prg *p) {
    if (g->has_spare) { g->has_spare = 0; return g->spare; }
    double u1 = prg_unif(p), u2 = prg_unif(p);
    /* refine small u1 (the tail of the radius) to 106-bit resolution, so that the discrete set of
       radii is fine enough for every sigma used here (sigma < 2^40) */
    if (u1 < 0x1p-20) u1 = (u1 * 0x1p53 - 1.0 + prg_unif(p)) * 0x1p-53;
    double rad = sqrt(-2.0 * log(u1)), th = 6.283185307179586476925 * u2;
    g->spare = rad * sin(th);
    g->has_spare = 1;
    return rad * cos(th);
}

int64_t gauss_sample(gauss_sampler *g, prg *p) {
    for (;;) {
        double y = g->sigma * normal01(g, p);
        int64_t x = (int64_t)llround(y);
        int64_t ax = x < 0 ? -x : x;
        double u = prg_unif(p);
        if (g->use_table) {
            if (ax > g->xmax) continue;          /* never happens: Box-Muller gives |y| < 8.6 sigma */
            if (u <= g->acc[ax]) return x;
        } else {
            if (u < g->fast) return x;
            if (u <= gauss_accept_prob(g->sigma, ax)) return x;
        }
    }
}

/* D_{Z,1} by an integer cumulative distribution table: CDT[x] = floor(2^64 Pr[|X| <= x]).
   Integer-only, hence identical on every platform (the randomness of a share is re-derived from
   its seed by the authority); statistical distance at most about 2^-60 per sample. */
static const uint64_t CDT1[9] = {7359186107371418017ULL, 16286330116675483855ULL, 18278245189139596012ULL,
                                 18441751535121735928ULL, 18446688998943340873ULL, 18446743849211842779ULL,
                                 18446744073372353484ULL, 18446744073709365182ULL, 18446744073709551578ULL};
static void gauss1_vec(prg *p, int64_t *out, size_t n) {
    uint8_t sb = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t r = prg_u64(p);
        int64_t x = 0;
        for (int j = 0; j < 9; j++) x += (r >= CDT1[j]);
        if ((i & 7) == 0) prg_bytes(p, &sb, 1);
        int neg = (sb >> (i & 7)) & 1;
        out[i] = (neg && x) ? -x : x;
    }
}

/* the Box-Muller spare is discarded at both ends, so that a vector depends only on the
   stream it is drawn from (needed when randomness is re-derived from a seed) */
void gauss_vec(gauss_sampler *g, prg *p, int64_t *out, size_t n) {
    if (g->sigma == 1.0) { gauss1_vec(p, out, n); return; }
    g->has_spare = 0;
    for (size_t i = 0; i < n; i++) out[i] = gauss_sample(g, p);
    g->has_spare = 0;
}

/* ---------------- challenges ---------------- */

void sample_challenge(challenge *c, prg *p) {      /* SampleInBall as in Dilithium */
    int8_t d[TV_N];
    memset(d, 0, sizeof d);
    uint8_t signs[8];
    prg_bytes(p, signs, 8);
    uint64_t sb = 0;
    for (int i = 7; i >= 0; i--) sb = (sb << 8) | signs[i];
    for (int i = TV_N - TV_CH_W; i < TV_N; i++) {
        unsigned j = prg_below(p, (unsigned)i + 1);
        d[i] = d[j];
        d[j] = (sb & 1) ? -1 : 1;
        sb >>= 1;
    }
    challenge_from_dense(c, d);
}

void challenge_dense(int8_t out[TV_N], const challenge *c) {
    memset(out, 0, TV_N);
    for (int k = 0; k < TV_CH_W; k++) out[c->pos[k]] = c->sgn[k];
}

int challenge_from_dense(challenge *c, const int8_t in[TV_N]) {
    int k = 0;
    for (int i = 0; i < TV_N; i++) {
        if (in[i] == 0) continue;
        if ((in[i] != 1 && in[i] != -1) || k == TV_CH_W) return -1;
        c->pos[k] = (uint8_t)i;
        c->sgn[k] = in[i];
        k++;
    }
    return k == TV_CH_W ? 0 : -1;
}

void sample_perm(perm_ch *pi, prg *p) {
    for (int i = 0; i < TV_N; i++) pi->s[i] = (uint8_t)i;
    for (int i = TV_N - 1; i > 0; i--) {
        unsigned j = prg_below(p, (unsigned)i + 1);
        uint8_t t = pi->s[i]; pi->s[i] = pi->s[j]; pi->s[j] = t;
    }
    uint8_t bits[8];
    prg_bytes(p, bits, 8);
    for (int i = 0; i < TV_CH_W; i++) pi->b[i] = (bits[i >> 3] >> (i & 7)) & 1;
}

/* g[s[p]] = f[p] * (-1)^{b[rank_f(p)]}, rank_f(p) = index of p among the nonzero positions of f */
void perm_apply(challenge *out, const perm_ch *pi, const challenge *f) {
    int8_t fd[TV_N], gd[TV_N];
    challenge_dense(fd, f);
    memset(gd, 0, sizeof gd);
    int rank = 0;
    for (int p = 0; p < TV_N; p++) {
        if (!fd[p]) continue;
        gd[pi->s[p]] = pi->b[rank++] ? (int8_t)-fd[p] : fd[p];
    }
    challenge_from_dense(out, gd);
}

void perm_apply_inv(challenge *out, const perm_ch *pi, const challenge *g) {
    int8_t gd[TV_N], fd[TV_N];
    challenge_dense(gd, g);
    memset(fd, 0, sizeof fd);
    int rank = 0;
    for (int p = 0; p < TV_N; p++) {
        int8_t v = gd[pi->s[p]];
        if (!v) continue;
        fd[p] = pi->b[rank++] ? (int8_t)-v : v;
    }
    challenge_from_dense(out, fd);
}

void ch_mul_int(int64_t *out, const challenge *c, const int64_t *v, size_t k) {
    for (size_t e = 0; e < k; e++) {
        const int64_t *a = v + e * TV_N;
        int64_t *r = out + e * TV_N;
        memset(r, 0, sizeof(int64_t) * TV_N);
        for (int t = 0; t < TV_CH_W; t++) {
            unsigned p = c->pos[t];
            int64_t s = c->sgn[t];
            for (unsigned i = 0; i < TV_N; i++) {
                unsigned idx = i + p;
                if (idx < TV_N) r[idx] += s * a[i];
                else r[idx - TV_N] -= s * a[i];
            }
        }
    }
}
