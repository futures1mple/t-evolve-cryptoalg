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

#include "cdt_tables.h"

/* |X| by a cumulative table T[0..len) of 192-bit entries: |X| = #{j : r >= T[j]} for r uniform in
   [0, 2^192). r is drawn lazily: first its top 16 bits u; the entries with top 16 bits < u are all
   <= r and those with top 16 bits > u are all > r, so only the entries whose top 16 bits equal u
   (usually none) need more bits of r. The stream consumption depends only on the stream itself. */
typedef struct { uint16_t *lt; } cdt_index;      /* lt[u] = #{j : top16(T[j]) < u}, u = 0..65536 */

static void cdt_index_build(cdt_index *ix, const uint64_t (*tab)[3], int len) {
    ix->lt = malloc(sizeof(uint16_t) * 65537);
    int j = 0;
    for (uint32_t u = 0; u <= 65536; u++) {
        while (j < len && (tab[j][0] >> 48) < u) j++;
        ix->lt[u] = (uint16_t)j;
    }
}

static uint64_t prg_bits48(prg *p) {
    uint8_t b[6];
    prg_bytes(p, b, 6);
    uint64_t x = 0;
    for (int i = 5; i >= 0; i--) x = (x << 8) | b[i];
    return x;
}

static int64_t cdt_abs(prg *p, const uint64_t (*tab)[3], const cdt_index *ix) {
    uint8_t b[2];
    prg_bytes(p, b, 2);
    uint32_t u = (uint32_t)b[0] | ((uint32_t)b[1] << 8);
    int lo = ix->lt[u], hi = ix->lt[u + 1];
    if (lo == hi) return lo;
    /* rare: complete r and compare with the entries in [lo, hi) */
    uint64_t w[3];
    w[0] = ((uint64_t)u << 48) | prg_bits48(p);
    int have = 1, x = lo;
    for (int j = lo; j < hi; j++) {
        int ge = 1;
        for (int k = 0; k < 3; k++) {
            if (have <= k) { w[k] = prg_u64(p); have = k + 1; }
            if (w[k] != tab[j][k]) { ge = w[k] > tab[j][k]; break; }
        }
        if (!ge) break;
        x = j + 1;
    }
    return x;
}

static cdt_index IX1, IX256;
static void cdt_ready(void) {
    if (!IX1.lt) cdt_index_build(&IX1, CDT_S1, CDT_S1_LEN);
    if (!IX256.lt) cdt_index_build(&IX256, CDT_S256, CDT_S256_LEN);
}

typedef struct { uint64_t bits; int left; } signbits;

static int64_t cdt_sample(prg *p, signbits *sb, const uint64_t (*tab)[3], const cdt_index *ix) {
    int64_t x = cdt_abs(p, tab, ix);
    if (x == 0) return 0;
    if (sb->left == 0) { sb->bits = prg_u64(p); sb->left = 64; }
    int neg = (int)(sb->bits & 1);
    sb->bits >>= 1; sb->left--;
    return neg ? -x : x;
}

/* eta = sqrt(ln(2 + 2^161)/pi) >= eta_eps(Z) for eps = 2^-160 (Gaussian parameter s = sqrt(2 pi) sigma).
   The convolution theorem needs s_in >= sqrt(2) max(a,b) eta, i.e. max(a,b)^2 <= (pi/eta^2) sigma_in^2.
   C2 = pi/eta^2 = pi^2 / ln(2 + 2^161), rounded down by a relative 1e-12 for safety. */
static double conv_c2(void) { return (9.8696044010893586188 / (161.0 * 0.69314718055994530942)) * (1.0 - 1e-12); }

static int64_t gcd64(int64_t a, int64_t b) { while (b) { int64_t t = a % b; a = b; b = t; } return a; }

/* the best (a,b) with gcd 1, b <= a <= zmax, a^2+b^2 >= need, minimizing a^2+b^2; returns 0 if none */
static int64_t best_pair(int64_t need, int64_t zmax, int32_t *pa, int32_t *pb) {
    int64_t best = 0;
    if (need < 2) need = 2;                      /* (1,1) is the smallest pair */
    int64_t a0 = (int64_t)floor(sqrt((double)need / 2.0));
    if (a0 < 1) a0 = 1;
    while (a0 > 1 && 2 * (a0 - 1) * (a0 - 1) >= need) a0--;
    for (int64_t a = a0; a <= zmax; a++) {
        if (best && a * a >= best) break;        /* a^2 + b^2 > a^2 >= best */
        int64_t rem = need - a * a, b;
        if (rem <= 1) b = 1;
        else {
            b = (int64_t)floor(sqrt((double)rem));
            while (b * b < rem) b++;
            while (b > 1 && (b - 1) * (b - 1) >= rem) b--;
        }
        for (; b <= a; b++) if (gcd64(a, b) == 1) break;
        if (b > a) continue;
        int64_t v = a * a + b * b;
        if (!best || v < best) { best = v; *pa = (int32_t)a; *pb = (int32_t)b; }
    }
    return best;
}

typedef struct { double target; gauss_sampler g; } plan_cache;
static plan_cache CACHE[16];
static int CACHE_N = 0;

static int plan(gauss_sampler *g, double target) {
    for (int i = 0; i < CACHE_N; i++) if (CACHE[i].target == target) { *g = CACHE[i].g; return 0; }
    memset(g, 0, sizeof *g);
    const double s0 = CDT_BASE_SIGMA, c2 = conv_c2();
    if (target <= s0) { g->levels = 0; g->sigma = s0; g->sigma2 = (uint64_t)(s0 * s0); goto done; }
    /* n_total = (a1^2+b1^2)(a2^2+b2^2) must be >= need = ceil((target/s0)^2) */
    double r = target / s0;
    int64_t need = (int64_t)ceil(r * r);
    int64_t zmax0 = (int64_t)floor(sqrt(c2 * s0 * s0));
    int32_t a, b;
    int64_t best = 0;
    int64_t n1 = best_pair(need, zmax0, &a, &b);
    if (n1) { best = n1; g->levels = 1; g->a[0] = a; g->b[0] = b; }
    for (int64_t a1 = 1; a1 <= zmax0; a1++)
        for (int64_t b1 = 1; b1 <= a1; b1++) {
            if (gcd64(a1, b1) != 1) continue;
            int64_t m1 = a1 * a1 + b1 * b1;
            int64_t zmax1 = (int64_t)floor(sqrt(c2 * s0 * s0 * (double)m1));
            int64_t need2 = (need + m1 - 1) / m1;
            if (best && m1 * 2 >= best) continue;          /* level 2 adds a factor >= 2 */
            int32_t a2, b2;
            int64_t m2 = best_pair(need2, zmax1, &a2, &b2);
            if (!m2) continue;
            if (!best || m1 * m2 < best) {
                best = m1 * m2; g->levels = 2;
                g->a[0] = (int32_t)a1; g->b[0] = (int32_t)b1; g->a[1] = a2; g->b[1] = b2;
            }
        }
    if (!best) return -1;
    g->sigma = s0 * sqrt((double)best);
    g->sigma2 = (uint64_t)(s0 * s0) * (uint64_t)best;
done:
    if (CACHE_N < 16) { CACHE[CACHE_N].target = target; CACHE[CACHE_N].g = *g; CACHE_N++; }
    return 0;
}

double gauss_round_sigma(double sigma) {
    if (sigma == 1.0) return 1.0;
    gauss_sampler g;
    if (plan(&g, sigma)) return -1.0;
    return g.sigma;
}

int gauss_init(gauss_sampler *g, double sigma) {
    cdt_ready();
    memset(g, 0, sizeof *g);
    if (sigma == 1.0) { g->sigma = 1.0; g->sigma2 = 1; g->levels = -1; return 0; }
    if (!(sigma > 1.0)) return -1;
    return plan(g, sigma);
}

void gauss_free(gauss_sampler *g) { (void)g; }

uint64_t gauss_sigma2(double sigma) {
    gauss_sampler g;
    if (gauss_init(&g, sigma)) return 0;
    return g.sigma2;
}

/* ---------------- exact rejection ---------------- */

typedef unsigned __int128 u128;

/* uniform integer in [0, m), 0 < m < 2^127, by rejection from the smallest power of two >= m */
static u128 uniform_below(prg *p, u128 m) {
    int bits = 0;
    while (bits < 128 && ((u128)1 << bits) < m) bits++;
    u128 mask = bits >= 128 ? ~(u128)0 : (((u128)1 << bits) - 1);
    for (;;) {
        u128 x = ((u128)prg_u64(p) << 64) | prg_u64(p);
        x &= mask;
        if (x < m) return x;
    }
}

/* 1 with probability num/den, 0 <= num <= den */
static int bern_frac(prg *p, u128 num, u128 den) { return uniform_below(p, den) < num; }

/* 1 with probability exp(-num/den), 0 <= num <= den [CKS20, Alg. 1] */
static int bern_exp01(prg *p, u128 num, u128 den) {
    u128 K = 1;
    while (bern_frac(p, num, den * K)) K++;
    return (int)(K & 1);
}

int bern_exp(prg *p, i128 num, i128 den) {
    if (num <= 0) return 1;
    u128 n = (u128)num, d = (u128)den;
    while (n > d) {                       /* exp(-g) = exp(-1) * exp(-(g-1)) */
        if (!bern_exp01(p, d, d)) return 0;
        n -= d;
    }
    return bern_exp01(p, n, d);
}

int reject_accept(prg *p, i128 zs, i128 ss, uint64_t sigma2, int64_t cn, int64_t cd) {
    /* R = (2 zs - ss)/(2 sigma2) + cn/cd = ((2 zs - ss) cd + 2 sigma2 cn) / (2 sigma2 cd); accept w.p. exp(-R) */
    i128 den = (i128)2 * (i128)sigma2 * (i128)cd;
    i128 num = ((i128)2 * zs - ss) * (i128)cd + (i128)2 * (i128)sigma2 * (i128)cn;
    return bern_exp(p, num, den);
}

static int64_t conv_sample(const gauss_sampler *g, prg *p, signbits *sb, int level) {
    if (level == 0) return cdt_sample(p, sb, CDT_S256, &IX256);
    int64_t x = conv_sample(g, p, sb, level - 1), y = conv_sample(g, p, sb, level - 1);
    return (int64_t)g->a[level - 1] * x + (int64_t)g->b[level - 1] * y;
}

int64_t gauss_sample(gauss_sampler *g, prg *p) {
    signbits sb = {0, 0};
    if (g->levels < 0) return cdt_sample(p, &sb, CDT_S1, &IX1);
    return conv_sample(g, p, &sb, g->levels);
}

void gauss_vec(gauss_sampler *g, prg *p, int64_t *out, size_t n) {
    signbits sb = {0, 0};
    if (g->levels < 0) { for (size_t i = 0; i < n; i++) out[i] = cdt_sample(p, &sb, CDT_S1, &IX1); return; }
    for (size_t i = 0; i < n; i++) out[i] = conv_sample(g, p, &sb, g->levels);
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
