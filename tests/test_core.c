/* Tests of the core: SHAKE known answers, ring arithmetic against a schoolbook reference,
 * challenge permutations, and the discrete Gaussian sampler. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "keccak.h"
#include "ring.h"
#include "sample.h"

static int fails = 0;
#define CHECK(cond, ...) do { if (!(cond)) { fails++; printf("FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

static void hex(char *out, const uint8_t *b, size_t n) { for (size_t i = 0; i < n; i++) sprintf(out + 2 * i, "%02x", b[i]); }

static void test_shake(void) {
    struct { int f; int len; const char *head, *tail; } kat[] = {
        {128, 0, "7f9c2ba4e88f827d616045507605853ed73b8093f6efbc88eb1a6eacfa66ef263cb1eea988004b93", "3ea3aeeb613a7f1b1de33fd75081f592305f2e45"},
        {128, 3, "5881092dd818bf5cf8a3ddb793fbcba74097d5c526a6d35f97b83351940f2cc844c50af32acd3f2c", "6bdb2e06a3eed543a38919b57ecbec737f4086be"},
        {128, 200, "0c4234ca1e31801ae606f8b8d8e0665c66f42a21d601c2681858a92c79ad5d69e143c3b1393dd894", "2e126fbe086095b9433e06a84f609a0c91793cc7"},
        {128, 1000, "36514c827683dd1b85c3304d8079021cdcce289a0257b8acf0317f7d8ce455d4b6ca10bddb068b43", "e8d9f0a09e6257b7a9f201ae6d8165de3750514c"},
        {256, 0, "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762fd75dc4ddd8c0f200", "ed96d477ff96390bf9a66d1368b208e21f7c10d0"},
        {256, 3, "483366601360a8771c6863080cc4114d8db44530f8f1e1ee4f94ea37e78b5739d5a15bef186a5386", "66caa7d8ddcbec7da52b42215c11d5f8ee57f341"},
        {256, 200, "4ee1ca03272b05d3bfb1e1c79a967f823b9fc5e4bb3987b1ba9e9cb5afb07a5ee3a07fbd457a9436", "01358ae64f3f0ccedfa05b29e84e1a11a635bfe7"},
        {256, 1000, "ea7e0f6345ce9236910251b203e1b228aa83bf8e2ea18178c1f1e9dc083524040c47d6476433c448", "514951c70d38d05890e828633d6f19f872cbafad"},
    };
    uint8_t msg[1000], out[300], out2[300];
    char h[700];
    for (size_t k = 0; k < sizeof kat / sizeof kat[0]; k++) {
        int len = kat[k].len;
        if (len == 3) memcpy(msg, "abc", 3);
        else if (len == 200) for (int i = 0; i < 200; i++) msg[i] = (uint8_t)i;
        else for (int i = 0; i < len; i++) msg[i] = (uint8_t)((i * 7) % 256);
        if (kat[k].f == 128) shake128(out, 300, msg, (size_t)len); else shake256(out, 300, msg, (size_t)len);
        hex(h, out, 300);
        CHECK(strncmp(h, kat[k].head, 80) == 0 && strcmp(h + 600 - 40, kat[k].tail) == 0, "SHAKE%d KAT len %d", kat[k].f, len);
        /* incremental absorb in odd chunks and squeeze in odd chunks */
        keccak_state st;
        if (kat[k].f == 128) shake128_init(&st); else shake256_init(&st);
        for (int i = 0; i < len;) { int c = 1 + (i * 13) % 97; if (i + c > len) c = len - i; keccak_absorb(&st, msg + i, (size_t)c); i += c; }
        keccak_finalize(&st);
        for (int i = 0; i < 300;) { int c = 1 + (i * 5) % 71; if (i + c > 300) c = 300 - i; keccak_squeeze(&st, out2 + i, (size_t)c); i += c; }
        CHECK(memcmp(out, out2, 300) == 0, "SHAKE%d incremental len %d", kat[k].f, len);
    }
}

static void test_ring(uint64_t q) {
    CHECK(ring_init(q) == 0, "ring_init(%llu)", (unsigned long long)q);
    uint8_t seed[32] = {1, 2, 3};
    prg p;
    prg_init(&p, 0, seed, 32);
    for (int t = 0; t < 20; t++) {
        poly a, b, r1, r2, bm;
        sample_uniform_poly(&a, &p);
        sample_uniform_poly(&b, &p);
        poly_mul_ref(&r1, &a, &b);
        poly an = a, bn = b;
        poly_ntt(&an);
        poly_ntt(&bn);
        bm = bn;
        poly_to_mont(&bm);
        poly_basemul_mont(&r2, &an, &bm);
        poly_invntt(&r2);
        CHECK(memcmp(&r1, &r2, sizeof r1) == 0, "NTT product q=%llu", (unsigned long long)q);
        poly back = an;
        poly_invntt(&back);
        CHECK(memcmp(&back, &a, sizeof a) == 0, "NTT roundtrip");
    }
    /* dot product of 30 terms with maximal coefficients q-1 (worst case for the accumulator) */
    {
        enum { L = 30 };
        poly A[L], y[L], r, ref, t;
        for (int j = 0; j < L; j++) for (int i = 0; i < TV_N; i++) { A[j].c[i] = q - 1 - (uint64_t)((i + j) % 3); y[j].c[i] = q - 1 - (uint64_t)((i * j) % 5); }
        poly_zero(&ref);
        for (int j = 0; j < L; j++) { poly_mul_ref(&t, &A[j], &y[j]); poly_add(&ref, &ref, &t); }
        for (int j = 0; j < L; j++) { poly_ntt(&A[j]); poly_to_mont(&A[j]); poly_ntt(&y[j]); }
        poly_dot_mont(&r, A, y, L);
        poly_invntt(&r);
        CHECK(memcmp(&r, &ref, sizeof r) == 0, "dot product worst case");
    }
    /* sparse challenge product */
    for (int t = 0; t < 10; t++) {
        challenge c;
        sample_challenge(&c, &p);
        int8_t cd[TV_N];
        challenge_dense(cd, &c);
        poly cp, a, r1, r2;
        for (int i = 0; i < TV_N; i++) cp.c[i] = mod_from_int(cd[i]);
        sample_uniform_poly(&a, &p);
        poly_mul_ref(&r1, &cp, &a);
        poly_mul_sparse(&r2, &a, c.pos, c.sgn, TV_CH_W);
        CHECK(memcmp(&r1, &r2, sizeof r1) == 0, "sparse product");
        int64_t v[TV_N], w[TV_N];
        for (int i = 0; i < TV_N; i++) v[i] = (int64_t)(prg_u64(&p) % 2001) - 1000;
        ch_mul_int(w, &c, v, 1);
        poly vp, wp;
        poly_from_int(&vp, v);
        poly_mul_ref(&wp, &cp, &vp);
        int ok = 1;
        for (int i = 0; i < TV_N; i++) ok &= (mod_from_int(w[i]) == wp.c[i]);
        CHECK(ok, "integer challenge product");
    }
}

static void test_challenges(void) {
    uint8_t seed[32] = {9};
    prg p;
    prg_init(&p, 1, seed, 32);
    for (int t = 0; t < 200; t++) {
        challenge f, g, h;
        perm_ch pi;
        sample_challenge(&f, &p);
        int8_t d[TV_N];
        challenge_dense(d, &f);
        int nz = 0;
        for (int i = 0; i < TV_N; i++) nz += d[i] != 0;
        CHECK(nz == TV_CH_W, "challenge weight");
        sample_perm(&pi, &p);
        perm_apply(&g, &pi, &f);
        perm_apply_inv(&h, &pi, &g);
        CHECK(memcmp(&f, &h, sizeof f) == 0, "perm inverse");
    }
}

/* probability of [lo, hi] under D_{Z,sigma}; exact sums for small sigma, erf otherwise */
static double dg_mass(double sigma, double Z, long lo, long hi) {
    if (sigma <= 3000) {
        double m = 0;
        for (long x = lo; x <= hi; x++) m += exp(-(double)x * x / (2 * sigma * sigma));
        return m / Z;
    }
    double c = sqrt(2.0) * sigma;
    return 0.5 * (erf((hi + 0.5) / c) - erf((lo - 0.5) / c));
}

/* chi-square test of gauss_vec against D_{Z,sigma_out} with 2B+1 bins of width w plus two tails */
static void test_gauss(double target, int n) {
    gauss_sampler g;
    CHECK(gauss_init(&g, target) == 0, "gauss_init %g", target);
    double sigma = g.sigma;
    CHECK(sigma >= target && (target < 20000 || sigma <= target * 1.001), "sigma_out %g for target %g", sigma, target);
    uint8_t seed[32] = {7};
    prg p;
    prg_init(&p, 2, seed, 32);
    long w = sigma < 8 ? 1 : (long)(sigma / 4);
    int B = (int)ceil(4 * sigma / w);
    long *cnt = calloc((size_t)(2 * B + 3), sizeof(long));
    double s1 = 0, s2 = 0;
    int64_t *vec = malloc(sizeof(int64_t) * 1000);
    for (int i = 0; i < n; i++) {
        if (i % 1000 == 0) gauss_vec(&g, &p, vec, 1000);
        int64_t x = vec[i % 1000];
        s1 += (double)x; s2 += (double)x * (double)x;
        long t = (long)x + w / 2, q = t >= 0 ? t / w : -((-t + w - 1) / w);   /* floor((x + w/2)/w) */
        int bin = q < -B ? 0 : q > B ? 2 * B + 2 : (int)(q + B + 1);
        cnt[bin]++;
    }
    double Z = 0;
    if (sigma <= 3000) for (long x = -(long)(40 * sigma); x <= (long)(40 * sigma); x++) Z += exp(-(double)x * x / (2 * sigma * sigma));
    double chi = 0, inner = 0;
    int dof = 0;
    for (int q = -B; q <= B; q++) {
        long lo = q * w - w / 2, hi = lo + w - 1;
        if (w == 1) lo = hi = q;
        double e = n * dg_mass(sigma, Z, lo, hi);
        inner += e;
        double o = (double)cnt[q + B + 1];
        if (e > 5) { chi += (o - e) * (o - e) / e; dof++; }
    }
    double et = (n - inner) / 2;
    if (et > 5) { chi += (cnt[0] - et) * (cnt[0] - et) / et + (cnt[2 * B + 2] - et) * (cnt[2 * B + 2] - et) / et; dof += 2; }
    dof -= 1;
    double var = s2 / n - (s1 / n) * (s1 / n);
    printf("  gauss target=%g sigma_out=%.3f levels=%d: n=%d var/sigma^2=%.4f chi2=%.1f dof=%d\n",
           target, sigma, g.levels, n, var / (sigma * sigma), chi, dof);
    CHECK(chi < dof + 5 * sqrt(2.0 * dof), "gauss chi-square sigma=%g", sigma);
    CHECK(fabs(var / (sigma * sigma) - 1) < 6 * sqrt(2.0 / n) + (sigma < 3 ? 0.01 : 0), "gauss variance sigma=%g", sigma);
    free(cnt);
    free(vec);
    gauss_free(&g);
}

/* the conditions of the convolution theorem hold for every plan used by the parameter sets */
static void test_gauss_plans(void) {
    const double c2 = 9.8696044010893586188 / log(2.0 + pow(2.0, 221));
    double targets[] = {300, 1000, 25113.9, 27000, 54000, 1.2e5, 2.8e6, 8.9e6, 3.0e7, 9.1e7, 1.8e8};
    for (unsigned i = 0; i < sizeof targets / sizeof *targets; i++) {
        gauss_sampler g;
        CHECK(gauss_init(&g, targets[i]) == 0, "plan %g", targets[i]);
        double s_in = 256;
        double prod = 1;
        for (int l = 0; l < g.levels; l++) {
            long a = g.a[l], b = g.b[l], m = a > b ? a : b, x = a, y = b;
            while (y) { long t = x % y; x = y; y = t; }
            CHECK(x == 1, "gcd level %d target %g", l, targets[i]);
            CHECK((double)m * m <= c2 * s_in * s_in, "smoothing condition level %d target %g", l, targets[i]);
            s_in *= sqrt((double)(a * a + b * b));
            prod *= (double)(a * a + b * b);
        }
        CHECK(fabs(s_in - g.sigma) < 1e-6 * g.sigma && g.sigma >= targets[i], "sigma_out target %g", targets[i]);
        printf("  plan target=%g: levels=%d a=(%d,%d) b=(%d,%d) sigma_out/target=%.5f\n", targets[i], g.levels,
               g.a[0], g.a[1], g.b[0], g.b[1], g.sigma / targets[i]);
    }
}

/* exact Bernoulli(exp(-num/den)) and the rejection rule: frequencies against exp(-g) (binomial, 5 sd) */
static void test_bern_exp(void) {
    uint8_t seed[32] = {5};
    prg p;
    prg_init(&p, 4, seed, 32);
    const long long nd[][2] = {{0, 1}, {3, 10}, {1, 1}, {27, 10}, {5, 1}, {1234567, 1000000}};
    for (unsigned i = 0; i < sizeof nd / sizeof *nd; i++) {
        long n = 400000, c = 0;
        for (long j = 0; j < n; j++) c += bern_exp(&p, nd[i][0], nd[i][1]);
        double pe = exp(-(double)nd[i][0] / nd[i][1]), sd = sqrt(n * pe * (1 - pe)) + 1e-9;
        printf("  bern_exp g=%g: %.5f (exp(-g) = %.5f)\n", (double)nd[i][0] / nd[i][1], (double)c / n, pe);
        CHECK(fabs(c - n * pe) <= 5 * sd + 1, "bern_exp frequency g=%g", (double)nd[i][0] / nd[i][1]);
    }
    /* reject_accept: accept w.p. min(1, exp((ss - 2 zs)/(2 s2) - cn/cd)) */
    long n = 400000, c = 0;
    for (long j = 0; j < n; j++) c += reject_accept(&p, 30, 20, 50, 265, 242);   /* (20-60)/100 - 265/242 */
    double pe = exp(-0.4 - 265.0 / 242.0), sd = sqrt(n * pe * (1 - pe));
    printf("  reject_accept: %.5f (expected %.5f)\n", (double)c / n, pe);
    CHECK(fabs(c - n * pe) <= 5 * sd + 1, "reject_accept frequency");
    c = 0;
    for (long j = 0; j < 1000; j++) c += reject_accept(&p, -1000, 0, 50, 265, 242);   /* R < 0: always accept */
    CHECK(c == 1000, "reject_accept with R < 0");
    /* saturating sums of squares (verifier side) never wrap around */
    i128 sq = 0;
    for (int j = 0; j < 4; j++) sq = sat_madd(sq, INT64_MAX, INT64_MAX);
    CHECK(sq == I128_MAX, "sat_madd saturates");
    CHECK(chk_madd(5, -3, 4, "test") == -7, "chk_madd");
    /* a large exponent: exp(-40.5) ~ 2.6e-18, so no acceptance in 10^5 trials, and many exp(-1) factors */
    c = 0;
    for (long j = 0; j < 100000; j++) c += bern_exp(&p, 81, 2);
    CHECK(c == 0, "bern_exp with gamma = 40.5");
}

/* known answer for seed-derived commitment randomness: must be identical on every platform,
   because the authority re-derives the voter's randomness from the seed */
static void test_seed_kat(void) {
    uint8_t s[32];
    for (int i = 0; i < 32; i++) s[i] = (uint8_t)i;
    prg p;
    prg_init(&p, 0x20, s, 32);
    gauss_sampler g;
    gauss_init(&g, 1.0);
    int64_t v[4352];
    gauss_vec(&g, &p, v, 4352);
    uint8_t b[4352], h[16];
    for (int i = 0; i < 4352; i++) b[i] = (uint8_t)v[i];
    shake256(h, 16, b, 4352);
    char hx[40];
    hex(hx, h, 16);
    CHECK(strcmp(hx, "117ced2534df93428c0c0c3f826bf947") == 0, "seed-derived randomness KAT (%s)", hx);
    gauss_free(&g);
}

int main(void) {
    test_shake();
    test_ring(2199023255633ULL);        /* 2^41 + 81 */
    test_ring(70368744177937ULL);       /* 2^46 + 273 */
    test_ring(562949953422097ULL);      /* 2^49 + 785 */
    test_challenges();
    test_seed_kat();
    test_bern_exp();
    test_gauss_plans();
    test_gauss(1.0, 4000000);
    test_gauss(256.0, 2000000);
    test_gauss(1000.0, 2000000);
    test_gauss(25113.9, 2000000);
    test_gauss(8.9e6, 2000000);
    printf("%s (%d failures)\n", fails ? "FAILED" : "all core tests passed", fails);
    return fails != 0;
}
