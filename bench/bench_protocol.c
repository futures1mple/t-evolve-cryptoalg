/* Benchmarks of T-EVOLVE and EVOLVE on the same core (single thread).
 *
 *   bench_protocol [-nv N_V] [-reps R] [-seed S] [-d d] [-q q] [-out file.csv] [-skip-agg]
 *
 * Measures, with medians and interquartile ranges over R repetitions:
 *   micro      NTT, product in R_q, matrix-vector product, Gaussian samplers, SHAKE256;
 *   ballot     Vote (with all attempts), Verify, Check_k per share, sizes after encoding,
 *              for T-EVOLVE (yes/no; one of two candidates) and for EVOLVE (yes/no);
 *   aggregation  the complete tally of N_V ballots by all n authorities: round 1 (tree and kappa
 *              committed attempts), round 2 (search for a non-aborting attempt), repeated rounds
 *              until t authorities succeed, VerAgg, Combine; per-block costs for extrapolation.
 * Every row of the CSV carries the machine, compiler, flags, commit and seed. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "agg.h"
#include "codec.h"
#include "evolve.h"
#include "shamir.h"
#include "tally.h"
#include "timer.h"

#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif
#ifndef BUILD_FLAGS
#define BUILD_FLAGS "unknown"
#endif

static char cpu_name[256] = "unknown";
static FILE *csv;
static unsigned long long g_seed = 1;

static void detect_cpu(void) {
#ifdef _WIN32
    const char *e = getenv("PROCESSOR_IDENTIFIER");
    if (e) snprintf(cpu_name, sizeof cpu_name, "%s", e);
#else
    FILE *f = fopen("/proc/cpuinfo", "r");
    if (!f) return;
    char line[512];
    while (fgets(line, sizeof line, f))
        if (!strncmp(line, "model name", 10)) {
            char *c = strchr(line, ':');
            if (c) { snprintf(cpu_name, sizeof cpu_name, "%s", c + 2); cpu_name[strcspn(cpu_name, "\n")] = 0; }
            break;
        }
    fclose(f);
#endif
    for (char *p = cpu_name; *p; p++) if (*p == ',') *p = ';';
}

static int cmpd(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return x < y ? -1 : x > y; }
typedef struct { double med, q1, q3, min; } st4;
static st4 stats(double *v, int n) {
    qsort(v, (size_t)n, sizeof(double), cmpd);
    st4 s = {v[n / 2], v[n / 4], v[(3 * n) / 4], v[0]};
    return s;
}

static void row(const char *group, const char *name, const char *config, double med, double q1, double q3, double mn, int reps, const char *unit) {
    printf("%-12s %-34s %-28s %14.4f %-6s (IQR %.4f-%.4f, n=%d)\n", group, name, config, med, unit, q1, q3, reps);
    fprintf(csv, "%s,%s,%s,%.6f,%.6f,%.6f,%.6f,%d,%s,%s,%s,%s,%s,%llu\n", group, name, config, med, q1, q3, mn, reps, unit,
            cpu_name, __VERSION__, BUILD_FLAGS, GIT_COMMIT, g_seed);
    fflush(csv);
}
static void row1(const char *group, const char *name, const char *config, double v, const char *unit) { row(group, name, config, v, v, v, v, 1, unit); }

static void seed_from(uint8_t s[32], unsigned long long a, unsigned long long b) {
    memset(s, 0, 32);
    for (int i = 0; i < 8; i++) { s[i] = (uint8_t)(a >> (8 * i)); s[8 + i] = (uint8_t)(b >> (8 * i)); }
}

/* ------------------------------------------------------------------ micro benchmarks */
static void bench_micro(const tv_pub *pub, int reps) {
    const tv_params *p = &pub->prm;
    char cfg[64];
    snprintf(cfg, sizeof cfg, "log2q=%.1f d=%d mu=%d", log2((double)RING.q), p->d, p->mu);
    double *v = malloc(sizeof(double) * (size_t)reps);
    uint8_t s[32];
    seed_from(s, g_seed, 11);
    prg g;
    prg_init(&g, 5, s, 32);
    poly a, b, r;
    sample_uniform_poly(&a, &g);
    sample_uniform_poly(&b, &g);
    const int inner = 200;
    for (int i = 0; i < reps; i++) { double t0 = now_ms(); for (int j = 0; j < inner; j++) poly_ntt(&a); v[i] = (now_ms() - t0) * 1000 / inner; }
    st4 x = stats(v, reps); row("micro", "ntt_forward", cfg, x.med, x.q1, x.q3, x.min, reps, "us");
    poly_ntt(&b); poly_to_mont(&b);
    for (int i = 0; i < reps; i++) { double t0 = now_ms(); for (int j = 0; j < inner; j++) poly_basemul_mont(&r, &a, &b); v[i] = (now_ms() - t0) * 1000 / inner; }
    x = stats(v, reps); row("micro", "product_ntt_domain", cfg, x.med, x.q1, x.q3, x.min, reps, "us");
    int64_t *vec = malloc(sizeof(int64_t) * (size_t)p->mu * TV_N);
    gauss_sampler g1, gJ;
    gauss_init(&g1, p->sigma);
    gauss_init(&gJ, p->sigma_J);
    gauss_vec(&g1, &g, vec, (size_t)p->mu * TV_N);
    poly *out = malloc(sizeof(poly) * p->rows);
    for (int i = 0; i < reps; i++) { double t0 = now_ms(); tv_commit(out, pub, out + p->d, vec); v[i] = (now_ms() - t0) * 1000; }
    x = stats(v, reps); row("micro", "commitment (C r + m)", cfg, x.med, x.q1, x.q3, x.min, reps, "us");
    for (int i = 0; i < reps; i++) { double t0 = now_ms(); gauss_vec(&g1, &g, vec, (size_t)p->mu * TV_N); v[i] = (now_ms() - t0) * 1000; }
    x = stats(v, reps); row("micro", "gauss_vector sigma=1 (mu*N)", cfg, x.med, x.q1, x.q3, x.min, reps, "us");
    for (int i = 0; i < reps; i++) { double t0 = now_ms(); gauss_vec(&gJ, &g, vec, (size_t)p->mu * TV_N); v[i] = (now_ms() - t0) * 1000; }
    x = stats(v, reps); row("micro", "gauss_vector sigma_J (mu*N)", cfg, x.med, x.q1, x.q3, x.min, reps, "us");
    uint8_t buf[1024], h[32];
    memset(buf, 1, sizeof buf);
    for (int i = 0; i < reps; i++) { double t0 = now_ms(); for (int j = 0; j < inner; j++) shake256(h, 32, buf, sizeof buf); v[i] = (now_ms() - t0) * 1000 / inner; }
    x = stats(v, reps); row("micro", "shake256_1KiB", cfg, x.med, x.q1, x.q3, x.min, reps, "us");
    gauss_free(&g1); gauss_free(&gJ);
    free(vec); free(out); free(v);
}

/* ------------------------------------------------------------------ ballots */
static void bench_tballot(int n, int t, int L, int w, int d, uint64_t q, int reps) {
    tv_params prm;
    tv_params_init(&prm, n, t, L, w, d, q);
    uint8_t s[32];
    seed_from(s, g_seed, 21);
    tv_pub pub;
    tv_pub_init(&pub, &prm, s);
    tv_akey *keys = calloc((size_t)n, sizeof(tv_akey));
    for (int k = 0; k < n; k++) seed_from(keys[k].key, g_seed, 100 + (unsigned)k);
    char cfg[96];
    snprintf(cfg, sizeof cfg, "T-EVOLVE n=%d t=%d L=%d w=%d d=%d log2q=%.1f", n, t, L, w, d, log2((double)RING.q));
    double *tv = malloc(sizeof(double) * reps), *tvf = malloc(sizeof(double) * reps), *tc = malloc(sizeof(double) * reps);
    double *att = malloc(sizeof(double) * reps), *by = malloc(sizeof(double) * reps), *bp = malloc(sizeof(double) * reps), *bcm = malloc(sizeof(double) * reps);
    tv_ballot b, b2;
    tv_voter_secret sec;
    tv_ballot_alloc(&b, &prm); tv_ballot_alloc(&b2, &prm); tv_secret_alloc(&sec, &prm);
    uint8_t *buf = malloc(tv_ballot_maxbytes(&prm));
    poly *mk = calloc((size_t)L, sizeof(poly));
    int64_t *rk = malloc(sizeof(int64_t) * (size_t)prm.mu * TV_N);
    int allok = 1;
    for (int i = 0; i < reps; i++) {
        int v[16] = {0};
        if (w == TV_W_FREE) v[0] = i & 1; else v[i % L] = 1;
        seed_from(s, g_seed, 1000 + (unsigned)i);
        tv_prove_stats st;
        double t0 = now_ms();
        tv_vote(&b, &sec, &st, &pub, keys, (uint64_t)i, v, s);
        tv[i] = now_ms() - t0;
        tv_sizes sz;
        size_t len = tv_ballot_encode(buf, tv_ballot_maxbytes(&prm), &pub, &b, &sz);
        t0 = now_ms();
        allok &= tv_ballot_decode(&b2, &pub, buf, len) == 0 && tv_verify_ballot(&pub, &b2);
        tvf[i] = now_ms() - t0;
        t0 = now_ms();
        allok &= tv_check_share(&pub, &keys[0], 1, &b2, mk, rk);
        tc[i] = now_ms() - t0;
        att[i] = st.attempts; by[i] = (double)len; bp[i] = (double)sz.proof; bcm[i] = (double)sz.commitments;
    }
    st4 x;
    x = stats(tv, reps); row("ballot", "vote (incl. all attempts)", cfg, x.med, x.q1, x.q3, x.min, reps, "ms");
    x = stats(tvf, reps); row("ballot", "decode+verify", cfg, x.med, x.q1, x.q3, x.min, reps, "ms");
    x = stats(tc, reps); row("ballot", "check_share (one authority)", cfg, x.med, x.q1, x.q3, x.min, reps, "ms");
    double ma = 0; for (int i = 0; i < reps; i++) ma += att[i];
    row1("ballot", "attempts (mean)", cfg, ma / reps, "count");
    x = stats(by, reps); row("ballot", "ballot size", cfg, x.med / 1024, x.q1 / 1024, x.q3 / 1024, x.min / 1024, reps, "KiB");
    x = stats(bp, reps); row("ballot", "  of which proof", cfg, x.med / 1024, x.q1 / 1024, x.q3 / 1024, x.min / 1024, reps, "KiB");
    x = stats(bcm, reps); row("ballot", "  of which commitments", cfg, x.med / 1024, x.q1 / 1024, x.q3 / 1024, x.min / 1024, reps, "KiB");
    row1("ballot", "all verified", cfg, allok, "bool");
    tv_ballot_free(&b); tv_ballot_free(&b2); tv_secret_free(&sec);
    free(buf); free(mk); free(rk); free(tv); free(tvf); free(tc); free(att); free(by); free(bp); free(bcm); free(keys);
    tv_pub_free(&pub);
}

static void bench_eballot(int n, int d, uint64_t q, int reps) {
    tv_params prm;
    tv_params_init(&prm, n, n, 1, TV_W_FREE, d, q);
    uint8_t s[32];
    seed_from(s, g_seed, 31);
    tv_pub pub;
    tv_pub_init(&pub, &prm, s);
    ev_params ep;
    ev_params_init(&ep, &prm);
    tv_akey *keys = calloc((size_t)n, sizeof(tv_akey));
    for (int k = 0; k < n; k++) seed_from(keys[k].key, g_seed, 200 + (unsigned)k);
    char cfg[96];
    snprintf(cfg, sizeof cfg, "EVOLVE n=%d L=1 d=%d log2q=%.1f", n, d, log2((double)RING.q));
    double *tv = malloc(sizeof(double) * reps), *tvf = malloc(sizeof(double) * reps), *att = malloc(sizeof(double) * reps), *by = malloc(sizeof(double) * reps);
    ev_ballot b;
    ev_ballot_alloc(&b, &prm);
    poly *sh = malloc(sizeof(poly) * n);
    int64_t *rnd = malloc(sizeof(int64_t) * (size_t)n * prm.mu * TV_N);
    int allok = 1;
    for (int i = 0; i < reps; i++) {
        seed_from(s, g_seed, 2000 + (unsigned)i);
        int a;
        double t0 = now_ms();
        ev_vote(&b, sh, rnd, &a, &pub, &ep, keys, (uint64_t)i, i & 1, s);
        tv[i] = now_ms() - t0;
        t0 = now_ms();
        allok &= ev_verify(&pub, &ep, &b);
        tvf[i] = now_ms() - t0;
        att[i] = a;
        by[i] = (double)ev_ballot_bytes(&pub, &ep, &b, NULL);
    }
    st4 x;
    x = stats(tv, reps); row("ballot", "vote (incl. all attempts)", cfg, x.med, x.q1, x.q3, x.min, reps, "ms");
    x = stats(tvf, reps); row("ballot", "verify", cfg, x.med, x.q1, x.q3, x.min, reps, "ms");
    double ma = 0; for (int i = 0; i < reps; i++) ma += att[i];
    row1("ballot", "attempts (mean)", cfg, ma / reps, "count");
    x = stats(by, reps); row("ballot", "ballot size", cfg, x.med / 1024, x.q1 / 1024, x.q3 / 1024, x.min / 1024, reps, "KiB");
    row1("ballot", "all verified", cfg, allok, "bool");
    ev_ballot_free(&b); free(sh); free(rnd); free(tv); free(tvf); free(att); free(by); free(keys);
    tv_pub_free(&pub);
}

/* ------------------------------------------------------------------ aggregation */
static void bench_agg(int NV, int n, int t, int d, uint64_t q, c2_mode mode) {
    tv_params prm;
    tv_params_init(&prm, n, t, 1, TV_W_FREE, d, q);
    uint8_t s[32];
    seed_from(s, g_seed, 41);
    tv_pub pub;
    tv_pub_init(&pub, &prm, s);
    char cfg[112];
    snprintf(cfg, sizeof cfg, "NV=%d n=%d t=%d L=1 d=%d log2q=%.1f C2=%s", NV, n, t, d, log2((double)RING.q), mode == C2_SIGNED ? "signed" : "binary");
    size_t vn = (size_t)prm.mu * TV_N;
    /* leaves: commitments to uniformly random shares with seed-derived randomness, as produced
       by accepted ballots (the ballot proofs do not enter the aggregation cost) */
    poly *lc = malloc(sizeof(poly) * (size_t)n * NV * prm.rows), *lm = malloc(sizeof(poly) * (size_t)n * NV);
    int64_t *lr = malloc(sizeof(int64_t) * (size_t)n * NV * vn);
    prg g;
    prg_init(&g, 6, s, 32);
    long expect = 0;
    int T[16];
    for (int i = 0; i < NV; i++) {
        tv_voter_secret sec;
        tv_secret_alloc(&sec, &prm);
        poly v;
        poly_zero(&v);
        v.c[0] = (uint64_t)(i % 3 == 0);
        expect += (long)v.c[0];
        shamir_share(sec.m, &v, &prm, &g);
        for (int k = 1; k <= n; k++) {
            uint8_t sd[32];
            prg_bytes(&g, sd, 32);
            size_t idx = (size_t)(k - 1) * NV + (size_t)i;
            tv_rand_from_seed(lr + idx * vn, &pub, (uint64_t)i, k, sd);
            lm[idx] = sec.m[k];
            tv_commit(lc + idx * prm.rows, &pub, &lm[idx], lr + idx * vn);
        }
        tv_secret_free(&sec);
    }
    agg_params ap;
    agg_params_init(&ap, &prm, NV, mode, 64);
    agg_params_print(&ap, &prm);
    agg_contrib *c = calloc((size_t)n, sizeof(agg_contrib));
    agg_state *st = calloc((size_t)n, sizeof(agg_state));
    int *succ = calloc((size_t)n, sizeof(int));
    double r1ms = 0, r2ms = 0, vms = 0, att_total = 0;
    size_t by1 = 0, by2 = 0, by_failed = 0;
    int rounds = 0, nsucc = 0;
    uint8_t bb1[TV_HASHBYTES];
    while (nsucc < t && rounds < 20) {
        rounds++;
        for (int k = 0; k < n; k++) {
            if (rounds > 1) { agg_state_free(&st[k]); agg_contrib_free(&c[k]); }
            uint8_t as[32];
            seed_from(as, g_seed, 5000 + (unsigned)(k + 100 * rounds));
            agg_stats sa;
            agg_round1(&st[k], &c[k], &sa, &pub, &ap, k + 1, lc + (size_t)k * NV * prm.rows, lm + (size_t)k * NV, lr + (size_t)k * NV * vn, NV, as);
            r1ms += sa.t_round1_ms;
            by1 += sa.bytes_round1;
        }
        agg_bb1_digest(bb1, &pub, c, n);
        nsucc = 0;
        for (int k = 0; k < n; k++) {
            agg_stats sa;
            succ[k] = agg_round2(&st[k], &c[k], &sa, &pub, &ap, bb1) == 0;
            nsucc += succ[k];
            r2ms += sa.t_round2_ms;
            att_total += sa.attempts_used;
            if (succ[k]) by2 += sa.bytes_round2;
        }
        if (nsucc < t) for (int k = 0; k < n; k++) if (succ[k]) by_failed += agg_contrib_bytes(&c[k], &pub, &ap, 2);
    }
    int nt = 0;
    const poly *V[16];
    int allok = 1;
    for (int k = 0; k < n; k++) {
        if (!succ[k]) continue;
        double t0 = now_ms();
        allok &= agg_verify(&pub, &ap, &c[k], lc + (size_t)k * NV * prm.rows, NV, bb1);
        vms += now_ms() - t0;
        if (nt < t) { T[nt] = k + 1; V[nt] = c[k].V; nt++; }
    }
    long cnt[1];
    double t0 = now_ms();
    allok &= nt == t && tv_combine(cnt, &pub, T, V, NV) == 0 && cnt[0] == expect;
    double combine_ms = now_ms() - t0;
    int part = n;   /* all n authorities took part in every round */
    row1("agg", "blocks E", cfg, ap.E, "count");
    row1("agg", "rounds until t succeeded", cfg, rounds, "count");
    row1("agg", "round 1 per authority per round", cfg, r1ms / (part * rounds) / 1000, "s");
    row1("agg", "round 2 per authority per round", cfg, r2ms / (part * rounds) / 1000, "s");
    row1("agg", "round 2 per attempt", cfg, r2ms / att_total / 1000, "s");
    row1("agg", "round 1 per block per attempt", cfg, r1ms / (part * rounds) / ap.kappa / ap.E / 1000, "s");
    row1("agg", "round 2 per block per attempt", cfg, r2ms / att_total / ap.E / 1000, "s");
    row1("agg", "VerAgg per authority", cfg, vms / (nsucc ? nsucc : 1) / 1000, "s");
    row1("agg", "VerAgg per block", cfg, vms / (nsucc ? nsucc : 1) / ap.E / 1000, "s");
    row1("agg", "Combine", cfg, combine_ms, "ms");
    row1("agg", "published round 1 per authority per round", cfg, (double)by1 / (part * rounds) / 1024, "KiB");
    row1("agg", "published round 2 per successful authority", cfg, (double)by2 / (nsucc ? nsucc : 1) / 1024, "KiB");
    row1("agg", "published per ballot per authority", cfg, ((double)by1 / (part * rounds) + (double)by2 / (nsucc ? nsucc : 1)) / NV / 1024, "KiB");
    row1("agg", "responses withheld in failed rounds (not published)", cfg, (double)by_failed / 1024, "KiB");
    row1("agg", "tally correct and all VerAgg ok", cfg, allok, "bool");
    for (int k = 0; k < n; k++) { agg_state_free(&st[k]); agg_contrib_free(&c[k]); }
    free(c); free(st); free(succ); free(lc); free(lm); free(lr);
    tv_pub_free(&pub);
}

int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IONBF, 0);
    int NV = 1000, reps = 30, d = 7, skip_agg = 0, only_agg = 0;
    uint64_t q = 4398046511953ULL;    /* 2^42 + 849: the set chosen by tools/chain.py for N_V = 10^4 */
    const char *out = "results/bench.csv";
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-nv") && i + 1 < argc) NV = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-reps") && i + 1 < argc) reps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-seed") && i + 1 < argc) g_seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-d") && i + 1 < argc) d = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-q") && i + 1 < argc) q = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-out") && i + 1 < argc) out = argv[++i];
        else if (!strcmp(argv[i], "-skip-agg")) skip_agg = 1;
        else if (!strcmp(argv[i], "-only-agg")) only_agg = 1;
    }
    detect_cpu();
    csv = fopen(out, "a");
    if (!csv) { fprintf(stderr, "cannot open %s (does the directory exist?)\n", out); return 1; }
    fseek(csv, 0, SEEK_END);
    if (ftell(csv) == 0) fprintf(csv, "group,name,config,median,q1,q3,min,reps,unit,cpu,compiler,flags,commit,seed\n");
    printf("T-EVOLVE benchmarks | cpu: %s | compiler: %s | flags: %s | commit: %s | seed: %llu\n", cpu_name, __VERSION__, BUILD_FLAGS, GIT_COMMIT, g_seed);
    if (!only_agg) {
        tv_params prm;
        if (tv_params_init(&prm, 4, 3, 1, TV_W_FREE, d, q)) { fprintf(stderr, "bad q\n"); return 1; }
        uint8_t s[32];
        seed_from(s, g_seed, 1);
        tv_pub pub;
        tv_pub_init(&pub, &prm, s);
        bench_micro(&pub, reps);
        tv_pub_free(&pub);
        bench_tballot(4, 3, 1, TV_W_FREE, d, q, reps);
        bench_tballot(4, 3, 2, 1, d, q, reps);
        bench_tballot(5, 3, 1, TV_W_FREE, d, q, reps);
        bench_eballot(4, d, q, reps);
    }
    if (!skip_agg) bench_agg(NV, 4, 3, d, q, C2_SIGNED);
    fclose(csv);
    return 0;
}
