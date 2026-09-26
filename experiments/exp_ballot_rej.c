/* Experiment for the ballot proof with a single rejection step.
 *
 * For every configuration, generates independent ballots (fresh voter randomness, fresh
 * challenges in every attempt, sigma_J fixed in advance by params.c, exact discrete Gaussian
 * sampler), runs the complete prover and verifier, and records:
 *   - the number of attempts (expected: M = exp(43669/39762) ~ 3.00, independent of L);
 *   - the ratio ||shift|| / T (the rejection lemma needs <= 1);
 *   - the response norm relative to B_J, and that every ballot verifies after encoding;
 *   - a chi-square goodness-of-fit test of the accepted responses against D_{sigma_J};
 *   - a two-sample test of witness independence: the OR responses r_{1,0} of branch 0 are
 *     produced by rejection sampling when the vote is 0 and sampled directly when it is 1.
 *
 * usage: exp_ballot_rej [ballots] [seed] [outdir]  */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ballot.h"
#include "codec.h"
#include "timer.h"

#define NBINS 64      /* bins of z/sigma on [-4,4], plus two tails */

static double Phi(double x) { return 0.5 * erfc(-x / sqrt(2.0)); }

typedef struct { double cnt[NBINS + 2]; double n; } hist;
static void hist_add(hist *h, const int64_t *x, size_t k, double sigma) {
    for (size_t i = 0; i < k; i++) {
        double u = (double)x[i] / sigma;
        int b;
        if (u < -4) b = 0;
        else if (u >= 4) b = NBINS + 1;
        else { b = 1 + (int)((u + 4) / 8 * NBINS); if (b > NBINS) b = NBINS; }
        h->cnt[b] += 1;
    }
    h->n += (double)k;
}
/* goodness of fit against N(0,1) (sigma_J > 2^14: the discretization is negligible) */
static double chi2_gof(const hist *h, int *dof) {
    double chi = 0;
    for (int b = 0; b < NBINS + 2; b++) {
        double lo = b == 0 ? -INFINITY : -4 + 8.0 * (b - 1) / NBINS;
        double hi = b == NBINS + 1 ? INFINITY : -4 + 8.0 * b / NBINS;
        double e = h->n * (Phi(hi) - Phi(lo));
        chi += (h->cnt[b] - e) * (h->cnt[b] - e) / e;
    }
    *dof = NBINS + 1;
    return chi;
}
/* homogeneity test of two histograms */
static double chi2_two(const hist *a, const hist *b, int *dof) {
    double chi = 0;
    int k = 0;
    for (int i = 0; i < NBINS + 2; i++) {
        double tot = a->cnt[i] + b->cnt[i];
        if (tot == 0) continue;
        double ea = tot * a->n / (a->n + b->n), eb = tot * b->n / (a->n + b->n);
        chi += (a->cnt[i] - ea) * (a->cnt[i] - ea) / ea + (b->cnt[i] - eb) * (b->cnt[i] - eb) / eb;
        k++;
    }
    *dof = k - 1;
    return chi;
}
/* upper tail of chi-square via the Wilson-Hilferty approximation */
static double chi2_pvalue(double x, int k) {
    double z = (pow(x / k, 1.0 / 3) - (1 - 2.0 / (9 * k))) / sqrt(2.0 / (9 * k));
    return 1 - Phi(z);
}

static void seed_from(uint8_t s[32], uint64_t a, uint64_t b) {
    memset(s, 0, 32);
    for (int i = 0; i < 8; i++) { s[i] = (uint8_t)(a >> (8 * i)); s[8 + i] = (uint8_t)(b >> (8 * i)); }
}

static int cmpd(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return x < y ? -1 : x > y; }

int main(int argc, char **argv) {
    int nb = argc > 1 ? atoi(argv[1]) : 1000;
    uint64_t seed = argc > 2 ? strtoull(argv[2], NULL, 10) : 1;
    const char *outdir = argc > 3 ? argv[3] : "results";
    struct { int n, t, L, w, d; uint64_t q; } cfg[] = {
        {4, 3, 1, TV_W_FREE, 7, 4398046510961ULL},
        {4, 3, 2, 1, 7, 4398046510961ULL},
        {4, 3, 5, 1, 7, 4398046510961ULL},
        {4, 3, 10, 1, 7, 4398046510961ULL},
        {5, 3, 1, TV_W_FREE, 7, 4398046510961ULL},
    };
    int ncfg = sizeof cfg / sizeof cfg[0];
    char path[512];
    snprintf(path, sizeof path, "%s/ballot_rej_summary.csv", outdir);
    FILE *fs = fopen(path, "w");
    if (!fs) { fprintf(stderr, "cannot write %s (create the directory first)\n", path); return 1; }
    fprintf(fs, "n,t,L,w,d,logq,sigma_J,log2_sigma_J,T,B_J,M,ballots,attempts_mean,attempts_sd,attempts_ci95,"
                "accept_rate,max_shift_over_T,max_norm_over_BJ,verified,bytes_mean,prove_ms_median,verify_ms_median,"
                "gof_chi2,gof_dof,gof_p,indep_chi2,indep_dof,indep_p,seed\n");
    printf("%-3s %-3s %-4s %-8s %-8s %-10s %-9s %-9s %-9s %-9s %-7s %-7s\n", "n", "L", "w", "log2 sJ", "M", "attempts",
           "+-95%", "maxS/T", "max/BJ", "verified", "gof p", "indep p");
    for (int ci = 0; ci < ncfg; ci++) {
        tv_params prm;
        if (tv_params_init(&prm, cfg[ci].n, cfg[ci].t, cfg[ci].L, cfg[ci].w, cfg[ci].d, cfg[ci].q)) return 1;
        uint8_t s[32];
        seed_from(s, seed, 0xC0FFEE + (uint64_t)ci);
        tv_pub pub;
        tv_pub_init(&pub, &prm, s);
        tv_akey *keys = calloc((size_t)prm.n, sizeof(tv_akey));
        for (int k = 0; k < prm.n; k++) seed_from(keys[k].key, seed, 1000 + (uint64_t)k);
        snprintf(path, sizeof path, "%s/ballot_rej_n%d_L%d.csv", outdir, prm.n, prm.L);
        FILE *fb = fopen(path, "w");
        fprintf(fb, "ballot,vote,attempts,shift_over_T,norm_over_BJ,bytes,prove_ms,verify_ms,verified\n");
        tv_ballot b, b2;
        tv_voter_secret sec;
        tv_ballot_alloc(&b, &prm);
        tv_ballot_alloc(&b2, &prm);
        tv_secret_alloc(&sec, &prm);
        uint8_t *buf = malloc(tv_ballot_maxbytes(&prm));
        size_t vn = (size_t)prm.mu * TV_N;
        hist gof = {{0}, 0}, br0_real = {{0}, 0}, br0_sim = {{0}, 0};
        double sa = 0, sa2 = 0, maxS = 0, maxN = 0, bytes = 0;
        double *tp = malloc(sizeof(double) * (size_t)nb), *tv = malloc(sizeof(double) * (size_t)nb);
        int verified = 0;
        long total_att = 0;
        for (int i = 0; i < nb; i++) {
            /* admissible vote, drawn from the seed */
            int v[16] = {0};
            uint64_t r = seed * 1000003ULL + (uint64_t)i * 7919ULL + (uint64_t)ci;
            if (prm.w == TV_W_FREE) { for (int a = 0; a < prm.L; a++) v[a] = (int)((r >> a) & 1); }
            else { int pos = (int)(r % (uint64_t)prm.L); v[pos] = 1; }
            seed_from(s, seed, ((uint64_t)ci << 32) | (uint64_t)i);
            tv_prove_stats st;
            double t0 = now_ms();
            if (tv_vote(&b, &sec, &st, &pub, keys, (uint64_t)i, v, s)) { fprintf(stderr, "vote failed\n"); return 1; }
            double t1 = now_ms();
            tv_sizes sz;
            size_t len = tv_ballot_encode(buf, tv_ballot_maxbytes(&prm), &pub, &b, &sz);
            int ok = tv_ballot_decode(&b2, &pub, buf, len) == 0;
            double t2 = now_ms();
            ok = ok && tv_verify_ballot(&pub, &b2);
            double t3 = now_ms();
            verified += ok;
            tp[i] = t1 - t0; tv[i] = t3 - t2;
            sa += st.attempts; sa2 += (double)st.attempts * st.attempts; total_att += st.attempts;
            if (st.max_shift_ratio > maxS) maxS = st.max_shift_ratio;
            double nr = st.resp_norm / prm.B_J;
            if (nr > maxN) maxN = nr;
            bytes += (double)len;
            /* accepted witness-dependent responses: z_0..z_n and the real OR branches */
            hist_add(&gof, b.z, (size_t)(prm.n + 1) * vn, prm.sigma_J);
            for (int a = 0; a < prm.L; a++) hist_add(&gof, b.or_r + ((size_t)a * 2 + (size_t)v[a]) * vn, vn, prm.sigma_J);
            /* branch 0 of candidate 1: real (vote 0) versus simulated (vote 1) */
            hist_add(v[0] == 0 ? &br0_real : &br0_sim, b.or_r, vn, prm.sigma_J);
            fprintf(fb, "%d,%d,%d,%.4f,%.4f,%zu,%.3f,%.3f,%d\n", i, v[0], st.attempts, st.max_shift_ratio, nr, len, t1 - t0, t3 - t2, ok);
        }
        fclose(fb);
        double mean = sa / nb, sd = sqrt(sa2 / nb - mean * mean), ci95 = 1.96 * sd / sqrt((double)nb);
        int dg, di;
        double cg = chi2_gof(&gof, &dg), cii = chi2_two(&br0_real, &br0_sim, &di);
        double pg = chi2_pvalue(cg, dg), pi_ = chi2_pvalue(cii, di);
        qsort(tp, (size_t)nb, sizeof(double), cmpd);
        qsort(tv, (size_t)nb, sizeof(double), cmpd);
        double M = exp(prm.logM);
        printf("%-3d %-3d %-4d %-8.2f %-8.3f %-10.3f %-9.3f %-9.3f %-9.3f %4d/%-4d %-7.3f %-7.3f\n", prm.n, prm.L, prm.w,
               log2(prm.sigma_J), M, mean, ci95, maxS, maxN, verified, nb, pg, pi_);
        fprintf(fs, "%d,%d,%d,%d,%d,%u,%.1f,%.3f,%.1f,%.1f,%.4f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.1f,%.3f,%.3f,%.1f,%d,%.4f,%.1f,%d,%.4f,%llu\n",
                prm.n, prm.t, prm.L, prm.w, prm.d, RING.logq, prm.sigma_J, log2(prm.sigma_J), prm.T, prm.B_J, M, nb, mean, sd, ci95,
                (double)nb / (double)total_att, maxS, maxN, verified, bytes / nb, tp[nb / 2], tv[nb / 2], cg, dg, pg, cii, di, pi_,
                (unsigned long long)seed);
        fflush(fs);
        free(tp); free(tv); free(buf); free(keys);
        tv_ballot_free(&b); tv_ballot_free(&b2); tv_secret_free(&sec);
        tv_pub_free(&pub);
    }
    fclose(fs);
    return 0;
}
