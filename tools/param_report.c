/* Prints, as one JSON object, every quantity of the parameter chain for one configuration:
 * witnesses -> bounds on their norms -> rejection parameters -> extraction bounds -> beta_SIS,
 * together with expected sizes. The same functions are used by the implementation, so the
 * numbers in the paper and the code cannot drift apart.
 *
 * usage: param_report NV n t L w d q mode logQ      (w = -1 for no weight; mode 0 binary, 1 signed) */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "agg.h"

/* expected Golomb-Rice length (bits) of D_sigma with parameter k = floor(log2 sigma) */
static double rice_bits(double sigma) {
    int k = sigma >= 1 ? (int)floor(log2(sigma)) : 0;
    double e_hi = 0, Z = 0;
    long lim = (long)ceil(12 * sigma) + 2;
    double step = sigma > 2000 ? sigma / 2000 : 1;   /* integrate on a grid for large sigma */
    for (double x = 0; x <= lim; x += step) {
        double pr = exp(-x * x / (2 * sigma * sigma)) * (x == 0 ? 1 : 2);
        Z += pr;
        e_hi += pr * floor(x / pow(2, k));
    }
    return 1.0 + k + 1.0 + e_hi / Z;
}

int main(int argc, char **argv) {
    if (argc < 10) { fprintf(stderr, "usage: %s NV n t L w d q mode logQ\n", argv[0]); return 1; }
    long NV = atol(argv[1]);
    int n = atoi(argv[2]), t = atoi(argv[3]), L = atoi(argv[4]), w = atoi(argv[5]), d = atoi(argv[6]);
    uint64_t q = strtoull(argv[7], NULL, 10);
    c2_mode mode = atoi(argv[8]) ? C2_SIGNED : C2_BINARY;
    double logQ = atof(argv[9]);
    tv_params p;
    if (tv_params_init(&p, n, t, L, w < 0 ? TV_W_FREE : w, d, q)) { fprintf(stderr, "bad parameters\n"); return 1; }
    agg_params ap;
    agg_params_init(&ap, &p, NV, mode, logQ);
    double N = TV_N, Nmu = N * p.mu;
    double eta = 120.0;                               /* ||c - c'||_1 <= 2*60 */
    double beta_h = 2.0 * p.sigma * sqrt(Nmu);
    double Bbar = 2.0 * p.B_J;                        /* extraction bound of the ballot proof */
    double slack = mode == C2_SIGNED ? 2.0 : 1.0;     /* slack of the authority openings */
    double B_A = 4.0 * ap.zinf_mult * ap.sigma2 * sqrt(Nmu);   /* ||S_bar||_inf <= 2*||Z2-Z2'||_inf <= 36 sigma2 */
    double betaBL_inf = (mode == C2_SIGNED ? pow(2, 2.5) : pow(2, 1.5)) * ap.k * ap.sigma1 + pow(2, 2.5) * sqrt(60.0 * N) * ap.sigma2;
    double betaBL = betaBL_inf * sqrt(Nmu);
    double t_root = slack * beta_h + B_A, t_tree = (ap.fanin + 2) * B_A;
    double t_ballot = slack * Bbar + eta * B_A, t_or = eta * (Bbar + Bbar);
    double beta = t_root;
    if (t_tree > beta) beta = t_tree;
    if (t_ballot > beta) beta = t_ballot;
    if (t_or > beta) beta = t_or;
    if (betaBL > beta) beta = betaBL;
    /* sizes */
    double logq = log2((double)q), qbits = (double)RING.logq;   /* packed encoding uses ceil(log2 q) bits */
    double com_bytes = (n + 1) * p.rows * N * qbits / 8;
    double proof_bytes = (32 + ((n + 1) + 2.0 * L) * Nmu * rice_bits(p.sigma_J) / 8 + L * 68);
    double ct_bytes = n * 32;
    double auth_r1 = (ap.ninner * p.rows * N * qbits + L * N * qbits) / 8 + ap.kappa * 32;
    double auth_r2 = 33 + ((double)(ap.nnodes + ap.ninner) * Nmu * rice_bits(ap.sigma1) + (double)ap.E * Nmu * ap.ell * rice_bits(ap.sigma2)) / 8 + Nmu * 3 / 8;
    printf("{\"NV\":%ld,\"n\":%d,\"t\":%d,\"L\":%d,\"w\":%d,\"d\":%d,\"q\":%llu,\"log2q\":%.3f,\"mode\":\"%s\",\"logQ\":%.0f,",
           NV, n, t, L, w, d, (unsigned long long)q, logq, mode == C2_SIGNED ? "signed" : "binary", logQ);
    printf("\"mu\":%d,\"E\":%d,\"E_open\":%d,\"E_zero\":%d,\"nodes\":%d,\"inner\":%d,", p.mu, ap.E, ap.E_open, ap.E_zero, ap.nnodes, ap.ninner);
    printf("\"T_ballot\":%.6g,\"sigma_J\":%.6g,\"B_J\":%.6g,\"Bbar\":%.6g,", p.T, p.sigma_J, p.B_J, Bbar);
    printf("\"T1\":%.6g,\"T2\":%.6g,\"sigma1\":%.6g,\"sigma2\":%.6g,\"B_A\":%.6g,", ap.T1, ap.T2, ap.sigma1, ap.sigma2, B_A);
    printf("\"beta_h\":%.6g,\"t_root\":%.6g,\"t_tree\":%.6g,\"t_ballot\":%.6g,\"t_or\":%.6g,\"beta_BL\":%.6g,\"beta_SIS\":%.6g,\"log2beta\":%.3f,",
           beta_h, t_root, t_tree, t_ballot, t_or, betaBL, beta, log2(beta));
    printf("\"ballot_bytes\":%.0f,\"com_bytes\":%.0f,\"proof_bytes\":%.0f,\"ct_bytes\":%.0f,\"auth_round1_bytes\":%.0f,\"auth_round2_bytes\":%.0f,\"auth_per_ballot_bytes\":%.1f}\n",
           com_bytes + proof_bytes + ct_bytes + 8, com_bytes, proof_bytes, ct_bytes, auth_r1, auth_r2, (auth_r1 + auth_r2) / NV);
    return 0;
}
