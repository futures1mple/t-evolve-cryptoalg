#include <math.h>
#include <string.h>
#include "params.h"

/* Rigorous bound on the norm of the shift s = (c r_0, ..., c r_n, f_1 r_0, ..., f_L r_0).
 *
 * The coefficients of r_j are i.i.d. D_{Z,sigma}, which is sigma-subgaussian
 * [Micciancio-Peikert 2012, Lemma 2.8]. For a fixed challenge tuple, s = Phi r is linear in r;
 * with Sigma = Phi^T Phi, the quadratic-form tail bound of Hsu, Kakade and Zhang (2012, Thm. 1)
 * gives Pr[ ||s||^2 > sigma^2 (tr + 2 sqrt(tr(Sigma^2) tau) + 2 ||Sigma|| tau) ] <= e^{-tau}.
 * For challenges with 60 coefficients in {-1,1}:  tr(Sigma) = 60 N mu (n+1+L) exactly,
 * ||Sigma|| <= 60^2 (1+L) (r_0 appears in 1+L blocks, ||Rot(c)|| <= ||c||_1 = 60), and
 * tr(Sigma^2) <= ||Sigma|| tr(Sigma). */
static double shift_bound(int N, int mu, int n, int L, double sigma, double tau) {
    double tr = 60.0 * N * mu * (n + 1 + L);
    double lam = 3600.0 * (1 + L);
    return sigma * sqrt(tr + 2.0 * sqrt(lam * tr * tau) + 2.0 * lam * tau);
}

int tv_params_init(tv_params *p, int n, int t, int L, int w, int d, uint64_t q) {
    memset(p, 0, sizeof *p);
    if (n < 1 || t < 1 || t > n || L < 1 || d < 1) return -1;
    if (w != TV_W_FREE && (w < 0 || w > L)) return -1;
    if ((uint64_t)n >= q) return -1;
    p->n = n; p->t = t; p->L = L; p->w = w;
    p->q = q; p->d = d;
    p->mu = 2 * d + L;
    p->rows = d + L;
    p->sigma = 1.0;
    p->alpha = 11.0;
    p->tau_bits = 128.0 + 40.0;        /* 2^-168 per attempt, hence 2^-128 over 2^40 attempts */
    p->T = shift_bound(TV_N, p->mu, n, L, p->sigma, p->tau_bits * log(2.0));
    p->sigma_J = p->alpha * p->T;
    p->logM = 12.0 / p->alpha + 1.0 / (2.0 * p->alpha * p->alpha);
    p->B_J = 2.0 * p->sigma_J * sqrt((double)tv_resp_coeffs(p));
    p->fanin = 30; p->k_blk = 500; p->ell = 517; p->kappa = 4;
    return ring_init(q);
}
