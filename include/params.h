/* Parameters of T-EVOLVE and the quantities derived from them (see docs/PARAMETERS.md). */
#ifndef TV_PARAMS_H
#define TV_PARAMS_H

#include <stdint.h>
#include "ring.h"

#define TV_W_FREE (-1)          /* no weight constraint (w = bottom) */
#define TV_SEEDBYTES 32
#define TV_HASHBYTES 32

typedef struct {
    /* election */
    int n, t;                   /* authorities, threshold */
    int L, w;                   /* candidates, weight (TV_W_FREE for none) */
    /* lattice */
    uint64_t q;
    int d;                      /* module rank */
    int mu;                     /* commitment width 2d+L */
    int rows;                   /* d+L */
    double sigma;               /* commitment randomness D_sigma */
    /* ballot proof (joint rejection) */
    double alpha;               /* sigma_J = alpha * T */
    double tau_bits;            /* failure exponent of the norm bound on the shift */
    double T;                   /* bound on the l2 norm of the shift */
    double sigma_J;
    double logM;                /* ln M, M = exp(12/alpha + 1/(2 alpha^2)) */
    double B_J;                 /* l2 bound on the concatenation of all responses */
    /* aggregation (Baum-Lyubashevsky 2017), filled by agg_params() */
    int fanin, k_blk, ell;
    int kappa;
} tv_params;

/* fills every derived field; returns 0 on success */
int tv_params_init(tv_params *p, int n, int t, int L, int w, int d, uint64_t q);

/* number of coefficients in the concatenation of all ballot-proof responses */
static inline long tv_resp_coeffs(const tv_params *p) { return (long)TV_N * p->mu * (p->n + 1 + 2 * p->L); }

#endif
