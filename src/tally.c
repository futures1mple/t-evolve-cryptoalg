#include <stdlib.h>
#include "shamir.h"
#include "tally.h"

int tv_combine(long *counts, const tv_pub *pub, const int *T, const poly *const *V, long nacc) {
    const tv_params *p = &pub->prm;
    uint64_t lam[64];
    if (p->t > 64) return -1;
    shamir_lagrange(lam, T, p->t);
    int ok = 0;
    for (int a = 0; a < p->L; a++) {
        poly acc, tmp;
        poly_zero(&acc);
        for (int i = 0; i < p->t; i++) { poly_scale(&tmp, &V[i][a], to_mont(lam[i])); poly_add(&acc, &acc, &tmp); }
        for (int j = 1; j < TV_N; j++) if (acc.c[j] != 0) ok = -1;
        if (acc.c[0] > (uint64_t)nacc) ok = -1;
        counts[a] = (long)acc.c[0];
    }
    return ok;
}
