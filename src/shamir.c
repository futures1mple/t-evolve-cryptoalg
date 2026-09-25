#include <stdlib.h>
#include "shamir.h"

void shamir_parity(uint64_t *H, const tv_params *prm) {
    int n = prm->n, rows = n + 1 - prm->t;
    for (int i = 0; i <= n; i++) {
        uint64_t u = 1;
        for (int l = 0; l <= n; l++) {
            if (l == i) continue;
            u = mod_mul(u, mod_sub((uint64_t)i, (uint64_t)l));   /* alpha_i - alpha_l = i - l */
        }
        u = mod_inv(u);
        uint64_t apow = 1;                                          /* alpha_i^j, with 0^0 = 1 */
        for (int j = 0; j < rows; j++) {
            H[j * (n + 1) + i] = mod_mul(u, apow);
            apow = mod_mul(apow, (uint64_t)i);
        }
    }
}

void shamir_share(poly *m, const poly *v, const tv_params *prm, prg *g) {
    int n = prm->n, t = prm->t, L = prm->L;
    poly *coef = malloc(sizeof(poly) * (size_t)((t > 1 ? t - 1 : 1) * L));
    for (int l = 0; l < (t - 1) * L; l++) sample_uniform_poly(&coef[l], g);
    for (int a = 0; a < L; a++) m[a] = v[a];
    for (int k = 1; k <= n; k++) {
        for (int a = 0; a < L; a++) {
            /* Horner: P(k) = v + k(a_1 + k(a_2 + ...)) */
            poly acc;
            poly_zero(&acc);
            for (int l = t - 1; l >= 1; l--) {
                poly_add(&acc, &acc, &coef[(l - 1) * L + a]);
                poly_scale(&acc, &acc, to_mont((uint64_t)k));
            }
            poly_add(&m[k * L + a], &acc, &v[a]);
        }
    }
    free(coef);
}

void shamir_lagrange(uint64_t *lam, const int *T, int t) {
    for (int i = 0; i < t; i++) {
        uint64_t num = 1, den = 1;
        for (int j = 0; j < t; j++) {
            if (j == i) continue;
            num = mod_mul(num, (uint64_t)T[j]);
            den = mod_mul(den, mod_sub((uint64_t)T[j], (uint64_t)T[i]));
        }
        lam[i] = mod_mul(num, mod_inv(den));
    }
}
