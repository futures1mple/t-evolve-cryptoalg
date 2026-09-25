#include <stdlib.h>
#include <string.h>
#include "bdlop.h"
#include "sample.h"
#include "shamir.h"

int tv_pub_init(tv_pub *pub, const tv_params *prm, const uint8_t seed[TV_SEEDBYTES]) {
    memset(pub, 0, sizeof *pub);
    pub->prm = *prm;
    memcpy(pub->seed, seed, TV_SEEDBYTES);
    if (ring_init(prm->q)) return -1;
    size_t nC = (size_t)prm->rows * prm->mu;
    pub->C = malloc(sizeof(poly) * nC);
    int hr = prm->n + 1 - prm->t, hc = prm->n + 1;
    pub->H = malloc(sizeof(uint64_t) * (size_t)(hr > 0 ? hr : 1) * hc);
    pub->H_m = malloc(sizeof(uint64_t) * (size_t)(hr > 0 ? hr : 1) * hc);
    if (!pub->C || !pub->H || !pub->H_m) return -1;
    prg g;
    prg_init(&g, 0x10, seed, TV_SEEDBYTES);          /* expansion of C from its seed */
    for (size_t i = 0; i < nC; i++) {
        sample_uniform_poly(&pub->C[i], &g);         /* uniform; directly interpreted in the NTT domain */
        poly_to_mont(&pub->C[i]);
    }
    pub->Bsum = malloc(sizeof(poly) * prm->mu);
    if (!pub->Bsum) return -1;
    for (int j = 0; j < prm->mu; j++) {
        poly_zero(&pub->Bsum[j]);
        for (int a = 0; a < prm->L; a++) poly_add(&pub->Bsum[j], &pub->Bsum[j], &pub->C[(size_t)(prm->d + a) * prm->mu + j]);
    }
    shamir_parity(pub->H, prm);
    for (int i = 0; i < hr * hc; i++) pub->H_m[i] = to_mont(pub->H[i]);
    /* digest of the parameters */
    keccak_state st;
    shake256_init(&st);
    int64_t vals[7] = {prm->n, prm->t, prm->L, prm->w, prm->d, (int64_t)prm->q, prm->mu};
    uint8_t le[56];
    for (int i = 0; i < 7; i++) for (int b = 0; b < 8; b++) le[8 * i + b] = (uint8_t)((uint64_t)vals[i] >> (8 * b));
    keccak_absorb(&st, (const uint8_t *)"T-EVOLVE/par", 12);
    keccak_absorb(&st, le, sizeof le);
    keccak_absorb(&st, seed, TV_SEEDBYTES);
    keccak_finalize(&st);
    keccak_squeeze(&st, pub->digest, TV_HASHBYTES);
    return 0;
}

void tv_pub_free(tv_pub *pub) {
    free(pub->C); free(pub->H); free(pub->H_m); free(pub->Bsum);
    memset(pub, 0, sizeof *pub);
}

void vec_ntt_from_int(poly *out, const int64_t *x, int k) {
    for (int j = 0; j < k; j++) {
        poly_from_int(&out[j], x + (size_t)j * TV_N);
        poly_ntt(&out[j]);
    }
}

void C_mul(poly *out, const tv_pub *pub, int row0, int nrows, const poly *xh) {
    for (int i = 0; i < nrows; i++) {
        poly_dot_mont(&out[i], C_row(pub, row0 + i), xh, (size_t)pub->prm.mu);
        poly_invntt(&out[i]);
    }
}

void tv_commit(poly *c, const tv_pub *pub, const poly *m, const int64_t *r) {
    const tv_params *p = &pub->prm;
    poly *rh = malloc(sizeof(poly) * p->mu);
    vec_ntt_from_int(rh, r, p->mu);
    C_mul(c, pub, 0, p->rows, rh);
    for (int a = 0; a < p->L; a++) poly_add(&c[p->d + a], &c[p->d + a], &m[a]);
    free(rh);
}
