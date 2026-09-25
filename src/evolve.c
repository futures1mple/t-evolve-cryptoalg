#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "codec.h"
#include "evolve.h"

#define VN(p) ((size_t)(p)->mu * TV_N)

void ev_params_init(ev_params *ep, const tv_params *p) {
    /* shift f r with r = sum of n independent D_sigma vectors: tr = 60 N mu n, ||Sigma|| <= 3600 n */
    double tr = 60.0 * TV_N * p->mu * p->n, lam = 3600.0 * p->n, tau = p->tau_bits * log(2.0);
    ep->T = p->sigma * sqrt(tr + 2.0 * sqrt(tr * lam * tau) + 2.0 * lam * tau);
    ep->sigma_OR = p->alpha * ep->T;
    ep->logM = p->logM;
    ep->B_OR = 2.0 * ep->sigma_OR * sqrt(2.0 * (double)VN(p));
}

void ev_ballot_alloc(ev_ballot *b, const tv_params *p) {
    memset(b, 0, sizeof *b);
    b->c = calloc((size_t)p->n * p->rows, sizeof(poly));
    b->e = calloc((size_t)p->n, TV_CT_BYTES);
    b->r0 = calloc(VN(p), sizeof(int64_t));
    b->r1 = calloc(VN(p), sizeof(int64_t));
}
void ev_ballot_free(ev_ballot *b) { free(b->c); free(b->e); free(b->r0); free(b->r1); memset(b, 0, sizeof *b); }

static void stub_enc(uint8_t out[TV_CT_BYTES], const tv_akey *key, uint64_t id, const uint8_t in[TV_CT_BYTES]) {
    uint8_t buf[40], pad[TV_CT_BYTES];
    memcpy(buf, key->key, 32);
    for (int i = 0; i < 8; i++) buf[32 + i] = (uint8_t)(id >> (8 * i));
    shake256(pad, TV_CT_BYTES, buf, sizeof buf);
    for (int i = 0; i < TV_CT_BYTES; i++) out[i] = in[i] ^ pad[i];
}

static void sum_com(poly *s, const tv_params *p, const poly *c) {
    for (int i = 0; i < p->rows; i++) {
        s[i] = c[i];
        for (int j = 1; j < p->n; j++) poly_add(&s[i], &s[i], &c[(size_t)j * p->rows + i]);
    }
}

static void st_digest(uint8_t out[64], const tv_pub *pub, const ev_ballot *b) {
    const tv_params *p = &pub->prm;
    keccak_state st;
    shake256_init(&st);
    keccak_absorb(&st, (const uint8_t *)"EVOLVE/bal-st", 13);
    keccak_absorb(&st, pub->digest, TV_HASHBYTES);
    uint8_t idb[8];
    for (int i = 0; i < 8; i++) idb[i] = (uint8_t)(b->id >> (8 * i));
    keccak_absorb(&st, idb, 8);
    absorb_polys(&st, b->c, (size_t)p->n * p->rows);
    keccak_absorb(&st, b->e, (size_t)p->n * TV_CT_BYTES);
    keccak_finalize(&st);
    keccak_squeeze(&st, out, 64);
}

static void ch_hash(uint8_t chat[TV_HASHBYTES], const tv_params *p, const uint8_t std[64], const poly *t0, const poly *t1) {
    keccak_state st;
    shake256_init(&st);
    keccak_absorb(&st, (const uint8_t *)"EVOLVE/bal-ch", 13);
    keccak_absorb(&st, std, 64);
    absorb_polys(&st, t0, (size_t)p->rows);
    absorb_polys(&st, t1, (size_t)p->rows);
    keccak_finalize(&st);
    keccak_squeeze(&st, chat, TV_HASHBYTES);
}

static void expand_pi(perm_ch *pi, const uint8_t chat[TV_HASHBYTES]) {
    prg g;
    prg_init(&g, 0x31, chat, TV_HASHBYTES);
    sample_perm(pi, &g);
}

/* t = C x + f (0;bit) - f c, x integer */
static void or_t(poly *t, const tv_pub *pub, const int64_t *x, const challenge *f, int bit, const poly *c) {
    const tv_params *p = &pub->prm;
    poly *xh = malloc(sizeof(poly) * p->mu), tmp;
    vec_ntt_from_int(xh, x, p->mu);
    C_mul(t, pub, 0, p->rows, xh);
    free(xh);
    if (!f) return;
    for (int i = 0; i < p->rows; i++) {
        poly_mul_sparse(&tmp, &c[i], f->pos, f->sgn, TV_CH_W);
        poly_sub(&t[i], &t[i], &tmp);
    }
    if (bit)
        for (int k = 0; k < TV_CH_W; k++) {
            uint64_t *cf = &t[p->d].c[f->pos[k]];
            *cf = f->sgn[k] > 0 ? mod_add(*cf, 1) : mod_sub(*cf, 1);
        }
}

int ev_vote(ev_ballot *b, poly *shares, int64_t *rnd, int *attempts, const tv_pub *pub, const ev_params *ep,
            const tv_akey *keys, uint64_t id, int vote, const uint8_t seed[32]) {
    const tv_params *p = &pub->prm;
    if (p->L != 1 || (vote != 0 && vote != 1)) return -1;
    int n = p->n, rows = p->rows;
    size_t vn = VN(p);
    prg g;
    prg_init(&g, 0x02, seed, 32);
    b->id = id;
    /* additive shares of the vote and their commitments; seeds encrypted */
    poly sum;
    poly_zero(&sum);
    for (int j = 0; j < n - 1; j++) { sample_uniform_poly(&shares[j], &g); poly_add(&sum, &sum, &shares[j]); }
    poly vp;
    poly_zero(&vp);
    vp.c[0] = (uint64_t)vote;
    poly_sub(&shares[n - 1], &vp, &sum);
    for (int j = 0; j < n; j++) {
        uint8_t s[32];
        prg_bytes(&g, s, 32);
        tv_rand_from_seed(rnd + (size_t)j * vn, pub, id, j + 1, s);
        tv_commit(&b->c[(size_t)j * rows], pub, &shares[j], rnd + (size_t)j * vn);
        stub_enc(b->e + (size_t)j * TV_CT_BYTES, &keys[j], id, s);
    }
    /* OR-proof on the sum commitment, randomness r = sum r_j */
    poly *c = malloc(sizeof(poly) * rows), *t0 = malloc(sizeof(poly) * rows), *t1 = malloc(sizeof(poly) * rows);
    sum_com(c, p, b->c);
    int64_t *r = calloc(vn, sizeof(int64_t)), *rho = malloc(sizeof(int64_t) * vn), *sh = malloc(sizeof(int64_t) * vn);
    for (int j = 0; j < n; j++) for (size_t i = 0; i < vn; i++) r[i] += rnd[(size_t)j * vn + i];
    uint8_t std[64];
    st_digest(std, pub, b);
    gauss_sampler gs;
    gauss_init(&gs, ep->sigma_OR);
    int m = vote, s = 1 - m, att = 0;
    int64_t *rs = s ? b->r1 : b->r0, *rm = m ? b->r1 : b->r0;
    for (;;) {
        att++;
        challenge fs, fm;
        gauss_vec(&gs, &g, rs, vn);
        sample_challenge(&fs, &g);
        or_t(s ? t1 : t0, pub, rs, &fs, s, c);
        gauss_vec(&gs, &g, rho, vn);
        or_t(m ? t1 : t0, pub, rho, NULL, 0, c);
        ch_hash(b->chat, p, std, t0, t1);
        perm_ch pi;
        expand_pi(&pi, b->chat);
        if (m == 1) { b->f0 = fs; perm_apply(&fm, &pi, &fs); }
        else { perm_apply_inv(&fm, &pi, &fs); b->f0 = fm; }
        ch_mul_int(sh, &fm, r, (size_t)p->mu);
        long double zv = 0, vv = 0;
        for (size_t i = 0; i < vn; i++) { rm[i] = rho[i] + sh[i]; zv += (long double)rm[i] * sh[i]; vv += (long double)sh[i] * sh[i]; }
        double lhs = (double)((-2.0L * zv + vv) / (2.0L * (long double)ep->sigma_OR * ep->sigma_OR)) - ep->logM;
        if (log(prg_unif(&g)) > lhs) continue;
        long double n2 = 0;
        for (size_t i = 0; i < vn; i++) n2 += (long double)b->r0[i] * b->r0[i] + (long double)b->r1[i] * b->r1[i];
        if (sqrtl(n2) > ep->B_OR) continue;
        break;
    }
    if (attempts) *attempts = att;
    gauss_free(&gs);
    free(c); free(t0); free(t1); free(r); free(rho); free(sh);
    return 0;
}

int ev_verify(const tv_pub *pub, const ev_params *ep, const ev_ballot *b) {
    const tv_params *p = &pub->prm;
    size_t vn = VN(p);
    long double n2 = 0;
    for (size_t i = 0; i < vn; i++) n2 += (long double)b->r0[i] * b->r0[i] + (long double)b->r1[i] * b->r1[i];
    if (!(sqrtl(n2) <= ep->B_OR)) return 0;
    int rows = p->rows;
    poly *c = malloc(sizeof(poly) * rows), *t0 = malloc(sizeof(poly) * rows), *t1 = malloc(sizeof(poly) * rows);
    sum_com(c, p, b->c);
    perm_ch pi;
    expand_pi(&pi, b->chat);
    challenge f1;
    perm_apply(&f1, &pi, &b->f0);
    or_t(t0, pub, b->r0, &b->f0, 0, c);
    or_t(t1, pub, b->r1, &f1, 1, c);
    uint8_t std[64], chat[TV_HASHBYTES];
    st_digest(std, pub, b);
    ch_hash(chat, p, std, t0, t1);
    free(c); free(t0); free(t1);
    return memcmp(chat, b->chat, TV_HASHBYTES) == 0;
}

size_t ev_ballot_bytes(const tv_pub *pub, const ev_params *ep, const ev_ballot *b, tv_sizes *sz) {
    const tv_params *p = &pub->prm;
    size_t vn = VN(p);
    int k = rice_param(ep->sigma_OR);
    bitw w;
    bw_init(&w, NULL, 0);
    bw_put(&w, 0, 32); bw_put(&w, 0, 32);
    size_t s0 = w.bytes;
    for (size_t i = 0; i < (size_t)p->n * p->rows; i++) pack_poly(&w, &b->c[i]);
    size_t s1 = w.bytes;
    for (size_t i = 0; i < (size_t)p->n * TV_CT_BYTES; i++) bw_put(&w, 0, 8);
    size_t s2 = w.bytes;
    for (int i = 0; i < TV_HASHBYTES; i++) bw_put(&w, 0, 8);
    for (int i = 0; i < TV_CH_W; i++) bw_put(&w, 0, 9);
    for (size_t i = 0; i < vn; i++) { rice_put(&w, b->r0[i], k); rice_put(&w, b->r1[i], k); }
    size_t tot = bw_finish(&w);
    if (sz) { sz->total = tot; sz->commitments = s1 - s0; sz->ciphertexts = s2 - s1; sz->proof = tot - s2; }
    return tot;
}
