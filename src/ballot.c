#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "ballot.h"
#include "codec.h"
#include "shamir.h"

#define VN(p) ((size_t)(p)->mu * TV_N)     /* coefficients of one randomness vector */

void tv_ballot_alloc(tv_ballot *b, const tv_params *p) {
    memset(b, 0, sizeof *b);
    b->c = calloc((size_t)(p->n + 1) * p->rows, sizeof(poly));
    b->e = calloc((size_t)p->n, TV_CT_BYTES);
    b->z = calloc((size_t)(p->n + 1) * VN(p), sizeof(int64_t));
    b->or_r = calloc((size_t)p->L * 2 * VN(p), sizeof(int64_t));
    b->f0 = calloc((size_t)p->L, sizeof(challenge));
}
void tv_ballot_free(tv_ballot *b) { free(b->c); free(b->e); free(b->z); free(b->or_r); free(b->f0); memset(b, 0, sizeof *b); }
void tv_secret_alloc(tv_voter_secret *s, const tv_params *p) {
    s->m = calloc((size_t)(p->n + 1) * p->L, sizeof(poly));
    s->r = calloc((size_t)(p->n + 1) * VN(p), sizeof(int64_t));
    s->seeds = calloc((size_t)p->n, 32);
}
void tv_secret_free(tv_voter_secret *s) { free(s->m); free(s->r); free(s->seeds); }

double tv_beta(const tv_params *p) { return 2.0 * p->sigma * sqrt((double)VN(p)); }
uint64_t tv_beta2(const tv_params *p) { return 4 * (uint64_t)VN(p); }      /* beta^2 for sigma = 1 */

/* r = SampleD_sigma(X(par, id, k, s)): the input is domain-separated by the parameters, the voter
   identity and the authority index, so that a seed reused in several ballots or shares yields
   independent randomness (needed by the bounds on the aggregation witnesses) */
void tv_rand_from_seed(int64_t *r, const tv_pub *pub, uint64_t id, int k, const uint8_t seed[32]) {
    const tv_params *p = &pub->prm;
    uint8_t in[TV_HASHBYTES + 12 + 32];
    memcpy(in, pub->digest, TV_HASHBYTES);
    for (int i = 0; i < 8; i++) in[TV_HASHBYTES + i] = (uint8_t)(id >> (8 * i));
    for (int i = 0; i < 4; i++) in[TV_HASHBYTES + 8 + i] = (uint8_t)((uint32_t)k >> (8 * i));
    memcpy(in + TV_HASHBYTES + 12, seed, 32);
    prg g;
    gauss_sampler gs;
    prg_init(&g, 0x20, in, sizeof in);
    gauss_init(&gs, p->sigma);
    gauss_vec(&gs, &g, r, VN(p));
    gauss_free(&gs);
}

/* placeholder PKE, NOT secure: e = s xor SHAKE256(key || id). It only fixes the data flow and
   the payload size; the cost of a real KEM (e.g. ML-KEM) is not included. */
static void pke_stub(uint8_t out[TV_CT_BYTES], const tv_akey *key, uint64_t id, const uint8_t in[TV_CT_BYTES]) {
    uint8_t buf[40], pad[TV_CT_BYTES];
    memcpy(buf, key->key, 32);
    for (int i = 0; i < 8; i++) buf[32 + i] = (uint8_t)(id >> (8 * i));
    shake256(pad, TV_CT_BYTES, buf, sizeof buf);
    for (int i = 0; i < TV_CT_BYTES; i++) out[i] = in[i] ^ pad[i];
}

/* ------------------------------------------------------------------ */

static void statement_digest(uint8_t out[64], const tv_pub *pub, const tv_ballot *b) {
    const tv_params *p = &pub->prm;
    keccak_state st;
    shake256_init(&st);
    keccak_absorb(&st, (const uint8_t *)"T-EVOLVE/bal-st", 15);
    keccak_absorb(&st, pub->digest, TV_HASHBYTES);
    uint8_t idb[8];
    for (int i = 0; i < 8; i++) idb[i] = (uint8_t)(b->id >> (8 * i));
    keccak_absorb(&st, idb, 8);
    absorb_polys(&st, b->c, (size_t)(p->n + 1) * p->rows);
    keccak_absorb(&st, b->e, (size_t)p->n * TV_CT_BYTES);
    keccak_finalize(&st);
    keccak_squeeze(&st, out, 64);
}

/* first messages, in the order they are hashed */
typedef struct {
    poly *w;      /* (n+1) * d */
    poly *v;      /* (n+1-t) * L */
    poly vw;      /* weight */
    poly *t;      /* L * 2 * (d+1) */
} firstmsg;

static void fm_alloc(firstmsg *f, const tv_params *p) {
    f->w = calloc((size_t)(p->n + 1) * p->d, sizeof(poly));
    f->v = calloc((size_t)(p->n + 1 - p->t) * p->L + 1, sizeof(poly));
    f->t = calloc((size_t)p->L * 2 * (p->d + 1), sizeof(poly));
    poly_zero(&f->vw);
}
static void fm_free(firstmsg *f) { free(f->w); free(f->v); free(f->t); }

static void fm_hash(uint8_t chat[TV_HASHBYTES], const tv_params *p, const uint8_t stdig[64], const firstmsg *f) {
    keccak_state st;
    shake256_init(&st);
    keccak_absorb(&st, (const uint8_t *)"T-EVOLVE/bal-ch", 15);
    keccak_absorb(&st, stdig, 64);
    absorb_polys(&st, f->w, (size_t)(p->n + 1) * p->d);
    absorb_polys(&st, f->v, (size_t)(p->n + 1 - p->t) * p->L);
    if (p->w != TV_W_FREE) absorb_polys(&st, &f->vw, 1);
    absorb_polys(&st, f->t, (size_t)p->L * 2 * (p->d + 1));
    keccak_finalize(&st);
    keccak_squeeze(&st, chat, TV_HASHBYTES);
}

/* c and pi_1..pi_L from the challenge digest */
static void expand_challenges(challenge *c, perm_ch *pi, int L, const uint8_t chat[TV_HASHBYTES]) {
    prg g;
    prg_init(&g, 0x30, chat, TV_HASHBYTES);
    sample_challenge(c, &g);
    for (int a = 0; a < L; a++) sample_perm(&pi[a], &g);
}

/* out (d+1 polys) = C_a x with C_a = (A ; b_a); xh in NTT form */
static void Ca_mul(poly *out, const tv_pub *pub, int a, const poly *xh) {
    C_mul(out, pub, 0, pub->prm.d, xh);
    C_mul(out + pub->prm.d, pub, pub->prm.d + a, 1, xh);
}

/* t (d+1 polys) -= f * (c_{0,1} ; c_{0,2,a}) and, if bit, row d += f */
static void or_adjust(poly *t, const tv_params *p, const tv_ballot *b, int a, const challenge *f, int bit) {
    poly tmp;
    for (int i = 0; i < p->d; i++) {
        poly_mul_sparse(&tmp, &b->c[i], f->pos, f->sgn, TV_CH_W);
        poly_sub(&t[i], &t[i], &tmp);
    }
    poly_mul_sparse(&tmp, &b->c[p->d + a], f->pos, f->sgn, TV_CH_W);
    poly_sub(&t[p->d], &t[p->d], &tmp);
    if (bit) {
        for (int k = 0; k < TV_CH_W; k++) {
            uint64_t *cf = &t[p->d].c[f->pos[k]];
            *cf = f->sgn[k] > 0 ? mod_add(*cf, 1) : mod_sub(*cf, 1);
        }
    }
}

/* v_l (L polys each) = B (sum_j H[l][j] x_j) for x in NTT form ((n+1) * mu) */
static void parity_first(poly *v, const tv_pub *pub, const poly *xh) {
    const tv_params *p = &pub->prm;
    int mu = p->mu;
    poly *comb = malloc(sizeof(poly) * mu), tmp;
    for (int l = 0; l <= p->n - p->t; l++) {
        for (int i = 0; i < mu; i++) poly_zero(&comb[i]);
        for (int j = 0; j <= p->n; j++) {
            uint64_t h = pub->H_m[l * (p->n + 1) + j];
            for (int i = 0; i < mu; i++) {
                poly_scale(&tmp, &xh[(size_t)j * mu + i], h);
                poly_add(&comb[i], &comb[i], &tmp);
            }
        }
        C_mul(&v[(size_t)l * p->L], pub, p->d, p->L, comb);
    }
    free(comb);
}

static i128 sqnorm(const int64_t *x, size_t n) {      /* exact */
    i128 s = 0;
    for (size_t i = 0; i < n; i++) s += (i128)x[i] * x[i];
    return s;
}

int tv_vote(tv_ballot *b, tv_voter_secret *sec, tv_prove_stats *stats, const tv_pub *pub,
            const tv_akey *keys, uint64_t id, const int *v, const uint8_t rnd_seed[32]) {
    const tv_params *p = &pub->prm;
    int n = p->n, L = p->L, d = p->d, mu = p->mu, rows = p->rows;
    size_t vn = VN(p);
    int wsum = 0;
    for (int a = 0; a < L; a++) { if (v[a] != 0 && v[a] != 1) return -1; wsum += v[a]; }
    if (p->w != TV_W_FREE && wsum != p->w) return -1;

    prg g;
    prg_init(&g, 0x01, rnd_seed, 32);
    gauss_sampler g1, gJ;
    gauss_init(&g1, p->sigma);
    gauss_init(&gJ, p->sigma_J);

    /* 1. shares */
    poly *vote = calloc((size_t)L, sizeof(poly));
    for (int a = 0; a < L; a++) vote[a].c[0] = (uint64_t)v[a];
    shamir_share(sec->m, vote, p, &g);
    free(vote);
    /* 2. randomness: r_0 fresh, r_k from seeds; commitments */
    b->id = id;
    gauss_vec(&g1, &g, sec->r, vn);
    for (int k = 1; k <= n; k++) {
        prg_bytes(&g, sec->seeds + (size_t)(k - 1) * 32, 32);
        tv_rand_from_seed(sec->r + (size_t)k * vn, pub, id, k, sec->seeds + (size_t)(k - 1) * 32);
    }
    for (int j = 0; j <= n; j++) tv_commit(&b->c[(size_t)j * rows], pub, &sec->m[(size_t)j * L], sec->r + (size_t)j * vn);
    /* 3. encryption of the seeds */
    for (int k = 1; k <= n; k++) pke_stub(b->e + (size_t)(k - 1) * TV_CT_BYTES, &keys[k - 1], id, sec->seeds + (size_t)(k - 1) * 32);

    /* 4. proof */
    uint8_t stdig[64];
    statement_digest(stdig, pub, b);
    size_t nz = (size_t)(n + 1) * vn;               /* z part */
    size_t nwit = nz + (size_t)L * vn;               /* witness-dependent responses */
    int64_t *y = malloc(sizeof(int64_t) * nz);
    int64_t *rho = malloc(sizeof(int64_t) * (size_t)L * vn);
    int64_t *zall = malloc(sizeof(int64_t) * nwit);  /* (z_0..z_n, r_{1,m_1}..r_{L,m_L}) */
    int64_t *shift = malloc(sizeof(int64_t) * nwit);
    poly *yh = malloc(sizeof(poly) * (size_t)(n + 1) * mu);
    poly *xh = malloc(sizeof(poly) * mu);
    poly *rows_tmp = malloc(sizeof(poly) * rows);
    challenge *fsim = malloc(sizeof(challenge) * L);
    perm_ch *pi = malloc(sizeof(perm_ch) * L);
    firstmsg fm;
    fm_alloc(&fm, p);
    int attempts = 0;
    double max_ratio = 0;
    for (;;) {
        attempts++;
        if (attempts > 1000) return -2;
        /* masks and first messages of the linear part */
        gauss_vec(&gJ, &g, y, nz);
        vec_ntt_from_int(yh, y, (n + 1) * mu);
        for (int j = 0; j <= n; j++) C_mul(&fm.w[(size_t)j * d], pub, 0, d, &yh[(size_t)j * mu]);
        parity_first(fm.v, pub, yh);
        if (p->w != TV_W_FREE) { poly_dot_mont(&fm.vw, pub->Bsum, yh, (size_t)mu); poly_invntt(&fm.vw); }
        /* OR first messages */
        for (int a = 0; a < L; a++) {
            int m = v[a], s = 1 - m;
            poly *ts = &fm.t[((size_t)a * 2 + s) * (d + 1)], *tm = &fm.t[((size_t)a * 2 + m) * (d + 1)];
            int64_t *rs = b->or_r + ((size_t)a * 2 + s) * vn;
            gauss_vec(&gJ, &g, rs, vn);                         /* simulated branch */
            sample_challenge(&fsim[a], &g);
            vec_ntt_from_int(xh, rs, mu);
            Ca_mul(ts, pub, a, xh);
            or_adjust(ts, p, b, a, &fsim[a], s);
            gauss_vec(&gJ, &g, rho + (size_t)a * vn, vn);       /* real branch */
            vec_ntt_from_int(xh, rho + (size_t)a * vn, mu);
            Ca_mul(tm, pub, a, xh);
        }
        fm_hash(b->chat, p, stdig, &fm);
        challenge c;
        expand_challenges(&c, pi, L, b->chat);
        /* responses and shift */
        ch_mul_int(shift, &c, sec->r, (size_t)(n + 1) * mu);
        for (size_t i = 0; i < nz; i++) zall[i] = y[i] + shift[i];
        for (int a = 0; a < L; a++) {
            int m = v[a];
            challenge fm_real;
            if (m == 1) { b->f0[a] = fsim[a]; perm_apply(&fm_real, &pi[a], &fsim[a]); }
            else { perm_apply_inv(&fm_real, &pi[a], &fsim[a]); b->f0[a] = fm_real; }
            int64_t *sh = shift + nz + (size_t)a * vn;
            ch_mul_int(sh, &fm_real, sec->r, (size_t)mu);             /* f_{a,m} r_0 */
            for (size_t i = 0; i < vn; i++) zall[nz + (size_t)a * vn + i] = rho[(size_t)a * vn + i] + sh[i];
        }
        /* single rejection step on all witness-dependent responses */
        i128 zv = 0, vv = 0;                      /* exact integers */
        for (size_t i = 0; i < nwit; i++) { zv += (i128)zall[i] * shift[i]; vv += (i128)shift[i] * shift[i]; }
        double ratio = sqrt((double)vv) / p->T;   /* statistics only */
        if (ratio > max_ratio) max_ratio = ratio;
        if (!reject_accept(&g, zv, vv, p->sigma_J2, p->logM_num, p->logM_den)) continue;
        /* assemble and check the norm bound (fails only with negligible probability) */
        memcpy(b->z, zall, sizeof(int64_t) * nz);
        for (int a = 0; a < L; a++) {
            int m = v[a];
            memcpy(b->or_r + ((size_t)a * 2 + m) * vn, zall + nz + (size_t)a * vn, sizeof(int64_t) * vn);
        }
        i128 nrm2 = sqnorm(b->z, nz) + sqnorm(b->or_r, (size_t)L * 2 * vn);
        if (nrm2 > (i128)p->B_J2) continue;
        if (stats) { stats->attempts = attempts; stats->max_shift_ratio = max_ratio; stats->resp_norm = sqrt((double)nrm2); }
        break;
    }
    fm_free(&fm);
    free(y); free(rho); free(zall); free(shift); free(yh); free(xh); free(rows_tmp); free(fsim); free(pi);
    gauss_free(&g1);
    gauss_free(&gJ);
    return 0;
}

int tv_verify_ballot(const tv_pub *pub, const tv_ballot *b) {
    const tv_params *p = &pub->prm;
    int n = p->n, L = p->L, d = p->d, mu = p->mu, rows = p->rows;
    size_t vn = VN(p), nz = (size_t)(n + 1) * vn;
    /* norm bound */
    i128 nrm2 = sqnorm(b->z, nz) + sqnorm(b->or_r, (size_t)L * 2 * vn);
    if (nrm2 > (i128)p->B_J2) return 0;
    challenge c;
    perm_ch *pi = malloc(sizeof(perm_ch) * L);
    expand_challenges(&c, pi, L, b->chat);
    firstmsg fm;
    fm_alloc(&fm, p);
    poly *zh = malloc(sizeof(poly) * (size_t)(n + 1) * mu), tmp;
    poly *xh = malloc(sizeof(poly) * mu);
    vec_ntt_from_int(zh, b->z, (n + 1) * mu);
    /* w_j = A z_j - c c_{j,1} */
    for (int j = 0; j <= n; j++) {
        C_mul(&fm.w[(size_t)j * d], pub, 0, d, &zh[(size_t)j * mu]);
        for (int i = 0; i < d; i++) {
            poly_mul_sparse(&tmp, &b->c[(size_t)j * rows + i], c.pos, c.sgn, TV_CH_W);
            poly_sub(&fm.w[(size_t)j * d + i], &fm.w[(size_t)j * d + i], &tmp);
        }
    }
    /* v_l = B (sum_j H_lj z_j) - c sum_j H_lj c_{j,2} */
    parity_first(fm.v, pub, zh);
    for (int l = 0; l <= n - p->t; l++) {
        for (int a = 0; a < L; a++) {
            poly acc, t2;
            poly_zero(&acc);
            for (int j = 0; j <= n; j++) {
                poly_scale(&t2, &b->c[(size_t)j * rows + d + a], pub->H_m[l * (n + 1) + j]);
                poly_add(&acc, &acc, &t2);
            }
            poly_mul_sparse(&tmp, &acc, c.pos, c.sgn, TV_CH_W);
            poly_sub(&fm.v[(size_t)l * L + a], &fm.v[(size_t)l * L + a], &tmp);
        }
    }
    /* v_w = Bsum z_0 - c sum_a c_{0,2,a} + c w */
    if (p->w != TV_W_FREE) {
        poly_dot_mont(&fm.vw, pub->Bsum, zh, (size_t)mu);
        poly_invntt(&fm.vw);
        poly acc;
        poly_zero(&acc);
        for (int a = 0; a < L; a++) poly_add(&acc, &acc, &b->c[d + a]);
        acc.c[0] = mod_sub(acc.c[0], (uint64_t)p->w);
        poly_mul_sparse(&tmp, &acc, c.pos, c.sgn, TV_CH_W);
        poly_sub(&fm.vw, &fm.vw, &tmp);
    }
    /* OR first messages */
    for (int a = 0; a < L; a++) {
        challenge f1;
        perm_apply(&f1, &pi[a], &b->f0[a]);
        const challenge *f[2] = {&b->f0[a], &f1};
        for (int bit = 0; bit < 2; bit++) {
            poly *t = &fm.t[((size_t)a * 2 + bit) * (d + 1)];
            vec_ntt_from_int(xh, b->or_r + ((size_t)a * 2 + bit) * vn, mu);
            Ca_mul(t, pub, a, xh);
            or_adjust(t, p, b, a, f[bit], bit);
        }
    }
    uint8_t stdig[64], chat[TV_HASHBYTES];
    statement_digest(stdig, pub, b);
    fm_hash(chat, p, stdig, &fm);
    int ok = memcmp(chat, b->chat, TV_HASHBYTES) == 0;
    fm_free(&fm);
    free(pi); free(zh); free(xh);
    return ok;
}

int tv_check_share(const tv_pub *pub, const tv_akey *key, int k, const tv_ballot *b, poly *m_k, int64_t *r_k) {
    const tv_params *p = &pub->prm;
    uint8_t seed[32];
    pke_stub(seed, key, b->id, b->e + (size_t)(k - 1) * TV_CT_BYTES);
    tv_rand_from_seed(r_k, pub, b->id, k, seed);
    if (sqnorm(r_k, VN(p)) > (i128)tv_beta2(p)) return 0;
    poly *cr = malloc(sizeof(poly) * p->rows);
    poly *rh = malloc(sizeof(poly) * p->mu);
    vec_ntt_from_int(rh, r_k, p->mu);
    C_mul(cr, pub, 0, p->rows, rh);
    int ok = 1;
    const poly *ck = &b->c[(size_t)k * p->rows];
    for (int i = 0; i < p->d; i++) ok &= memcmp(&cr[i], &ck[i], sizeof(poly)) == 0;
    for (int a = 0; a < p->L; a++) poly_sub(&m_k[a], &ck[p->d + a], &cr[p->d + a]);
    free(cr); free(rh);
    return ok;
}

/* ------------------------------------------------------------------ serialization */

size_t tv_ballot_maxbytes(const tv_params *p) {
    size_t vn = VN(p);
    size_t coms = (size_t)(p->n + 1) * p->rows * TV_N * 8;
    size_t resp = ((size_t)(p->n + 1) + 2 * (size_t)p->L) * vn * 16;   /* generous for Rice codes */
    return 64 + coms + (size_t)p->n * TV_CT_BYTES + resp + (size_t)p->L * 80;
}

static void put_ch(bitw *w, const challenge *c) {
    for (int k = 0; k < TV_CH_W; k++) bw_put(w, c->pos[k], 8);
    for (int k = 0; k < TV_CH_W; k++) bw_put(w, c->sgn[k] < 0, 1);
}
static int get_ch(bitr *r, challenge *c) {
    int last = -1;
    for (int k = 0; k < TV_CH_W; k++) {
        c->pos[k] = (uint8_t)br_get(r, 8);
        if ((int)c->pos[k] <= last) return -1;            /* canonical: strictly increasing */
        last = c->pos[k];
    }
    for (int k = 0; k < TV_CH_W; k++) c->sgn[k] = br_get(r, 1) ? -1 : 1;
    return 0;
}

size_t tv_ballot_encode(uint8_t *buf, size_t cap, const tv_pub *pub, const tv_ballot *b, tv_sizes *sz) {
    const tv_params *p = &pub->prm;
    size_t vn = VN(p);
    int k = rice_param(p->sigma_J);
    bitw w;
    bw_init(&w, buf, cap);
    bw_put(&w, b->id & 0xFFFFFFFFULL, 32);
    bw_put(&w, b->id >> 32, 32);
    size_t s0 = w.bytes;
    for (size_t i = 0; i < (size_t)(p->n + 1) * p->rows; i++) pack_poly(&w, &b->c[i]);
    size_t s1 = w.bytes;
    for (size_t i = 0; i < (size_t)p->n * TV_CT_BYTES; i++) bw_put(&w, b->e[i], 8);
    size_t s2 = w.bytes;
    for (int i = 0; i < TV_HASHBYTES; i++) bw_put(&w, b->chat[i], 8);
    for (size_t i = 0; i < (size_t)(p->n + 1) * vn; i++) rice_put(&w, b->z[i], k);
    for (int a = 0; a < p->L; a++) {
        put_ch(&w, &b->f0[a]);
        for (size_t i = 0; i < 2 * vn; i++) rice_put(&w, b->or_r[(size_t)a * 2 * vn + i], k);
    }
    size_t tot = bw_finish(&w);
    if (sz) { sz->total = tot; sz->commitments = s1 - s0; sz->ciphertexts = s2 - s1; sz->proof = tot - s2; }
    return w.err ? 0 : tot;
}

int tv_ballot_decode(tv_ballot *b, const tv_pub *pub, const uint8_t *buf, size_t len) {
    const tv_params *p = &pub->prm;
    size_t vn = VN(p);
    int k = rice_param(p->sigma_J);
    bitr r;
    br_init(&r, buf, len);
    b->id = br_get(&r, 32);
    b->id |= br_get(&r, 32) << 32;
    for (size_t i = 0; i < (size_t)(p->n + 1) * p->rows; i++) unpack_poly(&r, &b->c[i]);
    for (size_t i = 0; i < (size_t)p->n * TV_CT_BYTES; i++) b->e[i] = (uint8_t)br_get(&r, 8);
    for (int i = 0; i < TV_HASHBYTES; i++) b->chat[i] = (uint8_t)br_get(&r, 8);
    for (size_t i = 0; i < (size_t)(p->n + 1) * vn; i++) b->z[i] = rice_get(&r, k);
    for (int a = 0; a < p->L; a++) {
        if (get_ch(&r, &b->f0[a])) return -1;
        for (size_t i = 0; i < 2 * vn; i++) b->or_r[(size_t)a * 2 * vn + i] = rice_get(&r, k);
    }
    /* canonical encoding: every byte consumed and the padding bits of the last byte zero */
    if (r.pos != len || r.acc != 0) return -1;
    return r.err ? -1 : 0;
}
