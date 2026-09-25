#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "agg.h"
#include "codec.h"
#include "timer.h"

#define VN(p) ((size_t)(p)->mu * TV_N)

/* ------------------------------------------------------------------ tree shape and blocks */

int agg_tree_shape(agg_tree *t, int nleaves, int fanin) {
    memset(t, 0, sizeof *t);
    t->nleaves = nleaves;
    int cnt[64], lv = 0, tot = 0;
    cnt[0] = nleaves;
    tot = nleaves;
    while (cnt[lv] > 1 || lv == 0) {
        cnt[lv + 1] = (cnt[lv] + fanin - 1) / fanin;
        lv++;
        tot += cnt[lv];
    }
    t->nlevels = lv + 1;
    t->nnodes = tot;
    t->lvl_start = malloc(sizeof(int) * (size_t)t->nlevels);
    t->lvl_count = malloc(sizeof(int) * (size_t)t->nlevels);
    t->child_start = malloc(sizeof(int) * (size_t)tot);
    t->nchild = malloc(sizeof(int) * (size_t)tot);
    int s = 0;
    for (int l = 0; l < t->nlevels; l++) { t->lvl_start[l] = s; t->lvl_count[l] = cnt[l]; s += cnt[l]; }
    for (int u = 0; u < tot; u++) { t->child_start[u] = -1; t->nchild[u] = 0; }
    for (int l = 1; l < t->nlevels; l++)
        for (int i = 0; i < cnt[l]; i++) {
            int u = t->lvl_start[l] + i;
            int first = i * fanin, last = first + fanin;
            if (last > cnt[l - 1]) last = cnt[l - 1];
            t->child_start[u] = t->lvl_start[l - 1] + first;
            t->nchild[u] = last - first;
        }
    return 0;
}

int agg_make_blocks(agg_block **blocks, const agg_tree *t, int k) {
    int E = (t->nnodes + k - 1) / k;
    for (int l = 1; l < t->nlevels; l++) E += (t->lvl_count[l] + k - 1) / k;
    agg_block *b = malloc(sizeof(agg_block) * (size_t)E);
    int e = 0;
    for (int f = 0; f < t->nnodes; f += k) { b[e].zero = 0; b[e].level = -1; b[e].first = f; b[e].count = (t->nnodes - f < k) ? t->nnodes - f : k; e++; }
    for (int l = 1; l < t->nlevels; l++)
        for (int i = 0; i < t->lvl_count[l]; i += k) {
            b[e].zero = 1; b[e].level = l; b[e].first = t->lvl_start[l] + i;
            b[e].count = (t->lvl_count[l] - i < k) ? t->lvl_count[l] - i : k;
            e++;
        }
    *blocks = b;
    return E;
}

static void tree_free(agg_tree *t) {
    free(t->lvl_start); free(t->lvl_count); free(t->child_start); free(t->nchild);
    free(t->com); free(t->msg); free(t->rnd);
    memset(t, 0, sizeof *t);
}

/* ------------------------------------------------------------------ parameters */

/* HKZ tail bound for ||Phi x||^2, x sigma-subgaussian: tr, lam = ||Phi^T Phi||, tau in nats */
static double hkz(double tr, double lam, double tau) { return tr + 2.0 * sqrt(tr * lam * tau) + 2.0 * lam * tau; }

void agg_params_init(agg_params *ap, const tv_params *p, long nleaves, c2_mode mode, double logQ) {
    memset(ap, 0, sizeof *ap);
    ap->fanin = p->fanin; ap->k = p->k_blk; ap->ell = p->ell; ap->kappa = p->kappa;
    ap->mode = mode;
    ap->alpha = 12.0 / log(sqrt(3.0));
    ap->logM = 12.0 / ap->alpha + 1.0 / (2.0 * ap->alpha * ap->alpha);
    ap->zinf_mult = 9.0;
    ap->logQ = logQ;
    ap->tau = (128.0 + 40.0) * log(2.0);
    agg_tree t;
    agg_tree_shape(&t, (int)nleaves, ap->fanin);
    agg_block *blk;
    int E = agg_make_blocks(&blk, &t, ap->k);
    ap->E = E; ap->nnodes = t.nnodes; ap->ninner = t.nnodes - t.nleaves;
    for (int e = 0; e < E; e++) { if (blk[e].zero) ap->E_zero++; else ap->E_open++; }
    double Nmu = (double)TV_N * p->mu, s2 = p->sigma * p->sigma;
    double ncols = (double)t.nnodes + ap->ninner;
    double lnQ = logQ * log(2.0);
    double T1sq = 0, Fsum = 0, Fsq = 0, Fmax = 0, T2bin = 0;
    for (int e = 0; e < E; e++) {
        double Fe = 0;
        for (int i = 0; i < blk[e].count; i++) {
            int u = blk[e].first + i;
            double v, g;                                   /* variance multiplicity, grindable summands */
            if (!blk[e].zero) { v = 1; g = (u < t.nleaves) ? 1 : 0; }
            else { v = 1 + t.nchild[u]; g = (blk[e].level == 1) ? t.nchild[u] : 0; }
            double tc = ap->tau + log(ncols) + g * lnQ;
            /* c1 = H_agg(par, BB_1) can be steered by the adversary among its queries: one more ln Q */
            T1sq += s2 * v * hkz(60.0 * Nmu, 3600.0, tc + lnQ);    /* ||c1 s_col||^2 */
            Fe += s2 * v * hkz(Nmu, 1.0, tc);                      /* ||s_col||^2 */
        }
        Fsum += Fe; Fsq += Fe * Fe; if (Fe > Fmax) Fmax = Fe;
        /* binary C2: ||S C2|| <= ||S||_F s1(C2), s1(C2) <= sqrt(k l)/2 + s1(Rademacher)/2 (net bound, eps=0.1) */
        double ke = blk[e].count, eps = 0.1;
        double srad = sqrt(2.0 * ((ke + ap->ell) * log(1 + 2 / eps) + ap->tau + log((double)E))) / (1 - 2 * eps);
        double s1c = sqrt(ke * ap->ell) / 2 + srad / 2;
        T2bin += Fe * s1c * s1c;
    }
    ap->T1 = sqrt(T1sq);
    if (mode == C2_SIGNED) ap->T2 = sqrt(ap->ell * Fsum + 2.0 * sqrt(ap->ell * Fsq * ap->tau) + 2.0 * Fmax * ap->tau);
    else ap->T2 = sqrt(T2bin);
    ap->sigma1 = ap->alpha * ap->T1;
    ap->sigma2 = ap->alpha * ap->T2;
    free(blk);
    tree_free(&t);
}

void agg_params_print(const agg_params *ap, const tv_params *p) {
    (void)p;
    printf("aggregation: fanin=%d k=%d ell=%d kappa=%d C2=%s logQ=%.0f E=%d (open %d, zero %d) nodes=%d\n",
           ap->fanin, ap->k, ap->ell, ap->kappa, ap->mode == C2_SIGNED ? "signed" : "binary", ap->logQ,
           ap->E, ap->E_open, ap->E_zero, ap->nnodes);
    printf("  T1=2^%.2f T2=2^%.2f sigma1=2^%.2f sigma2=2^%.2f\n", log2(ap->T1), log2(ap->T2), log2(ap->sigma1), log2(ap->sigma2));
}

/* ------------------------------------------------------------------ helpers */

static void seed_ctx(uint8_t *buf, const uint8_t seed[32], int a, int b, int c) {
    memcpy(buf, seed, 32);
    int v[3] = {a, b, c};
    for (int j = 0; j < 3; j++) for (int i = 0; i < 4; i++) buf[32 + 4 * j + i] = (uint8_t)(v[j] >> (8 * i));
}

/* masks of block e in attempt with seed: Y1 (mu*N x count) and Y2 (mu*N x ell), column-major */
static void gen_masks(int64_t *Y1, int64_t *Y2, const uint8_t seed[32], int e, int count, size_t vn, int ell,
                      gauss_sampler *g1, gauss_sampler *g2) {
    uint8_t buf[44];
    prg g;
    seed_ctx(buf, seed, e, 1, 0);
    prg_init(&g, 0x40, buf, 44);
    gauss_vec(g1, &g, Y1, vn * (size_t)count);
    seed_ctx(buf, seed, e, 2, 0);
    prg_init(&g, 0x40, buf, 44);
    gauss_vec(g2, &g, Y2, vn * (size_t)ell);
}

static void first_challenge(challenge *c1, const uint8_t bb1[TV_HASHBYTES], int k, int iota, int e) {
    uint8_t buf[44];
    prg g;
    seed_ctx(buf, bb1, k, iota, e);
    prg_init(&g, 0x41, buf, 44);
    sample_challenge(c1, &g);
}

/* C2 = H(Z1, c1, context): count x ell entries, binary {0,1} or signed {-1,1} */
static void second_challenge(int8_t *C2, const int64_t *Z1, size_t nz1, const challenge *c1, const uint8_t bb1[TV_HASHBYTES],
                             int k, int iota, int e, int count, int ell, c2_mode mode) {
    keccak_state st;
    shake256_init(&st);
    uint8_t buf[44];
    seed_ctx(buf, bb1, k, iota, e);
    keccak_absorb(&st, (const uint8_t *)"T-EVOLVE/agg-H", 14);
    keccak_absorb(&st, buf, 44);
    keccak_absorb(&st, c1->pos, TV_CH_W);
    keccak_absorb(&st, (const uint8_t *)c1->sgn, TV_CH_W);
    uint8_t tmp[4096];
    size_t fill = 0;
    for (size_t i = 0; i < nz1; i++) {
        uint64_t v = (uint64_t)Z1[i];
        for (int b = 0; b < 8; b++) tmp[fill++] = (uint8_t)(v >> (8 * b));
        if (fill == sizeof tmp) { keccak_absorb(&st, tmp, fill); fill = 0; }
    }
    keccak_absorb(&st, tmp, fill);
    keccak_finalize(&st);
    size_t nbits = (size_t)count * ell, nbytes = (nbits + 7) / 8;
    uint8_t *bits = malloc(nbytes);
    keccak_squeeze(&st, bits, nbytes);
    for (size_t i = 0; i < nbits; i++) {
        int bit = (bits[i >> 3] >> (i & 7)) & 1;
        C2[i] = (int8_t)(mode == C2_SIGNED ? (bit ? -1 : 1) : bit);
    }
    free(bits);
}

/* witness column of node u in block type zero/open (mu*N ints) */
static void witness_col(int64_t *out, const agg_tree *t, size_t vn, int u, int zero) {
    const int8_t *r = t->rnd + (size_t)u * vn;
    for (size_t i = 0; i < vn; i++) out[i] = r[i];
    if (zero)
        for (int c = 0; c < t->nchild[u]; c++) {
            const int8_t *rc = t->rnd + (size_t)(t->child_start[u] + c) * vn;
            for (size_t i = 0; i < vn; i++) out[i] -= rc[i];
        }
}

/* statement column of node u: open -> first d rows of d_u; zero -> delta_u (rows) */
static void statement_col(poly *out, const poly *com_of_node(const void *, int), const void *ctx, const tv_params *p,
                          const agg_tree *t, int u, int zero) {
    int mrows = zero ? p->rows : p->d;
    const poly *du = com_of_node(ctx, u);
    for (int i = 0; i < mrows; i++) out[i] = du[i];
    if (zero)
        for (int c = 0; c < t->nchild[u]; c++) {
            const poly *dc = com_of_node(ctx, t->child_start[u] + c);
            for (int i = 0; i < mrows; i++) poly_sub(&out[i], &out[i], &dc[i]);
        }
}

/* out = M x for an integer column x (mu polys), M = rows 0..mrows-1 of C */
static void M_mul_int(poly *out, const tv_pub *pub, int mrows, const int64_t *x, poly *scratch) {
    vec_ntt_from_int(scratch, x, pub->prm.mu);
    C_mul(out, pub, 0, mrows, scratch);
}

/* absorb a set of polynomials into G */
static void g_absorb(keccak_state *st, const poly *w, size_t n) { absorb_polys(st, w, n); }

/* ------------------------------------------------------------------ round 1 */

typedef struct { const poly *com; int rows; } com_view;
static const poly *view_com(const void *ctx, int u) { const com_view *v = ctx; return v->com + (size_t)u * v->rows; }

/* first messages of one block for masks Y1, Y2, absorbed into G */
static void block_first_messages(keccak_state *G, const tv_pub *pub, const agg_block *b, const int64_t *Y1, const int64_t *Y2,
                                 int ell, poly *wbuf, poly *scratch) {
    const tv_params *p = &pub->prm;
    size_t vn = VN(p);
    int mrows = b->zero ? p->rows : p->d;
    for (int i = 0; i < b->count; i++) { M_mul_int(wbuf, pub, mrows, Y1 + (size_t)i * vn, scratch); g_absorb(G, wbuf, (size_t)mrows); }
    for (int j = 0; j < ell; j++) { M_mul_int(wbuf, pub, mrows, Y2 + (size_t)j * vn, scratch); g_absorb(G, wbuf, (size_t)mrows); }
}

static void put_le32(uint8_t *b, int v) { for (int i = 0; i < 4; i++) b[i] = (uint8_t)((uint32_t)v >> (8 * i)); }

static void g_init(keccak_state *G, int k, int iota, const uint8_t leaf_dig[TV_HASHBYTES]) {
    shake256_init(G);
    keccak_absorb(G, (const uint8_t *)"T-EVOLVE/agg-G", 14);
    uint8_t v[8];
    put_le32(v, k);
    put_le32(v + 4, iota);
    keccak_absorb(G, v, 8);
    keccak_absorb(G, leaf_dig, TV_HASHBYTES);
}

static void leaf_digest(uint8_t out[TV_HASHBYTES], const tv_pub *pub, const poly *leaf_com, int nleaves) {
    keccak_state st;
    shake256_init(&st);
    keccak_absorb(&st, (const uint8_t *)"T-EVOLVE/leaves", 15);
    uint8_t v[4];
    put_le32(v, nleaves);
    keccak_absorb(&st, v, 4);
    absorb_polys(&st, leaf_com, (size_t)nleaves * pub->prm.rows);
    keccak_finalize(&st);
    keccak_squeeze(&st, out, TV_HASHBYTES);
}

int agg_round1(agg_state *st, agg_contrib *out, agg_stats *stats, const tv_pub *pub, const agg_params *ap, int k,
               const poly *leaf_com, const poly *leaf_msg, const int64_t *leaf_rnd, int nleaves, const uint8_t seed[32]) {
    double t0 = now_ms();
    const tv_params *p = &pub->prm;
    size_t vn = VN(p);
    int rows = p->rows, L = p->L;
    memset(st, 0, sizeof *st);
    memset(out, 0, sizeof *out);
    prg_init(&st->rng, 0x42, seed, 32);
    agg_tree *t = &st->tree;
    agg_tree_shape(t, nleaves, ap->fanin);
    t->com = malloc(sizeof(poly) * (size_t)t->nnodes * rows);
    t->msg = malloc(sizeof(poly) * (size_t)t->nnodes * L);
    t->rnd = malloc((size_t)t->nnodes * vn);
    memcpy(t->com, leaf_com, sizeof(poly) * (size_t)nleaves * rows);
    memcpy(t->msg, leaf_msg, sizeof(poly) * (size_t)nleaves * L);
    for (size_t i = 0; i < (size_t)nleaves * vn; i++) {
        if (leaf_rnd[i] < -127 || leaf_rnd[i] > 127) return -1;
        t->rnd[i] = (int8_t)leaf_rnd[i];
    }
    gauss_sampler g1;
    gauss_init(&g1, p->sigma);
    int64_t *rho = malloc(sizeof(int64_t) * vn);
    for (int u = nleaves; u < t->nnodes; u++) {
        poly *m = t->msg + (size_t)u * L;
        for (int a = 0; a < L; a++) poly_zero(&m[a]);
        for (int c = 0; c < t->nchild[u]; c++)
            for (int a = 0; a < L; a++) poly_add(&m[a], &m[a], &t->msg[(size_t)(t->child_start[u] + c) * L + a]);
        gauss_vec(&g1, &st->rng, rho, vn);
        for (size_t i = 0; i < vn; i++) {
            if (rho[i] < -127 || rho[i] > 127) return -1;
            t->rnd[(size_t)u * vn + i] = (int8_t)rho[i];
        }
        tv_commit(t->com + (size_t)u * rows, pub, m, rho);
    }
    gauss_free(&g1);
    free(rho);
    int E = agg_make_blocks(&st->blocks, t, ap->k);
    if (E != ap->E) return -1;
    /* publish */
    out->k_auth = k;
    leaf_digest(out->leaf_digest, pub, leaf_com, nleaves);
    out->ninner = t->nnodes - nleaves;
    out->inner_com = malloc(sizeof(poly) * (size_t)out->ninner * rows);
    memcpy(out->inner_com, t->com + (size_t)nleaves * rows, sizeof(poly) * (size_t)out->ninner * rows);
    out->V = malloc(sizeof(poly) * L);
    memcpy(out->V, t->msg + (size_t)(t->nnodes - 1) * L, sizeof(poly) * L);
    out->h = malloc(sizeof(*out->h) * (size_t)ap->kappa);
    out->iota = -1;
    st->seed = malloc(sizeof(*st->seed) * (size_t)ap->kappa);
    st->salt = malloc(sizeof(*st->salt) * (size_t)ap->kappa);
    /* kappa attempts: first messages committed by h = G(W, salt) */
    gauss_sampler gs1, gs2;
    gauss_init(&gs1, ap->sigma1);
    gauss_init(&gs2, ap->sigma2);
    int64_t *Y1 = malloc(sizeof(int64_t) * vn * (size_t)ap->k);
    int64_t *Y2 = malloc(sizeof(int64_t) * vn * (size_t)ap->ell);
    poly *wbuf = malloc(sizeof(poly) * rows), *scratch = malloc(sizeof(poly) * p->mu);
    for (int it = 0; it < ap->kappa; it++) {
        prg_bytes(&st->rng, st->seed[it], 32);
        prg_bytes(&st->rng, st->salt[it], 32);
        keccak_state G;
        g_init(&G, k, it, out->leaf_digest);
        for (int e = 0; e < E; e++) {
            gen_masks(Y1, Y2, st->seed[it], e, st->blocks[e].count, vn, ap->ell, &gs1, &gs2);
            block_first_messages(&G, pub, &st->blocks[e], Y1, Y2, ap->ell, wbuf, scratch);
        }
        keccak_absorb(&G, st->salt[it], 32);
        keccak_finalize(&G);
        keccak_squeeze(&G, out->h[it], TV_HASHBYTES);
    }
    free(Y1); free(Y2); free(wbuf); free(scratch);
    gauss_free(&gs1); gauss_free(&gs2);
    if (stats) {
        memset(stats, 0, sizeof *stats);
        stats->t_round1_ms = now_ms() - t0;
        stats->bytes_round1 = agg_contrib_bytes(out, pub, ap, 1);
    }
    return 0;
}

void agg_bb1_digest(uint8_t out[TV_HASHBYTES], const tv_pub *pub, const agg_contrib *c, int nc) {
    keccak_state st;
    shake256_init(&st);
    keccak_absorb(&st, (const uint8_t *)"T-EVOLVE/BB1", 12);
    keccak_absorb(&st, pub->digest, TV_HASHBYTES);
    for (int i = 0; i < nc; i++) {
        uint8_t kb[4];
        put_le32(kb, c[i].k_auth);
        keccak_absorb(&st, kb, 4);
        keccak_absorb(&st, c[i].leaf_digest, TV_HASHBYTES);
        absorb_polys(&st, c[i].inner_com, (size_t)c[i].ninner * pub->prm.rows);
        absorb_polys(&st, c[i].V, (size_t)pub->prm.L);
        keccak_absorb(&st, (const uint8_t *)c[i].h, sizeof(*c[i].h) * (size_t)pub->prm.kappa);
    }
    keccak_finalize(&st);
    keccak_squeeze(&st, out, TV_HASHBYTES);
}

/* ------------------------------------------------------------------ round 2 */

static size_t z1_offset(const agg_block *blk, int e, size_t vn) {
    size_t off = 0;
    for (int i = 0; i < e; i++) off += (size_t)blk[i].count * vn;
    return off;
}

int agg_round2(agg_state *st, agg_contrib *out, agg_stats *stats, const tv_pub *pub, const agg_params *ap,
               const uint8_t bb1[TV_HASHBYTES]) {
    double t0 = now_ms();
    const tv_params *p = &pub->prm;
    size_t vn = VN(p);
    agg_tree *t = &st->tree;
    int E = ap->E, ell = ap->ell, k = out->k_auth;
    size_t nz1 = (size_t)t->nnodes * vn + (size_t)(t->nnodes - t->nleaves) * vn;   /* all columns of both proofs */
    free(out->Z1); free(out->Z2);
    out->Z1 = malloc(sizeof(int64_t) * nz1);
    out->Z2 = malloc(sizeof(int64_t) * (size_t)E * vn * ell);
    gauss_sampler gs1, gs2;
    gauss_init(&gs1, ap->sigma1);
    gauss_init(&gs2, ap->sigma2);
    int64_t *Y1 = malloc(sizeof(int64_t) * vn * (size_t)ap->k), *Y2 = malloc(sizeof(int64_t) * vn * (size_t)ell);
    int64_t *S = malloc(sizeof(int64_t) * vn * (size_t)ap->k), *B = malloc(sizeof(int64_t) * vn);
    int64_t *acc = malloc(sizeof(int64_t) * vn);
    int8_t *C2 = malloc((size_t)ap->k * ell);
    double lim_row = 2.0 * ap->k * ap->sigma1 * ap->sigma1;
    double lim_ent = ap->zinf_mult * ap->sigma2, lim_col = 2.0 * TV_N * ap->sigma2 * ap->sigma2;
    double *rowsq = malloc(sizeof(double) * vn);
    double r1max = 0, r2max = 0;
    int found = -1, used = 0;
    for (int it = 0; it < ap->kappa && found < 0; it++) {
        used++;
        long double zb1 = 0, bb1s = 0, zb2 = 0, bb2s = 0;
        int small = 1;
        for (int e = 0; e < E; e++) {
            const agg_block *b = &st->blocks[e];
            gen_masks(Y1, Y2, st->seed[it], e, b->count, vn, ell, &gs1, &gs2);
            challenge c1;
            first_challenge(&c1, bb1, k, it, e);
            int64_t *Z1 = out->Z1 + z1_offset(st->blocks, e, vn);
            for (size_t r = 0; r < vn; r++) rowsq[r] = 0;
            for (int i = 0; i < b->count; i++) {
                int64_t *s = S + (size_t)i * vn;
                witness_col(s, t, vn, b->first + i, b->zero);
                ch_mul_int(B, &c1, s, (size_t)p->mu);
                for (size_t r = 0; r < vn; r++) {
                    int64_t z = Y1[(size_t)i * vn + r] + B[r];
                    zb1 += (long double)z * B[r];
                    bb1s += (long double)B[r] * B[r];
                    Z1[(size_t)i * vn + r] = (int64_t)z;
                    rowsq[r] += (double)z * (double)z;
                }
            }
            for (size_t r = 0; r < vn; r++) if (rowsq[r] > lim_row) small = 0;
            second_challenge(C2, Z1, (size_t)b->count * vn, &c1, bb1, k, it, e, b->count, ell, ap->mode);
            int64_t *Z2 = out->Z2 + (size_t)e * vn * ell;
            for (int j = 0; j < ell; j++) {
                for (size_t r = 0; r < vn; r++) acc[r] = 0;
                for (int i = 0; i < b->count; i++) {
                    int8_t cij = C2[(size_t)i * ell + j];
                    if (!cij) continue;
                    const int64_t *s = S + (size_t)i * vn;
                    if (cij > 0) for (size_t r = 0; r < vn; r++) acc[r] += s[r];
                    else for (size_t r = 0; r < vn; r++) acc[r] -= s[r];
                }
                for (size_t r = 0; r < vn; r++) {
                    int64_t z = Y2[(size_t)j * vn + r] + acc[r];
                    zb2 += (long double)z * acc[r];
                    bb2s += (long double)acc[r] * acc[r];
                    if (fabs((double)z) > lim_ent) small = 0;
                    Z2[(size_t)j * vn + r] = (int64_t)z;
                }
                for (int pp = 0; pp < p->mu; pp++) {
                    double cs = 0;
                    for (int c = 0; c < TV_N; c++) { double z = Z2[(size_t)j * vn + (size_t)pp * TV_N + c]; cs += z * z; }
                    if (cs > lim_col) small = 0;
                }
            }
        }
        double r1 = sqrt((double)bb1s) / ap->T1, r2 = sqrt((double)bb2s) / ap->T2;
        if (r1 > r1max) r1max = r1;
        if (r2 > r2max) r2max = r2;
        double l1 = (double)((-2.0L * zb1 + bb1s) / (2.0L * (long double)ap->sigma1 * ap->sigma1)) - ap->logM;
        double l2 = (double)((-2.0L * zb2 + bb2s) / (2.0L * (long double)ap->sigma2 * ap->sigma2)) - ap->logM;
        int u1 = log(prg_unif(&st->rng)) <= l1, u2 = log(prg_unif(&st->rng)) <= l2;
        if (u1 && u2 && small) found = it;
    }
    out->iota = found;
    if (found < 0) {                   /* responses of rejected attempts must never be published */
        free(out->Z1); free(out->Z2);
        out->Z1 = out->Z2 = NULL;
    }
    if (found >= 0) {
        memcpy(out->salt, st->salt[found], 32);
        out->rho_root = malloc(sizeof(int64_t) * vn);
        const int8_t *rr = t->rnd + (size_t)(t->nnodes - 1) * vn;
        for (size_t i = 0; i < vn; i++) out->rho_root[i] = rr[i];
    }
    free(Y1); free(Y2); free(S); free(B); free(acc); free(C2); free(rowsq);
    gauss_free(&gs1); gauss_free(&gs2);
    if (stats) {
        stats->attempts_used = used;
        stats->shift1_ratio = r1max; stats->shift2_ratio = r2max;
        stats->t_round2_ms = now_ms() - t0;
        stats->bytes_round2 = found >= 0 ? agg_contrib_bytes(out, pub, ap, 2) : 1;
    }
    return found >= 0 ? 0 : 1;
}

/* ------------------------------------------------------------------ verification */

int agg_verify(const tv_pub *pub, const agg_params *ap, const agg_contrib *c, const poly *leaf_com, int nleaves,
               const uint8_t bb1[TV_HASHBYTES]) {
    const tv_params *p = &pub->prm;
    size_t vn = VN(p);
    int rows = p->rows, L = p->L, ell = ap->ell;
    if (c->iota < 0 || c->iota >= ap->kappa) return 0;
    uint8_t ld[TV_HASHBYTES];
    leaf_digest(ld, pub, leaf_com, nleaves);
    if (memcmp(ld, c->leaf_digest, TV_HASHBYTES) != 0) return 0;
    agg_tree t;
    agg_tree_shape(&t, nleaves, ap->fanin);
    if (t.nnodes - nleaves != c->ninner) { tree_free(&t); return 0; }
    agg_block *blk;
    int E = agg_make_blocks(&blk, &t, ap->k);
    /* all node commitments */
    poly *com = malloc(sizeof(poly) * (size_t)t.nnodes * rows);
    memcpy(com, leaf_com, sizeof(poly) * (size_t)nleaves * rows);
    memcpy(com + (size_t)nleaves * rows, c->inner_com, sizeof(poly) * (size_t)c->ninner * rows);
    com_view cv = {com, rows};
    int ok = (E == ap->E);
    keccak_state G;
    g_init(&G, c->k_auth, c->iota, c->leaf_digest);
    poly *T = malloc(sizeof(poly) * (size_t)ap->k * rows), *w = malloc(sizeof(poly) * rows), *scratch = malloc(sizeof(poly) * p->mu);
    int8_t *C2 = malloc((size_t)ap->k * ell);
    int64_t *col = malloc(sizeof(int64_t) * vn);
    int64_t *accq = malloc(sizeof(int64_t) * (size_t)rows * TV_N);
    double lim_row = 2.0 * ap->k * ap->sigma1 * ap->sigma1;
    double lim_ent = ap->zinf_mult * ap->sigma2, lim_col = 2.0 * TV_N * ap->sigma2 * ap->sigma2;
    for (int e = 0; e < E && ok; e++) {
        const agg_block *b = &blk[e];
        int mrows = b->zero ? rows : p->d;
        const int64_t *Z1 = c->Z1 + z1_offset(blk, e, vn);
        const int64_t *Z2 = c->Z2 + (size_t)e * vn * ell;
        /* IsSmall */
        for (size_t r = 0; r < vn && ok; r++) {
            double s = 0;
            for (int i = 0; i < b->count; i++) { double z = Z1[(size_t)i * vn + r]; s += z * z; }
            if (s > lim_row) ok = 0;
        }
        for (int j = 0; j < ell && ok; j++)
            for (int pp = 0; pp < p->mu; pp++) {
                double cs = 0;
                for (int q2 = 0; q2 < TV_N; q2++) {
                    double z = Z2[(size_t)j * vn + (size_t)pp * TV_N + q2];
                    if (fabs(z) > lim_ent) ok = 0;
                    cs += z * z;
                }
                if (cs > lim_col) ok = 0;
            }
        if (!ok) break;
        challenge c1;
        first_challenge(&c1, bb1, c->k_auth, c->iota, e);
        second_challenge(C2, Z1, (size_t)b->count * vn, &c1, bb1, c->k_auth, c->iota, e, b->count, ell, ap->mode);
        for (int i = 0; i < b->count; i++) statement_col(&T[(size_t)i * mrows], view_com, &cv, p, &t, b->first + i, b->zero);
        /* W1 = M Z1 - c1 T */
        for (int i = 0; i < b->count; i++) {
            for (size_t r = 0; r < vn; r++) col[r] = Z1[(size_t)i * vn + r];
            M_mul_int(w, pub, mrows, col, scratch);
            for (int rr = 0; rr < mrows; rr++) {
                poly tmp;
                poly_mul_sparse(&tmp, &T[(size_t)i * mrows + rr], c1.pos, c1.sgn, TV_CH_W);
                poly_sub(&w[rr], &w[rr], &tmp);
            }
            g_absorb(&G, w, (size_t)mrows);
        }
        /* W2 = M Z2 - T C2 */
        for (int j = 0; j < ell; j++) {
            for (size_t r = 0; r < vn; r++) col[r] = Z2[(size_t)j * vn + r];
            M_mul_int(w, pub, mrows, col, scratch);
            memset(accq, 0, sizeof(int64_t) * (size_t)mrows * TV_N);
            for (int i = 0; i < b->count; i++) {
                int8_t cij = C2[(size_t)i * ell + j];
                if (!cij) continue;
                for (int rr = 0; rr < mrows; rr++) {
                    const uint64_t *tc = T[(size_t)i * mrows + rr].c;
                    int64_t *a = accq + (size_t)rr * TV_N;
                    if (cij > 0) for (int q2 = 0; q2 < TV_N; q2++) a[q2] += (int64_t)tc[q2];
                    else for (int q2 = 0; q2 < TV_N; q2++) a[q2] -= (int64_t)tc[q2];
                }
            }
            for (int rr = 0; rr < mrows; rr++) {
                poly tc;
                for (int q2 = 0; q2 < TV_N; q2++) tc.c[q2] = mod_from_int(accq[(size_t)rr * TV_N + q2]);
                poly_sub(&w[rr], &w[rr], &tc);
            }
            g_absorb(&G, w, (size_t)mrows);
        }
    }
    if (ok) {
        uint8_t h[TV_HASHBYTES];
        keccak_absorb(&G, c->salt, 32);
        keccak_finalize(&G);
        keccak_squeeze(&G, h, TV_HASHBYTES);
        ok = memcmp(h, c->h[c->iota], TV_HASHBYTES) == 0;
    }
    /* root opening */
    if (ok) {
        double n2 = 0;
        for (size_t i = 0; i < vn; i++) n2 += (double)c->rho_root[i] * (double)c->rho_root[i];
        if (sqrt(n2) > tv_beta(p)) ok = 0;
        poly *cr = malloc(sizeof(poly) * rows);
        tv_commit(cr, pub, c->V, c->rho_root);
        if (memcmp(cr, com + (size_t)(t.nnodes - 1) * rows, sizeof(poly) * rows) != 0) ok = 0;
        free(cr);
    }
    (void)L;
    free(T); free(w); free(scratch); free(C2); free(col); free(accq); free(com); free(blk);
    tree_free(&t);
    return ok;
}

/* ------------------------------------------------------------------ sizes and cleanup */

size_t agg_contrib_bytes(const agg_contrib *c, const tv_pub *pub, const agg_params *ap, int round) {
    const tv_params *p = &pub->prm;
    size_t vn = VN(p);
    bitw w;
    bw_init(&w, NULL, 0);
    if (round == 1) {
        for (size_t i = 0; i < (size_t)c->ninner * p->rows; i++) pack_poly(&w, &c->inner_com[i]);
        for (int a = 0; a < p->L; a++) pack_poly(&w, &c->V[a]);
        for (int i = 0; i < ap->kappa * TV_HASHBYTES; i++) bw_put(&w, 0, 8);
        return bw_finish(&w);
    }
    bw_put(&w, (uint64_t)c->iota, 8);
    for (int i = 0; i < 32; i++) bw_put(&w, 0, 8);
    int k1 = rice_param(ap->sigma1), k2 = rice_param(ap->sigma2);
    size_t nz1 = (size_t)(ap->nnodes + ap->ninner) * vn;
    for (size_t i = 0; i < nz1; i++) rice_put(&w, c->Z1[i], k1);
    for (size_t i = 0; i < (size_t)ap->E * vn * ap->ell; i++) rice_put(&w, c->Z2[i], k2);
    for (size_t i = 0; i < vn; i++) rice_put(&w, c->rho_root[i], 0);
    return bw_finish(&w);
}

void agg_state_free(agg_state *st) { tree_free(&st->tree); free(st->blocks); free(st->seed); free(st->salt); memset(st, 0, sizeof *st); }
void agg_contrib_free(agg_contrib *c) { free(c->inner_com); free(c->V); free(c->h); free(c->Z1); free(c->Z2); free(c->rho_root); memset(c, 0, sizeof *c); }
