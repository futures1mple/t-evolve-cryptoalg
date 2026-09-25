#include <math.h>
#include <stdlib.h>
#include "codec.h"

void bw_init(bitw *w, uint8_t *buf, size_t cap) { w->buf = buf; w->cap = cap; w->bytes = 0; w->acc = 0; w->nacc = 0; w->err = 0; }

static void bw_flush8(bitw *w) {
    while (w->nacc >= 8) {
        if (w->buf) {
            if (w->bytes < w->cap) w->buf[w->bytes] = (uint8_t)w->acc; else w->err = 1;
        }
        w->bytes++;
        w->acc >>= 8;
        w->nacc -= 8;
    }
}

void bw_put(bitw *w, uint64_t v, int nbits) {
    if (nbits <= 0) return;
    v &= (nbits == 64) ? ~0ULL : ((1ULL << nbits) - 1);
    w->acc |= v << w->nacc;
    w->nacc += nbits;
    bw_flush8(w);
}

size_t bw_finish(bitw *w) {
    if (w->nacc > 0) { w->nacc = 8; bw_flush8(w); w->nacc = 0; }
    return w->bytes;
}

void br_init(bitr *r, const uint8_t *buf, size_t len) { r->buf = buf; r->len = len; r->pos = 0; r->acc = 0; r->nacc = 0; r->err = 0; }

uint64_t br_get(bitr *r, int nbits) {
    while (r->nacc < nbits) {
        uint64_t b = 0;
        if (r->pos < r->len) b = r->buf[r->pos]; else r->err = 1;
        r->pos++;
        r->acc |= b << r->nacc;
        r->nacc += 8;
    }
    uint64_t v = r->acc & ((nbits == 64) ? ~0ULL : ((1ULL << nbits) - 1));
    r->acc >>= nbits;
    r->nacc -= nbits;
    return v;
}

void pack_poly(bitw *w, const poly *a) { for (int i = 0; i < TV_N; i++) bw_put(w, a->c[i], (int)RING.logq); }
void unpack_poly(bitr *r, poly *a) {
    for (int i = 0; i < TV_N; i++) {
        a->c[i] = br_get(r, (int)RING.logq);
        if (a->c[i] >= RING.q) r->err = 1;
    }
}

void rice_put(bitw *w, int64_t x, int k) {
    uint64_t ax = (uint64_t)(x < 0 ? -x : x);
    bw_put(w, x < 0, 1);
    uint64_t hi = ax >> k;
    while (hi >= 32) { bw_put(w, 0xFFFFFFFFULL, 32); hi -= 32; }
    bw_put(w, (1ULL << hi) - 1, (int)hi);   /* hi ones */
    bw_put(w, 0, 1);                         /* terminator */
    bw_put(w, ax, k);
}

int64_t rice_get(bitr *r, int k) {
    int neg = (int)br_get(r, 1);
    uint64_t hi = 0;
    while (br_get(r, 1)) {
        hi++;
        if (r->err || hi > (1u << 20)) { r->err = 1; return 0; }
    }
    uint64_t ax = (hi << k) | br_get(r, k);
    if (neg && ax == 0) r->err = 1;          /* non-canonical zero */
    return neg ? -(int64_t)ax : (int64_t)ax;
}

int rice_param(double sigma) { int k = (int)floor(log2(sigma)); return k < 0 ? 0 : k; }

void absorb_polys(keccak_state *st, const poly *a, size_t k) {
    size_t cap = (size_t)TV_N * RING.logq / 8 + 16;
    uint8_t *buf = malloc(cap);
    for (size_t j = 0; j < k; j++) {
        bitw w;
        bw_init(&w, buf, cap);
        pack_poly(&w, &a[j]);
        size_t len = bw_finish(&w);
        keccak_absorb(st, buf, len);
    }
    free(buf);
}
