/* Straightforward Keccak-f[1600] (FIPS 202). Lanes are little-endian. */
#include <string.h>
#include "keccak.h"

static const uint64_t RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808AULL, 0x8000000080008000ULL,
    0x000000000000808BULL, 0x0000000080000001ULL, 0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008AULL, 0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000AULL,
    0x000000008000808BULL, 0x800000000000008BULL, 0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL, 0x000000000000800AULL, 0x800000008000000AULL,
    0x8000000080008081ULL, 0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL};
static const unsigned ROTC[24] = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14,
                                  27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44};
static const unsigned PILN[24] = {10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4,
                                  15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1};

#define ROTL64(x, y) (((x) << (y)) | ((x) >> (64 - (y))))

static void keccakf(uint64_t st[25]) {
    uint64_t t, bc[5];
    for (int r = 0; r < 24; r++) {
        for (int i = 0; i < 5; i++) bc[i] = st[i] ^ st[i + 5] ^ st[i + 10] ^ st[i + 15] ^ st[i + 20];
        for (int i = 0; i < 5; i++) {
            t = bc[(i + 4) % 5] ^ ROTL64(bc[(i + 1) % 5], 1);
            for (int j = 0; j < 25; j += 5) st[j + i] ^= t;
        }
        t = st[1];
        for (int i = 0; i < 24; i++) {
            unsigned j = PILN[i];
            bc[0] = st[j];
            st[j] = ROTL64(t, ROTC[i]);
            t = bc[0];
        }
        for (int j = 0; j < 25; j += 5) {
            for (int i = 0; i < 5; i++) bc[i] = st[j + i];
            for (int i = 0; i < 5; i++) st[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
        }
        st[0] ^= RC[r];
    }
}

static void xor_byte(uint64_t *s, unsigned idx, uint8_t b) { s[idx >> 3] ^= (uint64_t)b << (8 * (idx & 7)); }
static uint8_t get_byte(const uint64_t *s, unsigned idx) { return (uint8_t)(s[idx >> 3] >> (8 * (idx & 7))); }

static void init(keccak_state *st, unsigned rate) {
    memset(st->s, 0, sizeof st->s);
    st->rate = rate;
    st->pos = 0;
    st->squeezing = 0;
}
void shake128_init(keccak_state *st) { init(st, SHAKE128_RATE); }
void shake256_init(keccak_state *st) { init(st, SHAKE256_RATE); }

void keccak_absorb(keccak_state *st, const uint8_t *in, size_t len) {
    while (len > 0) {
        if (st->pos == 0 && len >= st->rate) {          /* fast path: whole block */
            for (unsigned i = 0; i < st->rate / 8; i++) {
                uint64_t w = 0;
                for (int b = 7; b >= 0; b--) w = (w << 8) | in[8 * i + b];
                st->s[i] ^= w;
            }
            keccakf(st->s);
            in += st->rate;
            len -= st->rate;
            continue;
        }
        xor_byte(st->s, st->pos++, *in++);
        len--;
        if (st->pos == st->rate) {
            keccakf(st->s);
            st->pos = 0;
        }
    }
}

void keccak_finalize(keccak_state *st) {
    xor_byte(st->s, st->pos, 0x1F);
    xor_byte(st->s, st->rate - 1, 0x80);
    keccakf(st->s);
    st->pos = 0;
    st->squeezing = 1;
}

void keccak_squeeze(keccak_state *st, uint8_t *out, size_t len) {
    while (len > 0) {
        if (st->pos == st->rate) {
            keccakf(st->s);
            st->pos = 0;
        }
        *out++ = get_byte(st->s, st->pos++);
        len--;
    }
}

void shake128(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen) {
    keccak_state st;
    shake128_init(&st);
    keccak_absorb(&st, in, inlen);
    keccak_finalize(&st);
    keccak_squeeze(&st, out, outlen);
}
void shake256(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen) {
    keccak_state st;
    shake256_init(&st);
    keccak_absorb(&st, in, inlen);
    keccak_finalize(&st);
    keccak_squeeze(&st, out, outlen);
}
