/* Keccak-f[1600] and the SHAKE128/SHAKE256 extendable-output functions (FIPS 202). */
#ifndef TV_KECCAK_H
#define TV_KECCAK_H

#include <stddef.h>
#include <stdint.h>

#define SHAKE128_RATE 168
#define SHAKE256_RATE 136

typedef struct {
    uint64_t s[25];
    unsigned rate;   /* in bytes */
    unsigned pos;    /* absorb: bytes absorbed in current block; squeeze: bytes consumed */
    int squeezing;
} keccak_state;

void shake128_init(keccak_state *st);
void shake256_init(keccak_state *st);
void keccak_absorb(keccak_state *st, const uint8_t *in, size_t len);
void keccak_finalize(keccak_state *st);            /* pads with the SHAKE domain byte 0x1F */
void keccak_squeeze(keccak_state *st, uint8_t *out, size_t len);

/* one-shot helpers */
void shake128(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);
void shake256(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);

#endif
