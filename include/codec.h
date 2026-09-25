/* Bit-level encoding used for sizes and for hashing:
 *   - elements of Z_q packed in ceil(log2 q) bits,
 *   - Gaussian integers Golomb-Rice coded (sign, unary high part, k low bits). */
#ifndef TV_CODEC_H
#define TV_CODEC_H

#include <stddef.h>
#include <stdint.h>
#include "ring.h"

typedef struct {
    uint8_t *buf;
    size_t cap, bytes;   /* capacity, bytes in use */
    uint64_t acc;
    int nacc;            /* bits in acc */
    int err;
} bitw;

typedef struct {
    const uint8_t *buf;
    size_t len, pos;
    uint64_t acc;
    int nacc;
    int err;
} bitr;

void bw_init(bitw *w, uint8_t *buf, size_t cap);
void bw_put(bitw *w, uint64_t v, int nbits);     /* nbits <= 57 */
size_t bw_finish(bitw *w);                        /* returns number of bytes */
void br_init(bitr *r, const uint8_t *buf, size_t len);
uint64_t br_get(bitr *r, int nbits);

void pack_poly(bitw *w, const poly *a);           /* logq bits per coefficient */
void unpack_poly(bitr *r, poly *a);               /* sets r->err if a value >= q */
void rice_put(bitw *w, int64_t x, int k);
int64_t rice_get(bitr *r, int k);
int rice_param(double sigma);                     /* floor(log2 sigma) */

/* hashing of polynomials: absorb their packed encoding */
#include "keccak.h"
void absorb_polys(keccak_state *st, const poly *a, size_t k);

#endif
