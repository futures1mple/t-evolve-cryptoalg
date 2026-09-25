/* Ballots of T-EVOLVE with the EVOLVE-style ballot proof (Section 6 of the paper), in the
 * variant where the randomness of share k is derived from a seed s_k and only s_k is
 * encrypted to authority k:  r_k = SampleD_sigma(X(par, id, k, s_k)), where X is SHAKE256 with
 * domain byte 0x20 over (parameter digest || id || k || s_k). Binding id and k into X keeps the
 * leaves of one aggregation column independent even if a voter reuses a seed. */
#ifndef TV_BALLOT_H
#define TV_BALLOT_H

#include "bdlop.h"
#include "sample.h"

#define TV_CT_BYTES 32      /* payload of e_k: the seed (see pke_stub below) */

typedef struct {
    uint64_t id;
    poly *c;                /* (n+1) * rows */
    uint8_t *e;             /* n * TV_CT_BYTES */
    uint8_t chat[TV_HASHBYTES];
    int64_t *z;             /* (n+1) * mu * N */
    int64_t *or_r;          /* L * 2 * mu * N, index (a*2+b) */
    challenge *f0;          /* L */
} tv_ballot;

typedef struct {            /* voter-side secrets, kept for tests and experiments */
    poly *m;                /* (n+1) * L shares, m_0 = vote */
    int64_t *r;             /* (n+1) * mu * N */
    uint8_t *seeds;         /* n * 32 */
} tv_voter_secret;

typedef struct {            /* statistics of one proof generation */
    int attempts;
    double max_shift_ratio; /* max over attempts of ||shift|| / T */
    double resp_norm;       /* l2 norm of the accepted concatenation */
} tv_prove_stats;

/* authority keys; the PKE is a placeholder (see README): e_k = s_k xor SHAKE256(key_k || id) */
typedef struct { uint8_t key[32]; } tv_akey;

void tv_ballot_alloc(tv_ballot *b, const tv_params *p);
void tv_ballot_free(tv_ballot *b);
void tv_secret_alloc(tv_voter_secret *s, const tv_params *p);
void tv_secret_free(tv_voter_secret *s);

/* vote: bits v[0..L-1] in {0,1}; rnd_seed drives all voter randomness */
int tv_vote(tv_ballot *b, tv_voter_secret *sec, tv_prove_stats *st, const tv_pub *pub,
            const tv_akey *keys, uint64_t id, const int *v, const uint8_t rnd_seed[32]);
/* returns 1 iff the ballot proof verifies */
int tv_verify_ballot(const tv_pub *pub, const tv_ballot *b);

/* authority k: decrypt s_k, re-derive r_k, check it opens c_k (first d rows) with ||r_k|| <= beta;
 * outputs the share m_k (L polys) and r_k. Returns 1 if the opening is valid. */
int tv_check_share(const tv_pub *pub, const tv_akey *key, int k, const tv_ballot *b, poly *m_k, int64_t *r_k);
double tv_beta(const tv_params *p);

/* derive r = SampleD_sigma(X(par, id, k, seed)) */
void tv_rand_from_seed(int64_t *r, const tv_pub *pub, uint64_t id, int k, const uint8_t seed[32]);

/* serialization; sizes in bytes of the three parts */
typedef struct { size_t total, commitments, ciphertexts, proof; } tv_sizes;
size_t tv_ballot_encode(uint8_t *buf, size_t cap, const tv_pub *pub, const tv_ballot *b, tv_sizes *sz);
int tv_ballot_decode(tv_ballot *b, const tv_pub *pub, const uint8_t *buf, size_t len);
size_t tv_ballot_maxbytes(const tv_params *p);

#endif
