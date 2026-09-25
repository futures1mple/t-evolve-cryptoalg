/* EVOLVE ballots (del Pino, Lyubashevsky, Neven, Seiler, CCS 2017) on the same core, for a
 * like-for-like comparison with T-EVOLVE: additive shares of a yes/no vote among n
 * authorities, one commitment per share, the seeds of the share randomness encrypted to the
 * authorities (as in the T-EVOLVE variant), and EVOLVE's OR-proof on the sum of the
 * commitments. The rejection parameter of the OR-proof is set by the same rigorous rule as in
 * T-EVOLVE (sigma = 11 T, T a tail bound on ||f r||). The authorities' aggregation is the same
 * tree and the same amortized proofs as in T-EVOLVE (agg.h), run by all n authorities. */
#ifndef TV_EVOLVE_H
#define TV_EVOLVE_H

#include "ballot.h"

typedef struct {
    double T, sigma_OR, B_OR, logM;
} ev_params;

typedef struct {
    uint64_t id;
    poly *c;                /* n * rows */
    uint8_t *e;             /* n * TV_CT_BYTES */
    uint8_t chat[TV_HASHBYTES];
    int64_t *r0, *r1;       /* mu*N each */
    challenge f0;
} ev_ballot;

void ev_params_init(ev_params *ep, const tv_params *p);
void ev_ballot_alloc(ev_ballot *b, const tv_params *p);
void ev_ballot_free(ev_ballot *b);
/* shares m_1..m_n (L=1) and randomness returned for the authorities' side */
int ev_vote(ev_ballot *b, poly *shares, int64_t *rnd, int *attempts, const tv_pub *pub, const ev_params *ep,
            const tv_akey *keys, uint64_t id, int vote, const uint8_t seed[32]);
int ev_verify(const tv_pub *pub, const ev_params *ep, const ev_ballot *b);
size_t ev_ballot_bytes(const tv_pub *pub, const ev_params *ep, const ev_ballot *b, tv_sizes *sz);

#endif
