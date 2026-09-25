/* Authority aggregation of T-EVOLVE (Section 6 of the paper): a tree of recommitments with
 * fan-in l over the accepted share commitments, the amortized proofs Pi_open (short r with
 * A r = d_{u,s,1} for every node, leaves included) and Pi_zero (short r with C r = delta_{u,s}
 * for every inner node) of Baum and Lyubashevsky (2017), split into blocks of at most k
 * equations, with one joint rejection step per stage over all blocks, kappa attempts
 * committed in round 1 by h = G(W, salt), and the index of the first non-aborting attempt
 * announced before its responses are revealed.
 *
 * The second challenge C2 is either binary (as in BL17) or signed ({-1,1}); see
 * docs/PARAMETERS.md for the effect on the bounds. */
#ifndef TV_AGG_H
#define TV_AGG_H

#include "ballot.h"

typedef enum { C2_BINARY = 0, C2_SIGNED = 1 } c2_mode;

typedef struct {
    int fanin, k, ell, kappa;
    c2_mode mode;
    double alpha, logM;          /* sigma = alpha * T, M = exp(12/alpha + 1/(2 alpha^2)) */
    double T1, T2;               /* bounds on the norms of the two shifts (all blocks) */
    double sigma1, sigma2;
    double zinf_mult;            /* entries of Z2 must be <= zinf_mult * sigma2 (9) */
    double logQ;                 /* grinding: log2 of the adversary's hash queries */
    double tau;                  /* failure exponent (nats) of every tail bound */
    int E, E_open, E_zero;
    int nnodes, ninner;
} agg_params;

typedef struct {
    int nleaves, nlevels, nnodes;
    int *lvl_start, *lvl_count;
    int *child_start, *nchild;   /* children of node u (inner nodes) */
    poly *com;                   /* nnodes * rows */
    poly *msg;                   /* nnodes * L */
    int8_t *rnd;                 /* nnodes * mu*N */
} agg_tree;

typedef struct { int zero; int level; int first; int count; } agg_block;

/* what an authority publishes */
typedef struct {
    int k_auth;                  /* authority index 1..n */
    int ninner;
    poly *inner_com;             /* commitments of the inner nodes (level >= 1), in node order */
    poly *V;                     /* L polys */
    uint8_t (*h)[TV_HASHBYTES];  /* kappa commitments to first messages */
    /* round 2 */
    int iota;                    /* index of the non-aborting attempt, -1 if none */
    uint8_t salt[32];
    uint8_t leaf_digest[TV_HASHBYTES]; /* hash of the leaf commitments (the accepted ballots) */
    int64_t *Z1;                 /* sum over blocks of mu*N*count, column-major per block */
    int64_t *Z2;                 /* E * mu*N*ell */
    int64_t *rho_root;           /* mu*N */
} agg_contrib;

typedef struct {                 /* statistics of one aggregation */
    int attempts_used;           /* attempts examined in round 2 */
    double shift1_ratio, shift2_ratio;   /* max ||shift|| / T over attempts */
    double t_round1_ms, t_round2_ms;
    size_t bytes_round1, bytes_round2;
} agg_stats;

typedef struct {                 /* authority-side state between the rounds */
    agg_tree tree;
    agg_block *blocks;
    uint8_t (*seed)[32];         /* kappa attempt seeds */
    uint8_t (*salt)[32];
    prg rng;
} agg_state;

void agg_params_init(agg_params *ap, const tv_params *p, long nleaves, c2_mode mode, double logQ);
void agg_params_print(const agg_params *ap, const tv_params *p);

/* round 1: build the tree over the accepted share commitments of authority k (leaves in board
   order, with messages and randomness from tv_check_share), prepare kappa attempts */
int agg_round1(agg_state *st, agg_contrib *out, agg_stats *stats, const tv_pub *pub, const agg_params *ap, int k,
               const poly *leaf_com, const poly *leaf_msg, const int64_t *leaf_rnd, int nleaves, const uint8_t seed[32]);
/* digest of BB_1: everything posted in round 1 */
void agg_bb1_digest(uint8_t out[TV_HASHBYTES], const tv_pub *pub, const agg_contrib *c, int nc);
/* round 2: find the first non-aborting attempt; fills iota (or -1) and, if found, the responses */
int agg_round2(agg_state *st, agg_contrib *out, agg_stats *stats, const tv_pub *pub, const agg_params *ap,
               const uint8_t bb1[TV_HASHBYTES]);
int agg_verify(const tv_pub *pub, const agg_params *ap, const agg_contrib *c, const poly *leaf_com, int nleaves,
               const uint8_t bb1[TV_HASHBYTES]);

void agg_state_free(agg_state *st);
void agg_contrib_free(agg_contrib *c);
size_t agg_contrib_bytes(const agg_contrib *c, const tv_pub *pub, const agg_params *ap, int round);

/* structure helpers (also used by the experiments) */
int agg_tree_shape(agg_tree *t, int nleaves, int fanin);   /* fills levels and children only */
int agg_make_blocks(agg_block **blocks, const agg_tree *t, int k);

#endif
