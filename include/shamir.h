/* Shamir sharing over R_q^L and the generalized Reed-Solomon parity checks. */
#ifndef TV_SHAMIR_H
#define TV_SHAMIR_H

#include "params.h"
#include "sample.h"

/* H in Z_q^{(n+1-t) x (n+1)}, H[j][i] = u_i alpha_i^j, alpha_0 = 0, alpha_k = k */
void shamir_parity(uint64_t *H, const tv_params *prm);
/* m[j*L + a] for j = 0..n: m_0 = v, m_k = P(k); coefficients a_l uniform from the stream */
void shamir_share(poly *m, const poly *v, const tv_params *prm, prg *g);
/* Lagrange coefficients at 0 for the index set T (size t, values in 1..n) */
void shamir_lagrange(uint64_t *lam, const int *T, int t);

#endif
