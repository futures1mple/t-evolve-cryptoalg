/* Combine and Verify of T-EVOLVE: interpolation of t column sums at 0. */
#ifndef TV_TALLY_H
#define TV_TALLY_H

#include "agg.h"

/* R = sum_{k in T} lambda_k V_k; returns 0 if every entry of R is an integer constant in
   [0, nacc] (the check of Verify), and writes the counts to counts[0..L-1] */
int tv_combine(long *counts, const tv_pub *pub, const int *T, const poly *const *V, long nacc);

#endif
