#!/usr/bin/env python3
"""Opening bound of EVOLVE (del Pino, Lyubashevsky, Neven, Seiler, CCS 2017) versus its modulus.

Part 1 evaluates the parameter constraints listed in Section 5 of the EVOLVE paper for the
parameters of its Table 2 (n = 256, q = 2^31 - 2^7 - 2^5 + 1, d = 7, sigma = 1, N_A = 4, l = 30)
and prints the resulting lower bound on the opening bound B_r.

Part 2 constructs, for a random commitment key of the same shape, a vector z != 0 mod q with
A z = 0 and B z != 0 by linear algebra alone, splits it as z = r - r', and checks that r and r'
are both below that lower bound: a commitment C r + (0, m) is then also opened by r' to the
message m + B z != m, so both openings are valid under the opening bound the constraints require.

Appendix A of the T-EVOLVE paper reports the output of this script.
usage: python3 scripts/evolve_opening_check.py [d] [seed]      (needs numpy)
"""
import math
import sys

import numpy as np

N = 256
q = 2**31 - 2**7 - 2**5 + 1
d = int(sys.argv[1]) if len(sys.argv) > 1 else 7
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
sigma, NA, l = 1.0, 4, 30
m = 2 * d + 1            # width of the commitment randomness (one message row)
lg = math.log2

# ---------------------------------------------------------------- part 1: EVOLVE, Section 5
w = math.sqrt(N * m)                                   # sqrt(n(2d+1)) in EVOLVE's notation
B_OR = 2 * NA * w * sigma                              # correctness of pi^V
B_OR_p = 44 * math.sqrt(60 * N * m) * B_OR             # zero knowledge of pi^V
B_A1 = 2 * w * sigma                                   # correctness of pi^A
B_A1_p = 2684 * N * math.sqrt(m) * B_A1                # zero knowledge of pi^A
B_A2 = 2 * (l + 1) * w * sigma                         # correctness of pi^A'
B_A2_p = 2684 * N * math.sqrt(m) * B_A2                # zero knowledge of pi^A'
cons = {
    "(3)  2 B'_OR": 2 * B_OR_p,
    "(4)  2 sqrt(60) N_A B'_Amo,1": 2 * math.sqrt(60) * NA * B_A1_p,
    "(7)  (l+1) B'_Amo,1": (l + 1) * B_A1_p,
    "(8)  B'_Amo,2": B_A2_p,
}
print(f"EVOLVE parameters: n={N}, q=2^31-2^7-2^5+1 (2^{lg(q):.2f}), d={d}, sigma={sigma}, N_A={NA}, l={l}")
print(f"  B_OR = 2^{lg(B_OR):.2f}, B'_OR = 2^{lg(B_OR_p):.2f}")
print(f"  B_Amo,1 = 2^{lg(B_A1):.2f}, B'_Amo,1 = 2^{lg(B_A1_p):.2f}")
print(f"  B_Amo,2 = 2^{lg(B_A2):.2f}, B'_Amo,2 = 2^{lg(B_A2_p):.2f}")
for k, v in cons.items():
    print(f"  consistency {k:30s} <= B_r  requires B_r >= 2^{lg(v):.2f}")
Br = max(cons.values())
print(f"  hence B_r >= 2^{lg(Br):.2f}, and the binding bound 2 B_r = 2^{lg(2 * Br):.2f} exceeds q = 2^{lg(q):.2f}")

# ---------------------------------------------------------------- part 2: two openings
rng = np.random.default_rng(seed)


def rot(a):
    """Negacyclic N x N matrix of the ring element a."""
    M = np.empty((N, N), dtype=np.int64)
    col = a.copy()
    for j in range(N):
        M[:, j] = col
        col = np.roll(col, 1)
        col[0] = (-col[0]) % q
    return M


Ar = rng.integers(0, q, size=(d, m, N), dtype=np.int64)
Br_ring = rng.integers(0, q, size=(1, m, N), dtype=np.int64)
A = np.block([[rot(Ar[i, j]) for j in range(m)] for i in range(d)])
B = np.block([[rot(Br_ring[0, j]) for j in range(m)]])
nz = d * N
A1, A2 = A[:, :nz].copy(), A[:, nz:]
z2 = np.zeros(A2.shape[1], dtype=np.int64)
z2[0] = 1
rhs = (-(A2 @ z2)) % q
M = np.concatenate([A1, rhs[:, None]], axis=1) % q
for c in range(nz):                                    # Gauss-Jordan elimination mod q
    p = c + np.nonzero(M[c:, c])[0][0]
    if p != c:
        M[[c, p]] = M[[p, c]]
    M[c] = (M[c] * pow(int(M[c, c]), q - 2, q)) % q
    col = M[:, c].copy()
    col[c] = 0
    rows = np.nonzero(col)[0]
    M[rows] = (M[rows] - (col[rows, None] * M[c][None, :]) % q) % q
z = np.concatenate([M[:, -1], z2])


def mulmod(X, v):
    lo, hi = v & 0xFFFF, v >> 16
    return ((X @ lo) % q + (((X @ hi) % q) * 65536) % q) % q


assert not np.any(mulmod(A, z)), "A z != 0"
Bz = mulmod(B, z)
zc = np.where(z > q // 2, z - q, z).astype(np.float64)
r = np.round(zc / 2)
rp = r - zc
print(f"kernel vector: A z = 0 mod q, B z != 0: {bool(np.any(Bz))}")
print(f"  ||z||_2 = 2^{lg(np.linalg.norm(zc)):.2f}  (expected about q sqrt(dN/12) = 2^{lg(q * math.sqrt(nz / 12)):.2f})")
print(f"  z = r - r': ||r||_2 = 2^{lg(np.linalg.norm(r)):.2f}, ||r'||_2 = 2^{lg(np.linalg.norm(rp)):.2f}")
ok = max(np.linalg.norm(r), np.linalg.norm(rp)) <= Br
print(f"  both below the required B_r = 2^{lg(Br):.2f}: {ok}  ->  two valid openings of one commitment: {ok and bool(np.any(Bz))}")
