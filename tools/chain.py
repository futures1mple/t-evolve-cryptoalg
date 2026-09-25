#!/usr/bin/env python3
"""Parameter chain of T-EVOLVE.

For every number of voters N_V and every ballot format, takes the bounds computed by
build/param_report (the same C code as the implementation), searches the smallest module rank d
and, for it, the smallest bit length logq and the largest prime q = 17 (mod 32) below 2^logq
with q > beta_SIS such that the lattice estimator
reports at least 2^128 for both M-LWE (hiding) and M-SIS (binding), and writes every estimator
result (all attacks, not only the minimum) to a JSON-lines log.

Requires SageMath and the lattice estimator (https://github.com/malb/lattice-estimator),
commit 53da598, importable as `estimator`.

usage: python3 tools/chain.py [out.jsonl] [mode: signed|binary] [logQ]
"""
import json, math, subprocess, sys, time

from sage.all import log, oo
from estimator import LWE, SIS, ND, RC

N = 256
OUT = sys.argv[1] if len(sys.argv) > 1 else "results/params_chain.jsonl"
MODE = sys.argv[2] if len(sys.argv) > 2 else "signed"
LOGQ = sys.argv[3] if len(sys.argv) > 3 else "64"
TARGET = 128


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_above(logq):
    q = 2 ** logq + 17
    while not is_prime(q):
        q += 32
    return q


def prime_below(logq):
    """Largest prime q = 17 (mod 32) below 2^logq, so that elements of Z_q take exactly logq bits."""
    q = 2 ** logq - 15
    while not is_prime(q):
        q -= 32
    return q


def report(NV, n, t, L, w, d, q):
    mode = "1" if MODE == "signed" else "0"
    out = subprocess.run(["./build/param_report", str(NV), str(n), str(t), str(L), str(w), str(d), str(q), mode, LOGQ],
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out)


def rop_bits(res):
    return {k: (float(log(v["rop"], 2)) if v["rop"] != oo else None) for k, v in res.items()}


def estimate(d, L, q, beta):
    t0 = time.time()
    lwe = LWE.estimate(LWE.Parameters(n=d * N, q=q, Xs=ND.DiscreteGaussian(1.0), Xe=ND.DiscreteGaussian(1.0), m=(d + L) * N),
                       red_cost_model=RC.MATZOV, quiet=True)
    if beta >= q:
        sis = None
        sis_bits = {"trivial (beta >= q)": 0.0}
    else:
        sis = SIS.estimate(SIS.Parameters(n=d * N, q=q, length_bound=beta, m=(2 * d + L) * N, norm=2), quiet=True)
        sis_bits = rop_bits(sis)
    lwe_bits = rop_bits(lwe)
    h = min(v for v in lwe_bits.values() if v is not None)
    b = min(v for v in sis_bits.values() if v is not None)
    return h, b, lwe_bits, sis_bits, time.time() - t0


def main():
    configs = [(NV, 4, 3, L, w) for NV in (10**4, 10**5, 10**6) for (L, w) in ((1, -1), (2, 1))]
    with open(OUT, "a") as f:
        f.write(json.dumps({"# chain": "T-EVOLVE", "mode": MODE, "logQ": LOGQ, "estimator": "53da598",
                            "lwe_cost_model": "MATZOV", "sis_norm": 2, "target": TARGET, "time": time.ctime()}) + "\n")
        for (NV, n, t, L, w) in configs:
            chosen = None
            for d in (7, 8, 9, 10):
                r0 = report(NV, n, t, L, w, d, prime_above(40))
                lq0 = math.floor(r0["log2beta"]) + 1
                for logq in range(lq0, lq0 + 4):
                    if logq > 50:
                        break
                    q = prime_below(logq)
                    r = report(NV, n, t, L, w, d, q)
                    if r["beta_SIS"] >= q:
                        continue
                    h, b, lb, sb, dt = estimate(d, L, q, r["beta_SIS"])
                    rec = dict(r, hiding_bits=h, binding_bits=b, lwe_attacks=lb, sis_attacks=sb, seconds=round(dt, 1))
                    f.write(json.dumps(rec) + "\n")
                    f.flush()
                    print(f"NV={NV} L={L} d={d} logq={logq} beta=2^{r['log2beta']:.1f} hiding=2^{h:.1f} binding=2^{b:.1f}", flush=True)
                    if h >= TARGET and b >= TARGET:
                        chosen = rec
                        break
                    if h < TARGET:
                        break
                if chosen:
                    break
            if chosen:
                f.write(json.dumps({"chosen": True, "NV": NV, "L": L, "d": chosen["d"], "q": chosen["q"],
                                    "hiding_bits": chosen["hiding_bits"], "binding_bits": chosen["binding_bits"]}) + "\n")
                f.flush()


if __name__ == "__main__":
    main()
