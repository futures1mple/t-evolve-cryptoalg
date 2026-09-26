#!/usr/bin/env python3
"""Compare two parameter-chain outputs (results/params_chain_signed.jsonl): for every numeric bound
present in both, print the largest change in log2 over all records, and the lines that differ in
(N_V, L, d, log q, hiding, binding)."""
import json, math, sys

def load(p):
    return [json.loads(l) for l in open(p) if l.strip()]

def flat(o, p=''):
    r = {}
    if isinstance(o, dict):
        for k, v in o.items():
            r.update(flat(v, p + '.' + k))
    elif isinstance(o, (int, float)) and not isinstance(o, bool):
        r[p] = o
    return r

a, b = load(sys.argv[1]), load(sys.argv[2])
mx = {}
for x, y in zip(a, b):
    fx, fy = flat(x), flat(y)
    for k in fx:
        if 'second' in k or 'time' in k:
            continue            # running time of the estimator, not a bound
        if k in fy and fx[k] > 0 and fy[k] > 0:
            d = abs(math.log2(fy[k]) - math.log2(fx[k]))
            mx[k] = max(mx.get(k, 0.0), d)
for k, v in sorted(mx.items(), key=lambda t: -t[1]):
    if v > 0:
        print(f"{k:40s} {v:.4f} bit")
