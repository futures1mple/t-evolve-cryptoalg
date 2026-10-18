#!/usr/bin/env python3
"""Hardness of the parameter sets of T-EVOLVE under four cost models of lattice reduction.

The parameter chain (tools/chain.py) targets 2^128 under the MATZOV cost model (classical). This
script re-estimates the chosen sets, the same sets under the quantum variant of the MATZOV model
(nearest-neighbour search costed in the depth x width metric of Albrecht, Gheorghiu, Postlethwaite
and Schanck, ASIACRYPT 2020), and under the core-SVP model of Alkim, Ducas, Poppelmann and Schwabe
(USENIX Security 2016), 2^(0.292 beta) classical and 2^(0.265 beta) quantum. For comparison it
estimates the parameters of Kyber-512 (ML-KEM-512) under the same models, the set for
N_V = 10^6 with module rank d = 8 at the same modulus, and the concrete instantiation of the encryption
with publicly verifiable opening (Appendix D of the paper): the ciphertext RLWE instance with ternary
secret and noise over N' + 256 samples (eo_2048_rlwe, eo_1024_rlwe) and NTRU key recovery for the
master key with f, g of standard deviation 1.17 sqrt(q'/2N') (eo_2048_ntru, estimator NTRU module,
circulant type; falls back to the LWE model if the NTRU module raises).

Requires SageMath and the lattice estimator (https://github.com/malb/lattice-estimator), commit
53da598, importable as `estimator` (or set ESTIMATOR_DIR). LWE: secret and error D_1, m = (d+L)N
samples, L = 1; the attacks arora-gb and bkw are skipped (never the minimum for these parameters).
SIS: Euclidean norm, m = (2d+L)N columns, beta_SIS of results/params_chain_signed.jsonl (for d = 8 from
build/param_report 1000000 4 3 1 -1 8 35184372088241 1 64).

usage: python3 tools/quantum_models.py SET MODEL > results/quantum_models/SET_MODEL.json
  SET   one of 1e4, 1e5, 1e6, 1e6_d8, kyber512, eo_2048_rlwe, eo_2048_ntru, eo_1024_rlwe
  MODEL one of MATZOV, MATZOV_q, ADPS16, ADPS16_q
Each run takes 1 to 20 minutes on one core.
"""
import json, os, sys, time

sys.path.insert(0, os.environ.get("ESTIMATOR_DIR", "lattice-estimator"))
from sage.all import log, oo
from estimator import LWE, NTRU, SIS, ND, RC
from estimator.reduction import ADPS16, MATZOV

N = 256
SETS = {  # d, q, beta_SIS (from results/params_chain_signed.jsonl; for d = 8 from param_report)
    "1e4": (7, 4398046510961, 3257500000000.0),
    "1e5": (7, 8796093022033, 8326850000000.0),
    "1e6": (7, 35184372088241, 24105000000000.0),
    "1e6_d8": (8, 35184372088241, 26660200000000.0),
}
EO = {"eo_2048_rlwe": (2048, 9123495937), "eo_2048_ntru": (2048, 9123495937), "eo_1024_rlwe": (1024, 2256363521)}
MODELS = {"MATZOV": RC.MATZOV, "MATZOV_q": MATZOV(nn="quantum"), "ADPS16": RC.ADPS16, "ADPS16_q": ADPS16(mode="quantum")}


def bits(res):
    return {k: (float(log(v["rop"], 2)) if v["rop"] != oo else None) for k, v in res.items()}


def betas(res):
    return {k: int(v["beta"]) for k, v in res.items() if "beta" in v and v["rop"] != oo}


def main():
    key, model = sys.argv[1], sys.argv[2]
    t0 = time.time()
    if key in EO:
        import math
        Np, qp = EO[key]
        if key.endswith("rlwe"):
            lwe_p = LWE.Parameters(n=Np, q=qp, Xs=ND.Uniform(-1, 1), Xe=ND.Uniform(-1, 1), m=Np + 256)
            lwe = LWE.estimate(lwe_p, red_cost_model=MODELS[model], deny_list=("arora-gb", "bkw"), quiet=True)
        else:
            sf = 1.17 * math.sqrt(qp / (2 * Np))
            ntru_p = NTRU.Parameters(n=Np, q=qp, Xs=ND.DiscreteGaussian(sf), Xe=ND.DiscreteGaussian(sf), ntru_type="circulant")
            try:
                lwe = NTRU.estimate(ntru_p, red_cost_model=MODELS[model], deny_list=("arora-gb", "bkw"), quiet=True)
            except Exception as ex:
                sys.stderr.write("NTRU.estimate failed (%r); LWE model of key recovery\n" % (ex,))
                lwe_p = LWE.Parameters(n=Np, q=qp, Xs=ND.DiscreteGaussian(sf), Xe=ND.DiscreteGaussian(sf), m=Np)
                lwe = LWE.estimate(lwe_p, red_cost_model=MODELS[model], deny_list=("arora-gb", "bkw"), quiet=True)
        b = bits(lwe)
        print(json.dumps({"set": key, "model": model, "estimator": "53da598", "N": Np, "q": qp, "attacks": b,
                          "min": min(v for v in b.values() if v), "seconds": round(time.time() - t0)}))
        return
    if key == "kyber512":
        lwe_p = LWE.Parameters(n=512, q=3329, Xs=ND.CenteredBinomial(3), Xe=ND.CenteredBinomial(3), m=512, tag="Kyber512")
        sis_p = None
    else:
        d, q, b = SETS[key]
        lwe_p = LWE.Parameters(n=d * N, q=q, Xs=ND.DiscreteGaussian(1.0), Xe=ND.DiscreteGaussian(1.0), m=(d + 1) * N)
        sis_p = SIS.Parameters(n=d * N, q=q, length_bound=b, m=(2 * d + 1) * N, norm=2) if b else None
    lwe = LWE.estimate(lwe_p, red_cost_model=MODELS[model], deny_list=("arora-gb", "bkw"), quiet=True)
    out = {"set": key, "model": model, "estimator": "53da598", "lwe": bits(lwe), "lwe_beta": betas(lwe),
           "hiding": min(v for v in bits(lwe).values() if v)}
    if sis_p is not None:
        sis = SIS.estimate(sis_p, red_cost_model=MODELS[model], quiet=True)
        out.update({"sis": bits(sis), "sis_beta": betas(sis), "binding": min(v for v in bits(sis).values() if v)})
    out["seconds"] = round(time.time() - t0)
    print(json.dumps(out))


if __name__ == "__main__":
    main()
