# Hiding (M-LWE) estimate for L = 2 at the moduli of Table 2; usage: python tools/hiding_L2.py <path-to-lattice-estimator> <q>  (SageMath, estimator commit 53da598)
import sys, json, time
sys.path.insert(0, sys.argv[1])
from sage.all import log, oo
from estimator import LWE, ND, RC
N=256; d=7; L=2; q=int(sys.argv[2])
t0=time.time()
lwe = LWE.estimate(LWE.Parameters(n=d*N, q=q, Xs=ND.DiscreteGaussian(1.0), Xe=ND.DiscreteGaussian(1.0), m=(d+L)*N), red_cost_model=RC.MATZOV, quiet=True)
lb={k:(float(log(v["rop"],2)) if v["rop"]!=oo else None) for k,v in lwe.items()}
print(json.dumps({"d":d,"L":L,"q":q,"hiding_bits":min(v for v in lb.values() if v),"lwe":lb,"seconds":time.time()-t0}))
