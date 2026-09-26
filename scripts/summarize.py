#!/usr/bin/env python3
"""Summarize one results folder (bench.csv, ballot_rej_summary.csv) as SUMMARY.md.

Measured values are reported as measured. For N_V = 10^5 and 10^6 the aggregation cost is
extrapolated from the measured per-block costs at the largest measured N_V, using the number of
blocks E of each tree; these rows are marked "extrapolated".

usage: python3 scripts/summarize.py results/<machine>/<date>
"""
import csv, glob, math, os, sys
from math import comb

FANIN, K = 30, 500


def blocks(nv):
    counts = [nv]
    while counts[-1] > 1 or len(counts) == 1:
        counts.append((counts[-1] + FANIN - 1) // FANIN)
    nodes = sum(counts)
    return (nodes + K - 1) // K + sum((c + K - 1) // K for c in counts[1:])


def fisher(ps):
    """Fisher's method: X = -2 sum ln p ~ chi^2 with 2k degrees of freedom (closed form for even dof)"""
    x = -2 * sum(math.log(max(p, 1e-300)) for p in ps) / 2
    k = len(ps)
    term, tot = 1.0, 1.0
    for i in range(1, k):
        term *= x / i
        tot += term
    return math.exp(-x) * tot


def rounds_model(p, kappa, n, t):
    """expected rounds and attempts when every authority prepares kappa attempts per round"""
    s = 1 - (1 - p) ** kappa
    ok_all = sum(comb(n, i) * s ** i * (1 - s) ** (n - i) for i in range(t, n + 1))
    ok_worst = s ** t
    return s, 1 / ok_all, 1 / ok_worst


def main(folder):
    rows = list(csv.DictReader(open(os.path.join(folder, "bench.csv"))))
    get = lambda g, name, cfgpart="": [r for r in rows if r["group"] == g and r["name"].strip() == name and cfgpart in r["config"]]
    out = []
    cpu = rows[0]["cpu"] if rows else "?"
    out.append(f"# Measurements\n\nMachine: {cpu}; compiler {rows[0]['compiler']}; flags `{rows[0]['flags']}`; commit {rows[0]['commit']}.\n")
    out.append("All times are medians over the repetitions (single thread).\n")
    out.append("## Ballots\n\n| Scheme | Vote (ms) | Verify (ms) | Check_k (ms) | Attempts | Size (KiB) |\n|---|---|---|---|---|---|")
    configs = sorted({r["config"] for r in rows if r["group"] == "ballot"})
    for c in configs:
        def v(name):
            x = [r for r in rows if r["group"] == "ballot" and r["name"].strip() == name and r["config"] == c]
            return float(x[0]["median"]) if x else float("nan")
        ver = v("decode+verify") if "T-EVOLVE" in c else v("verify")
        out.append(f"| {c} | {v('vote (incl. all attempts)'):.1f} | {ver:.1f} | {v('check_share (one authority)'):.2f} | "
                   f"{v('attempts (mean)'):.2f} | {v('ballot size'):.1f} |")
    agg = [r for r in rows if r["group"] == "agg"]
    if agg:
        out.append("\n## Aggregation (per authority)\n")
        out.append("| N_V | E | Round 1 (s) | Round 2 (s) | VerAgg (s) | Published per ballot (KiB) | Status |\n|---|---|---|---|---|---|---|")
        nvs = sorted({int(r["config"].split()[0].split("=")[1]) for r in agg})
        last = None
        for nv in nvs:
            f = lambda name: float([r for r in agg if r["name"].strip() == name and r["config"].startswith(f"NV={nv} ")][0]["median"])
            E = f("blocks E")
            last = dict(nv=nv, E=E, r1b=f("round 1 per block per attempt"), r2b=f("round 2 per block per attempt"),
                        vb=f("VerAgg per block"), kib=f("published per ballot per authority"))
            out.append(f"| {nv} | {E:.0f} | {f('round 1 per authority per round'):.1f} | {f('round 2 per authority per round'):.1f} | "
                       f"{f('VerAgg per authority'):.1f} | {last['kib']:.1f} | measured |")
        for nv in (10**5, 10**6):
            if last and nv > last["nv"]:
                E = blocks(nv)
                out.append(f"| {nv} | {E} | {4 * E * last['r1b']:.0f} | {3 * E * last['r2b']:.0f} | {E * last['vb']:.0f} | — | extrapolated |")
        s, r_all, r_worst = rounds_model(math.exp(-2 * 86929 / 157922), 4, 4, 3)
        out.append(f"\nRounds (kappa = 4, attempt success 1/M^2 = 0.333): an authority succeeds in a round with probability {s:.2f}; "
                   f"expected rounds {r_all:.2f} with 4 honest authorities, at most {r_worst:.2f} with 3. If fewer than t authorities "
                   f"announce an attempt, no responses are published. Misbehaviour cannot make a round fail that would otherwise succeed; "
                   f"an authority that announces and then posts no valid contribution only makes the honest authorities publish their "
                   f"responses in a round that fails anyway, and is excluded afterwards (at most n - t such rounds). "
                   f"The benchmark has no faulty authorities.\n")
    rej = os.path.join(folder, "ballot_rej_summary.csv")
    if os.path.exists(rej):
        out.append("## Ballot rejection experiment\n\n| n | L | log2 sigma_J | Ballots | Attempts (mean ± 95%) | max ‖s‖/T | Verified | GOF p | Independence p |\n|---|---|---|---|---|---|---|---|---|")
        for r in csv.DictReader(open(rej)):
            out.append(f"| {r['n']} | {r['L']} | {float(r['log2_sigma_J']):.2f} | {r['ballots']} | {float(r['attempts_mean']):.2f} ± {float(r['attempts_ci95']):.2f} | "
                       f"{float(r['max_shift_over_T']):.2f} | {r['verified']} | {float(r['gof_p']):.3f} | {float(r['indep_p']):.3f} |")
        out.append("\nExpected attempts: M = exp(43669/39762) = 3.00 for every L.")
        reps = sorted(glob.glob(os.path.join(folder, "replications", "*", "ballot_rej_summary.csv")))
        g, ind = [], []
        for f in [rej] + reps:
            for r in csv.DictReader(open(f)):
                g.append(float(r["gof_p"])); ind.append(float(r["indep_p"]))
        out.append(f"\nOver {len(g)} tests of each kind ({1 + len(reps)} runs with different seeds), Fisher's combined p-value is "
                   f"{fisher(g):.2f} for goodness of fit and {fisher(ind):.2f} for witness independence; the smallest single p-values "
                   f"are {min(g):.3f} and {min(ind):.3f} (Bonferroni threshold 0.05/{len(g)} = {0.05/len(g):.4f}).")
    open(os.path.join(folder, "SUMMARY.md"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "results/linux_cloud")
