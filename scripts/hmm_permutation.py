#!/usr/bin/env python3
"""Permutation test for the gene-order HMM of block_lore.py.

Statistic: lnL(HMM) - lnL(independent mixture), both fitted on the same
gene-ordered families (each with its own class weights). Under the null,
timing classes are not autocorrelated along gene order, so families are
exchangeable among ordered positions. Each permutation reassigns the ordered
families to the positions at random (segment number and lengths kept), refits
the HMM, and recomputes the statistic. The independent mixture does not depend
on order, so only the HMM is refitted. p = (1 + #{perm >= obs}) / (N + 1).

usage: hmm_permutation.py --grid r=dir ... --order order.tsv --perm 100
                          [--perm-start 0] [--seed 1] --out prefix
       hmm_permutation.py --merge prefix_shard*_hmm_perm.json --out prefix
Permutations run serially; split them into shards (--perm-start / --perm,
separate processes) and combine with --merge. Permutation k always uses seed
seed*100003 + k, so the merged result does not depend on the sharding.
Writes prefix_hmm_perm.json (observed fit, null statistics, p).
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import block_lore as bl

G = {}


def _init(LL, segs, pi0):
    G.update(LL=LL, segs=segs, pi0=pi0)


def _perm(seed):
    rng = np.random.default_rng(seed)
    flat = np.concatenate(G["segs"])
    perm = rng.permutation(flat)
    cuts = np.cumsum([len(s) for s in G["segs"]])[:-1]
    lh, _, rho = bl.fit_hmm(G["LL"], np.split(perm, cuts), G["pi0"])
    return lh, rho


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", nargs="+")
    ap.add_argument("--order")
    ap.add_argument("--perm", type=int, default=100)
    ap.add_argument("--perm-start", type=int, default=0)
    ap.add_argument("--merge", nargs="+")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.merge:
        sh = [json.load(open(f)) for f in a.merge]
        out = {k: v for k, v in sh[0].items() if k not in ("null_gains", "perm_start")}
        null = np.array([g for x in sh for g in x["null_gains"]])
        obs = out["observed_gain"]
        out.update(perm=int(len(null)), null_gain=dict(mean=float(null.mean()), sd=float(null.std()),
                   max=float(null.max()), q95=float(np.quantile(null, 0.95))),
                   p=float((1 + np.sum(null >= obs)) / (len(null) + 1)), null_gains=null.tolist())
        json.dump(out, open(a.out + "_hmm_perm.json", "w"), indent=1)
        print(json.dumps({k: v for k, v in out.items() if k != "null_gains"}, indent=1))
        return
    df, r = bl.load_grid(a.grid)
    fams = list(df.index); pos = {f: i for i, f in enumerate(fams)}
    o = pd.read_csv(a.order, sep="\t", header=None, names=["family", "chrom", "position"],
                    dtype={"family": str, "chrom": str}, keep_default_na=False)
    o = o[o.family.isin(fams)].sort_values(["chrom", "position"])
    segs = [np.array([pos[f] for f in g.family]) for _, g in o.groupby("chrom", sort=False)]
    segs = [s for s in segs if len(s)]
    LL = df.values
    flat = np.concatenate(segs)
    li, pii, _ = bl.em_mixture(LL[flat])
    pi0 = bl.em_mixture(LL)[1]
    lh, pih, rho = bl.fit_hmm(LL, segs, pi0)
    obs = lh - li
    _init(LL, segs, pi0)
    res = [_perm(a.seed * 100003 + k) for k in range(a.perm_start, a.perm_start + a.perm)]
    null = np.array([x[0] - li for x in res])
    out = dict(n_ordered=int(len(flat)), n_segments=len(segs), lnL_indep=li, lnL_hmm=lh, rho=rho,
               mean_run_length=1.0 / rho, observed_gain=obs, perm=a.perm,
               null_gain=dict(mean=float(null.mean()), sd=float(null.std()), max=float(null.max()),
                              q95=float(np.quantile(null, 0.95))),
               null_rho_median=float(np.median([x[1] for x in res])),
               p=float((1 + np.sum(null >= obs)) / (a.perm + 1)), perm_start=a.perm_start,
               null_gains=null.tolist())
    json.dump(out, open(a.out + "_hmm_perm.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "null_gains"}, indent=1))


if __name__ == "__main__":
    main()
