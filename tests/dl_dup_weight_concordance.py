#!/usr/bin/env python3
"""Regression gate for the UndatedDL duplication weighting (forward-process
weights: factor 2 on every D and WGD/LORe doubling split; upstream AleRax uses
the same factor for D since 8e34d87).

Same design as dtl_dup_weight_concordance.py, but for UndatedDL. Gene families
are forward-simulated (scripts/kalerax_sim.py) under two parameter sets, each
family with 3-6 genes is collapsed to its unrooted species-labelled shape, and
one representative per frequent shape is scored with kalerax at FIXED
parameters. For an exact likelihood (up to a parameter-independent per-tree
constant) X = freq(shape)/L(shape) must be the same under both sets.

Two pairs are tested:
  rates: D 0.15 vs 0.90 (L 0.60)  -> tests the D-split weight;
  wgd:   WGD on (A,B) with q 0.2 vs 0.95 (D 0.15, L 0.60) -> tests the WGD
         doubling-split weight (dupBracket).
The removed ALE weighting (each split once; kalerax <= 0a024d6) is discordant
because the number of duplication splits differs between reconciliations of
one gene tree (e.g. S versus D followed by two SL). At 8M families it was
rejected in both pairs (rates chi2 332.2/211, p 2e-7; wgd 600.6/281, p 3e-25)
while the forward weights were not (185.4/211, p 0.90; 270.1/281, p 0.67).

usage: dl_dup_weight_concordance.py <kalerax> <workdir> [N]
exit 0 = PASS (forward weights not rejected, chi2 p > 1e-3, in either pair)
Default N = 8,000,000 per parameter set (~10 min); smaller N lacks power.
"""
import os, subprocess, sys
from collections import Counter
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
from kalerax_sim import SpeciesTree, Params, simulate_families
from dtl_dup_weight_concordance import shape_key

SP = "(O:1,((A:1,B:1):1,(C:1,D:1):1):1);"
PAIRS = {
    "rates": {"a": dict(D=0.15, L=0.60, q=None), "b": dict(D=0.90, L=0.60, q=None)},
    "wgd": {"a": dict(D=0.15, L=0.60, q=0.2), "b": dict(D=0.15, L=0.60, q=0.95)},
}


def run_pair(kalerax, wd, N, name, pair, tree, sp):
    e_ab = [i for i in range(tree.n) if tree.leafset[i] == frozenset("AB")][0]
    X = {}
    for ci, (tag, c) in enumerate(pair.items()):
        wgds = [(e_ab, c["q"], 1.0)] if c["q"] is not None else []
        fams, _ = simulate_families(tree, Params(tree, c["D"], c["L"], 0.0, wgds=wgds), N,
                                    211 + 7 * ci + (13 if name == "wgd" else 0), keep_all=True)
        freq, rep = Counter(), {}
        for f in fams:
            if 3 <= f["n_genes"] <= 6:
                k = shape_key(f["newick"]); freq[k] += 1; rep.setdefault(k, f)
        top = [k for k, n in freq.most_common() if n >= 150]
        lines = ["[FAMILIES]"]
        for i, k in enumerate(top):
            f = rep[k]
            t = os.path.join(wd, f"{name}_{tag}_{i}.nw"); m = os.path.join(wd, f"{name}_{tag}_{i}.link")
            open(t, "w").write(f["newick"] + "\n")
            open(m, "w").write("".join(f"{g}\t{s}\n" for g, s in f["mapping"]))
            lines += [f"- s{i}", f"gene_tree = {t}", f"mapping = {m}"]
        ff = os.path.join(wd, f"{name}_{tag}.txt"); open(ff, "w").write("\n".join(lines) + "\n")
        for fw in (True,):
            out = os.path.join(wd, f"out_{name}_{tag}_{int(fw)}")
            cmd = [kalerax, "-f", ff, "-s", sp, "--rec-model", "UndatedDL",
                   "--gene-tree-rooting", "UNIFORM", "--species-tree-search", "SKIP", "--summary-only",
                   "--gene-tree-samples", "0", "--fix-rates", "--d", str(c["D"]), "--l", str(c["L"]),
                   "--min-covered-species", "1", "--skip-family-filtering", "--seed", "1", "-p", out]
            if c["q"] is not None:
                cmd += ["--wgd", f"A,B:{c['q']}", "--fix-wgd-q"]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            ll = dict(l.split() for l in open(os.path.join(out, "per_fam_likelihoods.txt")) if l.strip())
            for i, k in enumerate(top):
                if f"s{i}" in ll:
                    X.setdefault((fw, k), {})[tag] = (freq[k] / N / np.exp(float(ll[f"s{i}"])), freq[k])
    res = {}
    for fw in (True,):
        z = [np.log(v["b"][0] / v["a"][0]) / np.sqrt(1 / v["b"][1] + 1 / v["a"][1])
             for (f, k), v in X.items() if f == fw and len(v) == 2]
        res[fw] = (float(np.sum(np.square(z))), len(z))
    return res


if __name__ == "__main__":
    kalerax, wd = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 8000000
    os.makedirs(wd, exist_ok=True)
    sp = os.path.join(wd, "sp.nw"); open(sp, "w").write(SP + "\n")
    tree = SpeciesTree(SP)
    from scipy.stats import chi2 as chi2d
    ok = True
    for name, pair in PAIRS.items():
        c1, n1 = run_pair(kalerax, wd, N, name, pair, tree, sp)[True]
        p1 = chi2d.sf(c1, n1)
        print(f"  [{name}] forward-process weights: chi2 = {c1:.1f} on {n1} shapes (p = {p1:.2g})")
        ok &= n1 >= 15 and p1 > 1e-3
    print("PASS: UndatedDL duplication weights concordant with forward simulation" if ok
          else "FAIL: UndatedDL dup-weight concordance")
    sys.exit(0 if ok else 1)
