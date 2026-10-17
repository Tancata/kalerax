#!/usr/bin/env python3
"""Regression gate for the UndatedDTL duplication weighting (default = forward
process; --ale-dup-weights = ALE/AleRax convention).

Forward-simulates gene families (scripts/kalerax_sim.py) under two UndatedDTL
parameter sets that differ only in the transfer rate, collapses every family
with 3-6 genes to its unrooted species-labelled shape, and scores one
representative per frequent shape with kalerax at FIXED rates
(--fix-rates, per_fam_likelihoods.txt). For an exact likelihood (up to a
parameter-independent per-tree constant) the ratio X = freq(shape)/L(shape)
must be the same under both parameter sets. With the default (forward-process)
weighting this holds (chi2/df ~ 1); with --ale-dup-weights it does not, because a
duplication split is counted once while a transfer split sums both orderings.

usage: dtl_dup_weight_concordance.py <kalerax> <workdir> [N]
exit 0 = PASS (default concordant AND --ale-dup-weights discordant), 1 = FAIL
"""
import itertools, os, subprocess, sys
from collections import Counter, defaultdict
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from kalerax_sim import SpeciesTree, Params, simulate_families

SP = "(O:1,((A:1,B:1):1,(C:1,D:1):1):1);"
CONFIGS = {"lowT": (0.25, 0.35, 0.05), "highT": (0.25, 0.35, 0.30)}


def parse(nwk):
    s = nwk.strip().rstrip(";"); stack, clades, leaves, i = [[]], [], [], 0
    while i < len(s):
        c = s[i]
        if c == "(":
            stack.append([]); i += 1
        elif c == ",":
            i += 1
        elif c == ")":
            cl = frozenset().union(*stack.pop()); clades.append(cl); stack[-1].append(cl); i += 1
        else:
            j = i
            while j < len(s) and s[j] not in ",()":
                j += 1
            leaves.append(s[i:j]); stack[-1].append(frozenset([s[i:j]])); i = j
    return leaves, clades


def shape_key(nwk):
    leaves, clades = parse(nwk)
    allL = frozenset(leaves); ref = min(leaves)
    sps = {cl if ref not in cl else allL - cl for cl in clades}
    sps = {x for x in sps if 1 < len(x) < len(allL) - 1}
    groups = defaultdict(list)
    for l in leaves:
        groups[l.split("_")[0]].append(l)
    names = sorted(groups); best = None
    for perm in itertools.product(*[itertools.permutations(groups[s]) for s in names]):
        ren = {l: f"{s}{k}" for s, pl in zip(names, perm) for k, l in enumerate(pl)}
        A = frozenset(ren.values()); r0 = min(A); rep = []
        for side in sps:
            s2 = frozenset(ren[l] for l in side)
            rep.append(tuple(sorted(A - s2 if r0 in s2 else s2)))
        rep = tuple(sorted(rep))
        best = rep if best is None or rep < best else best
    return (tuple(sorted(Counter(l.split("_")[0] for l in leaves).items())), best)


def run(kalerax, wd, N):
    os.makedirs(wd, exist_ok=True)
    sp = os.path.join(wd, "sp.nw"); open(sp, "w").write(SP + "\n")
    tree = SpeciesTree(SP)
    X = {}
    for ci, (tag, (D, L, T)) in enumerate(CONFIGS.items()):
        fams, _ = simulate_families(tree, Params(tree, D, L, T), N, 101 + ci, keep_all=True)
        freq, rep = Counter(), {}
        for f in fams:
            if 3 <= f["n_genes"] <= 6:
                k = shape_key(f["newick"]); freq[k] += 1; rep.setdefault(k, f)
        top = [k for k, c in freq.most_common() if c >= 150]
        lines = ["[FAMILIES]"]
        for i, k in enumerate(top):
            f = rep[k]; t = os.path.join(wd, f"{tag}_{i}.nw"); m = os.path.join(wd, f"{tag}_{i}.link")
            open(t, "w").write(f["newick"] + "\n")
            open(m, "w").write("".join(f"{g}\t{s}\n" for g, s in f["mapping"]))
            lines += [f"- s{i}", f"gene_tree = {t}", f"mapping = {m}"]
        ff = os.path.join(wd, f"{tag}.txt"); open(ff, "w").write("\n".join(lines) + "\n")
        for flag in (False, True):
            out = os.path.join(wd, f"out_{tag}_{int(flag)}")
            cmd = [kalerax, "-f", ff, "-s", sp, "--rec-model", "UndatedDTL", "--transfer-constraint", "PARENTS",
                   "--gene-tree-rooting", "UNIFORM", "--species-tree-search", "SKIP", "--summary-only",
                   "--gene-tree-samples", "0", "--fix-rates", "--d", str(D), "--l", str(L), "--t", str(T),
                   "--min-covered-species", "1", "--skip-family-filtering", "--seed", "1", "-p", out]
            if not flag:
                cmd.append("--ale-dup-weights")
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            ll = dict(l.split() for l in open(os.path.join(out, "per_fam_likelihoods.txt")) if l.strip())
            for i, k in enumerate(top):
                if f"s{i}" in ll:
                    X.setdefault((flag, k), {})[tag] = (freq[k] / N / np.exp(float(ll[f"s{i}"])), freq[k])
    res = {}
    for flag in (False, True):
        z = [np.log(v["highT"][0] / v["lowT"][0]) / np.sqrt(1 / v["highT"][1] + 1 / v["lowT"][1])
             for (fl, k), v in X.items() if fl == flag and len(v) == 2]
        res[flag] = (float(np.sum(np.square(z))), len(z))
    return res


if __name__ == "__main__":
    kalerax, wd = sys.argv[1], sys.argv[2]
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 600000
    res = run(kalerax, wd, N)
    (c0, n0), (c1, n1) = res[False], res[True]
    print(f"  --ale-dup-weights:              chi2 = {c0:.1f} on {n0} shapes")
    print(f"  default (forward-process):      chi2 = {c1:.1f} on {n1} shapes")
    ok = n1 >= 20 and c1 / n1 < 1.8 and c0 / n0 > 4.0
    print("PASS: default duplication weights are concordant with forward simulation"
          if ok else "FAIL: dup-weight concordance")
    sys.exit(0 if ok else 1)
