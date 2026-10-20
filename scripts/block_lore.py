#!/usr/bin/env python3
"""Block-coupled LORe: do families in the same syntenic block rediploidize together?

Input: per-family log-likelihoods L_f(r_c) on a grid of LORe resolution values
r_c, from kalerax runs at otherwise fixed parameters (rates, q, sigmas held at
the free LORe fit; see DTL_LORE.md, "Fixed-parameter evaluation"), plus a
family -> block map (and, for the HMM, a gene order).

Models (all on the same L_f(r_c) matrix, so likelihoods are comparable):
  single    one genome-wide r: sum_f log L_f(r_hat) (r_hat column, or best grid r)
  indep     mixture over r classes, every family its own class draw:
            sum_f log sum_c pi_c L_f(r_c)
  block     mixture over r classes shared within a block:
            sum_b log sum_c pi_c prod_{f in b} L_f(r_c)
            (families without a block are singletons)
  hmm       hidden class along gene order; switch with prob rho to a class drawn
            from pi: A = (1 - rho) I + rho 1 pi^T  (block boundaries inferred)
Evidence for block coupling: lnL(block) - lnL(indep), calibrated by permuting
block labels across families (block sizes preserved). Rates/q/sigma are held
at the free-fit values, which is an approximation.

P(resolved before the WGD subtree's first speciation) = r on the WGD branch,
so a block's timing summary is E[r | block] = sum_c P(c | block) r_c.

usage:
  block_lore.py --grid 0=dir0 0.25=dir1 ... [--rhat 0.443=dirhat]
                --blocks fam2block.tsv [--order order.tsv] [--perm 1000]
                --out prefix
  fam2block.tsv: family<TAB>block  (blank/NA block = no block)
  order.tsv:     family<TAB>chrom<TAB>position  (reference-genome gene order)
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.optimize import minimize


# --------------------------------------------------------------------- input
def read_per_fam(path):
    d = pd.read_csv(os.path.join(path, "per_fam_likelihoods.txt") if os.path.isdir(path) else path,
                    sep=r"\s+", header=None, names=["family", "ll"], dtype={"family": str},
                    keep_default_na=False)
    return d.set_index("family")["ll"]


def load_grid(spec):
    """spec: list of 'r=path'. Returns (DataFrame families x r, r array)."""
    cols = {}
    for s in spec:
        r, p = s.split("=", 1)
        cols[float(r)] = read_per_fam(p)
    df = pd.DataFrame(cols).sort_index(axis=1)
    if df.isna().any().any():
        raise SystemExit("grid runs do not share the same families")
    return df, np.array(df.columns, float)


def block_index(families, fam2block):
    """List of integer index arrays, one per block; unblocked families are singletons."""
    pos = {f: i for i, f in enumerate(families)}
    groups, single = {}, []
    for f in families:
        b = fam2block.get(f)
        if b is None or b == "" or b == "NA":
            single.append(pos[f])
        else:
            groups.setdefault(b, []).append(pos[f])
    blocks = [np.array(v) for v in groups.values()] + [np.array([i]) for i in single]
    names = list(groups.keys()) + [None] * len(single)
    return blocks, names


# ------------------------------------------------------------------ mixtures
def block_sums(LL, blocks):
    return np.vstack([LL[b].sum(0) for b in blocks])


def em_mixture(B, tol=1e-10, max_iter=10000, pi0=None):
    """Maximise sum_b log sum_c pi_c exp(B[b,c]). Returns lnL, pi, posterior."""
    nb, C = B.shape
    pi = np.full(C, 1.0 / C) if pi0 is None else np.asarray(pi0, float)
    old = -np.inf
    for _ in range(max_iter):
        lp = np.log(np.maximum(pi, 1e-300)) + B
        z = logsumexp(lp, axis=1)
        ll = z.sum()
        post = np.exp(lp - z[:, None])
        pi = post.mean(0)
        if ll - old < tol:
            break
        old = ll
    lp = np.log(np.maximum(pi, 1e-300)) + B
    z = logsumexp(lp, axis=1)
    return z.sum(), pi, np.exp(lp - z[:, None])


def fit_block(LL, blocks):
    return em_mixture(block_sums(LL, blocks))


def permutation_test(LL, blocks, n_perm, seed=1, pi0=None):
    """Shuffle block membership among blocked families (block sizes kept)."""
    rng = np.random.default_rng(seed)
    multi = [b for b in blocks if len(b) > 1]
    single = [b for b in blocks if len(b) == 1]
    members = np.concatenate(multi) if multi else np.array([], int)
    sizes = np.cumsum([len(b) for b in multi])[:-1]
    Bs = block_sums(LL, single) if single else np.zeros((0, LL.shape[1]))
    out = np.empty(n_perm)
    for k in range(n_perm):
        perm = np.split(rng.permutation(members), sizes)
        B = np.vstack([block_sums(LL, perm), Bs])
        out[k] = em_mixture(B, tol=1e-8, pi0=pi0)[0]
    return out


# ----------------------------------------------------------------------- HMM
def hmm_forward(LL, segments, pi, rho):
    """lnL of the switching HMM; segments = list of index arrays in gene order
    (one per chromosome/scaffold), each started from pi."""
    lpi = np.log(np.maximum(pi, 1e-300))
    l1r = np.log1p(-rho)
    lr = np.log(rho)
    tot = 0.0
    for seg in segments:
        a = lpi + LL[seg[0]]
        for i in seg[1:]:
            za = logsumexp(a)
            a = np.logaddexp(l1r + a, lr + lpi + za) + LL[i]
        tot += logsumexp(a)
    return tot


def hmm_posterior(LL, segments, pi, rho):
    """Per-family posterior class probabilities (forward-backward)."""
    lpi = np.log(np.maximum(pi, 1e-300))
    l1r, lr = np.log1p(-rho), np.log(rho)
    post = np.zeros_like(LL)
    for seg in segments:
        n = len(seg)
        F = np.empty((n, LL.shape[1])); Bk = np.zeros_like(F)
        F[0] = lpi + LL[seg[0]]
        for k in range(1, n):
            F[k] = np.logaddexp(l1r + F[k - 1], lr + lpi + logsumexp(F[k - 1])) + LL[seg[k]]
        for k in range(n - 2, -1, -1):
            e = LL[seg[k + 1]] + Bk[k + 1]
            Bk[k] = np.logaddexp(l1r + e, lr + logsumexp(lpi + e))
        z = logsumexp(F[-1])
        post[seg] = np.exp(F + Bk - z)
    return post


def fit_hmm(LL, segments, pi0):
    C = LL.shape[1]
    def unpack(x):
        rho = 1.0 / (1.0 + np.exp(-x[0]))
        w = np.concatenate([[0.0], x[1:]]); pi = np.exp(w - logsumexp(w))
        return pi, rho
    x0 = np.concatenate([[np.log(0.1 / 0.9)], np.log(np.maximum(pi0[1:], 1e-6) / max(pi0[0], 1e-6))])
    res = minimize(lambda x: -hmm_forward(LL, segments, *unpack(x)), x0, method="Nelder-Mead",
                   options=dict(maxiter=4000, xatol=1e-5, fatol=1e-4))
    pi, rho = unpack(res.x)
    return -res.fun, pi, rho


# ----------------------------------------------------------------------- run
def analyse(df, r, fam2block=None, order=None, rhat=None, n_perm=1000, seed=1):
    """df: families x r grid of per-family lnL. rhat: Series of per-family lnL at the
    free-fit r (single-r model), optional. Returns (summary dict, per-block table,
    per-family table)."""
    fams = list(df.index)
    LL = df.values
    res = {"n_families": len(fams), "r_grid": list(map(float, r))}
    tot = LL.sum(0)
    res["single_grid"] = {"r": float(r[np.argmax(tot)]), "lnL": float(tot.max())}
    if rhat is not None:
        res["single_rhat"] = {"r": float(rhat[0]), "lnL": float(rhat[1].loc[fams].sum())}
    li, pii, posti = em_mixture(LL)
    res["indep"] = {"lnL": float(li), "pi": pii.tolist(), "df": len(r) - 1}
    fam_tab = pd.DataFrame({"family": fams, "E_r_indep": posti @ r})
    blk_tab = None
    if fam2block is not None:
        blocks, names = block_index(fams, fam2block)
        multi = [(b, nm) for b, nm in zip(blocks, names) if len(b) > 1]
        mb = [b for b, _ in multi]
        inb = np.concatenate(mb) if mb else np.array([], int)
        out = np.setdiff1d(np.arange(len(fams)), inb)
        # Blocked and unblocked families are different populations (block membership
        # is defined from the same synteny data that enters L_f), so each subset
        # gets its own class weights. The coupling test compares block-shared vs
        # independent mixtures on the blocked families only; the unblocked
        # families contribute the same independent-mixture term to both models.
        lu, piu, postu = em_mixture(LL[out]) if len(out) else (0.0, np.zeros(len(r)), None)
        lib, pib_i, postib = em_mixture(LL[inb])
        lb, pib, postb = em_mixture(block_sums(LL, mb))
        res["unblocked"] = {"n_families": int(len(out)), "indep_lnL": float(lu), "pi": piu.tolist()}
        res["block"] = {"n_blocks": len(mb), "n_blocked_families": int(len(inb)),
                        "indep_lnL_blocked": float(lib), "indep_pi_blocked": pib_i.tolist(),
                        "block_lnL_blocked": float(lb), "pi": pib.tolist(),
                        "lnL": float(lb + lu), "indep_split_lnL": float(lib + lu),
                        "gain_vs_indep": float(lb - lib)}
        if n_perm:
            perm = permutation_test(LL, mb, n_perm, seed, pi0=pib)
            res["block"]["perm"] = {"n": n_perm, "mean": float(perm.mean()), "sd": float(perm.std()),
                                    "max": float(perm.max()),
                                    "p": float((1 + np.sum(perm >= lb)) / (n_perm + 1)),
                                    "gain_vs_perm_mean": float(lb - perm.mean()),
                                    "z": float((lb - perm.mean()) / max(perm.std(), 1e-12))}
        rows = []
        for (b, nm), pb in zip(multi, postb):
            rows.append(dict(block=nm, n_families=len(b), E_r=float(pb @ r),
                             MAP_r=float(r[np.argmax(pb)]), P_MAP=float(pb.max()),
                             lnL_shared=float(logsumexp(np.log(np.maximum(pib, 1e-300)) + LL[b].sum(0))),
                             lnL_indep=float(logsumexp(np.log(np.maximum(pib_i, 1e-300)) + LL[b], axis=1).sum()),
                             **{f"P_r{rc:g}": float(x) for rc, x in zip(r, pb)}))
        blk_tab = pd.DataFrame(rows)
        f2E = {}
        for (b, _), pb in zip(multi, postb):
            for i in b:
                f2E[fams[i]] = float(pb @ r)
        fam_tab["E_r_block"] = fam_tab["family"].map(f2E)
        fam_tab["block"] = fam_tab["family"].map(lambda f: fam2block.get(f, ""))
    if order is not None:
        o = order[order["family"].isin(fams)].sort_values(["chrom", "position"])
        pos = {f: i for i, f in enumerate(fams)}
        segs = [np.array([pos[f] for f in g["family"]]) for _, g in o.groupby("chrom", sort=False)]
        segs = [s for s in segs if len(s)]
        inseg = set(np.concatenate(segs)) if segs else set()
        rest = [i for i in range(len(fams)) if i not in inseg]
        lh, pih, rho = fit_hmm(LL, segs, pii)
        # families without coordinates contribute an independent mixture with their
        # own class weights (as for unblocked families in the block model)
        lrest = float(em_mixture(LL[rest])[0]) if rest else 0.0
        li_seg = float(logsumexp(np.log(np.maximum(pii, 1e-300)) + LL[list(inseg)], axis=1).sum())
        res["hmm"] = {"lnL_ordered": float(lh), "lnL_total": float(lh + lrest), "pi": pih.tolist(),
                      "rho": float(rho), "mean_run_length": float(1 / rho),
                      "n_ordered_families": len(inseg), "n_segments": len(segs),
                      "indep_lnL_same_families": li_seg, "gain_vs_indep": float(lh - li_seg)}
        ph = hmm_posterior(LL, segs, pih, rho)
        fam_tab["E_r_hmm"] = np.where(np.isin(np.arange(len(fams)), list(inseg)), ph @ r, np.nan)
    return res, blk_tab, fam_tab


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", nargs="+", required=True)
    ap.add_argument("--rhat", default=None, help="r=dir of the run at the free-fit r")
    ap.add_argument("--blocks", default=None)
    ap.add_argument("--order", default=None)
    ap.add_argument("--perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    df, r = load_grid(a.grid)
    rhat = None
    if a.rhat:
        rv, p = a.rhat.split("=", 1); rhat = (float(rv), read_per_fam(p))
    f2b = None
    if a.blocks:
        t = pd.read_csv(a.blocks, sep="\t", header=None, names=["family", "block"], dtype=str,
                        keep_default_na=False)
        f2b = dict(zip(t["family"], t["block"]))
    order = None
    if a.order:
        order = pd.read_csv(a.order, sep="\t", header=None, names=["family", "chrom", "position"],
                            dtype={"family": str, "chrom": str}, keep_default_na=False)
    res, blk, fam = analyse(df, r, f2b, order, rhat, a.perm, a.seed)
    json.dump(res, open(a.out + "_summary.json", "w"), indent=1)
    if blk is not None:
        blk.to_csv(a.out + "_blocks.tsv", sep="\t", index=False)
    fam.to_csv(a.out + "_families.tsv", sep="\t", index=False)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
