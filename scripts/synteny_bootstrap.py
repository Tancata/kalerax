#!/usr/bin/env python3
"""Parametric-bootstrap calibration of the synteny-aware WGD likelihood-ratio test.

Why: with synteny, the WGD-vs-null LRT is non-standard. Under H0 (q = 0) sigma_W
is not identifiable and q sits on its boundary (a Davies-type problem). The
pair-level synteny emission is also a composite likelihood that ignores the
dependence of pairs within a collinear block, so the statistic can be
overconfident. This script simulates data under the fitted null model (no
WGD), refits the null and the WGD model exactly as for the real data, and
returns the bootstrap null distribution of 2*dlnL and a p-value.

Typical use (after fitting the real data):
  python3 scripts/synteny_bootstrap.py \
      --kalerax build/bin/kalerax --species-tree sp.nw \
      --null-run OBS_NULL_DIR --wgd-run OBS_WGD_DIR \
      --wgd Glycine_max,Glycine_falcata --min-sigmaW 0.3 \
      --reps 100 --workers 16 --out boot_n2 [--block-frac 0.5] [--families N]

  * --null-run / --wgd-run: output directories (or their alerax.log files) of
    the observed fits with --synteny-dir. The null log provides D, L, T (from
    model_parameters/model_parameters.txt or the log) and sigma_N; both logs
    provide the observed lnL.
  * The simulated families use the null rates, UNIFORM origination with
    survival conditioning and the kalerax >=4-species filter, which matches
    what the likelihood assumes. --families defaults to the number of
    families kalerax kept in the observed null run.
  * Synteny: every same-genome pair whose LCA is a duplication or transfer
    gets y = 1 with prob sigma_N. A fraction --block-frac of (gene-tree node,
    genome) groups share one draw, mimicking collinear SSD blocks
    (segmental duplications). --na-frac sets the fraction of NA genes.
    Choose --block-frac as a sensitivity parameter; 0 is the composite-
    likelihood assumption itself.
  * Gene trees are simulated without inference error (true trees). The
    observed data use UFBoot sets, so treat the bootstrap as approximate in
    that respect; see the gene-tree-error results in the simulation study.

--null-model aore calibrates the LORe-vs-AORe test instead: data are simulated
under the fitted AORe model (WGD with q-hat, r = 1; WGD-split pairs collinear
with sigma_W-hat, other pairs with sigma_N-hat, same --block-frac), and each
replicate is fitted with AORe (null) and LORe (alternative). Then --null-run is
the observed AORe fit and --wgd-run the observed LORe fit.

Gamma rate heterogeneity: if the observed --null-run was fitted with
--speciation-gamma-categories K (its model_parameters/wgd_parameters.txt has a
'gamma alpha' row), replicates are simulated with the same K categories at the
fitted alpha (one category per family, kalerax_sim.simulate_families) and every
replicate fit adds --speciation-gamma-categories K. --no-gamma disables this.

Each replicate's files are written, fitted and deleted before the next one
starts, so disk use stays at about `workers` replicates.
"""
import argparse, json, os, re, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kalerax_sim import SpeciesTree, Params, simulate_families, write_family_set, synteny_lines


def read_log(path):
    log = path if path.endswith(".log") else os.path.join(path, "alerax.log")
    txt = open(log).read()
    ll = float(re.findall(r"Final species tree likelihood: ll=(-?[0-9.eE+-]+)", txt)[-1])
    sig = re.findall(r"sigmaW=([0-9.eE+-]+) sigmaN=([0-9.eE+-]+)", txt)
    sW, sN = (float(sig[-1][0]), float(sig[-1][1])) if sig else (None, None)
    kept = re.findall(r"Remaining families: (\d+)", txt)
    qs = re.findall(r"q\[node \d+\]=([0-9.eE+-]+)", txt)
    rates = None
    mp = os.path.join(os.path.dirname(log), "model_parameters", "model_parameters.txt")
    if not os.path.exists(mp):
        mp = os.path.join(os.path.dirname(log), "model_parameters.txt")
    if os.path.exists(mp):
        rates = [float(x) for x in open(mp).readlines()[1].split()[1:4]]
    gamma = None  # (alpha, categories) from the full-precision parameter file
    wp = os.path.join(os.path.dirname(log), "model_parameters", "wgd_parameters.txt")
    if os.path.exists(wp):
        for l in open(wp):
            f = l.rstrip("\n").split("\t")
            if len(f) == 5 and f[0] == "gamma" and f[1] == "alpha":
                gamma = (float(f[3]), int(float(f[4])))
    return dict(lnL=ll, sigmaW=sW, sigmaN=sN, kept=int(kept[-1]) if kept else None, rates=rates,
                q=float(qs[-1]) if qs else None, gamma=gamma)


def fit(args, famfile, syndir, out, wgd, lore=None):
    cmd = [args.kalerax, "-f", famfile, "-s", args.species_tree, "--rec-model", "UndatedDTL",
           "--transfer-constraint", "PARENTS", "--gene-tree-rooting", "UNIFORM",
           "--species-tree-search", "SKIP", "--summary-only", "--gene-tree-samples", "0",
           "--seed", "1", "--synteny-dir", syndir, "-p", out]
    if wgd:
        cmd += ["--wgd", args.wgd]
        if args.min_sigmaW > 0:
            cmd += ["--synteny-min-sigmaW", str(args.min_sigmaW)]
        if (args.lore if lore is None else lore):
            cmd += ["--lore"]
    if args.gamma:
        cmd += ["--speciation-gamma-categories", str(args.gamma[1])]
    cmd += args.extra
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    r = read_log(out)
    pf = os.path.join(out, "per_fam_likelihoods.txt")
    if os.path.exists(pf):  # exact lnL (the log prints 6 significant figures)
        r["lnL"] = sum(float(l.split()[1]) for l in open(pf) if l.strip())
    q = re.findall(r"q\[node \d+\]=([0-9.eE+-]+)", open(os.path.join(out, "alerax.log")).read())
    r["q"] = float(q[-1]) if (wgd and q) else None
    return r


def replicate(args, tree, params, b):
    wd = os.path.join(args.out, f"rep{b}")
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    rng = np.random.default_rng(args.seed + 7919 * b)
    ga, gk = args.gamma if args.gamma else (None, 1)
    fams, _ = simulate_families(tree, params, args.families, args.seed + b,
                                origin_strategy="UNIFORM", keep_nodes=True,
                                gamma_alpha=ga, gamma_categories=gk)
    ff = write_family_set(fams, wd)
    syn = os.path.join(wd, "syn"); os.makedirs(syn)
    for i, f in enumerate(fams):
        lines = synteny_lines(f["root"], rng, sigma_w=args.sim_sigmaW, sigma_n=args.sigmaN,
                              na_frac=args.na_frac, block_frac=args.block_frac)
        if lines:
            open(os.path.join(syn, f"fam{i}.syn"), "w").write("\n".join(lines) + "\n")
    if args.null_model == "aore":  # AORe (null) vs LORe (alternative)
        r0 = fit(args, ff, syn, os.path.join(wd, "null"), True, lore=False)
        r1 = fit(args, ff, syn, os.path.join(wd, "wgd"), True, lore=True)
    else:  # no WGD (null) vs WGD (alternative)
        r0 = fit(args, ff, syn, os.path.join(wd, "null"), False)
        r1 = fit(args, ff, syn, os.path.join(wd, "wgd"), True)
    shutil.rmtree(wd, ignore_errors=True)
    return dict(rep=b, lnL_null=r0["lnL"], lnL_wgd=r1["lnL"], stat=2 * (r1["lnL"] - r0["lnL"]),
                q=r1["q"], sigmaW=r1["sigmaW"], sigmaN_null=r0["sigmaN"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kalerax", required=True)
    ap.add_argument("--species-tree", required=True)
    ap.add_argument("--null-run", required=True)
    ap.add_argument("--wgd-run", required=True)
    ap.add_argument("--wgd", required=True, help="LABEL[,LABEL2] as passed to kalerax --wgd")
    ap.add_argument("--min-sigmaW", type=float, default=0.3)
    ap.add_argument("--lore", action="store_true")
    ap.add_argument("--null-model", choices=["none", "aore"], default="none",
                    help="none: simulate without WGD, test WGD vs none (default); "
                         "aore: simulate under the fitted AORe model, test LORe vs AORe")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--families", type=int, default=None)
    ap.add_argument("--block-frac", type=float, default=0.0)
    ap.add_argument("--na-frac", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-gamma", action="store_true",
                    help="ignore a gamma fit in --null-run (simulate and fit without gamma)")
    ap.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                    help="extra kalerax arguments for every fit (put last)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    obs0, obs1 = read_log(args.null_run), read_log(args.wgd_run)
    if obs0["rates"] is None or obs0["sigmaN"] is None:
        sys.exit("could not read null rates / sigma_N from --null-run")
    args.sigmaN = obs0["sigmaN"]
    args.gamma = None if args.no_gamma else obs0["gamma"]
    if any(a == "--speciation-gamma-categories" for a in args.extra):
        sys.exit("pass gamma through the --null-run fit, not --extra")
    args.families = args.families or obs0["kept"]
    tree = SpeciesTree(open(args.species_tree).read())
    D, L, T = obs0["rates"]
    if args.null_model == "aore":
        if obs0["q"] is None or obs0["sigmaW"] is None:
            sys.exit("could not read q / sigma_W from the AORe --null-run")
        e = tree.lca(args.wgd.split(","))
        params = Params(tree, D, L, T, wgds=[(e, obs0["q"], 1.0)])
        args.sim_sigmaW = obs0["sigmaW"]
    else:
        params = Params(tree, D, L, T)
        args.sim_sigmaW = args.sigmaN
    obs_stat = 2 * (obs1["lnL"] - obs0["lnL"])
    print(f"observed 2dlnL = {obs_stat:.2f}; null rates D={D} L={L} T={T}, sigma_N={args.sigmaN}, "
          f"families={args.families}, block_frac={args.block_frac}, gamma={args.gamma}", flush=True)
    res = []
    outp = os.path.join(args.out, "bootstrap.jsonl")
    with ThreadPoolExecutor(max_workers=args.workers) as ex, open(outp, "a") as fh:
        for r in ex.map(lambda b: replicate(args, tree, params, b), range(args.reps)):
            res.append(r)
            fh.write(json.dumps(r) + "\n"); fh.flush()
    stats = np.array([r["stat"] for r in res])
    p = (1 + np.sum(stats >= obs_stat)) / (1 + len(stats))
    summ = dict(observed_stat=obs_stat, reps=len(stats), p_value=float(p),
                null_quantiles={k: float(np.quantile(stats, v)) for k, v in
                                [("50%", .5), ("90%", .9), ("95%", .95), ("99%", .99), ("max", 1.0)]},
                settings=dict(null_model=args.null_model, q_null=obs0["q"], sigmaW_sim=args.sim_sigmaW,
                              wgd=args.wgd, min_sigmaW=args.min_sigmaW, block_frac=args.block_frac,
                              na_frac=args.na_frac, families=args.families, sigmaN=args.sigmaN,
                              rates=dict(D=D, L=L, T=T),
                              gamma=dict(alpha=args.gamma[0], categories=args.gamma[1]) if args.gamma else None))
    json.dump(summ, open(os.path.join(args.out, "summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
