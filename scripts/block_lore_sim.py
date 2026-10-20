#!/usr/bin/env python3
"""Simulation check for block-coupled LORe (scripts/block_lore.py).

One replicate: simulate n_blocks syntenic blocks on the acipenseriform tree
(WGD on the Polyodon/Acipenser LCA branch) with block-shared resolution
classes (scripts/kalerax_sim.simulate_block_families), emission-level synteny,
then: free kalerax --lore fit -> fixed-parameter grid of per-family lnL at
r in GRID and at r-hat -> block / independent mixtures, permutation test,
HMM along the simulated gene order (blocks are contiguous). coupling=0 is the
no-block-structure null (permutation p must be ~ Uniform(0,1)).

usage: block_lore_sim.py <kalerax> <species_tree> <workdir> <seed> <coupling>
       [n_blocks=160] [mean_block_size=8] [n_perm=200]
Writes <workdir>/result.json and deletes the kalerax run directories.
"""
import json, os, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kalerax_sim as ks
import block_lore as bl

WGD = "Polyodon_spathula,Acipenser_ruthenus"
GRID = [0.0, 0.25, 0.5, 0.75, 1.0]
R_CLASSES, WEIGHTS = [1.0, 0.25], [0.45, 0.55]  # early (before split), late
BASE = ["--rec-model", "UndatedDTL", "--transfer-constraint", "PARENTS", "--gene-tree-rooting",
        "UNIFORM", "--species-tree-search", "SKIP", "--summary-only", "--gene-tree-samples", "0",
        "--seed", "1", "--synteny-min-sigmaW", "0.3"]


def read_par(out):
    par = {}
    for l in open(os.path.join(out, "model_parameters", "wgd_parameters.txt")):
        if l.startswith("#"):
            continue
        k, lab, _, v1, v2 = l.rstrip("\n").split("\t")
        par[(k, lab)] = (float(v1), None if v2 == "NA" else float(v2))
    return par


def fixed_extra(par, r):
    q = [v[0] for (k, _), v in par.items() if k == "wgd"][0]
    rate = {lab: v[0] for (k, lab), v in par.items() if k == "rate"}
    sw, sn = par[("synteny", "sigma")]
    return ["--fix-rates", "--d", repr(rate["D"]), "--l", repr(rate["L"]), "--t", repr(rate["T"]),
            "--wgd", f"{WGD}:{q!r}", "--fix-wgd-q", "--fix-synteny", "--synteny-params",
            f"{sw!r},{sn!r}", "--skip-thorough-rates", "--lore", "--lore-fixed-r", repr(r)]


def main():
    kx, spf, wd, seed, coupling = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5])
    nb = int(sys.argv[6]) if len(sys.argv) > 6 else 160
    mbs = float(sys.argv[7]) if len(sys.argv) > 7 else 8
    nperm = int(sys.argv[8]) if len(sys.argv) > 8 else 200
    os.makedirs(wd, exist_ok=True)
    spt = open(spf).read().strip()
    t = ks.SpeciesTree(spt); e = t.lca(WGD.split(","))
    fams, _ = ks.simulate_block_families(t, 0.0574, 0.0922, 0.0173, e, 0.9, R_CLASSES, WEIGHTS,
                                         n_blocks=nb, mean_block_size=mbs, seed=seed, coupling=coupling)
    data = os.path.join(wd, "data")
    ff = ks.write_family_set(fams, data)
    syn = os.path.join(data, "syn"); os.makedirs(syn, exist_ok=True)
    rng = np.random.default_rng(seed + 100000)
    for i, f in enumerate(fams):
        L = ks.synteny_lines(f["root"], rng, 0.95, 0.03, 0.1)
        if L:
            open(os.path.join(syn, f"fam{i}.syn"), "w").write("\n".join(L) + "\n")
    sp = os.path.join(wd, "sp.nw"); open(sp, "w").write(spt + "\n")

    def run(out, extra):
        cmd = [kx, "-f", ff, "-s", sp, "--synteny-dir", syn, "-p", out] + BASE + extra
        p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert p.returncode == 0, extra

    run(os.path.join(wd, "free"), ["--wgd", WGD, "--lore"])
    par = read_par(os.path.join(wd, "free")); rhat = par[[k for k in par if k[0] == "wgd"][0]][1]
    rs = GRID + [rhat]
    with ThreadPoolExecutor(len(rs)) as ex:
        list(ex.map(lambda r: run(os.path.join(wd, f"grid_{r:g}"), fixed_extra(par, r)), rs))
    df, r = bl.load_grid([f"{x:g}={wd}/grid_{x:g}" for x in GRID])
    rh = (rhat, bl.read_per_fam(f"{wd}/grid_{rhat:g}"))
    names = [f"fam{i}" for i in range(len(fams))]
    f2b = {n: f"b{f['block']}" for n, f in zip(names, fams)}
    order = pd.DataFrame({"family": names, "chrom": "1", "position": np.arange(len(names))})
    res, blk, fam = bl.analyse(df, r, f2b, order, rh, n_perm=nperm, seed=seed)
    truth = {f"b{f['block']}": R_CLASSES[f["block_class"]] for f in fams}
    blk["r_true_block"] = blk["block"].map(truth)
    r_of = {n: R_CLASSES[f["r_class"]] for n, f in zip(names, fams)}
    rt = fam["family"].map(r_of).values  # fam rows follow the grid index order
    res["sim"] = dict(seed=seed, coupling=coupling, n_blocks=nb, mean_block_size=mbs,
                      r_classes=R_CLASSES, weights=WEIGHTS, free_fit=dict(
                          q=par[[k for k in par if k[0] == "wgd"][0]][0], r=rhat,
                          D=par[("rate", "D")][0], L=par[("rate", "L")][0], T=par[("rate", "T")][0],
                          sigmaW=par[("synteny", "sigma")][0], sigmaN=par[("synteny", "sigma")][1]),
                      block_MAP_acc=float((blk["MAP_r"] == blk["r_true_block"]).mean()),
                      block_E_r_corr=float(np.corrcoef(blk["E_r"], blk["r_true_block"])[0, 1]),
                      fam_hmm_E_r_corr=float(np.corrcoef(fam["E_r_hmm"], rt)[0, 1]),
                      fam_indep_E_r_corr=float(np.corrcoef(fam["E_r_indep"], rt)[0, 1]),
                      fam_block_E_r_corr=float(np.corrcoef(fam["E_r_block"], rt)[0, 1]))
    json.dump(res, open(os.path.join(wd, "result.json"), "w"), indent=1)
    for d in os.listdir(wd):
        if d in ("free", "data") or d.startswith("grid_"):
            shutil.rmtree(os.path.join(wd, d))
    print(json.dumps(res["sim"]))


if __name__ == "__main__":
    main()
