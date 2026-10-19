#!/usr/bin/env python3
"""Regression gate for --synteny-min-sigmaW.

Simulates (scripts/kalerax_sim.py) an 11-taxon data set with a burst of
small-scale duplications on one internal branch and no WGD, plus emission-level
synteny (WGD-split pairs collinear w.p. 0.7, other duplicate pairs w.p. 0.03).
  - free sigma_W: the WGD model absorbs the burst (sigma_W -> ~sigma_N), and the
    LRT is large;
  - --synteny-min-sigmaW 0.3: the LRT drops by a large factor and q-hat becomes
    small. It need not reach exactly 0: with a small q the WGD can still be
    placed selectively on the few nodes whose SSD pairs happen to be
    collinear. Calibrate the remaining statistic with scripts/synteny_bootstrap.py.
A second data set with a real WGD (q = 0.3) checks that the bound does not
change the fit when sigma_W is identified (sigma_W-hat ~ 0.7 > bound).
usage: synteny_bound_regression.py <kalerax> <workdir> [families]
"""
import os, re, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from kalerax_sim import SpeciesTree, Params, simulate_families, write_family_set, synteny_lines

SP = ("(Passiflora_edulis:0.158589,(Casearia_decandra:0.138849,(Dianyuea_turbinata:0.145867,"
      "((Itoa_orientalis:0.0579734,(Idesia_polycarpa:0.0425043,(Populus_alba:0.0320969,"
      "Salix_arbutifolia:0.0484897):0.0206582):0.00545601):0.0175212,(Dovyalis_caffra:0.0489127,"
      "((Flacourtia_jangomas:0.0231842,Xylosma_longifolia:0.0194545):0.00992546,"
      "Scolopia_chinensis:0.0296868):0.0162224):0.015656):0.031979):0.0275982):0.158589);")
WGD = "Itoa_orientalis,Salix_arbutifolia"


def dataset(wd, n, seed, burst):
    t = SpeciesTree(SP); e = t.lca(WGD.split(","))
    p = Params(t, 0.11, 0.075, 0.043, wgds=[] if burst else [(e, 0.3, 1.0)],
               branch_rates={e: (0.333, 0.075, 0.043)} if burst else None)
    fams, _ = simulate_families(t, p, n, seed, origin_strategy="UNIFORM", keep_nodes=True)
    ff = write_family_set(fams, wd)
    os.makedirs(os.path.join(wd, "syn"), exist_ok=True)
    rng = np.random.default_rng(seed)
    for i, f in enumerate(fams):
        L = synteny_lines(f["root"], rng, 0.7, 0.03, 0.1)
        if L:
            open(os.path.join(wd, "syn", f"fam{i}.syn"), "w").write("\n".join(L) + "\n")
    return ff


def run(kx, ff, sp, out, extra):
    cmd = [kx, "-f", ff, "-s", sp, "--rec-model", "UndatedDTL", "--transfer-constraint", "PARENTS",
           "--gene-tree-rooting", "UNIFORM", "--species-tree-search", "SKIP", "--summary-only",
           "--gene-tree-samples", "0", "--seed", "1", "--synteny-dir",
           os.path.join(os.path.dirname(ff), "syn"), "-p", out] + extra
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    txt = open(os.path.join(out, "alerax.log")).read()
    ll = float(re.findall(r"Final species tree likelihood: ll=(-?[0-9.eE+-]+)", txt)[-1])
    q = re.findall(r"q\[node \d+\]=([0-9.eE+-]+)", txt)
    s = re.findall(r"sigmaW=([0-9.eE+-]+)", txt)
    return ll, (float(q[-1]) if q else 0.0), (float(s[-1]) if s else None)


def main():
    kx, wd = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 1200
    os.makedirs(wd, exist_ok=True)
    sp = os.path.join(wd, "sp.nw"); open(sp, "w").write(SP + "\n")
    ok = True
    # (1) SSD burst, no WGD
    ff = dataset(os.path.join(wd, "burst"), n, 11, True)
    l0, _, _ = run(kx, ff, sp, os.path.join(wd, "b_null"), [])
    lf, qf, sf = run(kx, ff, sp, os.path.join(wd, "b_free"), ["--wgd", WGD])
    lb, qb, sb = run(kx, ff, sp, os.path.join(wd, "b_bound"), ["--wgd", WGD, "--synteny-min-sigmaW", "0.3"])
    print(f"  SSD burst: free sigma_W 2dlnL={2*(lf-l0):.1f} q={qf:.3f} sigmaW={sf}; "
          f"bounded 2dlnL={2*(lb-l0):.1f} q={qb:.3f} sigmaW={sb}")
    if not (2 * (lf - l0) > 20 and 2 * (lb - l0) < 0.15 * 2 * (lf - l0) and qb < 0.06):
        print("  FAIL: bound does not remove the SSD-burst WGD"); ok = False
    # (2) real WGD
    ff = dataset(os.path.join(wd, "wgd"), n, 12, False)
    lf, qf, sf = run(kx, ff, sp, os.path.join(wd, "w_free"), ["--wgd", WGD])
    lb, qb, sb = run(kx, ff, sp, os.path.join(wd, "w_bound"), ["--wgd", WGD, "--synteny-min-sigmaW", "0.3"])
    print(f"  real WGD: free q={qf:.3f} sigmaW={sf}; bounded q={qb:.3f} sigmaW={sb}; dlnL={lb-lf:.2f}")
    if not (abs(lb - lf) < 1.0 and sb is not None and sb > 0.5 and qb > 0.2):
        print("  FAIL: bound changed a fit where sigma_W is identified"); ok = False
    print("PASS: --synteny-min-sigmaW" if ok else "FAIL: --synteny-min-sigmaW")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
