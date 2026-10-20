#!/usr/bin/env python3
"""Regression gate for --lore-fixed-r / --fix-wgd-q (fixed-parameter evaluation).

Simulates (scripts/kalerax_sim.py) an 11-taxon data set with one internal LORe
WGD (q = 0.4, r = 0.4) and emission-level synteny, then:
  (1) free fit: --wgd W --lore (+ synteny) -> full-precision parameters from
      model_parameters/wgd_parameters.txt;
  (2) fully fixed re-evaluation at those values (--fix-rates --d/--l/--t
      --wgd W:q --fix-wgd-q --lore --lore-fixed-r r --fix-synteny
      --synteny-params sW,sN --skip-thorough-rates) must reproduce the free-fit
      lnL (sum of per_fam_likelihoods.txt) to < 0.05;
  (3) r fixed at r-hat only (q, rates, sigmas re-fitted) must reproduce the
      free lnL to < 0.5 and the free q to < 0.02;
  (4) fully fixed at other r on a grid must not beat r-hat;
  (5) --lore-fixed-r 1 must equal the AORe model (no --lore) at the same q;
  (6) argument checks: --lore-fixed-r without --lore and r outside [0,1] fail.
usage: lore_fixed_r_regression.py <kalerax> <workdir> [families]
"""
import os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from kalerax_sim import SpeciesTree, Params, simulate_families, write_family_set, synteny_lines

SP = ("(Passiflora_edulis:0.158589,(Casearia_decandra:0.138849,(Dianyuea_turbinata:0.145867,"
      "((Itoa_orientalis:0.0579734,(Idesia_polycarpa:0.0425043,(Populus_alba:0.0320969,"
      "Salix_arbutifolia:0.0484897):0.0206582):0.00545601):0.0175212,(Dovyalis_caffra:0.0489127,"
      "((Flacourtia_jangomas:0.0231842,Xylosma_longifolia:0.0194545):0.00992546,"
      "Scolopia_chinensis:0.0296868):0.0162224):0.015656):0.031979):0.0275982):0.158589);")
WGD = "Itoa_orientalis,Salix_arbutifolia"


def dataset(wd, n, seed):
    t = SpeciesTree(SP); e = t.lca(WGD.split(","))
    p = Params(t, 0.11, 0.075, 0.043, wgds=[(e, 0.4, 0.4)])
    fams, _ = simulate_families(t, p, n, seed, origin_strategy="UNIFORM", keep_nodes=True)
    ff = write_family_set(fams, wd)
    os.makedirs(os.path.join(wd, "syn"), exist_ok=True)
    rng = np.random.default_rng(seed)
    for i, f in enumerate(fams):
        L = synteny_lines(f["root"], rng, 0.7, 0.03, 0.1)
        if L:
            open(os.path.join(wd, "syn", f"fam{i}.syn"), "w").write("\n".join(L) + "\n")
    return ff


def run(kx, ff, sp, out, extra, check=True):
    cmd = [kx, "-f", ff, "-s", sp, "--rec-model", "UndatedDTL", "--transfer-constraint", "PARENTS",
           "--gene-tree-rooting", "UNIFORM", "--species-tree-search", "SKIP", "--summary-only",
           "--gene-tree-samples", "0", "--seed", "1", "--synteny-dir",
           os.path.join(os.path.dirname(ff), "syn"), "-p", out] + extra
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not check:
        return r.returncode
    assert r.returncode == 0, f"kalerax failed: {' '.join(extra)}"
    ll = sum(float(l.split()[1]) for l in open(os.path.join(out, "per_fam_likelihoods.txt")) if l.strip())
    par = {}
    for l in open(os.path.join(out, "model_parameters", "wgd_parameters.txt")):
        if l.startswith("#"):
            continue
        k, lab, _, v1, v2 = l.rstrip("\n").split("\t")
        par[(k, lab)] = (float(v1), None if v2 == "NA" else float(v2))
    return ll, par


def fixed_args(par, r=None, lore=True):
    q = [v[0] for (k, _), v in par.items() if k == "wgd"][0]
    rate = {lab: v[0] for (k, lab), v in par.items() if k == "rate"}
    sw, sn = par[("synteny", "sigma")]
    a = ["--fix-rates", "--d", repr(rate["D"]), "--l", repr(rate["L"]), "--t", repr(rate["T"]),
         "--wgd", f"{WGD}:{q!r}", "--fix-wgd-q", "--fix-synteny", "--synteny-params", f"{sw!r},{sn!r}",
         "--skip-thorough-rates"]
    if lore:
        a += ["--lore", "--lore-fixed-r", repr(r)]
    return a


def main():
    kx, wd = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    os.makedirs(wd, exist_ok=True)
    sp = os.path.join(wd, "sp.nw"); open(sp, "w").write(SP + "\n")
    ff = dataset(os.path.join(wd, "data"), n, 21)
    ok = True
    # (1) free fit
    lf, pf = run(kx, ff, sp, os.path.join(wd, "free"), ["--wgd", WGD, "--lore"])
    qf, rf = [v for (k, _), v in pf.items() if k == "wgd"][0]
    print(f"  free fit: lnL={lf:.4f} q={qf:.4f} r={rf:.4f} sigma={pf[('synteny', 'sigma')]}")
    # (2) fully fixed at the free estimate
    lx, _ = run(kx, ff, sp, os.path.join(wd, "fixed"), fixed_args(pf, rf))
    print(f"  fully fixed at r-hat: lnL={lx:.4f} (diff {lx - lf:+.2e})")
    if abs(lx - lf) > 0.05:
        print("  FAIL: fixed-parameter evaluation does not reproduce the free fit"); ok = False
    # (3) r fixed only
    lr, pr = run(kx, ff, sp, os.path.join(wd, "rfixed"), ["--wgd", WGD, "--lore", "--lore-fixed-r", repr(rf)])
    qr, rr = [v for (k, _), v in pr.items() if k == "wgd"][0]
    print(f"  r fixed at r-hat, rest fitted: lnL={lr:.4f} (diff {lr - lf:+.3f}) q={qr:.4f} r={rr:.4f}")
    if abs(lr - lf) > 0.5 or abs(qr - qf) > 0.02 or abs(rr - rf) > 1e-12:
        print("  FAIL: --lore-fixed-r alone does not reproduce the free fit"); ok = False
    # (4) grid
    for r in (0.0, 0.25, 0.5, 0.75, 1.0):
        lg, pg = run(kx, ff, sp, os.path.join(wd, f"grid_{r}"), fixed_args(pf, r))
        rg = [v for (k, _), v in pg.items() if k == "wgd"][0][1]
        print(f"  grid r={r}: lnL={lg:.4f} (vs r-hat {lg - lx:+.3f}); reported r={rg}")
        if lg > lx + 1e-6 or abs(rg - r) > 1e-12:
            print("  FAIL: grid point beats r-hat or r was not held"); ok = False
        if r == 1.0:
            l1 = lg
    # (5) r = 1 == AORe at the same q
    la, _ = run(kx, ff, sp, os.path.join(wd, "aore"), fixed_args(pf, lore=False))
    print(f"  AORe fixed: lnL={la:.6f}; LORe r=1 fixed: {l1:.6f} (diff {l1 - la:+.2e})")
    if abs(l1 - la) > 1e-6 * abs(la):
        print("  FAIL: --lore-fixed-r 1 differs from AORe"); ok = False
    # (6) argument checks
    bad1 = run(kx, ff, sp, os.path.join(wd, "bad1"), ["--wgd", WGD, "--lore-fixed-r", "0.5"], check=False)
    bad2 = run(kx, ff, sp, os.path.join(wd, "bad2"), ["--wgd", WGD, "--lore", "--lore-fixed-r", "1.5"], check=False)
    if bad1 == 0 or bad2 == 0:
        print(f"  FAIL: invalid --lore-fixed-r accepted (exit {bad1}, {bad2})"); ok = False
    print("PASS: --lore-fixed-r / --fix-wgd-q" if ok else "FAIL: --lore-fixed-r / --fix-wgd-q")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
