"""Independent check of kalerax DTL duplication weighting with Zombi (Davin et al. 2020).

For each config: one Zombi species tree (pure birth, 12 extant species), K
independent Zombi Gm runs of M root families each (additive transfers to
contemporaneous lineages, no inversion/transposition/origination). True
observable D/T = events at internal nodes of the pruned gene trees. kalerax
(UndatedDTL, PARENTS, --origination ROOT) is fitted with default and
--ordered-dup-weights; exported D/T totals (-g 10, --summary-only) are compared
with the truth over the families kalerax kept.

usage: ZOMBI_DIR=... KALERAX=... python zombi_check.py <outdir> <tag> D T L [K M]
Results (October 2026): validation/zombi_dup_weights.md
"""
import os, re, shutil, subprocess, sys, json
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ZOMBI = os.environ.get("ZOMBI_DIR", os.path.join(HERE, "Zombi"))
PY = sys.executable
KALERAX = os.environ.get("KALERAX", os.path.join(HERE, "..", "build", "bin", "kalerax"))


def params(src, dst, subs):
    s = open(src).read()
    for a, b in subs:
        assert a in s, a
        s = s.replace(a, b)
    open(dst, "w").write(s)


def run_zombi(args, cwd_out):
    env = dict(os.environ, PYTHONPATH=ZOMBI)
    subprocess.run([PY, os.path.join(ZOMBI, "Zombi.py")] + args, env=env, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def truth_and_files(gdir, fams_dir, prefix):
    """Write <fam>.nw/.link for every family with >=1 extant gene; return
    {family: (nD, nT)} of observable events."""
    out = {}
    gt = os.path.join(gdir, "Gene_trees"); ev = os.path.join(gdir, "Gene_families")
    for f in os.listdir(gt):
        if not f.endswith("_prunedtree.nwk"):
            continue
        fam = f.split("_")[0]
        nwk = open(os.path.join(gt, f)).read().strip()
        if nwk in (";", "") or "(" not in nwk:
            continue
        evt = {}
        for l in open(os.path.join(ev, f"{fam}_events.tsv")).readlines()[1:]:
            t, e, nodes = l.rstrip("\n").split("\t")
            if e in ("D", "T", "S"):
                p = nodes.split(";")
                evt[p[1]] = e  # parent gene id -> event that split it
        internal = re.findall(r"\)([A-Za-z0-9]+)_(\d+)", nwk)
        nD = sum(1 for _, g in internal if evt.get(g) == "D")
        nT = sum(1 for _, g in internal if evt.get(g) == "T")
        leaves = re.findall(r"[(,]([A-Za-z0-9]+)_(\d+):", nwk)
        if len({sp for sp, _ in leaves}) < 4:   # kalerax drops these; CCPs need >=3 leaves
            continue
        name = f"{prefix}f{fam}"
        clean = re.sub(r"\)[A-Za-z0-9]+_\d+", ")", nwk)  # drop internal labels
        clean = re.sub(r"\)Root_\d+", ")", clean)
        open(os.path.join(fams_dir, name + ".nw"), "w").write(clean + "\n")
        open(os.path.join(fams_dir, name + ".link"), "w").write(
            "".join(f"{sp}_{g}\t{sp}\n" for sp, g in leaves))
        out[name] = (nD, nT, len(leaves))
    return out


def main():
    outdir, tag = os.path.abspath(sys.argv[1]), sys.argv[2]
    D, T, L = sys.argv[3:6]
    K = int(sys.argv[6]) if len(sys.argv) > 6 else 14
    M = int(sys.argv[7]) if len(sys.argv) > 7 else 300
    wd = os.path.join(outdir, tag); shutil.rmtree(wd, ignore_errors=True); os.makedirs(wd)
    # species tree
    params(os.path.join(ZOMBI, "Parameters/SpeciesTreeParameters.tsv"), os.path.join(wd, "T.tsv"),
           [("\nEXTINCTION\tf:0.2", "\nEXTINCTION\tf:0"), ("\nTOTAL_LINEAGES\t20", "\nTOTAL_LINEAGES\t12"),
            ("\nSEED 0", "\nSEED 11")])
    base = os.path.join(wd, "base"); run_zombi(["T", os.path.join(wd, "T.tsv"), base], wd)
    sp = open(os.path.join(base, "T", "ExtantTree.nwk")).read().strip()
    sp = re.sub(r"\)Root:[0-9.eE-]+;", ");", sp)
    open(os.path.join(wd, "species.nw"), "w").write(sp + "\n")
    fams_dir = os.path.join(wd, "fams"); os.makedirs(fams_dir)

    def one(k):
        rd = os.path.join(wd, f"z{k}"); shutil.copytree(os.path.join(base, "T"), os.path.join(rd, "T"))
        os.makedirs(os.path.join(rd, "G"), exist_ok=True)
        gp = os.path.join(rd, "G.tsv")
        params(os.path.join(ZOMBI, "Parameters/GenomeParameters.tsv"), gp,
               [("DUPLICATION f:1", f"DUPLICATION f:{D}"), ("TRANSFER f:2", f"TRANSFER f:{T}"),
                ("LOSS f:3", f"LOSS f:{L}"), ("INVERSION f:2", "INVERSION f:0"),
                ("TRANSPOSITION f:2", "TRANSPOSITION f:0"), ("ORIGINATION f:2", "ORIGINATION f:0"),
                ("REPLACEMENT_TRANSFER 1", "REPLACEMENT_TRANSFER 0"),
                ("INITIAL_GENOME_SIZE\t100", f"INITIAL_GENOME_SIZE\t{M}"),
                ("MIN_GENOME_SIZE\t10", "MIN_GENOME_SIZE\t0"), ("PROFILES 1", "PROFILES 0"),
                ("EVENTS_PER_BRANCH 1", "EVENTS_PER_BRANCH 0"), ("VERBOSE 1", "VERBOSE 0"),
                ("SEED 0", f"SEED {1000 + k}")])
        run_zombi(["Gm", gp, rd], rd)
        tr = truth_and_files(os.path.join(rd, "G"), fams_dir, f"z{k}")
        shutil.rmtree(rd, ignore_errors=True)
        return tr

    truth = {}
    with ThreadPoolExecutor(max_workers=14) as ex:
        for tr in ex.map(one, range(K)):
            truth.update(tr)
    lines = ["[FAMILIES]"]
    for name in sorted(truth):
        lines += [f"- {name}", f"gene_tree = {fams_dir}/{name}.nw", f"mapping = {fams_dir}/{name}.link"]
    ff = os.path.join(wd, "families.txt"); open(ff, "w").write("\n".join(lines) + "\n")

    def fit(ordered):
        o = os.path.join(wd, f"fit_{int(ordered)}")
        cmd = [KALERAX, "-f", ff, "-s", os.path.join(wd, "species.nw"), "--rec-model", "UndatedDTL",
               "--transfer-constraint", "PARENTS", "--gene-tree-rooting", "UNIFORM", "--origination", "ROOT",
               "--species-tree-search", "SKIP", "--summary-only", "--gene-tree-samples", "10", "--seed", "1", "-p", o]
        if not ordered:
            cmd.append("--ale-dup-weights")   # ALE/AleRax convention (new default = forward process)
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        kept = [l.split()[0] for l in open(os.path.join(o, "per_fam_likelihoods.txt")) if l.strip()]
        rows = [l.replace(",", " ").split() for l in open(os.path.join(o, "reconciliations", "totalSpeciesEventCounts.txt"))]
        hdr = rows[0]; iD, iT = hdr.index("duplications") , hdr.index("transfers")
        expD = sum(float(r[iD]) for r in rows[1:]); expT = sum(float(r[iT]) for r in rows[1:])
        mp = open(os.path.join(o, "model_parameters", "model_parameters.txt")).readlines()[1].split()[1:4]
        ll = re.findall(r"Final species tree likelihood: ll=(-?[0-9.]+)", open(os.path.join(o, "alerax.log")).read())[-1]
        res = dict(ordered=ordered, kept=len(kept), exported_D=expD, exported_T=expT,
                   true_D=sum(truth[k][0] for k in kept), true_T=sum(truth[k][1] for k in kept),
                   D_hat=float(mp[0]), L_hat=float(mp[1]), T_hat=float(mp[2]), lnL=float(ll))
        shutil.rmtree(o, ignore_errors=True)
        return res

    with ThreadPoolExecutor(max_workers=2) as ex:
        res = list(ex.map(fit, [False, True]))
    for r in res:
        r.update(tag=tag, D=float(D), T=float(T), L=float(L),
                 D_err_pct=100 * (r["exported_D"] - r["true_D"]) / r["true_D"],
                 T_err_pct=100 * (r["exported_T"] - r["true_T"]) / r["true_T"])
        print(json.dumps(r))
    json.dump(res, open(os.path.join(outdir, f"{tag}.json"), "w"), indent=1)
    shutil.rmtree(fams_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
