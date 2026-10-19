#!/bin/bash
# Regression test: synteny split weights far below the double range.
# A family with many collinear pairs whose LCA cannot be a declared WGD split
# (here 320 same-genome collinear pairs across a root duplication, WGD declared
# on another branch)
# has a synteny log-weight ~ 320*log(sigmaN/sigmaW) on that split. Before the
# fix this weight was stored as a double, underflowed to 0, and the family got
# ll=-inf even in high precision (seen on the Glycine dataset).
# Checks:
#  1. with scaling disabled (old behaviour) the run fails (test is sensitive)
#  2. null model: lnL(synteny) - lnL(no synteny) == 320*log(sigmaN) +
#     280*log(1-sigmaN) exactly (same-genome pairs only: 5 species x 8 x 8
#     listed collinear pairs; the 280 unlisted within-copy paralog pairs
#     default to y=0)
#  3. WGD declared: finite lnL, identical when the scaled path is forced on
#     every split (KALERAX_SYN_SCALE_THRESHOLD=-1e-9)
#  4. a large random family (184 genes, 30 trees, ~5% collinear pairs, UNIFORM
#     rooting, DTL; high-precision mode with competing CCP splits): lnL does
#     not depend on the scaling threshold. Before ScaledValue terms were fully
#     renormalised, unnormalised synteny-weighted terms won additions against
#     larger terms; kalerax 97ec8fe lost 3-7 log units here.
# usage: synteny_underflow_regression.sh [path/to/kalerax]
set -euo pipefail
BIN="${1:-$(dirname "$0")/../build/bin/kalerax}"
[ -x "$BIN" ] || { echo "FAIL: binary not found: $BIN"; exit 1; }
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
fail() { echo "FAIL: $*"; exit 1; }
d=$TMP; mkdir -p $d/g $d/m $d/syn
echo '(OG:2,((A1:1,A2:1):1,(B1:1,B2:1):1):1);' > $d/sp.nw
python3 - $d <<'E'
import sys
d = sys.argv[1]; K = 8; sp = ["OG", "A1", "A2", "B1", "B2"]
def cat(s, c):  # K recent paralogs of species s in copy c (caterpillar)
    g = [f"{s}_{c}{i}" for i in range(K)]
    t = g[0]
    for x in g[1:]: t = f"({t},{x})"
    return t, g
def copy(c):
    t = {s: cat(s, c) for s in sp}
    nw = f"({t['OG'][0]},(({t['A1'][0]},{t['A2'][0]}),({t['B1'][0]},{t['B2'][0]})))"
    return nw, [x for s in sp for x in t[s][1]]
(n1, g1), (n2, g2) = copy("a"), copy("b")
fams = ["[FAMILIES]"]
def add(name, nw, genes, pairs):
    open(f"{d}/g/{name}.nw", "w").write(nw + ";\n")
    open(f"{d}/m/{name}.link", "w").write("".join(f"{x}\t{x.split('_')[0]}\n" for x in genes))
    open(f"{d}/syn/{name}.syn", "w").write("".join(f"{a} {b} 1\n" for a, b in pairs))
    fams.append(f"- {name}\ngene_tree = {d}/g/{name}.nw\nmapping = {d}/m/{name}.link")
add("big", f"({n1},{n2})", g1 + g2, [(a, b) for a in g1 for b in g2])
for i in range(6):
    add(f"s{i}", "(OG_1,((A1_1,A2_1),(B1_1,B2_1)))", [f"{s}_1" for s in sp], [])
open(f"{d}/families.txt", "w").write("\n".join(fams) + "\n")
E
C=(-f $d/families.txt -s $d/sp.nw --gene-tree-rooting ROOTED --species-tree-search SKIP
   --seed 1 --summary-only --gene-tree-samples 0 --rec-model UndatedDL --fix-rates)
S=(--synteny-dir $d/syn --synteny-params 0.95,0.05 --fix-synteny)
lnl() { awk '$1=="big"{print $2}' $1/per_fam_likelihoods.txt; }
echo "binary: $BIN"
if KALERAX_SYN_SCALE_THRESHOLD=-1e300 "$BIN" "${C[@]}" "${S[@]}" -p $d/old > $d/old.log 2>&1 \
   && ! grep -q "ll=-inf" $d/old.log; then
  fail "with scaling disabled the family should underflow (test not sensitive)"
fi
echo "  PASS: unscaled weights underflow (old behaviour reproduced)"
"$BIN" "${C[@]}" -p $d/nosyn > $d/nosyn.log 2>&1 || fail "no-synteny run"
"$BIN" "${C[@]}" "${S[@]}" -p $d/syn_null > $d/syn_null.log 2>&1 || fail "null synteny run"
python3 - $(lnl $d/nosyn) $(lnl $d/syn_null) <<'E' || fail "null synteny offset"
import sys, math
a, b = map(float, sys.argv[1:3]); ex = 320 * math.log(0.05) + 280 * math.log(0.95)
print(f"  null: lnL(syn)-lnL(nosyn) = {b-a:.6f}, expected {ex:.6f}")
sys.exit(0 if abs(b - a - ex) < 0.01 else 1)  # 6 significant digits in per_fam_likelihoods
E
echo "  PASS: null synteny term exact"
"$BIN" "${C[@]}" "${S[@]}" --wgd B1,B2:0.5 -p $d/wgd > $d/wgd.log 2>&1 || fail "WGD synteny run"
KALERAX_SYN_SCALE_THRESHOLD=-1e-9 "$BIN" "${C[@]}" "${S[@]}" --wgd B1,B2:0.5 -p $d/wgdf > $d/wgdf.log 2>&1 \
  || fail "WGD synteny run, forced scaling"
python3 - $(lnl $d/wgd) $(lnl $d/wgdf) <<'E' || fail "forced scaling changed lnL"
import sys, math
a, b = map(float, sys.argv[1:3]); print(f"  WGD: lnL={a:.6f} forced-scaling lnL={b:.6f}")
sys.exit(0 if math.isfinite(a) and abs(a - b) < 1e-6 * abs(a) else 1)
E
# 4. threshold invariance on a large random family
r=$d/rnd; mkdir -p $r/g $r/m $r/syn
python3 - $r <<'E'
import sys, random
d = sys.argv[1]; random.seed(3)
cnt = {"OG": 6, "A1": 8, "A2": 10, "B1": 120, "B2": 40}
genes = [f"{s}_{i}" for s, n in cnt.items() for i in range(n)]
def rt():
    nodes = genes[:]
    while len(nodes) > 1:
        i = random.randrange(len(nodes))
        j = min(len(nodes) - 1, max(0, i + random.choice([-2, -1, 1, 2])))
        if i == j: continue
        a, b = nodes[i], nodes[j]; k = min(i, j)
        nodes = [x for t, x in enumerate(nodes) if t not in (i, j)]
        nodes.insert(k, f"({a},{b})")
    return nodes[0]
open(f"{d}/g/big.nw", "w").write("".join(rt() + ";\n" for _ in range(30)))
open(f"{d}/m/big.link", "w").write("".join(f"{g}\t{g.split('_')[0]}\n" for g in genes))
L = []
for s, n in cnt.items():
    for i in range(n):
        for j in range(i + 1, n):
            L.append(f"{s}_{i} {s}_{j} {1 if random.random() < 0.05 else 0}\n")
open(f"{d}/syn/big.syn", "w").write("".join(L))
open(f"{d}/families.txt", "w").write(
    f"[FAMILIES]\n- big\ngene_tree = {d}/g/big.nw\nmapping = {d}/m/big.link\n")
E
R=(-f $r/families.txt -s $d/sp.nw --gene-tree-rooting UNIFORM --species-tree-search SKIP
   --seed 1 --summary-only --gene-tree-samples 0 --rec-model UndatedDTL --fix-rates
   --synteny-dir $r/syn --synteny-params 0.8,0.0227 --fix-synteny --wgd B1,B2:0.5)
"$BIN" "${R[@]}" -p $r/def > $r/def.log 2>&1 || fail "random family, default threshold"
for t in -1e-9 -50; do
  KALERAX_SYN_SCALE_THRESHOLD=$t "$BIN" "${R[@]}" -p $r/t$t > $r/t$t.log 2>&1 \
    || fail "random family, threshold $t"
  python3 - $(lnl $r/def) $(lnl $r/t$t) $t <<'E' || fail "lnL depends on the scaling threshold"
import sys, math
a, b = map(float, sys.argv[1:3]); print(f"  random family: lnL default={a:.4f} threshold {sys.argv[3]}={b:.4f}")
sys.exit(0 if math.isfinite(a) and abs(a - b) < 1e-6 * abs(a) else 1)
E
done
echo "  PASS: lnL independent of the scaling threshold"
echo "PASS: synteny weights below the double range are exact"
