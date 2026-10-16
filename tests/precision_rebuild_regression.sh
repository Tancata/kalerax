#!/usr/bin/env bash
# Regression gate: a family rebuilt in high precision (double underflow,
# lnL < ~ -745) must keep its WGD. Before the fix, AleEvaluator::resetEvaluation
# re-created the model without re-applying q / r, so after resetAllPrecisions()
# the large family was silently evaluated (and exported) WITHOUT the WGD.
# usage: precision_rebuild_regression.sh [path/to/kalerax]
set -euo pipefail
BIN="${1:-$(dirname "$0")/../build/bin/kalerax}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT; export OMP_NUM_THREADS=1
W=$TMP; echo '(OG:2,((A1:1,A2:1):1,(B1:1,B2:1):1):1);' > $W/sp.nw
python3 - $W <<'E'
import random, sys
random.seed(5); W = sys.argv[1]; subs = []; leaves = []
for i in range(130):  # 130 WGD-structured subtrees -> lnL ~ -900 (underflows double)
    subs.append(f"(OG_{i},((A1_{i},A2_{i}),((B1_{i}a,B2_{i}a),(B1_{i}b,B2_{i}b))))")
    leaves += [f"OG_{i}", f"A1_{i}", f"A2_{i}", f"B1_{i}a", f"B2_{i}a", f"B1_{i}b", f"B2_{i}b"]
nodes = subs[:]
while len(nodes) > 1:
    a = nodes.pop(random.randrange(len(nodes))); b = nodes.pop(random.randrange(len(nodes))); nodes.append(f"({a},{b})")
open(f"{W}/big.nw", "w").write(nodes[0] + ";\n")
open(f"{W}/big.link", "w").write("".join(f"{l}\t{l.split('_')[0]}\n" for l in leaves))
open(f"{W}/s.nw", "w").write("(OG_1,((A1_1,A2_1),((B1_1,B2_1),(B1_2,B2_2))));\n")
open(f"{W}/s.link", "w").write("A1_1\tA1\nA2_1\tA2\nB1_1\tB1\nB1_2\tB1\nB2_1\tB2\nB2_2\tB2\nOG_1\tOG\n")
with open(f"{W}/fam.txt", "w") as f:
    f.write(f"[FAMILIES]\n- big\ngene_tree = {W}/big.nw\nmapping = {W}/big.link\n")
    for i in range(10): f.write(f"- s{i}\ngene_tree = {W}/s.nw\nmapping = {W}/s.link\n")
E
C=(-s $W/sp.nw -f $W/fam.txt --gene-tree-rooting UNIFORM --species-tree-search SKIP --seed 1 -g 3 --rec-model UndatedDL --fix-rates)
"$BIN" "${C[@]}" -p $W/nowgd > $W/nowgd.log 2>&1
"$BIN" "${C[@]}" --wgd B1,B2 -p $W/wgd > $W/wgd.log 2>&1
L0=$(awk '$1=="big"{print $2}' $W/nowgd/per_fam_likelihoods.txt)
L1=$(awk '$1=="big"{print $2}' $W/wgd/per_fam_likelihoods.txt)
LQ=$(grep -E "After WGD retention opt" $W/wgd.log | tail -1 | sed -E 's/.*ll=([-0-9.e+]+).*/\1/')
LF=$(grep -E "Final species tree likelihood" $W/wgd.log | tail -1 | sed -E 's/.*ll=//')
echo "  big family lnL: no WGD=$L0  +WGD=$L1 | ll after q-opt=$LQ  final ll=$LF"
awk -v l0=$L0 -v l1=$L1 -v lq=$LQ -v lf=$LF 'BEGIN{
  if (l1 - l0 < 50) { print "FAIL: high-precision family lost its WGD (lnL == no-WGD value)"; exit 1 }
  if (lf - lq < -1 || lf - lq > 1) { print "FAIL: final lnL inconsistent with the fitted model"; exit 1 }
  print "PASS: precision rebuild keeps the WGD (q, r, synteny re-applied)"
}'
