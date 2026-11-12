#!/usr/bin/env bash
# Regression: clade values above 1 are valid under the forward-process
# duplication weights (factor 2 per D / WGD doubling split) and must not be
# rejected.
#
# A terminal WGD with retention q ~ 1 and negligible loss makes a retained
# cherry score ~2*R*R ~ 2: the factor 2 is the labelling factor of the two
# exchangeable copies (cladeProbaBound in the undated models). b959574
# rejected such evaluations (lnL -1e15). For UndatedDTL and UndatedDL:
#   1. q=0.999: finite lnL, no "proba >" / "lh >" errors;
#   2. lnL(q=0.999) - lnL(q=0.5) == ln(0.999/0.5) for the one-split family
#      (with D, L, T ~ 0 the retained cherry contributes q * 2 R R);
#   3. lnL(q=0.999) matches the value validated in b67fd23 / 0a024d6
#      (= the old ALE-weighted value + ln 2).
# Usage: tests/forward_weight_bound_regression.sh <path/to/kalerax>
set -euo pipefail
KX=${1:?usage: $0 <kalerax binary>}
KX=$(cd "$(dirname "$KX")" && pwd)/$(basename "$KX")
WD=$(mktemp -d)
trap 'rm -rf "$WD"' EXIT
cd "$WD"
echo "((A:1,B:1)AB:1,(C:1,D:1)CD:1)R;" > sp.nwk
echo "(((a1,a2),b1),(c1,d1));" > fam1.trees
printf "A:a1;a2\nB:b1\nC:c1\nD:d1\n" > fam1.link
printf "[FAMILIES]\n- fam1\ngene_tree = %s/fam1.trees\nmapping = %s/fam1.link\n" "$WD" "$WD" > fams.txt
fail=0
for model in UndatedDTL UndatedDL; do
  t=""; [ $model = UndatedDTL ] && t="--t 1e-6"
  for q in 0.5 0.999; do
    "$KX" --families fams.txt --species-tree sp.nwk --species-tree-search SKIP \
      --gene-tree-samples 0 --rec-model $model --fix-rates --d 1e-6 --l 1e-6 $t \
      --skip-thorough-rates --wgd A:$q --fix-wgd-q --prefix "o_${model}_$q" \
      > "o_${model}_$q.out" 2>&1
    if grep -q -e 'proba >' -e 'lh >' "o_${model}_$q.out"; then
      echo "FAIL $model q=$q: bound error logged"; fail=1
    fi
  done
  python3 - "$model" "$(awk '{print $2}' o_${model}_0.5/per_fam_likelihoods.txt)" \
    "$(awk '{print $2}' o_${model}_0.999/per_fam_likelihoods.txt)" <<'PY' || fail=1
import sys, math
model, a, b = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
ref = {"UndatedDTL": -3.19969546682, "UndatedDL": -3.19968861685}[model]
ok = b > -1e10 and abs((b - a) - math.log(0.999 / 0.5)) < 1e-4 and abs(b - ref) < 1e-6
print(f"{model}: lnL q=0.5 {a:.6f}, q=0.999 {b:.6f}; diff {b-a:.6f} "
      f"(expected {math.log(0.999/0.5):.6f}); reference {ref:.6f}", "OK" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
done
if [ $fail -eq 0 ]; then echo "PASS forward_weight_bound_regression"; else exit 1; fi
