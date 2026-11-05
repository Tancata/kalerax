#!/usr/bin/env bash
# Regression: clade values above 1 are valid under forward-process duplication
# weighting (the default) and must not be rejected.
#
# A terminal WGD with retention q ~ 1 and negligible loss makes a retained
# cherry score ~2*R*R ~ 2: the factor 2 is the labelling factor of the two
# exchangeable copies (UndatedDTLMultiModel::cladeProbaBound). b959574 rejected
# such evaluations (lnL -1e15). Checks:
#   1. q=0.999, forward weighting: finite lnL, no "proba >" errors;
#   2. forward minus --ale-dup-weights lnL == ln 2 for the one-split family, at
#      q=0.5 and q=0.999 (the weighting factor, parameter-independent here);
#   3. a family with no duplication is identical under both weightings.
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
echo "((a1,b1),(c1,d1));" > fam2.trees
printf "A:a1\nB:b1\nC:c1\nD:d1\n" > fam2.link
printf "[FAMILIES]\n- fam1\ngene_tree = %s/fam1.trees\nmapping = %s/fam1.link\n- fam2\ngene_tree = %s/fam2.trees\nmapping = %s/fam2.link\n" \
  "$WD" "$WD" "$WD" "$WD" > fams.txt
run() { # tag, extra args...
  local tag=$1; shift
  "$KX" --families fams.txt --species-tree sp.nwk --species-tree-search SKIP \
    --gene-tree-samples 0 --rec-model UndatedDTL --fix-rates --d 1e-6 \
    --l 1e-6 --t 1e-6 --skip-thorough-rates --prefix "o_$tag" "$@" \
    > "o_$tag.out" 2>&1
}
ll() { awk -v f="$2" '$1==f{print $2}' "o_$1/per_fam_likelihoods.txt"; }
fail=0
for q in 0.5 0.999; do
  run "fw$q" --wgd A:$q --fix-wgd-q
  run "ale$q" --wgd A:$q --fix-wgd-q --ale-dup-weights
  for t in "fw$q" "ale$q"; do
    if grep -q 'proba >' "o_$t.out"; then echo "FAIL $t: proba bound error logged"; fail=1; fi
  done
  python3 - "$(ll fw$q fam1)" "$(ll ale$q fam1)" "$(ll fw$q fam2)" "$(ll ale$q fam2)" "$q" <<'PY' || fail=1
import sys, math
fw1, al1, fw2, al2 = map(float, sys.argv[1:5]); q = sys.argv[5]
ok = True
if fw1 < -1e10: print(f"FAIL q={q}: forward fam1 rejected ({fw1})"); ok = False
if abs((fw1 - al1) - math.log(2)) > 1e-6: print(f"FAIL q={q}: fw-ale = {fw1-al1}, expected ln2"); ok = False
if abs(fw2 - al2) > 1e-8: print(f"FAIL q={q}: no-dup family differs ({fw2} vs {al2})"); ok = False
print(f"q={q}: fam1 forward {fw1:.6f}, ALE {al1:.6f}, diff {fw1-al1:.6f} (ln2 {math.log(2):.6f}); fam2 {fw2:.6f} / {al2:.6f}", "OK" if ok else "")
sys.exit(0 if ok else 1)
PY
done
if [ $fail -eq 0 ]; then echo "PASS forward_weight_bound_regression"; else exit 1; fi
