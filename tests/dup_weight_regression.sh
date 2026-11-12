#!/usr/bin/env bash
# Regression gate for the duplication-split weighting. kalerax always uses
# the forward-process weights (factor 2 on every D and WGD/LORe doubling
# split, UndatedDTL and UndatedDL; upstream AleRax uses the same factor for D
# since 8e34d87). The ALE option was removed after 0a024d6.
#  (1) --ale-dup-weights is rejected with an explanation; --ordered-dup-weights
#      is accepted (no-op);
#  (2) a new run writes the checkpoint marker "2 dl_and_dtl";
#  (3) checkpoints that used other weights are refused on resume: no marker
#      (pre-October 2026), "1" (--ale-dup-weights) and, for UndatedDL, a bare
#      "2" (DL before the AleRax v1.4.3 merge);
#  (4) the weights are concordant with forward simulation
#      (dtl_dup_weight_concordance.py, dl_dup_weight_concordance.py; python3 +
#      numpy + scipy; the DL test takes ~10 min at its default 8M families).
# usage: dup_weight_regression.sh [kalerax] [N DTL families] [N DL families]
set -euo pipefail
BIN="${1:-$(dirname "$0")/../build/bin/kalerax}"
BIN=$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")
N="${2:-600000}"
NDL="${3:-8000000}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
D="$(cd "$(dirname "$0")" && pwd)"
E="$D/../validation/example-1"
[ -f "$E/families.txt" ] || (cd "$E" && python3 make_inputs.py > /dev/null)
sed -e "s|^gene_tree = .*/|gene_tree = $E/|" -e "s|^mapping = .*/|mapping = $E/mappings/|" \
    "$E/families.txt" > "$TMP/fam.txt"
C=(-f "$TMP/fam.txt" -s "$E/tree.nw" --gene-tree-rooting UNIFORM
   --species-tree-search SKIP --gene-tree-samples 0 --seed 1 --wgd ATHA,ATRI)
# (1)
set +e
"$BIN" "${C[@]}" --summary-only --rec-model UndatedDTL --ale-dup-weights -p "$TMP/ale" > "$TMP/ale.log" 2>&1; rc=$?
"$BIN" "${C[@]}" --summary-only --rec-model UndatedDTL --ordered-dup-weights -p "$TMP/ord" > /dev/null 2>&1; rc2=$?
set -e
if [ $rc -ne 0 ] && grep -q "ale-dup-weights was removed" "$TMP/ale.log" && [ $rc2 -eq 0 ]; then
  echo "  PASS: --ale-dup-weights rejected with an explanation; --ordered-dup-weights accepted"
else
  echo "  FAIL: option handling (ale rc=$rc, ordered rc=$rc2)"; exit 1
fi
# (2) (no --summary-only: it deletes the CCPs, which a resume needs)
for m in UndatedDTL UndatedDL; do
  "$BIN" "${C[@]}" --rec-model $m -p "$TMP/$m" > /dev/null 2>&1
  mk=$(cat "$TMP/$m/checkpoint/dup_weighting.txt")
  [ "$mk" = "2 dl_and_dtl" ] || { echo "  FAIL: $m checkpoint marker '$mk'"; exit 1; }
done
echo "  PASS: checkpoint marker '2 dl_and_dtl' written (UndatedDTL, UndatedDL)"
# (3) rewind to step 5 (line 2; line 1 is the completeness header) and resume
# with an old marker
resume_refused() { # model, marker content ('' = no marker)
  local m=$1 mk=$2 d="$TMP/r_$1_${2// /_}"
  # a fresh run per case: a resume must repeat the original command (same -p)
  rm -rf "$d"; "$BIN" "${C[@]}" --rec-model $m -p "$d" > /dev/null 2>&1
  if [ -z "$mk" ]; then rm "$d/checkpoint/dup_weighting.txt"; else echo "$mk" > "$d/checkpoint/dup_weighting.txt"; fi
  python3 -c "import sys;p=sys.argv[1];L=open(p).read().split('\n');L[1]='5';open(p,'w').write('\n'.join(L))" "$d/checkpoint/mainCheckpoint.txt"
  set +e; "$BIN" "${C[@]}" --rec-model $m -p "$d" > "$d.log" 2>&1; local rc=$?; set -e
  [ $rc -ne 0 ] && grep -q "cannot resume this checkpoint" "$d.log"
}
for case in "UndatedDTL|" "UndatedDTL|1" "UndatedDL|" "UndatedDL|1" "UndatedDL|2"; do
  m=${case%%|*}; mk=${case#*|}
  if resume_refused "$m" "$mk"; then
    echo "  PASS: $m checkpoint with marker '${mk:-none}' refused on resume"
  else
    echo "  FAIL: $m checkpoint with marker '${mk:-none}' was not refused"; exit 1
  fi
done
# (4)
python3 "$D/dtl_dup_weight_concordance.py" "$BIN" "$TMP/conc_dtl" "$N"
python3 "$D/dl_dup_weight_concordance.py" "$BIN" "$TMP/conc_dl" "$NDL"
echo "PASS: dup_weight_regression"
