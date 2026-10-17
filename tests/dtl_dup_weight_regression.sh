#!/usr/bin/env bash
# Regression gate for the duplication vs transfer split weighting.
#  (1) UndatedDL: --ale-dup-weights is a no-op (identical per-family lnL);
#  (2) labelling and checkpoint marker: a new UndatedDTL run records the
#      forward-process weighting; --ale-dup-weights records "ale"; a checkpoint
#      without a marker (created before the default switch) resumes with the
#      ALE weighting, so a run never mixes weightings;
#  (3) the default is concordant with forward simulation while
#      --ale-dup-weights is not (tests/dtl_dup_weight_concordance.py;
#      needs python3 + numpy).
# usage: dtl_dup_weight_regression.sh [path/to/kalerax] [N families per config]
set -euo pipefail
BIN="${1:-$(dirname "$0")/../build/bin/kalerax}"
N="${2:-600000}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
D="$(cd "$(dirname "$0")" && pwd)"
E="$D/../validation/example-1"
[ -f "$E/families.txt" ] || (cd "$E" && python3 make_inputs.py > /dev/null)
sed -e "s|^gene_tree = .*/|gene_tree = $E/|" -e "s|^mapping = .*/|mapping = $E/mappings/|" \
    "$E/families.txt" > "$TMP/fam.txt"
C=(-f "$TMP/fam.txt" -s "$E/tree.nw" --gene-tree-rooting UNIFORM
   --species-tree-search SKIP --summary-only --gene-tree-samples 0 --seed 1 --wgd ATHA,ATRI)
"$BIN" "${C[@]}" --rec-model UndatedDL -p "$TMP/dl0" > /dev/null 2>&1
"$BIN" "${C[@]}" --rec-model UndatedDL --ale-dup-weights -p "$TMP/dl1" > /dev/null 2>&1
if cmp -s <(sort "$TMP/dl0/per_fam_likelihoods.txt") <(sort "$TMP/dl1/per_fam_likelihoods.txt"); then
  echo "  PASS: UndatedDL per-family likelihoods identical with/without --ale-dup-weights"
else
  echo "  FAIL: UndatedDL changed by --ale-dup-weights"; exit 1
fi
# (no --summary-only for this run: it deletes the CCPs, which a resume needs)
CR=(-f "$TMP/fam.txt" -s "$E/tree.nw" --gene-tree-rooting UNIFORM
    --species-tree-search SKIP --gene-tree-samples 0 --seed 1 --wgd ATHA,ATRI)
"$BIN" "${CR[@]}" --rec-model UndatedDTL -p "$TMP/dtl" > /dev/null 2>&1
"$BIN" "${C[@]}" --rec-model UndatedDTL --ale-dup-weights -p "$TMP/dtla" > /dev/null 2>&1
w0=$(cut -d' ' -f1 "$TMP/dtl/model_parameters/dup_weighting.txt")
w1=$(cut -d' ' -f1 "$TMP/dtla/model_parameters/dup_weighting.txt")
m0=$(cat "$TMP/dtl/checkpoint/dup_weighting.txt")
if [ "$w0" = forward_process ] && [ "$w1" = ale ] && [ "$m0" = 2 ]; then
  echo "  PASS: weighting labelled (default=forward_process, --ale-dup-weights=ale) and checkpointed"
else
  echo "  FAIL: weighting label/marker ($w0, $w1, $m0)"; exit 1
fi
# simulate an interrupted pre-switch run: drop the marker and rewind the
# checkpoint to the final rate optimisation (step 5), then resume with the
# same command
rm "$TMP/dtl/checkpoint/dup_weighting.txt"
python3 - "$TMP/dtl/checkpoint/mainCheckpoint.txt" <<'PY'
import sys
p = sys.argv[1]; L = open(p).read().split("\n"); L[0] = "5"; open(p, "w").write("\n".join(L))
PY
"$BIN" "${CR[@]}" --rec-model UndatedDTL -p "$TMP/dtl" > "$TMP/resume.log" 2>&1
if grep -q "predates the default switch" "$TMP/resume.log" && \
   [ "$(cut -d' ' -f1 "$TMP/dtl/model_parameters/dup_weighting.txt")" = ale ]; then
  echo "  PASS: a checkpoint without a weighting marker resumes with the ALE weighting"
else
  echo "  FAIL: pre-switch checkpoint resume"; exit 1
fi
python3 "$D/dtl_dup_weight_concordance.py" "$BIN" "$TMP/conc" "$N"
