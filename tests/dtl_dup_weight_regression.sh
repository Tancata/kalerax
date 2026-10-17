#!/usr/bin/env bash
# Regression gate for --ordered-dup-weights.
#  (1) UndatedDL: the flag is a no-op (identical per-family likelihoods);
#  (2) UndatedDTL default (no flag) is unchanged by construction (weight 1);
#  (3) UndatedDTL with the flag is concordant with forward simulation while the
#      default is not (tests/dtl_dup_weight_concordance.py; needs python3+numpy).
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
C=(-f "$TMP/fam.txt" -s "$E/tree.nw" --rec-model UndatedDL --gene-tree-rooting UNIFORM
   --species-tree-search SKIP --summary-only --gene-tree-samples 0 --seed 1 --wgd ATHA,ATRI)
"$BIN" "${C[@]}" -p "$TMP/dl0" > /dev/null 2>&1
"$BIN" "${C[@]}" --ordered-dup-weights -p "$TMP/dl1" > /dev/null 2>&1
if cmp -s <(sort "$TMP/dl0/per_fam_likelihoods.txt") <(sort "$TMP/dl1/per_fam_likelihoods.txt"); then
  echo "  PASS: UndatedDL per-family likelihoods identical with/without the flag"
else
  echo "  FAIL: UndatedDL changed by --ordered-dup-weights"; exit 1
fi
python3 "$D/dtl_dup_weight_concordance.py" "$BIN" "$TMP/conc" "$N"
