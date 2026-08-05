#!/usr/bin/env bash
# Regression gate for checkpoint resume of WGD retention (q) and per-event LORe
# resolution (r).
#
# Before this fix, AleState (the struct serialized to mainCheckpoint.txt) never
# stored the WGD retention(s) or LORe resolution(s): they lived only in
# AleEvaluator, and declareWGDs() unconditionally reset them from the
# command-line starting values (q0, r=0.9) on EVERY run, including resumes. So
# a run resumed after the ModelRateOpt2 step (which fits q/r) but before
# Reconciliation would reconcile with the un-fitted starting q/r instead of the
# fitted ones -- silently wrong final output (wgdSummary.txt, resolution
# profiles). See DTL_LORE.md, "Caveat -- checkpoint resume does not restore
# q/r" (now resolved).
#
# This test:
#   1. Runs a synthetic two-WGD DL+LORe analysis to completion (dir "full"),
#      recording the fitted per-WGD q/r from wgdSummary.txt.
#   2. Clones that checkpoint into a second dir ("resume"), rolls the
#      checkpoint's step back from End to Reconciliation (simulating a kill
#      right after the last rate/retention optimization finished, before
#      reconciliation ran), and deletes the old reconciliation output.
#   3. Re-runs the EXACT same command against "resume" (checkCheckpointCmd
#      requires this). It must detect the checkpoint, skip straight to
#      Reconciliation, and reconcile using the FITTED q/r -- not q0=0.2 / r=0.9.
#   4. Asserts the "resume" wgdSummary.txt q/r match the "full" one.
#
# Usage: checkpoint_lore_regression.sh [path/to/kalerax]
set -euo pipefail

BIN="${1:-$(dirname "$0")/../build/bin/kalerax}"
[ -x "$BIN" ] || BIN="${1:-$(dirname "$0")/../build/bin/alerax}"
[ -x "$BIN" ] || { echo "FAIL: binary not found/executable: $BIN"; exit 1; }
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=1

fail() { echo "FAIL: $*"; exit 1; }

echo "binary: $BIN"

# --- synthetic two-clade input (same construction as multi_wgd_regression.sh,
# clade A = true LORe internal WGD, clade B = true AORe internal WGD) ---------
mkdir -p "$TMP/gtrees" "$TMP/maps"
SP="$TMP/species.nw"
echo '(OG:2,((A1:1,A2:1):1,(B1:1,B2:1):1):1);' > "$SP"
NA=24
NB=24
FAM="$TMP/families.txt"
echo '[FAMILIES]' > "$FAM"
for i in $(seq 1 "$NA"); do
  f="famA$i"
  echo "(OG_1,(((A1_1,A1_2),(A2_1,A2_2)),(B1_1,B2_1)));" > "$TMP/gtrees/$f.nw"
  printf "A1_1\tA1\nA1_2\tA1\nA2_1\tA2\nA2_2\tA2\nB1_1\tB1\nB2_1\tB2\nOG_1\tOG\n" \
    > "$TMP/maps/$f.link"
  printf -- "- %s\ngene_tree = %s/gtrees/%s.nw\nmapping = %s/maps/%s.link\n" \
    "$f" "$TMP" "$f" "$TMP" "$f" >> "$FAM"
done
for i in $(seq 1 "$NB"); do
  f="famB$i"
  echo "(OG_1,((A1_1,A2_1),((B1_1,B2_1),(B1_2,B2_2))));" > "$TMP/gtrees/$f.nw"
  printf "A1_1\tA1\nA2_1\tA2\nB1_1\tB1\nB1_2\tB1\nB2_1\tB2\nB2_2\tB2\nOG_1\tOG\n" \
    > "$TMP/maps/$f.link"
  printf -- "- %s\ngene_tree = %s/gtrees/%s.nw\nmapping = %s/maps/%s.link\n" \
    "$f" "$TMP" "$f" "$TMP" "$f" >> "$FAM"
done

# NOTE: all args (including -p) must be byte-identical between the initial run
# and the resume, since checkCheckpointCmd() aborts on any mismatch. So we
# reuse a single output dir throughout, and snapshot it (rather than rerun
# under a different -p) to keep a copy of the ground-truth fitted output.
out="$TMP/out"
common=(-f "$FAM" -s "$SP" --gene-tree-rooting UNIFORM
        --species-tree-search SKIP --fix-rates --seed 1 --rec-model UndatedDL
        -g 1 --wgd A1,A2 --wgd B1,B2 --lore -p "$out")

echo "=== initial run to completion ==="
"$BIN" "${common[@]}" > "$out.log" 2>&1 || fail "initial run exited non-zero"
grep -q "End of AleRax execution" "$out.log" || fail "initial run did not finish"

sumFull="$out/reconciliations/wgdSummary.txt"
[ -s "$sumFull" ] || fail "no wgdSummary.txt written by the initial run"
qaFull=$(awk '!/^#/{print $3; exit}' "$sumFull")
raFull=$(awk '!/^#/{print $4; exit}' "$sumFull")
qbFull=$(awk '!/^#/{c++; if(c==2){print $3; exit}}' "$sumFull")
rbFull=$(awk '!/^#/{c++; if(c==2){print $4; exit}}' "$sumFull")
echo "  fitted (full run):    qA=$qaFull rA=$raFull   qB=$qbFull rB=$rbFull"
awk -v ra="$raFull" -v rb="$rbFull" 'BEGIN{
  if (ra > 0.85) { print "FAIL: clade A r-hat not < 1 (LORe not recovered); test fixture is not exercising the bug"; exit 1 }
  if (rb < 0.95) { print "FAIL: clade B r-hat not ~ 1 (AORe not recovered); test fixture is not exercising the bug"; exit 1 }
}' || exit 1

# --- roll the checkpoint back to right before Reconciliation
# (AleStep::Reconciliation == 6), as if the process had been killed right
# after the last ModelRateOpt2 pass finished, then re-run the SAME command
# (same -p) so checkCheckpointCmd() accepts it as a resume ------------------
rm -rf "$out/reconciliations"
ckpt="$out/checkpoint/mainCheckpoint.txt"
[ -s "$ckpt" ] || fail "no checkpoint written by the initial run"
{ echo 6; tail -n +2 "$ckpt"; } > "$ckpt.new" && mv "$ckpt.new" "$ckpt"

echo "=== resume from checkpoint at step Reconciliation ==="
"$BIN" "${common[@]}" > "$out.resume.log" 2>&1 || fail "resumed run exited non-zero"
grep -q "Checkpoint detected" "$out.resume.log" || fail "resumed run did not detect the checkpoint"
grep -q "End of AleRax execution" "$out.resume.log" || fail "resumed run did not finish"

sumResume="$out/reconciliations/wgdSummary.txt"
[ -s "$sumResume" ] || fail "no wgdSummary.txt written by the resumed run"
qaResume=$(awk '!/^#/{print $3; exit}' "$sumResume")
raResume=$(awk '!/^#/{print $4; exit}' "$sumResume")
qbResume=$(awk '!/^#/{c++; if(c==2){print $3; exit}}' "$sumResume")
rbResume=$(awk '!/^#/{c++; if(c==2){print $4; exit}}' "$sumResume")
echo "  fitted (resumed run):  qA=$qaResume rA=$raResume   qB=$qbResume rB=$rbResume"

awk -v qa1="$qaFull" -v ra1="$raFull" -v qb1="$qbFull" -v rb1="$rbFull" \
    -v qa2="$qaResume" -v ra2="$raResume" -v qb2="$qbResume" -v rb2="$rbResume" 'BEGIN{
  tol = 1e-6
  if (qa1-qa2 > tol || qa2-qa1 > tol) { print "FAIL: resumed qA=" qa2 " != fitted qA=" qa1; exit 1 }
  if (ra1-ra2 > tol || ra2-ra1 > tol) { print "FAIL: resumed rA=" ra2 " != fitted rA=" ra1; exit 1 }
  if (qb1-qb2 > tol || qb2-qb1 > tol) { print "FAIL: resumed qB=" qb2 " != fitted qB=" qb1; exit 1 }
  if (rb1-rb2 > tol || rb2-rb1 > tol) { print "FAIL: resumed rB=" rb2 " != fitted rB=" rb1; exit 1 }
  print "PASS: checkpoint resume reconciles with the fitted q/r, not the starting values"
}'
