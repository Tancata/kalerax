#!/usr/bin/env bash
# Regression gate for the synteny (Layer 1) emission in kalerax.
#   1. --synteny-dir pointing at an EMPTY directory reproduces the no-synteny
#      per-family log-likelihoods exactly (weights 1, constant 0).
#   2. Discrimination (T3 design, AORe WGD on B1-B2): ohnolog-grouped families
#      whose ohnolog pairs are collinear stay WGD; families with the SAME
#      topology but non-collinear pairs move to SSD once synteny is supplied
#      (they are WGD-assigned without synteny). A tip-SSD background calibrates
#      the duplication rate and sigma_N.
#   3. Checkpoint resume at the Reconciliation step restores the FITTED
#      (sigma_W, sigma_N), not the starting values.
# usage: synteny_regression.sh [path/to/kalerax]
set -euo pipefail
BIN="${1:-$(dirname "$0")/../build/bin/kalerax}"
[ -x "$BIN" ] || { echo "FAIL: binary not found: $BIN"; exit 1; }
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=1
fail() { echo "FAIL: $*"; exit 1; }
d=$TMP; mkdir -p $d/g $d/m $d/syn $d/empty
echo '(OG:2,((A1:1,A2:1):1,(B1:1,B2:1):1):1);' > $d/sp.nw
echo '[FAMILIES]' > $d/families.txt
mk() { local f=$1 nw=$2 y=$3; shift 3
  echo "$nw" > $d/g/$f.nw
  echo "$nw" | tr -c 'A-Za-z0-9_\n' ' ' | tr ' ' '\n' | grep -E '_' | awk '{split($1,a,"_"); print $1"\t"a[1]}' > $d/m/$f.link
  printf -- "- %s\ngene_tree = %s/g/%s.nw\nmapping = %s/m/%s.link\n" $f $d $f $d $f >> $d/families.txt
  : > $d/syn/$f.syn; for p in "$@"; do echo "$p $y" >> $d/syn/$f.syn; done; }
nw="(OG_1,((A1_1,A2_1),((B1_1,B2_1),(B1_2,B2_2))));"
for i in $(seq 1 40); do mk W$i "$nw" 1 "B1_1 B1_2" "B2_1 B2_2"; done
for i in $(seq 1 16); do mk S$i "$nw" 0 "B1_1 B1_2" "B2_1 B2_2"; done
k=0
for tip in A1 A2 B1 B2; do for i in $(seq 1 12); do k=$((k+1)); y=0; [ $i -eq 1 ] && y=1
  case $tip in
    A1) t="(OG_1,(((A1_1,A1_2),A2_1),(B1_1,B2_1)));";; A2) t="(OG_1,((A1_1,(A2_1,A2_2)),(B1_1,B2_1)));";;
    B1) t="(OG_1,((A1_1,A2_1),((B1_1,B1_2),B2_1)));";; B2) t="(OG_1,((A1_1,A2_1),(B1_1,(B2_1,B2_2))));";;
  esac
  mk D$k "$t" $y "${tip}_1 ${tip}_2"; done; done
for i in $(seq 1 24); do mk C$i "(OG_1,((A1_1,A2_1),(B1_1,B2_1)));" 0; done
C=(-f $d/families.txt -s $d/sp.nw --gene-tree-rooting UNIFORM --species-tree-search SKIP
   --seed 1 -g 10 --rec-model UndatedDL --wgd B1,B2)
echo "binary: $BIN"
"$BIN" "${C[@]}" -p $d/nosyn > $d/nosyn.log 2>&1 || fail "no-synteny run"
"$BIN" "${C[@]}" --synteny-dir $d/empty -p $d/empty_out > $d/empty.log 2>&1 || fail "empty-synteny run"
cmp -s <(sort $d/nosyn/per_fam_likelihoods.txt) <(sort $d/empty_out/per_fam_likelihoods.txt) \
  || fail "empty synteny dir changed the per-family likelihoods"
echo "  PASS: empty synteny dir == no synteny (per-family lnL identical)"
out=$d/syn_out
"$BIN" "${C[@]}" --synteny-dir $d/syn -p $out > $d/syn.log 2>&1 || fail "synteny run"
grp() { # dir group -> "meanSSD meanWGD"
  python3 - "$1/reconciliations/summaries" "$2" <<'E'
import sys, glob, os
D, g = sys.argv[1], sys.argv[2]; s = w = n = 0.0
for f in glob.glob(os.path.join(D, g + "*_meanSpeciesEventCounts.txt")):
    rows = [l.strip().split(", ") for l in open(f)]; h = rows[0]
    iD, iW = h.index("duplications"), h.index("wgd_duplications"); n += 1
    for r in rows[1:]: s += float(r[iD]); w += float(r[iW])
print(f"{s/max(n,1):.3f} {w/max(n,1):.3f}")
E
}
read sN wN <<< "$(grp $d/nosyn S)"; read sW wW <<< "$(grp $out W)"; read sS wS <<< "$(grp $out S)"
echo "  no synteny: S-families WGD splits/fam=$wN | synteny: W WGD/fam=$wW, S SSD/fam=$sS"
awk -v wN=$wN -v wW=$wW -v sS=$sS 'BEGIN{
  if (wN < 0.8) { print "FAIL: without synteny the S families should be WGD-assigned"; exit 1 }
  if (wW < 0.9) { print "FAIL: collinear W families should stay WGD"; exit 1 }
  if (sS < 0.9) { print "FAIL: non-collinear S families should move to SSD"; exit 1 }
  print "  PASS: synteny separates collinear (WGD) from non-collinear (SSD) ohnolog-grouped families"
}' || exit 1
read fitW fitN <<< "$(awk '!/^#/{print $1, $2}' $out/reconciliations/syntenySummary.txt)"
rm -rf "$out/reconciliations"
ck="$out/checkpoint/mainCheckpoint.txt"
{ echo 6; tail -n +2 "$ck"; } > "$ck.new" && mv "$ck.new" "$ck"
"$BIN" "${C[@]}" --synteny-dir $d/syn -p $out > $d/resume.log 2>&1 || fail "resumed run"
grep -q "Checkpoint detected" $d/resume.log || fail "resume did not detect the checkpoint"
read resW resN <<< "$(awk '!/^#/{print $1, $2}' $out/reconciliations/syntenySummary.txt)"
echo "  fitted sigma (full run): $fitW $fitN | resumed: $resW $resN"
[ "$fitW" = "$resW" ] && [ "$fitN" = "$resN" ] || fail "resume did not restore the fitted synteny parameters"
echo "  PASS: checkpoint resume restores the fitted sigma_W, sigma_N"
echo "PASS: synteny regression"
