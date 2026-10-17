#!/usr/bin/env bash
# Exported-scenario consistency gate: the WGD/LORe-aware sampler must draw from
# the exact posterior. Compares, for DL/DTL x AORe/LORe, the exact posterior
# expectation E[# marked pairs at WGD splits] = d lnL / d ln(sigmaW) with its
# Monte-Carlo estimate from exported scenarios (z-test, |z| < 4). The last two
# configurations use the --ordered-dup-weights split weight (2) under DTL.
# usage: export_consistency.sh [path/to/build/bin]
set -euo pipefail
BINDIR="${1:-$(dirname "$0")/../build/bin}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
W=$TMP; echo '(OG:2,((A1:1,A2:1):1,(B1:1,B2:1):1):1);' > $W/sp.nw
cat > $W/g.nw <<'E'
(OG_1,((A1_1,A2_1),((B1_1,B2_1),(B1_2,B2_2))));
(OG_1,((A1_1,A2_1),((B1_1,B1_2),(B2_1,B2_2))));
(OG_1,((A1_1,A2_1),((B1_1,B2_2),(B1_2,B2_1))));
(OG_1,((A1_1,A2_1),((B1_1,B2_1),(B1_2,B2_2))));
E
printf "A1_1\tA1\nA2_1\tA2\nB1_1\tB1\nB1_2\tB1\nB2_1\tB2\nB2_2\tB2\nOG_1\tOG\n" > $W/map.link
printf "B1_1 B1_2 1\nB2_1 B2_2 1\n" > $W/syn.txt
status=0
for cfg in "DL 0.6 1.0 0.1 0.1" "DL 0.6 0.5 0.5 0.05" "DTL 0.6 1.0 0.1 0.1" "DTL 0.6 0.5 0.5 0.05" "DL 0.2 0.3 0.5 0.05" "DTL 0.6 1.0 0.1 0.1 2" "DTL 0.6 0.5 0.5 0.05 2"; do
  set -- $cfg
  out=$("$BINDIR/export_consistency" $W/sp.nw $W/g.nw $W/map.link B1 B2 $2 $3 $1 30000 $W/syn.txt $W/t.trees 17 $4 $5 ${6:-1} | head -1)
  exact=$(echo "$out" | sed -E 's/.*dlnL\/dln\(sigmaW\)=([0-9.eE+-]+).*/\1/')
  python3 - $W/t.trees "$exact" "$cfg" <<'E' || status=1
import sys, re, math
marked = {frozenset(("B1_1", "B1_2")), frozenset(("B2_1", "B2_2"))}
def parse(s):
    s = s.strip().rstrip(';'); pos = 0
    def node():
        nonlocal pos
        kids = []
        if s[pos] == '(':
            pos += 1; kids.append(node())
            while s[pos] == ',': pos += 1; kids.append(node())
            pos += 1
        m = re.match(r'[^,():;]*', s[pos:]); lab = m.group(0); pos += len(lab)
        if pos < len(s) and s[pos] == ':':
            pos += 1; m = re.match(r'[0-9.eE+-]*', s[pos:]); pos += len(m.group(0))
        return (lab, kids)
    return node()
def leaves(n): return [n[0]] if not n[1] else [x for k in n[1] for x in leaves(k)]
xs = []
for line in open(sys.argv[1]):
    if not line.strip(): continue
    st = [parse(line)]; c = 0
    while st:
        n = st.pop(); st.extend(n[1])
        if n[0] == 'W' and len(n[1]) == 2:
            L, R = set(leaves(n[1][0])), set(leaves(n[1][1]))
            c += sum(1 for p in marked if len(p & L) == 1 and len(p & R) == 1)
    xs.append(c)
m = sum(xs) / len(xs); sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1)); se = sd / math.sqrt(len(xs))
ex = float(sys.argv[2]); z = (m - ex) / se
ok = abs(z) < 4
print(f"  [{sys.argv[3]}] exact={ex:.4f} sampled={m:.4f} se={se:.4f} z={z:+.2f} {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
E
done
[ $status -eq 0 ] && echo "PASS: exported scenarios match the exact posterior (DL/DTL, AORe/LORe, synteny)"
exit $status
