"""Forward simulator matching kalerax's undated DL/DTL + WGD (q) + LORe (r) model.

Mirrors src/ale/Undated{DL,DTL}MultiModel.hpp @ 2853924:
  R lineage on branch e (at the top of e, below any WGD on e):
      P_D, P_L, P_T, P_S = (D, L, T, 1) / (1 + D + L + T)
      D -> two R lineages on e
      L -> nothing
      T -> one R lineage on e + one R lineage on d, d ~ Uniform{non-ancestors of e}
           (PARENTS constraint: d != e and d not an ancestor of e); lands raw
      S -> leaf: one observed gene; internal: enter(child) for both children
  enter(f): WGD on f -> U on f with prob q, else R on f; no WGD -> R on f
  U lineage on branch e (r = r_e, painted on the WGD subtree; 1 elsewhere):
      prob r   -> resolve: WGD split, two R lineages on e (commit recorded on e)
      prob 1-r -> leaf e: a single R lineage on e (tetrasomic tip)
                  internal: S with P_S/(P_S+P_L) -> U on both children (raw, no
                  further WGD transform); L with P_L/(P_S+P_L) -> nothing
  Origination: uniform over all species branches, as a raw R lineage on that
  branch; families are conditioned on survival (>=1 observed gene) by rejection.
  The kalerax family filter (<4 covered species dropped) is applied separately.
"""
from __future__ import annotations

import numpy as np


# ----------------------------------------------------------------- species tree
class SpeciesTree:
    def __init__(self, newick: str):
        self.parent, self.children, self.label, self.bl = [], [], [], []
        self._parse(newick.strip().rstrip(";"))
        n = len(self.parent)
        self.n = n
        self.root = [i for i in range(n) if self.parent[i] < 0][0]
        self.leaves = [i for i in range(n) if not self.children[i]]
        # leaf sets (for naming / matching kalerax nodes)
        self.leafset = [None] * n
        for i in self.postorder():
            if not self.children[i]:
                self.leafset[i] = frozenset([self.label[i]])
            else:
                self.leafset[i] = frozenset().union(*[self.leafset[c] for c in self.children[i]])
        # ancestors including self
        self.anc = []
        for i in range(n):
            a, j = set(), i
            while j >= 0:
                a.add(j)
                j = self.parent[j]
            self.anc.append(a)
        self.candidates = [np.array([d for d in range(n) if d not in self.anc[e]]) for e in range(n)]
        self.by_label = {self.label[i]: i for i in self.leaves}

    def _parse(self, s):
        stack, cur, i = [], None, 0
        def new(par):
            self.parent.append(par); self.children.append([]); self.label.append(""); self.bl.append(0.0)
            k = len(self.parent) - 1
            if par >= 0:
                self.children[par].append(k)
            return k
        cur = new(-1)
        while i < len(s):
            c = s[i]
            if c == "(":
                cur = new(cur); i += 1
            elif c == ",":
                cur = new(self.parent[cur]); i += 1
            elif c == ")":
                cur = self.parent[cur]; i += 1
            elif c == ":":
                j = i + 1
                while j < len(s) and s[j] not in ",()": j += 1
                self.bl[cur] = float(s[i + 1:j]); i = j
            else:
                j = i
                while j < len(s) and s[j] not in ",():": j += 1
                self.label[cur] = s[i:j]; i = j
        # the first created node is the root wrapper; parse yields it correctly

    def postorder(self):
        out, st = [], [(self.root, False)]
        while st:
            v, done = st.pop()
            if done or not self.children[v]:
                out.append(v)
            else:
                st.append((v, True))
                for c in reversed(self.children[v]):
                    st.append((c, False))
        return out

    def lca(self, labels):
        sets = [self.anc[self.by_label[l]] for l in labels]
        common = set.intersection(*sets)
        # deepest common ancestor = the one with smallest leafset
        return min(common, key=lambda k: len(self.leafset[k]))

    def subtree(self, e):
        return [k for k in range(self.n) if e in self.anc[k]]

    def name(self, e):
        return "|".join(sorted(self.leafset[e]))


# ----------------------------------------------------------------- parameters
class Params:
    """Per-branch event probabilities + WGDs.

    D, L, T : global rates (kalerax 'rates'); per-branch overrides via dict.
    wgds    : list of (branch_index, q, r). r painted over the WGD subtree,
              but pinned to 1 for a WGD on a terminal branch (as kalerax).
    """
    def __init__(self, tree: SpeciesTree, D, L, T=0.0, wgds=(), branch_rates=None):
        n = tree.n
        Dv, Lv, Tv = np.full(n, D, float), np.full(n, L, float), np.full(n, T, float)
        for e, (d, l, t) in (branch_rates or {}).items():
            Dv[e], Lv[e], Tv[e] = d, l, t
        s = 1.0 + Dv + Lv + Tv
        self.PD, self.PL, self.PT, self.PS = Dv / s, Lv / s, Tv / s, 1.0 / s
        self.q = np.zeros(n)
        self.r = np.ones(n)
        self.has_wgd = np.zeros(n, bool)
        for e, q, r in wgds:
            self.has_wgd[e] = True
            self.q[e] = q
            if tree.children[e]:  # internal WGD: paint subtree
                for k in tree.subtree(e):
                    self.r[k] = r
        self.dtl = T > 0 or any(t > 0 for (_, _, t) in (branch_rates or {}).values())


# ----------------------------------------------------------------- simulation
class Node:
    __slots__ = ("kids", "leaf", "event", "branch")
    def __init__(self, kids=None, leaf=None, event=None, branch=None):
        self.kids, self.leaf, self.event, self.branch = kids, leaf, event, branch


class FamilySim:
    """Simulate one family; returns gene tree root (None if extinct) + counts."""
    MAXGENES = 5000

    def __init__(self, tree: SpeciesTree, p: Params, rng: np.random.Generator):
        self.t, self.p, self.rng = tree, p, rng

    def run(self, origin):
        n = self.t.n
        self.cnt = {k: np.zeros(n) for k in ("SSD", "WGD", "T_out", "T_in", "L", "tetra_tip", "U_loss")}
        self.ngenes = 0
        root = self._R(origin)
        return root

    # iterative-safe recursion depth is small for realistic rates
    def _R(self, e):
        p, rng, t = self.p, self.rng, self.t
        u = rng.random()
        if u < p.PD[e]:
            self.cnt["SSD"][e] += 1
            return self._join(self._R(e), self._R(e), "D", e)
        u -= p.PD[e]
        if u < p.PL[e]:
            self.cnt["L"][e] += 1
            return None
        u -= p.PL[e]
        if u < p.PT[e]:
            cands = t.candidates[e]
            d = int(cands[rng.integers(len(cands))])
            self.cnt["T_out"][e] += 1
            self.cnt["T_in"][d] += 1
            return self._join(self._R(e), self._R(d), "T", e)
        # speciation
        if not t.children[e]:
            self.ngenes += 1
            if self.ngenes > self.MAXGENES:
                raise OverflowError
            return Node(leaf=t.label[e], branch=e)
        f, g = t.children[e]
        return self._join(self._enter(f), self._enter(g), "S", e)

    def _enter(self, f):
        if self.p.has_wgd[f] and self.rng.random() < self.p.q[f]:
            return self._U(f)
        return self._R(f)

    def _U(self, e):
        p, rng, t = self.p, self.rng, self.t
        if rng.random() < p.r[e]:
            self.cnt["WGD"][e] += 1
            return self._join(self._R(e), self._R(e), "W", e)
        if not t.children[e]:
            self.cnt["tetra_tip"][e] += 1
            return self._R(e)
        if rng.random() < p.PS[e] / (p.PS[e] + p.PL[e]):
            f, g = t.children[e]
            return self._join(self._U(f), self._U(g), "S", e)
        self.cnt["U_loss"][e] += 1
        return None

    @staticmethod
    def _join(a, b, ev, e):
        if a is None:
            return b
        if b is None:
            return a
        return Node(kids=(a, b), event=ev, branch=e)


def to_newick(node, counter):
    """Leaf names <species>_g<k>; returns newick string and mapping list."""
    mapping = []
    def rec(v):
        if v.kids is None:
            counter[0] += 1
            nm = f"{v.leaf}_g{counter[0]}"
            mapping.append((nm, v.leaf))
            return nm
        return "(" + rec(v.kids[0]) + "," + rec(v.kids[1]) + ")"
    # iterative fallback not needed for MAXGENES-limited trees
    s = rec(node) + ";"
    return s, mapping


def simulate_families(tree: SpeciesTree, p: Params, n_fam: int, seed: int,
                      min_species: int = 4, origin_strategy="UNIFORM",
                      keep_all=False, keep_nodes=False):
    """Simulate surviving families (rejection on extinction) and apply the
    kalerax coverage filter (>= min_species distinct species) unless keep_all.

    Returns list of dicts: newick, mapping, counts (per-branch arrays), origin,
    n_genes, n_species; plus a dict of filter statistics.
    """
    import sys
    sys.setrecursionlimit(100000)
    rng = np.random.default_rng(seed)
    sim = FamilySim(tree, p, rng)
    out, stats = [], dict(attempts=0, extinct=0, filtered=0, overflow=0)
    nodes = np.arange(tree.n) if origin_strategy == "UNIFORM" else np.array([tree.root])
    while len(out) < n_fam:
        stats["attempts"] += 1
        o = int(nodes[rng.integers(len(nodes))])
        try:
            root = sim.run(o)
        except OverflowError:
            stats["overflow"] += 1
            continue
        if root is None:
            stats["extinct"] += 1
            continue
        nwk, mp = to_newick(root, [0])
        nsp = len({s for _, s in mp})
        if not keep_all and nsp < min_species:
            stats["filtered"] += 1
            continue
        out.append(dict(newick=nwk, mapping=mp, origin=o, n_genes=len(mp), n_species=nsp,
                        counts={k: v.copy() for k, v in sim.cnt.items()},
                        **({"root": root} if keep_nodes else {})))
    return out, stats


def write_family_set(fams, outdir, prefix="fam"):
    """Write gene trees + mapping files + a kalerax families file. Gene trees
    for all families are written into ONE directory (2 files per family)."""
    import os
    os.makedirs(os.path.join(outdir, "trees"), exist_ok=True)
    lines = ["[FAMILIES]"]
    for i, f in enumerate(fams):
        nm = f"{prefix}{i}"
        tp = os.path.join(outdir, "trees", nm + ".nw")
        mp = os.path.join(outdir, "trees", nm + ".link")
        with open(tp, "w") as fh:
            fh.write(f["newick"] + "\n")
        with open(mp, "w") as fh:
            fh.write("".join(f"{g}\t{s}\n" for g, s in f["mapping"]))
        lines += [f"- {nm}", f"gene_tree = {tp}", f"mapping = {mp}"]
    ff = os.path.join(outdir, "families.txt")
    with open(ff, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return ff


# ----------------------------------------------------------- branch lengths
def species_depths(tree: SpeciesTree, root_stem=0.05):
    """top/bottom depth (subs/site) of every species branch"""
    top, bot = np.zeros(tree.n), np.zeros(tree.n)
    for v in reversed(tree.postorder()):  # preorder
        if tree.parent[v] < 0:
            top[v] = 0.0
            bot[v] = root_stem if tree.bl[v] == 0 else tree.bl[v]
        else:
            top[v] = bot[tree.parent[v]]
            bot[v] = top[v] + tree.bl[v]
    return top, bot


def to_newick_bl(node, tree, rng, top, bot, fam_rate=1.0, noise_sd=0.3, min_bl=1e-4):
    """Newick with branch lengths. Events on branch e are placed uniformly in
    (max(lower bound, top[e]), bot[e]); S nodes and leaves at bot[e]."""
    mapping = []
    counter = [0]

    def depth(v, lb):
        e = v.branch
        if v.kids is None or v.event == "S":
            return bot[e]
        lo = max(lb, top[e])
        return lo if lo >= bot[e] else rng.uniform(lo, bot[e])

    def rec(v, d_parent):
        d = depth(v, d_parent)
        bl = max(d - d_parent, 0.0) * fam_rate * np.exp(rng.normal(0, noise_sd) - noise_sd ** 2 / 2)
        bl = max(bl, min_bl)
        if v.kids is None:
            counter[0] += 1
            nm = f"{v.leaf}_g{counter[0]}"
            mapping.append((nm, v.leaf))
            return f"{nm}:{bl:.6f}"
        return "(" + rec(v.kids[0], d) + "," + rec(v.kids[1], d) + f"):{bl:.6f}"

    d0 = depth(node, 0.0)
    s = "(" + rec(node.kids[0], d0) + "," + rec(node.kids[1], d0) + ");" if node.kids else None
    return s, mapping


# ----------------------------------------------------------- synteny emission
def synteny_lines(root, rng, sigma_w=0.7, sigma_n=0.03, na_frac=0.0, block=False, block_frac=None):
    """Emission-level synteny observations for every same-genome gene pair.
    Leaf names follow to_newick() order (<species>_g<k>). The LCA event of each
    pair decides the emission: 'W' (WGD / LORe resolution split) -> y=1 w.p.
    sigma_w, otherwise (D, T) -> sigma_n. block=True draws one y per (LCA node,
    species), shared by all pairs of that node in that genome (block-level
    collinearity, violating pairwise independence). block_frac in [0,1]: the
    fraction of (node, species) groups that share one draw (block=True is
    block_frac=1). na_frac: each gene is NA with this probability."""
    if block_frac is None:
        block_frac = 1.0 if block else 0.0
    counter = [0]
    out = []

    def rec(v):
        if v.kids is None:
            counter[0] += 1
            return {v.leaf: [f"{v.leaf}_g{counter[0]}"]}
        a, b = rec(v.kids[0]), rec(v.kids[1])
        for s in set(a) & set(b):
            sig = sigma_w if v.event == "W" else sigma_n
            yb = int(rng.random() < sig)
            shared = rng.random() < block_frac
            for x in a[s]:
                for y in b[s]:
                    yy = yb if shared else int(rng.random() < sig)
                    out.append(f"{x} {y} {yy}")
        for s, l in b.items():
            a.setdefault(s, []).extend(l)
        return a

    allg = rec(root)
    if na_frac > 0:
        for s, gl in allg.items():
            for g in gl:
                if rng.random() < na_frac:
                    out.append(f"{g} NA")
    return out
