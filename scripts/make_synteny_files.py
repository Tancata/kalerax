#!/usr/bin/env python3
"""
Build kalerax synteny (.syn) files from within-genome collinearity results.

For every family in an AleRax/kalerax families file, every pair of genes from
the SAME species is written as
    geneA  geneB  1    # anchor pair of a collinear duplicated block
    geneA  geneB  0    # both genes placed, not an anchor pair
and genes whose position is unreliable (short/unplaced scaffolds, collapsed
regions) can be written as
    geneX  NA          # all pairs involving geneX are undetermined

Anchor input (one of):
  --mcscanx FILE      MCScanX *.collinearity output (all blocks, or only those
                      with >= --min-anchors anchors)
  --anchors FILE      plain whitespace-separated "geneA geneB" pairs (e.g. from
                      i-ADHoRe, WGDI, JCVI/MCscan, or a doubly-conserved-synteny
                      analysis against an unduplicated outgroup)

Do NOT pre-filter anchor pairs by a Ks window around the WGD peak: ohnologs in
lineage-specifically rediploidized (LORe) regions have systematically lower Ks,
so such a filter removes exactly the LORe signal.

Gene names must match the gene-tree leaf labels (use --rename to translate
annotation IDs to leaf labels). Mapping files may be in "gene species" (link)
or "species:gene1;gene2" format.

usage:
  make_synteny_files.py --families families.txt (--mcscanx F | --anchors F)
                        --out SYN_DIR [--placed-genes F] [--na-genes F]
                        [--rename F] [--min-anchors 5]
then run kalerax with --synteny-dir SYN_DIR.
"""
import argparse, itertools, os, re, sys
from collections import defaultdict


def read_families(path):
    fams, cur = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith("- "):
            cur = line[2:].strip(); fams[cur] = {}
        elif "=" in line and cur:
            k, v = [x.strip() for x in line.split("=", 1)]
            fams[cur][k] = v
    return fams


def read_mapping(path):
    gene2sp = {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if ":" in line and ";" in line or (":" in line and len(line.split()) == 1):
            sp, genes = line.split(":", 1)
            for g in genes.split(";"):
                if g:
                    gene2sp[g.strip()] = sp.strip()
        else:
            parts = line.split()
            if len(parts) >= 2:
                gene2sp[parts[0]] = parts[1]
    return gene2sp


def read_mcscanx(path, min_anchors):
    pairs, block, n = set(), [], 0
    def flush():
        if len(block) >= min_anchors:
            pairs.update(block)
    for line in open(path):
        if line.startswith("## Alignment"):
            flush(); block = []
            continue
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 3:
            a, b = f[1].strip(), f[2].strip()
            block.append(frozenset((a, b)))
    flush()
    return pairs


def read_pairs(path):
    pairs = set()
    for line in open(path):
        f = line.split()
        if len(f) >= 2 and not line.startswith("#"):
            pairs.add(frozenset((f[0], f[1])))
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--mcscanx"); g.add_argument("--anchors")
    ap.add_argument("--out", required=True)
    ap.add_argument("--placed-genes", help="genes with reliable positions; all others -> NA")
    ap.add_argument("--na-genes", help="genes to mark NA (e.g. on collapsed/unplaced sequence)")
    ap.add_argument("--rename", help="2-column map: annotation_id leaf_label")
    ap.add_argument("--min-anchors", type=int, default=5)
    a = ap.parse_args()

    ren = {}
    if a.rename:
        for line in open(a.rename):
            f = line.split()
            if len(f) >= 2:
                ren[f[0]] = f[1]
    pairs = read_mcscanx(a.mcscanx, a.min_anchors) if a.mcscanx else read_pairs(a.anchors)
    if ren:
        pairs = {frozenset(ren.get(x, x) for x in p) for p in pairs}
    placed = None
    if a.placed_genes:
        placed = {ren.get(l.split()[0], l.split()[0]) for l in open(a.placed_genes) if l.strip()}
    na = set()
    if a.na_genes:
        na = {ren.get(l.split()[0], l.split()[0]) for l in open(a.na_genes) if l.strip()}

    os.makedirs(a.out, exist_ok=True)
    tot = defaultdict(int)
    for fam, info in read_families(a.families).items():
        gene2sp = read_mapping(info["mapping"])
        bysp = defaultdict(list)
        for gname, sp in gene2sp.items():
            bysp[sp].append(gname)
        lines = []
        nagenes = sorted(gname for gname in gene2sp if gname in na or (placed is not None and gname not in placed))
        for gname in nagenes:
            lines.append(f"{gname} NA")
        nas = set(nagenes)
        for sp, genes in bysp.items():
            for x, y in itertools.combinations(sorted(genes), 2):
                if x in nas or y in nas:
                    continue
                obs = 1 if frozenset((x, y)) in pairs else 0
                lines.append(f"{x} {y} {obs}"); tot[obs] += 1
        tot["NAgenes"] += len(nagenes)
        with open(os.path.join(a.out, fam + ".syn"), "w") as fh:
            fh.write("\n".join(lines) + ("\n" if lines else ""))
    print(f"wrote {a.out}: collinear pairs={tot[1]} non-collinear={tot[0]} NA genes={tot['NAgenes']}", file=sys.stderr)


if __name__ == "__main__":
    main()
