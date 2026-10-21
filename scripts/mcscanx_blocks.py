#!/usr/bin/env python3
"""Family -> collinear-block maps and within-block gene order for block_lore.py.

Inputs
  --syn-dir DIR          kalerax synteny files (<family>.syn; lines 'geneA geneB 1|0|NA')
  --collinearity F [F..] MCScanX *.collinearity files, one per genome. The FIRST is
                         the reference genome; the others are used, in the order
                         given, only for the fallback map.
  --out PREFIX

A family belongs to the reference alignment ("## Alignment k") that contains its
collinear (y = 1) same-genome pair(s); if its pairs fall in several alignments,
the one holding most of them is used (ties: highest MCScanX score). Families
without a collinear pair in the reference genome are unblocked.

Outputs (tab-separated, no header unless stated)
  PREFIX_fam2block_alignment.tsv   family, block (reference alignment) or empty
  PREFIX_fam2block_chrpair.tsv     family, reference chromosome pair (coarser)
  PREFIX_fam2block_fallback.tsv    family, block; unblocked families take the
                                   alignment of the next genome that has one
  PREFIX_order_alignment.tsv       family, block, position (index of the family's
                                   first anchor pair along the alignment) -> HMM
  PREFIX_coverage.tsv              (header) one row per map: blocked families,
                                   blocks with >= 2 families, median/max size
Use the same synteny tool output that produced the .syn files.
"""
import argparse, collections, glob, os, re, statistics


def parse_collinearity(path):
    alns, pair2aln, cur = {}, {}, None
    for line in open(path):
        if line.startswith("## Alignment"):
            m = re.match(r"## Alignment (\d+): score=(\S+) e_value=(\S+) N=(\d+) (\S+)&(\S+) (\w+)", line)
            cur = int(m.group(1))
            alns[cur] = dict(score=float(m.group(2)), chrA=m.group(5), chrB=m.group(6), pairs=[])
        elif line.startswith("#") or not line.strip():
            continue
        else:
            f = line.rstrip("\n").split("\t")
            a, b = f[1].strip(), f[2].strip()
            alns[cur]["pairs"].append((a, b))
            pair2aln.setdefault(frozenset((a, b)), []).append(cur)
    return alns, pair2aln


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--syn-dir", required=True)
    ap.add_argument("--collinearity", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genomes = []
    for f in a.collinearity:
        label = os.path.basename(f).replace(".collinearity", "")
        alns, p2a = parse_collinearity(f)
        genomes.append((label, alns, p2a))

    fam_pairs = {}
    for f in sorted(glob.glob(os.path.join(a.syn_dir, "*.syn"))):
        fam = os.path.basename(f)[:-4]
        fam_pairs[fam] = [tuple(x[:2]) for x in (l.split() for l in open(f)) if len(x) == 3 and x[2] == "1"]

    genes = set()
    for _, alns, _ in genomes:
        for k in alns.values():
            for x, y in k["pairs"]:
                genes.update((x, y))
    n_pairs = 0  # collinear pairs whose genes belong to a supplied genome
    found = 0

    def majority(fam, gi):
        label, alns, p2a = genomes[gi]
        hits = [k for p in fam_pairs[fam] for k in p2a.get(frozenset(p), [])]
        if not hits:
            return None
        c = collections.Counter(hits)
        top = max(c.values())
        return max((k for k, v in c.items() if v == top), key=lambda k: alns[k]["score"])

    maps = {"alignment": {}, "chrpair": {}, "fallback": {}}
    order = []
    for fam in fam_pairs:
        for p in fam_pairs[fam]:
            if p[0] in genes or p[1] in genes:
                n_pairs += 1
                found += any(frozenset(p) in g[2] for g in genomes)
        ref = majority(fam, 0)
        label, alns, _ = genomes[0]
        if ref is not None:
            maps["alignment"][fam] = f"{label}:{ref}"
            maps["chrpair"][fam] = f"{label}:" + "&".join(sorted((alns[ref]["chrA"], alns[ref]["chrB"])))
            anchors = set(map(frozenset, fam_pairs[fam]))
            pos = min(k for k, p in enumerate(alns[ref]["pairs"]) if frozenset(p) in anchors)
            order.append((fam, maps["alignment"][fam], pos))
        fb = None
        for gi in range(len(genomes)):
            k = majority(fam, gi)
            if k is not None:
                fb = f"{genomes[gi][0]}:{k}"
                break
        if fb:
            maps["fallback"][fam] = fb

    rows = ["map\tblocked_families\tblocks_ge2\tmedian_size\tmax_size"]
    for name, m in maps.items():
        if name == "fallback" and len(genomes) == 1:
            continue
        with open(f"{a.out}_fam2block_{name}.tsv", "w") as fh:
            for fam in fam_pairs:
                fh.write(f"{fam}\t{m.get(fam, '')}\n")
        sizes = collections.Counter(m.values())
        multi = [s for s in sizes.values() if s > 1]
        rows.append(f"{name}\t{len(m)}\t{len(multi)}\t{statistics.median(sizes.values()) if sizes else 0}\t{max(sizes.values()) if sizes else 0}")
    with open(f"{a.out}_order_alignment.tsv", "w") as fh:
        for fam, b, p in order:
            fh.write(f"{fam}\t{b}\t{p}\n")
    open(f"{a.out}_coverage.tsv", "w").write("\n".join(rows) + "\n")
    print(f"{len(fam_pairs)} families; {n_pairs} collinear pairs in the supplied genomes, {found} of them are "
          f"anchor pairs in the collinearity files ({100.0 * found / max(n_pairs, 1):.1f}%; ~100% expected if "
          f"these files produced the .syn data)")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
