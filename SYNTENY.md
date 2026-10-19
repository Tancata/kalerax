# Synteny-informed WGD / LORe inference in kalerax (Layer 1)

From gene trees alone a WGD is not identifiable against a burst of small-scale
duplications (SSDs) on the same branch, and LORe (lineage-specific ohnolog
resolution) is not identifiable against independent SSDs in each descendant
lineage: both pairs of hypotheses produce the same gene-tree topologies. The
decisive evidence is positional. WGD-derived pairs (ohnologs) sit in collinear
duplicated blocks; most SSD pairs (tandem, dispersed, retrotransposed) do not.
`--synteny-dir` adds that information to the reconciliation likelihood.

## Model

For every pair of genes (g, h) of a family that lie **in the same genome**, an
observation y_gh is supplied:

| y | meaning |
|---|---|
| 1 | g and h are an anchor pair of a collinear duplicated block (paralogon), ideally doubly conserved against an unduplicated outgroup |
| 0 | both genes are placed on adequately assembled sequence, but they are not an anchor pair |
| NA | undetermined (short/unplaced scaffold, collapsed region, ...) |

In any reconciled gene tree each same-genome pair has exactly **one** lowest
common ancestor. That node is a WGD split (AORe doubling at the WGD branch, or
a LORe U→R resolution below it) or a non-WGD event (SSD, or occasionally S/T
with transfers). The emission is

    P(y = 1 | LCA is a WGD split) = sigma_W,   P(y = 1 | otherwise) = sigma_N,

and pairs are conditionally independent given the reconciliation. The
emission therefore factorises over gene-tree nodes, hence over CCP clade
splits, and is integrated **exactly** inside the existing amalgamated DP:

- every split k of clade γ gets the counts n₁(k), n₀(k) of same-species pairs
  across it with y = 1 and y = 0;
- each pair's emission is factored by m(y) = max(f_W(y), f_N(y)), so the
  non-WGD weight w_N(k) = ∏ f_N/m ≤ 1 is folded into the CCP split frequency
  (all S, D, T and U-speciation terms), the WGD weight w_W(k) = ∏ f_W/m ≤ 1
  replaces the frequency in the WGD doubling-bracket split term, and
  log ∏_pairs m(y) is added to the family log-likelihood.

All weights are ≤ 1, so every conditional likelihood stays a (sub-)probability.
Weights can be far below the double range: a split carrying hundreds of
collinear pairs that must be explained by a non-WGD event (e.g. ohnologs from an
older, undeclared WGD) has log-weight ≈ n₁·log(sigma_N/sigma_W), or ≈
n₀·log((1-sigma_W)/(1-sigma_N)) for many non-collinear pairs at a WGD split. A
log-weight below −600 is stored as an in-range mantissa times 2^(−256·s) and
the integer s is applied to the CLV term; the family then underflows in double
precision and is re-evaluated exactly in ScaledValue precision (before this,
such families failed with `ll=-inf`).
Because LORe resolution reuses the doubling bracket, a LORe-resolved ohnolog
pair automatically receives the WGD emission, while two parallel SSDs receive
the non-WGD emission.

`sigma_W` and `sigma_N` are two global parameters, estimated by maximum
likelihood together with the WGD retention(s) q (coarse grid + golden-section,
keeping the best point; q is re-refined afterwards because q and sigma_W trade
off). `sigma_N` is fitted in the null (no-WGD) model too, so likelihoods with
and without `--wgd` remain comparable for an LRT. The fitted values are
checkpointed.

**No synteny = bit-for-bit unchanged.** Families without a `.syn` file, or
with no determined pairs, get weights 1 and constant 0.

## Input

`--synteny-dir DIR` reads `DIR/<family>.syn` for every family (missing files
are fine). One observation per line:

    geneA geneB 1
    geneA geneB 0
    geneA geneB NA
    geneX NA          # all pairs involving geneX are undetermined
    # comment

Gene names are the gene-tree leaf labels. Same-genome pairs that are not
listed default to y = 0; pass `--synteny-unlisted-na` to make them NA instead
(e.g. if you list only anchor pairs and cannot vouch for the rest).

`scripts/make_synteny_files.py` writes the `.syn` files from MCScanX
`*.collinearity` output (`--mcscanx`, with `--min-anchors`) or from a plain
anchor-pair list (`--anchors`: i-ADHoRe, WGDI, JCVI/MCscan, or a
doubly-conserved-synteny analysis), with optional `--placed-genes`,
`--na-genes` and `--rename` (annotation ID → leaf label).

**Do not filter anchor pairs by a Ks window around the WGD peak.** LORe
ohnologs have systematically lower Ks than AORe ohnologs, so such a filter
removes the LORe signal. Prefer doubly conserved synteny against an
unduplicated outgroup to assign blocks to a WGD (this also separates nested
WGDs). Flag assembly-collapsed regions as NA: collapsed near-identical LORe
duplicates look like single-copy genes and cannot be repaired by the emission.

## Options

| option | meaning |
|---|---|
| `--synteny-dir DIR` | enable the synteny emission |
| `--synteny-params sigmaW,sigmaN` | starting values (default `0.5,0.05`) |
| `--fix-synteny` | keep sigma_W, sigma_N fixed at the given values |
| `--synteny-unlisted-na` | unlisted same-genome pairs are NA instead of y = 0 |
| `--synteny-min-sigmaW X` | fit sigma_W in [X, 1): a declared WGD must leave at least this fraction of its surviving ohnolog pairs collinear (default 0 = unconstrained) |

## Output

- `reconciliations/syntenySummary.txt`: fitted sigma_W, sigma_N, families with
  pairs, and the numbers of collinear / non-collinear pairs.
- The usual per-species event counts, where WGD ohnolog splits are counted
  separately (`wgd_duplications` column) from SSDs (`duplications`), and
  sampled gene trees label WGD splits `W` (NewickEvents) / `:WGD=Y` (NHX).

## Interpretation and caveats

- sigma_W is estimated from the data. When non-collinear duplicate pairs fit a
  WGD better than an SSD (for instance, two parallel SSDs are a priori much less
  likely than one LORe history), the ML solution can be "WGD with lower
  collinearity" (a lower sigma_W) instead of "SSD". If sigma_W is known
  independently (e.g. the anchor fraction among ohnologs of a well-characterised
  AORe region), fix it with `--fix-synteny --synteny-params`. This is most
  relevant for LORe.
- Pairs in the same block are treated as conditionally independent given the
  reconciliation (a composite likelihood). The planned Layer 2 (block-coupled
  resolution histories) addresses that dependence.

## Testing for a WGD with synteny (`--synteny-min-sigmaW`, `scripts/synteny_bootstrap.py`)

Synteny only discriminates when sigma_W is clearly larger than sigma_N. Each
collinear pair across a duplication node multiplies the WGD:SSD odds by
sigma_W/sigma_N, and each non-collinear pair by (1−sigma_W)/(1−sigma_N). With
sigma_W free and no real WGD, the WGD model can set sigma_W ≈ sigma_N, which
makes the emission identical under both hypotheses. In simulations an SSD
burst was then still called a WGD (2ΔlnL of several hundred per 2000 families).

1. **Make collinearity part of H1.** Use `--synteny-min-sigmaW 0.3` (or similar)
   in the WGD fit; the null fit is unchanged.
   - When a real WGD is present, sigma_W is identified (about 0.7 in the
     simulations) and the bound is not binding.
   - On an SSD burst the bound cut 2ΔlnL from about 190 to about 12–13 (1200
     families) or to 0 (2000 root-origin families).
   - The value matters little: fixed values from 0.4 to 0.95 behaved alike.
   - A small residual can remain, because with a small q the WGD can sit
     selectively on the few nodes whose SSD pairs happen to be collinear.
   - With several declared WGDs, sigma_W is shared. A secure WGD therefore
     pins it for the others, so testing a candidate conditional on known WGDs
     needs no external value.
2. **Calibrate the statistic by parametric bootstrap, not χ².** Under H0, q sits
   on its boundary and sigma_W is not identified. The pair-level emission is
   also a composite likelihood that ignores within-block dependence: in
   simulations with block-correlated collinearity among SSD pairs, 2ΔlnL
   stayed in the hundreds even with sigma_W ≥ 0.7, while q-hat was only 1–2%.
   `scripts/synteny_bootstrap.py`:
   - takes the observed null and WGD runs;
   - simulates data under the fitted null (rates and sigma_N; optional
     `--block-frac` for collinear SSD blocks, `--na-frac` for NA genes);
   - refits both models identically and reports the null distribution of
     2ΔlnL and a p-value.

   A tiny q-hat together with a large ΔlnL is a warning sign in its own right.
3. **Report a sensitivity profile:** 2ΔlnL for a few fixed sigma_W values
   (`--fix-synteny --synteny-params W,N`).

## Validation (bundled tests)

- `tests/synteny_underflow_regression.sh [kalerax]`: a family with 320
  same-genome collinear pairs across a root duplication (log-weight ≈ −942):
  fails with scaling disabled (old behaviour), gives the exact null synteny
  term 320·log sigma_N + 280·log(1−sigma_N), and gives a finite WGD likelihood
  that is unchanged when the scaled path is forced on every split.

- `tests/synteny_bound_regression.py <kalerax> <workdir> [families]`: simulated
  SSD burst — the bound removes most of the spurious WGD signal that a free
  sigma_W lets through; real WGD — the bound does not change the fit.

- `tests/synteny_regression.sh`: (1) an empty synteny dir reproduces the
  no-synteny per-family log-likelihoods exactly; (2) discrimination: on a toy
  with 40 collinear WGD families, 16 SSD families with the *same* ohnolog-grouped
  topology but non-collinear pairs, a 48-family tip-SSD background and 24
  single-copy families, the SSD families are WGD-assigned without synteny
  (0.95 WGD splits per family) and SSD-assigned with it (1.00 SSD per family),
  while the collinear families stay WGD (1.00); (3) checkpoint resume restores
  the fitted sigma.
- `tests/export_consistency.sh`: the exact identity
  ∂lnL/∂ln(sigma_W) = E_post[# collinear pairs at WGD splits] (when only y = 1
  pairs are marked) is checked against exported scenario samples for DL/DTL ×
  AORe/LORe, at sigma_W = sigma_N and at sigma_W = 0.5, sigma_N = 0.05
  (all |z| < 1 in the reference run).
