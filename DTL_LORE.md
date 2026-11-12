# DTL + LORe: lineage-specific rediploidization with transfers

This documents the extension of the **LORe** (lineage-specific ordered
rediploidization) model from the DL model to the **DTL** (duplication–transfer–
loss) model in `kalerax`. Before this, `--lore` was DL-only; it now works with
`--rec-model UndatedDTL` as well.

## What LORe models (recap)

After a WGD the two ohnolog copies need not diverge immediately. A gene lineage
carries a hidden state — **R** (resolved/disomic) or **U** (unresolved/
tetrasomic: a single lineage still holding a not-yet-diverged ohnolog pair). At
the WGD branch the retained duplicate enters **U** with probability `q`; a U
lineage passes through speciation **as a unit** and, on each branch, **resolves
with probability `r`** (resolution = ohnolog divergence, recorded on the branch
where it happens).

- `q = 0` → no WGD (the null);
- `r = 1` → resolution certain at the WGD branch → **AORe** (immediate, WHALE-like);
- `r < 1` → resolution can be **delayed** into descendant lineages → **LORe**.

`r` is fitted **per declared WGD** (one `r` per `--wgd`, pooled over families),
jointly with `q` by maximum likelihood.

**Normalisation.** A U lineage is vertical-only (no D, T or highway channel),
so its per-branch event probabilities are P_S/(P_S+P_L) and P_L/(P_S+P_L).
(Before October 2026 the unnormalised P_S, P_L were used, which made the U
process defective and biased r̂ towards 1; AORe and no-WGD likelihoods are
unaffected.)

## How the DTL port works

The DL LORe machinery is mirrored into `UndatedDTLMultiModel`:

- a per-clade unresolved CLV `_uu` (parallel to the resolved `_uq`), an
  unresolved extinction `_uEU`, and a per-branch resolution probability
  `_resolutionProbs`;
- the **WGD transform** the parent consumes becomes the R/U mixture
  `P_top = (1−q)·P + q·U` (and `E_top = (1−q)·E + q·EU`), where the doubling
  bracket `U`/`EU` is the resolved value when `r = 1`;
- a **U-aware backtrace** (`sampleResolutionCommits`) that reads out, as a
  posterior marginal under the fitted `r`, the species branch on which each
  ohnolog pair resolved. An always-on consistency guard (`check=true`) asserts
  that every sampled-term weight equals the exact inside-CLV (`_uu`, `_uqTop`).

### The one modelling choice: U × transfers

The U-state propagates **vertically only** — through speciation (both daughters
inherit U as a unit), speciation-with-unresolved-loss, and the duplication
bracket. **Transfers and highways act on resolved lineages only**: a tetrasomic
(U) lineage does not transfer, and a transferred copy is born **resolved**,
landing within the destination branch *below* any WGD there. This is the natural
extension of the existing AORe-DTL convention ("a transferred lineage lands below
the WGD, so transfer reads stay raw") and keeps the U recursion transfer-free.
Consequently `_uEU` has no transfer-extinction term, and the backtrace continues
a transferred copy with `btR` (resolved) at the destination.

This is a deliberate first-cut assumption, not a derived necessity; a model in
which a tetrasomic locus can itself transfer (carrying U) would be a strictly
richer — and much harder to identify — alternative.

## Reduction to AORe (the regression gate)

At `r = 1` every U term collapses to the resolved doubling bracket, so
`P_top` and `E_top` become exactly the AORe/WHALE transforms. The model is
therefore **bit-for-bit identical to the existing DTL+WGD (AORe) likelihood**
when `r = 1`, which is the nested null an LRT on `r` rests on.

**Validation performed:**

1. **Bit-for-bit r=1 == AORe.** On a 30-family salmonid set with fixed rates, the
   post-port DTL null and DTL+WGD joint log-likelihoods reproduce the pre-port
   (baseline) values exactly: `−1863.892400` (no WGD) and `−1854.694800`
   (+WGD, q̂=0.36303). Adding the U-state machinery does not perturb the AORe /
   no-WGD likelihood.
2. **Backtrace consistency.** With the `check=true` guard enabled, the U-cell and
   WGD R/U-coin sampling weights equal the inside CLVs at every cell across many
   families × samples (no assertion fired).
3. **Resolution profile localisation.** All inferred ohnolog-resolution events
   fall within the WGD clade; every branch outside it reads exactly 0.
4. **Bundled regression gate.** `tests/dtl_lore_regression.sh` runs DTL null,
   DTL+WGD (AORe), and DTL+WGD+`--lore` on `validation/example-1` and asserts the
   nesting `null ≤ AORe ≤ LORe` with `r̂ ∈ [0,1]`. On the moss (PPAT) terminal
   branch the WGD is detected and LORe correctly reverts to `r = 1` (a terminal
   branch has no descendant lineages, so there is nothing to delay).

## Using it

```sh
kalerax -f families.txt -s species_tree.nw --rec-model UndatedDTL \
  --gene-tree-rooting UNIFORM --species-tree-search SKIP \
  --wgd TAXONA,TAXONB --lore -p out
#   -> out/reconciliations/{totalSpeciesResolutionCounts,wgdSummary}.txt
```

## Caveats

- Inherits all the LORe caveats (identifiability lives in gene-tree topology;
  the resolution *branch* is a posterior readout under one global `r`, not a
  dated time; undated model infers *where*, not *when*).
- The U × transfer interaction is the assumption above; results should be read
  with that in mind. As a robustness check, compare DL+LORe and DTL+LORe: if the
  resolution-timing profile is stable across the two, transfers are not driving
  the LORe signal.
- Like the rest of the WGD/LORe code, this is a research prototype.

## Duplication-split weighting (forward process; fixed)

The undated recursion inherited from ALE/AleRax weighted a duplication split,
and the WGD / LORe doubling-bracket split, once per *unordered* pair of child
clades (`P_D · R(γ′) · R(γ″)`). Speciation and transfer splits, by contrast,
sum both orientations of the split (`P_T · [R_e(γ′) T(γ″) + R_e(γ″) T(γ′)]`).
Under the forward process a duplication or WGD produces two exchangeable
copies, either of which can found either subclade, so every reconciliation was
weighted by 2^-(number of duplication splits) relative to the model.

**kalerax weights every duplication split ×2, in UndatedDTL and UndatedDL.**
This covers the D event and the doubling-bracket split used by the WGD
transform, LORe resolution and the U-state sampler; the likelihood, the
reconciliation backtrace and the samplers use the same weights. The factor is
a fixed constant (`DUP_SPLIT_WEIGHT` in the undated models), and the D term is
written exactly as in upstream AleRax, which adopted the same factor in 8e34d87
(its check: summing the recursion over every labelled rooted gene tree on one
branch gives k!·P(N = k), where before it gave 2^-(k-1) times that).

**History.** UndatedDTL switched to these weights in October 2026 (3f6a5ff),
with `--ale-dup-weights` to restore the old numerics. UndatedDL followed in
November 2026 with the AleRax v1.4.3 merge (0a024d6). The option was then
removed: **0a024d6 is the last build with `--ale-dup-weights`**, so build that
commit (or earlier) to reproduce ALE-weighted results. `--ordered-dup-weights`
(the opt-in flag of 5add0d8) is still accepted as a no-op; `--ale-dup-weights`
now stops with an explanation.

An earlier version of this document said UndatedDL was unaffected because "the
factor is a constant per gene tree". That holds on a single branch but not on a
species tree, where one gene tree can be reconciled with different numbers of
duplication splits (e.g. S versus D followed by two SL). On 291 acipenseriform
families, the per-family difference between the two weightings changed by up to
0.53 lnL between D = 0.05 and D = 0.5, and by up to 0.22 between WGD q = 0.2 and
q = 0.9.

**Evidence**

- **Forward-simulation concordance, UndatedDTL** (`tests/dtl_dup_weight_concordance.py`).
  `scripts/kalerax_sim.py` simulates families under two parameter sets that
  differ only in T; for an exact likelihood the frequency/likelihood ratio of
  every gene-tree shape must be the same under both (χ² ≈ number of shapes).
  The forward weights pass; the ALE weights failed (χ²/df > 4).
- **Forward-simulation concordance, UndatedDL** (`tests/dl_dup_weight_concordance.py`,
  8,000,000 simulated families per parameter set):

  | pair | ALE weights (≤ 0a024d6) | forward weights |
  |---|---|---|
  | D 0.15 vs 0.90 (L 0.60) | χ² 332.2 / 211 shapes, p = 2e-7 | χ² 185.4 / 211, p = 0.90 |
  | WGD q 0.2 vs 0.95 | χ² 600.6 / 281 shapes, p = 3e-25 | χ² 270.1 / 281, p = 0.67 |

  The DL test needs this many families: at 600,000 neither weighting was
  rejected.
- **Independent check with Zombi** (Davín et al. 2020, a dated forward
  simulator with additive transfers to contemporaneous lineages; 12 species,
  4,200 simulated families per configuration of which 2,200–3,100 pass the
  ≥4-species filter; transfer rate varied). Exported D and T event totals were
  compared with Zombi's true observable events. The ALE weighting under-counted
  D by 11–17% and over-counted T by 13–42% in every configuration; the forward
  weights gave D within −2% to +4% and T within +2% to +13%. The numbers are in
  `validation/zombi_dup_weights.md`.
- **Effect in simulations on an 11-taxon tree** (2000 families per replicate,
  UndatedDTL). ALE weighting gave D̂ −5%, T̂ +13% and q̂ about −5%; exported WGD
  splits were about 8% low and transfers about 21% high. With the forward
  weights the biases are about −1%, +4% and −2%.

**Clade values above 1 are expected**

With the ×2 weights a clade value `P_e(γ)` is no longer bounded by 1. Every D or
WGD/LORe doubling split in a scenario carries the factor 2, which is the
labelling factor of the two exchangeable copies. A clade with n leaves has at
most n − 1 such splits, and its value under the old weights is a probability,
so a valid value satisfies `P_e(γ) ≤ 2^(n−1)`; a survival-conditioned family
likelihood likewise satisfies `≤ 2^(genes−1)`. Values above 1 are typical for a
well-retained, low-loss WGD: two retained copies of one gene score about
`2 · R · R ≈ 2`, because the doubling, unlike a D event, carries no `P_D` cost.
`computeProbability` and `computeLogLikelihood` therefore check against these
bounds (`cladeProbaBound`, `familyProbaBound`) rather than against 1; a value
beyond the bound marks the evaluation invalid (lnL −1e15). b959574 instead
rejected every value above 1, which kept optimisers out of the
high-q/low-loss region: on 600 snow carp families with 11 declared WGDs that
build ended 97.8 lnL below the correct optimum, with distorted q. Builds that
logged `proba > 1.0` errors but used the values were correct: their per-family
lnL is bit-identical to the fixed build. Test:
`tests/forward_weight_bound_regression.sh`.

**Comparability and checkpoints**

- Likelihoods under the two weightings are on different scales; never mix them
  in one LRT. **Results from kalerax UndatedDTL before 3f6a5ff (October 2026)
  and UndatedDL before 0a024d6 used the ALE weights** and should be re-run.
- New checkpoints record `2 dl_and_dtl` in `checkpoint/dup_weighting.txt`. A
  checkpoint that used other weights is refused on resume with an explanation:
  no marker (created before October 2026), `1` (`--ale-dup-weights`), or a bare
  `2` from an UndatedDL run (DL before the v1.4.3 merge still used the ALE
  weights). Restart such runs, or resume them with a build at 0a024d6 or earlier.

**Tests:** `tests/dup_weight_regression.sh` checks the option handling, the
checkpoint marker and the refused resumes, then runs both concordance gates.
`tests/export_consistency.sh` checks that sampled scenarios match the exact
posterior in DL and DTL.

## Resolved issue: exported scenarios ignored the WGD (October 2026)

The stock scenario sampler (`MultiModel::backtrace`) picked a speciation using
the WGD-transformed child value (`_uqTop`) but then sampled the child's event
from the raw value, so the q·U mass was never sampled. Ohnolog splits were
exported as SSDs (DL) or as D/T (DTL), and under LORe the species-grouped
ohnologs appeared as a spray of transfers. That mismatch was also the root
cause of the degenerate resample chains described below. The sampler now
descends through the WGD R/U coin and samples the U recursion
(`MultiModel::descend` / `backtraceU`, model hooks `sampleWGDCoin` /
`sampleUEvent`, shared with the resolution-profile sampler). WGD ohnolog splits
are exported as `isWGD` duplications (`wgd_duplications` column, `W` label), and
the backbone is sampled at the fitted r.

## Resolved issue: reconciliation-export stack overflow

The DTL+LORe **optimisation** (q̂, r̂, joint lnL, LRT) was always exact. The
**reconciliation export** (writing the resolution-branch profile and the sampled
gene trees) used to crash ~91 % of the way through a large run.

**Root cause.** The reconciliation backtrace recurses on the DL/TL
"nothing observed, resample" events. Under the DTL+WGD/LORe fitted state the
retention transform suppresses speciation probability at the WGD branch, so for
a few families a cell's terminating-event probability becomes vanishingly small
(precision-dependent) and the resample recursion goes essentially unbounded ->
**stack overflow (SIGSEGV)**. This was a latent property of the stock
amalgamated-likelihood scenario sampler (`MultiModel::backtrace`), not of the
LORe code, exposed by the DTL+LORe likelihood surface.

**Fix** (commit `f83197c`): the scenario backtrace now caps the recursion
*depth* (RAII guard; large-but-shallow gene trees are unaffected) and
discards-and-resamples a runaway scenario; the LORe `btR`/`btU` backtraces have a
per-cell resample cap and a hard per-sample step backstop. A few un-samplable
families contribute fewer (or no) sampled gene trees but are otherwise handled
gracefully; the resolution profile is still written. The bundled regression gate
still passes and the full salmonid DTL+LORe export now completes.

**Second fix — consensus-tree null-deref at the corrected high `q`.** Once the
retention optimiser was fixed (golden-section refinement; `q̂` rises from the
`s`-parameterisation plateau of 0.667 to ~0.78–0.93), the export crashed again,
but for an unrelated reason: at the higher fitted `q` a handful of families yield
a degenerate or empty set of sampled gene trees, and the per-family **consensus
tree** builder (`corax_utree_weight_consensus`) returns null / a malformed
newick. In a release build `NDEBUG` strips corax's internal `assert`s, so this
fell through to `export_newick(null)` (SIGSEGV) or an uncaught `LibpllException`
(SIGABRT). Fix: convert those stripped asserts into runtime guards that return an
empty consensus, wrap the per-tree parse so an unparseable sampled tree is
skipped, and have `AleOptimizer::saveGeneConsensusTree` skip writing when the
consensus is empty (`src/ale/AleOptimizer.cpp`, plus the GeneRaxCore submodule —
see `patches/generaxcore_consensus_export_guards.patch`). The recursion-depth cap
was also lowered (50000 → 8000) to stay well below the 64 MB stack guard page,
since the higher `q` deepens the degenerate resample chain. The resolution
profile is unaffected (it comes from a separate U-aware backtrace, not the
gene-tree consensus). The full salmonid DL+LORe and DTL+LORe exports now complete
at the corrected `q̂`.

**Operational note — gene-tree sample output is large.** A full `-g 30` export
writes per-family sampled gene trees to `reconciliations/all/` (~150 GB for the
16,760-family salmonid set). Only the small `totalSpeciesResolutionCounts.txt`
profile is needed for the LORe readout; delete `all/` after a run if disk is
tight.

**Fixed — checkpoint resume now restores q/r.** `AleState` serialises the WGD
retention(s) `q` and the per-event LORe resolution(s) `r` (by branch label, so
they survive a species-tree node reindex), and `declareWGDs()` no longer
re-declares them from the command-line starting values (q0, r=0.9) on a
resumed run. Resuming a `--wgd … --lore` run after the retention/resolution
optimisation has completed now reconciles with the *fitted* q/r, matching an
uninterrupted run bit-for-bit (`tests/checkpoint_lore_regression.sh`).

## Fixed-parameter evaluation (`--lore-fixed-r`, `--fix-wgd-q`)

- `--lore-fixed-r <r>` (with `--lore` or `--lore-wgd`) holds the resolution
  probability at `r` in [0,1] on every resolvable WGD. Non-target and terminal
  WGDs stay at r = 1. The AORe/LORe model choice and revert are skipped. q is
  still fitted unless `--fix-wgd-q` is given.
- `--fix-wgd-q` holds every WGD retention at its declared `--wgd LABEL:q0`.
  `--fix-rates` alone does not fix q.
- Every run writes `model_parameters/wgd_parameters.txt` at full precision: the
  global rates, q and r per WGD, and σ_W and σ_N. `per_fam_likelihoods.txt` is
  now written with 12 significant digits.
- To evaluate per-family likelihoods L_f(r) at fixed parameters, for example on
  an r grid for block-level models:

      --fix-rates --d D --l L --t T --wgd W:q --fix-wgd-q \
      --lore --lore-fixed-r r --fix-synteny --synteny-params sW,sN \
      --skip-thorough-rates

  Re-evaluating at the free-fit values reproduces the free-fit lnL exactly, and
  `--lore-fixed-r 1` equals the AORe model exactly
  (`tests/lore_fixed_r_regression.py`).
- Invalid command-line arguments now exit with status 1 (previously 0).
