# Whole-genome duplication (WGD) inference in kalerax

> **Naming.** The tool/binary described here is **`kalerax`**, a fork of
> [AleRax](https://github.com/BenoitMorel/AleRax). Throughout this report
> "kalerax" means this fork (the `kalerax` binary); "AleRax" means the upstream
> framework it builds on. (This fork was previously called *valerax*.)

This is a working report on **kalerax**, an extension to
[AleRax](https://github.com/BenoitMorel/AleRax)
that adds **whole-genome duplication (WGD) modelling** to its undated
gene-tree/species-tree reconciliation likelihood, with a per-branch retention
parameter `q`. The goal is to let kalerax test hypothetical WGDs (and estimate
how many duplicated copies were retained) the way
[WHALE](https://github.com/arzwa/Whale.jl) (Zwaenepoel & Van de Peer 2019) does,
while reusing AleRax's fast amalgamated-likelihood machinery and extending it to
the duplication–transfer–loss (DTL) model.

> **Changes (October 2026).**
> 1. *Export fix.* Exported reconciliations (event counts, sampled gene trees,
>    RecPhyloXML, transfers) previously ignored the WGD: the scenario sampler
>    chose a speciation with the WGD-transformed child value but sampled the
>    child from the raw value, so ohnolog splits were exported as SSDs (24 of 24
>    on a toy with q̂ = 0.999) and, under LORe/DTL, as spurious transfers. The
>    sampler now draws the WGD R/U choice and the U state with the exact inside
>    weights; WGD ohnolog splits are counted in a new `wgd_duplications` column
>    (the `duplications` column is SSD-only) and labelled `W` / `:WGD=Y` in
>    sampled gene trees. LORe backbones are sampled at the fitted r (the former
>    r = 1 workaround is gone). Validated against the exact posterior
>    (`tests/export_consistency.sh`).
> 2. *Normalised U state.* An unresolved (U) lineage has no duplication (DL) or
>    transfer/highway (DTL) channel, so its speciation and loss probabilities are
>    renormalised over P_S + P_L. Previously the U process lost (1−r)·P_D (DL) or
>    (1−r)·(P_D + P_T + highways) (DTL) probability per branch, which biased r̂
>    towards 1. AORe (r = 1) and no-WGD likelihoods are unchanged bit-for-bit.
> 3. *Precision rebuild.* A family whose likelihood underflows double precision
>    is rebuilt in high precision; this rebuild previously dropped the WGD q and
>    LORe r, so the largest families were silently evaluated (and exported)
>    without the WGD after `resetAllPrecisions()` (final lnL off by 323 units on
>    `tests/precision_rebuild_regression.sh`). The state is now re-applied.
> 4. *Optimiser.* The q/r golden-section refinement now evaluates the interval
>    endpoints and keeps the best point seen, so q̂ = 0 (no WGD) is reachable
>    exactly (AORe can no longer score below the null).
> 5. *Synteny.* `--synteny-dir` (see [`SYNTENY.md`](SYNTENY.md)).

> **Status:** research prototype. The model is implemented and validated against
> WHALE on a small published dataset. It is a *scaffold* — see *Limitations*.

---

## 1. What was implemented

### The model
At a WGD, a surviving lineage instantaneously **doubles**; the duplicate is
**retained with probability `q`** and lost with probability `1 − q` (the
asymmetric "1-or-2" retention model used by WHALE / Rabier et al. 2014).
Crucially, **`q = 0` recovers the ordinary no-WGD model**, which is the null
hypothesis a likelihood-ratio test rests on.

This gives two closed-form transforms applied at the top of a species branch
(`E` = extinction probability of the branch, `P` = clade likelihood at the
branch, `Σ_splits` = sum over conditional-clade splits):

```
E_top = (1 − q)·E + q·E²
P_top = (1 − q)·P + q·[ 2·E·P + Σ_splits ]
```

### No tree surgery
Rather than inserting a degree-2 node into the species tree (which fights
corax's node arrays), the WGD is attached to the **top of an existing branch**
and the transform is applied only to the values a branch's **parent** consumes.
The rule that keeps this correct:

- reads of a **child** branch (during speciation / speciation-loss) go through
  the transformed `…Top` arrays;
- reads of a branch's **own** events (duplication, loss, transfer, origination)
  stay raw, because the WGD sits *above* the branch while those events happen
  *within* it.

With no WGD declared, every `…Top` array is an exact copy of its raw
counterpart, so the likelihood is **bit-for-bit identical to stock AleRax**
(this is the primary regression test).

### Where it lives
- `src/ale/UndatedDLMultiModel.hpp` — DL model (duplication + loss).
- `src/ale/UndatedDTLMultiModel.hpp` — DTL model (adds transfers). Here a
  transferred lineage is taken to land *within* the destination branch, below
  any WGD there, so transfer/highway reads stay raw; only the vertical
  parent→child reads use the transformed arrays.
- `src/ale/MultiModel.hpp`, `AleEvaluator.*`, `AleArguments.*`, `ale.cpp` —
  plumbing: a `setWGD` broadcast, the CLI option, and a retention optimizer.

### Using it
```sh
kalerax -f families.txt -s species_tree.nw --rec-model UndatedDL \
       --gene-tree-rooting UNIFORM --species-tree-search SKIP \
       --wgd PPAT --wgd ATHA,ATRI -p output_dir
```
`--wgd LABEL` places a WGD on a terminal branch; `--wgd LABEL1,LABEL2` places it
on the branch leading to the LCA of two taxa; an optional `:q0` sets the
starting retention. The retention `q` of each declared WGD is then estimated by
maximum likelihood (reparameterised as `q = r/(1+r)`, `r ≥ 0`), in a dedicated
optimisation pass that runs even when D/L rates are fixed — so `q` can be
estimated against a fixed-rate background for a clean likelihood-ratio test.

### Controlling output volume (`--summary-only`)
A normal run writes the **full per-sample reconciliation** of every gene family
to `reconciliations/all/` — for each family, the sampled gene trees, a
RecPhyloXML file *per sample*, and per-sample event/transfer counts. On large
datasets (many families × many `-g` samples) this is the dominant cost and can
reach **tens of GB**, even though most analyses only consume the summaries.

`--summary-only` suppresses that bulk while keeping everything an aggregate
analysis needs:

| Kept | Dropped |
|---|---|
| `reconciliations/summaries/<fam>_*` (consensus tree, mean species-event / transfer / resolution counts) | `reconciliations/all/<fam>_*` (per-sample trees, RecPhyloXML, `.alerec`, per-sample event/transfer counts) |
| `reconciliations/{totalSpeciesEventCounts,totalTransfers,totalSpeciesResolutionCounts,totalSpeciesTetrasomyCounts}.txt` | `highways/candidate_tests/` (per-candidate per-family LL diffs) |
| `reconciliations/wgdSummary.txt`, `species_trees/`, `model_parameters/`, the gene-family origins, accepted highways, and all `*.log` files | the CCP binaries `ccps/*.ccp` (cleaned up, as with `--cleanup-ccp`) |

```sh
kalerax -f families.txt -s species_tree.nw --rec-model UndatedDL \
        --wgd ATHA,ATRI --lore --summary-only -p out
```
The per-family summaries are still produced exactly as in a full run — under the
hood the few per-sample files they are built from are written transiently and
deleted as soon as each family's summary is computed, so peak disk use stays at
roughly one family rather than the whole dataset. The numerical results
(likelihoods, `q̂`, `r̂`, the resolution profile) are identical; only the bulky
per-sample artefacts are absent.

### Tests
`tests/wgd_regression.cpp` (CMake target `wgd_regression`) checks:
1. **no-WGD == stock AleRax** likelihood, to the bit
   (`−58.689981412938` on the bundled `simulated_2` data), and
2. declaring a WGD and sweeping `q → 0` converges monotonically back to that
   value (confirming `q = 0` is the no-WGD null).

---

## 2. Benchmarking against WHALE

We validated against WHALE on its own **`example-1`** dataset (12 land-plant
gene families over 9 species; the inputs are staged here under
`validation/example-1/`). Both tools use the **same gene→species mapping** (the
species code is the label prefix before the first `_`) and the **same two
hypothetical WGDs** WHALE's documentation uses:

- **PPAT** — the *Physcomitrella* (moss) terminal branch;
- **ATHA–ATRI LCA** — the angiosperm crown branch.

As a sanity check, the WHALE install reproduces its documented reference
likelihoods exactly (e.g. `−592.0186` summed over the 12 families), so the
comparison is against a faithful WHALE.

**Maximum-likelihood retention estimates** (WHALE: constant-rate DLWGD, η fixed;
AleRax: UndatedDL). p-values use the boundary-correct ½χ²₀+½χ²₁ mixture, because
`q = 0` sits on the edge of the parameter space:

| tool | WGD branch | q̂ | 2·ΔLL | p | detect @0.05 |
|---|---|---|---|---|---|
| **AleRax** | PPAT | 0.31 | 3.31 | 0.034 | **yes** |
| **WHALE** | PPAT | 0.19 | 0.51 | 0.24 | no |
| **AleRax** | angiosperm | 0.00 | 0.00 | 1.0 | no |
| **WHALE** | angiosperm | 0.00 | 0.00 | 1.0 | no |

**What agrees:** both tools put the retention signal on the **same branch**
(PPAT), and both drive the **angiosperm `q̂` to exactly 0** — i.e. they agree
there is no evidence for the deeper WGD in this tiny dataset. The agreement on
the null is exact.

**What differs:** AleRax estimates a *stronger* moss WGD (higher `q̂`, ~6× larger
LRT) and would call it significant, where WHALE would not. This is expected and
not a bug — WHALE is a **dated** model (it uses the branch lengths and slices
time), whereas AleRax `UndatedDL` ignores branch lengths; the two also use
different root/origination priors and place the WGD at slightly different points
on the branch. WHALE's `q̂` is stable across its time discretisation, so the gap
is a genuine model difference, not a numerical artifact. WHALE's own docs treat
`example-1` as illustrative and caution against over-reading an LRT on so little
data — good context for both tools' verdicts.

**Bottom line:** the implementation is validated on the null (no false WGD) and
agrees with WHALE on *where* the signal is; it is simply more liberal about
*how strong* it is, consistent with being an undated model.

---

## 3. Effect of allowing transfers (DL vs DTL)

A natural worry: in the DTL model, could **horizontal transfers explain away** a
WGD — i.e. account for the extra gene copies that would otherwise be evidence
for genome doubling? We ported the WGD transform into the DTL model and re-ran
the same analysis (transfers with the realistic `PARENTS` constraint).

> Note (October 2026): the DTL row below was computed with the ALE/AleRax
> duplication weighting, which is now `--ale-dup-weights`. The default DTL
> weighting changed (see `DTL_LORE.md`), so a re-run with the current default
> gives different DTL likelihoods.

| model | null LL | +PPAT LL | q̂[PPAT] | 2·ΔLL | p | detect |
|---|---|---|---|---|---|---|
| **DL** (no transfers) | −362.4 | −360.7 | 0.31 | 3.31 | 0.034 | yes |
| **DTL** (transfers) | −330.2 | −327.9 | **0.34** | **4.71** | **0.015** | yes |

Two observations:

1. **DTL fits the data much better overall** (null −330 vs −362): transfers
   absorb ~32 log-units of gene-tree structure that the DL model otherwise
   charges to duplication and loss.
2. **The WGD signal survives — and slightly strengthens** (`q̂` 0.31 → 0.34,
   LRT 3.3 → 4.7). Allowing transfers did **not** explain away the moss WGD.

**Why** transfers can't absorb it: the duplicated *Physcomitrella* genes form
**moss–moss sister pairs** in the gene trees. A transfer would make a moss gene
sister to *another species'* gene, not to a second moss gene, so the ohnolog
topology is only explainable by duplication/WGD. Transfers soak up the *other*
incongruence in the dataset (hence the much better baseline fit) but leave the
ohnolog pairs to the WGD. The same four multi-copy families drive the signal in
both models, and the angiosperm WGD stays undetected (`q̂ = 0`) either way. The
finding is robust to the transfer constraint (`NONE` gives essentially the same
`q̂`).

---

## 4. Delayed rediploidization (LORe)

The WGD model above assumes the two ohnolog copies diverge **immediately** at the
WGD — *autopolyploid-like, ancestral rediploidization* (**AORe**). But after many
real WGDs the doubled genome stays tetrasomic for a while and the ohnologs
**rediploidize later, independently in different descendant lineages**
(*lineage-specific ordered rediploidization*, **LORe**; e.g. the salmonid Ss4R).
This fork adds an (opt-in, DL-only) LORe model.

### The model
A gene lineage carries a hidden state: **R** (resolved/disomic, ordinary) or
**U** (unresolved/tetrasomic — a single lineage holding a not-yet-diverged
ohnolog pair). At the WGD branch the retained duplicate enters the **U** state
(with probability `q`). A U lineage passes through speciation **as a unit** (both
daughters inherit U) and, on each branch, **resolves with probability `r`** — and
*resolution is the ohnolog divergence*, placed on the branch where it happens. So:

- `q = 0` → no WGD (the null);
- `r = 1` → resolution is certain at the WGD branch → recovers AORe / WHALE exactly;
- `r < 1` → resolution can be delayed to descendant lineages → LORe.

Concretely, the resolved extinction `E` and clade likelihood `R` gain unresolved
counterparts `EU`, `U`, with `EU = r·E² + (1−r)(…)` and
`U = r·dupBracket + (1−r)(speciate-U + speciate-loss-U)`, and the WGD branch mixes
them: `transform = (1−q)·R + q·U` (and likewise for extinction). At `r = 1` this is
**bit-for-bit** the AORe transform (a hard regression gate).

`r` is fitted **per declared WGD** — one `r` per `--wgd`, pooled across all
families (not per-family, not free per-branch *within* a WGD) — jointly with the
`q`'s by maximum likelihood. Each WGD paints its own subtree (the WGD branch plus
all descendants) with its own `r`; branches below no WGD stay at `r = 1` (inert
there, since no unresolved mass reaches them). This requires the WGD subtrees to
be **disjoint**: with multiple `--wgd` under `--lore`, if one WGD is ancestral to
another the run aborts (per-event `r` is unidentifiable on the shared branches —
the single U-channel per branch cannot record which WGD a U-lineage descends
from). A WGD on a **terminal** branch governs only its own branch (no descendant
lineage to delay resolution to), so its `r` is **pinned to 1 (AORe)** and excluded
from the free parameters. Because the LORe model *nests* AORe (`r = 1`), the fit
is guarded to never score below the AORe optimum: with no delayed-resolution
signal, `r̂ → 1` and the LRT ties. (Nested WGDs **without** `--lore` are fine — the
per-branch WGD transform is local, so stacked `q`'s compose and each is estimated
independently under AORe.)

**Testing LORe for one WGD only (`--lore-wgd`).** Use `--lore-wgd LABEL[,LABEL2]`
(same label syntax as `--wgd`, repeatable) to fit `r` for a **specific** declared
WGD while every other `--wgd` is pinned to AORe (`r = 1`). Only the named targets
have a free `r`, so the disjointness requirement applies **only among the
targets** — a LORe-tested WGD may now be **nested inside** an AORe-pinned WGD (an
AORe WGD resolves at its own branch and sends no unresolved mass downward, so the
shared branches carry only the inner target's `r`). This is the supported way to
ask "did *this* WGD rediploidize late?" without forcing every declared WGD to be
disjoint. `--lore-wgd` implies LORe and is mutually exclusive with bare `--lore`
(which fits every internal WGD). Non-target WGDs are logged as `(pinned;non-target)`,
terminal ones as `(pinned;terminal)`.

```sh
# Two WGDs, but only test the (nested) inner one for late rediploidization:
kalerax -f families.txt -s species_tree.nw --rec-model UndatedDL \
        --gene-tree-rooting UNIFORM --species-tree-search SKIP \
        --wgd OUTER1,OUTER2 --wgd INNER1,INNER2 --lore-wgd INNER1,INNER2 -p out
```
The fit log then reports a free `r` for the target and `r = 1` for the rest, e.g.
`... r[node 7]=0.59 r[node 6]=1(pinned;non-target)`, and every declared WGD's
fitted `q`, fitted `r`, and resolution status is written to
`out/reconciliations/wgdSummary.txt` (one row per `--wgd`). Pass `--lore-wgd`
more than once to fit several targets at once; the targets must still be
mutually disjoint (only the *non-target* WGDs may be nested freely).

### Reading out *where* rediploidization happened
The resolution **branch** is not a parameter — it is a **posterior marginal**
read out (under each WGD's fitted `r`) by a *U-aware backtrace*: the
reconciliation sampler is extended with the U state and a U→R **commit** event
that records the species branch where each ohnolog pair diverged. Every sampling
weight equals the exact inside-likelihood term (an always-on consistency guard
asserts this), so the commit frequencies are unbiased.

This is written to the reconciliation output:

- `reconciliations/totalSpeciesResolutionCounts.txt` — the **genome-wide
  rediploidization-timing profile**: expected number of ohnolog-divergence
  commits on each species branch, summed over families.
- `reconciliations/summaries/<fam>_meanResolutionCounts.txt` — per family.
- `reconciliations/wgdSummary.txt` — per declared WGD: branch, fitted `q`,
  fitted `r`, and the expected commits in its subtree (so the WGD location, its
  retention, and its AORe/LORe status are self-contained in one file).

Under **AORe** the commits sit on the WGD branch itself; under **LORe** they sit
on the descendant lineage branches (and the WGD branch reads ≈ 0). So the WGD
*doubling* branch is always the user-declared one (in `wgdSummary.txt`), while the
*resolution* profile shows the timing.

### Validation
On a controlled toy (WGD on the LCA branch, D/L fixed):

| gene tree | fitted `r̂` | resolution commits land on |
|---|---|---|
| species-grouped `((Aa,Ab),(Ba,Bb))` (true LORe) | **0.62** | the two lineage/terminal branches (≈0.95 each), WGD branch ≈0 |
| ohnolog-grouped `((Aa,Ba),(Ab,Bb))` (true AORe) | **1.0** | the WGD branch (≈1), terminals 0 |

So the highest-posterior commit branch matches the simulated truth in both cases,
under the fitted `r`, with no per-branch tuning. On `example-1` (no LORe
signal) `r̂ → 1` and the model ties AORe.

> **Use it:** add `--lore` to a `--wgd` run (DL model). Without `--lore`, `r = 1`
> (AORe); the resolution profile is still written and attributes the WGD-derived
> duplications to the WGD branch.

### LORe-specific caveats
- **Identifiability lives in gene-tree topology** (species-grouped vs
  ohnolog-grouped ohnologs). If the CCPs don't carry multi-species ohnolog
  structure, `r` trades off against loss/retention and `r̂` is not meaningful —
  check the input actually contains that signal before interpreting it.
- The resolution branch is a *posterior readout under the fitted `r`* (one `r`
  per declared WGD), not a free per-branch rate; and the timing is a resolution
  **branch**, not a dated time (the undated model infers *where*, not *when*).
- Implemented for both the DL and the DTL model (`--lore` accepts
  `--rec-model UndatedDL` or `UndatedDTL`). In the DTL model the unresolved
  (tetrasomic) state propagates **vertically only**: transfers act on resolved
  lineages and a transferred copy is born resolved (see [`DTL_LORE.md`](DTL_LORE.md)).

---

## 5. Limitations

- **Tiny dataset.** `example-1` is 12 families; significance values here are
  illustrative, not conclusions about plant evolution.
- **Undated model.** Unlike WHALE, AleRax `UndatedDL/DTL` ignore branch lengths,
  which is why absolute `q̂`/LRT differ from WHALE. Only relative comparisons
  (DL vs DTL, branch vs branch) should be trusted across the two tools.
- **Synteny is optional.** From gene trees alone, a WGD is **not identifiable**
  against a burst of small-scale duplications on the same branch (nor LORe
  against parallel SSDs). `--synteny-dir` adds a pair-level collinearity
  emission (Layer 1, [`SYNTENY.md`](SYNTENY.md)); block-coupled resolution
  histories (Layer 2) are a planned extension.
- The DTL WGD/LORe port is covered by empirical regression checks (no-WGD and
  r=1/AORe likelihoods reproduced bit-for-bit; backtrace consistency guard;
  `tests/dtl_lore_regression.sh`) rather than a dedicated C++ unit test.

---

## 6. Reproducing

```sh
# build (no MPI):
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DDISABLE_MPI=ON
cmake --build . --target kalerax wgd_regression lore_marginal export_consistency -j4

# WGD regression + LORe r=1 gate (no-WGD == stock, q->0, r=1 == WHALE):
./bin/wgd_regression <species_tree> <gene_trees> <mapping>

# LORe resolution-branch marginal (recovers the simulated commit branch):
./bin/lore_marginal <species_tree> <geneTree> <taxonA> <taxonB> <AORe|LORe>

# run a WGD + LORe analysis (DL *or* DTL model), writing the resolution profile:
./bin/kalerax -f families.txt -s species_tree.nw --rec-model UndatedDL \
  --gene-tree-rooting UNIFORM --species-tree-search SKIP \
  --wgd ATAXON,BTAXON --lore -p out
#   (--rec-model UndatedDTL also supported for --lore; see DTL_LORE.md)
#   -> out/reconciliations/{totalSpeciesResolutionCounts,wgdSummary}.txt

# test LORe for ONE WGD only (may be nested inside other, AORe-pinned WGDs):
./bin/kalerax -f families.txt -s species_tree.nw --rec-model UndatedDL \
  --gene-tree-rooting UNIFORM --species-tree-search SKIP \
  --wgd OUTER1,OUTER2 --wgd INNER1,INNER2 --lore-wgd INNER1,INNER2 -p out

# suppress the bulky per-family output (GB-scale); keep summaries + logs:
./bin/kalerax ... --wgd ... --lore --summary-only -p out

# DTL+LORe regression gate (null <= AORe <= LORe on validation/example-1,
# terminal (PPAT) and internal (ATHA,ATRI) WGD):
bash tests/dtl_lore_regression.sh ./bin/kalerax

# exported scenarios == exact posterior (needs the export_consistency target):
bash tests/export_consistency.sh ./bin
# precision rebuild keeps the WGD; synteny emission gates:
bash tests/precision_rebuild_regression.sh ./bin/kalerax
bash tests/synteny_regression.sh ./bin/kalerax

# Per-event r gate (disjoint WGDs -> distinct r-hat; terminal pinning; nested
# abort; --lore-wgd single-target incl. nested; --summary-only suppression) on a
# synthetic two-clade input, for UndatedDL AND UndatedDTL:
bash tests/multi_wgd_regression.sh ./bin/kalerax

# WHALE benchmark inputs and the exact run recipe:
#   validation/example-1/README.md
```

---

*This extension is a prototype for research use. The retention estimates should
be interpreted alongside the limitations above, and cross-checked against a
dated and/or synteny-aware method (e.g. WHALE) before drawing biological
conclusions.*
