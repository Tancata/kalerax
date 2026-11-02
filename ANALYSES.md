# kalerax analyses: a how-to

This guide covers the analyses kalerax supports, from a plain reconciliation
to block-level rediploidization timing. Each section gives the question it
answers, its inputs, a command, the outputs and how to read them. Model
details are in [WGD_REPORT.md](WGD_REPORT.md) (WGD, LORe),
[SYNTENY.md](SYNTENY.md) (synteny) and [DTL_LORE.md](DTL_LORE.md) (DTL and
fixed-parameter evaluation).

| Analysis | Question | Needs | Section |
|---|---|---|---|
| Reconciliation | D/T/L rates and events | gene trees | [A](#a-reconciliation-no-wgd) |
| WGD test | Did a WGD occur on branch *w*? | gene trees | [B](#b-wgd-test-aore) |
| WGD placement | On which branch? | gene trees | [C](#c-wgd-placement-scan) |
| LORe | Was rediploidization delayed past speciation? Where? | gene trees | [D](#d-delayed-rediploidization-lore) |
| Synteny | Use collinearity to separate WGD from SSDs | + `.syn` files | [E](#e-adding-synteny) |
| Fixed-parameter evaluation | Per-family likelihood at chosen parameters | a fit | [F](#f-fixed-parameter-evaluation) |
| Regional timing | Did syntenic blocks rediploidize together, and when? | + collinear blocks | [G](#g-regional-block-coupled-resolution-timing) |
| Simulation | Calibration, power, null distributions | species tree | [H](#h-simulation-and-calibration) |

**Is synteny usable without block structure?** Yes. Synteny (section E)
needs only the per-family `.syn` files. Each same-genome gene pair is scored
as collinear or not, and that enters every family's likelihood through
σ_W and σ_N. Block IDs are used **only** by the optional, post hoc
regional-timing analysis (section G), which runs on top of kalerax output.
The `.syn` files can come from any collinearity tool (MCScanX, i-ADHoRe,
WGDI, JCVI) via `scripts/make_synteny_files.py`.

---

## Common inputs and settings

- **Families file** (`-f`): an AleRax families file listing each family's
  gene-tree sample (for example 1000 UFBoot trees) and gene→species mapping.
  `scripts/build_family_file.py trees_dir mappings_dir families.txt` builds
  one. Use absolute paths, or run from the directory the relative paths refer
  to.
- **Species tree** (`-s`): rooted Newick with internal node labels. A WGD is
  named by one taxon (terminal branch) or two taxa (the branch to their LCA).
- **Typical settings** for WGD work:

```
COMMON="--rec-model UndatedDTL --transfer-constraint PARENTS \
        --gene-tree-rooting UNIFORM --species-tree-search SKIP \
        --summary-only --seed 1"
```

  `--rec-model UndatedDL` is faster and appropriate when transfers are
  implausible. `--species-tree-search SKIP` keeps the topology fixed, which a
  WGD analysis requires.
- **Output volume.** kalerax writes several files per family.
  `--summary-only` keeps the summaries, parameters and logs.
  `--gene-tree-samples 0` additionally skips reconciliation sampling; that is
  fast, but then the posterior outputs listed under D (`wgdSummary.txt`,
  resolution profiles, `syntenySummary.txt`) are **not** written.
- **Reading results.** The log's "Final species tree likelihood" has only 6
  significant figures. Sum `per_fam_likelihoods.txt` for the exact lnL.
  `model_parameters/wgd_parameters.txt` holds the fitted rates, q, r, σ_W and
  σ_N at full precision.
- **Parallel runs.** Many single-process runs in parallel work well. Give
  every run its own output directory; reusing one is unsafe.

---

## A. Reconciliation (no WGD)

The baseline AleRax analysis, and the null model for every WGD test.

```
kalerax -f families.txt -s species.nw $COMMON --gene-tree-samples 0 -p out_null
```

Outputs: `per_fam_likelihoods.txt`, `model_parameters/model_parameters.txt`
(D, L, T) and `alerax.log`.

## B. WGD test (AORe)

Did a WGD occur at the top of branch *w*? Under AORe the extra copy is
retained with probability q and the ohnologs diverge immediately.

```
kalerax -f families.txt -s species.nw $COMMON --gene-tree-samples 0 \
        --wgd Glycine_max,Glycine_falcata -p out_aore
```

- `--wgd LABEL[,LABEL2][:q0]` is repeatable; `:q0` sets the starting q
  (default 0.2).
- **Test:** 2ΔlnL = 2(lnL_AORe − lnL_null), against ½χ²₀ + ½χ²₁ (q = 0 is a
  boundary). The 5% critical value is 2.71.
- Gene trees alone cannot distinguish a WGD from a burst of small-scale
  duplications on the same branch. To make that distinction, add synteny
  (E).

## C. WGD placement scan

Which branch carries the WGD? Fit B once per candidate branch and compare the
lnL values; all models have the same number of parameters. Adjacent branches
often both score well when the true event lies between them or the signal is
spread. A nested two-WGD fit (two `--wgd` flags) tests whether both are
needed.

## D. Delayed rediploidization (LORe)

Were the ohnologs still tetrasomic at the first speciation after the WGD?
r is the per-branch probability that an unresolved pair rediploidizes. r = 1
is AORe; small r means resolution mostly in descendant lineages (gene trees
then group ohnologs by species).

```
kalerax -f families.txt -s species.nw $COMMON --gene-tree-samples 100 \
        --wgd Polyodon_spathula,Acipenser_ruthenus --lore -p out_lore
```

- `--lore` fits r for every declared internal-branch WGD; their subtrees must
  be disjoint. To test one WGD nested inside others, use
  `--lore-wgd LABEL[,LABEL2]`: only that WGD gets a free r and the others are
  held at AORe.
- **Test of delay:** 2(lnL_LORe − lnL_AORe), against ½χ²₀ + ½χ²₁ (r = 1 is a
  boundary).
- **Where resolution happened** (requires `--gene-tree-samples` > 0):
  - `reconciliations/totalSpeciesResolutionCounts.txt`: expected U→R
    resolutions per species branch.
  - `reconciliations/wgdSummary.txt`: per WGD, q, r and the expected numbers
    of commits and of tetrasomic pairs.
  - Per family: `scripts/summarize_rediploidization.py --summaries-dir
    out_lore/reconciliations/summaries --wgd-branch N2 --out-prefix lore`.
- r is a single genome-wide value. To ask whether timing varies by region,
  use G.

## E. Adding synteny

Gene order breaks two ambiguities that gene trees can't resolve: WGD versus
an SSD burst, and LORe versus independent SSDs in each lineage. Each
same-genome gene pair of a family is observed as collinear (1), not
collinear (0) or undetermined (NA). Pairs separated by a WGD split are
collinear with probability σ_W, all other pairs with σ_N. Both are global
and fitted.

**1. Build the `.syn` files** from within-genome collinearity, one tool output
per genome:

```
python scripts/make_synteny_files.py --families families.txt \
       --mcscanx Glycine_max.collinearity --out syn/ \
       [--placed-genes placed.txt] [--na-genes collapsed.txt] [--rename ids.tsv]
```

Format of `<family>.syn`: `geneA geneB 1|0|NA` per same-genome pair, or
`gene NA`. Gene names must match the gene-tree leaf labels. **Do not filter
anchor pairs by a Ks window around the WGD peak.** Ohnologs in LORe regions
have lower Ks, so such a filter removes exactly the LORe signal.

**2. Fit with synteny.** Add the same flags to the null, AORe and LORe fits:

```
SYN="--synteny-dir syn/ --synteny-min-sigmaW 0.3"
kalerax ... $COMMON $SYN -p out_null_syn
kalerax ... $COMMON $SYN --wgd A,B -p out_aore_syn
kalerax ... $COMMON $SYN --wgd A,B --lore -p out_lore_syn
```

- `--synteny-min-sigmaW 0.3` requires a declared WGD to leave at least 30%
  of its surviving ohnolog pairs collinear. Without it, a WGD model fitted to
  an SSD burst can set σ_W ≈ σ_N and win on gene trees alone.
- `--synteny-params sW,sN` sets starting values; `--fix-synteny` holds them.
  `--synteny-unlisted-na` treats unlisted pairs as NA instead of 0.
- Families without a `.syn` file contribute their plain gene-tree
  likelihood. Synteny needs no block IDs.

**3. Test the WGD.** With synteny, the χ² reference is invalid: σ_W is not
identified under the null, and pairs from one block are not independent.
Calibrate by parametric bootstrap under the fitted null:

```
python scripts/synteny_bootstrap.py --kalerax build/bin/kalerax \
       --species-tree species.nw --null-run out_null_syn --wgd-run out_aore_syn \
       --wgd A,B --min-sigmaW 0.3 --reps 100 --workers 16 --out boot \
       [--block-frac 0.5] [--na-frac 0.1] [--lore] --extra <kalerax flags>
```

`--block-frac` makes simulated collinearity correlated within blocks, a
conservative null.

**4. Test LORe against AORe with synteny.** Add `--null-model aore`: data are
simulated under the fitted AORe + synteny model (q-hat, r = 1, sigma_W-hat,
sigma_N-hat), and each replicate is fitted with AORe (null) and LORe
(alternative). Pass the observed AORe fit as `--null-run` and the observed
LORe fit as `--wgd-run`.

## F. Fixed-parameter evaluation

This computes every family's likelihood at chosen parameter values, with no
optimisation. Uses: likelihood profiles, re-scoring a fit and the per-family
grids for G. Fixed values are read back from a fit's
`model_parameters/wgd_parameters.txt`.

```
kalerax ... $COMMON --gene-tree-samples 0 $SYN \
        --fix-rates --d D --l L --t T \
        --wgd A,B:q --fix-wgd-q \
        --fix-synteny --synteny-params sW,sN \
        --lore --lore-fixed-r r --skip-thorough-rates -p out_r
```

- `--fix-rates` does not fix q; `--fix-wgd-q` holds q at the `:q` value.
  `--lore-fixed-r` holds r on every resolvable WGD (r = 1 is exactly AORe).
- Any subset can be fixed. For example, `--lore-fixed-r r` alone gives the
  profile likelihood of r, with everything else re-fitted.
- Re-evaluating a fit at its own estimates reproduces its lnL
  (`tests/lore_fixed_r_regression.py`). One evaluation costs about the CCP
  build time of the data set.

## G. Regional (block-coupled) resolution timing

Did families in the same collinear block rediploidize together, and did each
block resolve before or after the first speciation below the WGD? This is a
mixture model over per-family likelihoods from F. Block IDs come from the
same collinearity output used for the `.syn` files.

**1. Per-family likelihood grid.** Fit LORe with synteny (E). Then run F at
a grid of r values that includes r̂ and is dense around it, for example
`0 0.1 0.2 0.3 0.35 0.4 r̂ 0.5 0.6 0.7 0.75 0.8 0.9 1`; if r̂ ≈ 1, use
`… 0.9 0.95 0.98 0.99 r̂ 1`. Optionally run a q grid at r = r̂ (for example
`0.85 0.9 q̂ 0.98 0.995`) as a control. Keep each run's
`per_fam_likelihoods.txt` and delete the rest of the run directory.

**2. Block maps and gene order** from MCScanX, with the reference genome
first:

```
python scripts/mcscanx_blocks.py --syn-dir syn/ \
       --collinearity Polyodon_spathula.collinearity Huso_huso.collinearity ... \
       --out blocks
# -> blocks_fam2block_{alignment,chrpair,fallback}.tsv, blocks_order_alignment.tsv,
#    blocks_coverage.tsv; prints the fraction of .syn collinear pairs found
#    (should be ~100%)
```

A family's block is the reference-genome MCScanX alignment containing its
collinear pair(s). The `chrpair` map is coarser (reference chromosome pair).
The `fallback` map assigns otherwise-unblocked families to other genomes'
alignments; it gives broader coverage but is open to ascertainment effects.

**3. Mixture models, coupling test and HMM.**

```
python scripts/block_lore.py \
       --grid 0=out_r0 0.1=out_r0.1 ... 1=out_r1 --rhat 0.4426=out_rhat \
       --blocks blocks_fam2block_alignment.tsv --order blocks_order_alignment.tsv \
       --perm 1000 --out aln
# q control: pass the q-grid runs as --grid q=dir (labels are the q values), no --order
```

Outputs:
- `aln_summary.json`:
  - `single_*`: one genome-wide r.
  - `indep`, `unblocked`: independent mixtures.
  - `block.gain_vs_indep`: block-shared minus independent, on blocked
    families.
  - `block.perm`: permutation null (mean, sd, max, p).
  - `hmm`: ρ, mean run length, and its gain over independent on the same
    families.
- `aln_blocks.tsv`: per block, the number of families, E[r | block] (the
  probability of resolution before the first speciation below the WGD), the
  MAP class, and its lnL contribution.
- `aln_families.tsv`: per-family E[r] under each model.

**4. Calibrate the HMM** (optional). The HMM gain over the independent
mixture is calibrated by permuting families among the gene-ordered positions:

```
python scripts/hmm_permutation.py --grid 0=out_r0 ... --order blocks_order_alignment.tsv \
       --perm 50 --perm-start 0 --seed 1 --out hmm_shard0     # run shards in parallel
python scripts/hmm_permutation.py --merge hmm_shard*_hmm_perm.json --out hmm
```

**How to read it:**
- **Coupling:** evidence for coupling is `block.gain_vs_indep` far above
  the permutation maximum, with a small p. The raw gain can be negative when
  timing is not shared; only the permutation-calibrated comparison counts.
- **Two populations:** blocked and unblocked families get separate class
  weights, because block membership depends on the same synteny data.
- **Use block timings only after a positive test.** When blocks are not
  coupled, forcing shared timing makes per-family estimates worse.
- **The q control** shows whether block structure is in timing or in
  retention. It has low power when the q likelihoods are flat.
- **Caveats:**
  - The permutation null treats families as exchangeable, so any regionally
    autocorrelated property that affects L_f(r) (local loss, assembly
    collapse, rearrangement) also registers as coupling.
  - Rates, q and σ are held at the genome-wide fit.
  - E[r] values are grid classes, so they give relative timing.

Worked results: Acipenseriformes shows strong coupling (+1,030 lnL,
p ≤ 0.001); *Glycine* is a negative control (no coupling, r̂ ≈ 1).

## H. Simulation and calibration

- `scripts/kalerax_sim.py` is a Python module implementing the model's
  forward process: DL/DTL, WGD q, LORe r, emission-level synteny with
  optional NA and within-block correlation. Main entry points:
  `SpeciesTree`, `Params`, `simulate_families`, `simulate_block_families`
  (block-shared resolution classes with coupling κ), `write_family_set` and
  `synteny_lines`. Simulations with a fixed seed are reproducible.
- `scripts/block_lore_sim.py <kalerax> <species_tree> <workdir> <seed> <κ>
  [n_blocks] [mean_block_size] [n_perm]` runs one end-to-end replicate of G
  (simulate, free fit, grid, mixtures, permutations, HMM) and writes
  `result.json`. Run many seeds at κ = 0 to check calibration and at κ > 0
  for power.
- `scripts/synteny_bootstrap.py` provides the null distribution for the
  synteny WGD test (E3).
- Regression tests are in `tests/`. `lore_fixed_r_regression.py` and
  `synteny_bound_regression.py` take `<kalerax> <workdir>`; the shell tests
  take the binary path and need a writable `mktemp`.

## I. Rate heterogeneity among families (gamma)

`--speciation-gamma-categories K` (inherited from AleRax) gives each family one
of K equiprobable rate categories. In category c the speciation weight is the
gamma scaler s_c (mean 1, shape alpha) and D, L, T, S are normalised jointly,
so D, L and T are scaled together by 1/s_c; the D:L:T ratios are fixed. The
per-family likelihood is the mean over categories, each conditioned on survival.

Without it, LORe can absorb rate heterogeneity: a WGD with r -> 0 and q > 0
is a duplication/transfer-free lineage class, i.e. a two-class rate mixture.
Fit gamma before interpreting r or q.

- The fitted alpha is written to `model_parameters/wgd_parameters.txt`
  (`gamma alpha NA <alpha> <K>`).
- `--speciation-gamma-alpha <a>` sets the starting alpha; with
  `--fix-gamma-alpha` it is held. Fixed-parameter grids (section F) under gamma
  add `--speciation-gamma-categories K --speciation-gamma-alpha <a>
  --fix-gamma-alpha`; this re-evaluates a free fit exactly.
- Optimisation: alpha and the overall rate lie on a strong ridge. kalerax first
  moves alpha with D, L, T rescaled to keep the mean event intensity
  E_c[1/s_c] fixed, then refines alpha alone (Brent tolerance 1e-3; upstream
  0.1). `--converge-rounds n` repeats the thorough pass until lnL gains < 0.1.
- Simulation: `kalerax_sim.simulate_families(..., gamma_alpha=a,
  gamma_categories=K)` draws one category per family and redoes rejections
  within it. `synteny_bootstrap.py` picks up gamma from the `--null-run`
  parameter file automatically (`--no-gamma` to disable); `block_lore_sim.py`
  takes `gamma_alpha` as its 9th argument.

## Practical notes

- **Runtime:** a full LORe + synteny fit of about 18k families with 1000
  trees each takes about 1–1.5 h on one core. A fixed-parameter evaluation
  takes 10–35 min, mostly CCP construction.
- **Invalid arguments** exit with status 1.
- **MPI** is supported, as in AleRax, but the WGD/LORe/synteny additions have
  mostly been run single-process.
- **Upstream AleRax options** (highways, per-species rates, species-tree
  search) still work, but WGD analyses assume a fixed species tree.
