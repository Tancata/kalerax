# Independent check of the DTL duplication weighting with Zombi

This note asks whether kalerax's UndatedDTL recovers the true numbers of
duplications and transfers under each duplication weighting, using gene
families simulated by an independent forward simulator. The simulator is Zombi
(Davín et al. 2020, *Bioinformatics* 36:1286; github.com/AADavin/Zombi @ 8db13ee).

## Setup (`validation/zombi_check.py`)

Run as `ZOMBI_DIR=/path/to/Zombi KALERAX=build/bin/kalerax python3 validation/zombi_check.py OUTDIR TAG D T L [K M]`
(needs ete3, networkx, pyvolve and numpy for Zombi).


- **Species tree:** a Zombi `T`-mode pure-birth tree (no extinction, so no
  transfers from unsampled lineages) with 12 extant species.
- **Gene families:** simulated in Zombi `Gm` mode as 14 independent runs of 300
  root families. Fixed per-family rates of D, T and L per unit time; additive
  transfers to a uniformly chosen contemporaneous lineage; no inversions,
  transpositions or originations.
- **Truth:** the D and T splits at internal nodes of Zombi's pruned (extant)
  gene trees, i.e. the observable events.
- **Fits:** kalerax UndatedDTL with PARENTS, UNIFORM gene-tree rooting and ROOT
  origination, `--summary-only -g 10`, run with the default weighting and with
  `--ale-dup-weights`. Exported duplication and transfer totals
  (`totalSpeciesEventCounts.txt`) are compared with the truth over the families
  kalerax kept (≥4 species).

Zombi's process is dated and its transfers go only to contemporaneous lineages,
whereas kalerax's undated model is an approximation. Some residual error is
therefore expected under both weightings. The question is which weighting
recovers the D vs T balance.

## Results

| Zombi config (D/T/L per time unit) | families kept | weighting | exported D vs true | exported T vs true |
|---|---|---|---|---|
| lowT (0.8/0.2/1.0) | 2196 | ALE (--ale-dup-weights) | -11.2% (15287 vs 17220) | +42.2% (6085 vs 4278) |
| lowT (0.8/0.2/1.0) | 2196 | forward-process (default) | -2.2% (16849 vs 17220) | +13.4% (4852 vs 4278) |
| midT (0.8/0.5/1.0) | 2598 | ALE (--ale-dup-weights) | -15.5% (23064 vs 27310) | +25.6% (20953 vs 16687) |
| midT (0.8/0.5/1.0) | 2598 | forward-process (default) | -1.9% (26798 vs 27310) | +8.0% (18024 vs 16687) |
| highT (0.8/1.0/1.0) | 2981 | ALE (--ale-dup-weights) | -17.2% (48123 vs 58118) | +13.4% (81517 vs 71858) |
| highT (0.8/1.0/1.0) | 2981 | forward-process (default) | +2.0% (59260 vs 58118) | +2.0% (73285 vs 71858) |
| lowD (0.3/0.5/0.6) | 3084 | ALE (--ale-dup-weights) | -16.8% (8696 vs 10447) | +14.5% (19299 vs 16858) |
| lowD (0.3/0.5/0.6) | 3084 | forward-process (default) | +3.6% (10825 vs 10447) | +2.6% (17300 vs 16858) |

The ALE/AleRax weighting under-counts duplications (by 11–17%) and over-counts
transfers (by 13–42%) in every configuration, as expected when D splits receive
half the weight of T splits. The default (forward-process) weighting removes
most of this: D is within −2% to +4%, and T is within +2% to +13%. The remaining
T excess, largest at low transfer rates, is plausibly the dated-vs-undated
mismatch. Likelihood values differ between the two weightings by construction (different
per-tree constants), so they are not a model-comparison criterion. Raw numbers are in `validation/zombi_dup_weights.csv`.
