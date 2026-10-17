

# kalerax — a WGD/LORe-focused fork of AleRax

> **This is a fork of [AleRax](https://github.com/BenoitMorel/AleRax)** (the
> binary is named `kalerax`, to keep it distinct from upstream AleRax). It adds
> support for modelling **whole-genome duplications (WGD)** in AleRax's undated
> reconciliation models: a WGD can be declared on any species branch with a
> per-branch retention parameter `q`, which is then estimated by maximum
> likelihood — letting you test WGD hypotheses the way
> [WHALE](https://github.com/arzwa/Whale.jl) does, but on AleRax's amalgamated
> DL and DTL likelihoods. It also models **lineage-specific rediploidization
> (LORe)** via `--lore` — delayed ohnolog divergence with a fitted resolution
> probability `r` **per declared WGD** (disjoint WGD subtrees; a WGD on a
> terminal branch has its `r` pinned to 1/AORe) — for **both** the DL and the
> DTL model. `--lore-wgd LABEL[,LABEL2]` instead tests LORe for one specific WGD
> (which may be **nested inside** other, AORe-pinned WGDs), and `--summary-only`
> suppresses the bulky per-family output (keeping the summaries, `wgdSummary.txt`
> and logs) for large runs.
>
> **DTL weighting.** Under `UndatedDTL`, duplication and WGD splits are weighted
> as in the forward process relative to transfers (default since October 2026).
> `--ale-dup-weights` restores the ALE/AleRax numerics. DTL results from earlier
> versions are not comparable; see [`DTL_LORE.md`](DTL_LORE.md).
>
> **Synteny (`--synteny-dir`).** Per-family collinearity observations for
> same-genome gene pairs are integrated into the reconciliation likelihood, so
> WGD can be separated from bursts of small-scale duplication and LORe from
> parallel SSDs — see [`SYNTENY.md`](SYNTENY.md) and
> `scripts/make_synteny_files.py`.
>
> **Event export.** Sampled reconciliations are drawn from the fitted WGD/LORe
> model: WGD ohnolog splits are counted in a separate `wgd_duplications` column
> of the per-species event counts (the `duplications` column holds SSDs only)
> and are labelled `W` / `:WGD=Y` in the sampled gene trees.
>
> 📄 **Start here: [`WGD_REPORT.md`](WGD_REPORT.md)** — a short, readable report
> covering the implementation, a benchmark against WHALE, and what allowing
> horizontal transfers (DL vs DTL) does to the evidence for WGD;
> [`DTL_LORE.md`](DTL_LORE.md) documents the DTL+LORe model;
> [`SYNTENY.md`](SYNTENY.md) documents the synteny emission.
>
> The WGD/LORe code lives mainly in `src/ale/UndatedDLMultiModel.hpp` and
> `src/ale/UndatedDTLMultiModel.hpp`, exposed via the `--wgd` and `--lore`
> command-line options; cross-validation inputs and recipes are under
> `validation/`. This is a
> research prototype — see the report's *Limitations*. Everything below is the
> upstream AleRax documentation.
>
> **Inspiration & acknowledgement.** The WGD retention model and the validation
> approach are directly inspired by **WHALE**
> ([arzwa/Whale.jl](https://github.com/arzwa/Whale.jl)), which implements a
> duplication–loss–WGD (DLWGD) model on the same amalgamated-likelihood
> framework. If you use the WGD features here, please also cite the WHALE paper:
>
> > Zwaenepoel, A. and Van de Peer, Y., 2019. *Inference of Ancient Whole-Genome
> > Duplications and the Evolution of Gene Duplication and Loss Rates.*
> > Molecular Biology and Evolution, 36(7), pp. 1384–1404.
> > [doi:10.1093/molbev/msz088](https://doi.org/10.1093/molbev/msz088)

---

# AleRax  

AleRax is a parallel tool for species tree - gene tree inference and reconciliation under gene duplication, loss, and HGT. For each gene family, it takes as input a gene tree distribution (typically inferred with Bayesian inference tools such as MrBayes, PhyloBayes, etc.). AleRax can perform the following operations:
* Species tree inference
* Species tree rooting 
* Reconciled gene tree sampling 
* Model parameter estimation (e.g. DTL event probabilities) 
* Statistical test of different species tree hypotheses (you'll need to instal consel)
  
We are also working on the following features:
* Relative order of speciation event (relative dating) from HGT constraints
* Inference of highways of transfers (pairs of species that exchanged many genes via HGT)

When using AleRax, please cite: [https://academic.oup.com/bioinformatics/article/40/4/btae162/7633408](https://academic.oup.com/bioinformatics/article/40/4/btae162/7633408)

## Requirement

* A Linux or MacOS environnement
* gcc 5.0 or > 
* CMake 3.6 or >
* MPI (required if you want to use parallelization)

## Installation 


To download AleRax, please use git,  and clone with --recursive!!!

```
git clone --recursive https://github.com/BenoitMorel/AleRax
```

For some reason, cloning coraxlib (one of the dependencies) sometimes fails:
```
Cloning into '/home/benoit/github/AleRax/ext/GeneRaxCore/ext/coraxlib'...
fatal: unable to access 'https://codeberg.org/Exelixis-Lab/coraxlib.git/': server certificate verification failed. CAfile: none CRLfile: none
```
You can skip the SSL certificate verification with:
```
git config --global http.sslVerify "false"
```
You can also have a look at this [post](https://forum.gitlab.com/t/gitlab-runner-server-certificate-verification-failed/59450/8) for a more detailed explanation of the problem.


To build the sources:
```
./install.sh
```

The generated executable is located here:
```
build/bin/kalerax
```

To copy the executable to your PATH, such that you can call kalerax from anywhere:
```
cd build
sudo make install
```

## Updating to the last stable version

In case we've made some changes since the last time you updated your repository, you can get and install the most recent version of AleRax with:

```
./gitpull.sh # does a bit more than git pull
./install.sh
```

Note that `git pull` does not update the submodules used by AleRax, hence the need for `./gitpull.sh`

## Running

See the wiki (https://github.com/BenoitMorel/AleRax/wiki)

## Issues and questions

If you encounter any issue, please report it!! I'm always happy to help.
For questions, issues or feedback, please post on the GeneRax (even if it's AleRax) google group: https://groups.google.com/g/generaxusers.
When reporting an issue, please send me at least the command line you ran, the logs file, and any file that might be relevant (e.g. the species tree file if AleRax failed to read the species tree...)


