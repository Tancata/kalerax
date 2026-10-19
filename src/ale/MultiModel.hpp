#pragma once
#include <cstdlib>

#include <IO/HighwayCandidateParser.hpp>
#include <ccp/ConditionalClades.hpp>
#include <likelihoods/reconciliation_models/BaseReconciliationModel.hpp>

#include "Synteny.hpp"

#include <cmath>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

/**
 *  Temp storage for a clade event sampled in
 *  the reconciliation space
 */
template <class REAL> struct ReconciliationCell {
  Scenario::Event event;
  REAL maxProba;
  double blLeft;
  double blRight;
};

/**
 *  Base class for reconciliation models that take as input
 *  gene tree distributions (represented by conditional clade
 *  probabilities)
 */
// multiply a CLV term by 2^(-256*s) (synteny split-weight scaler, see
// MultiModelInterface::_synNScale); exact for ScaledValue, plain for double
//
// ScaledValue::scale() renormalises by a single 2^256 step and ScaledValue
// addition keeps only the operand with the smaller scaler, which is correct
// only for normalised mantissas in [2^-256, 1]. A synteny-weighted term (CLV *
// CLV * rate * weight, weight down to e^-177) can sit several steps below that
// after the single scale() at the call site; left unnormalised it would win
// additions against genuinely larger terms and silently drop them. So first
// normalise fully (extra scale() calls are no-ops on normalised values), then
// apply the 2^-256 factors one step at a time.
template <class REAL> inline void applySynScale(REAL &t, int s) {
  for (int i = 0; i < 4; ++i) {
    scale(t);
  }
  for (int i = 0; i < s; ++i) {
    t *= JS_SCALE_THRESHOLD;
    scale(t);
  }
}

class MultiModelInterface : public BaseReconciliationModel {
public:
  MultiModelInterface(PLLRootedTree &speciesTree,
                      const GeneSpeciesMapping &geneSpeciesMapping,
                      const RecModelInfo &info, const std::string &ccpFile)
      : BaseReconciliationModel(speciesTree, geneSpeciesMapping, info),
        _memorySavings(info.memorySavings) {
    _ccp.unserialize(ccpFile);
    mapGenesToSpecies();
  }

  virtual ~MultiModelInterface() {}

  const ConditionalClades &getCCP() const { return _ccp; }

  virtual void setAlpha(double alpha) = 0;
  virtual void setRates(const RatesVector &rates) = 0;
  virtual void setHighways(const std::vector<Highway> &highways) {
    (void)(highways);
  }

  // Declare a whole-genome duplication at the top of a species branch, with
  // retention probability q. No-op for models that do not implement WGD.
  virtual void setWGD(unsigned int /*speciesNode*/, double /*q*/) {}

  // Weight applied to every duplication split (SSD D event and the WGD/LORe
  // doubling-bracket split) in models that also have transfers. 1 (default)
  // = the ALE/AleRax undated convention, which counts an unordered D split
  // once while a T split sums both stay/transfer orderings; 2 = weights that
  // match the forward process (two exchangeable copies) relative to transfers
  // (--ordered-dup-weights). No-op for models without transfers (UndatedDL),
  // where the factor is a per-gene-tree constant.
  virtual void setDupSplitWeight(double /*w*/) {}

  // Set the (global) LORe resolution probability r in [0,1]. r=1 recovers the
  // immediate-resolution (WHALE/AORe) WGD model. No-op for models without LORe.
  virtual void setResolutionProb(double /*r*/) {}

  // Set a *per-species-branch* LORe resolution vector (index == species
  // node_index). Used for per-event r: each declared WGD paints its own subtree
  // with its own r, every other branch stays at r=1 (AORe). No-op without LORe.
  virtual void setResolutionProbVector(const std::vector<double> & /*perBranchR*/) {}

  // Sample `samples` resolution histories and accumulate, per species branch,
  // the number of U->R commit (WGD resolution) events. CLVs must be current.
  // No-op (empty output) for models without LORe. See WGD_LORE_marginal.md.
  virtual void sampleResolutionCommits(unsigned int /*samples*/,
                                       std::vector<double> &commitCounts,
                                       std::vector<double> &tetraCounts) {
    commitCounts.clear();
    tetraCounts.clear();
  }


  // --- Synteny (Layer 1) -------------------------------------------------
  // Observation per same-genome gene pair (g,h): y=1 (anchor pair of a
  // collinear duplicated block), y=0 (both placed, not collinear), or NA.
  // Every same-genome pair has exactly one LCA in a reconciled gene tree; the
  // pair's emission is Bernoulli(sigma_W) if that LCA is a WGD split (AORe
  // doubling or LORe U->R resolution) and Bernoulli(sigma_N) otherwise. Hence
  //   P(Y | rec) = C(sigma_N) * prod_{WGD splits} rho(split),
  //   C = prod_pairs f_N(y),  rho = prod_{pairs across split} f_W(y)/f_N(y),
  // Numerically, each pair's emission is factored by m = max(f_W, f_N), so
  //   P(Y | rec) = Cm * prod_{non-WGD splits} wN(split) * prod_{WGD} wW(split)
  //   Cm = prod_pairs m(y),  wX(split) = prod_{pairs across} f_X(y)/m(y) <= 1.
  // wN is folded into the CCP split frequencies (all S/D/T/U-speciation
  // terms), wW replaces the frequency in the WGD doubling-bracket split term,
  // and log Cm is added to the family log-likelihood. All weights are <= 1,
  // so every CLV stays a (sub-)probability. No synteny => weights 1, Cm = 1.
  using SyntenyInput = ::SyntenyInput;
  // Count, for every clade split, the same-species cross pairs with y=1 / y=0.
  void setSynteny(const SyntenyInput &in);
  // Recompute rho and log C for the given emission parameters.
  void setSyntenyParams(double sigmaW, double sigmaN);
  bool hasSynteny() const { return _synActive; }
  // total observed pairs (y=1, y=0) in this family, for logging
  void getSyntenyTotals(unsigned int &n1, unsigned int &n0) const {
    n1 = _synTot1;
    n0 = _synTot0;
  }

  virtual double computeLogLikelihood() = 0;

  virtual bool inferMLScenario(Scenario &scenario) = 0;

  virtual bool
  sampleReconciliations(unsigned int samples,
                        std::vector<std::shared_ptr<Scenario>> &scenarios) = 0;

protected:
  // should be called on changing the model rates (alpha/DTLO/
  // highways), as the LLs computed earlier become irrelevant
  void resetCache() { _llCache.clear(); }

private:
  virtual void mapGenesToSpecies();

protected:
  const bool _memorySavings;
  ConditionalClades _ccp;
  std::unordered_map<size_t, double> _llCache;

  // --- Synteny (Layer 1) state ---
  bool _synActive = false;
  std::vector<std::vector<unsigned int>> _synN1; // per cid, per split
  std::vector<std::vector<unsigned int>> _synN0;
  std::vector<std::vector<double>> _synOrigFreq; // unweighted ccp frequencies
  std::vector<std::vector<double>> _synWFreq;    // WGD-split frequency * wW
  double _synLogConst = 0.0;                     // log Cm
  unsigned int _synTot1 = 0;
  unsigned int _synTot0 = 0;
  // frequency to use for split k of clade cid in the WGD doubling bracket
  double synWgdFreq(CID cid, size_t k, double freq) const {
    return _synActive ? _synWFreq[cid][k] : freq;
  }
  // Synteny split weights can be far below the double range (a split carrying
  // hundreds of collinear pairs that must be explained by a non-WGD event has
  // log-weight ~ n1*log(sigmaN/sigmaW)). Such weights are stored as an in-range
  // mantissa (in the split frequency / _synWFreq) times 2^(-256*s); s is the
  // integer stored here and applied to the CLV term with applySynScale, which
  // is exact for ScaledValue. In double precision the term underflows to 0 as
  // before, which triggers the high-precision rebuild of the family.
  std::vector<std::vector<int>> _synNScale; // per cid, per split (non-WGD)
  std::vector<std::vector<int>> _synWScale; // per cid, per split (WGD)
  int synScaleN(CID cid, size_t k) const {
    return _synNScale.empty() ? 0 : _synNScale[cid][k];
  }
  int synScaleN(CID cid, const CladeSplit &cs) const {
    return _synNScale.empty()
               ? 0
               : _synNScale[cid][&cs - _ccp.getCladeSplits(cid).data()];
  }
  int synScaleW(CID cid, size_t k) const {
    return _synWScale.empty() ? 0 : _synWScale[cid][k];
  }
  void restoreSyntenyFrequencies() {
    for (size_t cid = 0; cid < _synOrigFreq.size(); ++cid) {
      auto &splits = _ccp.getCladeSplitsMutable(cid);
      for (size_t k = 0; k < splits.size(); ++k) {
        splits[k].frequency = _synOrigFreq[cid][k];
      }
    }
    _synNScale.clear();
    _synWScale.clear();
  }
};

inline void MultiModelInterface::mapGenesToSpecies() {
  // fill _geneToSpecies, _speciesCoverage and _numberOfCoveredSpecies
  this->_geneToSpecies.clear();
  this->_speciesCoverage.resize(this->_speciesTree.getLeafNumber(), 0);
  this->_numberOfCoveredSpecies = 0;
  const auto &cidToLeaves = _ccp.getCidToLeaves();
  for (const auto &p : cidToLeaves) {
    const auto cid = p.first;
    const auto &geneName = p.second;
    const auto &speciesName = this->_geneNameToSpeciesName[geneName];
    const auto spid = this->_speciesNameToId[speciesName];
    this->_geneToSpecies[cid] = spid;
    if (!this->_speciesCoverage[spid]) {
      this->_numberOfCoveredSpecies++;
    }
    this->_speciesCoverage[spid]++;
  }
  // build the pruned species tree representation based
  // on the species coverage
  onSpeciesTreeChange(nullptr);
}


inline void MultiModelInterface::setSynteny(const SyntenyInput &in) {
  restoreSyntenyFrequencies(); // idempotent if synteny was set before
  const auto nClades = _ccp.getCladesNumber();
  // leaves per clade, bucketed by species: cid -> {species -> leaf cids}
  std::vector<std::unordered_map<unsigned int, std::vector<CID>>> bySpecies(
      nClades);
  const auto &cidToLeaves = _ccp.getCidToLeaves();
  for (CID cid = 0; cid < nClades; ++cid) {
    if (_ccp.isLeaf(cid)) {
      bySpecies[cid][this->_geneToSpecies[cid]].push_back(cid);
    } else {
      const auto &splits = _ccp.getCladeSplits(cid);
      assert(!splits.empty());
      const auto &cs = splits[0]; // any split partitions the same leaf set
      for (auto *side : {&bySpecies[cs.left], &bySpecies[cs.right]}) {
        for (const auto &kv : *side) {
          auto &dst = bySpecies[cid][kv.first];
          dst.insert(dst.end(), kv.second.begin(), kv.second.end());
        }
      }
    }
  }
  auto obs = [&](CID a, CID b) -> int { // 1, 0, or -1 (NA)
    const auto &ga = cidToLeaves.at(a);
    const auto &gb = cidToLeaves.at(b);
    if (in.naGenes.count(ga) || in.naGenes.count(gb)) {
      return -1;
    }
    auto key = syntenyKey(ga, gb);
    if (in.naPairs.count(key)) {
      return -1;
    }
    auto it = in.pairObs.find(key);
    if (it != in.pairObs.end()) {
      return it->second;
    }
    return in.unlistedNA ? -1 : 0;
  };
  _synN1.assign(nClades, {});
  _synN0.assign(nClades, {});
  _synTot1 = 0;
  _synTot0 = 0;
  for (CID cid = 0; cid < nClades; ++cid) {
    const auto &splits = _ccp.getCladeSplits(cid);
    _synN1[cid].assign(splits.size(), 0);
    _synN0[cid].assign(splits.size(), 0);
    for (size_t k = 0; k < splits.size(); ++k) {
      const auto &L = bySpecies[splits[k].left];
      const auto &R = bySpecies[splits[k].right];
      for (const auto &kv : L) {
        auto it = R.find(kv.first);
        if (it == R.end()) {
          continue;
        }
        for (auto a : kv.second) {
          for (auto b : it->second) {
            int y = obs(a, b);
            if (y == 1) {
              _synN1[cid][k]++;
            } else if (y == 0) {
              _synN0[cid][k]++;
            }
          }
        }
      }
    }
  }
  // family totals: every same-species leaf pair, counted once (root clade's
  // leaf set is the whole family)
  const CID root = nClades - 1;
  for (const auto &kv : bySpecies[root]) {
    const auto &v = kv.second;
    for (size_t i = 0; i < v.size(); ++i) {
      for (size_t j = i + 1; j < v.size(); ++j) {
        int y = obs(v[i], v[j]);
        if (y == 1) {
          _synTot1++;
        } else if (y == 0) {
          _synTot0++;
        }
      }
    }
  }
  _synActive = (_synTot1 + _synTot0) > 0;
  _synOrigFreq.assign(nClades, {});
  for (CID cid = 0; cid < nClades; ++cid) {
    for (const auto &cs : _ccp.getCladeSplits(cid)) {
      _synOrigFreq[cid].push_back(cs.frequency);
    }
  }
  if (!_synActive) {
    _synN1.clear();
    _synN0.clear();
    _synWFreq.clear();
    _synOrigFreq.clear();
    _synLogConst = 0.0;
  }
}

inline void MultiModelInterface::setSyntenyParams(double sigmaW,
                                                  double sigmaN) {
  if (!_synActive) {
    return;
  }
  const double lw1 = std::log(sigmaW), lw0 = std::log(1.0 - sigmaW);
  const double ln1 = std::log(sigmaN), ln0 = std::log(1.0 - sigmaN);
  const double lm1 = std::max(lw1, ln1), lm0 = std::max(lw0, ln0);
  _synLogConst = _synTot1 * lm1 + _synTot0 * lm0;
  _synWFreq.assign(_synN1.size(), {});
  _synNScale.assign(_synN1.size(), {});
  _synWScale.assign(_synN1.size(), {});
  bool anyScaled = false;
  // split a log-weight lw <= 0 into exp(lw + s*L) * 2^(-256*s), L = 256 ln 2
  // threshold below which a log-weight is carried with an explicit scaler
  // (a weight below ~1e-230 times a ScaledValue product would underflow in
  // the double mantissa); KALERAX_SYN_SCALE_THRESHOLD overrides it (testing
  // only: results must not depend on it)
  static const double thr = [] {
    const char *v = std::getenv("KALERAX_SYN_SCALE_THRESHOLD");
    return v ? std::atof(v) : -150.0;
  }();
  auto splitLog = [&anyScaled](double lw, int &s) {
    const double L = 256.0 * std::log(2.0);
    s = 0;
    if (lw < thr) {
      // largest s with lw + s*L <= 0: mantissa in (e^-L, 1] (a mantissa
      // above 1 could not be renormalised downwards by scale())
      s = static_cast<int>(std::floor(-lw / L));
      anyScaled = true;
    }
    return std::exp(lw + s * L);
  };
  for (size_t cid = 0; cid < _synN1.size(); ++cid) {
    auto &splits = _ccp.getCladeSplitsMutable(cid);
    _synWFreq[cid].resize(splits.size());
    _synNScale[cid].resize(splits.size());
    _synWScale[cid].resize(splits.size());
    for (size_t k = 0; k < splits.size(); ++k) {
      const double n1 = _synN1[cid][k], n0 = _synN0[cid][k];
      const double lwN = n1 * (ln1 - lm1) + n0 * (ln0 - lm0); // <= 0
      const double lwW = n1 * (lw1 - lm1) + n0 * (lw0 - lm0); // <= 0
      splits[k].frequency =
          _synOrigFreq[cid][k] * splitLog(lwN, _synNScale[cid][k]);
      _synWFreq[cid][k] =
          _synOrigFreq[cid][k] * splitLog(lwW, _synWScale[cid][k]);
    }
  }
  if (!anyScaled) {
    _synNScale.clear();
    _synWScale.clear();
  }
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

/**
 *  Implements all the methods that require the REAL template
 *  and are not dependent on modelling HGT
 */
template <class REAL> class MultiModel : public MultiModelInterface {
public:
  MultiModel(PLLRootedTree &speciesTree,
             const GeneSpeciesMapping &geneSpeciesMapping,
             const RecModelInfo &info, const std::string &ccpFile)
      : MultiModelInterface(speciesTree, geneSpeciesMapping, info, ccpFile) {}

  virtual ~MultiModel() {}

  /**
   *  Compute the LL of the current species tree
   */
  virtual double computeLogLikelihood();

  /**
   *  Find the ML reconciled gene tree and store it in
   *  the scenario variable
   */
  virtual bool inferMLScenario(Scenario &scenario);

  /**
   *  Sample reconciled gene trees and add them to
   *  the scenarios vector
   */
  virtual bool
  sampleReconciliations(unsigned int samples,
                        std::vector<std::shared_ptr<Scenario>> &scenarios);

protected:
  // functions to work with CLVs
  virtual void allocateMemory() = 0;
  virtual void deallocateMemory() = 0;
  virtual void updateCLV(CID cid) = 0;

  /**
   *  Probability of an arbitrary gene family to originate
   *  somewhere in the species tree and survive (at least as
   *  a single copy on a single leaf).
   *  All possible reconciliations of the current family are
   *  included by this scenario, thus it is the conditioning
   *  probability
   */
  virtual double getLikelihoodFactor(unsigned int category) = 0;

  /**
   *  Probability of the current family to evolve starting
   *  from a given species branch.
   *  Sum over all scenarios stemming from the species branch
   */
  virtual REAL getRootCladeLikelihood(corax_rnode_t *speciesNode,
                                      unsigned int category) = 0;

  /**
   *  Compute the CLV value for a given cid (clade id) and a given
   *  species node. Return true in case of success.
   *  If recCell is set, sample the next event in the reconciliation
   *  space
   */
  virtual bool
  computeProbability(CID cid, corax_rnode_t *speciesNode, unsigned int category,
                     REAL &proba,
                     ReconciliationCell<REAL> *recCell = nullptr) = 0;


  // --- WGD / LORe-aware scenario sampling (fix for the export backbone) ---
  // A parent's S/SL term consumes the WGD-transformed child value
  //   R_top = (1-q) R + q U,
  // so the sampler must choose R vs U with exactly those weights before
  // sampling the child's events (sampleWGDCoin), and must be able to sample
  // the U recursion itself (sampleUEvent). Models without WGD keep defaults.
  struct UEventSample {
    enum Type {
      RESOLVE_SPLIT,  // U->R with both copies observed: WGD ohnolog split
      RESOLVE_SINGLE, // U->R, one copy lost (2*E*R term): no gene node
      SPECIATE,       // unresolved speciation, both daughters stay U
      SPECIATE_LOSS,  // unresolved speciation, one daughter lost
      LEAF            // still unresolved at an extant tip: seen as 1 gene
    } type = LEAF;
    CID left = 0, right = 0;     // RESOLVE_SPLIT / SPECIATE gene clades
    double blLeft = 0.0, blRight = 0.0;
    bool leftToLeftSpecies = true; // SPECIATE: left clade -> left species
    corax_rnode_t *dest = nullptr; // SPECIATE_LOSS: surviving species
    corax_rnode_t *lost = nullptr; // SPECIATE_LOSS: lost species
  };
  // true => the lineage enters the U state at the top of branch sp
  virtual bool sampleWGDCoin(CID /*cid*/, corax_rnode_t * /*sp*/,
                             unsigned int /*category*/) {
    return false;
  }
  virtual bool sampleUEvent(CID /*cid*/, corax_rnode_t * /*sp*/,
                            unsigned int /*category*/, UEventSample & /*ev*/) {
    return false;
  }

  /**
   *  Return a uniform random REAL from the [0,max) interval
   */
  REAL getRandom(REAL max) {
    auto rand = max * Random::getProba();
    scale(rand);
    if (rand == max) { // may occur due to rounding issues
      rand = REAL();
    }
    return rand;
  }

private:
  /**
   *  Return an undated or dated species tree hash
   */
  virtual size_t getHash() = 0;

  /**
   *  Number of mixture categories for the per-branch
   *  probabilities of elementary events
   */
  unsigned int getGammaCatNumber() { return this->_info.gammaCategories; }

  /**
   *  Fill CLV for each observed gene clade
   */
  void updateCLVs();

  /**
   *  Jointly sample a category and a species branch based on
   *  the probability of the gene tree to be rooted in it
   */
  corax_rnode_t *sampleOriginationSpecies(unsigned int &category);

  /**
   *  Wrapper for the backtrace function
   */
  bool computeScenario(Scenario &scenario, bool stochastic);

  /**
   *  Recursively sample a reconciled gene tree
   */
  bool backtrace(CID cid, corax_rnode_t *speciesNode, corax_unode_t *geneNode,
                 unsigned int category, Scenario &scenario, bool stochastic);
  // descend into a child species branch through its WGD (R/U coin)
  bool descend(CID cid, corax_rnode_t *speciesNode, corax_unode_t *geneNode,
               unsigned int category, Scenario &scenario, bool stochastic);
  // sample the unresolved (U) state at (cid, speciesNode)
  bool backtraceU(CID cid, corax_rnode_t *speciesNode,
                  corax_unode_t *geneNode, unsigned int category,
                  Scenario &scenario, bool stochastic);

  // Backtrace recursion-depth guard. The reconciliation backtrace recurses on
  // the DL/TL "nothing observed, resample" events; at a degenerate cell (where
  // the terminating events have vanishing probability) this recursion can go
  // arbitrarily deep and overflow the stack (observed under DTL+WGD/LORe, where
  // the retention transform suppresses speciation at the WGD branch). The
  // counter is reset per scenario; once it exceeds the cap, backtrace returns
  // false (the scenario is discarded and resampled). A legitimate scenario
  // never approaches the cap.
  size_t _btScenarioSteps = 0;
  static const size_t kBtScenarioCap = 8000; // max recursion depth (kept well
  // below the ~64MB stack's overflow threshold: at the fitted high-q retention
  // the degenerate resample chain recurses deeper, and 50000 frames overflowed
  // the stack guard page *before* this cap fired. Real scenarios never approach
  // a few thousand; only degenerate resample chains do.
};

template <class REAL> double MultiModel<REAL>::computeLogLikelihood() {
  if (this->_ccp.skip()) {
    return 0.0;
  }
  auto hash = getHash();
  // check if LL is already computed for this species tree
  auto cacheIt = this->_llCache.find(hash);
  if (cacheIt != this->_llCache.end()) {
    // too many collisions for dated getHash(),
    // delete the condition line upon fixing
    if (!this->_info.isDated())
      return cacheIt->second + this->_synLogConst;
  }
  if (this->_memorySavings) {
    // initialize CLVs, as they have not been created with the model
    allocateMemory();
  }
  // compute CLVs based on recomputed per-branch probas
  updateCLVs();
  // compute the species tree LL
  REAL res = REAL();
  std::vector<REAL> categoryLikelihoods(getGammaCatNumber(), REAL());
  for (unsigned int c = 0; c < getGammaCatNumber(); ++c) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      categoryLikelihoods[c] += getRootCladeLikelihood(speciesNode, c);
    }
    // condition on survival
    categoryLikelihoods[c] /=
        getLikelihoodFactor(c); // no need to scale here as factor <= 1.0
  }
  // sum over the categories
  res = std::accumulate(categoryLikelihoods.begin(), categoryLikelihoods.end(),
                        REAL());
  // normalize by the number of categories
  res /= static_cast<double>(getGammaCatNumber());
  scale(res);
  auto ll = getLog(res);
  // remember the computed LL for this species tree
  this->_llCache[hash] = ll;
  if (this->_memorySavings) {
    // delete CLVs to use less RAM
    deallocateMemory();
  }
  return ll + this->_synLogConst; // synteny emission constant (0 if none)
}

template <class REAL>
bool MultiModel<REAL>::inferMLScenario(Scenario &scenario) {
  assert(false); // not implemented for MultiModel
  return computeScenario(scenario, false);
}

template <class REAL>
bool MultiModel<REAL>::sampleReconciliations(
    unsigned int samples, std::vector<std::shared_ptr<Scenario>> &scenarios) {
  if (this->_memorySavings) {
    // initialize CLVs, as they have not been created with the model
    allocateMemory();
  }
  // compute CLVs based on recomputed per-branch probas
  updateCLVs();
  // sample and save scenarios
  bool ok = true;
  for (unsigned int i = 0; i < samples; ++i) {
    // A scenario can be discarded if it hits the backtrace depth guard (a
    // degenerate DL/TL resample chain); retry with a fresh random path and keep
    // only a *complete* (uncapped) scenario. A capped backtrace yields an
    // incomplete gene tree (a DL/TL-trapped lineage never reaches an observed
    // leaf), so it must not be emitted. A few pathological families (degenerate
    // likelihood under the DTL+WGD/LORe state) hit the guard on every attempt
    // and contribute fewer (or no) sampled gene trees; the run still terminates,
    // and the LORe resolution profile is produced by a separate backtrace.
    for (unsigned int attempt = 0; attempt < 8; ++attempt) {
      auto scenario = std::make_shared<Scenario>();
      if (computeScenario(*scenario, true)) {
        scenarios.push_back(scenario);
        break;
      }
    }
  }
  if (this->_memorySavings) {
    // delete CLVs to use less RAM
    deallocateMemory();
  }
  return ok;
}

template <class REAL> void MultiModel<REAL>::updateCLVs() {
  // recompute the per-branch probas of elementary events
  this->beforeComputeCLVs();
  // perform updateCLV() for each observed gene clade:
  // postorder gene clade traversal is granted
  for (CID cid = 0; cid < this->_ccp.getCladesNumber(); ++cid) {
    updateCLV(cid);
  }
  this->_allSpeciesNodesInvalid = false;
  this->_invalidatedSpeciesNodes.clear();
}

template <class REAL>
corax_rnode_t *
MultiModel<REAL>::sampleOriginationSpecies(unsigned int &category) {
  auto totalLikelihood = REAL();
  for (unsigned int c = 0; c < getGammaCatNumber(); ++c) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      totalLikelihood += getRootCladeLikelihood(speciesNode, c);
    }
  }
  auto toSample = getRandom(totalLikelihood);
  auto sumLikelihood = REAL();
  for (unsigned int c = 0; c < getGammaCatNumber(); ++c) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      sumLikelihood += getRootCladeLikelihood(speciesNode, c);
      if (sumLikelihood > toSample) {
        category = c;
        return speciesNode;
      }
    }
  }
  assert(false);
  return nullptr;
}

template <class REAL>
bool MultiModel<REAL>::computeScenario(Scenario &scenario, bool stochastic) {
  assert(stochastic);
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  // sample rate category and origination species
  unsigned int category = 0;
  auto originationSpecies = sampleOriginationSpecies(category);
  // init scenario
  scenario.setSpeciesTree(&this->_speciesTree);
  auto geneRoot = scenario.generateGeneRoot();
  scenario.setGeneRoot(geneRoot);
  auto virtualRootIndex = 2 * this->_ccp.getLeafNumber();
  scenario.setVirtualRootIndex(virtualRootIndex);
  // run sampling recursion
  _btScenarioSteps = 0; // reset the per-scenario depth guard
  auto ok = backtrace(rootCID, originationSpecies, geneRoot, category, scenario,
                      stochastic);
  return ok;
}

template <class REAL>
bool MultiModel<REAL>::backtrace(CID cid, corax_rnode_t *speciesNode,
                                 corax_unode_t *geneNode, unsigned int category,
                                 Scenario &scenario, bool stochastic) {
  // Bound the recursion DEPTH (not total calls): a degenerate DL/TL resample
  // chain recurses without bound and overflows the stack, while large but
  // shallow gene trees are unaffected. RAII keeps _btScenarioSteps equal to the
  // current depth (incremented on entry, decremented on every return path).
  struct DepthGuard {
    size_t &d;
    DepthGuard(size_t &x) : d(x) { ++d; }
    ~DepthGuard() { --d; }
  } depthGuard(_btScenarioSteps);
  if (_btScenarioSteps > kBtScenarioCap) {
    return false; // too deep -> discard this scenario and resample
  }
  auto c = category;
  REAL proba;
  // compute the probability of the clade on the species branch
  // to obtain the sampling probability from it further
  if (!computeProbability(cid, speciesNode, c, proba)) {
    return false;
  }
  // create a recCell, set its sampling probability and sample
  // an event of the clade on the species branch
  ReconciliationCell<REAL> recCell;
  recCell.maxProba = getRandom(proba);
  if (!computeProbability(cid, speciesNode, c, proba, &recCell)) {
    return false;
  }
  // fill the recCell's event fields:
  // - the current node index of the reconciled gene tree
  // - the current species node index
  // - the type of the previous event in the scenario
  if (scenario.getGeneNodeBuffer().size() == 1) {
    recCell.event.geneNode = scenario.getVirtualRootIndex();
  } else {
    recCell.event.geneNode = geneNode->node_index;
  }
  recCell.event.speciesNode = speciesNode->node_index;
  recCell.event.previousType = scenario.getLastEventType();
  // overwrite the recCell's event child node indices:
  // replace the cids of the sampled child clades with
  // the accordingly built child nodes of the reconciled gene tree
  CID leftCid = 0;
  CID rightCid = 0;
  corax_unode_t *leftGeneNode = nullptr;
  corax_unode_t *rightGeneNode = nullptr;
  if (recCell.event.type == ReconciliationEventType::EVENT_S ||
      recCell.event.type == ReconciliationEventType::EVENT_D ||
      recCell.event.type == ReconciliationEventType::EVENT_T) {
    // save the original clade indices
    leftCid = recCell.event.leftGeneIndex;
    rightCid = recCell.event.rightGeneIndex;
    // make child nodes for the scenario
    scenario.generateGeneChildren(geneNode, leftGeneNode, rightGeneNode);
    leftGeneNode->length = recCell.blLeft;
    rightGeneNode->length = recCell.blRight;
    recCell.event.leftGeneIndex = leftGeneNode->node_index;
    recCell.event.rightGeneIndex = rightGeneNode->node_index;
  }
  // add the overwritten recCell's event to the scenario
  bool addEvent = true;
  if (recCell.event.type == ReconciliationEventType::EVENT_DL ||
      (recCell.event.type == ReconciliationEventType::EVENT_TL &&
       recCell.event.pllDestSpeciesNode == nullptr)) {
    addEvent = false;
  }
  if (addEvent) {
    scenario.addEvent(recCell.event);
  }
  // recursion to sample for the child clades or to resample
  bool ok = true;
  std::string geneLabel;
  switch (recCell.event.type) {
  case ReconciliationEventType::EVENT_S:
    scenario.setLastEventType(recCell.event.type);
    ok &= descend(leftCid, this->getSpeciesLeft(speciesNode), leftGeneNode, c,
                  scenario, stochastic);
    scenario.setLastEventType(recCell.event.type);
    ok &= descend(rightCid, this->getSpeciesRight(speciesNode), rightGeneNode,
                  c, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_D:
    scenario.setLastEventType(recCell.event.type);
    ok &=
        backtrace(leftCid, speciesNode, leftGeneNode, c, scenario, stochastic);
    scenario.setLastEventType(recCell.event.type);
    ok &= backtrace(rightCid, speciesNode, rightGeneNode, c, scenario,
                    stochastic);
    break;
  case ReconciliationEventType::EVENT_T:
    // source species
    ok &=
        backtrace(leftCid, speciesNode, leftGeneNode, c, scenario, stochastic);
    // receiving species
    scenario.setLastEventType(ReconciliationEventType::EVENT_T);
    ok &= backtrace(rightCid, recCell.event.pllDestSpeciesNode, rightGeneNode,
                    c, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_SL:
    scenario.setLastEventType(recCell.event.type);
    ok &= descend(cid, recCell.event.pllDestSpeciesNode, geneNode, c,
                  scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_DL:
    // the gene got duplicated, but one copy was lost, we resample again
    ok &= backtrace(cid, speciesNode, geneNode, c, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_TL:
    if (recCell.event.pllDestSpeciesNode == nullptr) {
      // the gene was lost in the receiving species, we resample again
      ok &= backtrace(cid, speciesNode, geneNode, c, scenario, stochastic);
    } else {
      scenario.setLastEventType(ReconciliationEventType::EVENT_TL);
      ok &= backtrace(cid, recCell.event.pllDestSpeciesNode, geneNode, c,
                      scenario, stochastic);
    }
    break;
  case ReconciliationEventType::EVENT_None:
    // on leaf speciation
    geneLabel = recCell.event.label;
    geneNode->label = new char[geneLabel.size() + 1];
    memcpy(geneNode->label, geneLabel.c_str(), geneLabel.size() + 1);
    break;
  default:
    ok = false;
  }
  return ok;
}

template <class REAL>
bool MultiModel<REAL>::descend(CID cid, corax_rnode_t *speciesNode,
                               corax_unode_t *geneNode, unsigned int category,
                               Scenario &scenario, bool stochastic) {
  if (sampleWGDCoin(cid, speciesNode, category)) {
    return backtraceU(cid, speciesNode, geneNode, category, scenario,
                      stochastic);
  }
  return backtrace(cid, speciesNode, geneNode, category, scenario, stochastic);
}

template <class REAL>
bool MultiModel<REAL>::backtraceU(CID cid, corax_rnode_t *speciesNode,
                                  corax_unode_t *geneNode,
                                  unsigned int category, Scenario &scenario,
                                  bool stochastic) {
  struct DepthGuard {
    size_t &d;
    DepthGuard(size_t &x) : d(x) { ++d; }
    ~DepthGuard() { --d; }
  } depthGuard(_btScenarioSteps);
  if (_btScenarioSteps > kBtScenarioCap) {
    return false;
  }
  UEventSample ev;
  if (!sampleUEvent(cid, speciesNode, category, ev)) {
    return false;
  }
  auto makeEvent = [&](ReconciliationEventType type) {
    Scenario::Event event;
    event.type = type;
    event.geneNode = (scenario.getGeneNodeBuffer().size() == 1)
                         ? scenario.getVirtualRootIndex()
                         : geneNode->node_index;
    event.speciesNode = speciesNode->node_index;
    event.previousType = scenario.getLastEventType();
    return event;
  };
  bool ok = true;
  switch (ev.type) {
  case UEventSample::RESOLVE_SPLIT: {
    // the WGD ohnolog split: a duplication node flagged isWGD
    auto event = makeEvent(ReconciliationEventType::EVENT_D);
    event.isWGD = true;
    corax_unode_t *l = nullptr, *r = nullptr;
    scenario.generateGeneChildren(geneNode, l, r);
    l->length = ev.blLeft;
    r->length = ev.blRight;
    event.leftGeneIndex = l->node_index;
    event.rightGeneIndex = r->node_index;
    scenario.addEvent(event);
    scenario.setLastEventType(ReconciliationEventType::EVENT_D);
    ok &= backtrace(ev.left, speciesNode, l, category, scenario, stochastic);
    scenario.setLastEventType(ReconciliationEventType::EVENT_D);
    ok &= backtrace(ev.right, speciesNode, r, category, scenario, stochastic);
    break;
  }
  case UEventSample::RESOLVE_SINGLE:
    // one ohnolog lost: nothing observed, continue resolved on this branch
    ok &= backtrace(cid, speciesNode, geneNode, category, scenario, stochastic);
    break;
  case UEventSample::SPECIATE: {
    auto event = makeEvent(ReconciliationEventType::EVENT_S);
    corax_unode_t *l = nullptr, *r = nullptr;
    scenario.generateGeneChildren(geneNode, l, r);
    // l goes to the left species, r to the right species
    CID cl = ev.leftToLeftSpecies ? ev.left : ev.right;
    CID cr = ev.leftToLeftSpecies ? ev.right : ev.left;
    l->length = ev.leftToLeftSpecies ? ev.blLeft : ev.blRight;
    r->length = ev.leftToLeftSpecies ? ev.blRight : ev.blLeft;
    event.leftGeneIndex = l->node_index;
    event.rightGeneIndex = r->node_index;
    scenario.addEvent(event);
    scenario.setLastEventType(ReconciliationEventType::EVENT_S);
    ok &= backtraceU(cl, this->getSpeciesLeft(speciesNode), l, category,
                     scenario, stochastic);
    scenario.setLastEventType(ReconciliationEventType::EVENT_S);
    ok &= backtraceU(cr, this->getSpeciesRight(speciesNode), r, category,
                     scenario, stochastic);
    break;
  }
  case UEventSample::SPECIATE_LOSS: {
    auto event = makeEvent(ReconciliationEventType::EVENT_SL);
    event.lostSpeciesNode = ev.lost->node_index;
    event.pllDestSpeciesNode = ev.dest;
    event.pllLostSpeciesNode = ev.lost;
    scenario.addEvent(event);
    scenario.setLastEventType(ReconciliationEventType::EVENT_SL);
    ok &= backtraceU(cid, ev.dest, geneNode, category, scenario, stochastic);
    break;
  }
  case UEventSample::LEAF:
    // (1-r) R at a tip: the unresolved locus is observed as a single gene
    ok &= backtrace(cid, speciesNode, geneNode, category, scenario, stochastic);
    break;
  }
  return ok;
}
