#pragma once

#include <cstdlib> // std::abort for the always-on LORe consistency guard

#include "MultiModel.hpp"

/**
 *  Implements all HGT modelling-dependent functions of
 *  the MultiModel class in the context of the UndatedDL
 *  model
 */
template <class REAL> class UndatedDLMultiModel : public MultiModel<REAL> {
public:
  UndatedDLMultiModel(PLLRootedTree &speciesTree,
                      const GeneSpeciesMapping &geneSpeciesMapping,
                      const RecModelInfo &info, const std::string &ccpFile);

  virtual ~UndatedDLMultiModel() {}

  virtual void setAlpha(double alpha);
  virtual void setRates(const RatesVector &rates);

  // --- WGD extension ---
  // Declare a WGD at the top of species branch e with retention prob q.
  void setWGD(unsigned int e, double q) override {
    _hasWGD[e] = 1;
    _q[e] = q;
    this->invalidateAllSpeciesNodes();
    this->resetCache();
  }

  // --- LORe extension ---
  // Set the resolution probability r in [0,1] on EVERY species branch (the
  // phase-1 global-r behaviour). r = 1 recovers the WHALE/AORe model.
  void setResolutionProb(double r) override {
    std::fill(_resolutionProbs.begin(), _resolutionProbs.end(), r);
    this->invalidateAllSpeciesNodes();
    this->resetCache();
  }

  // Per-event r: assign a full per-branch resolution vector at once (index ==
  // species node_index). Defensive on size: branches not covered stay at r=1.
  // One cache reset per call.
  void setResolutionProbVector(const std::vector<double> &perBranchR) override {
    std::fill(_resolutionProbs.begin(), _resolutionProbs.end(), 1.0);
    const size_t n = std::min(perBranchR.size(), _resolutionProbs.size());
    for (size_t e = 0; e < n; ++e) {
      _resolutionProbs[e] = perBranchR[e];
    }
    this->invalidateAllSpeciesNodes();
    this->resetCache();
  }

  // Set the resolution probability on a single species branch e (phase 1.5:
  // branch-specific r, used to place / recover lineage-specific resolution).
  void setResolutionProbBranch(unsigned int e, double r) {
    _resolutionProbs[e] = r;
    this->invalidateAllSpeciesNodes();
    this->resetCache();
  }
  double getResolutionProbBranch(unsigned int e) const {
    return _resolutionProbs[e];
  }

  // --- LORe resolution-branch marginal (WGD_LORE_marginal.md, phase 1b) ---
  // Sample `samples` resolution histories under the CURRENT (fitted) global r
  // via a U-aware backtrace, and accumulate per species branch the expected
  // number of U->R commit (resolution / ohnolog-divergence) events. CLVs must
  // be current: call computeLogLikelihood() first. With `check`, asserts the
  // STEP 5.2 consistency (sampled-term weights sum to the inside CLV) at every
  // U cell and at the WGD R-vs-U coin.
  void sampleResolutionCommits(unsigned int samples,
                               std::vector<double> &commitCounts,
                               std::vector<double> &tetraCounts, bool check);

  // Interface override (used by the reconciliation pipeline): no consistency
  // check, current CLVs assumed.
  void sampleResolutionCommits(unsigned int samples,
                               std::vector<double> &commitCounts,
                               std::vector<double> &tetraCounts) override {
    sampleResolutionCommits(samples, commitCounts, tetraCounts, false);
  }

protected:
  using UEventSample = typename MultiModel<REAL>::UEventSample;
  // scenario-sampler hooks (MultiModel::descend / backtraceU)
  bool sampleWGDCoin(CID cid, corax_rnode_t *sp, unsigned int c) override;
  bool sampleUEvent(CID cid, corax_rnode_t *sp, unsigned int c,
                    UEventSample &ev) override {
    return sampleUEventImpl(cid, sp, c, ev, false);
  }

private:
  // one draw from the U recursion (shared by the scenario sampler and btU)
  bool sampleUEventImpl(CID cid, corax_rnode_t *sp, unsigned int c,
                        UEventSample &ev, bool check);
  // U-aware backtrace helpers for the resolution-branch marginal. They reuse
  // computeProbability() to sample the resolved (R-state) events, and add the
  // unresolved (U-state) recursion and the WGD R/U coin on top.
  corax_rnode_t *sampleOriginationR(unsigned int &category);
  void btR(CID cid, corax_rnode_t *sp, unsigned int c,
           std::vector<unsigned int> &commits, bool check);
  void btU(CID cid, corax_rnode_t *sp, unsigned int c,
           std::vector<unsigned int> &commits, bool check);
  void btDescend(CID cid, corax_rnode_t *sp, unsigned int c,
                 std::vector<unsigned int> &commits, bool check);

  unsigned int _gammaCatNumber;
  std::vector<double> _gammaScalers;
  RatesVector _dlRates;
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _OP; // Origination probability, per species branch
  std::vector<REAL> _uE;   // Extinction probability, per species branch
  OriginationStrategy _originationStrategy;
  // Leaf species at which a U (unresolved/tetrasomic) lineage reaches an extant
  // tip without committing -> still tetrasomic today. Filled by btU, read+cleared
  // per sample in sampleResolutionCommits (parallel to the resolved `commits`).
  std::vector<unsigned int> _tetra;

  // Element e of the gene clade's DLCLV stores the probability of the clade,
  // given the clade is mapped to the species branch e.
  // In the paper: Pi_{e,gamma} of a clade gamma for each branch e
  using DLCLV = std::vector<REAL>;
  // vector of DLCLVs for all observed clades
  std::vector<DLCLV> _dlclvs;

  // --- WGD extension ---
  std::vector<char> _hasWGD;     // per species branch: WGD at top of branch?
  std::vector<double> _q;        // per species branch: retention prob (if WGD)
  std::vector<REAL> _uEtop;      // extinction seen by PARENT, per (e,cat)
  std::vector<DLCLV> _dlclvsTop; // clade CLV seen by PARENT, per clade per (e,cat)

  // --- LORe (delayed rediploidization) extension ---
  std::vector<DLCLV> _uclvs;       // unresolved (tetrasomic) CLV, per clade per (e,cat)
  std::vector<REAL> _uEU;          // unresolved extinction, per (e,cat)
  std::vector<double> _resolutionProbs; // per species branch: r in [0,1]; 1==AORe

  // Doubling bracket shared by the WGD and LORe transforms:
  //   2*E*R + sum_splits ccp * R_left * R_right   (no _PD, no q/r weighting).
  // Reads R (=_dlclvs) and E (=_uE) at node-cat ec; _dlclvs[cid][ec] must
  // already hold R for this clade (set just before this is called in updateCLV).
  REAL dupBracket(unsigned int cid, unsigned int ec) {
    REAL split = REAL();
    const auto &splits = this->_ccp.getCladeSplits(cid);
    for (size_t k = 0; k < splits.size(); ++k) {
      const auto &cs = splits[k];
      // WGD split: synteny-weighted frequency (== ccp frequency without synteny)
      REAL t = _dlclvs[cs.left][ec] * _dlclvs[cs.right][ec] *
               this->synWgdFreq(cid, k, cs.frequency);
      scale(t);
      split += t;
    }
    REAL out = _uE[ec] * _dlclvs[cid][ec] * 2.0 + split;
    scale(out);
    return out;
  }

  // Normaliser of the unresolved (U) event probabilities. A tetrasomic U
  // lineage has no duplication channel, so its per-branch events (speciation,
  // loss) are renormalised over PS + PL to sum to 1; without this the U
  // process is defective (mass (1-r)*PD lost per branch) and r is biased to 1.
  double uDen(unsigned int ec) const { return _PS[ec] + _PL[ec]; }
  // functions to work with CLVs
  virtual void allocateMemory();
  virtual void deallocateMemory();
  virtual void updateCLV(CID cid);

  // functions to work with probabilities
  virtual void recomputeSpeciesProbabilities();
  virtual double getLikelihoodFactor(unsigned int category);
  virtual REAL getRootCladeLikelihood(corax_rnode_t *speciesNode,
                                      unsigned int category);
  virtual bool computeProbability(CID cid, corax_rnode_t *speciesNode,
                                  unsigned int category, REAL &proba,
                                  ReconciliationCell<REAL> *recCell = nullptr);

  // functions to work with _llCache
  virtual size_t getHash() { return this->getSpeciesTreeHash(); }
};

/**
 *  Constructor
 */
template <class REAL>
UndatedDLMultiModel<REAL>::UndatedDLMultiModel(
    PLLRootedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMapping,
    const RecModelInfo &info, const std::string &ccpFile)
    : MultiModel<REAL>(speciesTree, geneSpeciesMapping, info, ccpFile),
      _gammaCatNumber(info.gammaCategories),
      _gammaScalers(_gammaCatNumber, 1.0),
      _PD(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 0.2),
      _PL(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 0.2),
      _PS(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 1.0),
      _OP(this->getAllSpeciesNodeNumber(),
          1.0 / static_cast<double>(this->getAllSpeciesNodeNumber())),
      _uE(this->getAllSpeciesNodeNumber() * _gammaCatNumber, REAL()),
      _originationStrategy(info.originationStrategy) {
  auto N = this->getAllSpeciesNodeNumber();
  // --- WGD extension: per-branch state, no WGD by default ---
  auto Nc = N * _gammaCatNumber;
  _hasWGD.assign(N, 0);
  _q.assign(N, 1.0);
  _uEtop.assign(Nc, REAL());
  // --- LORe extension: unresolved extinction (CLV is in allocateMemory) ---
  _uEU.assign(Nc, REAL());
  _resolutionProbs.assign(N, 1.0); // default r=1 everywhere == AORe/WHALE
  // _dlclvsTop and _uclvs are allocated alongside _dlclvs in allocateMemory()
  // set gamma scalers with the default alpha
  setAlpha(1.0);
  // set all DLO rates to the default value
  _dlRates.resize(this->_info.modelFreeParameters(),
                  std::vector<double>(N, 0.2));
  // initialize DLCLVs if needed
  if (!this->_memorySavings) {
    allocateMemory();
  }
}

template <class REAL> void UndatedDLMultiModel<REAL>::setAlpha(double alpha) {
  corax_compute_gamma_cats(alpha, _gammaScalers.size(), &_gammaScalers[0],
                           CORAX_GAMMA_RATES_MEAN);
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

template <class REAL>
void UndatedDLMultiModel<REAL>::setRates(const RatesVector &rates) {
  assert(rates.size() == this->_info.modelFreeParameters());
  _dlRates = rates;
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

/**
 *  Allocate memory to the CLVs
 */
template <class REAL> void UndatedDLMultiModel<REAL>::allocateMemory() {
  DLCLV nullCLV(this->getAllSpeciesNodeNumber() * _gammaCatNumber, REAL());
  _dlclvs = std::vector<DLCLV>(this->_ccp.getCladesNumber(), nullCLV);
  _dlclvsTop = std::vector<DLCLV>(this->_ccp.getCladesNumber(), nullCLV);
  _uclvs = std::vector<DLCLV>(this->_ccp.getCladesNumber(), nullCLV);
}

/**
 *  Free memory allocated to the CLVs
 */
template <class REAL> void UndatedDLMultiModel<REAL>::deallocateMemory() {
  _dlclvs = std::vector<DLCLV>();
  _dlclvsTop = std::vector<DLCLV>();
  _uclvs = std::vector<DLCLV>();
}

/**
 *  Compute the CLV for a given clade
 */
template <class REAL> void UndatedDLMultiModel<REAL>::updateCLV(CID cid) {
  auto &uq = _dlclvs[cid];
  auto &uu = _uclvs[cid];
  std::fill(uq.begin(), uq.end(), REAL());
  std::fill(uu.begin(), uu.end(), REAL());
  bool ok;
  for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
    // postorder species tree traversal is granted
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      auto e = speciesNode->node_index;
      auto ec = e * _gammaCatNumber + c;
      const double r = _resolutionProbs[e]; // per-branch resolution prob
      REAL p = REAL();
      ok = computeProbability(cid, speciesNode, c, p);
      assert(ok);
      uq[ec] = p; // resolved CLV R

      // --- LORe: unresolved (tetrasomic) CLV U (STEP 4) ---
      //   U = r * dupBracket(cid,ec)                 // resolve here -> doubling
      //     + (1-r) * [ speciation, both daughters inherit U (both orderings)
      //                 + speciation with unresolved loss of one daughter ]
      //   leaf node: U = r*dupBracket + (1-r)*R  (never resolved -> single gene)
      // U carries NO duplication channel in phase 1 (a still-tetrasomic locus
      // is not independently duplicated). r=1 => U == dupBracket => WHALE.
      REAL uval = dupBracket(cid, ec);
      uval *= r;
      scale(uval);
      if (r < 1.0) {
        REAL rest = REAL();
        if (this->getSpeciesLeft(speciesNode)) {
          auto f = this->getSpeciesLeft(speciesNode)->node_index;
          auto g = this->getSpeciesRight(speciesNode)->node_index;
          auto fc = f * _gammaCatNumber + c;
          auto gc = g * _gammaCatNumber + c;
          // speciation: both daughters inherit U (both child orderings)
          for (const auto &cs : this->_ccp.getCladeSplits(cid)) {
            REAL t = (_uclvs[cs.left][fc] * _uclvs[cs.right][gc] +
                      _uclvs[cs.right][fc] * _uclvs[cs.left][gc]) *
                     (_PS[ec] / uDen(ec) * cs.frequency);
            scale(t);
            rest += t;
          }
          // speciation + unresolved loss of one daughter (SL acting on U)
          REAL sl = (uu[fc] * _uEU[gc] + uu[gc] * _uEU[fc]) * (_PS[ec] / uDen(ec));
          scale(sl);
          rest += sl;
        } else {
          // leaf species node: never resolved, seen as the single leaf gene
          rest = p;
        }
        rest *= (1.0 - r);
        scale(rest);
        uval += rest;
        scale(uval);
      }
      uu[ec] = uval;

      // --- WGD transform: the retained fraction enters the U state (STEP 5) ---
      auto &uqTop = _dlclvsTop[cid];
      if (_hasWGD[e]) {
        double q = _q[e];
        REAL keep = p * (1.0 - q); // (1-q) R
        REAL inj = uu[ec] * q;     // q U  (== q*dupBracket when r=1 -> WHALE)
        scale(keep);
        scale(inj);
        uqTop[ec] = keep + inj;
      } else {
        uqTop[ec] = p;
      }
    }
  }
}

/**
 *  Compute the per species branch probabilities of
 *  the elementary events of clade evolution
 */
template <class REAL>
void UndatedDLMultiModel<REAL>::recomputeSpeciesProbabilities() {
  auto allSpeciesNumber = this->getAllSpeciesNodeNumber();
  // recompute _PD, _PL, _PS
  auto &dupRates = _dlRates[0];
  auto &lossRates = _dlRates[1];
  assert(allSpeciesNumber == dupRates.size());
  assert(allSpeciesNumber == lossRates.size());
  std::fill(_PD.begin(), _PD.end(), 0.0);
  std::fill(_PL.begin(), _PL.end(), 0.0);
  std::fill(_PS.begin(), _PS.end(), 0.0);
  for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      auto e = speciesNode->node_index;
      auto ec = e * _gammaCatNumber + c;
      _PD[ec] = dupRates[e];
      _PL[ec] = lossRates[e];
      _PS[ec] = _gammaScalers[c];
      if (this->_info.noDup) {
        assert(false);
        _PD[ec] = 0.0;
      }
      auto sum = _PD[ec] + _PL[ec] + _PS[ec];
      _PD[ec] /= sum;
      _PL[ec] /= sum;
      _PS[ec] /= sum;
    }
  }
  // recompute _OP
  std::vector<corax_rnode_t *> speciesNodesBuffer;
  std::vector<corax_rnode_t *> *possibleOriginationSpeciesNodes = nullptr;
  switch (_originationStrategy) {
  case OriginationStrategy::UNIFORM:
  case OriginationStrategy::OPTIMIZE:
    possibleOriginationSpeciesNodes = &(this->getPrunedSpeciesNodes());
    break;
  case OriginationStrategy::ROOT:
    speciesNodesBuffer.push_back(this->getPrunedRoot());
    possibleOriginationSpeciesNodes = &speciesNodesBuffer;
    break;
  case OriginationStrategy::LCA:
    speciesNodesBuffer.push_back(this->getCoveredSpeciesLCA());
    possibleOriginationSpeciesNodes = &speciesNodesBuffer;
    break;
  }
  std::fill(_OP.begin(), _OP.end(), 0.0);
  if (_originationStrategy == OriginationStrategy::OPTIMIZE) {
    auto &oriRates = _dlRates[2];
    assert(allSpeciesNumber == oriRates.size());
    double sum = 0.0;
    for (auto speciesNode : *possibleOriginationSpeciesNodes) {
      auto e = speciesNode->node_index;
      sum += oriRates[e];
    }
    for (auto speciesNode : *possibleOriginationSpeciesNodes) {
      auto e = speciesNode->node_index;
      _OP[e] = oriRates[e] / sum;
    }
  } else {
    double sum = static_cast<double>(possibleOriginationSpeciesNodes->size());
    for (auto speciesNode : *possibleOriginationSpeciesNodes) {
      auto e = speciesNode->node_index;
      _OP[e] = 1.0 / sum;
    }
  }
  // recompute _uE (and _uEtop, the post-WGD extinction seen by the parent)
  std::fill(_uE.begin(), _uE.end(), REAL());
  std::fill(_uEtop.begin(), _uEtop.end(), REAL());
  // iterate several times to resolve _uE probas with
  // fixed point optimization
  unsigned int maxIt = 4;
  for (unsigned int it = 0; it < maxIt; ++it) {
    for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
      // postorder species tree traversal is granted
      for (auto speciesNode : this->getPrunedSpeciesNodes()) {
        auto e = speciesNode->node_index;
        auto ec = e * _gammaCatNumber + c;
        REAL temp;
        REAL proba = REAL();
        // L scenario
        temp = REAL(_PL[ec]);
        scale(temp);
        proba += temp;
        // S scenario
        if (this->getSpeciesLeft(speciesNode)) {
          // internal branch
          auto f = this->getSpeciesLeft(speciesNode)->node_index;
          auto g = this->getSpeciesRight(speciesNode)->node_index;
          auto fc = f * _gammaCatNumber + c;
          auto gc = g * _gammaCatNumber + c;
          // read children's POST-WGD extinction (_uEtop)
          temp = _uEtop[fc] * _uEtop[gc] * _PS[ec]; // SEE scenario
        } else {
          // terminal branch
          temp = REAL(_PS[ec] * this->_fm[e]); // S but not observed scenario
        }
        scale(temp);
        proba += temp;
        // DEE scenario (own branch, below WGD -> raw _uE)
        temp = _uE[ec] * _uE[ec] * _PD[ec];
        scale(temp);
        proba += temp;
        assert(proba < REAL(1.000001));
        _uE[ec] = proba;
        // --- LORe: unresolved (tetrasomic) extinction EU (STEP 3) ---
        //   EU = r*E^2 + (1-r)*( _PL + _PS*EU[fc]*EU[gc] )    (internal)
        //   EU = r*E^2 + (1-r)*E                              (leaf)
        // EU has no self-dependence, so it is a single postorder pass riding on
        // the E fixed point. r=1 => EU == E^2 (the WHALE doubling term).
        {
          const double r = _resolutionProbs[e];
          REAL E = _uE[ec];
          REAL eu = E * E;
          eu *= r;
          scale(eu);
          if (r < 1.0) {
            REAL rest = REAL();
            if (this->getSpeciesLeft(speciesNode)) {
              auto f = this->getSpeciesLeft(speciesNode)->node_index;
              auto g = this->getSpeciesRight(speciesNode)->node_index;
              auto fc = f * _gammaCatNumber + c;
              auto gc = g * _gammaCatNumber + c;
              REAL lossU = REAL(_PL[ec] / uDen(ec)); // lost while still unresolved
              scale(lossU);
              rest += lossU;
              REAL specU = _uEU[fc] * _uEU[gc] * (_PS[ec] / uDen(ec)); // speciate, both U
              scale(specU);
              rest += specU;
            } else {
              rest = E; // leaf: single-lineage non-observation prob (via _fm)
            }
            rest *= (1.0 - r);
            scale(rest);
            eu += rest;
            scale(eu);
          }
          _uEU[ec] = eu;
        }
        // --- WGD transform: the retained fraction enters the U state (STEP 5) ---
        // E_top = (1-q) E + q EU  (== (1-q)E + q E^2 when r=1 -> WHALE/AORe).
        if (_hasWGD[e]) {
          double q = _q[e];
          REAL E = _uE[ec];
          REAL keep = E * (1.0 - q); // (1-q) E
          REAL inj = _uEU[ec] * q;   // q EU
          scale(keep);
          scale(inj);
          _uEtop[ec] = keep + inj;
        } else {
          _uEtop[ec] = _uE[ec];
        }
      }
    }
  } // end of iteration
}

/**
 *  Correction factor to the species tree likelihood,
 *  because we condition on survival
 */
template <class REAL>
double UndatedDLMultiModel<REAL>::getLikelihoodFactor(unsigned int category) {
  double factor(0.0);
  auto c = category;
  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    auto ec = e * _gammaCatNumber + c;
    factor += (1.0 - _uE[ec]) * _OP[e];
  }
  return factor;
}

/**
 *  Probability of the current family to evolve starting
 *  from a given species branch
 */
template <class REAL>
REAL UndatedDLMultiModel<REAL>::getRootCladeLikelihood(
    corax_rnode_t *speciesNode, unsigned int category) {
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto c = category;
  auto e = speciesNode->node_index;
  auto ec = e * _gammaCatNumber + c;
  REAL likelihood = _dlclvs[rootCID][ec] * _OP[e];
  scale(likelihood);
  return likelihood;
}

/**
 *  Compute the CLV value for a given cid (clade id) and a given
 *  species node and write it to the proba variable
 */
template <class REAL>
bool UndatedDLMultiModel<REAL>::computeProbability(
    CID cid, corax_rnode_t *speciesNode, unsigned int category, REAL &proba,
    ReconciliationCell<REAL> *recCell) {
  proba = REAL();
  REAL temp;
  bool isGeneLeaf = this->_ccp.isLeaf(cid);
  bool isSpeciesLeaf = !this->getSpeciesLeft(speciesNode);
  auto c = category;
  auto e = speciesNode->node_index;
  auto ec = e * _gammaCatNumber + c;
  unsigned int f = 0;
  unsigned int g = 0;
  unsigned int fc = 0;
  unsigned int gc = 0;
  if (!isSpeciesLeaf) {
    f = this->getSpeciesLeft(speciesNode)->node_index;
    g = this->getSpeciesRight(speciesNode)->node_index;
    fc = f * _gammaCatNumber + c;
    gc = g * _gammaCatNumber + c;
  }
  REAL maxProba = REAL();
  if (recCell) {
    maxProba = recCell->maxProba;
  }
  // S events on terminal species branches can happen
  // for terminal gene nodes only:
  if (isGeneLeaf) {
    // - S event on a terminal species branch (only for compatible genes and
    // species)
    if (isSpeciesLeaf && this->_geneToSpecies[cid] == e) {
      temp = REAL(_PS[ec]);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_None;
        recCell->event.label = this->_ccp.getLeafLabel(cid);
        return true;
      }
    }
  }
  // S events on internal species branches and D events can happen
  // for ancestral gene nodes only:
  for (const auto &cladeSplit : this->_ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left;
    auto cidRight = cladeSplit.right;
    auto freq = cladeSplit.frequency;
    // - S event on an internal species branch
    if (!isSpeciesLeaf) {
      temp = _dlclvsTop[cidLeft][fc] * _dlclvsTop[cidRight][gc] * (_PS[ec] * freq);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidLeft;
        recCell->event.rightGeneIndex = cidRight;
        recCell->blLeft = cladeSplit.blLeft;
        recCell->blRight = cladeSplit.blRight;
        return true;
      }
      temp = _dlclvsTop[cidRight][fc] * _dlclvsTop[cidLeft][gc] * (_PS[ec] * freq);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidRight;
        recCell->event.rightGeneIndex = cidLeft;
        recCell->blLeft = cladeSplit.blRight;
        recCell->blRight = cladeSplit.blLeft;
        return true;
      }
    }
    // - D event
    temp = _dlclvs[cidLeft][ec] * _dlclvs[cidRight][ec] * (_PD[ec] * freq);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_D;
      recCell->event.leftGeneIndex = cidLeft;
      recCell->event.rightGeneIndex = cidRight;
      recCell->blLeft = cladeSplit.blLeft;
      recCell->blRight = cladeSplit.blRight;
      return true;
    }
  }
  // SL events and DL events can happen
  // for any of gene nodes:
  // - SL event (only on an internal species branch)
  if (!isSpeciesLeaf) {
    temp = _dlclvsTop[cid][fc] * (_uEtop[gc] * _PS[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.lostSpeciesNode = g;
      recCell->event.pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      recCell->event.pllLostSpeciesNode = this->getSpeciesRight(speciesNode);
      return true;
    }
    temp = _dlclvsTop[cid][gc] * (_uEtop[fc] * _PS[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.lostSpeciesNode = f;
      recCell->event.pllDestSpeciesNode = this->getSpeciesRight(speciesNode);
      recCell->event.pllLostSpeciesNode = this->getSpeciesLeft(speciesNode);
      return true;
    }
  }
  if (!this->_info.noDL) {
    // - DL event
    temp = proba / (1.0 - (_uE[ec] * _PD[ec] * 2.0));
    scale(temp);
    proba = temp;
    if (recCell && proba > maxProba) {
      // in fact, nothing happens, we'll have to resample
      recCell->event.type = ReconciliationEventType::EVENT_DL;
      return true;
    }
  }
  if (recCell) {
    Logger::error << "error: proba=" << proba << ", maxProba=" << maxProba
                  << " (proba < maxProba)" << std::endl;
    return false; // we haven't sampled any event, this should not happen
  }
  if (proba > REAL(1.0)) {
    Logger::error << "error: proba=" << proba << " (proba > 1.0)" << std::endl;
    return false;
  }
  return true;
}

/**
 *  LORe resolution-branch marginal (WGD_LORE_marginal.md)
 *
 *  A U-aware backtrace. computeProbability() already samples the resolved
 *  (R-state) events S/D/SL/DL/None proportionally to their inside contribution.
 *  On top of that we add:
 *   - the WGD R-vs-U coin when a lineage is inherited across the WGD branch
 *     (STEP 2), and
 *   - the unresolved (U-state) recursion (STEP 3), which emits a U->R commit
 *     (records the species branch) when resolution fires.
 *  All sampling weights are the exact terms of the inside recursion (same _PS,
 *  ccp frequencies, q/r prefactors), so the sampled marginals are unbiased; the
 *  `check` flag asserts weight==inside at every U cell and the coin (STEP 5.2).
 */

// helper: |log a - log b| small (robust to ScaledValue scaling); skips null
template <class REAL> static inline bool loreInsideClose(const REAL &a,
                                                         const REAL &b) {
  double la = getLog(a);
  double lb = getLog(b);
  if (!std::isfinite(la) && !std::isfinite(lb)) {
    return true; // both effectively zero
  }
  return std::fabs(la - lb) < 1e-6;
}

template <class REAL>
corax_rnode_t *
UndatedDLMultiModel<REAL>::sampleOriginationR(unsigned int &category) {
  REAL total = REAL();
  for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
    for (auto sp : this->getPrunedSpeciesNodes()) {
      total += getRootCladeLikelihood(sp, c);
    }
  }
  REAL toSample = this->getRandom(total);
  REAL acc = REAL();
  corax_rnode_t *last = nullptr;
  for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
    for (auto sp : this->getPrunedSpeciesNodes()) {
      acc += getRootCladeLikelihood(sp, c);
      last = sp;
      if (acc > toSample) {
        category = c;
        return sp;
      }
    }
  }
  category = 0;
  return last;
}

template <class REAL>
void UndatedDLMultiModel<REAL>::sampleResolutionCommits(
    unsigned int samples, std::vector<double> &commitCounts,
    std::vector<double> &tetraCounts, bool check) {
  commitCounts.assign(this->getAllSpeciesNodeNumber(), 0.0);
  tetraCounts.assign(this->getAllSpeciesNodeNumber(), 0.0);
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  for (unsigned int s = 0; s < samples; ++s) {
    unsigned int cat = 0;
    auto origin = sampleOriginationR(cat);
    // A lineage that ORIGINATES anywhere starts resolved (the inside root term
    // uses raw _dlclvs, not the WGD-transformed _dlclvsTop).
    std::vector<unsigned int> commits;
    _tetra.clear();
    btR(rootCID, origin, cat, commits, check);
    for (auto b : commits) {
      commitCounts[b] += 1.0;
    }
    for (auto b : _tetra) {
      tetraCounts[b] += 1.0; // leaf lineages still tetrasomic at the tip
    }
  }
}

// WGD R-vs-U coin at the top of branch sp, with the exact weights of the
// parent's consumed value R_top = (1-q) R + q U. Returns true => enter U.
template <class REAL>
bool UndatedDLMultiModel<REAL>::sampleWGDCoin(CID cid, corax_rnode_t *sp,
                                              unsigned int c) {
  auto e = sp->node_index;
  if (!_hasWGD[e]) {
    return false;
  }
  auto ec = e * _gammaCatNumber + c;
  double q = _q[e];
  REAL wR = _dlclvs[cid][ec] * (1.0 - q);
  REAL wU = _uclvs[cid][ec] * q;
  scale(wR);
  scale(wU);
  REAL total = wR + wU;
  REAL toSample = this->getRandom(total);
  return !(wR > toSample);
}

// Descend into a child species node; sample the WGD R-vs-U coin if it is the
// WGD branch (STEP 2), otherwise continue resolved.
template <class REAL>
void UndatedDLMultiModel<REAL>::btDescend(CID cid, corax_rnode_t *sp,
                                          unsigned int c,
                                          std::vector<unsigned int> &commits,
                                          bool check) {
  auto e = sp->node_index;
  if (check && _hasWGD[e]) {
    auto ec = e * _gammaCatNumber + c;
    REAL wR = _dlclvs[cid][ec] * (1.0 - _q[e]);
    REAL wU = _uclvs[cid][ec] * _q[e];
    scale(wR);
    scale(wU);
    if (!loreInsideClose(REAL(wR + wU), _dlclvsTop[cid][ec])) {
      std::cerr << "LORe STEP5.2 FAIL: WGD coin weights != _dlclvsTop at cid="
                << cid << " e=" << e << std::endl;
      std::abort();
    }
  }
  if (sampleWGDCoin(cid, sp, c)) {
    btU(cid, sp, c, commits, check);
  } else {
    btR(cid, sp, c, commits, check);
  }
}

// Resolved-state backtrace: reuse computeProbability to sample the R event,
// then recurse (applying the WGD coin when descending into the WGD branch).
template <class REAL>
void UndatedDLMultiModel<REAL>::btR(CID cid, corax_rnode_t *sp, unsigned int c,
                                    std::vector<unsigned int> &commits,
                                    bool check) {
  REAL proba = REAL();
  if (!computeProbability(cid, sp, c, proba)) {
    return;
  }
  ReconciliationCell<REAL> recCell;
  recCell.maxProba = this->getRandom(proba);
  REAL tmp = REAL();
  if (!computeProbability(cid, sp, c, tmp, &recCell)) {
    return;
  }
  switch (recCell.event.type) {
  case ReconciliationEventType::EVENT_None:
    return; // observed as a leaf gene
  case ReconciliationEventType::EVENT_S:
    btDescend(recCell.event.leftGeneIndex, this->getSpeciesLeft(sp), c, commits,
              check);
    btDescend(recCell.event.rightGeneIndex, this->getSpeciesRight(sp), c,
              commits, check);
    return;
  case ReconciliationEventType::EVENT_D:
    // duplication within the branch (below the WGD): children stay resolved
    btR(recCell.event.leftGeneIndex, sp, c, commits, check);
    btR(recCell.event.rightGeneIndex, sp, c, commits, check);
    return;
  case ReconciliationEventType::EVENT_SL:
    btDescend(cid, recCell.event.pllDestSpeciesNode, c, commits, check);
    return;
  case ReconciliationEventType::EVENT_DL:
    btR(cid, sp, c, commits, check); // duplication+loss -> resample this cell
    return;
  default:
    return;
  }
}

// One draw from the unresolved (U) recursion at (cid, sp). Every weight is the
// exact inside term of updateCLV (normalised U events, synteny rho on the
// WGD split), so the sampled scenario/resolution marginals are exact.
template <class REAL>
bool UndatedDLMultiModel<REAL>::sampleUEventImpl(CID cid, corax_rnode_t *sp,
                                                 unsigned int c,
                                                 UEventSample &ev,
                                                 bool check) {
  auto e = sp->node_index;
  auto ec = e * _gammaCatNumber + c;
  double r = _resolutionProbs[e];
  bool isLeaf = !this->getSpeciesLeft(sp);
  unsigned int fc = 0, gc = 0;
  if (!isLeaf) {
    fc = this->getSpeciesLeft(sp)->node_index * _gammaCatNumber + c;
    gc = this->getSpeciesRight(sp)->node_index * _gammaCatNumber + c;
  }
  const auto &splits = this->_ccp.getCladeSplits(cid);
  const double psU = _PS[ec] / uDen(ec);
  REAL wResolve = dupBracket(cid, ec);
  wResolve *= r;
  scale(wResolve);
  REAL total = wResolve;
  if (isLeaf) {
    REAL wLeaf = _dlclvs[cid][ec] * (1.0 - r);
    scale(wLeaf);
    total = total + wLeaf;
  } else {
    for (const auto &cs : splits) {
      double w0 = (1.0 - r) * psU * cs.frequency;
      REAL s1 = _uclvs[cs.left][fc] * _uclvs[cs.right][gc] * w0;
      REAL s2 = _uclvs[cs.right][fc] * _uclvs[cs.left][gc] * w0;
      scale(s1);
      scale(s2);
      total = total + s1 + s2;
    }
    double w1 = (1.0 - r) * psU;
    REAL sl1 = _uclvs[cid][fc] * _uEU[gc] * w1;
    REAL sl2 = _uclvs[cid][gc] * _uEU[fc] * w1;
    scale(sl1);
    scale(sl2);
    total = total + sl1 + sl2;
  }
  if (check && !loreInsideClose(total, _uclvs[cid][ec])) {
    std::cerr << "LORe STEP5.2 FAIL: U-cell weights != _uclvs at cid=" << cid
              << " e=" << e << std::endl;
    std::abort();
  }
  if (!(total > REAL())) {
    return false;
  }
  REAL toSample = this->getRandom(total);
  REAL acc = wResolve;
  if (wResolve > toSample) {
    // resolve on branch e: sub-sample within 2*E*R + sum_k rho_k ccp_k R_l R_r
    REAL dTotal = dupBracket(cid, ec);
    REAL dSample = this->getRandom(dTotal);
    REAL dAcc = REAL();
    for (size_t k = 0; k < splits.size(); ++k) {
      const auto &cs = splits[k];
      REAL t = _dlclvs[cs.left][ec] * _dlclvs[cs.right][ec] *
               this->synWgdFreq(cid, k, cs.frequency);
      scale(t);
      dAcc = dAcc + t;
      if (dAcc > dSample) {
        ev.type = UEventSample::RESOLVE_SPLIT;
        ev.left = cs.left;
        ev.right = cs.right;
        ev.blLeft = cs.blLeft;
        ev.blRight = cs.blRight;
        return true;
      }
    }
    ev.type = UEventSample::RESOLVE_SINGLE; // 2*E*R: one copy lost
    return true;
  }
  if (isLeaf) {
    ev.type = UEventSample::LEAF;
    return true;
  }
  for (const auto &cs : splits) {
    double w0 = (1.0 - r) * psU * cs.frequency;
    REAL s1 = _uclvs[cs.left][fc] * _uclvs[cs.right][gc] * w0;
    scale(s1);
    acc = acc + s1;
    if (acc > toSample) {
      ev.type = UEventSample::SPECIATE;
      ev.left = cs.left;
      ev.right = cs.right;
      ev.blLeft = cs.blLeft;
      ev.blRight = cs.blRight;
      ev.leftToLeftSpecies = true;
      return true;
    }
    REAL s2 = _uclvs[cs.right][fc] * _uclvs[cs.left][gc] * w0;
    scale(s2);
    acc = acc + s2;
    if (acc > toSample) {
      ev.type = UEventSample::SPECIATE;
      ev.left = cs.left;
      ev.right = cs.right;
      ev.blLeft = cs.blLeft;
      ev.blRight = cs.blRight;
      ev.leftToLeftSpecies = false;
      return true;
    }
  }
  double w1 = (1.0 - r) * psU;
  REAL sl1 = _uclvs[cid][fc] * _uEU[gc] * w1;
  scale(sl1);
  acc = acc + sl1;
  ev.type = UEventSample::SPECIATE_LOSS;
  if (acc > toSample) {
    ev.dest = this->getSpeciesLeft(sp); // keep f, lose g
    ev.lost = this->getSpeciesRight(sp);
  } else {
    ev.dest = this->getSpeciesRight(sp); // keep g, lose f
    ev.lost = this->getSpeciesLeft(sp);
  }
  return true;
}

// Unresolved-state backtrace for the resolution profile: one U draw, recording
// a U->R commit (species branch) when resolution fires.
template <class REAL>
void UndatedDLMultiModel<REAL>::btU(CID cid, corax_rnode_t *sp, unsigned int c,
                                    std::vector<unsigned int> &commits,
                                    bool check) {
  UEventSample ev;
  if (!sampleUEventImpl(cid, sp, c, ev, check)) {
    return;
  }
  auto e = sp->node_index;
  switch (ev.type) {
  case UEventSample::RESOLVE_SPLIT:
    commits.push_back(e);
    btR(ev.left, sp, c, commits, check);
    btR(ev.right, sp, c, commits, check);
    return;
  case UEventSample::RESOLVE_SINGLE:
    commits.push_back(e);
    btR(cid, sp, c, commits, check);
    return;
  case UEventSample::LEAF:
    _tetra.push_back(e); // still tetrasomic at the tip
    return;
  case UEventSample::SPECIATE:
    if (ev.leftToLeftSpecies) {
      btU(ev.left, this->getSpeciesLeft(sp), c, commits, check);
      btU(ev.right, this->getSpeciesRight(sp), c, commits, check);
    } else {
      btU(ev.right, this->getSpeciesLeft(sp), c, commits, check);
      btU(ev.left, this->getSpeciesRight(sp), c, commits, check);
    }
    return;
  case UEventSample::SPECIATE_LOSS:
    btU(cid, ev.dest, c, commits, check);
    return;
  }
}
