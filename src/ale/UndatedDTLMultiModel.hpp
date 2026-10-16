#pragma once

#include <cstdlib> // std::abort for the always-on LORe consistency guard

#include <trees/DatedTree.hpp>

#include "MultiModel.hpp"

/**
 *  Highway with the proba scaled per rate category
 */
struct WeightedHighway {
  Highway highway;
  std::vector<double> proba;
};

/**
 *  Implements all HGT modelling-dependent functions of
 *  the MultiModel class in the context of the UndatedDTL
 *  model
 */
template <class REAL> class UndatedDTLMultiModel : public MultiModel<REAL> {
public:
  UndatedDTLMultiModel(DatedTree &speciesTree,
                       const GeneSpeciesMapping &geneSpeciesMapping,
                       const RecModelInfo &info, const std::string &ccpFile);

  virtual ~UndatedDTLMultiModel() {}

  virtual void setAlpha(double alpha);
  virtual void setRates(const RatesVector &rates);
  virtual void setHighways(const std::vector<Highway> &highways);

  // --- WGD extension ---
  // Declare a WGD at the top of species branch e with retention prob q.
  void setWGD(unsigned int e, double q) override {
    _hasWGD[e] = 1;
    _q[e] = q;
    this->invalidateAllSpeciesNodes();
    this->resetCache();
  }

  // --- LORe (delayed rediploidization) extension ---
  // Set the global resolution probability r in [0,1] on EVERY species branch.
  // r = 1 recovers the WHALE/AORe model bit-for-bit (the nested null). The
  // U-state propagates only vertically (speciation / duplication-bracket /
  // loss); transfers and highways act on resolved lineages and read raw
  // (post-WGD-below) quantities, exactly as in the AORe DTL model.
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
  void setResolutionProbBranch(unsigned int e, double r) {
    _resolutionProbs[e] = r;
    this->invalidateAllSpeciesNodes();
    this->resetCache();
  }
  double getResolutionProbBranch(unsigned int e) const {
    return _resolutionProbs[e];
  }

  // --- LORe resolution-branch marginal ---
  // Sample `samples` resolution histories under the CURRENT (fitted) global r
  // via a U-aware backtrace, accumulating per species branch the expected
  // number of U->R commit (ohnolog-divergence) events. CLVs must be current.
  void sampleResolutionCommits(unsigned int samples,
                               std::vector<double> &commitCounts,
                               std::vector<double> &tetraCounts, bool check);
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
  bool sampleUEventImpl(CID cid, corax_rnode_t *sp, unsigned int c,
                        UEventSample &ev, bool check);
  // U-aware backtrace helpers for the resolution-branch marginal. They reuse
  // computeProbability() to sample the resolved (R-state) events (S/D/T/SL/DL/
  // TL/highways), and add the unresolved (U-state, vertical-only) recursion and
  // the WGD R/U coin on top.
  corax_rnode_t *sampleOriginationR(unsigned int &category);
  void btR(CID cid, corax_rnode_t *sp, unsigned int c,
           std::vector<unsigned int> &commits, bool check);
  void btU(CID cid, corax_rnode_t *sp, unsigned int c,
           std::vector<unsigned int> &commits, bool check);
  void btDescend(CID cid, corax_rnode_t *sp, unsigned int c,
                 std::vector<unsigned int> &commits, bool check);

  DatedTree &_datedTree;
  unsigned int _gammaCatNumber;
  std::vector<double> _gammaScalers;
  RatesVector _dtlRates;
  // Highways, per species branch
  std::vector<std::vector<WeightedHighway>> _highways;
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PT; // Transfer probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _OP; // Origination probability, per species branch
  std::vector<REAL> _uE;   // Extinction probability, per species branch
  std::vector<REAL> _tE; // Transfer-extinction probability, per species branch
  // --- WGD extension ---
  std::vector<char> _hasWGD; // per species branch: WGD at top of branch?
  std::vector<double> _q;    // per species branch: retention prob (if WGD)
  std::vector<REAL> _uEtop;  // extinction seen by the PARENT, per (e,cat)
  // --- LORe (delayed rediploidization) extension ---
  std::vector<REAL> _uEU;          // unresolved (tetrasomic) extinction, per (e,cat)
  std::vector<double> _resolutionProbs; // per species branch: r in [0,1]; 1==AORe
  // resolution-backtrace: number of lineages abandoned at the resample cap (a
  // degenerate cell where only DL/TL "nothing observed" events have mass).
  size_t _btCapHits = 0;
  // hard per-sample backstop: total backtrace steps (across btR/btU/btDescend).
  // Reset before each sample; once exceeded, every call returns immediately so
  // the recursion unwinds and the (pathological) sample is abandoned. Bounds any
  // runaway, whether an iterative loop or an inter-function re-entry.
  size_t _btSteps = 0;
  static const size_t kBtStepCap = 50000;
  // Leaf species at which a U (unresolved/tetrasomic) lineage reaches an extant
  // tip without ever committing -> still tetrasomic in that species today.
  // Filled by btU during a resolution-commit sample; read+cleared per sample in
  // sampleResolutionCommits (parallel to the resolved-state `commits`).
  std::vector<unsigned int> _tetra;
  OriginationStrategy _originationStrategy;
  TransferConstaint _transferConstraint;
  // Forbidden transfers, per species branch
  std::vector<std::vector<corax_rnode_t *>> _transferForbiddenSpeciesNodes;
  // Allowed transfers, per species branch
  std::vector<std::vector<corax_rnode_t *>> _transferCandidateSpeciesNodes;

  struct DTLCLV {
    // Element e of the gene clade's _uq stores the probability of the clade,
    // given the clade is mapped to the species branch e.
    // In the paper: Pi_{e,gamma} of a clade gamma for each branch e
    std::vector<REAL> _uq;
    // Element e of the gene clade's _tq stores the probability of the clade
    // after its transfer from the species branch e to some other branch.
    // In the paper: \bar{Pi}_{e,gamma} of a clade gamma for each branch e
    std::vector<REAL> _tq;
    // Post-WGD clade CLV seen by the PARENT: equals _uq on branches without a
    // WGD, and the WHALE/LORe retention transform of _uq on branches that carry
    // one ( (1-q) _uq + q _uu ; _uu == doubling bracket when r=1 -> WHALE).
    std::vector<REAL> _uqTop;
    // Unresolved (tetrasomic) clade CLV U: a single lineage still holding a
    // not-yet-diverged ohnolog pair. Propagated vertically only (no transfers).
    std::vector<REAL> _uu;
    DTLCLV() : _uq(0), _tq(0), _uqTop(0), _uu(0) {}
    DTLCLV(unsigned int speciesNumber, unsigned int gammaCategories)
        : _uq(speciesNumber * gammaCategories, REAL()),
          _tq(speciesNumber * gammaCategories, REAL()),
          _uqTop(speciesNumber * gammaCategories, REAL()),
          _uu(speciesNumber * gammaCategories, REAL()) {}
  };
  // vector of DTLCLVs for all observed clades
  std::vector<DTLCLV> _dtlclvs;

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

  // functions to deal with transfers
  void updateTransferCandidates();
  bool sampleTransferEvent(CID cid, corax_rnode_t *srcSpeciesNode,
                           unsigned int category, Scenario::Event &event);

  double getTransferWeightNorm(unsigned int speciesNodeIndex) const {
    auto e = speciesNodeIndex;
    return static_cast<double>(_transferCandidateSpeciesNodes[e].size());
  }

  // Doubling bracket shared by the WGD and LORe transforms:
  //   2*E*R + sum_splits ccp * R_left * R_right   (no _PD, no _PT, no q/r).
  // Reads R (=_uq) and E (=_uE) at node-cat ec; _uq[cid][ec] must already hold R
  // for this clade. Transfers do not enter (the U-state is vertical-only).
  REAL dupBracket(unsigned int cid, unsigned int ec) {
    REAL split = REAL();
    const auto &splits = this->_ccp.getCladeSplits(cid);
    for (size_t k = 0; k < splits.size(); ++k) {
      const auto &cs = splits[k];
      // WGD split: synteny-weighted frequency (== ccp frequency without synteny)
      REAL t = _dtlclvs[cs.left]._uq[ec] * _dtlclvs[cs.right]._uq[ec] *
               this->synWgdFreq(cid, k, cs.frequency);
      scale(t);
      split += t;
    }
    REAL out = _uE[ec] * _dtlclvs[cid]._uq[ec] * 2.0 + split;
    scale(out);
    return out;
  }

  // Normaliser of the unresolved (U) event probabilities: a tetrasomic U
  // lineage is vertical-only (no D, T or highway channel), so its speciation
  // and loss probabilities are renormalised over PS + PL to sum to 1.
  double uDen(unsigned int ec) const { return _PS[ec] + _PL[ec]; }

  // TEMP debug guard for the resolution backtrace: validate (sp, cid) before any
  // pointer deref / array index, aborting with a clear message instead of
  // segfaulting. Used to localise the residual DTL export crash.
  void btGuard(const char *where, CID cid, corax_rnode_t *sp) {
    bool bad = false;
    if (!sp) {
      std::cerr << "[BTDBG] NULL sp in " << where << " cid=" << cid << std::endl;
      bad = true;
    } else if (sp->node_index >= this->getAllSpeciesNodeNumber()) {
      std::cerr << "[BTDBG] OOB node_index=" << sp->node_index
                << " (N=" << this->getAllSpeciesNodeNumber() << ") in " << where
                << " cid=" << cid << std::endl;
      bad = true;
    }
    if (cid >= _dtlclvs.size()) {
      std::cerr << "[BTDBG] OOB cid=" << cid << " (nClades=" << _dtlclvs.size()
                << ") in " << where << std::endl;
      bad = true;
    }
    if (bad) {
      std::cerr.flush();
      std::abort();
    }
  }

  // functions to work with _llCache
  virtual size_t getHash() {
    auto hash = this->getSpeciesTreeHash();
    return (_transferConstraint == TransferConstaint::RELDATED)
               ? _datedTree.getOrderingHash(hash)
               : hash;
  }
};

/**
 *  Constructor
 */
template <class REAL>
UndatedDTLMultiModel<REAL>::UndatedDTLMultiModel(
    DatedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMapping,
    const RecModelInfo &info, const std::string &ccpFile)
    : MultiModel<REAL>(speciesTree.getRootedTree(), geneSpeciesMapping, info,
                       ccpFile),
      _datedTree(speciesTree), _gammaCatNumber(info.gammaCategories),
      _gammaScalers(_gammaCatNumber, 1.0),
      _PD(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 0.2),
      _PL(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 0.2),
      _PT(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 0.1),
      _PS(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 1.0),
      _OP(this->getAllSpeciesNodeNumber(),
          1.0 / static_cast<double>(this->getAllSpeciesNodeNumber())),
      _uE(this->getAllSpeciesNodeNumber() * _gammaCatNumber, REAL()),
      _tE(this->getAllSpeciesNodeNumber() * _gammaCatNumber, REAL()),
      _originationStrategy(info.originationStrategy),
      _transferConstraint(info.transferConstraint) {
  auto N = this->getAllSpeciesNodeNumber();
  // --- WGD extension: per-branch state, no WGD by default ---
  _hasWGD.assign(N, 0);
  _q.assign(N, 1.0);
  _uEtop.assign(N * _gammaCatNumber, REAL());
  // --- LORe extension: unresolved extinction + per-branch resolution prob ---
  _uEU.assign(N * _gammaCatNumber, REAL());
  _resolutionProbs.assign(N, 1.0); // default r=1 everywhere == AORe/WHALE
  // set gamma scalers with the default alpha
  setAlpha(1.0);
  // set all DTLO rates to the default value
  _dtlRates.resize(this->_info.modelFreeParameters(),
                   std::vector<double>(N, 0.2));
  // initialize highways and per-species transfer candidates
  _highways.resize(N);
  _transferForbiddenSpeciesNodes.resize(N);
  _transferCandidateSpeciesNodes.resize(N);
  // initialize DTLCLVs if needed
  if (!this->_memorySavings) {
    allocateMemory();
  }
}

template <class REAL> void UndatedDTLMultiModel<REAL>::setAlpha(double alpha) {
  corax_compute_gamma_cats(alpha, _gammaScalers.size(), &_gammaScalers[0],
                           CORAX_GAMMA_RATES_MEAN);
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::setRates(const RatesVector &rates) {
  assert(rates.size() == this->_info.modelFreeParameters());
  _dtlRates = rates;
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::setHighways(
    const std::vector<Highway> &highways) {
  for (auto &speciesWeightedHighways : _highways) {
    speciesWeightedHighways.clear();
  }
  for (auto highway : highways) {
    WeightedHighway hp;
    hp.highway = highway;
    // map the highway to the pruned species tree
    if (this->prunedMode()) {
      hp.highway.src = this->_speciesToPrunedNode[highway.src->node_index];
      hp.highway.dest = this->_speciesToPrunedNode[highway.dest->node_index];
    } else {
      hp.highway.src = highway.src;
      hp.highway.dest = highway.dest;
    }
    // in the pruned mode a highway from/to a species not covered by the gene
    // family may have the src/dest species absent or look like a self-transfer
    if (!hp.highway.src || !hp.highway.dest ||
        (hp.highway.src == hp.highway.dest)) {
      continue; // this highway should not affect this gene family
    }
    // this value will be normalized later on, separately per category
    hp.proba.resize(_gammaCatNumber, highway.proba);
    // update _highways
    _highways[hp.highway.src->node_index].push_back(hp);
  }
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

/**
 *  Allocate memory to the CLVs
 */
template <class REAL> void UndatedDTLMultiModel<REAL>::allocateMemory() {
  DTLCLV nullCLV(this->getAllSpeciesNodeNumber(), _gammaCatNumber);
  _dtlclvs = std::vector<DTLCLV>(this->_ccp.getCladesNumber(), nullCLV);
}

/**
 *  Free memory allocated to the CLVs
 */
template <class REAL> void UndatedDTLMultiModel<REAL>::deallocateMemory() {
  _dtlclvs = std::vector<DTLCLV>();
}

/**
 *  Compute the CLV for a given clade
 */
template <class REAL> void UndatedDTLMultiModel<REAL>::updateCLV(CID cid) {
  auto &uq = _dtlclvs[cid]._uq;
  auto &tq = _dtlclvs[cid]._tq;
  auto &uqTop = _dtlclvs[cid]._uqTop;
  auto &uu = _dtlclvs[cid]._uu;
  std::fill(uq.begin(), uq.end(), REAL());
  std::fill(tq.begin(), tq.end(), REAL());
  std::fill(uqTop.begin(), uqTop.end(), REAL());
  std::fill(uu.begin(), uu.end(), REAL());
  // iterate several times to resolve the DL and TL terms with
  // fixed point optimization: not needed if we don't model TL
  unsigned int maxIt = this->_info.noTL ? 1 : 4;
  for (unsigned int it = 0; it < maxIt; ++it) {
    bool ok;
    for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
      REAL transferSum = REAL();
      // postorder species tree traversal is granted
      for (auto speciesNode : this->getPrunedSpeciesNodes()) {
        auto e = speciesNode->node_index;
        auto ec = e * _gammaCatNumber + c;
        REAL p = REAL();
        ok = computeProbability(cid, speciesNode, c, p);
        assert(ok);
        transferSum += p;
        uq[ec] = p; // raw P_e(gamma)

        // --- LORe: unresolved (tetrasomic) CLV U (vertical-only) ---
        //   U = r * dupBracket(cid,ec)                 // resolve here -> doubling
        //     + (1-r) * [ speciation, both daughters inherit U (both orderings)
        //                 + speciation with unresolved loss of one daughter ]
        //   leaf node: U = r*dupBracket + (1-r)*R  (never resolved -> single gene)
        // No transfer channel: a still-tetrasomic locus is not transferred and
        // is not a transfer destination. r=1 => U == dupBracket => WHALE/AORe.
        {
          const double r = _resolutionProbs[e];
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
                REAL t = (_dtlclvs[cs.left]._uu[fc] * _dtlclvs[cs.right]._uu[gc] +
                          _dtlclvs[cs.right]._uu[fc] * _dtlclvs[cs.left]._uu[gc]) *
                         (_PS[ec] / uDen(ec) * cs.frequency);
                scale(t);
                rest += t;
              }
              // speciation + unresolved loss of one daughter (SL acting on U)
              REAL sl = (uu[fc] * _uEU[gc] + uu[gc] * _uEU[fc]) * (_PS[ec] / uDen(ec));
              scale(sl);
              rest += sl;
            } else {
              rest = p; // leaf: never resolved, seen as the single leaf gene
            }
            rest *= (1.0 - r);
            scale(rest);
            uval += rest;
            scale(uval);
          }
          uu[ec] = uval;
        }

        // --- clade CLV seen by the parent (WHALE/LORe retention transform) ---
        // P_top = (1-q) P + q U ; q=0 recovers no-WGD, and U==dupBracket==
        // 2 E P + split when r=1, so this reduces to the AORe/WHALE transform.
        // Uses raw own-branch quantities; transfer reads stay raw because a
        // transferred lineage lands below the WGD.
        if (_hasWGD[e]) {
          double q = _q[e];
          REAL keep = p * (1.0 - q); // (1-q) P
          REAL inj = uu[ec] * q;     // q U
          scale(keep);
          scale(inj);
          uqTop[ec] = keep + inj;
        } else {
          uqTop[ec] = p;
        }
      }
      // now that we've got the clade proba on every species branch,
      // we can compute the clade transfer probas
      for (auto speciesNode : this->getPrunedSpeciesNodes()) {
        auto e = speciesNode->node_index;
        auto ec = e * _gammaCatNumber + c;
        REAL tp = REAL();
        for (auto destSpeciesNode : _transferForbiddenSpeciesNodes[e]) {
          auto d = destSpeciesNode->node_index;
          auto dc = d * _gammaCatNumber + c;
          tp += uq[dc];
        }
        tp = transferSum - tp; // sum over the candidates, but faster
        tp /= getTransferWeightNorm(e);
        scale(tp);
        tq[ec] = tp;
      }
    }
  }
}

/**
 *  Update the list of the allowed transfer-receiving species
 *  for each species branch
 */
template <class REAL>
void UndatedDTLMultiModel<REAL>::updateTransferCandidates() {
  std::function<bool(unsigned int, unsigned int)> canTransfer;
  switch (_transferConstraint) {
  case TransferConstaint::NONE:
    // include all the species nodes except for the self
    canTransfer = [](unsigned int e, unsigned int d) { return d != e; };
    break;
  case TransferConstaint::PARENTS:
    // include all the species nodes except for the parents
    canTransfer = [this](unsigned int e, unsigned int d) {
      return !this->_speciesTree.isAncestorOf(d, e);
    };
    break;
  case TransferConstaint::RELDATED:
    // include all the species nodes younger than the self's parent
    canTransfer = [this](unsigned int e, unsigned int d) {
      return this->_datedTree.canTransferUnderRelDated(e, d);
    };
    break;
  }
  // clear for all species nodes
  for (unsigned int i = 0; i < this->getAllSpeciesNodeNumber(); ++i) {
    _transferForbiddenSpeciesNodes[i].clear();
    _transferCandidateSpeciesNodes[i].clear();
  }
  // fill for the current pruned species nodes
  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    for (auto destSpeciesNode : this->getPrunedSpeciesNodes()) {
      auto d = destSpeciesNode->node_index;
      if (canTransfer(e, d)) {
        _transferCandidateSpeciesNodes[e].push_back(destSpeciesNode);
      } else {
        _transferForbiddenSpeciesNodes[e].push_back(destSpeciesNode);
      }
    }
  }
}

/**
 *  Compute the per species branch probabilities of
 *  the elementary events of clade evolution
 */
template <class REAL>
void UndatedDTLMultiModel<REAL>::recomputeSpeciesProbabilities() {
  auto allSpeciesNumber = this->getAllSpeciesNodeNumber();
  // recompute _PD, _PL, _PT, _PS and highway.proba
  auto &dupRates = _dtlRates[0];
  auto &lossRates = _dtlRates[1];
  auto &transferRates = _dtlRates[2];
  assert(allSpeciesNumber == dupRates.size());
  assert(allSpeciesNumber == lossRates.size());
  assert(allSpeciesNumber == transferRates.size());
  std::fill(_PD.begin(), _PD.end(), 0.0);
  std::fill(_PL.begin(), _PL.end(), 0.0);
  std::fill(_PT.begin(), _PT.end(), 0.0);
  std::fill(_PS.begin(), _PS.end(), 0.0);
  for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      auto e = speciesNode->node_index;
      auto ec = e * _gammaCatNumber + c;
      _PD[ec] = dupRates[e];
      _PL[ec] = lossRates[e];
      _PT[ec] = transferRates[e];
      _PS[ec] = _gammaScalers[c];
      if (this->_info.noDup) {
        assert(!this->_info.noTL);
        _PD[ec] = 0.0;
      }
      auto sum = _PD[ec] + _PL[ec] + _PT[ec] + _PS[ec];
      for (const auto &highway : _highways[e]) {
        sum += highway.highway.proba;
      }
      _PD[ec] /= sum;
      _PL[ec] /= sum;
      _PT[ec] /= sum;
      _PS[ec] /= sum;
      for (auto &highway : _highways[e]) {
        assert(highway.highway.proba >= 0.0);
        // highway proba, normalized per category
        highway.proba[c] = highway.highway.proba / sum;
        assert(highway.proba[c] < 1.0);
      }
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
    auto &oriRates = _dtlRates[3];
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
  // recompute _uE and _tE
  updateTransferCandidates();
  std::fill(_uE.begin(), _uE.end(), REAL());
  std::fill(_tE.begin(), _tE.end(), REAL());
  std::fill(_uEtop.begin(), _uEtop.end(), REAL());
  std::fill(_uEU.begin(), _uEU.end(), REAL());
  // iterate several times to resolve _uE and _tE probas with
  // fixed point optimization
  unsigned int maxIt = 4;
  for (unsigned int it = 0; it < maxIt; ++it) {
    for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
      REAL extinctionSum = REAL();
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
        // TEE scenario
        temp = _uE[ec] * _tE[ec] * _PT[ec];
        scale(temp);
        proba += temp;
        // highway TEE scenario
        for (const auto &highway : _highways[e]) {
          auto d = highway.highway.dest->node_index;
          auto dc = d * _gammaCatNumber + c;
          temp = _uE[ec] * _uE[dc] * highway.proba[c];
          scale(temp);
          proba += temp;
        }
        assert(proba < REAL(1.000001));
        extinctionSum += proba;
        _uE[ec] = proba;
        // --- LORe: unresolved (tetrasomic) extinction EU (vertical-only) ---
        //   EU = r*E^2 + (1-r)*( _PL + _PS*EU[fc]*EU[gc] )    (internal)
        //   EU = r*E^2 + (1-r)*E                              (leaf)
        // No transfer-extinction channel: a still-tetrasomic lineage does not
        // transfer. EU has no self-dependence, so it rides on the E fixed point
        // in a single postorder pass. r=1 => EU == E^2 (the WHALE doubling term).
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
        // --- post-WGD extinction seen by the parent (WHALE/LORe form) ---
        // E_top = (1-q) E + q EU  (== (1-q)E + q E^2 when r=1 -> WHALE/AORe).
        // The transfer extinction _tE keeps reading the raw _uE (a transferred
        // lineage lands within the destination branch, below its WGD).
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
      // now that we've got extinction probas for every species branch,
      // we can compute transfer-extinction probas
      for (auto speciesNode : this->getPrunedSpeciesNodes()) {
        auto e = speciesNode->node_index;
        auto ec = e * _gammaCatNumber + c;
        REAL tproba = REAL();
        for (auto destSpeciesNode : _transferForbiddenSpeciesNodes[e]) {
          auto d = destSpeciesNode->node_index;
          auto dc = d * _gammaCatNumber + c;
          tproba += _uE[dc];
        }
        tproba = extinctionSum - tproba; // sum over the candidates, but faster
        tproba /= getTransferWeightNorm(e);
        scale(tproba);
        assert(tproba < REAL(1.000001));
        _tE[ec] = tproba;
      }
    }
  } // end of iteration
}

/**
 *  Correction factor to the species tree likelihood,
 *  because we condition on survival
 */
template <class REAL>
double UndatedDTLMultiModel<REAL>::getLikelihoodFactor(unsigned int category) {
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
REAL UndatedDTLMultiModel<REAL>::getRootCladeLikelihood(
    corax_rnode_t *speciesNode, unsigned int category) {
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto c = category;
  auto e = speciesNode->node_index;
  auto ec = e * _gammaCatNumber + c;
  REAL likelihood = _dtlclvs[rootCID]._uq[ec] * _OP[e];
  scale(likelihood);
  return likelihood;
}

/**
 *  Sample a transfer destination species branch and
 *  write it into the given recCell event
 */
template <class REAL>
bool UndatedDTLMultiModel<REAL>::sampleTransferEvent(
    unsigned int cid, corax_rnode_t *srcSpeciesNode, unsigned int category,
    Scenario::Event &event) {
  auto c = category;
  auto e = srcSpeciesNode->node_index;
  auto ec = e * _gammaCatNumber + c;
  auto survivingTransferSum = _dtlclvs[cid]._tq[ec] * getTransferWeightNorm(e);
  auto toSample = this->getRandom(survivingTransferSum);
  auto sumProba = REAL();
  for (auto destSpeciesNode : _transferCandidateSpeciesNodes[e]) {
    auto d = destSpeciesNode->node_index;
    auto dc = d * _gammaCatNumber + c;
    sumProba += _dtlclvs[cid]._uq[dc];
    if (sumProba > toSample) {
      event.destSpeciesNode = d;
      event.pllDestSpeciesNode = destSpeciesNode;
      return true;
    }
  }
  return false;
}

/**
 *  Compute the CLV value for a given cid (clade id) and a given
 *  species node and write it to the proba variable
 */
template <class REAL>
bool UndatedDTLMultiModel<REAL>::computeProbability(
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
  // S events on internal species branches, D events and T events can happen
  // for ancestral gene nodes only:
  for (const auto &cladeSplit : this->_ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left;
    auto cidRight = cladeSplit.right;
    auto freq = cladeSplit.frequency;
    // - S event on an internal species branch
    if (!isSpeciesLeaf) {
      temp = _dtlclvs[cidLeft]._uqTop[fc] * _dtlclvs[cidRight]._uqTop[gc] *
             (_PS[ec] * freq);
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
      temp = _dtlclvs[cidRight]._uqTop[fc] * _dtlclvs[cidLeft]._uqTop[gc] *
             (_PS[ec] * freq);
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
    temp = _dtlclvs[cidLeft]._uq[ec] * _dtlclvs[cidRight]._uq[ec] *
           (_PD[ec] * freq);
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
    // - T event
    temp = _dtlclvs[cidLeft]._uq[ec] *
           (_dtlclvs[cidRight]._tq[ec] * (_PT[ec] * freq));
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      if (!sampleTransferEvent(cidRight, speciesNode, c, recCell->event)) {
        return false;
      }
      recCell->event.leftGeneIndex = cidLeft;
      recCell->event.rightGeneIndex = cidRight;
      recCell->blLeft = cladeSplit.blLeft;
      recCell->blRight = cladeSplit.blRight;
      return true;
    }
    temp = _dtlclvs[cidRight]._uq[ec] *
           (_dtlclvs[cidLeft]._tq[ec] * (_PT[ec] * freq));
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      if (!sampleTransferEvent(cidLeft, speciesNode, c, recCell->event)) {
        return false;
      }
      recCell->event.leftGeneIndex = cidRight;
      recCell->event.rightGeneIndex = cidLeft;
      recCell->blLeft = cladeSplit.blRight;
      recCell->blRight = cladeSplit.blLeft;
      return true;
    }
    // - highway T event
    for (const auto &highway : _highways[e]) {
      auto d = highway.highway.dest->node_index;
      auto dc = d * _gammaCatNumber + c;
      temp = _dtlclvs[cidLeft]._uq[ec] *
             (_dtlclvs[cidRight]._uq[dc] * (highway.proba[c] * freq));
      scale(temp);
      proba += temp;
      if (proba > REAL(1.0)) {
        std::cerr << "error " << _dtlclvs[cidLeft]._uq[ec] << " "
                  << _dtlclvs[cidRight]._uq[dc] << " " << highway.proba[c]
                  << " " << freq << std::endl;
      }
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_T;
        recCell->event.destSpeciesNode = d;
        recCell->event.pllDestSpeciesNode = highway.highway.dest;
        recCell->event.leftGeneIndex = cidLeft;
        recCell->event.rightGeneIndex = cidRight;
        recCell->blLeft = cladeSplit.blLeft;
        recCell->blRight = cladeSplit.blRight;
        return true;
      }
      temp = _dtlclvs[cidRight]._uq[ec] *
             (_dtlclvs[cidLeft]._uq[dc] * (highway.proba[c] * freq));
      scale(temp);
      proba += temp;
      if (proba > REAL(1.0)) {
        std::cerr << "error " << _dtlclvs[cidRight]._uq[ec] << " "
                  << _dtlclvs[cidLeft]._uq[dc] << " " << highway.proba[c] << " "
                  << freq << std::endl;
      }
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_T;
        recCell->event.destSpeciesNode = d;
        recCell->event.pllDestSpeciesNode = highway.highway.dest;
        recCell->event.leftGeneIndex = cidRight;
        recCell->event.rightGeneIndex = cidLeft;
        recCell->blLeft = cladeSplit.blRight;
        recCell->blRight = cladeSplit.blLeft;
        return true;
      }
    }
  }
  // SL events, DL events and TL events can happen
  // for any of gene nodes:
  // - SL event (only on an internal species branch)
  if (!isSpeciesLeaf) {
    temp = _dtlclvs[cid]._uqTop[fc] * (_uEtop[gc] * _PS[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.lostSpeciesNode = g;
      recCell->event.pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      recCell->event.pllLostSpeciesNode = this->getSpeciesRight(speciesNode);
      return true;
    }
    temp = _dtlclvs[cid]._uqTop[gc] * (_uEtop[fc] * _PS[ec]);
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
  if (this->_info.noTL) { // no TL events allowed
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
  } else { // TL events allowed
    if (!this->_info.noDL) {
      // - DL event
      temp = _dtlclvs[cid]._uq[ec] * (_uE[ec] * _PD[ec] * 2.0);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        // in fact, nothing happens, we'll have to resample
        recCell->event.type = ReconciliationEventType::EVENT_DL;
        return true;
      }
    }
    // - TL event
    // we transfer, but the gene gets extinct in the receiving species
    temp = _dtlclvs[cid]._uq[ec] * (_tE[ec] * _PT[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      // in fact, nothing happens, we'll have to resample
      recCell->event.type = ReconciliationEventType::EVENT_TL;
      recCell->event.pllDestSpeciesNode = nullptr;
      return true;
    }
    // we transfer, and the gene gets extinct in the sending species
    temp = _dtlclvs[cid]._tq[ec] * (_uE[ec] * _PT[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_TL;
      if (!sampleTransferEvent(cid, speciesNode, c, recCell->event)) {
        return false;
      }
      return true;
    }
    // - highway TL event
    for (const auto &highway : _highways[e]) {
      auto d = highway.highway.dest->node_index;
      auto dc = d * _gammaCatNumber + c;
      // we transfer, but the gene gets extinct in the receiving species
      temp = _dtlclvs[cid]._uq[ec] * (_uE[dc] * highway.proba[c]);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        // in fact, nothing happens, we'll have to resample
        recCell->event.type = ReconciliationEventType::EVENT_TL;
        recCell->event.pllDestSpeciesNode = nullptr;
        return true;
      }
      // we transfer, and the gene gets extinct in the sending species
      temp = _dtlclvs[cid]._uq[dc] * (_uE[ec] * highway.proba[c]);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_TL;
        recCell->event.destSpeciesNode = d;
        recCell->event.pllDestSpeciesNode = highway.highway.dest;
        return true;
      }
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
 *  LORe resolution-branch marginal (DTL port)
 *
 *  A U-aware backtrace. computeProbability() already samples the resolved
 *  (R-state) events (S/D/T/SL/DL/TL and highways) proportionally to their inside
 *  contribution. On top of that we add the WGD R-vs-U coin when a lineage is
 *  inherited across the WGD branch, and the unresolved (U-state, vertical-only)
 *  recursion, which emits a U->R commit (records the species branch) when
 *  resolution fires. Transfers move resolved lineages only (a transferred copy
 *  lands resolved, below any WGD at its destination), so they never enter the U
 *  recursion. All sampling weights are exact terms of the inside recursion, so
 *  the marginals are unbiased; `check` asserts weight==inside at every U cell
 *  and at the WGD R/U coin.
 */

// helper: |log a - log b| small (robust to ScaledValue scaling); skips null.
// Distinct name from the DL helper to avoid a redefinition when both model
// headers are included in the same translation unit.
template <class REAL>
static inline bool dtlLoreInsideClose(const REAL &a, const REAL &b) {
  double la = getLog(a);
  double lb = getLog(b);
  if (!std::isfinite(la) && !std::isfinite(lb)) {
    return true; // both effectively zero
  }
  return std::fabs(la - lb) < 1e-6;
}

template <class REAL>
corax_rnode_t *
UndatedDTLMultiModel<REAL>::sampleOriginationR(unsigned int &category) {
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
void UndatedDTLMultiModel<REAL>::sampleResolutionCommits(
    unsigned int samples, std::vector<double> &commitCounts,
    std::vector<double> &tetraCounts, bool check) {
  commitCounts.assign(this->getAllSpeciesNodeNumber(), 0.0);
  tetraCounts.assign(this->getAllSpeciesNodeNumber(), 0.0);
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  for (unsigned int s = 0; s < samples; ++s) {
    unsigned int cat = 0;
    auto origin = sampleOriginationR(cat);
    // A lineage that ORIGINATES anywhere starts resolved (the root term uses
    // raw _uq, not the WGD-transformed _uqTop).
    std::vector<unsigned int> commits;
    _tetra.clear();
    _btSteps = 0; // reset the hard per-sample backstop
    btR(rootCID, origin, cat, commits, check);
    for (auto b : commits) {
      commitCounts[b] += 1.0;
    }
    for (auto b : _tetra) {
      tetraCounts[b] += 1.0; // leaf lineages still tetrasomic at the tip
    }
  }
}

// Descend into a child species node; sample the WGD R-vs-U coin if it is the
// WGD branch, otherwise continue resolved.
template <class REAL>
bool UndatedDTLMultiModel<REAL>::sampleWGDCoin(CID cid, corax_rnode_t *sp,
                                               unsigned int c) {
  auto e = sp->node_index;
  if (!_hasWGD[e]) {
    return false;
  }
  auto ec = e * _gammaCatNumber + c;
  double q = _q[e];
  REAL wR = _dtlclvs[cid]._uq[ec] * (1.0 - q); // stay R
  REAL wU = _dtlclvs[cid]._uu[ec] * q;         // become U
  scale(wR);
  scale(wU);
  REAL total = wR + wU;
  REAL toSample = this->getRandom(total);
  return !(wR > toSample);
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::btDescend(CID cid, corax_rnode_t *sp,
                                           unsigned int c,
                                           std::vector<unsigned int> &commits,
                                           bool check) {
  if (++_btSteps > kBtStepCap) {
    return; // hard per-sample backstop
  }
  btGuard("btDescend", cid, sp);
  auto e = sp->node_index;
  if (check && _hasWGD[e]) {
    auto ec = e * _gammaCatNumber + c;
    REAL wR = _dtlclvs[cid]._uq[ec] * (1.0 - _q[e]);
    REAL wU = _dtlclvs[cid]._uu[ec] * _q[e];
    scale(wR);
    scale(wU);
    if (!dtlLoreInsideClose(REAL(wR + wU), _dtlclvs[cid]._uqTop[ec])) {
      std::cerr << "LORe DTL coin weights != _uqTop at cid=" << cid
                << " e=" << e << std::endl;
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
//
// The "nothing observed, resample" events (DL, and TL) do not consume a gene-
// clade split, so they are handled by LOOPING (re-sampling the same cell, or
// the same gene at the transfer destination) rather than recursing — recursing
// on them is unbounded and overflows the stack on large families (the DL+TL
// self-loop mass is substantial under the DTL model). The branching events
// (S/D/T/SL) each consume a clade split and/or descend the species tree, so
// they stay recursive (bounded by the gene-tree / species-tree size).
template <class REAL>
void UndatedDTLMultiModel<REAL>::btR(CID cid, corax_rnode_t *sp, unsigned int c,
                                     std::vector<unsigned int> &commits,
                                     bool check) {
  // Cap on consecutive "nothing observed" (DL/TL) resamples. At a degenerate
  // cell where the terminating events (None/S/D/T/SL) all have zero mass, the
  // resample loop cannot escape (the stock recursive sampler overflows the
  // stack here -> SIGSEGV; the iterative form would spin). Such a lineage is an
  // unobservable DL/TL chain of negligible probability, so abandon it after the
  // cap (record no further commits). Any cell with escape probability > ~1e-4
  // escapes well within the cap, so the resolution marginals are unbiased.
  size_t resamples = 0;
  const size_t kResampleCap = 1000;
  while (true) {
    if (++_btSteps > kBtStepCap) {
      return; // hard per-sample backstop
    }
    if (++resamples > kResampleCap) {
      ++_btCapHits;
      return;
    }
    btGuard("btR.loop", cid, sp);
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
      btDescend(recCell.event.leftGeneIndex, this->getSpeciesLeft(sp), c,
                commits, check);
      btDescend(recCell.event.rightGeneIndex, this->getSpeciesRight(sp), c,
                commits, check);
      return;
    case ReconciliationEventType::EVENT_D:
      // duplication within the branch (below the WGD): children stay resolved
      btR(recCell.event.leftGeneIndex, sp, c, commits, check);
      btR(recCell.event.rightGeneIndex, sp, c, commits, check);
      return;
    case ReconciliationEventType::EVENT_T:
      // the kept copy (leftGeneIndex) continues here; the transferred copy
      // (rightGeneIndex) lands resolved within the destination branch, below
      // any WGD there -> btR (not btDescend) at the destination.
      btR(recCell.event.leftGeneIndex, sp, c, commits, check);
      btR(recCell.event.rightGeneIndex, recCell.event.pllDestSpeciesNode, c,
          commits, check);
      return;
    case ReconciliationEventType::EVENT_SL: {
      // surviving lineage continues in the non-lost child species. This is a
      // single-lineage tail move (no gene split consumed) that can chain down
      // the species tree, so loop rather than recurse. At a WGD child the R/U
      // coin must fire (may switch to the U state) -> delegate to btDescend.
      corax_rnode_t *child = recCell.event.pllDestSpeciesNode;
      if (!child) {
        std::cerr << "[BTDBG] SL null pllDestSpeciesNode cid=" << cid
                  << " e=" << sp->node_index << std::endl;
        std::cerr.flush();
        std::abort();
      }
      if (_hasWGD[child->node_index]) {
        btDescend(cid, child, c, commits, check);
        return;
      }
      sp = child; // stay resolved, continue down the species tree
      continue;
    }
    case ReconciliationEventType::EVENT_DL:
      continue; // duplication+loss -> nothing observed; resample this cell
    case ReconciliationEventType::EVENT_TL:
      // nothing observed; the single surviving lineage continues. If it died in
      // the sender (pllDestSpeciesNode set) it continues at the destination;
      // otherwise it continues here. Either way, resample (loop, do not recurse).
      if (recCell.event.pllDestSpeciesNode != nullptr) {
        sp = recCell.event.pllDestSpeciesNode;
      }
      continue;
    default:
      return;
    }
  }
}

// One draw from the unresolved (U, vertical-only) recursion at (cid, sp); every
// weight is the exact inside term of updateCLV (normalised U events, synteny
// rho on the WGD split).
template <class REAL>
bool UndatedDTLMultiModel<REAL>::sampleUEventImpl(CID cid, corax_rnode_t *sp,
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
    REAL wLeaf = _dtlclvs[cid]._uq[ec] * (1.0 - r);
    scale(wLeaf);
    total = total + wLeaf;
  } else {
    for (const auto &cs : splits) {
      double w0 = (1.0 - r) * psU * cs.frequency;
      REAL s1 = _dtlclvs[cs.left]._uu[fc] * _dtlclvs[cs.right]._uu[gc] * w0;
      REAL s2 = _dtlclvs[cs.right]._uu[fc] * _dtlclvs[cs.left]._uu[gc] * w0;
      scale(s1);
      scale(s2);
      total = total + s1 + s2;
    }
    double w1 = (1.0 - r) * psU;
    REAL sl1 = _dtlclvs[cid]._uu[fc] * _uEU[gc] * w1;
    REAL sl2 = _dtlclvs[cid]._uu[gc] * _uEU[fc] * w1;
    scale(sl1);
    scale(sl2);
    total = total + sl1 + sl2;
  }
  if (check && !dtlLoreInsideClose(total, _dtlclvs[cid]._uu[ec])) {
    std::cerr << "LORe DTL U-cell weights != _uu at cid=" << cid << " e=" << e
              << std::endl;
    std::abort();
  }
  if (!(total > REAL())) {
    return false;
  }
  REAL toSample = this->getRandom(total);
  REAL acc = wResolve;
  if (wResolve > toSample) {
    REAL dTotal = dupBracket(cid, ec);
    REAL dSample = this->getRandom(dTotal);
    REAL dAcc = REAL();
    for (size_t k = 0; k < splits.size(); ++k) {
      const auto &cs = splits[k];
      REAL t = _dtlclvs[cs.left]._uq[ec] * _dtlclvs[cs.right]._uq[ec] *
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
    ev.type = UEventSample::RESOLVE_SINGLE;
    return true;
  }
  if (isLeaf) {
    ev.type = UEventSample::LEAF;
    return true;
  }
  for (const auto &cs : splits) {
    double w0 = (1.0 - r) * psU * cs.frequency;
    REAL s1 = _dtlclvs[cs.left]._uu[fc] * _dtlclvs[cs.right]._uu[gc] * w0;
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
    REAL s2 = _dtlclvs[cs.right]._uu[fc] * _dtlclvs[cs.left]._uu[gc] * w0;
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
  REAL sl1 = _dtlclvs[cid]._uu[fc] * _uEU[gc] * w1;
  scale(sl1);
  acc = acc + sl1;
  ev.type = UEventSample::SPECIATE_LOSS;
  if (acc > toSample) {
    ev.dest = this->getSpeciesLeft(sp);
    ev.lost = this->getSpeciesRight(sp);
  } else {
    ev.dest = this->getSpeciesRight(sp);
    ev.lost = this->getSpeciesLeft(sp);
  }
  return true;
}

// Unresolved-state backtrace for the resolution profile (vertical-only). The
// single-lineage SL-U tail moves loop rather than recurse.
template <class REAL>
void UndatedDTLMultiModel<REAL>::btU(CID cid, corax_rnode_t *sp, unsigned int c,
                                     std::vector<unsigned int> &commits,
                                     bool check) {
  size_t uIters = 0;
  while (true) {
    if (++_btSteps > kBtStepCap) {
      return; // hard per-sample backstop
    }
    if (++uIters > 5000) {
      ++_btCapHits;
      return;
    }
    btGuard("btU.loop", cid, sp);
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
      _tetra.push_back(e);
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
      sp = ev.dest; // loop, do not recurse
      continue;
    }
  }
}
