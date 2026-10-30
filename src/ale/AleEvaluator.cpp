#include <cstdlib>
#include <functional>
#include "AleOptimizer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cmath>
#include <IO/FileSystem.hpp>

#include <IO/Logger.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <parallelization/ParallelContext.hpp>
#include <search/SpeciesTransferSearch.hpp>

#include "UndatedDLMultiModel.hpp"
#include "UndatedDTLMultiModel.hpp"

static MultiEvaluationPtr
createModel(SpeciesTree &speciesTree, const FamilyInfo &family,
            const RecModelInfo &info, const double alpha,
            const AleModelParameters &modelParameters,
            const std::vector<Highway> &highways, bool highPrecision) {
  MultiEvaluationPtr model;
  GeneSpeciesMapping mapping;
  mapping.fill(family.mappingFile, family.startingGeneTree);
  switch (info.model) {
  case RecModel::UndatedDL:
    if (highPrecision) {
      model = std::make_shared<UndatedDLMultiModel<ScaledValue>>(
          speciesTree.getTree(), mapping, info, family.ccpFile);
    } else {
      model = std::make_shared<UndatedDLMultiModel<double>>(
          speciesTree.getTree(), mapping, info, family.ccpFile);
    }
    break;
  case RecModel::UndatedDTL:
    if (highPrecision) {
      model = std::make_shared<UndatedDTLMultiModel<ScaledValue>>(
          speciesTree.getDatedTree(), mapping, info, family.ccpFile);
    } else {
      model = std::make_shared<UndatedDTLMultiModel<double>>(
          speciesTree.getDatedTree(), mapping, info, family.ccpFile);
    }
    break;
  default:
    assert(false);
  }
  model->setAlpha(alpha);
  RatesVector rates;
  modelParameters.getRateVector(rates);
  model->setRates(rates);
  model->setHighways(highways);
  return model;
}

AleEvaluator::AleEvaluator(
    AleOptimizer &optimizer, SpeciesTree &speciesTree, const RecModelInfo &info,
    const ModelParametrization &modelParametrization,
    const std::string &optimizationClassFile, double &mixtureAlpha,
    std::vector<Highway> &transferHighways,
    std::vector<AleModelParameters> &perLocalFamilyModelParams,
    std::vector<std::string> &wgdBranchLabels,
    std::vector<double> &wgdRetentions, std::vector<double> &wgdResolutions,
    bool &optimizeResolution, bool optimizeRates, bool optimizeVerbose,
    const Families &families, const PerCoreGeneTrees &geneTrees)
    : _optimizer(optimizer), _speciesTree(speciesTree), _info(info),
      _optimizationClasses(_speciesTree.getTree(), modelParametrization,
                           optimizationClassFile, _info),
      _mixtureAlpha(mixtureAlpha), _transferHighways(transferHighways),
      _modelParameters(perLocalFamilyModelParams),
      _optimizeRates(optimizeRates), _optimizeVerbose(optimizeVerbose),
      _families(families), _geneTrees(geneTrees),
      _wgdBranchLabels(wgdBranchLabels), _wgdRetentions(wgdRetentions),
      _wgdResolutionsState(wgdResolutions),
      _optimizeResolution(optimizeResolution),
      _highPrecisions(getLocalFamilyNumber(), -1) {
  Logger::timed << "Initializing ccps and evaluators..." << std::endl;
  _evaluations.resize(getLocalFamilyNumber());
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    resetEvaluation(i, false);
  }
  ParallelContext::barrier();
  // Restore any WGD / LORe state (e.g. from a checkpoint): _wgdBranchLabels /
  // _wgdRetentions / _wgdResolutionsState / _optimizeResolution may already be
  // populated at this point (AleState::unserialize runs before the evaluator
  // is constructed), but the newly-built per-family evaluations above don't
  // know about it yet. This restores q on each evaluation; the per-branch
  // resolution vectors (which need _loreTargets / buildWGDStructure(), only
  // available once declareWGDs() runs post-construction) are pushed there.
  assert(_wgdBranchLabels.size() == _wgdRetentions.size());
  assert(_wgdBranchLabels.size() == _wgdResolutionsState.size());
  if (!_wgdBranchLabels.empty()) {
    auto labelToNode = _speciesTree.getTree().getLabelToNode(false);
    for (unsigned int i = 0; i < _wgdBranchLabels.size(); ++i) {
      auto it = labelToNode.find(_wgdBranchLabels[i]);
      assert(it != labelToNode.end());
      auto speciesNode = it->second->node_index;
      _wgdNodes.push_back(speciesNode);
      _wgdQ.push_back(_wgdRetentions[i]);
      _wgdResolution.push_back(_wgdResolutionsState[i]);
      for (auto &evaluation : _evaluations) {
        evaluation->setWGD(speciesNode, _wgdRetentions[i]);
      }
    }
  }
  unsigned int totalCladesNumber = 0;
  unsigned int worstFamilyCladesNumber = 0;
  for (const auto &evaluation : _evaluations) {
    auto ccpSize = evaluation->getCCP().getCladesNumber();
    totalCladesNumber += ccpSize;
    worstFamilyCladesNumber = std::max(worstFamilyCladesNumber, ccpSize);
  }
  ParallelContext::barrier();
  ParallelContext::sumUInt(totalCladesNumber);
  ParallelContext::maxUInt(worstFamilyCladesNumber);
  double perCoreCladesNumber =
      double(totalCladesNumber) / double(ParallelContext::getSize());
  double loadBalancing = std::min(1.0, double(perCoreCladesNumber) /
                                           double(worstFamilyCladesNumber));
  unsigned int effectiveFamiliesNumber =
      totalCladesNumber / worstFamilyCladesNumber;
  Logger::timed << "Initializing ccps finished" << std::endl;
  Logger::timed << "Total number of clades: " << totalCladesNumber << std::endl;
  Logger::timed << "Load balancing: " << loadBalancing << std::endl;
  Logger::timed << "Recommended maximum number of cores: "
                << effectiveFamiliesNumber << std::endl;
}

void AleEvaluator::resetEvaluation(unsigned int i, bool highPrecision) {
  // testing only: evaluate every family in ScaledValue precision
  static const bool forceHigh = std::getenv("KALERAX_FORCE_HIGH_PRECISION");
  highPrecision = highPrecision || forceHigh;
  auto famIndex = _geneTrees.getTrees()[i].familyIndex;
  _evaluations[i] =
      createModel(_speciesTree, _families[famIndex], _info, _mixtureAlpha,
                  _modelParameters[i], _transferHighways, highPrecision);
  // A rebuilt evaluation (precision switch) must carry the current WGD q,
  // per-branch LORe r and synteny state; createModel() knows none of them.
  reapplyFamilyState(i);
  auto ll = _evaluations[i]->computeLogLikelihood();
  if (highPrecision) {
    _highPrecisions[i] = 1;
  } else {
    _highPrecisions[i] = -1;
    if (!std::isnormal(ll)) {
      resetEvaluation(i, true);
    }
  }
}

void AleEvaluator::setDupSplitWeight(double w) {
  _dupSplitWeight = w;
  for (auto &evaluation : _evaluations) {
    evaluation->setDupSplitWeight(w);
  }
}

void AleEvaluator::reapplyFamilyState(unsigned int i) {
  auto &ev = _evaluations[i];
  if (_dupSplitWeight != 1.0) {
    ev->setDupSplitWeight(_dupSplitWeight);
  }
  for (unsigned int j = 0; j < _wgdNodes.size() && j < _wgdQ.size(); ++j) {
    ev->setWGD(_wgdNodes[j], _wgdQ[j]);
  }
  if (!_perBranchR.empty()) {
    ev->setResolutionProbVector(_perBranchR);
  }
  if (i < _synInputs.size() && _synInputs[i]) {
    ev->setSynteny(*_synInputs[i]);
    ev->setSyntenyParams(_sigmaW, _sigmaN);
  }
}


// Robust 1-D maximisation on [lo, hi]: coarse grid (log-spaced when lo > 0
// and the range spans > 1 decade), then golden-section between the best
// grid point's neighbours. Leaves the model at the best value ever evaluated
// (including `start`), so it never returns something worse than it started.
static double robustMax1D(const std::function<double(double)> &eval, double lo,
                          double hi, double start, int gridN, int iters) {
  double bestX = start, bestF = eval(start);
  std::vector<double> xs(gridN);
  const bool logScale = (lo > 0.0) && (hi / lo > 10.0);
  for (int i = 0; i < gridN; ++i) {
    double t = (gridN == 1) ? 0.5 : double(i) / double(gridN - 1);
    xs[i] = logScale ? std::exp(std::log(lo) + t * (std::log(hi) - std::log(lo)))
                     : lo + t * (hi - lo);
  }
  int bi = -1;
  for (int i = 0; i < gridN; ++i) {
    double f = eval(xs[i]);
    if (f > bestF) {
      bestF = f;
      bestX = xs[i];
      bi = i;
    }
  }
  double a = (bi > 0) ? xs[bi - 1] : (bi == 0 ? lo : std::max(lo, bestX * 0.5));
  double b = (bi >= 0 && bi < gridN - 1) ? xs[bi + 1]
                                          : (bi == gridN - 1 ? hi : std::min(hi, bestX * 1.5 + 1e-9));
  const double gr = 0.6180339887498949;
  double c = b - gr * (b - a), d = a + gr * (b - a);
  double fc = eval(c), fd = eval(d);
  for (int it = 0; it < iters; ++it) {
    if (fc > fd) {
      b = d; d = c; fd = fc; c = b - gr * (b - a); fc = eval(c);
    } else {
      a = c; c = d; fc = fd; d = a + gr * (b - a); fd = eval(d);
    }
  }
  if (fc > bestF) { bestF = fc; bestX = c; }
  if (fd > bestF) { bestF = fd; bestX = d; }
  eval(bestX);
  return bestX;
}

void AleEvaluator::loadSynteny(const std::string &dir, bool unlistedNA,
                               double sigmaW0, double sigmaN0, bool fix,
                               std::vector<double> &state) {
  _synState = &state;
  _fixSynteny = fix;
  if (state.size() == 2) { // resumed: fitted values from the checkpoint
    _sigmaW = state[0];
    _sigmaN = state[1];
  } else {
    _sigmaW = sigmaW0;
    _sigmaN = sigmaN0;
  }
  if (!(_sigmaW > 0.0 && _sigmaW < 1.0 && _sigmaN > 0.0 && _sigmaN < 1.0)) {
    Logger::error << "Error: synteny parameters must lie in (0,1): sigmaW="
                  << _sigmaW << " sigmaN=" << _sigmaN << std::endl;
    ParallelContext::abort(1);
  }
  if (_sigmaWMin > 0.0) {
    if (!(_sigmaWMin < 1.0)) {
      Logger::error << "Error: --synteny-min-sigmaW must lie in [0,1)" << std::endl;
      ParallelContext::abort(1);
    }
    if (_fixSynteny && _sigmaW < _sigmaWMin) {
      Logger::error << "Error: --fix-synteny sigmaW=" << _sigmaW
                    << " is below --synteny-min-sigmaW " << _sigmaWMin << std::endl;
      ParallelContext::abort(1);
    }
    if (_sigmaW < _sigmaWMin) {
      _sigmaW = _sigmaWMin; // start inside the constrained region
    }
    Logger::timed << "[Synteny] sigma_W constrained to [" << _sigmaWMin
                  << ", 1) (--synteny-min-sigmaW)" << std::endl;
  }
  _synInputs.assign(getLocalFamilyNumber(), nullptr);
  unsigned int filesRead = 0, badLines = 0;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto path =
        FileSystem::joinPaths(dir, _geneTrees.getTrees()[i].name + ".syn");
    std::ifstream is(path);
    if (!is.good()) {
      continue;
    }
    auto in = std::make_shared<SyntenyInput>();
    in->unlistedNA = unlistedNA;
    std::string line;
    while (std::getline(is, line)) {
      if (line.empty() || line[0] == '#') {
        continue;
      }
      std::istringstream iss(line);
      std::string a, b, c;
      iss >> a >> b >> c;
      if (a.empty() || b.empty()) {
        continue;
      }
      if (c.empty()) { // "gene NA": all pairs of this gene undetermined
        if (b == "NA" || b == "na") {
          in->naGenes.insert(a);
        } else {
          badLines++;
        }
        continue;
      }
      auto key = syntenyKey(a, b);
      if (c == "1") {
        in->pairObs[key] = 1;
      } else if (c == "0") {
        in->pairObs[key] = 0;
      } else if (c == "NA" || c == "na") {
        in->naPairs.insert(key);
      } else {
        badLines++;
      }
    }
    _synInputs[i] = in;
    filesRead++;
    _evaluations[i]->setSynteny(*in);
    _evaluations[i]->setSyntenyParams(_sigmaW, _sigmaN);
  }
  ParallelContext::sumUInt(filesRead);
  ParallelContext::sumUInt(badLines);
  _synLoaded = (filesRead > 0);
  unsigned int n1 = 0, n0 = 0, fams = 0;
  getSyntenyTotals(n1, n0, fams);
  Logger::timed << "[Synteny] read " << filesRead << " .syn files from " << dir
                << "; families with observed pairs: " << fams
                << "; collinear pairs (y=1): " << n1
                << "; non-collinear (y=0): " << n0
                << "; sigmaW=" << _sigmaW << " sigmaN=" << _sigmaN
                << (fix ? " (fixed)" : "") << std::endl;
  if (badLines) {
    Logger::info << "[Synteny] warning: " << badLines
                 << " unparsable line(s) ignored" << std::endl;
  }
  if (_synState) {
    *_synState = {_sigmaW, _sigmaN};
  }
}

void AleEvaluator::setSyntenyParams(double sigmaW, double sigmaN) {
  _sigmaW = sigmaW;
  _sigmaN = sigmaN;
  for (unsigned int i = 0; i < _evaluations.size(); ++i) {
    if (i < _synInputs.size() && _synInputs[i]) {
      _evaluations[i]->setSyntenyParams(sigmaW, sigmaN);
    }
  }
  if (_synState) {
    *_synState = {_sigmaW, _sigmaN};
  }
}

void AleEvaluator::getSyntenyTotals(unsigned int &n1, unsigned int &n0,
                                    unsigned int &families) {
  n1 = n0 = families = 0;
  for (unsigned int i = 0; i < _evaluations.size(); ++i) {
    if (i < _synInputs.size() && _synInputs[i] &&
        _evaluations[i]->hasSynteny()) {
      unsigned int a = 0, b = 0;
      _evaluations[i]->getSyntenyTotals(a, b);
      n1 += a;
      n0 += b;
      families++;
    }
  }
  ParallelContext::sumUInt(n1);
  ParallelContext::sumUInt(n0);
  ParallelContext::sumUInt(families);
}

double AleEvaluator::optimizeSynteny(bool thorough) {
  const double lo = 1e-4, hi = 1.0 - 1e-4;
  const int grid = thorough ? 12 : 8, iters = thorough ? 24 : 14;
  // sigma_N: every family with observations (null and WGD models alike)
  robustMax1D(
      [&](double x) { setSyntenyParams(_sigmaW, x); return computeLikelihood(); },
      lo, hi, _sigmaN, grid, iters);
  // sigma_W: only identifiable when a WGD is declared
  if (!_wgdNodes.empty()) {
    robustMax1D(
        [&](double x) { setSyntenyParams(x, _sigmaN); return computeLikelihood(); },
        std::max(lo, _sigmaWMin), hi, std::max(_sigmaW, _sigmaWMin), grid, iters);
    // q and sigma_W trade off: re-refine each retention once
    // (unless --fix-wgd-q holds the retentions at their declared q0)
    for (unsigned int j = 0; j < _wgdNodes.size() && !_fixWgdQ; ++j) {
      auto node = _wgdNodes[j];
      robustMax1D(
          [&](double q) { setWGD(node, q); return computeLikelihood(); }, 0.0,
          0.999, _wgdQ[j], thorough ? 10 : 6, iters);
    }
  }
  auto ll = computeLikelihood();
  Logger::timed << "[Synteny]   After synteny opt, ll=" << ll
                << " sigmaW=" << _sigmaW << " sigmaN=" << _sigmaN;
  for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
    Logger::info << " q[node " << _wgdNodes[j] << "]=" << _wgdQ[j];
  }
  Logger::info << std::endl;
  return ll;
}

void AleEvaluator::resetAllPrecisions() {
  auto llBefore = computeLikelihoodFast();
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    if (_highPrecisions[i] != -1) {
      resetEvaluation(i, false);
    }
  }
  auto llAfter = computeLikelihoodFast();
  if (std::fabs(llBefore - llAfter) > 0.1) {
    Logger::info << "Likelihood changed after lowering the precision: "
                 << std::endl;
    Logger::info << "Before: ll=" << llBefore << std::endl;
    Logger::info << "After:  ll=" << llAfter << std::endl;
  }
}

void AleEvaluator::printHighPrecisionCount() {
  unsigned int high = 0;
  unsigned int low = 0;
  for (auto v : _highPrecisions) {
    if (v >= 0) {
      high++;
    } else {
      low++;
    }
  }
  ParallelContext::barrier();
  ParallelContext::sumUInt(high);
  ParallelContext::sumUInt(low);
  Logger::info << "Families in low-precision mode: " << low
               << ", in high-precision mode: " << high << std::endl;
}

double AleEvaluator::computeLikelihoodFast() {
  // printHighPrecisionCount();
  return computeLikelihood();
}

double AleEvaluator::computeLikelihood(PerFamLL *perFamLL) {
  if (perFamLL) {
    perFamLL->clear();
  }
  double sumLL = 0.0;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto ll = computeFamilyLikelihood(i);
    if (perFamLL) {
      perFamLL->push_back(ll);
    }
    sumLL += ll;
  }
  ParallelContext::barrier();
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}

double AleEvaluator::computeFamilyLikelihood(unsigned int i) {
  auto ll = _evaluations[i]->computeLogLikelihood();
  if (_highPrecisions[i] == -1 && !std::isnormal(ll)) {
    // We are in the low precision mode (we use double)
    // and it's not accurate enough, switch to the high
    // precision mode and recompute the ll
    resetEvaluation(i, true);
    ll = _evaluations[i]->computeLogLikelihood();
  }
  if (!std::isnormal(ll)) {
    // Bad ll even in the high precision mode!
    Logger::error << "Error: ll=" << ll << " for family "
                  << _geneTrees.getTrees()[i].name << std::endl;
    assert(false);
  }
  /*
  if (_highPrecisions[i] >= 0 && _highPrecisions[i] % 20 == 0) {
    // We are in the high precision mode, we now check if we can
    // switch to the low precision mode to make computations faster
    resetEvaluation(i, false);
  }
  */
  if (_highPrecisions[i] >= 0) {
    _highPrecisions[i]++;
  }
  return ll;
}

void AleEvaluator::setAlpha(double alpha) {
  _mixtureAlpha = alpha;
  for (auto &evaluation : _evaluations) {
    evaluation->setAlpha(_mixtureAlpha);
  }
}

void AleEvaluator::setFamilyParameters(unsigned int family,
                                       const Parameters &parameters) {
  RatesVector rateVector;
  _modelParameters[family].setParameters(parameters);
  _modelParameters[family].getRateVector(rateVector);
  _evaluations[family]->setRates(rateVector);
}

void AleEvaluator::setWGD(unsigned int speciesNode, double q) {
  auto it = std::find(_wgdNodes.begin(), _wgdNodes.end(), speciesNode);
  if (it == _wgdNodes.end()) {
    _wgdNodes.push_back(speciesNode);
    _wgdQ.push_back(q);
    _wgdResolution.push_back(1.0); // keep the parallel vectors aligned
  } else {
    _wgdQ[std::distance(_wgdNodes.begin(), it)] = q;
  }
  for (auto &evaluation : _evaluations) {
    evaluation->setWGD(speciesNode, q);
  }
  // keep the checkpointed (label-keyed) mirrors in sync
  auto label = _speciesTree.getTree().getNode(speciesNode)->label;
  auto labelIt =
      std::find(_wgdBranchLabels.begin(), _wgdBranchLabels.end(), label);
  if (labelIt == _wgdBranchLabels.end()) {
    _wgdBranchLabels.push_back(label);
    _wgdRetentions.push_back(q);
    _wgdResolutionsState.push_back(1.0); // mirrors _wgdResolution above
  } else {
    _wgdRetentions[std::distance(_wgdBranchLabels.begin(), labelIt)] = q;
  }
}

void AleEvaluator::setResolutionProb(double r) {
  _resolutionProb = r;
  // keep the per-event r state consistent on the AORe baseline / revert paths
  std::fill(_wgdResolution.begin(), _wgdResolution.end(), r);
  std::fill(_wgdResolutionsState.begin(), _wgdResolutionsState.end(), r);
  _perBranchR.assign(_speciesTree.getTree().getNodeNumber(), r);
  for (auto &evaluation : _evaluations) {
    evaluation->setResolutionProb(r);
  }
}

void AleEvaluator::buildWGDStructure() {
  auto &tree = _speciesTree.getTree(); // non-const: isAncestorOf is non-const
  const unsigned int n = _wgdNodes.size();
  _wgdSubtreeBranches.assign(n, {});
  _wgdResolvable.assign(n, 0);
  // When --lore-wgd named specific targets, only those WGDs are resolvable;
  // otherwise (bare --lore) every internal-branch WGD is.
  const bool restricted = !_loreTargets.empty();
  for (unsigned int j = 0; j < n; ++j) {
    const unsigned int w = _wgdNodes[j];
    for (auto node : tree.getNodes()) {
      const unsigned int e = node->node_index;
      if (e == w || tree.isAncestorOf(w, e)) {
        _wgdSubtreeBranches[j].push_back(e);
      }
    }
    // internal WGD (subtree has descendants) => r is a candidate free parameter;
    // a terminal WGD governs only its own branch => r unidentifiable, pinned.
    bool internal = (_wgdSubtreeBranches[j].size() > 1);
    bool selected =
        !restricted || (std::find(_loreTargets.begin(), _loreTargets.end(), w) !=
                        _loreTargets.end());
    _wgdResolvable[j] = (internal && selected) ? 1 : 0;
  }
  // Per-event r is only identifiable when no species branch lies below two
  // *free-r* WGDs (the U-mass recursion carries a single r per branch). A pinned
  // (AORe) WGD resolves at its own branch and sends no U-mass downward, so it
  // may be nested with a free-r WGD; only the free-r set must be disjoint.
  if (_optimizeResolution) {
    std::vector<int> owner(tree.getNodeNumber(), -1); // branch -> WGD idx, -1=none
    for (unsigned int j = 0; j < n; ++j) {
      if (!_wgdResolvable[j]) {
        continue;
      }
      for (auto e : _wgdSubtreeBranches[j]) {
        if (owner[e] >= 0) {
          Logger::error
              << "Error: the LORe-tested WGDs must have DISJOINT subtrees. "
              << "Per-event resolution r is only identifiable when no species "
                 "branch lies below two fitted WGDs, but the WGDs at nodes "
              << _wgdNodes[owner[e]] << " and " << _wgdNodes[j]
              << " are nested (branch " << e << " is below both). "
              << "Use --lore-wgd to fit r for a single WGD (which may be nested "
                 "inside other, AORe-pinned WGDs), or re-run with disjoint "
                 "LORe-tested WGDs." << std::endl;
          ParallelContext::abort(1);
        }
        owner[e] = static_cast<int>(j);
      }
    }
  }
}

void AleEvaluator::setWGDResolutions(const std::vector<double> &rPerWGD) {
  _wgdResolution = rPerWGD;
  // keep the checkpointed (label-keyed) mirror in sync: same size/order as
  // _wgdResolution, since _wgdBranchLabels grows in lockstep with _wgdNodes
  assert(_wgdResolutionsState.size() == rPerWGD.size());
  _wgdResolutionsState = rPerWGD;
  // legacy scalar: first resolvable WGD's r (or 1.0)
  _resolutionProb = 1.0;
  for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
    if (j < _wgdResolvable.size() && _wgdResolvable[j]) {
      _resolutionProb = rPerWGD[j];
      break;
    }
  }
  // paint a per-branch vector: 1.0 everywhere, then each WGD's subtree with its
  // r. On a branch shared by nested WGDs the deepest (most derived) WGD must
  // win: an enclosing AORe WGD resolves at its own branch and sends no U-mass
  // down, so the branch effectively carries only the inner WGD's resolution. We
  // therefore paint ancestors first (largest subtree first); descendants, whose
  // subtree is strictly smaller, overwrite them. Disjoint WGDs never overlap, so
  // this is order-independent for them (bit-for-bit identical to before).
  std::vector<unsigned int> order(_wgdNodes.size());
  for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
    order[j] = j;
  }
  std::sort(order.begin(), order.end(), [&](unsigned int a, unsigned int b) {
    return _wgdSubtreeBranches[a].size() > _wgdSubtreeBranches[b].size();
  });
  std::vector<double> perBranch(_speciesTree.getTree().getNodeNumber(), 1.0);
  for (auto j : order) {
    for (auto e : _wgdSubtreeBranches[j]) {
      perBranch[e] = rPerWGD[j];
    }
  }
  _perBranchR = perBranch;
  for (auto &evaluation : _evaluations) {
    evaluation->setResolutionProbVector(perBranch);
  }
}

/**
 *  Optimizes a set of DTL parameters that are shared among gene families
 */
class DTLGlobalParametersOptimizer : public FunctionToOptimize {
public:
  DTLGlobalParametersOptimizer(AleEvaluator &evaluator)
      : _evaluator(evaluator) {}

  void setParameters(Parameters &parameters) {
    parameters.ensurePositivity();
    auto fullParameters =
        _evaluator.getOptimizationClasses().getFullParameters(parameters);
    for (unsigned int i = 0; i < _evaluator.getLocalFamilyNumber(); ++i) {
      _evaluator.setFamilyParameters(i, fullParameters);
    }
  }

  virtual double evaluate(Parameters &parameters) {
    setParameters(parameters);
    auto res = _evaluator.computeLikelihood();
    parameters.setScore(res);
    return res;
  }

private:
  AleEvaluator &_evaluator;
};

/**
 *  Optimizes DTL parameters for a given family
 */
class DTLFamilyParametersOptimizer : public FunctionToOptimize {
public:
  DTLFamilyParametersOptimizer(AleEvaluator &evaluator, unsigned int family)
      : _evaluator(evaluator), _family(family) {}

  void setParameters(Parameters &parameters) {
    parameters.ensurePositivity();
    auto fullParameters =
        _evaluator.getOptimizationClasses().getFullParameters(parameters);
    _evaluator.setFamilyParameters(_family, fullParameters);
  }

  virtual double evaluate(Parameters &parameters) {
    setParameters(parameters);
    auto res = _evaluator.computeFamilyLikelihood(_family);
    parameters.setScore(res);
    return res;
  }

private:
  AleEvaluator &_evaluator;
  unsigned int _family;
};

/**
 *  Fast-path optimizer for WGD retention probabilities only (one free value
 *  per declared WGD branch), reusing DTLOptimizer::optimizeParameters.
 *
 *  The L-BFGS-B wrapper only guarantees positivity of the free parameters
 *  (r >= 0), but a retention probability must live in [0, 1]. We therefore
 *  optimize an unconstrained r >= 0 and map it through q = r / (1 + r), which
 *  is a smooth bijection from [0, inf) to [0, 1). r = 0 gives q = 0, i.e. the
 *  no-WGD null hypothesis.
 */
class WGDRetentionOptimizer : public FunctionToOptimize {
public:
  WGDRetentionOptimizer(AleEvaluator &evaluator,
                        const std::vector<unsigned int> &wgdNodes,
                        bool optimizeResolution,
                        const std::vector<char> &rFree) // aligned with wgdNodes
      : _evaluator(evaluator), _nodes(wgdNodes),
        _optimizeResolution(optimizeResolution), _rFree(rFree) {}

  // Parameter layout: [ s_q[0..n-1] ]  (+ [ s_r ] for each j with _rFree[j], in
  // node order). Each free s >= 0 maps to a probability via p = s / (1 + s), in
  // [0, 1). Terminal WGDs are excluded from the r block (their r is pinned to
  // AORe).
  void setParameters(Parameters &parameters) {
    parameters.ensurePositivity(); // optimizer guarantees s >= 0
    const unsigned int n = _nodes.size();
    for (unsigned int j = 0; j < n; ++j) {
      double s = parameters[j];
      double q = s / (1.0 + s); // map s in [0, inf) -> q in [0, 1)
      _evaluator.setWGD(_nodes[j], q);
    }
    if (_optimizeResolution) {
      std::vector<double> rVec(n, 1.0);
      unsigned int k = 0;
      for (unsigned int j = 0; j < n; ++j) {
        if (j < _rFree.size() && _rFree[j]) {
          double s = parameters[n + k];
          rVec[j] = s / (1.0 + s); // map s in [0, inf) -> r in [0, 1)
          ++k;
        }
      }
      _evaluator.setWGDResolutions(rVec);
    }
  }

  virtual double evaluate(Parameters &parameters) {
    setParameters(parameters);
    auto res = _evaluator.computeLikelihood();
    parameters.setScore(res);
    return res;
  }

private:
  AleEvaluator &_evaluator;
  std::vector<unsigned int> _nodes;
  bool _optimizeResolution;
  std::vector<char> _rFree;
};

double AleEvaluator::optimizeModelRates(bool thorough) {
  auto ll = computeLikelihood();
  if (_optimizeRates) {
    Logger::timed << "[Species search] Optimizing model rates ";
    OptimizationSettings settings;
    settings.listeners.push_back(&_optimizer);
    settings.verbose = _optimizeVerbose;
    settings.strategy = _info.recOpt;
    settings.lineSearchMinImprovement = std::max(0.1, -ll / 10000.0);
    if (!thorough) {
      Logger::info << "(light), ll=" << ll << std::endl;
      settings.minAlpha = 0.01;
      settings.startingAlpha = 0.5;
    } else {
      Logger::info << "(thorough), ll=" << ll << std::endl;
      if (-ll < 100.0) {
        settings.lineSearchMinImprovement = 0.01;
      }
      settings.minAlpha = 0.005;
      settings.startingAlpha = 0.01;
    }
    settings.optimizationMinImprovement = settings.lineSearchMinImprovement;
    settings.factr = LBFGSBPrecision::MEDIUM;
    if (_info.perFamilyRates) {
      Logger::timed << "[Species search]   Free parameters: "
                    << _optimizationClasses.getFreeParameters() << " per family"
                    << std::endl;
      _optimizer.enableCheckpoints(false); // to avoid MPI issues
      for (unsigned int family = 0; family < getLocalFamilyNumber(); ++family) {
        DTLFamilyParametersOptimizer function(*this, family);
        auto categorizedParameters =
            getOptimizationClasses().getCompressedParameters(
                _modelParameters[family].getParameters());
        auto bestParameters = DTLOptimizer::optimizeParameters(
            function, categorizedParameters, settings);
        // set the found best parameters to the family
        function.setParameters(bestParameters);
      }
      _optimizer.enableCheckpoints(true); // to avoid MPI issues
      ll = computeLikelihood();
    } else {
      Logger::timed << "[Species search]   Free parameters: "
                    << _optimizationClasses.getFreeParameters() << std::endl;
      DTLGlobalParametersOptimizer function(*this);
      auto categorizedParameters =
          getOptimizationClasses().getCompressedParameters(
              _modelParameters[0].getParameters());
      ParallelContext::barrier();
      auto bestParameters = DTLOptimizer::optimizeParameters(
          function, categorizedParameters, settings);
      // set the found best parameters to all families
      function.setParameters(bestParameters);
      ll = computeLikelihood();
    }
    Logger::timed << "[Species search]   After model rate opt, ll=" << ll
                  << std::endl;
    //            << ", rates=" << _modelParameters << std::endl;
  }
  // Fast-path WGD retention optimization: a separate pass after the DTL-rate
  // optimization, one free retention per declared WGD branch. Guarded only by
  // whether any WGD is declared (independent of _optimizeRates, so q can be
  // estimated even with D/L fixed, as the WHALE cross-validation needs).
  if (!_wgdNodes.empty()) {
    // Build the per-WGD subtree sets and resolvable mask first, so the masks /
    // subtrees are ready and the disjointness check (under --lore) fires before
    // any fitting.
    buildWGDStructure();
    unsigned int nFreeR = 0;
    for (auto c : _wgdResolvable) {
      nFreeR += (c ? 1u : 0u);
    }
    Logger::timed << "[Species search] Optimizing WGD retention(s) on "
                  << _wgdNodes.size() << " branch(es)";
    if (_optimizeResolution) {
      Logger::info << " + per-event LORe r (" << nFreeR << " free, "
                   << (_wgdNodes.size() - nFreeR)
                   << " pinned:terminal/non-target)";
    }
    Logger::info << ", ll=" << ll << std::endl;
    OptimizationSettings settings;
    settings.listeners.push_back(&_optimizer);
    settings.verbose = _optimizeVerbose;
    settings.strategy = _info.recOpt;
    settings.lineSearchMinImprovement = std::max(0.1, -ll / 10000.0);
    settings.optimizationMinImprovement = settings.lineSearchMinImprovement;
    settings.factr = LBFGSBPrecision::MEDIUM;
    if (thorough) {
      settings.minAlpha = 0.005;
      settings.startingAlpha = 0.01;
    } else {
      settings.minAlpha = 0.01;
      settings.startingAlpha = 0.5;
    }
    // Starting s-vector for the WGD retentions from the current per-branch q
    // (q = s/(1+s), so s = q/(1-q)).
    auto startQ = [&]() {
      std::vector<double> s(_wgdNodes.size());
      for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
        double q = std::min(std::max(_wgdQ[j], 0.0), 0.999999);
        s[j] = q / (1.0 - q);
      }
      return s;
    };
    auto logQ = [&]() {
      for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
        Logger::info << " q[node " << _wgdNodes[j] << "]=" << _wgdQ[j];
      }
    };
    auto logR = [&]() {
      for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
        Logger::info << " r[node " << _wgdNodes[j] << "]=" << _wgdResolution[j];
        if (j >= _wgdResolvable.size() || !_wgdResolvable[j]) {
          // terminal-branch WGDs are unidentifiable; internal WGDs that are not
          // a --lore-wgd target are pinned to AORe by choice.
          bool terminal =
              (j < _wgdSubtreeBranches.size() && _wgdSubtreeBranches[j].size() <= 1);
          Logger::info << (terminal ? "(pinned;terminal)" : "(pinned;non-target)");
        }
      }
    };
    // Robust 1-D refinement (golden-section) of a single retention/resolution
    // parameter in [lo, hi]. The retention likelihood is unimodal in q (and in
    // r) but very flat near the q->1 boundary, where the gradient method in the
    // s=q/(1+s) chart under-climbs and stalls (e.g. at q=2/3, s=2). Golden-
    // section finds the true maximum regardless of the gradient scaling.
    // eval(x) sets the parameter on the model and returns the joint lnL.
    auto goldenRefine = [](const std::function<double(double)> &eval, double lo,
                           double hi, int iters) -> double {
      const double gr = 0.6180339887498949;
      // Track the best point ever evaluated, including both endpoints: the
      // optimum is often ON the boundary (q = 0: no WGD; r = 0), which the
      // interior golden-section points never reach exactly.
      double bestX = lo, bestF = eval(lo);
      {
        double fh = eval(hi);
        if (fh > bestF) { bestF = fh; bestX = hi; }
      }
      double a = lo, b = hi;
      double c = b - gr * (b - a), d = a + gr * (b - a);
      double fc = eval(c), fd = eval(d);
      for (int it = 0; it < iters; ++it) {
        if (fc > fd) {
          b = d; d = c; fd = fc; c = b - gr * (b - a); fc = eval(c);
        } else {
          a = c; c = d; fc = fd; d = a + gr * (b - a); fd = eval(d);
        }
      }
      if (fc > bestF) { bestF = fc; bestX = c; }
      if (fd > bestF) { bestF = fd; bestX = d; }
      eval(bestX); // leave the model at the best value
      return bestX;
    };

    const bool rFixed = _optimizeResolution && _fixedR >= 0.0;
    if (rFixed || _fixWgdQ) {
      // Constrained WGD pass (--lore-fixed-r and/or --fix-wgd-q): no AORe/LORe
      // model choice, no revert. r is set to the fixed value (resolvable WGDs;
      // others pinned to 1), or kept free and golden-refined with q held; q is
      // held at its declared value, or fitted with r held.
      const unsigned int n = _wgdNodes.size();
      std::vector<double> rVec(n, 1.0);
      if (_optimizeResolution) {
        for (unsigned int j = 0; j < n; ++j) {
          if (_wgdResolvable[j]) {
            double cur = (j < _wgdResolution.size()) ? _wgdResolution[j] : 0.9;
            rVec[j] = rFixed ? _fixedR : std::min(cur, 0.9);
          }
        }
        setWGDResolutions(rVec);
      } else {
        setResolutionProb(1.0);
      }
      if (!_fixWgdQ) {
        // q only (the optimizer leaves the per-event r state untouched)
        Parameters startQOnly(startQ());
        WGDRetentionOptimizer qFun(*this, _wgdNodes, false, {});
        ParallelContext::barrier();
        auto bestQ = DTLOptimizer::optimizeParameters(qFun, startQOnly, settings);
        qFun.setParameters(bestQ);
        if (thorough) {
          for (unsigned int j = 0; j < n; ++j) {
            auto node = _wgdNodes[j];
            goldenRefine(
                [&](double q) { setWGD(node, q); return computeLikelihood(); },
                0.0, 0.999, 16);
          }
        }
      }
      if (_optimizeResolution && !rFixed) {
        // r only, q held (golden-section on each resolvable WGD's r)
        for (unsigned int j = 0; j < n; ++j) {
          if (!_wgdResolvable[j]) {
            continue;
          }
          goldenRefine(
              [&](double r) {
                rVec[j] = r;
                setWGDResolutions(rVec);
                return computeLikelihood();
              },
              0.0, 0.999, thorough ? 24 : 16);
        }
      }
      ll = computeLikelihood();
      Logger::timed << "[Species search]   After constrained WGD pass ("
                    << (_fixWgdQ ? "q fixed" : "q fitted") << ", "
                    << (!_optimizeResolution ? "AORe"
                                             : (rFixed ? "r fixed" : "r fitted"))
                    << "), ll=" << ll;
      logQ();
      logR();
      Logger::info << std::endl;
    } else {
    // (1) AORe optimum: WGD retentions only, with r = 1 (the nested null).
    setResolutionProb(1.0);
    Parameters startA(startQ());
    WGDRetentionOptimizer aoreFun(*this, _wgdNodes, false, {});
    ParallelContext::barrier();
    auto bestA = DTLOptimizer::optimizeParameters(aoreFun, startA, settings);
    aoreFun.setParameters(bestA); // applies the gradient AORe q-hat (-> _wgdQ)
    // Golden-section refinement of each retention (r=1): the gradient stalls on
    // the flat high-q tail (in s-space), so refine in q-space directly. Done in
    // the thorough pass only (the light pass just warm-starts it).
    if (thorough) {
      for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
        auto node = _wgdNodes[j];
        goldenRefine(
            [&](double q) { setWGD(node, q); return computeLikelihood(); }, 0.0,
            0.999, 16);
      }
    }
    double llAore = computeLikelihood();
    std::vector<double> aoreQ = _wgdQ; // refined AORe retentions (for revert)

    if (!_optimizeResolution) {
      ll = llAore;
      Logger::timed << "[Species search]   After WGD retention opt, ll=" << ll;
      logQ();
      Logger::info << std::endl;
    } else {
      // (2) LORe optimum: retentions + a *per-event* resolution r, one free r
      //     per *resolvable* (internal-branch) WGD, seeded at the refined AORe
      //     q-hat and r = 0.9 (near-AORe). Terminal-branch WGDs keep r pinned to
      //     1 (AORe). r = 1 is unreachable in the s/(1+s) chart, so we keep
      //     whichever of AORe / LORe is better.
      const std::vector<char> &rFree = _wgdResolvable;

      // start vector: q-hat warm start, then r=0.9 for each resolvable WGD
      std::vector<double> startL = startQ();
      for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
        if (rFree[j]) {
          startL.push_back(0.9 / (1.0 - 0.9)); // r = 0.9
        }
      }

      // seed the model: 0.9 on resolvable WGDs, 1.0 (pinned) on terminal ones
      std::vector<double> r0(_wgdNodes.size(), 1.0);
      for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
        if (rFree[j]) {
          r0[j] = 0.9;
        }
      }
      setWGDResolutions(r0);

      Parameters startLore(startL);
      WGDRetentionOptimizer loreFun(*this, _wgdNodes, true, rFree);
      ParallelContext::barrier();
      auto bestL = DTLOptimizer::optimizeParameters(loreFun, startLore, settings);
      loreFun.setParameters(bestL);

      // Coordinate-wise golden-section refinement of the retentions and the
      // per-event r (thorough only). setWGD changes only q (leaves the per-event
      // r state intact); the pinned (!rFree[j]) entries of rVec stay 1.0.
      if (thorough) {
        std::vector<double> rVec = _wgdResolution; // current per-event r (size n)
        for (int sweep = 0; sweep < 2; ++sweep) {
          for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
            auto node = _wgdNodes[j];
            goldenRefine(
                [&](double q) { setWGD(node, q); return computeLikelihood(); },
                0.0, 0.999, 16);
          }
          for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
            if (!rFree[j]) {
              continue;
            }
            goldenRefine(
                [&](double r) {
                  rVec[j] = r;
                  setWGDResolutions(rVec);
                  return computeLikelihood();
                },
                0.0, 0.999, 16);
          }
        }
      }
      double llLore = computeLikelihood();

      if (llLore >= llAore) {
        ll = llLore; // delayed resolution improves the fit -> keep LORe
        Logger::timed << "[Species search]   After WGD+LORe opt, ll=" << ll;
        logQ();
        logR();
        Logger::info << std::endl;
      } else {
        // LORe could not beat AORe -> revert to the refined AORe q-hat, r = 1.
        for (unsigned int j = 0; j < _wgdNodes.size(); ++j) {
          setWGD(_wgdNodes[j], aoreQ[j]);
        }
        setResolutionProb(1.0); // also fills _wgdResolution with 1.0
        ll = computeLikelihood();
        Logger::timed << "[Species search]   After WGD+LORe opt, ll=" << ll
                      << " (LORe did not improve; reverted to AORe r=1)";
        logQ();
        logR();
        Logger::info << std::endl;
      }
    }
    } // end unconstrained WGD pass
  }
  // synteny emission parameters (Layer 1)
  if (_synLoaded && !_fixSynteny) {
    optimizeSynteny(thorough);
  }
  ll = optimizeGammaRates();
  resetAllPrecisions();
  ll = computeLikelihood();
  return ll;
}

static double callback(void *p, double x) {
  auto *evaluator = (AleEvaluator *)p;
  evaluator->setAlpha(x);
  auto ll = evaluator->computeLikelihood();
  return -ll;
}

// kalerax: mean per-family event multiplier E_c[1/s_c] of the gamma categories.
// In category c, D/L/T are effectively scaled by 1/s_c (speciation weight s_c),
// so alpha and the overall rate lie on a strong ridge.
static double meanEventMultiplier(double alpha, unsigned int ncat) {
  std::vector<double> sc(ncat);
  corax_compute_gamma_cats(alpha, ncat, &sc[0], CORAX_GAMMA_RATES_MEAN);
  double m = 0.0;
  for (auto v : sc) {
    m += 1.0 / v;
  }
  return m / static_cast<double>(ncat);
}

struct CompensatedAlphaCtx {
  AleEvaluator *evaluator;
  std::vector<Parameters> base;  // per local family, at the starting alpha
  double m0;
  unsigned int ncat;
};

// alpha step that keeps the mean event intensity fixed: rates are rescaled by
// m(alpha0)/m(alpha) so the move runs along the ridge rather than across it
static double compensatedCallback(void *p, double x) {
  auto *ctx = (CompensatedAlphaCtx *)p;
  double f = ctx->m0 / meanEventMultiplier(x, ctx->ncat);
  for (unsigned int i = 0; i < ctx->base.size(); ++i) {
    ctx->evaluator->setFamilyParameters(i, ctx->base[i] * f);
  }
  ctx->evaluator->setAlpha(x);
  return -ctx->evaluator->computeLikelihood();
}

double AleEvaluator::optimizeGammaRates() {
  auto ll = computeLikelihood();
  if (_info.gammaCategories == 1 || _fixGammaAlpha) {
    return ll;
  }
  Logger::timed << "[Species search] Optimizing gamma categories" << std::endl;
  double minAlpha = CORAX_OPT_MIN_ALPHA;
  double maxAlpha = CORAX_OPT_MAX_ALPHA;
  double startingAlpha = _mixtureAlpha;
  // kalerax: 0.1 (upstream) is far too coarse; the lnL can be sharply peaked
  // in alpha (tens of lnL units per 0.1), and Brent then stalls off-optimum
  double tolerance = 1e-3;
  double fx = -ll;
  double f2x = 1.0;
  ParallelContext::barrier();
  double alpha = startingAlpha;
  if (_optimizeRates && !_info.perFamilyRates) {
    // rate-compensated step (kalerax), then a plain alpha step at the new rates
    CompensatedAlphaCtx ctx;
    ctx.evaluator = this;
    ctx.ncat = _info.gammaCategories;
    ctx.m0 = meanEventMultiplier(startingAlpha, ctx.ncat);
    for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
      ctx.base.push_back(_modelParameters[i].getParameters());
    }
    double fc = -ll;
    double f2c = 1.0;
    double ac = corax_opt_minimize_brent(minAlpha, startingAlpha, maxAlpha,
                                         tolerance, &fc, &f2c, (void *)&ctx,
                                         &compensatedCallback);
    if (-fc > ll) {
      compensatedCallback((void *)&ctx, ac);
      alpha = ac;
      fx = fc;
    } else {  // no gain: restore the starting point exactly
      for (unsigned int i = 0; i < ctx.base.size(); ++i) {
        setFamilyParameters(i, ctx.base[i]);
      }
      setAlpha(startingAlpha);
    }
    startingAlpha = alpha;
  }
  alpha = corax_opt_minimize_brent(minAlpha, startingAlpha, maxAlpha, tolerance,
                                   &fx, &f2x, (void *)this, &callback);
  ll = -fx;
  setAlpha(alpha);
  std::vector<double> categories(_info.gammaCategories);
  corax_compute_gamma_cats(alpha, categories.size(), &categories[0],
                           CORAX_GAMMA_RATES_MEAN);
  Logger::timed << "[Species search]   After gamma cat  opt, ll=" << ll
                << " alpha=" << alpha << std::endl;
  if (_optimizeVerbose) {
    Logger::info << "alpha=" << alpha << std::endl;
    Logger::info << "speciation rate categories: ";
    for (auto c : categories) {
      Logger::info << c << " ";
    }
    Logger::info << std::endl;
  }
  return ll;
}

void AleEvaluator::onSpeciesDatesChange() {
  for (auto &evaluation : _evaluations) {
    evaluation->onSpeciesDatesChange();
  }
}

void AleEvaluator::onSpeciesTreeChange(
    const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
  for (auto &evaluation : _evaluations) {
    evaluation->onSpeciesTreeChange(nodesToInvalidate);
  }
}

void AleEvaluator::sampleFamilyResolutionCommits(
    unsigned int i, unsigned int samples, std::vector<double> &commitCounts,
    std::vector<double> &tetraCounts) {
  assert(i < getLocalFamilyNumber());
  // ensure the family's CLVs reflect the current (fitted) rates and r
  computeFamilyLikelihood(i);
  getEvaluation(i).sampleResolutionCommits(samples, commitCounts, tetraCounts);
}

void AleEvaluator::sampleFamilyScenarios(
    unsigned int i, unsigned int samples,
    std::vector<std::shared_ptr<Scenario>> &scenarios) {
  assert(i < getLocalFamilyNumber());
  scenarios.clear();
  bool ok = getEvaluation(i).sampleReconciliations(samples, scenarios);
  if (_highPrecisions[i] == -1 && !ok) {
    // We are in the low precision mode (we use double)
    // and it's not accurate enough, switch to the high
    // precision mode and resample
    scenarios.clear();
    resetEvaluation(i, true);
    ok = getEvaluation(i).sampleReconciliations(samples, scenarios);
  }
  if (!ok) {
    // Couldn't sample even in the high precision mode!
    Logger::error << "Error: cannot sample reconciliations for family "
                  << _geneTrees.getTrees()[i].name << std::endl;
    assert(false);
  }
}

void AleEvaluator::getTransferInformation(
    SpeciesTree &speciesTree, TransferFrequencies &transferFrequencies,
    PerSpeciesEvents &perSpeciesEvents,
    PerCorePotentialTransfers &potentialTransfers) {
  // this is duplicated code from Routines...
  const auto labelToId = speciesTree.getTree().getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getTree().getDeterministicIdToLabel();
  const unsigned int labelsNumber = idToLabel.size();
  transferFrequencies.count =
      MatrixUint(labelsNumber, VectorUint(labelsNumber, 0));
  transferFrequencies.idToLabel = idToLabel;
  perSpeciesEvents = PerSpeciesEvents(speciesTree.getTree().getNodeNumber());
  auto infoCopy = _info;
  infoCopy.model = RecModel::UndatedDTL;
  infoCopy.originationStrategy = OriginationStrategy::UNIFORM;
  infoCopy.transferConstraint = TransferConstaint::PARENTS;
  for (const auto &geneTree : _geneTrees.getTrees()) {
    const auto &family = _families[geneTree.familyIndex];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    UndatedDTLMultiModel<ScaledValue> evaluation(
        speciesTree.getDatedTree(), mapping, infoCopy, family.ccpFile);
    std::vector<std::shared_ptr<Scenario>> scenarios;
    // Warning:
    // Using Random::getProba() in the sampling function makes
    // the random state inconsistent between the MPI ranks.
    // Call ParallelContext::makeRandConsistent() right after
    // all MPI ranks passed the loop
    bool ok = evaluation.sampleReconciliations(1, scenarios);
    assert(ok);
    assert(scenarios.size() == 1);
    auto &scenario = *scenarios[0];
    scenario.countTransfers(labelToId, transferFrequencies.count);
    scenario.gatherReconciliationStatistics(perSpeciesEvents);
    potentialTransfers.addScenario(scenario);
  }
  ParallelContext::barrier();
  ParallelContext::makeRandConsistent();
  for (unsigned int i = 0; i < labelsNumber; ++i) {
    ParallelContext::sumVectorUInt(transferFrequencies.count[i]);
  }
  perSpeciesEvents.parallelSum();
  assert(ParallelContext::isRandConsistent());
}

void AleEvaluator::addHighway(const Highway &highway) {
  _transferHighways.push_back(highway);
  for (auto &evaluation : _evaluations) {
    evaluation->setHighways(_transferHighways);
  }
}

void AleEvaluator::removeHighway() {
  _transferHighways.pop_back();
  for (auto &evaluation : _evaluations) {
    evaluation->setHighways(_transferHighways);
  }
}

void AleEvaluator::saveSnapshotPerFamilyLL() {
  std::vector<double> localLikelihoods;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto ll = computeFamilyLikelihood(i);
    localLikelihoods.push_back(ll);
  }
  ParallelContext::barrier();
  ParallelContext::concatenateHeterogeneousDoubleVectors(localLikelihoods,
                                                         _snapshotPerFamilyLL);
  assert(_snapshotPerFamilyLL.size() == _families.size());
}

void AleEvaluator::savePerFamilyLikelihoodDiff(const std::string &outputFile) {
  std::vector<unsigned int> localIndices;
  std::vector<double> localLikelihoods;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    auto ll = computeFamilyLikelihood(i);
    localIndices.push_back(famIndex);
    localLikelihoods.push_back(ll);
  }
  ParallelContext::barrier();
  std::vector<unsigned int> indices;
  std::vector<double> likelihoods;
  ParallelContext::concatenateHeterogeneousUIntVectors(localIndices, indices);
  ParallelContext::concatenateHeterogeneousDoubleVectors(localLikelihoods,
                                                         likelihoods);
  assert(indices.size() == _snapshotPerFamilyLL.size());
  if (ParallelContext::getRank() == 0) {
    std::vector<ScoredFamily> scoredFamilies;
    for (unsigned int i = 0; i < indices.size(); ++i) {
      const auto &family = _families[indices[i]];
      auto ll = likelihoods[i];
      scoredFamilies.push_back(
          ScoredFamily(family.name, ll - _snapshotPerFamilyLL[i]));
    }
    std::sort(scoredFamilies.begin(), scoredFamilies.end());
    std::ofstream os(outputFile);
    os << "fam, llDiff" << std::endl;
    for (const auto &sf : scoredFamilies) {
      os << sf.familyName << ", " << sf.score << std::endl;
    }
    os.close();
  }
  ParallelContext::barrier();
}

unsigned int AleEvaluator::getInputTreesNumber() const {
  unsigned int totalInputTrees = 0;
  for (const auto &evaluation : _evaluations) {
    totalInputTrees += evaluation->getCCP().getInputTreesNumber();
  }
  ParallelContext::barrier();
  ParallelContext::sumUInt(totalInputTrees);
  return totalInputTrees;
}
