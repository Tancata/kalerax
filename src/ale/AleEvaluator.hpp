#pragma once

#include <algorithm>
#include <iterator>
#include <memory>

#include <IO/Families.hpp>
#include <IO/HighwayCandidateParser.hpp>
#include <search/SpeciesSearchCommon.hpp>
#include <util/types.hpp>

#include "AleModelParameters.hpp"
#include "OptimizationClasses.hpp"

class AleOptimizer;
class SpeciesTree;
class RecModelInfo;
class PerCoreGeneTrees;
class MultiModelInterface;
class Scenario;

using MultiEvaluation = MultiModelInterface;
using MultiEvaluationPtr = std::shared_ptr<MultiEvaluation>;
using PerCoreMultiEvaluations = std::vector<MultiEvaluationPtr>;

struct ScoredFamily {
  std::string familyName;
  double score;
  ScoredFamily() : score(0.0) {}
  ScoredFamily(const std::string &familyName, double score = 0.0)
      : familyName(familyName), score(score) {}
  bool operator<(const ScoredFamily &other) const {
    if (score == other.score) {
      return familyName < other.familyName;
    }
    return score < other.score;
  }
};

/**
 *  Implements the main functions to evaluate the species tree likelihood
 *  and to sample reconciliations
 */
class AleEvaluator : public SpeciesTreeLikelihoodEvaluatorInterface {
public:
  /**
   *  Constructor
   *  @param optimizer
   *  @param speciesTree Reference to the AleState species tree
   *  @param info Description of the reconciliation model
   *  @param modelParametrization Describes how parameters should be optimized
   *  @param optimizationClassFile Path to the file defining model parameter
   *                               optimization classes
   *  @param mixtureAlpha Reference to the AleState model alpha
   *  @param transferHighways Reference to the AleState transfer highways
   *  @param perLocalFamilyModelParams Reference to the AleState model params
   *  @param wgdBranchLabels Reference to the AleState WGD branch labels
   *  @param wgdRetentions Reference to the AleState WGD retentions (aligned
   *                        with wgdBranchLabels)
   *  @param wgdResolutions Reference to the AleState per-event LORe
   *                         resolutions (aligned with wgdBranchLabels)
   *  @param optimizeResolution Reference to the AleState LORe optimize flag
   *  @param optimizeRates If false, model parameter optimization is skipped
   *  @param optimizeVerbose If true, the optimization routines print more logs
   *  @param families List of all gene families
   *  @param geneTrees List of the gene families assigned to the MPI rank
   */
  AleEvaluator(AleOptimizer &optimizer, SpeciesTree &speciesTree,
               const RecModelInfo &info,
               const ModelParametrization &modelParametrization,
               const std::string &optimizationClassFile, double &mixtureAlpha,
               std::vector<Highway> &transferHighways,
               std::vector<AleModelParameters> &perLocalFamilyModelParams,
               std::vector<std::string> &wgdBranchLabels,
               std::vector<double> &wgdRetentions,
               std::vector<double> &wgdResolutions, bool &optimizeResolution,
               bool optimizeRates, bool optimizeVerbose,
               const Families &families, const PerCoreGeneTrees &geneTrees);

  /**
   *  Destructor
   */
  virtual ~AleEvaluator() {}

  /**
   *  See description in the parent class.
   *  We don't need to implement those functions in AleRax
   */
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}

  /**
   *  Compute the species tree likelihood
   *
   *  @param perFamLL Structure to store the per-family likelihoods (if set)
   *  @return The log likelihood
   */
  virtual double computeLikelihood(PerFamLL *perFamLL = nullptr);

  /**
   *  Fast approximated version of computeLikelihood, used by some
   *  of the optimization heuristics to speedup the search, unless
   *  providesFastLikelihoodImpl returns false (and so far it returns false)
   */
  virtual double computeLikelihoodFast();
  virtual bool providesFastLikelihoodImpl() const { return false; }

  /**
   *  Compute the likelihood of family i (using the _geneTrees indexing)
   */
  double computeFamilyLikelihood(unsigned int i);

  /**
   *  Functions to set the recmodel parameters
   */
  void setAlpha(double alpha);
  void setFamilyParameters(unsigned int family, const Parameters &parameters);

  /**
   *  Declare a whole-genome duplication at the top of a species branch (by
   *  node index) with retention probability q, and broadcast it to every local
   *  family evaluation. The branch is registered so the retention optimizer can
   *  later sweep its q. Calling this again for the same branch just updates q.
   */
  void setWGD(unsigned int speciesNode, double q);

  /**
   *  The species branches (node indices) on which a WGD has been declared.
   */
  const std::vector<unsigned int> &getWGDNodes() const { return _wgdNodes; }

  /**
   *  The current (fitted) retention probability q of the WGD on species branch
   *  `speciesNode` (1.0 if no WGD is declared there).
   */
  double getWGDRetention(unsigned int speciesNode) const {
    auto it = std::find(_wgdNodes.begin(), _wgdNodes.end(), speciesNode);
    return (it == _wgdNodes.end())
               ? 1.0
               : _wgdQ[std::distance(_wgdNodes.begin(), it)];
  }

  /**
   *  Legacy scalar accessor: now returns the LORe resolution r of the first
   *  resolvable (internal-branch) WGD, or 1.0 if there is none. Prefer the
   *  per-event getWGDResolution() / getWGDResolvableMask() below.
   */
  double getResolutionProb() const { return _resolutionProb; }

  /**
   *  Per-event LORe resolution r of the WGD on branch `speciesNode` (1.0 if
   *  none, or if pinned to AORe because the WGD sits on a terminal branch).
   */
  double getWGDResolution(unsigned int speciesNode) const {
    auto it = std::find(_wgdNodes.begin(), _wgdNodes.end(), speciesNode);
    return (it == _wgdNodes.end())
               ? 1.0
               : _wgdResolution[std::distance(_wgdNodes.begin(), it)];
  }

  /**
   *  Per declared WGD (aligned with getWGDNodes()): is its r a *free* parameter?
   *  False for terminal-branch WGDs (r pinned to AORe). Valid after
   *  buildWGDStructure().
   */
  const std::vector<char> &getWGDResolvableMask() const { return _wgdResolvable; }

  /**
   *  Set per-event resolution probs (aligned with getWGDNodes()). Each WGD
   *  paints its own subtree; branches below no WGD stay at r=1. Requires
   *  disjoint WGD subtrees (enforced in buildWGDStructure under --lore).
   */
  void setWGDResolutions(const std::vector<double> &rPerWGD);

  /**
   *  Build the per-WGD subtree branch sets and the resolvable mask (terminal
   *  WGDs => r pinned to AORe). Under --lore with more than one WGD, abort if
   *  the WGD subtrees overlap (nested WGDs => per-event r is unidentifiable).
   */
  void buildWGDStructure();

  /**
   *  Set the (global, phase-1) LORe resolution probability r in [0,1] and
   *  broadcast it to every local family evaluation. r=1 == WHALE/AORe.
   */
  void setResolutionProb(double r);

  /**
   *  Enable/disable optimisation of the LORe resolution probability r. When
   *  enabled, r is fitted jointly with the WGD retentions in the fast-path pass.
   */
  void enableResolutionOptimization(bool on) { _optimizeResolution = on; }
  bool optimizesResolution() const { return _optimizeResolution; }

  /**
   *  Restrict LORe to a specific declared WGD (by node index): only the
   *  registered targets get a free resolution r; every other WGD is pinned to
   *  AORe (r=1). When no target is registered, every internal-branch WGD is
   *  resolvable (the default --lore behaviour). Targets must sit on declared
   *  WGD branches and be mutually disjoint, but may be nested inside (pinned)
   *  WGDs. Consumed by buildWGDStructure().
   */
  void addResolutionTarget(unsigned int speciesNode) {
    _loreTargets.push_back(speciesNode);
  }
  const std::vector<unsigned int> &getResolutionTargets() const {
    return _loreTargets;
  }

  /**
   *  Suppress the bulky per-family/per-sample output (--summary-only). Read by
   *  AleOptimizer::reconcile() and the highway candidate tests.
   */
  void setSummaryOnly(bool on) { _summaryOnly = on; }
  bool isSummaryOnly() const { return _summaryOnly; }

  /**
   *  Optimize the model rates
   */
  virtual double optimizeModelRates(bool thorough);

  /**
   *  Update node information after modifying the species tree
   */
  virtual void onSpeciesDatesChange();
  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);

  /**
   *  Sample reconciliations for family i (using the _geneTrees indexing)
   */
  void sampleFamilyScenarios(unsigned int i, unsigned int samples,
                             std::vector<std::shared_ptr<Scenario>> &scenarios);

  /**
   *  Sample `samples` LORe resolution histories for family i and accumulate,
   *  per species branch, the number of U->R commit (WGD resolution) events.
   *  Empty for models / families without a WGD. CLVs are (re)computed first.
   */
  void sampleFamilyResolutionCommits(unsigned int i, unsigned int samples,
                                     std::vector<double> &commitCounts,
                                     std::vector<double> &tetraCounts);

  /**
   *  Sample reconciliations and return the HGT frequencies and
   *  per-species events
   *
   *  Obtain possible transfers to be used in transfer-guided
   *  species tree search, relative species dating and
   *  transfer highway inference
   */
  virtual void
  getTransferInformation(SpeciesTree &speciesTree,
                         TransferFrequencies &frequencies,
                         PerSpeciesEvents &perSpeciesEvents,
                         PerCorePotentialTransfers &potentialTransfers);

  /**
   *  Functions to handle highways
   */
  void addHighway(const Highway &highway);
  void removeHighway();
  void saveSnapshotPerFamilyLL();
  void savePerFamilyLikelihoodDiff(const std::string &outputFile);
  unsigned int getInputTreesNumber() const;

  /**
   *  Is the recmodel dated?
   */
  virtual bool isDated() const { return _info.isDated(); }

  /**
   *  Are we in the prune species tree mode?
   */
  virtual bool pruneSpeciesTree() const { return _info.pruneSpeciesTree; }

  /**
   *  Should the optimization routines print verbose logs?
   */
  virtual bool isVerbose() const { return _optimizeVerbose; }

  /**
   *  Accessors
   */
  // general info
  unsigned int getLocalFamilyNumber() const {
    return _geneTrees.getTrees().size();
  }
  const RecModelInfo &getRecModelInfo() const { return _info; }
  // model parameter optimization classes
  const OptimizationClasses &getOptimizationClasses() const {
    return _optimizationClasses;
  }
  // evaluation of the i-th family
  MultiEvaluation &getEvaluation(unsigned int i) { return *_evaluations[i]; }

private:
  // optimize the alpha parameter for speciation probability gamma categories
  double optimizeGammaRates();
  // set precision of the evaluation of the i-th family
  void resetEvaluation(unsigned int i, bool highPrecision);
  // try to set all evaluations to low precision
  void resetAllPrecisions();
  // print the number of evaluations having high precision
  void printHighPrecisionCount();

private:
  AleOptimizer &_optimizer;
  SpeciesTree &_speciesTree;
  const RecModelInfo &_info;
  // _optimizationClasses["D"][5] is the index of the optimization class
  // of the duplication probability of species branch with index 5.
  // All parameters that belong to the same optimization class share
  // the same value
  const OptimizationClasses _optimizationClasses;
  double &_mixtureAlpha;
  std::vector<Highway> &_transferHighways;
  std::vector<AleModelParameters> &_modelParameters;
  bool _optimizeRates;
  bool _optimizeVerbose;
  // global families
  const Families &_families;
  // local families
  const PerCoreGeneTrees &_geneTrees;
  // evaluations for local families
  PerCoreMultiEvaluations _evaluations;
  PerCoreMultiEvaluations _approxEvaluations;
  // species branches (node indices) carrying a declared WGD, and the current
  // retention probability of each (aligned with _wgdNodes; used as the warm
  // start for the retention optimizer). Rebuilt from _wgdBranchLabels /
  // _wgdRetentions / _wgdResolutionsState on construction (e.g. from a
  // checkpoint), since node indices are not stable across a species tree
  // reload.
  std::vector<unsigned int> _wgdNodes;
  std::vector<double> _wgdQ;
  // LORe (delayed rediploidization): per-event resolution probability r (one per
  // declared WGD, aligned with _wgdNodes; default 1.0 == AORe), the species
  // branches each WGD governs (its subtree), and per WGD whether r is a free
  // parameter (false for terminal-branch WGDs: r pinned to AORe). Built lazily
  // by buildWGDStructure(). _resolutionProb is the legacy scalar (r of the first
  // resolvable WGD).
  std::vector<double> _wgdResolution;
  std::vector<std::vector<unsigned int>> _wgdSubtreeBranches;
  std::vector<char> _wgdResolvable;
  // node indices of the WGDs whose r is fitted under --lore-wgd (empty => every
  // internal-branch WGD is resolvable, i.e. the default --lore behaviour).
  // Populated post-construction by declareWGDs() from the (checkpoint-command-
  // guaranteed-identical) CLI args, so it is not itself checkpointed.
  std::vector<unsigned int> _loreTargets;
  // AleState-backed, label-keyed mirrors of _wgdNodes / _wgdQ / _wgdResolution:
  // this is what actually gets checkpointed, so every setWGD() /
  // setWGDResolutions() call keeps it in sync.
  std::vector<std::string> &_wgdBranchLabels;
  std::vector<double> &_wgdRetentions;
  std::vector<double> &_wgdResolutionsState;
  double _resolutionProb = 1.0;
  // Reference into AleState so that it is checkpointed like the other model
  // parameters.
  bool &_optimizeResolution;
  // suppress the bulky per-family output (--summary-only)
  bool _summaryOnly = false;
  std::vector<int> _highPrecisions;
  // LL buffer for global families
  std::vector<double> _snapshotPerFamilyLL;
};
