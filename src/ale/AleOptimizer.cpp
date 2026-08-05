#include <cstdio>
#include "AleOptimizer.hpp"

#include <algorithm>
#include <cassert>
#include <memory>

#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <IO/ParallelOfstream.hpp>
#include <maths/Random.hpp>
#include <parallelization/ParallelContext.hpp>
#include <search/DatedSpeciesTreeSearch.hpp>
#include <search/SpeciesRootSearch.hpp>
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>
#include <util/Paths.hpp>
#include <util/Scenario.hpp>

AleOptimizer::AleOptimizer(const std::string &speciesTreeFile,
                           const Families &families, const RecModelInfo &info,
                           const ModelParametrization &modelParametrization,
                           const std::string &optimizationClassFile,
                           const Parameters &startingRates, bool optimizeRates,
                           bool optimizeVerbose, const std::string &outputDir)
    : _families(families),               // global families
      _geneTrees(families, false, true), // init local families
      _info(info), _rootLikelihoods(getLocalFamilyNumber()),
      _outputDir(outputDir), _enableCheckpoints(true) {
  // init the current state (from checkpoint or de novo)
  _state.speciesTree = std::make_unique<SpeciesTree>(
      speciesTreeFile, true, getRecModelInfo().isDated());
  if (checkpointExists()) {
    Logger::info << "Loading checkpoint information..." << std::endl;
    for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
      _state.localFamilyNames.push_back(_geneTrees.getTrees()[i].name);
    }
    _state.unserialize(getCheckpointDir(_outputDir));
  } else {
    for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
      _state.localFamilyNames.push_back(_geneTrees.getTrees()[i].name);
      _state.perLocalFamilyModelParams.push_back(AleModelParameters(
          startingRates, getSpeciesTree().getTree().getNodeNumber()));
    }
  }
  // init the species tree interface
  _speciesTreeSearchState = std::make_unique<SpeciesSearchState>(
      getSpeciesTree(), speciesTreeFile, getLocalFamilyNumber());
  _speciesTreeSearchState->addListener(this);
  getSpeciesTree().addListener(this);
  ParallelContext::barrier();
  // init the evaluator
  _evaluator = std::make_unique<AleEvaluator>(
      *this, getSpeciesTree(), getRecModelInfo(), modelParametrization,
      optimizationClassFile, getMixtureAlpha(), getTransferHighways(),
      getModelParameters(), getWGDBranchLabels(), getWGDRetentions(),
      getWGDResolutions(), getOptimizeResolution(), optimizeRates,
      optimizeVerbose, _families, _geneTrees);
  Logger::timed << "Initial ll=" << getEvaluator().computeLikelihood()
                << std::endl;
}

void AleOptimizer::randomizeRoot() {
  assert(getCurrentStep() == AleStep::SpeciesTreeOpt);
  auto N = getSpeciesTree().getTree().getInnerNodeNumber();
  for (unsigned int i = 0; i < N; ++i) {
    auto direction = Random::getInt() % 4;
    if (SpeciesTreeOperator::canChangeRoot(getSpeciesTree(), direction)) {
      SpeciesTreeOperator::changeRoot(getSpeciesTree(), direction);
    }
  }
}

void AleOptimizer::optimize() {
  assert(getCurrentStep() == AleStep::SpeciesTreeOpt);
  PerFamLL initialPerFamLL;
  auto ll = getEvaluator().computeLikelihood(&initialPerFamLL);
  _speciesTreeSearchState->bestLL = ll;
  _speciesTreeSearchState->farFromPlausible = true;
  _speciesTreeSearchState->betterTreeCallback(ll, initialPerFamLL);
  // Initial root search
  Logger::info << std::endl;
  rootSearch(3);
  // Alternate transfer search and normal SPR search,
  // until one cannot find a better tree.
  // Run each at least once
  size_t hash1 = getSpeciesTree().getHash();
  size_t hash2 = 0;
  unsigned int index = 0;
  do {
    if (index++ % 2 == 0) {
      Logger::info << std::endl;
      transferSearch();
    } else {
      Logger::info << std::endl;
      sprSearch(1);
    }
    if (!_speciesTreeSearchState->farFromPlausible) {
      Logger::info << std::endl;
      rootSearch(3);
    }
    hash2 = getSpeciesTree().getHash();
    std::swap(hash1, hash2);
  } while (hash1 != hash2);
  // Final root search
  Logger::info << std::endl;
  rootSearch(-1);
}

void AleOptimizer::reroot() {
  assert(getCurrentStep() == AleStep::SpeciesTreeOpt);
  PerFamLL initialPerFamLL;
  auto ll = getEvaluator().computeLikelihood(&initialPerFamLL);
  _speciesTreeSearchState->bestLL = ll;
  _speciesTreeSearchState->farFromPlausible = true;
  _speciesTreeSearchState->betterTreeCallback(ll, initialPerFamLL);
  // Optimize rates and run a local root search.
  // If the species tree is dated, run the second
  // local root search with thorough date optimization
  Logger::info << std::endl;
  optimizeModelRates(true);
  if (_info.isDated()) {
    Logger::info << std::endl;
    Logger::timed << "First root search, non-thorough dating" << std::endl;
    rootSearch(5, false);
    Logger::info << std::endl;
    Logger::timed << "Second root search, thorough dating" << std::endl;
    rootSearch(2, true);
  } else {
    Logger::info << std::endl;
    rootSearch(5);
  }
}

double AleOptimizer::rootSearch(unsigned int maxDepth, bool thorough) {
  _rootLikelihoods.reset();
  if (thorough) {
    _speciesTreeSearchState->farFromPlausible = false;
  }
  SpeciesRootSearch::rootSearch(getSpeciesTree(), getEvaluator(),
                                *_speciesTreeSearchState, maxDepth,
                                &_rootLikelihoods);
  auto ll = _speciesTreeSearchState->bestLL;
  saveCheckpoint();
  saveSpeciesRootSupports();
  return ll;
}

double AleOptimizer::transferSearch() {
  SpeciesTransferSearch::transferSearch(getSpeciesTree(), getEvaluator(),
                                        *_speciesTreeSearchState);
  auto ll = _speciesTreeSearchState->bestLL;
  Logger::timed << "After normal search: LL=" << ll << std::endl;
  saveCheckpoint();
  return ll;
}

double AleOptimizer::sprSearch(unsigned int radius) {
  SpeciesSPRSearch::SPRSearch(getSpeciesTree(), getEvaluator(),
                              *_speciesTreeSearchState, radius);
  auto ll = _speciesTreeSearchState->bestLL;
  Logger::timed << "After normal search: LL=" << ll << std::endl;
  saveCheckpoint();
  saveSpeciesBranchSupports();
  return ll;
}

double AleOptimizer::optimizeModelRates(bool thorough) {
  getEvaluator().optimizeModelRates(thorough);
  PerFamLL perFamLL;
  auto ll = getEvaluator().computeLikelihood(&perFamLL);
  _speciesTreeSearchState->betterLikelihoodCallback(ll, perFamLL);
  saveCheckpoint();
  return ll;
}

void AleOptimizer::optimizeDates(bool thorough) {
  assert(getCurrentStep() == AleStep::RelDating);
  if (!_info.isDated()) {
    return;
  }
  auto initialLL = _speciesTreeSearchState->bestLL;
  unsigned int datingsNumber =
      500; // number of independent searches from random datings
  unsigned int backupsNumber =
      20; // max number of best datings to continue with
  Logger::timed << "[Species search] Inferring speciation order for "
                << datingsNumber << " random datings" << std::endl;
  auto scoredBackups = DatedSpeciesTreeSearch::getBestDatingsFromReconciliation(
      getSpeciesTree(), getEvaluator(), datingsNumber, backupsNumber);
  Logger::timed << "[Species search] Best " << backupsNumber
                << " dating likelihoods from random datings:" << std::endl;
  for (const auto &sb : scoredBackups) {
    Logger::info << sb.score << " ";
  }
  Logger::info << std::endl;
  auto bestLL = scoredBackups[0].score;
  auto llDiff = bestLL - initialLL;
  if (llDiff > 0.0) {
    Logger::timed << "[Species search] Accept improvement from random datings: "
                  << "bestLL=" << bestLL << ", llDiff=" << llDiff << std::endl;
    SpeciesTreeOperator::restoreDates(getSpeciesTree(),
                                      scoredBackups[0].backup);
  } else {
    Logger::timed << "[Species search] No improvement from random datings: "
                  << "bestLL=" << bestLL << ", llDiff=" << llDiff << std::endl;
  }
  DatedSpeciesTreeSearch::optimizeDates(getSpeciesTree(), getEvaluator(),
                                        *_speciesTreeSearchState, thorough);
  saveCheckpoint();
  saveSpeciesTree();
  Logger::timed << "Species tree with ordered speciations: "
                << _speciesTreeSearchState->pathToBestSpeciesTree << std::endl;
}

void AleOptimizer::onSpeciesDatesChange() {
  getEvaluator().onSpeciesDatesChange();
}

void AleOptimizer::onSpeciesTreeChange(
    const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
  getEvaluator().onSpeciesTreeChange(nodesToInvalidate);
}

void AleOptimizer::betterTreeCallback() { saveCheckpoint(); }

void AleOptimizer::onBetterParametersFoundCallback() {
  if (_enableCheckpoints) {
    saveCheckpoint();
  }
}

void AleOptimizer::saveSpeciesTree() {
  auto out = _speciesTreeSearchState->pathToBestSpeciesTree;
  getSpeciesTree().saveToFile(out, true);
  // Logger::info << "save species tree to " << out << std::endl;
  ParallelContext::barrier();
}

void AleOptimizer::saveSpeciesBranchSupports() {
  auto outKH =
      Paths::getSpeciesTreeFile(_outputDir, "species_tree_support_kh.newick");
  _speciesTreeSearchState->saveSpeciesTreeKH(outKH);
  // Logger::info << "save KH branch support tree to " << outKH << std::endl;
  auto outBP =
      Paths::getSpeciesTreeFile(_outputDir, "species_tree_support_bp.newick");
  _speciesTreeSearchState->saveSpeciesTreeBP(outBP);
  // Logger::info << "save BP branch support tree to " << outBP << std::endl;
}

void AleOptimizer::saveSpeciesRootSupports() {
  assert(!_rootLikelihoods.isEmpty());
  auto newick = getSpeciesTree().getTree().getNewickString();
  PLLRootedTree treeLLR(newick, false);
  _rootLikelihoods.fillTree(treeLLR);
  auto outLLR =
      Paths::getSpeciesTreeFile(_outputDir, "species_tree_root_llr.newick");
  treeLLR.save(outLLR);
  // Logger::info << "save LLR root support tree to " << outLLR << std::endl;
  PLLRootedTree treeRELL(newick, false);
  _rootLikelihoods.fillTreeBootstraps(treeRELL);
  auto outRELL =
      Paths::getSpeciesTreeFile(_outputDir, "species_tree_root_rell.newick");
  treeRELL.save(outRELL);
  // Logger::info << "save RELL root support tree to " << outRELL << std::endl;
}

void AleOptimizer::saveRatesAndLL() {
  // save per-family likelihoods
  auto perFamilyLikelihoodPath =
      FileSystem::joinPaths(_outputDir, "per_fam_likelihoods.txt");
  std::vector<unsigned int> localIndices;
  std::vector<double> localLikelihoods;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    auto ll = _evaluator->computeFamilyLikelihood(i);
    localIndices.push_back(famIndex);
    localLikelihoods.push_back(ll);
  }
  ParallelContext::barrier();
  std::vector<unsigned int> indices;
  std::vector<double> likelihoods;
  ParallelContext::concatenateHeterogeneousUIntVectors(localIndices, indices);
  ParallelContext::concatenateHeterogeneousDoubleVectors(localLikelihoods,
                                                         likelihoods);
  assert(indices.size() == _families.size());
  if (ParallelContext::getRank() == 0) {
    std::vector<ScoredFamily> familiesAndLLs;
    for (unsigned int i = 0; i < indices.size(); ++i) {
      familiesAndLLs.push_back(
          ScoredFamily(_families[indices[i]].name, likelihoods[i]));
    }
    std::sort(familiesAndLLs.begin(), familiesAndLLs.end());
    std::ofstream llOs(perFamilyLikelihoodPath);
    for (const auto &sf : familiesAndLLs) {
      llOs << sf.familyName << " " << sf.score << std::endl;
    }
    llOs.close();
  }
  // Logger::info << "save per-family likelihoods to " <<
  // perFamilyLikelihoodPath << std::endl;
  //  save the DTL rates
  auto parameterNames = _info.getParamTypes();
  auto ratesDir = FileSystem::joinPaths(_outputDir, "model_parameters");
  FileSystem::mkdir(ratesDir, true);
  ParallelContext::barrier();
  if (_info.perFamilyRates) {
    for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
      auto familyRatesPath = FileSystem::joinPaths(
          ratesDir, _geneTrees.getTrees()[i].name + "_rates.txt");
      std::ofstream ratesOs(familyRatesPath);
      for (const auto &name : parameterNames) {
        if (name != parameterNames[0])
          ratesOs << " ";
        ratesOs << name;
      }
      ratesOs << std::endl;
      const auto &parameters = getModelParameters()[i];
      assert(parameters.getParamTypeNumber() == parameterNames.size());
      for (unsigned int rate = 0; rate < parameterNames.size(); ++rate) {
        if (rate != 0)
          ratesOs << " ";
        ratesOs << parameters.getParameter(0, rate);
      }
      ratesOs << std::endl;
      ratesOs.close();
    }
  } else {
    auto globalRatesPath =
        FileSystem::joinPaths(ratesDir, "model_parameters.txt");
    ParallelOfstream ratesOs(globalRatesPath, true);
    ratesOs << "node";
    for (const auto &name : parameterNames) {
      ratesOs << " " << name;
    }
    ratesOs << std::endl;
    const auto &parameters = getModelParameters()[0];
    assert(parameters.getParamTypeNumber() == parameterNames.size());
    for (auto node : getSpeciesTree().getTree().getNodes()) {
      ratesOs << node->label;
      for (unsigned int rate = 0; rate < parameterNames.size(); ++rate) {
        ratesOs << " " << parameters.getParameter(node->node_index, rate);
      }
      ratesOs << std::endl;
    }
    ratesOs.close();
  }
  // Logger::info << "save model rates to " << ratesDir << std::endl;
  ParallelContext::barrier();
}

void AleOptimizer::saveGeneConsensusTree(const std::string &geneSampleFile,
                                         const std::string &outputFile) {
  std::ifstream is(geneSampleFile);
  std::string line;
  std::vector<std::string> rtreeStrs;
  while (std::getline(is, line)) {
    if (line.size() > 2) {
      rtreeStrs.push_back(line);
    }
  }
  is.close();
  if (rtreeStrs.empty()) {
    return; // no sampled gene trees for this family (degenerate at high q)
  }
  auto consRtreeStr = PLLRootedTree::buildConsensusTree(rtreeStrs, 0.50001);
  if (consRtreeStr.empty()) {
    return; // corax could not form a consensus; skip rather than crash
  }
  std::ofstream os(outputFile);
  os << consRtreeStr << std::endl;
  os.close();
}

void AleOptimizer::saveFamiliesTakingHighway(
    const Highway &highway, const VectorDouble &perFamilyTransfers,
    const std::string &directory) {
  std::vector<unsigned int> localIndices;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    localIndices.push_back(famIndex);
  }
  ParallelContext::barrier();
  std::vector<unsigned int> indices;
  std::vector<double> transfers;
  ParallelContext::concatenateHeterogeneousUIntVectors(localIndices, indices);
  ParallelContext::concatenateHeterogeneousDoubleVectors(perFamilyTransfers,
                                                         transfers);
  assert(indices.size() == _families.size());
  if (ParallelContext::getRank() == 0) {
    std::vector<ScoredFamily> scoredFamilies;
    for (unsigned int i = 0; i < indices.size(); ++i) {
      scoredFamilies.push_back(
          ScoredFamily(_families[indices[i]].name, transfers[i]));
    }
    std::sort(scoredFamilies.rbegin(), scoredFamilies.rend());
    auto outputFile = FileSystem::joinPaths(
        directory, std::string("highway_") + highway.src->label + "_to_" +
                       highway.dest->label + ".txt");
    std::ofstream os(outputFile);
    os << "fam, transfers" << std::endl;
    for (const auto &sf : scoredFamilies) {
      if (sf.score == 0.0) {
        break;
      }
      os << sf.familyName << ", " << sf.score << std::endl;
    }
    os.close();
  }
  ParallelContext::barrier();
}

void AleOptimizer::reconcile(unsigned int samples) {
  assert(getCurrentStep() == AleStep::Reconciliation);
  if (samples == 0) {
    return;
  }
  Logger::timed << "[Reconciliation] Sampling reconciled gene trees ("
                << samples << " samples per gene family)" << std::endl;
  // Initiating directories
  auto recDir = FileSystem::joinPaths(_outputDir, "reconciliations");
  FileSystem::mkdir(recDir, true);
  auto allRecDir = FileSystem::joinPaths(recDir, "all");
  FileSystem::mkdir(allRecDir, true);
  auto summariesDir = FileSystem::joinPaths(recDir, "summaries");
  FileSystem::mkdir(summariesDir, true);
  auto originsDir = FileSystem::joinPaths(recDir, "origins");
  FileSystem::mkdir(originsDir, true);
  auto highwaysOutputDir = FileSystem::joinPaths(_outputDir, "highways");
  auto highwayFamiliesOutputDir =
      FileSystem::joinPaths(highwaysOutputDir, "accepted_transfers");
  if (FileSystem::dirExists(highwaysOutputDir)) {
    FileSystem::mkdir(highwayFamiliesOutputDir, true);
  }
  ParallelContext::barrier();
  Logger::timed << "[Reconciliation] Inferring reconciliations" << std::endl;
  const auto &localFamilies = _geneTrees.getTrees();
  std::vector<std::string> summaryPerSpeciesEventCountsFiles;
  std::vector<std::string> summaryTransferFiles;
  std::vector<std::shared_ptr<Scenario>> allScenarios;
  MatrixDouble perHighwayPerFamTransfers;
  const auto &highways = getTransferHighways();
  if (highways.size() && localFamilies.size()) {
    assert(FileSystem::dirExists(highwayFamiliesOutputDir));
    perHighwayPerFamTransfers =
        MatrixDouble(highways.size(), VectorDouble(localFamilies.size(), 0.0));
  }
  // LORe: per-species WGD resolution-event (U->R commit) profile, accumulated
  // over families. For the DL or DTL model with a declared WGD, full (unpruned)
  // species tree, and CLVs kept in RAM (no memory savings). See
  // WGD_LORE_marginal.md STEP 4.
  const auto &recInfo = _evaluator->getRecModelInfo();
  bool doResolution = !_evaluator->getWGDNodes().empty() &&
                      (recInfo.model == RecModel::UndatedDL ||
                       recInfo.model == RecModel::UndatedDTL) &&
                      !recInfo.pruneSpeciesTree && !recInfo.memorySavings;
  std::vector<double> totalResolution(getSpeciesTree().getTree().getNodeNumber(),
                                      0.0);
  // Per-(leaf)-species expected number of ohnolog pairs still TETRASOMIC at the
  // tips (a U lineage that reaches an extant species without ever committing).
  std::vector<double> totalTetrasomy(getSpeciesTree().getTree().getNodeNumber(),
                                     0.0);
  // The gene-tree event backbone (D/T/L/copies + sampled trees) is sampled at
  // r=1 (AORe-resolved). Under LORe (r<1) the WGD-doubling / U-state mass would
  // otherwise leak into the *U-unaware* scenario sampler as a spray of spurious
  // transfers around the WGD branch (the model treats the U-state as
  // vertical-only, transfers acting on resolved lineages, but the scenario
  // backtrace does not). The LORe resolution profile is sampled separately,
  // below, at the fitted r (PASS B). sampleReconciliations() calls updateCLVs(),
  // so setting r here is enough to make the backbone use the resolved transform.
  // Save the fitted per-event r (one per declared WGD) so we can restore it
  // after sampling the backbone at r=1. LORe is active if any WGD has r<1.
  std::vector<double> fittedR;
  for (auto w : _evaluator->getWGDNodes()) {
    fittedR.push_back(_evaluator->getWGDResolution(w));
  }
  bool wgdActiveLORe = false;
  for (auto r : fittedR) {
    if (r < 1.0) {
      wgdActiveLORe = true;
      break;
    }
  }
  if (wgdActiveLORe) {
    _evaluator->setResolutionProb(1.0);
  }
  for (unsigned int i = 0; i < localFamilies.size(); ++i) {
    std::vector<std::string> perSpeciesEventCountsFiles;
    std::vector<std::string> transferFiles;
    std::vector<std::shared_ptr<Scenario>> scenarios;
    // Warning:
    // Using Random::getProba() in the sampling function makes
    // the random state inconsistent between the MPI ranks.
    // Call ParallelContext::makeRandConsistent() right after
    // all MPI ranks passed the loop
    _evaluator->sampleFamilyScenarios(i, samples, scenarios);
    allScenarios.insert(allScenarios.end(), scenarios.begin(), scenarios.end());
    // Usually scenarios.size() == samples. A few families are not samplable
    // under the DTL+WGD/LORe state (the scenario backtrace hits the depth guard
    // on every attempt; see MultiModel::sampleReconciliations) and yield fewer
    // (or no) gene-tree samples. We iterate over what was sampled; the LORe
    // resolution profile below is computed by a separate, robust backtrace.
    unsigned int sampledScenarios = scenarios.size();
    // (LORe resolution profile is sampled in PASS B below, at the fitted r.)
    // writing in the reconciliations/all/ dir. Under --summary-only we suppress
    // the bulky per-sample output (RecPhyloXML, .alerec, per-sample total event
    // counts) entirely, and write only what the kept per-family summaries are
    // built from (the sampled newick for the consensus, and the per-species
    // event / transfer counts for the means); those transient files are deleted
    // at the end of this family's iteration so peak disk stays at one family.
    const bool summaryOnly = _evaluator->isSummaryOnly();
    auto geneTreesPath = FileSystem::joinPaths(
        allRecDir, localFamilies[i].name + "_samples.newick");
    ParallelOfstream geneTreesOs(geneTreesPath, false);
    auto geneTreesAlePath = FileSystem::joinPaths(
        allRecDir, localFamilies[i].name + "_samples.alerec");
    std::unique_ptr<ParallelOfstream> geneTreesAleOs;
    if (!summaryOnly) {
      geneTreesAleOs.reset(new ParallelOfstream(geneTreesAlePath, false));
    }
    for (unsigned int sample = 0; sample < sampledScenarios; ++sample) {
      auto perSpeciesEventCountsFile = FileSystem::joinPaths(
          allRecDir, localFamilies[i].name + "_speciesEventCounts_" +
                         std::to_string(sample) + ".txt");
      auto transferFile = FileSystem::joinPaths(
          allRecDir, localFamilies[i].name + "_transfers_" +
                         std::to_string(sample) + ".txt");
      perSpeciesEventCountsFiles.push_back(perSpeciesEventCountsFile);
      transferFiles.push_back(transferFile);
      auto &scenario = *scenarios[sample];
      scenario.saveReconciliation(geneTreesOs,
                                  ReconciliationFormat::NewickEvents);
      if (!summaryOnly) {
        auto geneTreeXMLPath = FileSystem::joinPaths(
            allRecDir, localFamilies[i].name + "_sample_" +
                           std::to_string(sample) + ".xml");
        auto eventCountsFile = FileSystem::joinPaths(
            allRecDir, localFamilies[i].name + "_eventCounts_" +
                           std::to_string(sample) + ".txt");
        scenario.saveReconciliation(*geneTreesAleOs, ReconciliationFormat::ALE);
        scenario.saveReconciliation(geneTreeXMLPath,
                                    ReconciliationFormat::RecPhyloXML, false);
        scenario.saveEventsCounts(eventCountsFile, false);
      }
      scenario.savePerSpeciesEventsCounts(perSpeciesEventCountsFile, false);
      scenario.saveTransfers(transferFile, false);
      for (unsigned int hi = 0; hi < highways.size(); ++hi) {
        perHighwayPerFamTransfers[hi][i] += scenario.countTransfer(
            highways[hi].src->label, highways[hi].dest->label);
      }
    }
    for (unsigned int hi = 0; hi < highways.size(); ++hi) {
      perHighwayPerFamTransfers[hi][i] /= static_cast<double>(samples);
    }
    geneTreesOs.close();
    if (geneTreesAleOs) {
      geneTreesAleOs->close();
    }
    // writing in the reconciliations/summaries/ dir (skip for the rare families
    // with no sampled gene tree -- see the sampledScenarios note above)
    if (sampledScenarios > 0) {
      auto consensusFile = FileSystem::joinPaths(
          summariesDir, localFamilies[i].name + "_consensus_50.newick");
      saveGeneConsensusTree(geneTreesPath, consensusFile);
      auto perSpeciesEventCountsFile = FileSystem::joinPaths(
          summariesDir, localFamilies[i].name + "_meanSpeciesEventCounts.txt");
      Scenario::mergePerSpeciesEventCounts(
          getSpeciesTree().getTree(), perSpeciesEventCountsFile,
          perSpeciesEventCountsFiles, false, true);
      summaryPerSpeciesEventCountsFiles.push_back(perSpeciesEventCountsFile);
      auto transferFile = FileSystem::joinPaths(
          summariesDir, localFamilies[i].name + "_meanTransfers.txt");
      Scenario::mergeTransfers(getSpeciesTree().getTree(), transferFile,
                               transferFiles, false, true);
      summaryTransferFiles.push_back(transferFile);
    }
    // --summary-only: the per-family summaries are now built, so drop the
    // transient per-sample files (each rank wrote its own local families).
    if (summaryOnly) {
      std::remove(geneTreesPath.c_str());
      for (const auto &f : perSpeciesEventCountsFiles) {
        std::remove(f.c_str());
      }
      for (const auto &f : transferFiles) {
        std::remove(f.c_str());
      }
    }
  }
  ParallelContext::barrier();
  ParallelContext::makeRandConsistent();
  // Restore the fitted per-event r (the event backbone above was sampled at
  // r=1). setWGDResolutions repaints each WGD's subtree with its own r.
  if (wgdActiveLORe) {
    _evaluator->setWGDResolutions(fittedR);
  }
  // PASS B: LORe resolution-branch profile, sampled at the fitted r by the
  // U-aware backtrace. sampleResolutionCommits() does NOT refresh the CLVs, so
  // rebuild them at the fitted r first (the backbone pass left them at r=1).
  if (doResolution) {
    if (wgdActiveLORe) {
      _evaluator->computeLikelihood();
    }
    for (unsigned int i = 0; i < localFamilies.size(); ++i) {
      std::vector<double> famCommits;
      std::vector<double> famTetra;
      _evaluator->sampleFamilyResolutionCommits(i, samples, famCommits,
                                                famTetra);
      auto famResFile = FileSystem::joinPaths(
          summariesDir, localFamilies[i].name + "_meanResolutionCounts.txt");
      ParallelOfstream famOs(famResFile, false);
      famOs << "#species_branch\texpected_WGD_resolution_events"
               "\texpected_tetrasomic_at_tip"
            << std::endl;
      for (auto node : getSpeciesTree().getTree().getNodes()) {
        double expected =
            famCommits.empty()
                ? 0.0
                : famCommits[node->node_index] / static_cast<double>(samples);
        double tetra =
            famTetra.empty()
                ? 0.0
                : famTetra[node->node_index] / static_cast<double>(samples);
        totalResolution[node->node_index] += expected;
        totalTetrasomy[node->node_index] += tetra;
        famOs << (node->label ? node->label : "NA") << "\t" << expected << "\t"
              << tetra << std::endl;
      }
      famOs.close();
    }
    ParallelContext::barrier();
    ParallelContext::makeRandConsistent();
  }
  Logger::timed << "[Reconciliation] Exporting reconciliation summaries"
                << std::endl;
  // export total per-branch event counts
  auto totalPerSpeciesEventCountsFile =
      FileSystem::joinPaths(recDir, "totalSpeciesEventCounts.txt");
  Scenario::mergePerSpeciesEventCounts(
      getSpeciesTree().getTree(), totalPerSpeciesEventCountsFile,
      summaryPerSpeciesEventCountsFiles, true, false);
  // export origins
  Scenario::saveOriginsGlobal(getSpeciesTree().getTree(), allScenarios, samples,
                              originsDir);
  // export total pairwise transfer counts
  auto totalTransferFile = FileSystem::joinPaths(recDir, "totalTransfers.txt");
  Scenario::mergeTransfers(getSpeciesTree().getTree(), totalTransferFile,
                           summaryTransferFiles, true, false);
  // export the total per-branch WGD resolution-event (U->R commit) profile:
  // the expected number of ohnolog-divergence (rediploidization) commits on
  // each species branch, summed over all gene families, under the fitted r.
  if (doResolution) {
    ParallelContext::sumVectorDouble(totalResolution);
    auto totalResolutionFile =
        FileSystem::joinPaths(recDir, "totalSpeciesResolutionCounts.txt");
    ParallelOfstream os(totalResolutionFile, true); // master rank only
    os << "#species_branch\texpected_WGD_resolution_events" << std::endl;
    for (auto node : getSpeciesTree().getTree().getNodes()) {
      os << (node->label ? node->label : "NA") << "\t"
         << totalResolution[node->node_index] << std::endl;
    }
    os.close();
    Logger::timed << "[Reconciliation] WGD resolution-branch profile: "
                  << totalResolutionFile << std::endl;
    // export the total per-tip tetrasomy profile: the expected number of ohnolog
    // pairs that reach an extant species still unresolved (tetrasomic today),
    // summed over all gene families, under the fitted r.
    ParallelContext::sumVectorDouble(totalTetrasomy);
    auto totalTetraFile =
        FileSystem::joinPaths(recDir, "totalSpeciesTetrasomyCounts.txt");
    ParallelOfstream tos(totalTetraFile, true); // master rank only
    tos << "#species_branch\texpected_tetrasomic_pairs" << std::endl;
    for (auto node : getSpeciesTree().getTree().getNodes()) {
      tos << (node->label ? node->label : "NA") << "\t"
          << totalTetrasomy[node->node_index] << std::endl;
    }
    tos.close();
    Logger::timed << "[Reconciliation] tetrasomy (still-unresolved) profile: "
                  << totalTetraFile << std::endl;
    // Per-WGD summary: declared branch, fitted retention q and resolution r, and
    // the expected number of rediploidization commits within the WGD's subtree.
    // This makes the WGD location self-contained even under LORe (where the
    // resolution profile above points at descendant branches, not the WGD one).
    auto &stree = getSpeciesTree().getTree();
    auto wgdSummaryFile = FileSystem::joinPaths(recDir, "wgdSummary.txt");
    ParallelOfstream wos(wgdSummaryFile, true);
    wos << "#wgd_branch\tnode_index\tq_retention\tr_resolution"
        << "\texpected_resolution_commits\texpected_tetrasomic_pairs"
        << std::endl;
    for (auto w : _evaluator->getWGDNodes()) {
      double subtreeCommits = 0.0;
      double subtreeTetra = 0.0;
      for (auto node : stree.getNodes()) {
        auto e = node->node_index;
        if (e == w || stree.isAncestorOf(w, e)) {
          subtreeCommits += totalResolution[e];
          subtreeTetra += totalTetrasomy[e];
        }
      }
      auto *wnode = stree.getNode(w);
      wos << (wnode->label ? wnode->label : "NA") << "\t" << w << "\t"
          << _evaluator->getWGDRetention(w) << "\t"
          << _evaluator->getWGDResolution(w) << "\t" // per-event r (was rHat)
          << subtreeCommits << "\t" << subtreeTetra << std::endl;
    }
    wos.close();
    Logger::timed << "[Reconciliation] WGD summary (branch, q, r, commits): "
                  << wgdSummaryFile << std::endl;
  }
  // export highway information
  for (unsigned int hi = 0; hi < highways.size(); ++hi) {
    saveFamiliesTakingHighway(highways[hi], perHighwayPerFamTransfers[hi],
                              highwayFamiliesOutputDir);
  }
  assert(ParallelContext::isRandConsistent());
  Logger::timed << "Reconciliation output directory: " << recDir << std::endl;
}

void AleOptimizer::saveBestHighways(
    const std::vector<ScoredHighway> &scoredHighways,
    const std::string &outputFile) {
  // Logger::info << "save highways to " << outputFile << std::endl;
  assert(scoredHighways.size());
  ParallelOfstream os(outputFile, true);
  os << "src, dest, proba, llDiff" << std::endl;
  for (const auto &scoredHighway : scoredHighways) {
    os << scoredHighway.highway.src->label << ", "
       << scoredHighway.highway.dest->label << ", "
       << scoredHighway.highway.proba << ", " << scoredHighway.score
       << std::endl;
  }
  os.close();
  ParallelContext::barrier();
}

void AleOptimizer::inferHighways(const std::string &highwayCandidateFile,
                                 unsigned int highwayCandidatesStep1,
                                 unsigned int highwayCandidatesStep2) {
  // let's infer highways of transfers!
  assert(getCurrentStep() == AleStep::Highways);
  auto highwaysOutputDir = getHighwaysOutputDir();
  // Step 1: select initial candidates
  auto candidateHighwayOutput =
      FileSystem::joinPaths(highwaysOutputDir, "candidate_highways.txt");
  std::vector<ScoredHighway> candidateHighways;
  if (highwayCandidateFile.size()) {
    // the user sets the candidates
    auto highways = HighwayCandidateParser::parse(highwayCandidateFile,
                                                  getSpeciesTree().getTree());
    Highways::getSortedCandidatesFromList(*this, highways, candidateHighways);
  } else {
    // automatically search for candidates
    Highways::getCandidateHighways(*this, candidateHighways,
                                   highwayCandidatesStep1);
  }
  if (!candidateHighways.size()) {
    Logger::timed << "No candidate highways found!" << std::endl;
    return;
  }
  FileSystem::mkdir(highwaysOutputDir, true);
  saveBestHighways(candidateHighways, candidateHighwayOutput);
  // Step 2: candidate filtering. We add each highway candidate individually,
  // set a small highway probability and keep the highway if the likelihood
  // improves. We also sort the highways by likelihood and keep the
  // int(highwayCandidatesStep2) best of them
  std::vector<ScoredHighway> filteredHighways;
  Highways::filterCandidateHighways(*this, candidateHighways, filteredHighways,
                                    highwayCandidatesStep2);
  if (!filteredHighways.size()) {
    Logger::timed << "No candidate highways passed the filtering!" << std::endl;
    return;
  }
  // Step 3: optimize all the filtered highways together
  auto acceptedHighwayOutput =
      FileSystem::joinPaths(highwaysOutputDir, "accepted_highways.txt");
  std::vector<ScoredHighway> acceptedHighways;
  Highways::optimizeAllHighways(*this, filteredHighways, acceptedHighways,
                                true);
  assert(acceptedHighways.size()); // there must be some, since we can get here
  saveBestHighways(acceptedHighways, acceptedHighwayOutput);
  Logger::timed << "Highway output directory: " << highwaysOutputDir
                << std::endl;
}
