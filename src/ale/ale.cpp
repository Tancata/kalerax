#include <cassert>
#include <cstdio>
#include <set>
#include <unordered_map>

#include <DistanceMethods/MiniNJ.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/FileSystem.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <IO/ParallelOfstream.hpp>
#include <ccp/ConditionalClades.hpp>
#include <maths/Parameters.hpp>
#include <maths/Random.hpp>
#include <parallelization/ParallelContext.hpp>
#include <trees/SpeciesTree.hpp>
#include <util/Paths.hpp>
#include <util/RecModelInfo.hpp>
#include <util/enums.hpp>
#include <util/types.hpp>
#include <util/utils.hpp>

#include "AleArguments.hpp"
#include "AleOptimizer.hpp"
#include "TrimFamilies.hpp"

const char *version = "kalerax v1.4.1 (WGD/LORe fork of AleRax)";

/**
 *  Create AleRax top directories
 */
void initAleRaxDirectories(const AleArguments &args,
                           const std::string &ccpDir) {
  FileSystem::mkdir(args.output, true);
  FileSystem::mkdir(args.output + "/species_trees", true);
  FileSystem::mkdir(ccpDir, true);
  ParallelContext::barrier();
}

/**
 *  Print AleRax version, the command line and a summary of
 *  the parameters
 */
void printInitialMessage(const AleArguments &args) {
  Logger::initFileOutput(FileSystem::joinPaths(args.output, "alerax"));
  Logger::timed << version << std::endl;
  args.printCommand();
  args.printSummary();
}

/**
 *  Check that there are enough gene families
 */
void checkEnoughFamilies(const Families &families) {
  if (families.size() == 0) {
    Logger::info << "\nError: No valid families, aborting\n" << std::endl;
    ParallelContext::abort(0);
  }
  if (families.size() < ParallelContext::getSize()) {
    Logger::info << "\nError: More MPI ranks set than valid families, "
                 << "aborting\n"
                 << std::endl;
    ParallelContext::abort(0);
  }
}

/**
 *  Check the validity of each gene family:
 *  - check that the gene tree files exist
 *  - if set, check that the mapping files exist
 */
void filterInvalidFamilies(const AleArguments &args, Families &families) {
  if (args.skipFamilyFiltering) {
    return;
  }
  Logger::timed << "Checking families..." << std::endl;
  Families validFamilies;
  for (const auto &family : families) {
    // check that the gene tree distribution file exists
    std::ifstream isTrees(family.startingGeneTree);
    if (!isTrees || isTrees.peek() == std::ifstream::traits_type::eof()) {
      Logger::info << "Excluding family " << family.name
                   << ": cannot open the gene tree distribution file!"
                   << std::endl;
      continue;
    }
    if (family.mappingFile.size()) {
      // check that the mapping file exists if it is set
      std::ifstream isMap(family.mappingFile);
      if (!isMap || isMap.peek() == std::ifstream::traits_type::eof()) {
        Logger::info << "Excluding family " << family.name
                     << ": cannot open the mapping file!" << std::endl;
        continue;
      }
    }
    validFamilies.push_back(family);
  }
  families = validFamilies;
  checkEnoughFamilies(families);
  ParallelContext::barrier();
}

/**
 *  Read the gene tree distribution files (or the .ale files), compute
 *  the corresponding CCP objects and serialize them to a binary file.
 *  Each MPI rank will run this function on a subset of the families
 *
 *  @param families The gene families to process
 *  @param ccpRooting Mode to decide how to deal with rooted gene trees
 *  @param sampleFrequency Number used to subsample the gene trees (if set to 2,
 *  we only use half)
 *  @param ccpDimensionFile Output file with the per-family ccp sizes
 */
void generateCCPs(const AleArguments &args, const Families &families,
                  const std::string &ccpDimensionFile) {
  ParallelContext::barrier();
  Logger::timed << "Generating ccp files..." << std::endl;
  auto N = families.size();
  std::vector<unsigned int> localFamilyIndices;
  std::vector<unsigned int> localTreeSizes;
  std::vector<unsigned int> localCcpSizes;
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N);
       ++i) {
    // read and seralize the ccps
    ConditionalClades ccp(families[i].startingGeneTree,
                          families[i].likelihoodFile, args.ccpRooting,
                          args.sampleFrequency);
    if (!ccp.isValid()) {
      Logger::error << "Error: invalid CCP for family " << families[i].name
                    << std::endl;
      ParallelContext::abort(1);
    }
    ccp.serialize(families[i].ccpFile);
    // record the dimensions for the subsequent sorting
    localFamilyIndices.push_back(i);
    localTreeSizes.push_back(ccp.getLeafNumber());
    localCcpSizes.push_back(ccp.getCladesNumber());
  }
  ParallelContext::barrier();
  // gather the family ccp dimensions from all MPI ranks
  std::vector<unsigned int> familyIndices;
  std::vector<unsigned int> treeSizes;
  std::vector<unsigned int> ccpSizes;
  ParallelContext::concatenateHeterogeneousUIntVectors(localFamilyIndices,
                                                       familyIndices);
  ParallelContext::concatenateHeterogeneousUIntVectors(localTreeSizes,
                                                       treeSizes);
  ParallelContext::concatenateHeterogeneousUIntVectors(localCcpSizes, ccpSizes);
  assert(familyIndices.size() == families.size());
  // output the families and their cpp sizes, from the largest to the smallest
  auto sortedIndices = sort_indices_descending<unsigned int>(ccpSizes);
  ParallelOfstream os(ccpDimensionFile, true);
  os << "fam, leaves, ccps" << std::endl;
  for (unsigned int i = 0; i < familyIndices.size(); ++i) {
    auto j = sortedIndices[i];
    os << families[familyIndices[j]].name << ", " << treeSizes[j] << ", "
       << ccpSizes[j] << std::endl;
  }
  os.close();
  ParallelContext::barrier();
}

/**
 *  Remove the serialized CCP files to free some disk space
 */
void cleanupCCPs(const Families &families) {
  ParallelContext::barrier();
  Logger::timed << "Cleaning up ccp files..." << std::endl;
  for (const auto &family : families) {
    std::remove(family.ccpFile.c_str());
  }
  ParallelContext::barrier();
}

/**
 *  Exclude families that:
 *  - cover too few species
 *  - are too large (optional)
 *  - have too much uncertainty (optional)
 */
void trimFamilies(const AleArguments &args, Families &families) {
  Logger::timed << "Families: " << families.size() << std::endl;
  auto minSpecies = args.minCoveredSpecies;
  auto trimRatio = args.trimFamilyRatio;
  auto maxCladeSplitRatio = args.maxCladeSplitRatio;
  if (minSpecies > 0) {
    Logger::timed << "Trimming families covering less than " << minSpecies
                  << " species..." << std::endl;
    TrimFamilies::trimMinSpeciesCoverage(families, minSpecies);
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
  if (trimRatio > 0.0) {
    Logger::timed << "Trimming (at least) " << trimRatio * 100.0
                  << "\% top families "
                  << "sorted by the number of observed clades..." << std::endl;
    TrimFamilies::trimHighCladesNumber(families, (1.0 - trimRatio));
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
  if (maxCladeSplitRatio >= 0.0) {
    Logger::timed
        << "Trimming families with a ratio (observed clades)/(gene nodes) > "
        << maxCladeSplitRatio << "..." << std::endl;
    TrimFamilies::trimCladeSplitRatio(families, maxCladeSplitRatio);
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
  checkEnoughFamilies(families);
  ParallelContext::barrier();
}

/**
 *  Compute an initial species tree from
 *  the given gene families and save it to a file
 */
void computeInitialSpeciesTree(const Families &families,
                               const SpeciesTreeAlgorithm &algo,
                               const std::string &outputFile) {
  switch (algo) {
  case SpeciesTreeAlgorithm::MiniNJ:
    MiniNJ::runMiniNJ(families)->save(outputFile);
    break;
  case SpeciesTreeAlgorithm::Random:
    std::make_unique<SpeciesTree>(families)->getTree().save(outputFile);
    break;
  case SpeciesTreeAlgorithm::User:
    assert(false);
    break;
  default:
    Logger::info << "\nError: The specified species tree building algorithm "
                 << "is not implemented in AleRax\n"
                 << std::endl;
    ParallelContext::abort(0);
  }
}

/**
 *  Generate or read the starting species tree and
 *  use it to initialize the current species tree
 *
 *  The tree can be either:
 *  - random
 *  - user-defined
 *  - computed with MiniNJ
 *
 *  @param args The program arguments
 *  @param families The gene family descriptions
 */
void initStartingSpeciesTree(AleArguments &args, const Families &families) {
  Logger::timed << "Initializing starting species tree..." << std::endl;
  // set the starting tree path
  auto startingSpeciesTree =
      Paths::getSpeciesTreeFile(args.output, "starting_species_tree.newick");
  // obtain the initial tree and save it as the starting tree
  if (args.speciesTreeAlgorithm == SpeciesTreeAlgorithm::User) {
    // ensure the initial tree file is readable
    if (ParallelContext::getRank() == 0) {
      try {
        SpeciesTree reader(args.speciesTree, true, false);
      } catch (const std::exception &e) {
        Logger::error << "Error while trying to parse the species tree:"
                      << std::endl;
        Logger::info << e.what() << std::endl;
        ParallelContext::abort(153);
      }
    }
    ParallelContext::barrier();
    // read the initial tree, provide it with internal node labels
    // and save it as the starting tree
    PLLRootedTree::labelRootedTree(args.speciesTree, startingSpeciesTree);
  } else {
    // compute the initial tree and save it as the starting tree
    computeInitialSpeciesTree(families, args.speciesTreeAlgorithm,
                              startingSpeciesTree);
  }
  ParallelContext::barrier();
  // set the current tree path for the further steps
  args.speciesTree =
      Paths::getSpeciesTreeFile(args.output, "inferred_species_tree.newick");
  // read the starting tree and save it as the current tree
  if (ParallelContext::getRank() == 0) {
    bool isDated = (args.transferConstraint == TransferConstaint::RELDATED);
    SpeciesTree currentTree(startingSpeciesTree, true, isDated);
    currentTree.getTree().save(args.speciesTree);
  }
  ParallelContext::barrier();
}

/**
 *  Initialize the current species tree from the checkpoint
 */
void initSpeciesTreeFromCheckpoint(AleArguments &args,
                                   const std::string &ckpDir) {
  Logger::timed << "Initializing species tree from the checkpoint..."
                << std::endl;
  // set the current tree path for the further steps
  args.speciesTree =
      Paths::getSpeciesTreeFile(args.output, "inferred_species_tree.newick");
  // read the checkpointed tree and save it as the current tree
  auto newick = AleState::readCheckpointSpeciesTree(ckpDir);
  ParallelOfstream os(args.speciesTree, true);
  os << newick << std::endl;
  os.close();
  ParallelContext::barrier();
}

/**
 *  Check that a family conditional clade probability object
 *  is valid: check that the mapping gene <-> species is bijective
 */
bool checkCCP(const FamilyInfo &family, const ConditionalClades &ccp,
              const std::unordered_set<std::string> &allSpecies) {
  GeneSpeciesMapping mapping;
  mapping.fill(family.mappingFile, family.startingGeneTree);
  auto m = mapping.getMap();
  for (const auto &p : ccp.getCidToLeaves()) {
    auto gene = p.second;
    auto it = m.find(gene);
    if (it == m.end()) {
      Logger::error << "Error: gene " << gene << " from family " << family.name
                    << " is not mapped to any species" << std::endl;
      return false;
    }
    auto species = it->second;
    if (allSpecies.find(species) == allSpecies.end()) {
      Logger::error << "Error: gene " << gene << " from family " << family.name
                    << " is mapped to species " << species
                    << ", but this species is not in the species tree"
                    << std::endl;
      return false;
    }
  }
  return true;
}

/**
 *  Run checkCCP on all gene families
 */
void checkCCPAndSpeciesTree(const Families &families,
                            const std::string &speciesTreePath) {
  Logger::timed << "Checking that ccps and mappings are valid..." << std::endl;
  PLLRootedTree speciesTree(speciesTreePath);
  auto labels = speciesTree.getLabels(true);
  auto N = families.size();
  bool ok = true;
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N);
       ++i) {
    ConditionalClades ccp;
    ccp.unserialize(families[i].ccpFile);
    ok &= checkCCP(families[i], ccp, labels);
  }
  ParallelContext::barrier();
  ParallelContext::parallelAnd(ok);
  if (!ok) {
    Logger::info << "\nProblems with CCPs, aborting." << std::endl;
    Logger::info << "Please check the files mentioned above\n" << std::endl;
    ParallelContext::abort(0);
  }
}

RecModelInfo buildRecModelInfo(const AleArguments &args) {
  return RecModelInfo(
      ArgumentsHelper::strToRecModel(args.reconciliationModelStr), args.recOpt,
      (args.modelParametrization ==
       ModelParametrization::PER_FAMILY), // per family rates
      args.gammaCategories, args.originationStrategy, args.pruneSpeciesTree,
      false, // rooted gene tree (option specific to GeneRax)
      false, // force gene tree root (option specific to GeneRax)
      false, // mad rooting (option specific to GeneRax)
      -1.0,  // branch length threshold
      args.transferConstraint, args.noDup, args.noDL, args.noTL,
      args.fractionMissingFile, args.memorySavings);
}

Parameters buildStartingRates(const AleArguments &args,
                              const RecModelInfo &info) {
  Parameters res(info.modelFreeParameters());
  double average = 0.0;
  switch (info.model) {
  case RecModel::UndatedDL:
    res[0] = args.d;
    res[1] = args.l;
    average = (args.d + args.l) / 2.0;
    break;
  case RecModel::UndatedDTL:
    res[0] = args.d;
    res[1] = args.l;
    res[2] = args.t;
    average = (args.d + args.l + args.t) / 3.0;
    break;
  default:
    assert(false);
  }
  if (args.originationStrategy == OriginationStrategy::OPTIMIZE) {
    // The initial value of the origination probabilities doesn't matter,
    // because they will be normalized during the likelihood computation.
    // We set it to the average of the other parameters such that all parameters
    // are in the same range (to facilitate gradient optimization)
    res[res.dimensions() - 1] = average;
  }
  return res;
}

void initCheckpoint(const AleArguments &args, const Families &families,
                    const std::string &ckpDir) {
  FileSystem::mkdir(ckpDir, true);
  ParallelContext::barrier();
  AleState::writeCheckpointCmd(args.getCommand(), ckpDir);
  AleState::writeCheckpointFamilies(families, ckpDir);
}

void runSpeciesTreeSearch(const AleArguments &args,
                          AleOptimizer &speciesTreeOptimizer) {
  Logger::info << std::endl;
  Logger::timed << "Start the species tree optimization..." << std::endl;
  if (args.randomSpeciesRoot) {
    Logger::info << std::endl;
    Logger::timed << "Random root position!" << std::endl;
    speciesTreeOptimizer.randomizeRoot();
  }
  // Species tree search or root search
  switch (args.speciesSearchStrategy) {
  case SpeciesSearchStrategy::HYBRID:
    speciesTreeOptimizer.optimize();
    break;
  case SpeciesSearchStrategy::REROOT:
    speciesTreeOptimizer.reroot();
    break;
  case SpeciesSearchStrategy::SKIP:
    Logger::info << std::endl;
    Logger::timed << "Optimization skipped!" << std::endl;
    break;
  default:
    assert(false); // not implemented yet
    break;
  }
  Logger::info << std::endl;
  Logger::timed << "End of the species tree optimization" << std::endl;
  speciesTreeOptimizer.saveSpeciesTree();
  Logger::timed << "Final species tree topology: " << args.speciesTree
                << std::endl;
}

void runFirstFinalRateOptimization(AleOptimizer &speciesTreeOptimizer) {
  Logger::info << std::endl;
  Logger::timed << "Final model rate optimization, non-thorough..."
                << std::endl;
  speciesTreeOptimizer.optimizeModelRates(false);
}

void runDateOptimization(const AleArguments &args,
                         AleOptimizer &speciesTreeOptimizer) {
  if (!args.inferSpeciationOrders) {
    return;
  }
  Logger::info << std::endl;
  Logger::timed << "Speciation order inference..." << std::endl;
  speciesTreeOptimizer.optimizeDates(true);
}

void runTransferHighwayInference(const AleArguments &args,
                                 AleOptimizer &speciesTreeOptimizer) {
  if (!args.highways) {
    return;
  }
  Logger::info << std::endl;
  Logger::timed << "Transfer highway inference..." << std::endl;
  speciesTreeOptimizer.inferHighways(args.highwayCandidateFile,
                                     args.highwayCandidatesStep1,
                                     args.highwayCandidatesStep2);
}

void runSecondFinalRateOptimization(const AleArguments &args,
                                    AleOptimizer &speciesTreeOptimizer) {
  if (!args.skipThoroughRates) {
    Logger::info << std::endl;
    Logger::timed << "Final model rate optimization, thorough..." << std::endl;
    speciesTreeOptimizer.optimizeModelRates(true);
  }
  speciesTreeOptimizer.saveRatesAndLL();
  auto finalLL = speciesTreeOptimizer.getEvaluator().computeLikelihood();
  Logger::timed << "Final species tree likelihood: ll=" << finalLL << std::endl;
}

void runReconciliationInference(const AleArguments &args,
                                AleOptimizer &speciesTreeOptimizer) {
  if (!args.geneTreeSamples) {
    return;
  }
  Logger::info << std::endl;
  Logger::timed << "Reconciling gene trees with the species tree..."
                << std::endl;
  speciesTreeOptimizer.reconcile(args.geneTreeSamples);
}

/**
 *  Resolve every --wgd declaration to a species-tree node index and broadcast
 *  it to the evaluator. A single label selects that terminal branch; two labels
 *  select the branch leading to their LCA.
 *
 *  NOTE: node indices are tied to the current species tree topology, so a WGD
 *  declared here is only meaningful while the topology is fixed. --wgd is
 *  therefore intended for use with --species-tree-search SKIP.
 *
 *  On a resumed run, the WGD retention(s) and per-event LORe resolution(s)
 *  have already been restored from the checkpoint (by the AleEvaluator
 *  constructor), reflecting whatever was fitted before the run was
 *  interrupted. We must NOT re-declare them here from the command-line
 *  starting values (args.wgds[*].q0, r=0.9), or we would silently discard
 *  that progress. checkCheckpointCmd already guarantees --wgd/--lore/
 *  --lore-wgd are unchanged between the original and resumed runs, so
 *  re-registering the (purely config, non-fitted) --lore-wgd targets and
 *  rebuilding the WGD subtree structure from them is safe to redo.
 */
// Resolve a WGD declaration (1 or 2 taxon labels) to a species-tree node. A
// single label selects that terminal branch; two labels select the branch
// leading to their LCA. `flag` names the option for error messages.
static corax_rnode_t *
resolveWGDNode(PLLRootedTree &speciesTree,
               const std::unordered_map<std::string, corax_rnode_t *> &labelToNode,
               const WGDDeclaration &wgd, const std::string &flag) {
  if (wgd.labels.size() == 1) {
    auto it = labelToNode.find(wgd.labels[0]);
    if (it == labelToNode.end()) {
      Logger::error << "Error: " << flag << " label '" << wgd.labels[0]
                    << "' was not found in the species tree" << std::endl;
      ParallelContext::abort(1);
    }
    return it->second;
  } else if (wgd.labels.size() == 2) {
    auto it1 = labelToNode.find(wgd.labels[0]);
    auto it2 = labelToNode.find(wgd.labels[1]);
    if (it1 == labelToNode.end() || it2 == labelToNode.end()) {
      Logger::error << "Error: " << flag << " labels '" << wgd.labels[0]
                    << "' and '" << wgd.labels[1]
                    << "' were not both found in the species tree" << std::endl;
      ParallelContext::abort(1);
    }
    return speciesTree.getLCA(it1->second, it2->second);
  }
  Logger::error << "Error: " << flag << " expects 1 or 2 taxon labels"
                << std::endl;
  ParallelContext::abort(1);
  return nullptr;
}

void declareWGDs(const AleArguments &args, AleOptimizer &optimizer,
                 bool checkpointDetected) {
  if (args.wgds.empty()) {
    return;
  }
  if (checkpointDetected) {
    auto &speciesTree = optimizer.getSpeciesTree().getTree();
    auto labelToNode = speciesTree.getLabelToNode(false);
    auto &evaluator = optimizer.getEvaluator();
    // --lore-wgd targets are pure config (derived from args, not fitted), so
    // it's safe -- and necessary -- to re-register them before rebuilding the
    // WGD subtree structure below.
    for (const auto &wgd : args.loreWgds) {
      auto *node = resolveWGDNode(speciesTree, labelToNode, wgd, "--lore-wgd");
      evaluator.addResolutionTarget(node->node_index);
    }
    if (!evaluator.getWGDNodes().empty()) {
      evaluator.buildWGDStructure();
      if (evaluator.optimizesResolution()) {
        // re-apply the fitted per-event r (already restored into
        // _wgdResolution by the AleEvaluator constructor) to every
        // evaluation's per-branch resolution vector, which needed the
        // subtree structure we just (re)built above
        std::vector<double> restoredR;
        for (auto node : evaluator.getWGDNodes()) {
          restoredR.push_back(evaluator.getWGDResolution(node));
        }
        evaluator.setWGDResolutions(restoredR);
      }
    }
    for (auto node : evaluator.getWGDNodes()) {
      Logger::timed << "Restored WGD from checkpoint at node " << node
                    << ", retention q=" << evaluator.getWGDRetention(node)
                    << std::endl;
    }
    if (evaluator.optimizesResolution()) {
      for (auto node : evaluator.getWGDNodes()) {
        Logger::timed << "Restored LORe resolution from checkpoint at node "
                      << node << ", r=" << evaluator.getWGDResolution(node)
                      << std::endl;
      }
    }
    return;
  }
  if (args.speciesSearchStrategy != SpeciesSearchStrategy::SKIP &&
      args.speciesSearchStrategy != SpeciesSearchStrategy::EVAL) {
    Logger::info << "Warning: --wgd is being used together with a species tree "
                 << "search; the WGD branch is resolved on the starting tree "
                 << "and may no longer be correct if the topology changes. "
                 << "Use --species-tree-search SKIP for WGD analyses."
                 << std::endl;
  }
  auto &speciesTree = optimizer.getSpeciesTree().getTree();
  auto labelToNode = speciesTree.getLabelToNode(false);
  auto &evaluator = optimizer.getEvaluator();
  std::set<unsigned int> declaredNodes;
  for (const auto &wgd : args.wgds) {
    auto *node = resolveWGDNode(speciesTree, labelToNode, wgd, "--wgd");
    evaluator.setWGD(node->node_index, wgd.q0);
    declaredNodes.insert(node->node_index);
    Logger::timed << "Declared a WGD at the top of the branch leading to node "
                  << node->node_index << " (label: "
                  << (node->label ? node->label : "<internal>")
                  << "), starting retention q=" << wgd.q0 << std::endl;
  }
  // --lore-wgd: fit r only for the named WGD(s); register them as resolution
  // targets so every other declared WGD is pinned to AORe (r=1).
  for (const auto &wgd : args.loreWgds) {
    auto *node = resolveWGDNode(speciesTree, labelToNode, wgd, "--lore-wgd");
    if (declaredNodes.find(node->node_index) == declaredNodes.end()) {
      Logger::error << "Error: --lore-wgd node " << node->node_index
                    << " (label: " << (node->label ? node->label : "<internal>")
                    << ") is not a declared WGD; add a matching --wgd"
                    << std::endl;
      ParallelContext::abort(1);
    }
    evaluator.addResolutionTarget(node->node_index);
    Logger::timed << "LORe target: will fit the resolution r of the WGD on the "
                  << "branch leading to node " << node->node_index << " (label: "
                  << (node->label ? node->label : "<internal>") << ")"
                  << std::endl;
  }
  // LORe: enable joint estimation of a per-event resolution probability r (start
  // near-AORe at r=0.9). With bare --lore every internal-branch WGD is fitted;
  // with --lore-wgd only the registered targets are. r=1 is the WHALE/AORe
  // limit; terminal-branch WGDs (and, under --lore-wgd, non-targets) keep r=1.
  if (args.lore || !args.loreWgds.empty()) {
    evaluator.setResolutionProb(0.9);
    evaluator.enableResolutionOptimization(true);
    if (args.loreWgds.empty()) {
      Logger::timed << "LORe enabled: will estimate a per-event resolution "
                    << "probability r per declared WGD (starting r=0.9)"
                    << std::endl;
    } else {
      Logger::timed << "LORe enabled (--lore-wgd): will estimate r for "
                    << args.loreWgds.size()
                    << " targeted WGD(s); all others pinned to AORe r=1 "
                    << "(starting r=0.9)" << std::endl;
    }
  }
}

/**
 *  Main function of AleRax once the arguments have been parsed
 *
 *  @param args The program arguments
 */
void run(AleArguments &args) {
  Random::setSeed(static_cast<unsigned int>(args.seed));
  // initializing the output dir, ccp files, species tree, model info, etc.
  bool checkpointDetected = AleOptimizer::checkpointExists(args.output);
  auto ckpDir = AleOptimizer::getCheckpointDir(args.output);
  if (checkpointDetected) {
    Logger::info << "Checkpoint detected" << std::endl;
    AleState::checkCheckpointCmd(args.getCommand(), ckpDir);
  }
  auto ccpDir = FileSystem::joinPaths(args.output, "ccps");
  if (!checkpointDetected) {
    initAleRaxDirectories(args, ccpDir);
  }
  printInitialMessage(args);
  auto ccpDimensionFile = FileSystem::joinPaths(ccpDir, "ccpdim.txt");
  auto families = FamiliesFileParser::parseFamiliesFile(args.families);
  for (auto &family : families) {
    family.ccpFile = FileSystem::joinPaths(ccpDir, family.name + ".ccp");
  }
  if (!checkpointDetected) {
    filterInvalidFamilies(args, families);
    generateCCPs(args, families, ccpDimensionFile);
    trimFamilies(args, families);
    initStartingSpeciesTree(args, families);
  } else {
    AleState::filterCheckpointFamilies(families, ckpDir);
    initSpeciesTreeFromCheckpoint(args, ckpDir);
  }
  checkCCPAndSpeciesTree(families, args.speciesTree);
  auto coverageFile =
      FileSystem::joinPaths(args.output, "perSpeciesCoverage.txt");
  auto fractionMissingFile =
      FileSystem::joinPaths(args.output, "perSpeciesMissing.txt");
  Family::printStats(families, args.speciesTree, coverageFile,
                     fractionMissingFile);
  // initializing the optimizer
  auto info = buildRecModelInfo(args);
  auto startingRates = buildStartingRates(args, info);
  AleOptimizer speciesTreeOptimizer(
      args.speciesTree, families, info, args.modelParametrization,
      args.optimizationClassFile, startingRates, !args.fixRates,
      args.optVerbose, args.output);
  speciesTreeOptimizer.getEvaluator().setSummaryOnly(args.summaryOnly);
  // declare any command-line WGDs on the (now constructed) species tree, or
  // restore those already recovered from the checkpoint on a resumed run
  declareWGDs(args, speciesTreeOptimizer, checkpointDetected);
  if (!checkpointDetected) {
    initCheckpoint(args, families, ckpDir);
    speciesTreeOptimizer.setCurrentStep(AleStep::SpeciesTreeOpt);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // species tree search
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::SpeciesTreeOpt) {
    runSpeciesTreeSearch(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::ModelRateOpt1);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // model parameter optimization on the final species tree topology
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::ModelRateOpt1) {
    runFirstFinalRateOptimization(speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::RelDating);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // relative species dating on the final species tree topology
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::RelDating) {
    runDateOptimization(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::Highways);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // highway search on the final species tree topology
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::Highways) {
    runTransferHighwayInference(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::ModelRateOpt2);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // final model parameter optimization
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::ModelRateOpt2) {
    runSecondFinalRateOptimization(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::Reconciliation);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // sampling reconciled gene trees
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::Reconciliation) {
    runReconciliationInference(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::End);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // end of the run
  Logger::info << std::endl;
  if (args.cleanupCCP || args.summaryOnly) {
    cleanupCCPs(families);
  }
  args.printWarning();
  Logger::timed << "End of AleRax execution" << std::endl;
}

/**
 *  This main will be called by each process and is the entry
 *  point to AleRax. It inits the parallel context, the logger,
 *  reads the arguments, runs AleRax and closes everything.
 *
 *  If comm is set to nullptr, MPI_Init will be called to initiate the
 *  MPI context
 */
int alerax_main(int argc, char **argv, void *comm) {
  ParallelContext::init(comm);
  Logger::init();
  Logger::timed << version << std::endl;
  AleArguments args(argc, argv);
  args.checkValid();
  run(args);
  Logger::close();
  ParallelContext::finalize();
  return 0;
}

int main(int argc, char **argv) {
#ifdef WITH_MPI
  // the null communicator signals that we need to call MPI_Init
  return alerax_main(argc, argv, nullptr);
#else
  // this fake communicator signals that we're not using MPI
  int noMPIComm = -1;
  return alerax_main(argc, argv, &noMPIComm);
#endif
}
