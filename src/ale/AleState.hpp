#pragma once

#include <memory>

#include <IO/Families.hpp>
#include <IO/HighwayCandidateParser.hpp>
#include <trees/SpeciesTree.hpp>
#include <util/types.hpp>

#include "AleModelParameters.hpp"

/**
 *  Represents the current step in the AleRax pipeline,
 *  to restart from the right place after a checkpoint.
 *  The values should be ordered and contiguous
 */
enum class AleStep {
  Init = 0,
  SpeciesTreeOpt = 1,
  ModelRateOpt1 = 2,
  RelDating = 3,
  Highways = 4,
  ModelRateOpt2 = 5,
  Reconciliation = 6,
  End = 7
};

/**
 *  Stores all information about the current state of an AleRax run:
 *  - the current step in the pipeline
 *  - the current species tree
 *  - the current model parameters
 *
 *  Implements functions to handle checkpoints
 */
struct AleState {
  /**
   *  Constructor
   */
  AleState()
      : currentStep(AleStep::Init), mixtureAlpha(1.0),
        optimizeResolution(false) {}

  /**
   *  Dump the current run arguments to the checkpoint directory
   */
  static void writeCheckpointCmd(const std::string &currentCmd,
                                 const std::string &checkpointDir);

  /**
   *  Make sure the checkpoint used the same arguments as
   *  the current run
   */
  static void checkCheckpointCmd(const std::string &currentCmd,
                                 const std::string &checkpointDir);

  /**
   *  Dump the list of the accepted family names to the checkpoint
   *  directory
   */
  static void writeCheckpointFamilies(const Families &families,
                                      const std::string &checkpointDir);

  /**
   *  Retain only the families listed in the checkpoint
   */
  static void filterCheckpointFamilies(Families &families,
                                       const std::string &checkpointDir);

  /**
   *  Dump the current state to the checkpoint directory
   */
  void serialize(const std::string &checkpointDir) const;

  /**
   *  Load the current state from the checkpoint directory
   */
  void unserialize(const std::string &checkpointDir);

  /**
   *  Read the species tree newick from the checkpoint directory
   */
  static std::string
  readCheckpointSpeciesTree(const std::string &checkpointDir);

  // the running step
  AleStep currentStep;
  // the current species tree
  std::unique_ptr<SpeciesTree> speciesTree;
  // the current model parameters
  double mixtureAlpha;
  std::vector<Highway> transferHighways;
  std::vector<AleModelParameters> perLocalFamilyModelParams;
  // the names of the families, to map the model parameters to their
  // respective families
  std::vector<std::string> localFamilyNames;
  // WGDs: the species branch (by label, since node indices are not stable
  // across a species tree reload) carrying each declared WGD, its current
  // (possibly fitted) retention probability q, and its current (possibly
  // fitted) per-event LORe resolution probability r (1.0 == AORe / not
  // resolvable). All three vectors are aligned.
  std::vector<std::string> wgdBranchLabels;
  std::vector<double> wgdRetentions;
  std::vector<double> wgdResolutions;
  // LORe: whether the resolution probability(ies) are being jointly
  // optimized with the WGD retentions (mirrors --lore / --lore-wgd)
  bool optimizeResolution;
};
