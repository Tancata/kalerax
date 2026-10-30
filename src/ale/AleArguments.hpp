#pragma once

#include <string>
#include <vector>

#include <IO/ArgumentsHelper.hpp>

/**
 *  A whole-genome duplication declared on the command line.
 *  The WGD is placed at the top of the branch leading to the node identified by
 *  `labels`: a single taxon label selects that terminal branch; two taxon
 *  labels select the branch leading to their LCA (mirrors WHALE's
 *  insertnode-at-LCA). `q0` is the starting retention probability.
 */
struct WGDDeclaration {
  std::vector<std::string> labels; // 1 or 2 taxon labels
  double q0;                       // starting retention probability
};

/**
 *  Parses and stores the program arguments
 */
class AleArguments {
public:
  /**
   *  Parse the arguments from main()
   */
  AleArguments(int iargc, char **iargv);

  /**
   *  Print AleRax' help message
   */
  void printHelp() const;

  /**
   *  Return the command line used to call the AleRax' executable
   */
  const std::string getCommand() const;

  /**
   *  Print the command line used to call the AleRax' executable
   */
  void printCommand() const;

  /**
   *  Print a user-friendly summary of the most important
   *  parameters set by the user
   */
  void printSummary() const;

  /**
   *  Print warnings if some of the arguments used are valid, but
   *  still might result in undesired behaviour
   */
  void printWarning() const;

  /**
   *  Check that the arguments are compatible with each other.
   *  If not, terminates the program with an explicit error message
   */
  void checkValid() const;

public:
  // the parameters of the main() function
  int argc;
  char **argv;

  // input data
  std::string families;
  std::string speciesTree;

  // whole-genome duplications declared via --wgd
  std::vector<WGDDeclaration> wgds;
  // LORe: estimate a (global) delayed-rediploidization resolution prob r
  bool lore;
  // LORe restricted to specific WGD(s): each entry names (by 1 or 2 taxon
  // labels, like --wgd) a declared WGD whose resolution r is fitted; every
  // other declared WGD is pinned to AORe (r=1). Unlike bare --lore this allows
  // the LORe-tested WGD to be nested within (AORe) WGDs. q0 is ignored here.
  std::vector<WGDDeclaration> loreWgds;

  // model
  std::string reconciliationModelStr;
  TransferConstaint transferConstraint;
  bool noDup;
  bool noDL;
  bool noTL;
  bool pruneSpeciesTree;
  std::string fractionMissingFile;
  CCPRooting ccpRooting;
  OriginationStrategy originationStrategy;
  bool memorySavings;
  double d;
  double l;
  double t;

  // search
  SpeciesTreeAlgorithm speciesTreeAlgorithm;
  SpeciesSearchStrategy speciesSearchStrategy;
  bool inferSpeciationOrders;
  ModelParametrization modelParametrization;
  std::string optimizationClassFile;
  unsigned int gammaCategories;
  RecOpt recOpt;
  bool fixRates;
  // LORe: hold the per-event resolution r fixed at this value on every
  // resolvable WGD (--lore-fixed-r; <0 = fit r as usual)
  double loreFixedR;
  // hold every WGD retention q fixed at its --wgd LABEL:q0 value
  bool fixWgdQ;
  double gammaAlpha;     // --speciation-gamma-alpha (<=0: default start 1.0)
  bool fixGammaAlpha;    // --fix-gamma-alpha
  unsigned int convergeRounds;  // --converge-rounds (max thorough passes)
  bool skipThoroughRates;

  // highways
  bool highways;
  std::string highwayCandidateFile;
  unsigned int highwayCandidatesStep1;
  unsigned int highwayCandidatesStep2;

  // trimming
  bool skipFamilyFiltering;
  unsigned int minCoveredSpecies;
  double trimFamilyRatio;
  double maxCladeSplitRatio;
  unsigned int sampleFrequency;

  // output
  std::string output;
  unsigned int geneTreeSamples;
  bool cleanupCCP;
  // suppress the bulky per-family/per-sample output (reconciliations/all/,
  // RecPhyloXML, highway candidate tests, CCP binaries); keep the summaries,
  // global totals, wgdSummary, species trees, model parameters and logs.
  bool summaryOnly;
  bool orderedDupWeights; // deprecated no-op (forward-process weights are the default)
  bool aleDupWeights;     // --ale-dup-weights: restore ALE/AleRax D vs T weighting

  // synteny (Layer 1): per-family <dir>/<family>.syn pair observations
  std::string syntenyDir;
  double sigmaW0;
  double sigmaN0;
  bool fixSynteny;
  double sigmaWMin; // --synteny-min-sigmaW (0 = unconstrained)
  bool syntenyUnlistedNA;

  // random seed
  int seed;

  // experimental
  bool randomSpeciesRoot;
  bool optVerbose;
};
