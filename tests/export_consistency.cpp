// Export-consistency test for the WGD/LORe-aware scenario sampler.
//
// Identity used (independent of the sampler): give the selected same-genome
// pairs y=1 and every other pair NA. Then a WGD split carrying n marked cross
// pairs is weighted by rho = (sigmaW/sigmaN)^n, so
//   d lnL / d ln(sigmaW) |_{sigmaW = sigmaN} = E_post[ sum over WGD splits of n ].
// The left side is computed by central finite differences of the likelihood;
// the right side is estimated from scenarios drawn by the exported sampler
// (MultiModel::sampleReconciliations), whose WGD splits are labelled "W" in the
// NewickEvents output written here. Agreement validates that the exported
// scenarios are posterior draws under the fitted q and r (and the synteny
// weights) -- the property that was violated before the export fix.
//
// usage: export_consistency <species.nw> <genetrees.nw> <mapping.link>
//          <taxonA> <taxonB> <q> <r> <DL|DTL> <samples> <synfile> <out.trees>
//          [seed] [sigmaW] [sigmaN] [dupSplitWeight (1 = default, 2 = --ordered-dup-weights)]

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>

#include <IO/GeneSpeciesMapping.hpp>
#include <IO/ParallelOfstream.hpp>
#include <ccp/ConditionalClades.hpp>
#include <trees/DatedTree.hpp>
#include <trees/PLLRootedTree.hpp>
#include <util/RecModelInfo.hpp>
#include <util/Scenario.hpp>
#include <util/enums.hpp>
#include <maths/Random.hpp>

#include "../src/ale/UndatedDLMultiModel.hpp"
#include "../src/ale/UndatedDTLMultiModel.hpp"

static RecModelInfo makeInfo(RecModel m) {
  return RecModelInfo(m, RecOpt::LBFGSB, false, 1, OriginationStrategy::UNIFORM,
                      false, false, false, false, -1.0, TransferConstaint::PARENTS,
                      false, false, false, "", false);
}

static SyntenyInput readSyn(const std::string &f) {
  SyntenyInput in;
  in.unlistedNA = true;
  std::ifstream is(f);
  std::string a, b, c;
  while (is >> a >> b >> c) {
    if (c == "1") in.pairObs[syntenyKey(a, b)] = 1;
    else if (c == "0") in.pairObs[syntenyKey(a, b)] = 0;
  }
  return in;
}

template <class MODEL>
static int run(MODEL &model, PLLRootedTree &speciesTree, const std::string &tA,
               const std::string &tB, double q, double r, unsigned int samples,
               const SyntenyInput &syn, const std::string &out, double sW,
               double sN, double dupW) {
  auto labelToNode = speciesTree.getLabelToNode(false);
  auto *w = speciesTree.getLCA(labelToNode[tA], labelToNode[tB]);
  model.setWGD(w->node_index, q);
  std::vector<double> rv(speciesTree.getNodeNumber(), 1.0);
  for (auto node : speciesTree.getNodes()) {
    if (node == w || speciesTree.isAncestorOf(w->node_index, node->node_index)) {
      rv[node->node_index] = r;
    }
  }
  model.setResolutionProbVector(rv);
  model.setDupSplitWeight(dupW); // --ordered-dup-weights (DTL only)
  model.setSynteny(syn);
  // with only y=1 pairs marked (others NA) the identity holds for any
  // sigmaW >= sigmaN; sigmaW == sigmaN is the plain (unweighted) model
  const double h = 1e-4;
  model.setSyntenyParams(sW * std::exp(h), sN);
  double lp = model.computeLogLikelihood();
  model.setSyntenyParams(sW * std::exp(-h), sN);
  double lm = model.computeLogLikelihood();
  model.setSyntenyParams(sW, sN); // sample at this point
  double l0 = model.computeLogLikelihood();
  double deriv = (lp - lm) / (2.0 * h);
  std::printf("sigmaW=%g sigmaN=%g lnL=%.10f  dlnL/dln(sigmaW)=%.6f\n", sW, sN, l0, deriv);
  std::vector<std::shared_ptr<Scenario>> scenarios;
  model.sampleReconciliations(samples, scenarios);
  {
    ParallelOfstream os(out, false);
    for (auto &sc : scenarios) {
      sc->saveReconciliation(os, ReconciliationFormat::NewickEvents);
    }
  }
  std::printf("sampled=%zu\n", scenarios.size());
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 12) {
    std::fprintf(stderr, "usage: see header\n");
    return 2;
  }
  Random::setSeed(argc > 12 ? atoi(argv[12]) : 42);
  double sN = argc > 14 ? atof(argv[14]) : 0.1;
  double sW = argc > 13 ? atof(argv[13]) : sN;
  double dupW = argc > 15 ? atof(argv[15]) : 1.0;
  if (sW < sN) { std::fprintf(stderr, "need sigmaW >= sigmaN\n"); return 2; }
  PLLRootedTree speciesTree(argv[1], true);
  GeneSpeciesMapping mapping;
  mapping.fill(argv[3], argv[2]);
  ConditionalClades ccp(argv[2], "", CCPRooting::UNIFORM);
  std::string ccpFile = std::string(argv[11]) + ".ccp";
  ccp.serialize(ccpFile);
  double q = atof(argv[6]), r = atof(argv[7]);
  unsigned int samples = atoi(argv[9]);
  auto syn = readSyn(argv[10]);
  std::string m = argv[8];
  int rc = 0;
  if (m == "DL") {
    UndatedDLMultiModel<double> model(speciesTree, mapping,
                                      makeInfo(RecModel::UndatedDL), ccpFile);
    rc = run(model, speciesTree, argv[4], argv[5], q, r, samples, syn, argv[11], sW, sN, dupW);
  } else {
    DatedTree dated(speciesTree, false);
    UndatedDTLMultiModel<double> model(dated, mapping,
                                       makeInfo(RecModel::UndatedDTL), ccpFile);
    rc = run(model, speciesTree, argv[4], argv[5], q, r, samples, syn, argv[11], sW, sN, dupW);
  }
  std::remove(ccpFile.c_str());
  return rc;
}
