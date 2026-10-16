#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>

/**
 *  Synteny observations for one gene family (kalerax Layer 1).
 *  For same-genome gene pairs: y=1 (anchor pair in a collinear duplicated
 *  block), y=0 (both genes placed, not collinear) or NA (undetermined).
 */
struct SyntenyInput {
  // ordered-pair key (see syntenyKey) -> y in {0,1}
  std::unordered_map<std::string, int> pairObs;
  std::unordered_set<std::string> naGenes; // all pairs involving these: NA
  std::unordered_set<std::string> naPairs; // explicit NA pairs
  bool unlistedNA = false; // unlisted same-species pairs: NA (else y = 0)
};

inline std::string syntenyKey(const std::string &a, const std::string &b) {
  return (a < b) ? (a + '\t' + b) : (b + '\t' + a);
}
