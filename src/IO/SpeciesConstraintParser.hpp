#pragma once

#include <string>
#include <vector>

class PLLRootedTree;

class SpeciesConstraintParser {
public:
  SpeciesConstraintParser() = delete;
  static std::vector<std::string> parse(const std::string &constraintFile,
                                        const PLLRootedTree &speciesTree);
};
