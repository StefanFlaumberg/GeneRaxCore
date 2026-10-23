#include "SpeciesConstraintParser.hpp"

#include <cassert>
#include <fstream>
#include <sstream>

#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <trees/PLLRootedTree.hpp>

std::vector<std::string>
SpeciesConstraintParser::parse(const std::string &constraintFile,
                               const PLLRootedTree &speciesTree) {
  std::vector<std::string> speciesToFixedCladeLabel;
  speciesToFixedCladeLabel.resize(speciesTree.getNodeNumber());
  auto labelToNode = speciesTree.getLabelToNode(false);
  std::ifstream is(constraintFile);
  if (!is || is.peek() == std::ifstream::traits_type::eof()) {
    Logger::error << "Error: cannot open species tree constraint file "
                  << constraintFile << std::endl;
    ParallelContext::abort(1);
  }
  std::string line;
  unsigned int lineNumber = 0;
  while (std::getline(is, line)) {
    lineNumber++;
    std::istringstream iss(line);
    std::string speciesLabel;
    std::string bufferStr;
    if (!(iss >> speciesLabel)) {
      Logger::error << "Error: cannot read species label at line " << lineNumber
                    << " in file " << constraintFile << std::endl;
      ParallelContext::abort(1);
    }
    auto it = labelToNode.find(speciesLabel);
    if (it == labelToNode.end()) {
      Logger::error << "Error: label " << speciesLabel << " in file "
                    << constraintFile << " is not in the species tree"
                    << std::endl;
      ParallelContext::abort(1);
    } else {
      auto node = it->second;
      if (node->left) { // leaves can be ignored
        speciesToFixedCladeLabel[node->node_index] = speciesLabel;
      }
    }
    if (iss >> bufferStr) {
      Logger::error << "Error: line " << lineNumber << " in file "
                    << constraintFile << " contains an extra field"
                    << std::endl;
      ParallelContext::abort(1);
    }
  }
  is.close();
  return speciesToFixedCladeLabel;
}
