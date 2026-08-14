#pragma once

#include <cstddef>
#include <string>

#include "../tree/tree.hpp"

namespace mlipper {
namespace treeio {

std::string write_tree_to_newick_string(const TreeBuildResult& tree);

size_t write_tree_to_newick_file(
    const TreeBuildResult& tree,
    const std::string& path,
    double collapse_internal_epsilon = -1.0);

void print_tree_structure(const TreeBuildResult& tree);
} // namespace treeio

} // namespace mlipper
