#pragma once

#include <string>

struct TreeBuildResult;

namespace mlipper {
namespace treeio {

std::string format_newick_taxon_name(const std::string& name);
// Reject cycles, shared/disconnected nodes, inconsistent parent links, and
// invalid branch lengths before any recursive serializer walks the topology.
void validate_tree_for_output(const TreeBuildResult& tree);
std::string write_tree_to_newick_string(const TreeBuildResult& tree);

void write_tree_to_newick_file(
    const TreeBuildResult& tree,
    const std::string& path,
    double collapse_internal_epsilon = -1.0);

} // namespace treeio

} // namespace mlipper
