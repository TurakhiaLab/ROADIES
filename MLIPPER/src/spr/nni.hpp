#pragma once

#include <vector>

#include "tree/tree.hpp"

struct PruneInfo;

namespace mlipper::nni {

// Return the two legal NNI alternatives created by the supplied prune. Edges
// use MLIPPER's child-endpoint representation; the unchanged attachment edge
// and edges inside the detached subtree are excluded.
std::vector<int> candidate_regraft_edges(
    const TreeBuildResult& pruned_tree,
    const PruneInfo& prune_info);

} // namespace mlipper::nni
