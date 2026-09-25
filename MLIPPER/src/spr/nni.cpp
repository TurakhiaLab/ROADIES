#include "nni.hpp"

#include "local_spr_internal.hpp"

namespace mlipper::nni {

std::vector<int> candidate_regraft_edges(
    const TreeBuildResult& pruned_tree,
    const PruneInfo& prune_info)
{
    const int central_parent = prune_info.grandparent_id;
    if (central_parent < 0 ||
        central_parent == pruned_tree.root_id ||
        central_parent >= static_cast<int>(pruned_tree.nodes.size())) {
        return {};
    }

    const TreeNode& central =
        pruned_tree.nodes[static_cast<size_t>(central_parent)];
    const int uncle_id = central.left == prune_info.sibling_id
        ? central.right
        : (central.right == prune_info.sibling_id ? central.left : -1);
    if (uncle_id < 0) {
        return {};
    }

    // The sibling edge reconstructs the baseline topology and supplies its
    // likelihood. The uncle edge is the NNI alternative for this prune root;
    // pruning the central node's other child evaluates the second alternative.
    return {prune_info.sibling_id, uncle_id};
}

} // namespace mlipper::nni
