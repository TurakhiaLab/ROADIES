#pragma once

#include <utility>
#include <vector>

#include "tree.hpp"

// Recompute deterministic left-before-right preorder and postorder lists after
// a topology edit. This helper assumes parent/child links already form a valid
// rooted tree; it does not detect cycles, shared children, or disconnected
// nodes. An invalid root produces empty traversals without changing node IDs.
inline void rebuild_traversals(TreeBuildResult& tree)
{
    tree.preorder.clear();
    tree.postorder.clear();
    if (tree.root_id < 0 ||
        tree.root_id >= static_cast<int>(tree.nodes.size())) {
        return;
    }

    std::vector<int> stack;
    stack.push_back(tree.root_id);
    while (!stack.empty()) {
        const int node_id = stack.back();
        stack.pop_back();
        tree.preorder.push_back(node_id);

        const TreeNode& node = tree.nodes[static_cast<size_t>(node_id)];
        if (!node.is_tip) {
            if (node.right >= 0) {
                stack.push_back(node.right);
            }
            if (node.left >= 0) {
                stack.push_back(node.left);
            }
        }
    }

    std::vector<std::pair<int, bool>> postorder_stack;
    postorder_stack.emplace_back(tree.root_id, false);
    while (!postorder_stack.empty()) {
        const auto [node_id, expanded] = postorder_stack.back();
        postorder_stack.pop_back();
        if (node_id < 0) {
            continue;
        }

        const TreeNode& node = tree.nodes[static_cast<size_t>(node_id)];
        if (expanded) {
            tree.postorder.push_back(node_id);
            continue;
        }

        postorder_stack.emplace_back(node_id, true);
        if (!node.is_tip) {
            if (node.right >= 0) {
                postorder_stack.emplace_back(node.right, false);
            }
            if (node.left >= 0) {
                postorder_stack.emplace_back(node.left, false);
            }
        }
    }
}

enum class SequentialBranchTraversalStepKind {
    OptimizeChildEdge,
    RefreshParentUpward,
};

struct SequentialBranchTraversalStep {
    SequentialBranchTraversalStepKind kind;
    int node_id = -1;
};

inline void append_sequential_branch_traversal_steps(
    const TreeBuildResult& tree,
    int parent_id,
    std::vector<SequentialBranchTraversalStep>& steps)
{
    // A child edge is optimized while its current outside message is valid.
    // Descendants follow, then the parent's upward message is refreshed so the
    // caller can continue without a full-tree CLV rebuild after every edge.
    const TreeNode& parent = tree.nodes.at(static_cast<size_t>(parent_id));
    if (parent.is_tip) return;
    const int children[2] = {parent.left, parent.right};
    for (int child_id : children) {
        if (child_id < 0) continue;
        steps.push_back({
            SequentialBranchTraversalStepKind::OptimizeChildEdge,
            child_id});
        if (!tree.nodes.at(static_cast<size_t>(child_id)).is_tip) {
            append_sequential_branch_traversal_steps(tree, child_id, steps);
        }
    }
    steps.push_back({
        SequentialBranchTraversalStepKind::RefreshParentUpward,
        parent_id});
}

inline std::vector<SequentialBranchTraversalStep>
build_sequential_branch_traversal_steps(const TreeBuildResult& tree)
{
    std::vector<SequentialBranchTraversalStep> steps;
    if (tree.root_id < 0 ||
        tree.root_id >= static_cast<int>(tree.nodes.size())) {
        return steps;
    }
    steps.reserve(tree.nodes.size() * 2);
    append_sequential_branch_traversal_steps(tree, tree.root_id, steps);
    return steps;
}
