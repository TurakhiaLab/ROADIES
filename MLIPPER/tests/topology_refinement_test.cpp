#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "spr/local_spr_internal.hpp"
#include "spr/nni.hpp"
#include "tree/divide_and_conquer.hpp"
#include "tree/tree_topology_utils.hpp"

namespace {

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

TreeBuildResult make_test_tree()
{
    TreeBuildResult tree;
    tree.nodes.resize(9);
    tree.root_id = 8;
    for (int node_id = 0; node_id < 9; ++node_id) {
        tree.nodes[static_cast<size_t>(node_id)].id = node_id;
        tree.nodes[static_cast<size_t>(node_id)].stable_node_label = node_id;
        tree.nodes[static_cast<size_t>(node_id)].branch_length_to_parent =
            fp_t(0.1);
    }

    for (int tip_id = 0; tip_id < 5; ++tip_id) {
        TreeNode& tip = tree.nodes[static_cast<size_t>(tip_id)];
        tip.is_tip = true;
        tip.name = std::string(1, static_cast<char>('A' + tip_id));
        tree.tip_node_by_name.emplace(tip.name, tip_id);
    }

    tree.nodes[0].parent = 5;
    tree.nodes[1].parent = 5;
    tree.nodes[5].parent = 6;
    tree.nodes[5].left = 0;
    tree.nodes[5].right = 1;

    tree.nodes[2].parent = 6;
    tree.nodes[6].parent = 8;
    tree.nodes[6].left = 5;
    tree.nodes[6].right = 2;

    tree.nodes[3].parent = 7;
    tree.nodes[4].parent = 7;
    tree.nodes[7].parent = 8;
    tree.nodes[7].left = 3;
    tree.nodes[7].right = 4;

    tree.nodes[8].parent = -1;
    tree.nodes[8].left = 6;
    tree.nodes[8].right = 7;
    tree.nodes[8].branch_length_to_parent = fp_t(0);
    rebuild_traversals(tree);
    return tree;
}

void test_anchor_parent_refresh()
{
    TreeBuildResult tree = make_test_tree();
    const std::vector<LocalSPRInsertionAnchor> anchors{{5, 0}};

    tree.nodes[5].left = 2;
    tree.nodes[2].parent = 5;
    tree.nodes[6].right = 0;
    tree.nodes[0].parent = 6;
    rebuild_traversals(tree);

    const auto refreshed = refresh_local_spr_anchor_parents(tree, anchors);
    require(refreshed.size() == 1, "anchor refresh changed the anchor count");
    require(refreshed[0].endpoint_a == 6,
            "anchor refresh did not follow the tip's current parent");
    require(refreshed[0].endpoint_b == 0,
            "anchor refresh changed the stable tip id");
}

void test_nni_candidate_edges()
{
    TreeBuildResult tree = make_test_tree();
    PruneInfo prune_info;
    require(prune_subtree_for_spr(tree, 0, prune_info),
            "failed to prune the NNI test tip");
    require(prune_info.sibling_id == 1, "wrong NNI sibling");
    require(prune_info.grandparent_id == 6, "wrong NNI central edge");

    const std::vector<int> candidates =
        mlipper::nni::candidate_regraft_edges(tree, prune_info);
    require(candidates == std::vector<int>({1, 2}),
            "NNI regraft edges do not contain the baseline and alternative");
}

void test_nni_rejects_root_adjacent_and_malformed_edges()
{
    TreeBuildResult root_adjacent_tree = make_test_tree();
    PruneInfo root_adjacent_prune;
    require(prune_subtree_for_spr(
                root_adjacent_tree, 5, root_adjacent_prune),
            "failed to prepare root-adjacent NNI test");
    require(
        mlipper::nni::candidate_regraft_edges(
            root_adjacent_tree, root_adjacent_prune).empty(),
        "NNI accepted a central edge adjacent to the rooted representation");

    TreeBuildResult malformed_tree = make_test_tree();
    PruneInfo malformed_prune;
    require(prune_subtree_for_spr(malformed_tree, 0, malformed_prune),
            "failed to prepare malformed NNI test");
    malformed_prune.sibling_id = 4;
    require(
        mlipper::nni::candidate_regraft_edges(
            malformed_tree, malformed_prune).empty(),
        "NNI accepted a sibling outside the central node");
}

void test_spr_prune_and_regraft()
{
    TreeBuildResult tree = make_test_tree();
    PruneInfo prune_info;
    require(prune_subtree_for_spr(tree, 0, prune_info),
            "failed to prune the SPR test tip");

    regraft_subtree_for_spr(
        tree,
        prune_info,
        2,
        0.2,
        0.04,
        OPT_BRANCH_LEN_MIN,
        OPT_BRANCH_LEN_MIN);
    require(tree.nodes[0].parent == 5,
            "regrafted tip does not use the recycled internal node");
    require(tree.nodes[2].parent == 5,
            "regraft target is not attached below the recycled internal node");
    require(tree.nodes[5].parent == 6,
            "recycled internal node is attached to the wrong parent");
    require(tree.nodes[6].right == 5,
            "target edge was not replaced by the recycled internal node");
    require(tree.postorder.size() == tree.nodes.size(),
            "SPR regraft did not rebuild the postorder traversal");
    require(tree.preorder.size() == tree.nodes.size(),
            "SPR regraft did not rebuild the preorder traversal");
    require(std::isfinite(static_cast<double>(
                tree.nodes[0].branch_length_to_parent)) &&
            tree.nodes[0].branch_length_to_parent > fp_t(0),
            "SPR regraft produced an invalid pendant branch length");
}

void test_spr_candidate_legality()
{
    const TreeBuildResult tree = make_test_tree();
    LocalSPRRepairUnit unit;
    unit.envelope_mask.assign(tree.nodes.size(), 1);

    LocalSPRCandidateMove candidate;
    candidate.prune_root_id = 0;
    candidate.old_parent_id = 5;
    candidate.regraft_child_id = 2;
    candidate.regraft_parent_id = 6;
    require(local_spr_candidate_still_legal(tree, unit, candidate),
            "valid SPR candidate was rejected");

    candidate.old_parent_id = 6;
    require(!local_spr_candidate_still_legal(tree, unit, candidate),
            "candidate with a stale prune parent was accepted");
    candidate.old_parent_id = 5;

    unit.envelope_mask[2] = 0;
    require(!local_spr_candidate_still_legal(tree, unit, candidate),
            "candidate outside the repair envelope was accepted");
}

void test_divide_and_conquer_subtree_plan()
{
    const TreeBuildResult tree = make_test_tree();
    const mlipper::divide_and_conquer::SubtreePlan plan =
        mlipper::divide_and_conquer::build_subtree_plan_from_anchors(
            tree,
            {{5, 6}},
            4);
    require(plan.selected_tip_ids.size() == 4,
            "D&C subtree planner did not honor the tip budget");
    require(plan.local_to_global == plan.subtree_node_ids,
            "D&C subtree node mapping is not one-to-one");

    const TreeBuildResult subtree =
        mlipper::divide_and_conquer::build_subtree_from_plan(tree, plan);
    require(subtree.nodes.size() == plan.local_to_global.size(),
            "D&C subtree and local-to-global mapping sizes differ");
    require(subtree.preorder.size() == subtree.nodes.size(),
            "D&C subtree traversal is incomplete");
}

void test_nni_sector_plan_contains_only_the_core_skeleton()
{
    const TreeBuildResult tree = make_test_tree();
    const mlipper::divide_and_conquer::SubtreePlan one_edge =
        mlipper::divide_and_conquer::build_nni_sector_plan(
            tree,
            {{5, 6}});
    require(one_edge.subtree_node_ids == std::vector<int>({5, 6}),
            "NNI sector retained nodes outside a single owned core edge");
    require(one_edge.selected_tip_ids.empty(),
            "NNI sector unexpectedly expanded toward real tips");
    require(one_edge.global_to_local.empty(),
            "compact NNI sector allocated a full-tree node index");
    const TreeBuildResult one_edge_subtree =
        mlipper::divide_and_conquer::build_subtree_from_plan(
            tree,
            one_edge);
    require(one_edge_subtree.nodes.size() == 2,
            "compact NNI subtree contains nodes outside its core edge");

    const mlipper::divide_and_conquer::SubtreePlan connected =
        mlipper::divide_and_conquer::build_nni_sector_plan(
            tree,
            {{5, 6}, {7, 8}});
    require(connected.subtree_node_ids == std::vector<int>({5, 6, 7, 8}),
            "NNI sector did not retain the minimal connector between core edges");
}

void test_nni_repair_unit_covers_the_compact_sector()
{
    const TreeBuildResult tree = make_test_tree();
    const auto units = build_nni_repair_units(tree, {5});
    require(units.size() == 1,
            "NNI central edges should share one deterministic repair unit");
    require(units[0].envelope_nodes.size() == tree.nodes.size(),
            "NNI repair unit does not cover the compact sector");
    require(units[0].envelope_mask.size() == tree.nodes.size(),
            "NNI repair mask has the wrong size");
}

void test_cuda_error_preserves_runtime_code()
{
    try {
        throw CudaRuntimeError(cudaErrorMemoryAllocation, "test OOM");
    } catch (const CudaRuntimeError& error) {
        require(error.code() == cudaErrorMemoryAllocation,
                "CUDA exception discarded its runtime error code");
        return;
    }
    throw std::runtime_error("typed CUDA exception was not caught");
}

void test_sequential_branch_traversal_order()
{
    const auto steps =
        build_sequential_branch_traversal_steps(make_test_tree());
    const std::vector<int> expected_nodes{
        6, 5, 0, 1, 5, 2, 6, 7, 3, 4, 7, 8};
    require(steps.size() == expected_nodes.size(),
            "sequential branch traversal has the wrong step count");
    for (size_t index = 0; index < steps.size(); ++index) {
        require(steps[index].node_id == expected_nodes[index],
                "sequential branch traversal changed DFS order");
    }
    require(
        steps[4].kind ==
            SequentialBranchTraversalStepKind::RefreshParentUpward &&
        steps[6].kind ==
            SequentialBranchTraversalStepKind::RefreshParentUpward &&
        steps[10].kind ==
            SequentialBranchTraversalStepKind::RefreshParentUpward &&
        steps[11].kind ==
            SequentialBranchTraversalStepKind::RefreshParentUpward,
        "sequential branch traversal misplaced an upward refresh");
}

} // namespace

int main()
{
    try {
        test_anchor_parent_refresh();
        test_nni_candidate_edges();
        test_nni_rejects_root_adjacent_and_malformed_edges();
        test_spr_prune_and_regraft();
        test_spr_candidate_legality();
        test_divide_and_conquer_subtree_plan();
        test_nni_sector_plan_contains_only_the_core_skeleton();
        test_nni_repair_unit_covers_the_compact_sector();
        test_cuda_error_preserves_runtime_code();
        test_sequential_branch_traversal_order();
        std::cout << "topology_refinement_test: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "topology_refinement_test: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
