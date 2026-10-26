#include "local_spr_internal.hpp"
#include "tree/tree_topology_utils.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include "mlipper_session.hpp"
#include "util/checked_size.hpp"
#include "util/mlipper_util.h"

void local_spr_assert(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(
            "Local SPR subtree assertion failed: " + message);
    }
}

struct LocalSPRResolvedAnchorEndpoints {
    int anchor_id = -1;
    int tip_id = -1;
};
bool local_spr_tip_id_is_valid(
    const TreeBuildResult& tree,
    int tip_id);
LocalSPRResolvedAnchorEndpoints resolve_local_spr_anchor_from_placement(
    const TreeBuildResult& tree,
    const mlipper::PlacementResult& placement);
LocalSPRResolvedAnchorEndpoints resolve_local_spr_anchor_endpoints(
    const TreeBuildResult& tree,
    const LocalSPRInsertionAnchor& anchor);

void collect_subtree_node_ids(
    const TreeBuildResult& tree,
    int node_id,
    std::vector<int>& node_ids)
{
    if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) return;
    node_ids.clear();
    std::vector<int> stack;
    stack.push_back(node_id);
    while (!stack.empty()) {
        const int cur = stack.back();
        stack.pop_back();
        if (cur < 0 || cur >= static_cast<int>(tree.nodes.size())) continue;
        node_ids.push_back(cur);
        const TreeNode& node = tree.nodes[(size_t)cur];
        if (node.right >= 0) stack.push_back(node.right);
        if (node.left >= 0) stack.push_back(node.left);
    }
}

bool subtree_fully_inside_mask(
    const TreeBuildResult& tree,
    int node_id,
    const std::vector<char>& mask,
    std::vector<int>* node_ids_out)
{
    if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) return false;
    if (mask.empty()) return false;

    std::vector<int> node_ids;
    std::vector<int>* sink = node_ids_out ? node_ids_out : &node_ids;
    collect_subtree_node_ids(tree, node_id, *sink);
    for (int cur : *sink) {
        if (cur < 0 || cur >= static_cast<int>(mask.size()) || !mask[(size_t)cur]) {
            return false;
        }
    }
    return true;
}

void rebuild_host_topology_from_tree_local(
    const TreeBuildResult& tree,
    HostPacking& host)
{
    host.postorder = tree.postorder;
    host.preorder = tree.preorder;
    host.parent.resize(tree.nodes.size(), -1);
    host.left.resize(tree.nodes.size(), -1);
    host.right.resize(tree.nodes.size(), -1);
    host.is_tip.resize(tree.nodes.size(), 0);
    host.blen.resize(tree.nodes.size(), fp_t(0));
    for (size_t node_idx = 0; node_idx < tree.nodes.size(); ++node_idx) {
        const TreeNode& node = tree.nodes[node_idx];
        host.parent[node_idx] = node.parent;
        host.left[node_idx] = node.left;
        host.right[node_idx] = node.right;
        host.is_tip[node_idx] = node.is_tip ? 1 : 0;
        host.blen[node_idx] = node.branch_length_to_parent;
    }
}

static inline void append_selected_downward_op(
    std::vector<NodeOpInfo>& ops,
    int parent_id,
    int left_id,
    int right_id,
    bool left_is_tip,
    bool right_is_tip,
    const std::vector<int>& node_to_tip,
    uint8_t dir_tag)
{
    NodeOpInfo op{};
    op.parent_id = parent_id;
    op.left_id = left_id;
    op.right_id = right_id;
    op.left_tip_index = left_is_tip ? node_to_tip[static_cast<size_t>(left_id)] : -1;
    op.right_tip_index = right_is_tip ? node_to_tip[static_cast<size_t>(right_id)] : -1;
    op.clv_pool = static_cast<uint8_t>(CLV_POOL_DOWN);
    op.dir_tag = dir_tag;

    const bool target_is_left = (dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const bool target_is_tip = target_is_left ? left_is_tip : right_is_tip;
    const bool sibling_is_tip = target_is_left ? right_is_tip : left_is_tip;
    if (target_is_tip && sibling_is_tip) {
        op.op_type = static_cast<int>(OP_DOWN_TIP_TIP);
    } else if (target_is_tip) {
        op.op_type = static_cast<int>(OP_DOWN_TIP_INNER);
    } else if (sibling_is_tip) {
        op.op_type = static_cast<int>(OP_DOWN_INNER_TIP);
    } else {
        op.op_type = static_cast<int>(OP_DOWN_INNER_INNER);
    }
    ops.push_back(op);
}

static void append_downward_op_for_child(
    const TreeBuildResult& tree,
    const std::vector<int>& node_to_tip,
    int child_id,
    std::vector<NodeOpInfo>& host_ops,
    const char* mismatch_message)
{
    if (child_id < 0 || child_id >= static_cast<int>(tree.nodes.size())) return;
    const int parent_id = tree.nodes[static_cast<size_t>(child_id)].parent;
    if (parent_id < 0 || parent_id >= static_cast<int>(tree.nodes.size())) return;

    const TreeNode& parent = tree.nodes[static_cast<size_t>(parent_id)];
    const int left_id = parent.left;
    const int right_id = parent.right;
    if (left_id < 0 || right_id < 0) return;

    const bool left_is_tip = tree.nodes[static_cast<size_t>(left_id)].is_tip;
    const bool right_is_tip = tree.nodes[static_cast<size_t>(right_id)].is_tip;
    if (left_id == child_id) {
        append_selected_downward_op(
            host_ops,
            parent_id,
            left_id,
            right_id,
            left_is_tip,
            right_is_tip,
            node_to_tip,
            static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    } else if (right_id == child_id) {
        append_selected_downward_op(
            host_ops,
            parent_id,
            left_id,
            right_id,
            left_is_tip,
            right_is_tip,
            node_to_tip,
            static_cast<uint8_t>(CLV_DIR_DOWN_RIGHT));
    } else {
        throw std::runtime_error(
            "Local SPR subtree assertion failed: " +
            std::string(mismatch_message));
    }
}

std::vector<int> build_local_spr_node_to_tip(
    const TreeBuildResult& tree,
    const HostPacking& host)
{
    std::vector<int> node_to_tip(tree.nodes.size(), -1);
    for (int tip_idx = 0; tip_idx < static_cast<int>(host.tip_node_ids.size()); ++tip_idx) {
        const int node_id = host.tip_node_ids[static_cast<size_t>(tip_idx)];
        if (node_id >= 0 && node_id < static_cast<int>(node_to_tip.size())) {
            node_to_tip[static_cast<size_t>(node_id)] = tip_idx;
        }
    }
    return node_to_tip;
}

void build_selected_downward_ops(
    const TreeBuildResult& tree,
    const std::vector<int>& node_to_tip,
    const std::vector<int>& target_child_ids,
    std::vector<NodeOpInfo>& host_ops)
{
    host_ops.clear();
    host_ops.reserve(target_child_ids.size());
    std::unordered_set<int> seen;
    for (int target_child_id : target_child_ids) {
        if (!seen.insert(target_child_id).second) continue;
        append_downward_op_for_child(
            tree,
            node_to_tip,
            target_child_id,
            host_ops,
            "target child does not hang under its recorded parent");
    }
}

void build_required_downward_update_ops(
    const TreeBuildResult& tree,
    const std::vector<int>& node_to_tip,
    const std::vector<int>& target_child_ids,
    std::vector<NodeOpInfo>& host_ops)
{
    std::vector<char> closure_mask(tree.nodes.size(), 0);
    for (int target_child_id : target_child_ids) {
        if (target_child_id < 0 || target_child_id >= static_cast<int>(tree.nodes.size())) continue;
        for (int node_id = target_child_id;
             node_id >= 0 && node_id < static_cast<int>(tree.nodes.size());
             node_id = tree.nodes[static_cast<size_t>(node_id)].parent) {
            const int parent_id = tree.nodes[static_cast<size_t>(node_id)].parent;
            if (parent_id < 0) break;
            closure_mask[static_cast<size_t>(node_id)] = 1;
        }
    }

    host_ops.clear();
    host_ops.reserve(mlipper::util::checked_product(
        "local SPR downward-operation reserve",
        target_child_ids.size(),
        size_t{4}));
    for (int node_id : tree.preorder) {
        if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) continue;
        if (!closure_mask[static_cast<size_t>(node_id)]) continue;
        append_downward_op_for_child(
            tree,
            node_to_tip,
            node_id,
            host_ops,
            "closure node does not hang under its recorded parent");
    }
}

struct LocalSPRNodeOpEqual {
    bool operator()(const NodeOpInfo& lhs, const NodeOpInfo& rhs) const {
        return lhs.parent_id == rhs.parent_id &&
               lhs.left_id == rhs.left_id &&
               lhs.right_id == rhs.right_id &&
               lhs.left_tip_index == rhs.left_tip_index &&
               lhs.right_tip_index == rhs.right_tip_index &&
               lhs.op_type == rhs.op_type &&
               lhs.clv_pool == rhs.clv_pool &&
               lhs.dir_tag == rhs.dir_tag;
    }
};

struct LocalSPRNodeOpHash {
    size_t operator()(const NodeOpInfo& key) const {
        size_t h = static_cast<size_t>(key.parent_id + 0x9e3779b9);
        auto mix = [&](size_t value) {
            h ^= value + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };
        mix(static_cast<size_t>(key.left_id));
        mix(static_cast<size_t>(key.right_id));
        mix(static_cast<size_t>(key.left_tip_index));
        mix(static_cast<size_t>(key.right_tip_index));
        mix(static_cast<size_t>(key.op_type));
        mix(static_cast<size_t>(key.clv_pool));
        mix(static_cast<size_t>(key.dir_tag));
        return h;
    }
};

std::vector<NodeOpInfo> filter_local_spr_new_ops(
    const std::vector<NodeOpInfo>& ops,
    const std::vector<NodeOpInfo>& already_computed_ops)
{
    std::unordered_set<NodeOpInfo, LocalSPRNodeOpHash, LocalSPRNodeOpEqual> seen;
    const size_t prior_capacity = mlipper::util::checked_product(
        "local SPR operation-set reserve",
        already_computed_ops.size(),
        size_t{2});
    seen.reserve(mlipper::util::checked_add_size(
        prior_capacity,
        ops.size(),
        "local SPR operation-set reserve"));
    for (const NodeOpInfo& op : already_computed_ops) {
        seen.insert(op);
    }

    std::vector<NodeOpInfo> filtered_ops;
    filtered_ops.reserve(ops.size());
    for (const NodeOpInfo& op : ops) {
        if (!seen.insert(op).second) {
            continue;
        }
        filtered_ops.push_back(op);
    }
    return filtered_ops;
}

std::vector<char> build_local_spr_subtree_mask(
    const TreeBuildResult& tree,
    const std::vector<int>& subtree_nodes)
{
    std::vector<char> subtree_mask(tree.nodes.size(), 0);
    for (int node_id : subtree_nodes) {
        if (node_id >= 0 && node_id < static_cast<int>(subtree_mask.size())) {
            subtree_mask[static_cast<size_t>(node_id)] = 1;
        }
    }
    return subtree_mask;
}

std::vector<int> filter_local_spr_candidate_edges(
    const TreeBuildResult& pruned_tree,
    const std::vector<char>& envelope_mask,
    const std::vector<char>& subtree_mask,
    const std::vector<int>& edge_candidates)
{
    if (envelope_mask.size() != pruned_tree.nodes.size() ||
        subtree_mask.size() != pruned_tree.nodes.size()) {
        throw std::invalid_argument(
            "Local SPR candidate masks must match the pruned tree.");
    }
    std::vector<int> legal_candidate_edges;
    legal_candidate_edges.reserve(edge_candidates.size());
    for (int edge_child : edge_candidates) {
        if (edge_child < 0 ||
            edge_child >= static_cast<int>(pruned_tree.nodes.size())) {
            continue;
        }
        const int edge_parent = pruned_tree.nodes[static_cast<size_t>(edge_child)].parent;
        if (edge_parent < 0) continue;
        if (!envelope_mask[static_cast<size_t>(edge_child)] ||
            !envelope_mask[static_cast<size_t>(edge_parent)]) {
            continue;
        }
        if (subtree_mask[static_cast<size_t>(edge_child)] ||
            subtree_mask[static_cast<size_t>(edge_parent)]) {
            continue;
        }
        legal_candidate_edges.push_back(edge_child);
    }
    return legal_candidate_edges;
}

HostPacking build_local_spr_tree_host_packing(
    const TopologyRefinementRunContext& ctx,
    const TreeBuildResult& tree)
{
    TreeBuildResult packing_tree = tree;
    for (int node_id : ctx.directional_boundary_node_ids) {
        if (node_id < 0 ||
            node_id >= static_cast<int>(packing_tree.nodes.size())) {
            throw std::runtime_error(
                "build_local_spr_tree_host_packing: virtual boundary id is out of range.");
        }
        packing_tree.nodes[static_cast<size_t>(node_id)].is_tip = true;
    }
    HostPacking host_pack = pack_host_arrays_from_tree_and_msa(
        packing_tree,
        ctx.current_names,
        ctx.current_rows,
        ctx.sites,
        ctx.states);
    for (int node_id : ctx.directional_boundary_node_ids) {
        host_pack.is_tip[static_cast<size_t>(node_id)] = 0;
    }
    host_pack.pattern_weights = ctx.pattern_weights_arg;
    fill_pmats_in_host_packing(
        tree,
        host_pack,
        ctx.state.eig,
        ctx.rate_multipliers,
        ctx.states,
        ctx.rate_cats);
    return host_pack;
}

void local_spr_assert_candidate_legal(
    const TreeBuildResult& tree,
    const LocalSPRRepairUnit& unit,
    int prune_root_id,
    const std::vector<int>& subtree_nodes,
    const std::vector<int>& legal_edges,
    int regraft_child_id)
{
    local_spr_assert(
        prune_root_id >= 0 && prune_root_id < static_cast<int>(tree.nodes.size()),
        "prune root out of range");
    local_spr_assert(
        prune_root_id != tree.root_id,
        "prune root cannot be the tree root");
    local_spr_assert(
        tree.nodes[static_cast<size_t>(prune_root_id)].parent >= 0,
        "prune root must have a parent");

    std::vector<int> subtree_nodes_check;
    local_spr_assert(
        subtree_fully_inside_mask(tree, prune_root_id, unit.envelope_mask, &subtree_nodes_check),
        "committed subtree is not fully contained in its repair envelope");
    local_spr_assert(
        subtree_nodes_check == subtree_nodes,
        "subtree node set drifted between enumeration and commit");

    local_spr_assert(
        regraft_child_id >= 0 && regraft_child_id < static_cast<int>(tree.nodes.size()),
        "regraft child is out of range");
    const int regraft_parent_id = tree.nodes[static_cast<size_t>(regraft_child_id)].parent;
    local_spr_assert(regraft_parent_id >= 0, "regraft edge has no parent");
    local_spr_assert(
        unit.envelope_mask[static_cast<size_t>(regraft_child_id)] &&
        unit.envelope_mask[static_cast<size_t>(regraft_parent_id)],
        "regraft edge escaped the repair envelope");

    const std::vector<char> subtree_mask =
        build_local_spr_subtree_mask(tree, subtree_nodes);
    local_spr_assert(
        !subtree_mask[static_cast<size_t>(regraft_child_id)] &&
        !subtree_mask[static_cast<size_t>(regraft_parent_id)],
        "self-regraft detected: target edge lies inside pruned subtree");
    local_spr_assert(
        std::find(legal_edges.begin(), legal_edges.end(), regraft_child_id) != legal_edges.end(),
        "selected regraft edge is not legal under the bounded-radius search");
}

std::vector<int> bfs_distances(const TreeBuildResult& tree, int start) {
    const int node_count = static_cast<int>(tree.nodes.size());
    std::vector<int> dist((size_t)node_count, -1);
    if (start < 0 || start >= node_count) return dist;
    std::vector<int> queue;
    queue.reserve((size_t)node_count);
    dist[(size_t)start] = 0;
    queue.push_back(start);
    for (size_t qi = 0; qi < queue.size(); ++qi) {
        const int cur = queue[qi];
        const int next_dist = dist[(size_t)cur] + 1;
        const TreeNode& node = tree.nodes[cur];
        const int neighbors[3] = {node.parent, node.left, node.right};
        for (int neighbor : neighbors) {
            if (neighbor < 0 || neighbor >= node_count) continue;
            if (dist[(size_t)neighbor] >= 0) continue;
            dist[(size_t)neighbor] = next_dist;
            queue.push_back(neighbor);
        }
    }
    return dist;
}

std::vector<LocalSPRInsertionAnchor> build_local_spr_insertion_anchors(
    const TreeBuildResult& tree,
    const std::vector<mlipper::PlacementResult>& committed_placements)
{
    std::vector<LocalSPRInsertionAnchor> anchors;
    anchors.reserve(committed_placements.size());
    for (const mlipper::PlacementResult& placement :
         committed_placements) {
        const LocalSPRResolvedAnchorEndpoints resolved_anchor =
            resolve_local_spr_anchor_from_placement(tree, placement);
        anchors.push_back(LocalSPRInsertionAnchor{
            resolved_anchor.anchor_id,
            resolved_anchor.tip_id,
        });
    }
    return anchors;
}

std::vector<LocalSPRInsertionAnchor> refresh_local_spr_anchor_parents(
    const TreeBuildResult& tree,
    const std::vector<LocalSPRInsertionAnchor>& anchors)
{
    std::vector<LocalSPRInsertionAnchor> refreshed;
    refreshed.reserve(anchors.size());
    for (const LocalSPRInsertionAnchor& anchor : anchors) {
        if (anchor.endpoint_a < 0 || anchor.endpoint_b < 0 ||
            anchor.endpoint_a >= static_cast<int>(tree.nodes.size()) ||
            anchor.endpoint_b >= static_cast<int>(tree.nodes.size())) {
            throw std::runtime_error(
                "Local SPR anchor endpoint is out of range.");
        }
        const bool endpoint_a_is_tip =
            tree.nodes[static_cast<size_t>(anchor.endpoint_a)].is_tip;
        const bool endpoint_b_is_tip =
            tree.nodes[static_cast<size_t>(anchor.endpoint_b)].is_tip;
        if (endpoint_a_is_tip == endpoint_b_is_tip) {
            throw std::runtime_error(
                "Local SPR anchor must identify exactly one tip.");
        }
        const int tip_id = endpoint_a_is_tip
            ? anchor.endpoint_a
            : anchor.endpoint_b;
        const int parent_id =
            tree.nodes[static_cast<size_t>(tip_id)].parent;
        if (parent_id < 0) {
            throw std::runtime_error(
                "Local SPR anchor tip has no current parent.");
        }
        refreshed.push_back({parent_id, tip_id});
    }
    return refreshed;
}

// Cluster nearby insertion anchors with a disjoint set, then form each unit's
// search envelope as the union of nodes within envelope_radius of any anchor
// edge. Clustering controls which insertions compete as one repair problem;
// the envelope independently limits the legal prune/regraft neighborhood.
std::vector<LocalSPRRepairUnit> build_local_spr_repair_units(
    const TreeBuildResult& tree,
    const std::vector<LocalSPRInsertionAnchor>& anchors,
    int cluster_threshold,
    int envelope_radius)
{
    std::vector<LocalSPRRepairUnit> repair_units;
    if (tree.nodes.empty()) {
        return repair_units;
    }

    const int anchor_count = static_cast<int>(anchors.size());
    if (anchor_count == 0) {
        return repair_units;
    }

    repair_units.reserve(anchors.size());
    cluster_threshold = std::max(0, cluster_threshold);
    const int effective_radius = std::max(envelope_radius, 1);
    IntDisjointSet dsu(anchor_count);
    std::vector<char> valid_anchor(static_cast<size_t>(anchor_count), 0);
    std::vector<LocalSPRResolvedAnchorEndpoints> resolved_anchors(
        static_cast<size_t>(anchor_count));
    std::vector<std::vector<int>> anchor_center_distances(
        static_cast<size_t>(anchor_count));
    std::vector<std::vector<int>> anchor_edge_distances(
        static_cast<size_t>(anchor_count));

    for (int anchor_idx = 0; anchor_idx < anchor_count; ++anchor_idx) {
        const LocalSPRInsertionAnchor& anchor =
            anchors[static_cast<size_t>(anchor_idx)];
        if (anchor.endpoint_a < 0 ||
            anchor.endpoint_a >= static_cast<int>(tree.nodes.size()) ||
            anchor.endpoint_b < 0 ||
            anchor.endpoint_b >= static_cast<int>(tree.nodes.size())) {
            continue;
        }
        const LocalSPRResolvedAnchorEndpoints resolved_anchor =
            resolve_local_spr_anchor_endpoints(tree, anchor);
        resolved_anchors[static_cast<size_t>(anchor_idx)] = resolved_anchor;
        valid_anchor[static_cast<size_t>(anchor_idx)] = 1;

        anchor_center_distances[static_cast<size_t>(anchor_idx)] =
            bfs_distances(tree, resolved_anchor.anchor_id);

        std::vector<char> source_mask(tree.nodes.size(), 0);
        source_mask[static_cast<size_t>(anchor.endpoint_a)] = 1;
        source_mask[static_cast<size_t>(anchor.endpoint_b)] = 1;
        anchor_edge_distances[static_cast<size_t>(anchor_idx)] =
            mlipper::divide_and_conquer::multi_source_bfs_distances(
                tree,
                source_mask);
    }

    for (int i = 0; i < anchor_count; ++i) {
        if (valid_anchor[static_cast<size_t>(i)] == 0) {
            continue;
        }
        for (int j = i + 1; j < anchor_count; ++j) {
            if (valid_anchor[static_cast<size_t>(j)] == 0) {
                continue;
            }
            const int other_anchor_id =
                resolved_anchors[static_cast<size_t>(j)].anchor_id;
            const int distance =
                anchor_center_distances[static_cast<size_t>(i)]
                    [static_cast<size_t>(other_anchor_id)];
            if (distance >= 0 && distance <= cluster_threshold) {
                dsu.unite(i, j);
            }
        }
    }

    std::vector<int> root_to_unit(static_cast<size_t>(anchor_count), -1);
    for (int anchor_idx = 0; anchor_idx < anchor_count; ++anchor_idx) {
        if (valid_anchor[static_cast<size_t>(anchor_idx)] == 0) {
            continue;
        }
        const int root = dsu.find(anchor_idx);
        int unit_index = root_to_unit[static_cast<size_t>(root)];
        if (unit_index < 0) {
            unit_index = static_cast<int>(repair_units.size());
            root_to_unit[static_cast<size_t>(root)] = unit_index;
            LocalSPRRepairUnit unit;
            unit.unit_id = unit_index;
            unit.envelope_mask.assign(tree.nodes.size(), 0);
            unit.envelope_nodes.reserve(tree.nodes.size());
            repair_units.push_back(std::move(unit));
        }

        LocalSPRRepairUnit& unit =
            repair_units[static_cast<size_t>(unit_index)];
        unit.anchors.push_back(anchors[static_cast<size_t>(anchor_idx)]);
        unit.anchor_indices.push_back(anchor_idx);
    }

    for (LocalSPRRepairUnit& unit : repair_units) {
        for (int anchor_idx : unit.anchor_indices) {
            const std::vector<int>& node_distances =
                anchor_edge_distances[static_cast<size_t>(anchor_idx)];
            if (node_distances.empty()) {
                continue;
            }
            for (int node_id = 0;
                 node_id < static_cast<int>(tree.nodes.size());
                 ++node_id) {
                const int distance = node_distances[static_cast<size_t>(node_id)];
                if (distance < 0 || distance > effective_radius) {
                    continue;
                }
                unit.envelope_mask[static_cast<size_t>(node_id)] = 1;
            }
        }
        for (int node_id = 0;
             node_id < static_cast<int>(tree.nodes.size());
             ++node_id) {
            if (unit.envelope_mask[static_cast<size_t>(node_id)] != 0) {
                unit.envelope_nodes.push_back(node_id);
            }
        }
    }

    repair_units.erase(
        std::remove_if(
            repair_units.begin(),
            repair_units.end(),
            [](const LocalSPRRepairUnit& unit) {
                return unit.envelope_nodes.empty();
            }),
        repair_units.end());
    for (size_t unit_idx = 0; unit_idx < repair_units.size(); ++unit_idx) {
        repair_units[unit_idx].unit_id = static_cast<int>(unit_idx);
    }

    return repair_units;
}

std::vector<LocalSPRRepairUnit> build_nni_repair_units(
    const TreeBuildResult& tree,
    const std::unordered_set<int>& central_edge_child_ids)
{
    if (central_edge_child_ids.empty()) {
        return {};
    }

    LocalSPRRepairUnit unit;
    unit.unit_id = 0;
    unit.envelope_mask.assign(tree.nodes.size(), 1);
    unit.envelope_nodes.reserve(tree.nodes.size());
    for (int node_id = 0; node_id < static_cast<int>(tree.nodes.size()); ++node_id) {
        unit.envelope_nodes.push_back(node_id);
    }
    return {std::move(unit)};
}

bool local_spr_candidate_better(
    const LocalSPRCandidateMove& lhs,
    const LocalSPRCandidateMove& rhs)
{
    if (lhs.approx_gain != rhs.approx_gain) {
        return lhs.approx_gain > rhs.approx_gain;
    }
    if (lhs.prune_root_id != rhs.prune_root_id) {
        return lhs.prune_root_id < rhs.prune_root_id;
    }
    if (lhs.regraft_child_id != rhs.regraft_child_id) {
        return lhs.regraft_child_id < rhs.regraft_child_id;
    }
    if (lhs.regraft_parent_id != rhs.regraft_parent_id) {
        return lhs.regraft_parent_id < rhs.regraft_parent_id;
    }
    if (lhs.repair_unit_id != rhs.repair_unit_id) {
        return lhs.repair_unit_id < rhs.repair_unit_id;
    }
    return lhs.old_parent_id < rhs.old_parent_id;
}

void keep_local_spr_topk(
    std::vector<LocalSPRCandidateMove>& topk,
    LocalSPRCandidateMove candidate,
    int topk_limit)
{
    topk.push_back(std::move(candidate));
    std::sort(
        topk.begin(),
        topk.end(),
        local_spr_candidate_better);
    if (static_cast<int>(topk.size()) > topk_limit) {
        topk.resize((size_t)topk_limit);
    }
}

std::vector<int> select_local_spr_seed_edges(
    const RawPlacementResult& placement_result,
    const TreeBuildResult& tree,
    const std::vector<int>& center_dist,
    double baseline_logL,
    int seed_limit,
    int required_edge_radius)
{
    std::vector<std::pair<double, int>> ranked_edges;
    ranked_edges.reserve(placement_result.top_placements.size());
    for (const RawPlacementResult::RankedPlacement& placement :
         placement_result.top_placements) {
        const int edge_child_id = placement.target_id;
        if (edge_child_id < 0 || edge_child_id >= static_cast<int>(tree.nodes.size())) {
            continue;
        }
        const int edge_parent_id = tree.nodes[(size_t)edge_child_id].parent;
        if (edge_parent_id < 0) {
            continue;
        }
        if (required_edge_radius >= 0) {
            const int edge_child_dist = center_dist[(size_t)edge_child_id];
            const int edge_parent_dist = center_dist[(size_t)edge_parent_id];
            int edge_radius = -1;
            if (edge_child_dist < 0) {
                edge_radius = edge_parent_dist;
            } else if (edge_parent_dist < 0) {
                edge_radius = edge_child_dist;
            } else {
                edge_radius = std::min(edge_child_dist, edge_parent_dist);
            }
            if (edge_radius != required_edge_radius) {
                continue;
            }
        }
        const double approx_gain = placement.loglikelihood - baseline_logL;
        if (!(approx_gain > 0.0)) {
            continue;
        }
        ranked_edges.emplace_back(approx_gain, edge_child_id);
    }

    std::sort(
        ranked_edges.begin(),
        ranked_edges.end(),
        [](const std::pair<double, int>& lhs, const std::pair<double, int>& rhs) {
            if (lhs.first == rhs.first) return lhs.second < rhs.second;
            return lhs.first > rhs.first;
        });

    std::vector<int> seed_edges;
    seed_edges.reserve(ranked_edges.size());
    std::unordered_set<int> seen;
    for (const auto& [gain, edge_child_id] : ranked_edges) {
        (void)gain;
        if (!seen.insert(edge_child_id).second) continue;
        seed_edges.push_back(edge_child_id);
        if (seed_limit > 0 && static_cast<int>(seed_edges.size()) >= seed_limit) {
            break;
        }
    }
    return seed_edges;
}

// Greedily choose mutually independent moves from an already ranked list.
// At most one move is selected per repair unit, and selected prune subtrees or
// regraft paths may not overlap. This permits later application without one
// move invalidating another move selected in the same round.
std::vector<LocalSPRCandidateMove> select_local_spr_candidates(
    const std::vector<LocalSPRCandidateMove>& ranked_candidates,
    int node_count)
{
    std::vector<LocalSPRCandidateMove> selected;
    std::vector<char> selected_subtree_mask((size_t)std::max(0, node_count), 0);
    std::vector<char> selected_path_mask((size_t)std::max(0, node_count), 0);
    std::unordered_set<int> used_units;

    for (const LocalSPRCandidateMove& candidate : ranked_candidates) {
        if (!used_units.insert(candidate.repair_unit_id).second) {
            continue;
        }

        bool conflict = false;
        for (int node_id : candidate.subtree_nodes) {
            if (node_id >= 0 &&
                node_id < node_count &&
                selected_subtree_mask[(size_t)node_id]) {
                conflict = true;
                break;
            }
        }
        if (conflict) {
            used_units.erase(candidate.repair_unit_id);
            continue;
        }
        if (candidate.regraft_child_id >= 0 &&
            candidate.regraft_child_id < node_count &&
            selected_subtree_mask[(size_t)candidate.regraft_child_id]) {
            used_units.erase(candidate.repair_unit_id);
            continue;
        }
        if (candidate.regraft_parent_id >= 0 &&
            candidate.regraft_parent_id < node_count &&
            selected_subtree_mask[(size_t)candidate.regraft_parent_id]) {
            used_units.erase(candidate.repair_unit_id);
            continue;
        }
        for (int node_id : candidate.regraft_path_nodes) {
            if (node_id >= 0 &&
                node_id < node_count &&
                selected_path_mask[(size_t)node_id]) {
                conflict = true;
                break;
            }
        }
        if (conflict) {
            used_units.erase(candidate.repair_unit_id);
            continue;
        }

        selected.push_back(candidate);
        for (int node_id : candidate.subtree_nodes) {
            if (node_id >= 0 && node_id < node_count) {
                selected_subtree_mask[(size_t)node_id] = 1;
            }
        }
        for (int node_id : candidate.regraft_path_nodes) {
            if (node_id >= 0 && node_id < node_count) {
                selected_path_mask[(size_t)node_id] = 1;
            }
        }
    }

    return selected;
}

// Revalidate cached node relationships after earlier accepted moves. Candidate
// IDs are deliberately not trusted across topology mutations; a mismatch is a
// normal stale-candidate result rather than an invariant failure.
bool local_spr_candidate_still_legal(
    const TreeBuildResult& tree,
    const LocalSPRRepairUnit& unit,
    const LocalSPRCandidateMove& candidate,
    std::vector<int>* current_subtree_nodes_out)
{
    if (candidate.prune_root_id < 0 ||
        candidate.prune_root_id >= static_cast<int>(tree.nodes.size()) ||
        candidate.prune_root_id == tree.root_id) {
        return false;
    }
    const int current_parent =
        tree.nodes[(size_t)candidate.prune_root_id].parent;
    if (current_parent < 0 ||
        (candidate.old_parent_id >= 0 && current_parent != candidate.old_parent_id)) {
        return false;
    }

    std::vector<int> current_subtree_nodes_storage;
    std::vector<int>* current_subtree_nodes =
        current_subtree_nodes_out != nullptr
            ? current_subtree_nodes_out
            : &current_subtree_nodes_storage;
    if (!subtree_fully_inside_mask(
            tree,
            candidate.prune_root_id,
            unit.envelope_mask,
            current_subtree_nodes)) {
        return false;
    }

    if (candidate.regraft_child_id < 0 ||
        candidate.regraft_child_id >= static_cast<int>(tree.nodes.size())) {
        return false;
    }
    const int current_regraft_parent =
        tree.nodes[(size_t)candidate.regraft_child_id].parent;
    if (current_regraft_parent < 0 ||
        (candidate.regraft_parent_id >= 0 &&
         current_regraft_parent != candidate.regraft_parent_id)) {
        return false;
    }
    if (!unit.envelope_mask[(size_t)candidate.regraft_child_id] ||
        !unit.envelope_mask[(size_t)current_regraft_parent]) {
        return false;
    }

    const std::vector<char> current_subtree_mask =
        build_local_spr_subtree_mask(tree, *current_subtree_nodes);
    if (current_subtree_mask[(size_t)candidate.regraft_child_id] ||
        current_subtree_mask[(size_t)current_regraft_parent]) {
        return false;
    }

    return true;
}

// Detach pruned_id and suppress its former parent, preserving the represented
// unrooted edge length by adding the two adjacent branches. Validation happens
// before mutation so false guarantees that tree is unchanged.
bool prune_subtree_for_spr(
    TreeBuildResult& tree,
    int pruned_id,
    PruneInfo& info)
{
    if (pruned_id < 0 ||
        pruned_id >= static_cast<int>(tree.nodes.size())) {
        return false;
    }
    if (pruned_id == tree.root_id) return false;
    const int parent_id = tree.nodes[pruned_id].parent;
    if (parent_id < 0 ||
        parent_id >= static_cast<int>(tree.nodes.size())) {
        return false;
    }
    const TreeNode& parent = tree.nodes[parent_id];
    const int sibling_id = parent.left == pruned_id
        ? parent.right
        : (parent.right == pruned_id ? parent.left : -1);
    if (sibling_id < 0 ||
        sibling_id >= static_cast<int>(tree.nodes.size())) {
        return false;
    }

    const int grandparent_id = parent.parent;
    // Complete every structural check before the first mutation. Callers use
    // false to skip an invalid candidate and expect the input tree unchanged.
    if (grandparent_id >= 0) {
        if (grandparent_id >= static_cast<int>(tree.nodes.size())) {
            return false;
        }
        const TreeNode& grandparent = tree.nodes[grandparent_id];
        if (grandparent.left != parent_id &&
            grandparent.right != parent_id) {
            return false;
        }
    }

    info.pruned_id = pruned_id;
    info.free_internal_id = parent_id;
    info.sibling_id = sibling_id;
    info.grandparent_id = grandparent_id;
    tree.nodes[pruned_id].parent = -1;

    if (grandparent_id >= 0) {
        TreeNode& grandparent = tree.nodes[grandparent_id];
        if (grandparent.left == parent_id) {
            grandparent.left = sibling_id;
        } else if (grandparent.right == parent_id) {
            grandparent.right = sibling_id;
        }
        TreeNode& sibling = tree.nodes[sibling_id];
        sibling.parent = grandparent_id;
        sibling.branch_length_to_parent += parent.branch_length_to_parent;
    } else {
        tree.root_id = sibling_id;
        TreeNode& sibling = tree.nodes[sibling_id];
        sibling.parent = -1;
        sibling.branch_length_to_parent = fp_t(0);
    }

    TreeNode& free_internal = tree.nodes[parent_id];
    free_internal.left = -1;
    free_internal.right = -1;
    free_internal.parent = -1;
    free_internal.is_tip = false;
    free_internal.name.clear();
    free_internal.branch_length_to_parent = fp_t(0);
    rebuild_traversals(tree);
    return true;
}

std::vector<char> build_unit_anchor_skeleton_mask(
    const TreeBuildResult& tree,
    const std::vector<LocalSPRInsertionAnchor>& anchors)
{
    std::vector<char> mask(tree.nodes.size(), 0);
    if (anchors.empty()) return mask;

    int root_anchor_id = -1;
    for (const LocalSPRInsertionAnchor& anchor : anchors) {
        if (anchor.endpoint_a < 0 ||
            anchor.endpoint_a >= static_cast<int>(tree.nodes.size()) ||
            anchor.endpoint_b < 0 ||
            anchor.endpoint_b >= static_cast<int>(tree.nodes.size())) {
            continue;
        }
        const LocalSPRResolvedAnchorEndpoints resolved_anchor =
            resolve_local_spr_anchor_endpoints(tree, anchor);
        if (root_anchor_id < 0) {
            root_anchor_id = resolved_anchor.anchor_id;
        }

        const std::vector<int> path_nodes =
            mlipper::divide_and_conquer::build_tree_path_nodes(
                tree,
                resolved_anchor.anchor_id,
                resolved_anchor.tip_id);
        for (int node_id : path_nodes) {
            if (node_id < 0 || node_id >= static_cast<int>(mask.size())) continue;
            mask[(size_t)node_id] = 1;
        }

        if (root_anchor_id >= 0 &&
            root_anchor_id != resolved_anchor.anchor_id) {
            const std::vector<int> anchor_bridge_nodes =
                mlipper::divide_and_conquer::build_tree_path_nodes(
                    tree,
                    root_anchor_id,
                    resolved_anchor.anchor_id);
            for (int node_id : anchor_bridge_nodes) {
                if (node_id < 0 || node_id >= static_cast<int>(mask.size())) continue;
                mask[(size_t)node_id] = 1;
            }
        }
    }
    return mask;
}

bool local_spr_tip_id_is_valid(
    const TreeBuildResult& tree,
    int tip_id)
{
    return tip_id >= 0 &&
           tip_id < static_cast<int>(tree.nodes.size()) &&
           tree.nodes[static_cast<size_t>(tip_id)].is_tip;
}

int local_spr_find_node_id_by_stable_label(
    const TreeBuildResult& tree,
    int stable_node_label)
{
    if (stable_node_label < 0) {
        return -1;
    }
    for (const TreeNode& node : tree.nodes) {
        if (node.stable_node_label == stable_node_label) {
            return node.id;
        }
    }
    return -1;
}

LocalSPRResolvedAnchorEndpoints resolve_local_spr_anchor_from_placement(
    const TreeBuildResult& tree,
    const mlipper::PlacementResult& placement)
{
    const bool has_tip_name = !placement.query_name.empty();
    int tip_id = -1;
    if (has_tip_name) {
        const auto tip_it = tree.tip_node_by_name.find(placement.query_name);
        if (tip_it != tree.tip_node_by_name.end() &&
            local_spr_tip_id_is_valid(tree, tip_it->second)) {
            tip_id = tip_it->second;
        }
    }

    if (tip_id < 0) {
        const int resolved_tip_id =
            local_spr_find_node_id_by_stable_label(tree, placement.query_tip_label);
        if (local_spr_tip_id_is_valid(tree, resolved_tip_id)) {
            tip_id = resolved_tip_id;
        }
    }

    if (tip_id < 0) {
        if (!has_tip_name) {
            throw std::runtime_error(
                "Local SPR committed placement is missing a resolvable tip identity.");
        }
        throw std::runtime_error(
            "Local SPR could not find committed tip in tree: " +
            placement.query_name);
    }

    int anchor_id = -1;
    const int current_parent_id =
        tree.nodes[static_cast<size_t>(tip_id)].parent;
    if (current_parent_id >= 0 &&
        current_parent_id < static_cast<int>(tree.nodes.size()) &&
        tree.nodes[static_cast<size_t>(current_parent_id)].stable_node_label ==
            placement.attachment_node_label) {
        anchor_id = current_parent_id;
    }

    if (anchor_id < 0) {
        anchor_id = current_parent_id;
    }
    if (anchor_id < 0 || anchor_id >= static_cast<int>(tree.nodes.size())) {
        throw std::runtime_error(
            "Local SPR committed tip is missing a valid anchor parent: " +
            placement.query_name);
    }
    return LocalSPRResolvedAnchorEndpoints{anchor_id, tip_id};
}

LocalSPRResolvedAnchorEndpoints resolve_local_spr_anchor_endpoints(
    const TreeBuildResult& tree,
    const LocalSPRInsertionAnchor& anchor)
{
    if (anchor.endpoint_a < 0 ||
        anchor.endpoint_a >= static_cast<int>(tree.nodes.size()) ||
        anchor.endpoint_b < 0 ||
        anchor.endpoint_b >= static_cast<int>(tree.nodes.size())) {
        throw std::runtime_error(
            "Local SPR anchor endpoint is out of range.");
    }

    const bool a_is_tip =
        tree.nodes[static_cast<size_t>(anchor.endpoint_a)].is_tip;
    const bool b_is_tip =
        tree.nodes[static_cast<size_t>(anchor.endpoint_b)].is_tip;
    if (a_is_tip &&
        tree.nodes[static_cast<size_t>(anchor.endpoint_a)].parent ==
            anchor.endpoint_b) {
        return LocalSPRResolvedAnchorEndpoints{
            anchor.endpoint_b,
            anchor.endpoint_a,
        };
    }
    if (b_is_tip &&
        tree.nodes[static_cast<size_t>(anchor.endpoint_b)].parent ==
            anchor.endpoint_a) {
        return LocalSPRResolvedAnchorEndpoints{
            anchor.endpoint_a,
            anchor.endpoint_b,
        };
    }

    throw std::runtime_error(
        "Local SPR anchor endpoints must resolve to a committed tip and its current parent.");
}

int edge_distance_from_distances(const std::vector<int>& dist, int child_id, int parent_id) {
    if (child_id < 0 || parent_id < 0) return -1;
    const int d_child = dist[(size_t)child_id];
    const int d_parent = dist[(size_t)parent_id];
    if (d_child < 0) return d_parent;
    if (d_parent < 0) return d_child;
    return std::min(d_child, d_parent);
}

std::vector<int> collect_candidate_edges(
    const TreeBuildResult& tree,
    int endpoint_a,
    int endpoint_b,
    int radius,
    int exclude_a,
    int exclude_b)
{
    std::vector<int> edges;
    if (radius < 0) return edges;
    const int node_count = static_cast<int>(tree.nodes.size());
    std::vector<int> dist_a = bfs_distances(tree, endpoint_a);
    std::vector<int> dist_b;
    if (endpoint_b >= 0 && endpoint_b != endpoint_a) {
        dist_b = bfs_distances(tree, endpoint_b);
    }
    for (int node_id = 0; node_id < node_count; ++node_id) {
        const int parent_id = tree.nodes[node_id].parent;
        if (parent_id < 0) continue;
        if (node_id == exclude_a || node_id == exclude_b) continue;
        if (parent_id == exclude_a || parent_id == exclude_b) continue;
        int best = edge_distance_from_distances(dist_a, node_id, parent_id);
        if (!dist_b.empty()) {
            const int alt = edge_distance_from_distances(dist_b, node_id, parent_id);
            if (best < 0 || (alt >= 0 && alt < best)) best = alt;
        }
        if (best >= 0 && best <= radius) {
            edges.push_back(node_id);
        }
    }
    return edges;
}

std::vector<int> compute_center_distances(
    const TreeBuildResult& tree,
    int endpoint_a,
    int endpoint_b)
{
    std::vector<int> best = bfs_distances(tree, endpoint_a);
    if (endpoint_b >= 0 && endpoint_b != endpoint_a) {
        const std::vector<int> alt = bfs_distances(tree, endpoint_b);
        if (best.size() != alt.size()) {
            throw std::runtime_error("compute_center_distances produced mismatched BFS arrays.");
        }
        for (size_t i = 0; i < best.size(); ++i) {
            if (best[i] < 0) {
                best[i] = alt[i];
            } else if (alt[i] >= 0 && alt[i] < best[i]) {
                best[i] = alt[i];
            }
        }
    }
    return best;
}

std::vector<int> collect_outward_edges_from_seed(
    const TreeBuildResult& tree,
    const std::vector<int>& center_dist,
    int seed_edge_child_id,
    int min_radius_exclusive,
    int max_radius,
    int exclude_a,
    int exclude_b)
{
    std::vector<int> edges;
    if (max_radius <= min_radius_exclusive) return edges;
    if (seed_edge_child_id < 0 ||
        seed_edge_child_id >= static_cast<int>(tree.nodes.size())) {
        return edges;
    }

    const int seed_parent_id = tree.nodes[(size_t)seed_edge_child_id].parent;
    if (seed_parent_id < 0 ||
        seed_parent_id >= static_cast<int>(tree.nodes.size())) {
        return edges;
    }

    const int child_dist = center_dist[(size_t)seed_edge_child_id];
    const int parent_dist = center_dist[(size_t)seed_parent_id];
    if (child_dist < 0 || parent_dist < 0 || child_dist == parent_dist) {
        return edges;
    }

    struct FrontierState {
        int node_id = -1, prev_id = -1;
    };

    std::vector<FrontierState> queue;
    queue.reserve(tree.nodes.size());
    std::vector<char> visited(tree.nodes.size(), 0);
    std::vector<char> edge_seen(tree.nodes.size(), 0);

    auto push_start = [&](int node_id, int prev_id) {
        if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) return;
        if (visited[(size_t)node_id]) return;
        visited[(size_t)node_id] = 1;
        queue.push_back(FrontierState{node_id, prev_id});
    };

    if (child_dist > parent_dist) {
        push_start(seed_edge_child_id, seed_parent_id);
    } else {
        push_start(seed_parent_id, seed_edge_child_id);
    }

    auto maybe_record_edge = [&](int node_a, int node_b) {
        int edge_child_id = -1;
        int edge_parent_id = -1;
        if (tree.nodes[(size_t)node_a].parent == node_b) {
            edge_child_id = node_a;
            edge_parent_id = node_b;
        } else if (tree.nodes[(size_t)node_b].parent == node_a) {
            edge_child_id = node_b;
            edge_parent_id = node_a;
        } else {
            return;
        }

        if (edge_child_id == exclude_a || edge_child_id == exclude_b) return;
        if (edge_parent_id == exclude_a || edge_parent_id == exclude_b) return;
        if (edge_seen[(size_t)edge_child_id]) return;

        const int edge_radius =
            edge_distance_from_distances(center_dist, edge_child_id, edge_parent_id);
        if (edge_radius <= min_radius_exclusive || edge_radius > max_radius) return;

        edge_seen[(size_t)edge_child_id] = 1;
        edges.push_back(edge_child_id);
    };

    for (size_t qi = 0; qi < queue.size(); ++qi) {
        const FrontierState state = queue[qi];
        const int cur = state.node_id;
        const int cur_dist = center_dist[(size_t)cur];
        const TreeNode& node = tree.nodes[(size_t)cur];
        const int neighbors[3] = {node.parent, node.left, node.right};
        for (int neighbor : neighbors) {
            if (neighbor < 0 || neighbor == state.prev_id) continue;
            if (neighbor >= static_cast<int>(tree.nodes.size())) continue;
            const int neighbor_dist = center_dist[(size_t)neighbor];
            if (neighbor_dist < 0 || neighbor_dist <= cur_dist) continue;

            maybe_record_edge(cur, neighbor);
            if (!visited[(size_t)neighbor]) {
                visited[(size_t)neighbor] = 1;
                queue.push_back(FrontierState{neighbor, cur});
            }
        }
    }

    return edges;
}

void regraft_subtree_for_spr(
    TreeBuildResult& tree,
    const PruneInfo& info,
    int target_child_id,
    double pendant_length,
    double proximal_length,
    double pruned_branch_min,
    double target_branch_min)
{
    const size_t node_count = tree.nodes.size();
    if (target_child_id < 0 ||
        static_cast<size_t>(target_child_id) >= node_count ||
        info.free_internal_id < 0 ||
        static_cast<size_t>(info.free_internal_id) >= node_count ||
        info.pruned_id < 0 ||
        static_cast<size_t>(info.pruned_id) >= node_count) {
        throw std::invalid_argument("SPR regraft node ID is out of range.");
    }
    const int parent_id = tree.nodes[target_child_id].parent;
    if (parent_id < 0 || static_cast<size_t>(parent_id) >= node_count) {
        throw std::runtime_error("SPR regraft target has no parent edge.");
    }
    const int free_internal_id = info.free_internal_id;
    const int pruned_id = info.pruned_id;

    TreeNode& parent = tree.nodes[parent_id];
    if (parent.left == target_child_id) {
        parent.left = free_internal_id;
    } else if (parent.right == target_child_id) {
        parent.right = free_internal_id;
    } else {
        throw std::runtime_error("SPR regraft target not found under its parent.");
    }

    TreeNode& internal = tree.nodes[free_internal_id];
    internal.parent = parent_id;
    internal.left = target_child_id;
    internal.right = pruned_id;
    internal.is_tip = false;
    internal.name.clear();

    TreeNode& target_child = tree.nodes[target_child_id];
    const fp_t total_length = target_child.branch_length_to_parent;
    const double total = static_cast<double>(total_length);
    double proximal = 0.0;
    double distal = 0.0;
    normalize_split_branch_lengths(
        total,
        proximal_length,
        scalar_max(TOPOLOGY_INTERNAL_BRANCH_LEN_MIN, target_branch_min),
        proximal,
        distal);
    internal.branch_length_to_parent = static_cast<fp_t>(proximal);
    target_child.branch_length_to_parent = static_cast<fp_t>(distal);
    target_child.parent = free_internal_id;

    TreeNode& pruned = tree.nodes[pruned_id];
    pruned.parent = free_internal_id;
    pruned.branch_length_to_parent = static_cast<fp_t>(sanitize_branch_length(
        pendant_length,
        pruned_branch_min));
    rebuild_traversals(tree);
}
