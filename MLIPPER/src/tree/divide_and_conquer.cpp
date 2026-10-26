#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "divide_and_conquer.hpp"
#include "tree_topology_utils.hpp"
#include "util/checked_size.hpp"

namespace mlipper {
namespace divide_and_conquer {
namespace {

const std::unordered_map<std::string, size_t>& cached_msa_row_index(
    const std::vector<std::string>& names,
    const char* context)
{
    // Session-owned alignment names are append-only. Cache identity and the
    // indexed prefix so repeated sector extraction only indexes newly appended
    // query names instead of rebuilding the full map.
    struct Cache {
        const std::vector<std::string>* owner = nullptr;
        const std::string* data = nullptr;
        size_t indexed_size = 0;
        std::string first_name;
        std::string last_name;
        std::unordered_map<std::string, size_t> row_by_name;
    };
    thread_local Cache cache;
    const bool identity_changed =
        cache.owner != &names ||
        cache.data != names.data() ||
        names.size() < cache.indexed_size ||
        (cache.indexed_size > 0 &&
         (!names.empty() &&
          (names.front() != cache.first_name ||
           names[cache.indexed_size - 1] != cache.last_name)));
    if (identity_changed) {
        cache = Cache{};
        cache.owner = &names;
        cache.data = names.data();
        cache.row_by_name.reserve(names.size() * 2);
    }
    for (size_t row = cache.indexed_size; row < names.size(); ++row) {
        const auto [_, inserted] =
            cache.row_by_name.emplace(names[row], row);
        if (!inserted) {
            throw std::runtime_error(
                std::string(context) + ": duplicate alignment tip name.");
        }
    }
    cache.indexed_size = names.size();
    cache.data = names.data();
    cache.first_name = names.empty() ? std::string{} : names.front();
    cache.last_name = names.empty() ? std::string{} : names.back();
    return cache.row_by_name;
}

namespace neighborhood {

void validate_node_id(
    const TreeBuildResult& tree,
    int node_id,
    const char* context)
{
    if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) {
        throw std::invalid_argument(
            std::string(context) + ": node id is out of range.");
    }
}

std::vector<int> build_tree_path_nodes(
    const TreeBuildResult& tree,
    int start,
    int end)
{
    validate_node_id(tree, start, "build_tree_path_nodes");
    validate_node_id(tree, end, "build_tree_path_nodes");

    std::vector<int> path_a;
    std::unordered_map<int, int> pos_a;
    int current_id = start;
    while (current_id >= 0) {
        pos_a[current_id] = static_cast<int>(path_a.size());
        path_a.push_back(current_id);
        current_id = tree.nodes[static_cast<size_t>(current_id)].parent;
    }

    std::vector<int> path_b;
    int lca = -1;
    current_id = end;
    while (current_id >= 0) {
        const auto it = pos_a.find(current_id);
        if (it != pos_a.end()) {
            lca = current_id;
            break;
        }
        path_b.push_back(current_id);
        current_id = tree.nodes[static_cast<size_t>(current_id)].parent;
    }

    if (lca < 0) {
        throw std::runtime_error(
            "build_tree_path_nodes: failed to resolve an LCA for the requested endpoints.");
    }

    std::vector<int> out;
    const int lca_pos = pos_a[lca];
    out.insert(out.end(), path_a.begin(), path_a.begin() + lca_pos + 1);
    std::reverse(path_b.begin(), path_b.end());
    out.insert(out.end(), path_b.begin(), path_b.end());
    return out;
}

void try_enqueue_bfs_node(
    int node_id,
    const TreeBuildResult& tree,
    std::vector<char>& visited,
    std::vector<int>& bfs_queue,
    std::vector<int>& predecessor,
    int from_node_id)
{
    if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) {
        return;
    }
    if (visited[static_cast<size_t>(node_id)] != 0) {
        return;
    }
    visited[static_cast<size_t>(node_id)] = 1;
    predecessor[static_cast<size_t>(node_id)] = from_node_id;
    bfs_queue.push_back(node_id);
}

std::vector<int> multi_source_bfs_distances(
    const TreeBuildResult& tree,
    const std::vector<char>& source_mask)
{
    if (source_mask.size() != tree.nodes.size()) {
        throw std::invalid_argument(
            "multi_source_bfs_distances: source_mask size must match tree.nodes.size().");
    }

    std::vector<int> distances(tree.nodes.size(), -1);
    std::vector<int> queue;
    queue.reserve(tree.nodes.size());
    for (int node_id = 0;
         node_id < static_cast<int>(tree.nodes.size());
         ++node_id) {
        if (source_mask[static_cast<size_t>(node_id)] == 0) {
            continue;
        }
        distances[static_cast<size_t>(node_id)] = 0;
        queue.push_back(node_id);
    }

    for (size_t queue_idx = 0; queue_idx < queue.size(); ++queue_idx) {
        const int node_id = queue[queue_idx];
        const int next_distance = distances[static_cast<size_t>(node_id)] + 1;
        const TreeNode& node = tree.nodes[static_cast<size_t>(node_id)];
        const int neighbors[3] = {
            node.parent,
            node.left,
            node.right,
        };
        for (int neighbor_id : neighbors) {
            if (neighbor_id < 0 ||
                neighbor_id >= static_cast<int>(tree.nodes.size())) {
                continue;
            }
            if (distances[static_cast<size_t>(neighbor_id)] >= 0) {
                continue;
            }
            distances[static_cast<size_t>(neighbor_id)] = next_distance;
            queue.push_back(neighbor_id);
        }
    }
    return distances;
}

std::vector<char> build_skeleton_mask_from_anchors(
    const TreeBuildResult& tree,
    const std::vector<TreeEdgeEndpoints>& anchors)
{
    std::vector<char> skeleton_mask(tree.nodes.size(), 0);
    for (const TreeEdgeEndpoints& anchor : anchors) {
        const std::vector<int> path_nodes =
            build_tree_path_nodes(tree, anchor.endpoint_a, anchor.endpoint_b);
        for (int node_id : path_nodes) {
            skeleton_mask[static_cast<size_t>(node_id)] = 1;
        }
    }
    return skeleton_mask;
}

void connect_skeleton_components(
    const TreeBuildResult& tree,
    std::vector<char>& skeleton_mask)
{
    std::vector<char> visited(tree.nodes.size(), 0);
    std::vector<int> component_representatives;
    component_representatives.reserve(tree.nodes.size());

    for (int start_id = 0; start_id < static_cast<int>(tree.nodes.size()); ++start_id) {
        if (skeleton_mask[static_cast<size_t>(start_id)] == 0 ||
            visited[static_cast<size_t>(start_id)] != 0) {
            continue;
        }

        component_representatives.push_back(start_id);
        std::vector<int> queue{start_id};
        visited[static_cast<size_t>(start_id)] = 1;
        for (size_t queue_idx = 0; queue_idx < queue.size(); ++queue_idx) {
            const int node_id = queue[queue_idx];
            const TreeNode& node = tree.nodes[static_cast<size_t>(node_id)];
            const int neighbors[3] = {
                node.parent,
                node.left,
                node.right,
            };
            for (int neighbor_id : neighbors) {
                if (neighbor_id < 0 ||
                    neighbor_id >= static_cast<int>(tree.nodes.size()) ||
                    skeleton_mask[static_cast<size_t>(neighbor_id)] == 0 ||
                    visited[static_cast<size_t>(neighbor_id)] != 0) {
                    continue;
                }
                visited[static_cast<size_t>(neighbor_id)] = 1;
                queue.push_back(neighbor_id);
            }
        }
    }

    if (component_representatives.size() <= 1) {
        return;
    }

    const int hub_id = component_representatives.front();
    for (size_t comp_idx = 1;
         comp_idx < component_representatives.size();
         ++comp_idx) {
        const std::vector<int> connector_path =
            build_tree_path_nodes(
                tree,
                hub_id,
                component_representatives[comp_idx]);
        for (int node_id : connector_path) {
            skeleton_mask[static_cast<size_t>(node_id)] = 1;
        }
    }
}

std::vector<int> build_subtree_node_ids_from_selected_tips(
    const TreeBuildResult& tree,
    const std::vector<char>& skeleton_mask,
    const std::vector<int>& predecessor,
    const std::vector<int>& selected_tip_ids)
{
    std::vector<char> subtree_mask = skeleton_mask;
    for (int tip_id : selected_tip_ids) {
        validate_node_id(
            tree,
            tip_id,
            "build_subtree_node_ids_from_selected_tips");
        int current_id = tip_id;
        while (current_id >= 0) {
            subtree_mask[static_cast<size_t>(current_id)] = 1;
            if (skeleton_mask[static_cast<size_t>(current_id)] != 0) {
                break;
            }
            current_id = predecessor[static_cast<size_t>(current_id)];
        }
        if (current_id < 0) {
            throw std::runtime_error(
                "build_subtree_node_ids_from_selected_tips: selected tip "
                "does not connect to the anchor skeleton through BFS "
                "predecessors");
        }
    }

    std::vector<int> subtree_node_ids;
    subtree_node_ids.reserve(tree.nodes.size());
    for (int node_id = 0;
         node_id < static_cast<int>(tree.nodes.size());
         ++node_id) {
        if (subtree_mask[static_cast<size_t>(node_id)] != 0) {
            subtree_node_ids.push_back(node_id);
        }
    }
    return subtree_node_ids;
}

} // namespace neighborhood

} // namespace

std::vector<int> build_tree_path_nodes(
    const TreeBuildResult& tree,
    int start,
    int end)
{
    return neighborhood::build_tree_path_nodes(tree, start, end);
}

std::vector<int> multi_source_bfs_distances(
    const TreeBuildResult& tree,
    const std::vector<char>& source_mask)
{
    return neighborhood::multi_source_bfs_distances(tree, source_mask);
}

// Planning keeps global IDs until the complete connected region is known. The
// dense path favors small trees; the sparse path avoids repeatedly scanning a
// whole large backbone for a small sector. Both preserve the same tip-budget
// and complete-BFS-level semantics.
SubtreePlan build_subtree_plan_from_anchors(
    const TreeBuildResult& tree,
    const std::vector<TreeEdgeEndpoints>& anchors,
    int tip_budget)
{
    if (tree.nodes.empty()) {
        throw std::invalid_argument(
            "build_subtree_plan_from_anchors: tree must not be empty.");
    }
    if (anchors.empty()) {
        throw std::invalid_argument(
            "build_subtree_plan_from_anchors: anchors must not be empty.");
    }
    if (tip_budget <= 0) {
        throw std::invalid_argument(
            "build_subtree_plan_from_anchors: tip_budget must be > 0.");
    }
    // Dense masks are fastest for small trees, but scanning and clearing a
    // 400k-node array for every 200-tip extraction dominates large-backbone
    // D&C.  The sparse path visits only the anchor skeleton and the local BFS
    // frontier, and intentionally leaves global_to_local empty so subtree
    // construction can use its sparse index as well.
    if (tree.nodes.size() > 10000) {
        const int skeleton_origin = anchors.front().endpoint_a;
        neighborhood::validate_node_id(
            tree, skeleton_origin, "build_subtree_plan_from_anchors");
        std::unordered_set<int> skeleton_set;
        skeleton_set.reserve(anchors.size() * 16 + 16);
        for (const TreeEdgeEndpoints& anchor : anchors) {
            neighborhood::validate_node_id(
                tree, anchor.endpoint_a, "build_subtree_plan_from_anchors");
            neighborhood::validate_node_id(
                tree, anchor.endpoint_b, "build_subtree_plan_from_anchors");
            const TreeNode& endpoint_a =
                tree.nodes[static_cast<size_t>(anchor.endpoint_a)];
            const TreeNode& endpoint_b =
                tree.nodes[static_cast<size_t>(anchor.endpoint_b)];
            if (endpoint_a.parent != anchor.endpoint_b &&
                endpoint_b.parent != anchor.endpoint_a) {
                throw std::invalid_argument(
                    "build_subtree_plan_from_anchors: anchor endpoints must be adjacent.");
            }
            for (int endpoint : {anchor.endpoint_a, anchor.endpoint_b}) {
                const std::vector<int> path =
                    neighborhood::build_tree_path_nodes(
                        tree, skeleton_origin, endpoint);
                skeleton_set.insert(path.begin(), path.end());
            }
        }

        SubtreePlan sparse_plan;
        sparse_plan.skeleton_nodes.assign(
            skeleton_set.begin(), skeleton_set.end());
        std::sort(
            sparse_plan.skeleton_nodes.begin(),
            sparse_plan.skeleton_nodes.end());
        std::vector<int> bfs_queue = sparse_plan.skeleton_nodes;
        std::unordered_map<int, int> predecessor;
        predecessor.reserve(
            static_cast<size_t>(tip_budget) * 8 +
            sparse_plan.skeleton_nodes.size() * 2);
        for (int node_id : sparse_plan.skeleton_nodes) {
            predecessor.emplace(node_id, -1);
        }
        auto enqueue = [&](int node_id, int from_node_id) {
            if (node_id < 0 ||
                node_id >= static_cast<int>(tree.nodes.size()) ||
                predecessor.count(node_id) != 0) {
                return;
            }
            predecessor.emplace(node_id, from_node_id);
            bfs_queue.push_back(node_id);
        };

        size_t level_begin = 0;
        while (level_begin < bfs_queue.size() &&
               sparse_plan.selected_tip_ids.size() <
                   static_cast<size_t>(tip_budget)) {
            const size_t level_end = bfs_queue.size();
            std::vector<int> level_tip_ids;
            for (size_t queue_idx = level_begin;
                 queue_idx < level_end;
                 ++queue_idx) {
                const int node_id = bfs_queue[queue_idx];
                const TreeNode& node =
                    tree.nodes[static_cast<size_t>(node_id)];
                if (node.is_tip) {
                    level_tip_ids.push_back(node_id);
                }
                enqueue(node.parent, node_id);
                enqueue(node.left, node_id);
                enqueue(node.right, node_id);
            }
            std::sort(level_tip_ids.begin(), level_tip_ids.end());
            const size_t remaining_budget =
                static_cast<size_t>(tip_budget) -
                sparse_plan.selected_tip_ids.size();
            if (level_tip_ids.size() > remaining_budget) {
                level_tip_ids.resize(remaining_budget);
            }
            sparse_plan.selected_tip_ids.insert(
                sparse_plan.selected_tip_ids.end(),
                level_tip_ids.begin(),
                level_tip_ids.end());
            level_begin = level_end;
        }
        if (sparse_plan.selected_tip_ids.empty()) {
            throw std::runtime_error(
                "build_subtree_plan_from_anchors: failed to collect any tip nodes.");
        }

        std::unordered_set<int> subtree_set = skeleton_set;
        subtree_set.reserve(
            skeleton_set.size() +
            sparse_plan.selected_tip_ids.size() * 8);
        for (int tip_id : sparse_plan.selected_tip_ids) {
            int node_id = tip_id;
            while (node_id >= 0 && subtree_set.insert(node_id).second) {
                const auto predecessor_it = predecessor.find(node_id);
                if (predecessor_it == predecessor.end()) {
                    throw std::runtime_error(
                        "build_subtree_plan_from_anchors: sparse predecessor chain is incomplete.");
                }
                node_id = predecessor_it->second;
            }
        }
        sparse_plan.subtree_node_ids.assign(
            subtree_set.begin(), subtree_set.end());
        std::sort(
            sparse_plan.subtree_node_ids.begin(),
            sparse_plan.subtree_node_ids.end());
        sparse_plan.local_to_global = sparse_plan.subtree_node_ids;
        return sparse_plan;
    }

    SubtreePlan plan;
    std::vector<char> skeleton_mask =
        neighborhood::build_skeleton_mask_from_anchors(
            tree,
            anchors);
    neighborhood::connect_skeleton_components(
        tree,
        skeleton_mask);
    plan.skeleton_nodes.reserve(tree.nodes.size());
    for (int node_id = 0;
         node_id < static_cast<int>(tree.nodes.size());
         ++node_id) {
        if (skeleton_mask[static_cast<size_t>(node_id)] != 0) {
            plan.skeleton_nodes.push_back(node_id);
        }
    }
    if (plan.skeleton_nodes.empty()) {
        throw std::runtime_error(
            "build_subtree_plan_from_anchors: anchor skeleton is empty.");
    }

    std::vector<char> visited(tree.nodes.size(), 0);
    std::vector<int> predecessor(tree.nodes.size(), -1);
    std::vector<int> bfs_queue;
    bfs_queue.reserve(tree.nodes.size());
    for (int node_id : plan.skeleton_nodes) {
        neighborhood::try_enqueue_bfs_node(
            node_id,
            tree,
            visited,
            bfs_queue,
            predecessor,
            -1);
    }

    size_t level_begin = 0;
    while (level_begin < bfs_queue.size() && plan.selected_tip_ids.size() < static_cast<size_t>(tip_budget)) {
        const size_t level_end = bfs_queue.size();
        std::vector<int> level_tip_ids;
        for (size_t queue_idx = level_begin;
             queue_idx < level_end;
             ++queue_idx) {
            const int node_id = bfs_queue[queue_idx];
            const TreeNode& node = tree.nodes[static_cast<size_t>(node_id)];
            if (node.is_tip) {
                level_tip_ids.push_back(node_id);
            }

            const int neighbors[3] = {
                node.parent,
                node.left,
                node.right,
            };
            for (int neighbor_id : neighbors) {
                neighborhood::try_enqueue_bfs_node(
                    neighbor_id,
                    tree,
                    visited,
                    bfs_queue,
                    predecessor,
                    node_id);
            }
        }

        std::sort(level_tip_ids.begin(), level_tip_ids.end());
        const size_t remaining_budget =
            static_cast<size_t>(tip_budget) -
            plan.selected_tip_ids.size();
        if (level_tip_ids.size() > remaining_budget) {
            level_tip_ids.resize(remaining_budget);
        }
        plan.selected_tip_ids.insert(
            plan.selected_tip_ids.end(),
            level_tip_ids.begin(),
            level_tip_ids.end());
        level_begin = level_end;
    }

    if (plan.selected_tip_ids.empty()) {
        throw std::runtime_error(
            "build_subtree_plan_from_anchors: failed to collect any tip nodes.");
    }

    plan.subtree_node_ids =
        neighborhood::build_subtree_node_ids_from_selected_tips(
            tree,
            skeleton_mask,
            predecessor,
            plan.selected_tip_ids);
    plan.local_to_global = plan.subtree_node_ids;
    plan.global_to_local.assign(tree.nodes.size(), -1);
    for (size_t local_id = 0;
         local_id < plan.local_to_global.size();
         ++local_id) {
        const int global_id = plan.local_to_global[local_id];
        neighborhood::validate_node_id(
            tree,
            global_id,
            "build_subtree_plan_from_anchors");
        plan.global_to_local[static_cast<size_t>(global_id)] =
            static_cast<int>(local_id);
    }
    return plan;
}

SubtreePlan build_nni_sector_plan(
    const TreeBuildResult& tree,
    const std::vector<TreeEdgeEndpoints>& central_edges)
{
    if (tree.nodes.empty()) {
        throw std::invalid_argument(
            "build_nni_sector_plan: tree must not be empty.");
    }
    if (central_edges.empty()) {
        throw std::invalid_argument(
            "build_nni_sector_plan: central_edges must not be empty.");
    }

    for (const TreeEdgeEndpoints& edge : central_edges) {
        neighborhood::validate_node_id(
            tree, edge.endpoint_a, "build_nni_sector_plan");
        neighborhood::validate_node_id(
            tree, edge.endpoint_b, "build_nni_sector_plan");
        const TreeNode& a = tree.nodes[static_cast<size_t>(edge.endpoint_a)];
        const TreeNode& b = tree.nodes[static_cast<size_t>(edge.endpoint_b)];
        if (a.parent != b.id && b.parent != a.id) {
            throw std::invalid_argument(
                "build_nni_sector_plan: central-edge endpoints must be adjacent.");
        }
        if (a.is_tip || b.is_tip) {
            throw std::invalid_argument(
                "build_nni_sector_plan: a central NNI edge must connect two internal nodes.");
        }
    }

    std::vector<char> skeleton_mask =
        neighborhood::build_skeleton_mask_from_anchors(tree, central_edges);
    neighborhood::connect_skeleton_components(tree, skeleton_mask);

    SubtreePlan plan;
    plan.skeleton_nodes.reserve(central_edges.size() * 2);
    for (int node_id = 0;
         node_id < static_cast<int>(tree.nodes.size());
         ++node_id) {
        if (skeleton_mask[static_cast<size_t>(node_id)] == 0) {
            continue;
        }
        plan.skeleton_nodes.push_back(node_id);
        plan.subtree_node_ids.push_back(node_id);
    }
    if (plan.subtree_node_ids.empty()) {
        throw std::runtime_error(
            "build_nni_sector_plan: central-edge skeleton is empty.");
    }
    plan.local_to_global = plan.subtree_node_ids;
    // Keep this index sparse. Compact NNI sectors are intentionally much
    // smaller than the full tree, and build_subtree_from_plan() already
    // constructs the required local lookup from local_to_global.
    return plan;
}

TreeBuildResult build_subtree_from_plan(
    const TreeBuildResult& tree,
    const SubtreePlan& plan)
{
    if (plan.subtree_node_ids.empty()) {
        throw std::invalid_argument(
            "build_subtree_from_plan: subtree_node_ids must not be empty.");
    }

    // Plans may carry a full global-to-local table or only local_to_global for
    // large sparse sectors. Materialization normalizes either form to compact
    // local node IDs without changing the full-tree IDs stored in the plan.
    const bool dense_index =
        plan.global_to_local.size() == tree.nodes.size();
    std::unordered_map<int, int> sparse_old_to_new;
    if (!dense_index) {
        sparse_old_to_new.reserve(plan.subtree_node_ids.size() * 2);
        for (size_t new_id = 0;
             new_id < plan.subtree_node_ids.size();
             ++new_id) {
            const int old_id = plan.subtree_node_ids[new_id];
            if (old_id < 0 || old_id >= static_cast<int>(tree.nodes.size())) {
                throw std::runtime_error(
                    "build_subtree_from_plan: subtree node id is out of range.");
            }
            sparse_old_to_new.emplace(old_id, static_cast<int>(new_id));
        }
    }
    auto local_id = [&](int old_id) {
        if (old_id < 0) {
            return -1;
        }
        if (dense_index) {
            return plan.global_to_local[static_cast<size_t>(old_id)];
        }
        const auto found = sparse_old_to_new.find(old_id);
        return found == sparse_old_to_new.end() ? -1 : found->second;
    };

    TreeBuildResult subtree;
    subtree.nodes.reserve(plan.subtree_node_ids.size());
    for (int old_id : plan.subtree_node_ids) {
        if (old_id < 0 || old_id >= static_cast<int>(tree.nodes.size())) {
            throw std::runtime_error(
                "build_subtree_from_plan: subtree node id is out of range.");
        }
        TreeNode new_node = tree.nodes[static_cast<size_t>(old_id)];
        new_node.id = local_id(old_id);
        subtree.nodes.push_back(std::move(new_node));
    }

    int root_count = 0;
    subtree.tip_node_by_name.clear();
    subtree.tip_node_by_name.reserve(subtree.nodes.size());
    for (TreeNode& node : subtree.nodes) {
        node.parent = local_id(node.parent);
        node.left = local_id(node.left);
        node.right = local_id(node.right);

        if (node.parent < 0) {
            subtree.root_id = node.id;
            ++root_count;
        }
        if (node.is_tip) {
            if (node.name.empty()) {
                throw std::runtime_error(
                    "build_subtree_from_plan: subtree tip is missing a name.");
            }
            const auto [_, inserted] =
                subtree.tip_node_by_name.emplace(node.name, node.id);
            if (!inserted) {
                throw std::runtime_error(
                    "build_subtree_from_plan: duplicate tip name in subtree.");
            }
        }
    }

    if (root_count != 1) {
        throw std::runtime_error(
            "build_subtree_from_plan: subtree must have exactly one root.");
    }

    rebuild_traversals(subtree);
    return subtree;
}

enum class DirectionalBoundaryMode {
    PreserveRealTips,
    MovableNNIComponents,
};

static int append_directional_boundary(
    PreparedDirectionalSubtree& prepared,
    std::vector<int>& local_to_global,
    int& boundary_serial,
    int parent_local_id,
    fp_t branch_length,
    int outside_global_id,
    DirectionalBoundaryMode mode)
{
    TreeNode boundary;
    boundary.id = static_cast<int>(prepared.subtree.nodes.size());
    boundary.stable_node_label = -1;
    boundary.is_tip = true;
    boundary.parent = parent_local_id;
    boundary.branch_length_to_parent = branch_length;
    boundary.name = "__mlipper_directional_boundary_" +
        std::to_string(boundary_serial++) + "__";
    prepared.subtree.nodes.push_back(boundary);
    local_to_global.push_back(
        mode == DirectionalBoundaryMode::MovableNNIComponents
            ? outside_global_id
            : -1);
    prepared.subtree.tip_node_by_name.emplace(boundary.name, boundary.id);
    prepared.boundary_ports.push_back(DirectionalBoundaryPort{
        boundary.id,
        outside_global_id,
    });
    return boundary.id;
}

// Convert a connected global plan into a self-contained local likelihood tree.
// Every cut edge gains a virtual tip. Its sequence is structural padding while
// boundary_ports instruct the GPU setup to replace that tip contribution with
// the directional CLV copied from the omitted full-tree component.
static PreparedDirectionalSubtree prepare_directional_subtree_with_boundaries(
    const TreeBuildResult& tree,
    const SubtreePlan& plan,
    DirectionalBoundaryMode boundary_mode,
    const std::vector<std::string>& full_msa_tip_names,
    const std::vector<std::string>& full_msa_rows,
    size_t sites,
    int states,
    const EigResult& eig,
    const std::vector<double>& rate_multipliers,
    int rate_cats)
{
    if (full_msa_tip_names.size() != full_msa_rows.size()) {
        throw std::invalid_argument(
            "prepare_directional_subtree_with_boundaries: alignment size mismatch.");
    }

    PreparedDirectionalSubtree out;
    out.subtree = build_subtree_from_plan(tree, plan);
    std::vector<int> local_to_global = plan.local_to_global;

    int boundary_serial = 0;
    // Complete every selected internal node with an exact message for each
    // omitted child component.  The original child branch remains explicit,
    // so the installed upward CLV is transformed by the normal local PMAT.
    const size_t selected_node_count = plan.local_to_global.size();
    for (size_t local_id = 0; local_id < selected_node_count; ++local_id) {
        if (out.subtree.nodes[local_id].is_tip) {
            continue;
        }
        const int global_id = plan.local_to_global[local_id];
        const TreeNode& global_node = tree.nodes[static_cast<size_t>(global_id)];
        const int global_children[2] = {global_node.left, global_node.right};
        for (int slot = 0; slot < 2; ++slot) {
            const int global_child = global_children[slot];
            const int local_child = slot == 0
                ? out.subtree.nodes[local_id].left
                : out.subtree.nodes[local_id].right;
            if (global_child < 0 || local_child >= 0) {
                continue;
            }
            if (boundary_mode == DirectionalBoundaryMode::PreserveRealTips &&
                tree.nodes[static_cast<size_t>(global_child)].is_tip) {
                TreeNode local_tip =
                    tree.nodes[static_cast<size_t>(global_child)];
                local_tip.id =
                    static_cast<int>(out.subtree.nodes.size());
                local_tip.parent = static_cast<int>(local_id);
                out.subtree.nodes.push_back(local_tip);
                out.subtree.tip_node_by_name.emplace(
                    local_tip.name, local_tip.id);
                local_to_global.push_back(global_child);
                if (slot == 0) {
                    out.subtree.nodes[local_id].left = local_tip.id;
                } else {
                    out.subtree.nodes[local_id].right = local_tip.id;
                }
                continue;
            }
            const int boundary_id = append_directional_boundary(
                out,
                local_to_global,
                boundary_serial,
                static_cast<int>(local_id),
                tree.nodes[static_cast<size_t>(global_child)]
                    .branch_length_to_parent,
                global_child,
                boundary_mode);
            if (slot == 0) {
                out.subtree.nodes[local_id].left = boundary_id;
            } else {
                out.subtree.nodes[local_id].right = boundary_id;
            }
        }
    }

    // Keep the original full-tree rooting by connecting the selected
    // component to the global root through a narrow ancestor spine.  Every
    // omitted sibling component is represented by an upward message.  This
    // avoids treating a parent-side directional value as an ordinary child
    // CLV at a synthetic root, and only adds O(tree depth) explicit nodes.
    int spine_local_child = out.subtree.root_id;
    int spine_global_child =
        local_to_global[static_cast<size_t>(spine_local_child)];
    int spine_global_parent =
        tree.nodes[static_cast<size_t>(spine_global_child)].parent;
    while (spine_global_parent >= 0) {
        const TreeNode& global_parent =
            tree.nodes[static_cast<size_t>(spine_global_parent)];
        const bool child_is_left =
            global_parent.left == spine_global_child;
        if (!child_is_left && global_parent.right != spine_global_child) {
            throw std::runtime_error(
                "prepare_directional_subtree_with_boundaries: ancestor spine is inconsistent.");
        }
        const int sibling_global =
            child_is_left ? global_parent.right : global_parent.left;
        if (sibling_global < 0) {
            throw std::runtime_error(
                "prepare_directional_subtree_with_boundaries: ancestor is missing a sibling child.");
        }

        TreeNode local_parent = global_parent;
        local_parent.id =
            static_cast<int>(out.subtree.nodes.size());
        local_parent.parent = -1;
        local_parent.left = -1;
        local_parent.right = -1;
        out.subtree.nodes.push_back(local_parent);
        local_to_global.push_back(spine_global_parent);

        int sibling_local = -1;
        if (boundary_mode == DirectionalBoundaryMode::PreserveRealTips &&
            tree.nodes[static_cast<size_t>(sibling_global)].is_tip) {
            TreeNode local_tip =
                tree.nodes[static_cast<size_t>(sibling_global)];
            local_tip.id =
                static_cast<int>(out.subtree.nodes.size());
            local_tip.parent = local_parent.id;
            out.subtree.nodes.push_back(local_tip);
            local_to_global.push_back(sibling_global);
            out.subtree.tip_node_by_name.emplace(
                local_tip.name, local_tip.id);
            sibling_local = local_tip.id;
        } else {
            sibling_local = append_directional_boundary(
                out,
                local_to_global,
                boundary_serial,
                local_parent.id,
                tree.nodes[static_cast<size_t>(sibling_global)]
                    .branch_length_to_parent,
                sibling_global,
                boundary_mode);
        }

        out.subtree.nodes[static_cast<size_t>(spine_local_child)].parent =
            local_parent.id;
        if (child_is_left) {
            out.subtree.nodes[static_cast<size_t>(local_parent.id)].left =
                spine_local_child;
            out.subtree.nodes[static_cast<size_t>(local_parent.id)].right =
                sibling_local;
        } else {
            out.subtree.nodes[static_cast<size_t>(local_parent.id)].left =
                sibling_local;
            out.subtree.nodes[static_cast<size_t>(local_parent.id)].right =
                spine_local_child;
        }
        spine_local_child = local_parent.id;
        spine_global_child = spine_global_parent;
        spine_global_parent = global_parent.parent;
    }
    out.subtree.root_id = spine_local_child;
    rebuild_traversals(out.subtree);

    // Directional boundary completion makes every retained internal node a
    // proper binary node, while preserving every retained full-tree edge as
    // an independently scored placement candidate.
    for (const TreeNode& node : out.subtree.nodes) {
        if (node.is_tip) {
            continue;
        }
        const int child_count =
            static_cast<int>(node.left >= 0) +
            static_cast<int>(node.right >= 0);
        if (child_count != 2) {
            throw std::runtime_error(
                "prepare_directional_subtree_with_boundaries: "
                "directional completion left a non-binary internal node.");
        }
    }
    out.local_to_global = std::move(local_to_global);
    TreeBuildResult& subtree = out.subtree;
    const auto& full_row_by_name = cached_msa_row_index(
        full_msa_tip_names,
        "prepare_directional_subtree_with_boundaries");
    for (int node_id : subtree.postorder) {
        const TreeNode& node = subtree.nodes[static_cast<size_t>(node_id)];
        if (!node.is_tip) {
            continue;
        }
        out.alignment.names.push_back(node.name);
        const auto row_it = full_row_by_name.find(node.name);
        if (row_it == full_row_by_name.end()) {
            out.alignment.sequences.emplace_back(sites, 'N');
        } else {
            out.alignment.sequences.push_back(
                full_msa_rows[row_it->second]);
        }
    }
    out.alignment.sites = sites;
    out.host_packing = pack_host_arrays_from_tree_and_msa(
        subtree,
        out.alignment.names,
        out.alignment.sequences,
        sites,
        states);
    out.host_packing.pattern_weights.assign(sites, 1u);

    // Keep the allocated tip storage, but make traversal kernels treat every
    // boundary as an already-computed internal CLV.
    for (const DirectionalBoundaryPort& port : out.boundary_ports) {
        subtree.nodes[static_cast<size_t>(port.local_node_id)].is_tip = false;
        out.host_packing.is_tip[static_cast<size_t>(port.local_node_id)] = 0;
    }
    rebuild_traversals(subtree);
    out.host_packing.preorder = subtree.preorder;
    out.host_packing.postorder = subtree.postorder;
    fill_pmats_in_host_packing(
        subtree,
        out.host_packing,
        eig,
        rate_multipliers,
        states,
        rate_cats);
    return out;
}

PreparedDirectionalSubtree prepare_anchor_subtree_with_directional_boundaries(
    const TreeBuildResult& tree,
    const std::vector<TreeEdgeEndpoints>& anchors,
    int tip_budget,
    const std::vector<std::string>& full_msa_tip_names,
    const std::vector<std::string>& full_msa_rows,
    size_t sites,
    int states,
    const EigResult& eig,
    const std::vector<double>& rate_multipliers,
    int rate_cats)
{
    const SubtreePlan plan =
        build_subtree_plan_from_anchors(tree, anchors, tip_budget);
    return prepare_directional_subtree_with_boundaries(
        tree,
        plan,
        DirectionalBoundaryMode::PreserveRealTips,
        full_msa_tip_names,
        full_msa_rows,
        sites,
        states,
        eig,
        rate_multipliers,
        rate_cats);
}

PreparedDirectionalSubtree prepare_nni_sector_with_directional_boundaries(
    const TreeBuildResult& tree,
    const std::vector<TreeEdgeEndpoints>& central_edges,
    const std::vector<std::string>& full_msa_tip_names,
    const std::vector<std::string>& full_msa_rows,
    size_t sites,
    int states,
    const EigResult& eig,
    const std::vector<double>& rate_multipliers,
    int rate_cats)
{
    const SubtreePlan plan = build_nni_sector_plan(tree, central_edges);
    return prepare_directional_subtree_with_boundaries(
        tree,
        plan,
        DirectionalBoundaryMode::MovableNNIComponents,
        full_msa_tip_names,
        full_msa_rows,
        sites,
        states,
        eig,
        rate_multipliers,
        rate_cats);
}

PlacementQueryBatch build_query_batch_from_placement_queries(
    const std::vector<mlipper::SequenceRecord>& placement_queries,
    size_t sites,
    int states)
{
    if (states != 4 && states != 5) {
        throw std::invalid_argument(
            "build_query_batch_from_placement_queries supports 4- or 5-state encoding");
    }
    PlacementQueryBatch batch;
    batch.count = placement_queries.size();
    if (batch.empty()) {
        return batch;
    }
    batch.query_chars.resize(
        mlipper::util::checked_mul_size(
            batch.count, sites, "build_query_batch_from_placement_queries"),
        states == 4 ? static_cast<uint8_t>(15) : static_cast<uint8_t>(4));

    for (size_t query_idx = 0; query_idx < placement_queries.size(); ++query_idx) {
        const mlipper::SequenceRecord& query = placement_queries[query_idx];
        if (query.sequence.size() != sites) {
            throw std::runtime_error(
                "build_query_batch_from_placement_queries: query sequence length mismatch (sequence_sites=" +
                std::to_string(query.sequence.size()) + ", tree_sites=" +
                std::to_string(sites) + ").");
        }
        for (size_t site_idx = 0; site_idx < sites; ++site_idx) {
            batch.query_chars[query_idx * sites + site_idx] =
                (states == 4)
                    ? encode_state_DNA4_mask(query.sequence[site_idx])
                    : encode_state_DNA5(query.sequence[site_idx]);
        }
    }

    return batch;
}

} // namespace divide_and_conquer
} // namespace mlipper
