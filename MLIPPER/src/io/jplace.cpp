#include "io/jplace.hpp"
#include "io/tree_newick.hpp"
#include "placement/placement.cuh"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct JplacePlacementRow {
    int edge_num = -1;
    double likelihood = 0.0;
    double like_weight_ratio = 1.0;
    double distal_length = 0.0;
    double pendant_length = 0.0;
};

struct JplacePlacementRecord {
    std::string query_name;
    std::vector<JplacePlacementRow> rows;
};

struct JplaceTreeExport {
    std::string tree;
    std::vector<int> edge_num_by_node;
};

bool has_artificial_root_edge(
    const TreeBuildResult& tree,
    int child_id)
{
    if (child_id < 0 || child_id >= static_cast<int>(tree.nodes.size())) {
        return false;
    }
    const TreeNode& child = tree.nodes[child_id];
    return !child.is_tip &&
        child.parent == tree.root_id &&
        std::abs(child.branch_length_to_parent) <= 1e-15;
}

std::string json_escape_string(const std::string& input) {
    std::ostringstream out;
    for (unsigned char ch : input) {
        switch (ch) {
            case '\"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (ch < 0x20) {
                    out << "\\u"
                        << std::hex << std::setw(4) << std::setfill('0')
                        << static_cast<int>(ch)
                        << std::dec << std::setfill(' ');
                } else {
                    out << static_cast<char>(ch);
                }
                break;
        }
    }
    return out.str();
}

void append_jplace_row(
    const TreeBuildResult& tree,
    const JplaceTreeExport& tree_export,
    JplacePlacementRecord& record,
    int target_id,
    double loglikelihood,
    double like_weight_ratio,
    double distal_length,
    double pendant_length)
{
    if (target_id < 0 || target_id >= static_cast<int>(tree.nodes.size())) {
        throw std::runtime_error("Cannot write jplace: invalid placement target id.");
    }
    const TreeNode& target = tree.nodes[target_id];
    if (target.parent < 0) {
        throw std::runtime_error("Cannot write jplace: placement target is the root.");
    }
    const int edge_num = tree_export.edge_num_by_node[target_id];
    if (edge_num < 0) {
        // The zero-length helper edge introduced while rooting an unrooted
        // binary tree is intentionally flattened out of the jplace tree.
        // Placement candidates on that helper edge therefore have no jplace
        // edge number and must not be exported.
        if (has_artificial_root_edge(tree, target_id)) {
            return;
        }
        throw std::runtime_error("Cannot write jplace: target edge has no edge number.");
    }
    if (!std::isfinite(loglikelihood) ||
        !std::isfinite(like_weight_ratio) || like_weight_ratio < 0.0 ||
        !std::isfinite(distal_length) || distal_length < 0.0 ||
        !std::isfinite(pendant_length) || pendant_length < 0.0) {
        throw std::runtime_error(
            "Cannot write jplace: placement row contains invalid numeric values.");
    }

    JplacePlacementRow row;
    row.edge_num = edge_num;
    row.likelihood = loglikelihood;
    row.like_weight_ratio = like_weight_ratio;
    row.distal_length = distal_length;
    row.pendant_length = pendant_length;
    record.rows.push_back(row);
}

// Serialize recursively and assign edge numbers after each child subtree. The
// node-to-edge table is the authoritative mapping later used for placement
// rows; callers must never assume a TreeBuildResult node ID is a jplace edge.
std::string emit_jplace_tree_node(
    const TreeBuildResult& tree,
    JplaceTreeExport& tree_export,
    int& next_edge_num,
    int node_id,
    bool suppress_parent_edge)
{
    if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) {
        throw std::runtime_error("build_jplace_tree: invalid node id.");
    }

    const TreeNode& node = tree.nodes[node_id];
    std::ostringstream out;
    if (node.is_tip) {
        out << mlipper::treeio::format_newick_taxon_name(
            node.name.empty()
                ? ("tip_" + std::to_string(node.id))
                : node.name);
    } else {
        out << "("
            << emit_jplace_tree_node(
                tree, tree_export, next_edge_num, node.left, false)
            << ","
            << emit_jplace_tree_node(
                tree, tree_export, next_edge_num, node.right, false)
            << ")";
    }

    if (node.parent >= 0 && !suppress_parent_edge) {
        const int edge_num = next_edge_num++;
        tree_export.edge_num_by_node[node_id] = edge_num;
        out << ":" << std::setprecision(17) << node.branch_length_to_parent
            << "{" << edge_num << "}";
    }
    return out.str();
}

// Flatten only the zero-length helper edge introduced when an unrooted binary
// tree was represented with a degree-two root. This keeps the exported edge set
// aligned with the biological unrooted tree used for placement interpretation.
JplaceTreeExport build_jplace_tree_export(const TreeBuildResult& tree) {
    mlipper::treeio::validate_tree_for_output(tree);

    JplaceTreeExport result;
    result.edge_num_by_node.assign(tree.nodes.size(), -1);
    int next_edge_num = 0;

    const TreeNode& root = tree.nodes[tree.root_id];
    const bool flatten_left = has_artificial_root_edge(tree, root.left);
    const bool flatten_right = has_artificial_root_edge(tree, root.right);

    std::ostringstream out;
    if (!root.is_tip && flatten_left != flatten_right) {
        const int flat_child = flatten_left ? root.left : root.right;
        const int other_child = (flat_child == root.left) ? root.right : root.left;
        out << "("
            << emit_jplace_tree_node(
                tree, result, next_edge_num, flat_child, true)
            << ","
            << emit_jplace_tree_node(
                tree, result, next_edge_num, other_child, false)
            << ")";
    } else {
        out << emit_jplace_tree_node(
            tree, result, next_edge_num, tree.root_id, false);
    }

    result.tree = out.str() + ";";
    return result;
}

std::vector<JplacePlacementRecord> build_jplace_records(
    const TreeBuildResult& tree,
    const JplaceTreeExport& tree_export,
    const std::vector<::RawPlacementResult>& placement_results,
    const std::vector<std::string>& query_names)
{
    if (placement_results.size() != query_names.size()) {
        throw std::runtime_error(
            "Cannot write jplace: placement result and query name counts differ.");
    }

    std::vector<JplacePlacementRecord> records;
    records.reserve(placement_results.size());

    for (size_t i = 0; i < placement_results.size(); ++i) {
        const ::RawPlacementResult& placement = placement_results[i];
        JplacePlacementRecord record;
        record.query_name = query_names[i].empty()
            ? ("query_" + std::to_string(i))
            : query_names[i];

        if (!placement.top_placements.empty()) {
            for (const ::RawPlacementResult::RankedPlacement& candidate :
                 placement.top_placements) {
                append_jplace_row(
                    tree,
                    tree_export,
                    record,
                    candidate.target_id,
                    candidate.loglikelihood,
                    candidate.like_weight_ratio,
                    candidate.distal_length,
                    candidate.pendant_length);
            }
        } else {
            append_jplace_row(
                tree,
                tree_export,
                record,
                placement.target_id,
                placement.loglikelihood,
                1.0,
                placement.distal_length,
                placement.pendant_length);
        }

        if (record.rows.empty()) {
            throw std::runtime_error("Could not export any placement rows for jplace.");
        }
        double weight_sum = 0.0;
        for (const JplacePlacementRow& row : record.rows) {
            weight_sum += row.like_weight_ratio;
        }
        if (!std::isfinite(weight_sum) || weight_sum <= 0.0) {
            throw std::runtime_error(
                "Cannot write jplace: placement weights have no positive mass.");
        }
        // Filtering the artificial root edge removes probability mass. Restore
        // the per-query jplace invariant after all non-exportable rows are gone.
        for (JplacePlacementRow& row : record.rows) {
            row.like_weight_ratio /= weight_sum;
        }
        records.push_back(std::move(record));
    }

    return records;
}

} // namespace

namespace mlipper {
namespace jplaceio {

void write_jplace(
    const std::string& output_path,
    const TreeBuildResult& tree,
    const std::vector<::RawPlacementResult>& placement_results,
    const std::vector<std::string>& query_names,
    const std::string& invocation)
{
    const JplaceTreeExport tree_export = build_jplace_tree_export(tree);
    const std::vector<JplacePlacementRecord> placements =
        build_jplace_records(
            tree, tree_export, placement_results, query_names);

    const std::filesystem::path path(output_path);
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Cannot open jplace output: " + output_path);
    }

    out << "{\n";
    out << "  \"tree\": \"" << json_escape_string(tree_export.tree) << "\",\n";
    out << "  \"placements\": [\n";
    for (size_t i = 0; i < placements.size(); ++i) {
        const JplacePlacementRecord& rec = placements[i];
        out << "    {\n";
        out << "      \"p\": [";
        for (size_t row_idx = 0; row_idx < rec.rows.size(); ++row_idx) {
            const JplacePlacementRow& row = rec.rows[row_idx];
            if (row_idx == 0) {
                out << "[";
            } else {
                out << ", [";
            }
            out << row.edge_num << ", "
                << std::setprecision(17) << row.likelihood << ", "
                << std::setprecision(17) << row.like_weight_ratio << ", "
                << std::setprecision(17) << row.distal_length << ", "
                << std::setprecision(17) << row.pendant_length << "]";
        }
        out << "],\n";
        out << "      \"n\": [\"" << json_escape_string(rec.query_name) << "\"]\n";
        out << "    }";
        if (i + 1 < placements.size()) out << ",";
        out << "\n";
    }
    out << "  ],\n";
    out << "  \"metadata\": {\n";
    out << "    \"invocation\": \"" << json_escape_string(invocation) << "\"\n";
    out << "  },\n";
    out << "  \"version\": 3,\n";
    out << "  \"fields\": [\"edge_num\", \"likelihood\", "
           "\"like_weight_ratio\", \"distal_length\", \"pendant_length\"]\n";
    out << "}\n";
    if (!out) {
        throw std::runtime_error("Failed while writing jplace output: " + output_path);
    }
}

} // namespace jplaceio
} // namespace mlipper
