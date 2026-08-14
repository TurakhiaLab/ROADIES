#include "io/jplace.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

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
    const mlipper::jplaceio::JplaceTreeExport& tree_export,
    mlipper::jplaceio::JplacePlacementRecord& record,
    int target_id,
    double loglikelihood,
    double like_weight_ratio,
    double proximal_length,
    double pendant_length)
{
    if (target_id < 0 || target_id >= static_cast<int>(tree.nodes.size())) {
        return;
    }
    const TreeNode& target = tree.nodes[target_id];
    if (target.parent < 0) {
        return;
    }
    const int edge_num = tree_export.edge_num_by_node[target_id];
    if (edge_num < 0) {
        return;
    }

    mlipper::jplaceio::JplacePlacementRow row;
    row.edge_num = edge_num;
    row.likelihood = loglikelihood;
    row.like_weight_ratio = like_weight_ratio;
    // PlacementResult stores the jplace distal coordinate here.
    row.distal_length = proximal_length;
    row.pendant_length = pendant_length;
    record.rows.push_back(std::move(row));
}

} // namespace

namespace mlipper {
namespace jplaceio {

JplaceTreeExport build_jplace_tree_export(const TreeBuildResult& tree) {
    if (tree.root_id < 0 || tree.root_id >= static_cast<int>(tree.nodes.size())) {
        throw std::runtime_error("build_jplace_tree: invalid root_id.");
    }

    JplaceTreeExport result;
    result.edge_num_by_node.assign(tree.nodes.size(), -1);
    int next_edge_num = 0;

    std::function<std::string(int, bool)> emit_node = [&](int node_id, bool suppress_parent_edge) -> std::string {
        if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) {
            throw std::runtime_error("build_jplace_tree: invalid node id.");
        }

        const TreeNode& node = tree.nodes[node_id];
        std::ostringstream out;
        if (node.is_tip) {
            out << node.name;
        } else {
            out << "(" << emit_node(node.left, false) << "," << emit_node(node.right, false) << ")";
        }

        if (node.parent >= 0 && !suppress_parent_edge) {
            const int edge_num = next_edge_num++;
            result.edge_num_by_node[node.id] = edge_num;
            out << ":" << std::setprecision(17) << node.branch_length_to_parent
                << "{" << edge_num << "}";
        }
        return out.str();
    };

    const TreeNode& root = tree.nodes[tree.root_id];
    const auto should_flatten_root_child = [&](int child_id) -> bool {
        if (child_id < 0 || child_id >= static_cast<int>(tree.nodes.size())) return false;
        const TreeNode& child = tree.nodes[child_id];
        return !child.is_tip && child.parent == tree.root_id && std::abs(child.branch_length_to_parent) <= 1e-15;
    };

    std::ostringstream out;
    if (!root.is_tip && should_flatten_root_child(root.left) != should_flatten_root_child(root.right)) {
        const int flat_child = should_flatten_root_child(root.left) ? root.left : root.right;
        const int other_child = (flat_child == root.left) ? root.right : root.left;
        out << "(" << emit_node(flat_child, true) << "," << emit_node(other_child, false) << ")";
    } else {
        out << emit_node(tree.root_id, false);
    }

    result.tree = out.str() + ";";
    return result;
}

std::vector<JplacePlacementRecord> build_jplace_records(
    const TreeBuildResult& tree,
    const JplaceTreeExport& tree_export,
    const std::vector<PlacementResult>& placement_results,
    const std::vector<NewPlacementQuery>& placement_queries)
{
    const size_t export_count = std::min(placement_results.size(), placement_queries.size());
    std::vector<JplacePlacementRecord> records;
    records.reserve(export_count);

    for (size_t i = 0; i < export_count; ++i) {
        const PlacementResult& placement = placement_results[i];
        JplacePlacementRecord record;
        record.query_name = placement_queries[i].msa_name.empty()
            ? ("query_" + std::to_string(i))
            : placement_queries[i].msa_name;

        if (!placement.top_placements.empty()) {
            for (const PlacementResult::RankedPlacement& candidate : placement.top_placements) {
                append_jplace_row(
                    tree,
                    tree_export,
                    record,
                    candidate.target_id,
                    candidate.loglikelihood,
                    candidate.like_weight_ratio,
                    candidate.proximal_length,
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
                placement.proximal_length,
                placement.pendant_length);
        }

        if (record.rows.empty()) {
            throw std::runtime_error("Could not export any placement rows for jplace.");
        }
        records.push_back(std::move(record));
    }

    return records;
}

void write_jplace(
    const std::string& out_path,
    const std::string& tree_string,
    const std::vector<JplacePlacementRecord>& placements,
    const std::string& invocation)
{
    std::filesystem::path output_path(out_path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream out(out_path);
    if (!out) {
        throw std::runtime_error("Cannot open jplace output: " + out_path);
    }

    out << "{\n";
    out << "  \"tree\": \"" << json_escape_string(tree_string) << "\",\n";
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
    out << "  \"fields\": [\"edge_num\", \"likelihood\", \"like_weight_ratio\", \"distal_length\", \"pendant_length\"]\n";
    out << "}\n";
}

} // namespace jplaceio
} // namespace mlipper
