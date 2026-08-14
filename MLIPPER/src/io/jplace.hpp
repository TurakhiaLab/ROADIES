#pragma once

#include <string>
#include <vector>

#include "../placement/placement.cuh"

namespace mlipper {
namespace jplaceio {

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

JplaceTreeExport build_jplace_tree_export(const TreeBuildResult& tree);

std::vector<JplacePlacementRecord> build_jplace_records(
    const TreeBuildResult& tree,
    const JplaceTreeExport& tree_export,
    const std::vector<PlacementResult>& placement_results,
    const std::vector<NewPlacementQuery>& placement_queries);

void write_jplace(
    const std::string& out_path,
    const std::string& tree_string,
    const std::vector<JplacePlacementRecord>& placements,
    const std::string& invocation);

} // namespace jplaceio
} // namespace mlipper
