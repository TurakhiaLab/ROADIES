#pragma once

#include <string>
#include <vector>

#include "tree/tree.hpp"

struct LocalSprBatchRunContext {
    BuildToGpuResult& res;
    PlacementOpBuffer& placement_ops;
    cudaStream_t stream = nullptr;
    const std::vector<std::string>& inserted_names;
    std::vector<std::string>& current_names;
    std::vector<std::string>& current_rows;
    std::string& current_tree_newick;
    const std::vector<unsigned>& pattern_weights_arg;
    const std::vector<double>& rate_weights;
    const std::vector<double>& rate_multipliers;
    const std::vector<double>& pi;
    size_t sites = 0;
    int states = 0;
    int rate_cats = 0;
    bool per_rate_scaling = false;
    bool profile_batch_timing = false;
    bool local_spr_fast = false;
    int local_spr_radius = 4;
    int local_spr_cluster_threshold = 0;
    int local_spr_topk_per_unit = 0;
    bool local_spr_dynamic_validation_conflicts = false;
    int local_spr_rounds = 0;
};

void run_local_spr_batch_refinement(LocalSprBatchRunContext& ctx);
