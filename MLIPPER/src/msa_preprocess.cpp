#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <cstdint>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include "msa_preprocess.hpp"

struct ColumnInfo {
    size_t id;
    std::string pattern;

    bool operator<(const ColumnInfo& other) const { return pattern < other.pattern; }
    bool operator==(const ColumnInfo& other) const { return pattern == other.pattern; }
};

bool repetitive_column_compression_enabled() {
    const char* disable = std::getenv("MLIPPER_DISABLE_REPETITIVE_COLUMNS");
    if (disable && disable[0] && std::string(disable) != "0") {
        return false;
    }

    const char* legacy_enable = std::getenv("MLIPPER_ENABLE_REPETITIVE_COLUMNS");
    if (legacy_enable) {
        return legacy_enable[0] && std::string(legacy_enable) != "0";
    }
    return true;
}

void remove_repetitive_columns(
    std::vector<std::string>& rows,
    std::vector<NewPlacementQuery>& queries,
    std::vector<unsigned>& pattern_weights,
    size_t& sites
) {
    if (sites == 0) return;
    if (pattern_weights.size() != sites) {
        throw std::runtime_error("pattern_weights size mismatch in remove_repetitive_columns");
    }
    size_t total_sequences = rows.size() + queries.size();
    std::vector<ColumnInfo> columns(sites);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, sites), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t col = r.begin(); col != r.end(); col++) {
            std::string pattern;
            pattern.reserve(total_sequences);
            for (const auto& row : rows) pattern += row[col];
            for (const auto& q : queries) pattern += q.msa[col];

            columns[col].pattern = std::move(pattern);
            columns[col].id      = col;
        }
    });

    tbb::parallel_sort(columns.begin(), columns.end());

    std::vector<size_t> unique_id(sites);
    std::vector<unsigned> new_pattern_weights(sites, 0u);
    unique_id[0] = columns[0].id;
    new_pattern_weights[0] = pattern_weights[columns[0].id];
    size_t filtered_sites = 1;
    
    for(size_t i = 1; i < sites; i++) {
        if(!(columns[i] == columns[i-1])) {
            unique_id[filtered_sites] = columns[i].id;
            new_pattern_weights[filtered_sites] = pattern_weights[columns[i].id];
            filtered_sites++;
        } else {
            new_pattern_weights[filtered_sites - 1] += pattern_weights[columns[i].id];
        }
    }

    std::vector<size_t> order(filtered_sites);
    for (size_t i = 0; i < filtered_sites; ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return unique_id[a] < unique_id[b];
    });

    std::vector<size_t> sorted_unique_id(filtered_sites);
    std::vector<unsigned> sorted_pattern_weights(filtered_sites, 0u);
    for (size_t i = 0; i < filtered_sites; ++i) {
        sorted_unique_id[i] = unique_id[order[i]];
        sorted_pattern_weights[i] = new_pattern_weights[order[i]];
    }

    std::vector<std::string> new_rows(rows.size());

    tbb::parallel_for(tbb::blocked_range<size_t>(0, rows.size()), [&](const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i < r.end(); i++) {
            std::string pattern;
            pattern.reserve(filtered_sites);
            for(size_t j = 0; j < filtered_sites; j++) {
                size_t col = sorted_unique_id[j];
                pattern += rows[i][col];
            }
            new_rows[i] = std::move(pattern);
        }
    });

    tbb::parallel_for(tbb::blocked_range<size_t>(0, queries.size()), [&](const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i < r.end(); i++) {
            std::string pattern;
            pattern.reserve(filtered_sites);
            for(size_t j = 0; j < filtered_sites; j++) {
                size_t col = sorted_unique_id[j];
                pattern += queries[i].msa[col];
            }
            queries[i].msa = std::move(pattern);
        }
    });

    rows = std::move(new_rows);
    pattern_weights = std::move(sorted_pattern_weights);
    sites = filtered_sites;
}
