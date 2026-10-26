#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include "util/msa_preprocess.hpp"

// Keep reference and query columns coupled throughout compression: placement
// likelihood is preserved only when a duplicate pattern matches on both sides.
PreprocessedAlignments preprocess_alignments(
    const parse::Alignment& tree_alignment,
    const parse::Alignment& query_alignment)
{
    PreprocessedAlignments result;
    result.tree_alignment = tree_alignment;
    result.query_alignment = query_alignment;

    std::vector<mlipper::SequenceRecord> placement_queries =
        build_placement_query(
            query_alignment.names,
            query_alignment.sequences);
    size_t sites = tree_alignment.sites;
    result.pattern_weights.assign(sites, 1u);
    remove_repetitive_columns(
        result.tree_alignment.sequences,
        placement_queries,
        result.pattern_weights,
        sites);
    if (sites == 0) {
        throw std::runtime_error(
            "All columns were removed after repetitive-column compression.");
    }

    result.tree_alignment.sites = sites;
    result.query_alignment.sequences.clear();
    result.query_alignment.sequences.reserve(placement_queries.size());
    for (const mlipper::SequenceRecord& query : placement_queries) {
        result.query_alignment.sequences.push_back(query.sequence);
    }
    result.query_alignment.sites = sites;
    return result;
}

struct ColumnInfo {
    size_t id;
    std::string pattern;

    bool operator<(const ColumnInfo& other) const { return pattern < other.pattern; }
    bool operator==(const ColumnInfo& other) const { return pattern == other.pattern; }
};

namespace {

void validate_alignment_widths(
    const std::vector<std::string>& rows,
    const std::vector<mlipper::SequenceRecord>& queries,
    size_t sites)
{
    for (const std::string& row : rows) {
        if (row.size() != sites) {
            throw std::runtime_error(
                "Reference sequence length does not match alignment sites.");
        }
    }
    for (const mlipper::SequenceRecord& query : queries) {
        if (query.sequence.size() != sites) {
            throw std::runtime_error(
                "Query sequence length does not match alignment sites.");
        }
    }
}

} // namespace

void remove_repetitive_columns(
    std::vector<std::string>& rows,
    std::vector<mlipper::SequenceRecord>& queries,
    std::vector<unsigned>& pattern_weights,
    size_t& sites
) {
    if (sites == 0) return;
    if (pattern_weights.size() != sites) {
        throw std::runtime_error(
            "pattern_weights size mismatch in remove_repetitive_columns");
    }
    validate_alignment_widths(rows, queries, sites);
    if (queries.size() > std::numeric_limits<size_t>::max() - rows.size()) {
        throw std::length_error("Combined alignment sequence count overflow.");
    }
    const size_t total_sequences = rows.size() + queries.size();
    std::vector<ColumnInfo> columns(sites);
    // A site is duplicate only when both its reference and query states match.
    // Including queries here prevents compression from changing placement
    // likelihoods for sites that look identical in the reference alone.
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, sites),
        [&](const tbb::blocked_range<size_t>& range) {
            for (size_t col = range.begin(); col != range.end(); ++col) {
                std::string pattern;
                pattern.reserve(total_sequences);
                for (const auto& row : rows) {
                    pattern += static_cast<char>(
                        encode_state_DNA4_mask(row[col]));
                }
                for (const auto& query : queries) {
                    pattern += static_cast<char>(
                        encode_state_DNA4_mask(query.sequence[col]));
                }

                columns[col].pattern = std::move(pattern);
                columns[col].id = col;
            }
    });

    tbb::parallel_sort(columns.begin(), columns.end());

    std::vector<size_t> unique_id(sites);
    std::vector<unsigned> new_pattern_weights(sites, 0u);
    unique_id[0] = columns[0].id;
    new_pattern_weights[0] = pattern_weights[columns[0].id];
    size_t filtered_sites = 1;

    for (size_t i = 1; i < sites; ++i) {
        if (!(columns[i] == columns[i - 1])) {
            unique_id[filtered_sites] = columns[i].id;
            new_pattern_weights[filtered_sites] = pattern_weights[columns[i].id];
            ++filtered_sites;
        } else {
            const unsigned weight = pattern_weights[columns[i].id];
            if (weight > std::numeric_limits<unsigned>::max() -
                    new_pattern_weights[filtered_sites - 1]) {
                throw std::length_error("Compressed pattern weight overflow.");
            }
            new_pattern_weights[filtered_sites - 1] += weight;
        }
    }

    // Sorting groups duplicates but changes column order. Restore the retained
    // representatives to input order before rebuilding every sequence.
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

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, rows.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i < range.end(); ++i) {
                std::string pattern;
                pattern.reserve(filtered_sites);
                for (size_t j = 0; j < filtered_sites; ++j) {
                    const size_t col = sorted_unique_id[j];
                    pattern += rows[i][col];
                }
                new_rows[i] = std::move(pattern);
            }
    });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, queries.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i < range.end(); ++i) {
                std::string pattern;
                pattern.reserve(filtered_sites);
                for (size_t j = 0; j < filtered_sites; ++j) {
                    const size_t col = sorted_unique_id[j];
                    pattern += queries[i].sequence[col];
                }
                queries[i].sequence = std::move(pattern);
            }
    });

    rows = std::move(new_rows);
    pattern_weights = std::move(sorted_pattern_weights);
    sites = filtered_sites;
}
