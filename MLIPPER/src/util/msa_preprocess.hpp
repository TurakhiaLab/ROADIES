#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "io/parse_file.hpp"
#include "tree/tree.hpp"

struct PreprocessedAlignments {
    // Both alignments are rewritten to the same retained-pattern order.
    parse::Alignment tree_alignment;
    parse::Alignment query_alignment;
    // Multiplicity of each retained combined reference/query site pattern.
    std::vector<unsigned> pattern_weights;
};

// Collapse columns only when their states match across the combined reference
// and query rows. The returned weights preserve the original site likelihood.
PreprocessedAlignments preprocess_alignments(
    const parse::Alignment& tree_alignment,
    const parse::Alignment& query_alignment = {});

// In-place lower-level form used after query rows have become SequenceRecords.
// sites and pattern_weights are rewritten to the compressed column count.
void remove_repetitive_columns(
    std::vector<std::string>& rows,
    std::vector<mlipper::SequenceRecord>& queries,
    std::vector<unsigned>& pattern_weights,
    size_t& sites);
