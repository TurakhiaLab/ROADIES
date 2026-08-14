#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "tree/tree.hpp"

bool repetitive_column_compression_enabled();

void remove_repetitive_columns(
    std::vector<std::string>& rows,
    std::vector<NewPlacementQuery>& queries,
    std::vector<unsigned>& pattern_weights,
    size_t& sites);
