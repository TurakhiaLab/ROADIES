#pragma once

#include <string>

#include "io/parse_file.hpp"

namespace mlipper::workflow {

enum class DipperTreeMode {
    NJPlacement,
    DivideAndConquer,
};

// Converts one validated DNA alignment into a Newick starting tree. DIPPER's
// temporary CPU/GPU state is confined to the call; the returned string owns no
// adapter resources.
std::string buildDipperStartingTree(
    const parse::Alignment& alignment,
    DipperTreeMode tree_mode);

} // namespace mlipper::workflow
