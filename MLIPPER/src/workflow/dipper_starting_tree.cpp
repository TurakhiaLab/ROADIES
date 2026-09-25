#include "dipper_starting_tree.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

// dipper_api.hpp currently uses the stream manipulators and unordered_set above
// without including their standard headers itself; keep those prerequisites here.
#include "dipper_api.hpp"
namespace mlipper::workflow {
namespace {

MashPlacement::Param makeDipperTreeParams()
{
    MashPlacement::Param params(
        /*kmerSize=*/21,
        /*sketchSize=*/1000,
        /*threshold=*/0,
        /*distanceType=*/2,
        /*in=*/"m",
        /*out=*/"t");
    params.range = {0, -1};
    params.isProtein = false;
    params.distanceType = 2;
    return params;
}

} // namespace

std::string buildDipperStartingTree(
    const parse::Alignment& alignment,
    DipperTreeMode tree_mode)
{
    if (alignment.names.empty() ||
        alignment.names.size() != alignment.sequences.size() ||
        alignment.sites == 0) {
        throw std::invalid_argument(
            "DIPPER starting tree requires a non-empty aligned data set");
    }
    for (const std::string& sequence : alignment.sequences) {
        if (sequence.size() != alignment.sites) {
            throw std::invalid_argument(
                "DIPPER sequence length does not match alignment sites");
        }
    }

    // DipperSession owns all adapter-local resources. The caller's active GPU
    // reservation remains in force, but no DIPPER pointer escapes this scope.
    DipperSession dipper_session;
    dipper_session.loadSequences(alignment.sequences, alignment.names);
    MashPlacement::Param params = makeDipperTreeParams();

    const char* tree_mode_name = nullptr;
    switch (tree_mode) {
        case DipperTreeMode::NJPlacement:
            tree_mode_name = "nj-placement";
            dipper_session.initializeGPU(params);
            dipper_session.placeAllSequences();
            break;
        case DipperTreeMode::DivideAndConquer:
            tree_mode_name = "divide-and-conquer";
            dipper_session.initializeDivideAndConquerGPU(params);
            dipper_session.buildDivideAndConquerTree();
            break;
        default:
            throw std::logic_error("unsupported DIPPER starting-tree mode");
    }

    std::cout << "DIPPER starting tree completed: mode="
              << tree_mode_name
              << " tips=" << alignment.names.size() << "\n";
    std::string newick = dipper_session.exportTreeNewickString(true);
    if (newick.empty()) {
        throw std::runtime_error("DIPPER returned an empty starting tree");
    }
    return newick;
}

} // namespace mlipper::workflow
