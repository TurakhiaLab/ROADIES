#pragma once

#include <string>
#include <vector>

struct RawPlacementResult;
struct TreeBuildResult;

namespace mlipper {
namespace jplaceio {

// Serialize placements against the unchanged input tree using jplace edge
// numbers derived from tree. placement_results and query_names must describe
// the same query order. Throws on inconsistent topology or output I/O failure.
void write_jplace(
    const std::string& output_path,
    const TreeBuildResult& tree,
    const std::vector<::RawPlacementResult>& placement_results,
    const std::vector<std::string>& query_names,
    const std::string& invocation);

} // namespace jplaceio
} // namespace mlipper
