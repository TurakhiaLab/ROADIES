#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "io/jplace.hpp"
#include "io/parse_file.hpp"
#include "io/tree_newick.hpp"
#include "placement/placement.cuh"
#include "tree/tree.hpp"

namespace {

template <typename Function>
void expect_runtime_error(Function&& function)
{
    try {
        function();
    } catch (const std::runtime_error&) {
        return;
    }
    assert(false && "Expected operation to fail");
}

std::filesystem::path write_temp_file(
    const std::string& filename,
    const std::string& contents)
{
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / filename;
    std::ofstream output(path);
    output << contents;
    assert(output.good());
    return path;
}

TreeBuildResult make_test_tree()
{
    TreeBuildResult tree;
    tree.nodes.resize(3);
    tree.root_id = 2;

    tree.nodes[0].id = 0;
    tree.nodes[0].parent = 2;
    tree.nodes[0].is_tip = true;
    tree.nodes[0].name = "A B";
    tree.nodes[0].branch_length_to_parent = fp_t(0.1);

    tree.nodes[1].id = 1;
    tree.nodes[1].parent = 2;
    tree.nodes[1].is_tip = true;
    tree.nodes[1].name = "O'Brien";
    tree.nodes[1].branch_length_to_parent = fp_t(0.2);

    tree.nodes[2].id = 2;
    tree.nodes[2].parent = -1;
    tree.nodes[2].left = 0;
    tree.nodes[2].right = 1;
    tree.nodes[2].is_tip = false;
    return tree;
}

TreeBuildResult make_artificial_root_tree()
{
    TreeBuildResult tree;
    tree.nodes.resize(5);
    tree.root_id = 4;

    for (int id = 0; id < 3; ++id) {
        tree.nodes[id].id = id;
        tree.nodes[id].is_tip = true;
        tree.nodes[id].name = std::string(1, static_cast<char>('A' + id));
        tree.nodes[id].branch_length_to_parent = fp_t(0.1);
    }
    tree.nodes[0].parent = 3;
    tree.nodes[1].parent = 3;
    tree.nodes[2].parent = 4;

    tree.nodes[3].id = 3;
    tree.nodes[3].parent = 4;
    tree.nodes[3].left = 0;
    tree.nodes[3].right = 1;
    tree.nodes[3].is_tip = false;
    tree.nodes[3].branch_length_to_parent = fp_t(0.0);

    tree.nodes[4].id = 4;
    tree.nodes[4].parent = -1;
    tree.nodes[4].left = 3;
    tree.nodes[4].right = 2;
    tree.nodes[4].is_tip = false;
    return tree;
}

} // namespace

int main()
{
    const std::filesystem::path fasta_path = write_temp_file(
        "mlipper_tree_io_test.fa",
        ">A\nACGTACGT\n>B\nTGCATGCA\n");
    const parse::Alignment fasta = parse::read_alignment_file(fasta_path.string());
    std::filesystem::remove(fasta_path);
    assert(fasta.names == std::vector<std::string>({"A", "B"}));
    assert(fasta.sites == 8);

    const std::filesystem::path sequential_path = write_temp_file(
        "mlipper_tree_io_test_sequential.phy",
        "2 8\nA ACGTACGT\nB TGCATGCA\n");
    const parse::Alignment sequential =
        parse::read_alignment_file(sequential_path.string());
    std::filesystem::remove(sequential_path);
    assert(sequential.names == std::vector<std::string>({"A", "B"}));
    assert(sequential.sequences[0] == "ACGTACGT");

    const std::filesystem::path interleaved_path = write_temp_file(
        "mlipper_tree_io_test_interleaved.phy",
        "2 8\nA ACGT\nB TGCA\n\nACGT\nTGCA\n");
    const parse::Alignment interleaved =
        parse::read_alignment_file(interleaved_path.string());
    std::filesystem::remove(interleaved_path);
    assert(interleaved.names == std::vector<std::string>({"A", "B"}));
    assert(interleaved.sequences[1] == "TGCATGCA");

    const TreeBuildResult tree = make_test_tree();
    const std::string ordinary =
        mlipper::treeio::write_tree_to_newick_string(tree);
    assert(ordinary.find("'A B'") != std::string::npos);
    assert(ordinary.find("'O''Brien'") != std::string::npos);

    RawPlacementResult placement;
    placement.target_id = 0;
    placement.loglikelihood = -1.0;
    placement.distal_length = 0.05;
    placement.pendant_length = 0.1;
    const std::filesystem::path jplace_path =
        std::filesystem::temp_directory_path() /
        "mlipper_tree_io_test.jplace";
    mlipper::jplaceio::write_jplace(
        jplace_path.string(),
        tree,
        std::vector<RawPlacementResult>{placement},
        std::vector<std::string>{"query \"one\""},
        "tree_io_test");
    std::ifstream jplace_file(jplace_path);
    const std::string jplace{
        std::istreambuf_iterator<char>(jplace_file),
        std::istreambuf_iterator<char>()};
    std::filesystem::remove(jplace_path);
    assert(jplace.find("'A B'") != std::string::npos);
    assert(jplace.find("'O''Brien'") != std::string::npos);
    assert(jplace.find("\"distal_length\"") != std::string::npos);
    assert(jplace.find("0.050000000000000003") != std::string::npos);
    assert(jplace.find("query \\\"one\\\"") != std::string::npos);

    const TreeBuildResult artificial_root_tree = make_artificial_root_tree();
    RawPlacementResult artificial_root_placement;
    artificial_root_placement.top_placements.resize(2);
    artificial_root_placement.top_placements[0].target_id = 3;
    artificial_root_placement.top_placements[0].loglikelihood = -1.0;
    artificial_root_placement.top_placements[0].like_weight_ratio = 0.2;
    artificial_root_placement.top_placements[0].distal_length = 0.0;
    artificial_root_placement.top_placements[0].pendant_length = 0.1;
    artificial_root_placement.top_placements[1].target_id = 0;
    artificial_root_placement.top_placements[1].loglikelihood = -2.0;
    artificial_root_placement.top_placements[1].like_weight_ratio = 0.8;
    artificial_root_placement.top_placements[1].distal_length = 0.05;
    artificial_root_placement.top_placements[1].pendant_length = 0.1;
    mlipper::jplaceio::write_jplace(
        jplace_path.string(),
        artificial_root_tree,
        std::vector<RawPlacementResult>{artificial_root_placement},
        std::vector<std::string>{"query"},
        "tree_io_test");
    std::ifstream artificial_jplace_file(jplace_path);
    const std::string artificial_jplace{
        std::istreambuf_iterator<char>(artificial_jplace_file),
        std::istreambuf_iterator<char>()};
    assert(artificial_jplace.find("\"p\": [[0, -2, ") !=
           std::string::npos);
    assert(artificial_jplace.find("], [") == std::string::npos);
    assert(artificial_jplace.find("[[0, -2, 1, ") !=
           std::string::npos);

    TreeBuildResult invalid_tree = tree;
    invalid_tree.nodes[0].branch_length_to_parent =
        std::numeric_limits<fp_t>::quiet_NaN();
    expect_runtime_error([&invalid_tree] {
        (void)mlipper::treeio::write_tree_to_newick_string(invalid_tree);
    });

    TreeBuildResult cyclic_tree = tree;
    cyclic_tree.nodes[2].left = 2;
    expect_runtime_error([&cyclic_tree] {
        (void)mlipper::treeio::write_tree_to_newick_string(cyclic_tree);
    });

    placement.loglikelihood = std::numeric_limits<double>::quiet_NaN();
    expect_runtime_error([&] {
        mlipper::jplaceio::write_jplace(
            jplace_path.string(),
            tree,
            std::vector<RawPlacementResult>{placement},
            std::vector<std::string>{"query"},
            "tree_io_test");
    });
    std::filesystem::remove(jplace_path);
    return 0;
}
