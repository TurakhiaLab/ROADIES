#include <algorithm>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "io/parse_file.hpp"
#include "tree.hpp"

static void destroy_pll_rtree(pll_rtree_t* tree)
{
    if (tree != nullptr) {
        pll_rtree_destroy(tree, nullptr);
    }
}

static void destroy_pll_utree(pll_utree_t* tree)
{
    if (tree != nullptr) {
        pll_utree_destroy(tree, nullptr);
    }
}

using PllRootedTreePtr =
    std::unique_ptr<pll_rtree_t, decltype(&destroy_pll_rtree)>;
using PllUnrootedTreePtr =
    std::unique_ptr<pll_utree_t, decltype(&destroy_pll_utree)>;
using PllNewickStringPtr = std::unique_ptr<char, decltype(&std::free)>;

static PllRootedTreePtr parse_rooted_tree_with_pll(
    const std::string& newick_text)
{
    // Keep every libpll allocation under RAII while trying both accepted input
    // shapes. The exported string uses malloc and therefore requires free.
    PllRootedTreePtr rooted(
        pll_rtree_parse_newick_string(newick_text.c_str()),
        &destroy_pll_rtree);
    if (rooted) {
        return rooted;
    }

    PllUnrootedTreePtr unrooted(
        pll_utree_parse_newick_string(newick_text.c_str()),
        &destroy_pll_utree);
    if (!unrooted) {
        throw std::runtime_error(
            "libpll could not parse the input as rooted or unrooted Newick.");
    }
    if (!unrooted->binary) {
        throw std::runtime_error(
            "MLIPPER requires a binary Newick tree; libpll parsed a multifurcating tree.");
    }

    PllNewickStringPtr rooted_newick(
        pll_utree_export_newick_rooted(unrooted->vroot, 0.0),
        &std::free);
    if (!rooted_newick) {
        throw std::runtime_error(
            "libpll failed to root the unrooted Newick tree.");
    }

    rooted.reset(pll_rtree_parse_newick_string(rooted_newick.get()));
    if (!rooted) {
        throw std::runtime_error(
            "libpll failed to parse its rooted Newick export.");
    }
    return rooted;
}

static std::string preview_names(std::vector<std::string> names, size_t limit = 5) {
    std::sort(names.begin(), names.end());
    std::ostringstream oss;
    const size_t count = std::min(limit, names.size());
    for (size_t idx = 0; idx < count; ++idx) {
        if (idx) oss << ", ";
        oss << names[idx];
    }
    if (names.size() > limit) {
        oss << " ... (" << names.size() << " total)";
    }
    return oss.str();
}

std::vector<double> build_gtr_q_matrix(
    int states,
    const parse::ModelConfig& model,
    const std::vector<double>& pi)
{
    if (states != 4) {
        throw std::runtime_error("build_gtr_q_matrix currently supports only 4-state DNA.");
    }
    if (model.rates.size() != 6) {
        throw std::runtime_error("build_gtr_q_matrix requires exactly 6 GTR rates.");
    }
    if (pi.size() != 4) {
        throw std::runtime_error("build_gtr_q_matrix requires exactly 4 equilibrium frequencies.");
    }

    std::vector<double> q_matrix(states * states, 0.0);

    auto set_pair = [&](int row, int col, double rate) {
        q_matrix[row * states + col] = rate * pi[col];
        q_matrix[col * states + row] = rate * pi[row];
    };

    set_pair(0, 1, model.rates[0]);
    set_pair(0, 2, model.rates[1]);
    set_pair(0, 3, model.rates[2]);
    set_pair(1, 2, model.rates[3]);
    set_pair(1, 3, model.rates[4]);
    set_pair(2, 3, model.rates[5]);

    for (int row = 0; row < states; ++row) {
        double row_sum = 0.0;
        for (int col = 0; col < states; ++col) {
            if (row != col) row_sum += q_matrix[row * states + col];
        }
        q_matrix[row * states + row] = -row_sum;
    }

    // Normalize Q so one branch-length unit is one expected substitution.
    double mu = 0.0;
    for (int row = 0; row < states; ++row) {
        mu -= pi[row] * q_matrix[row * states + row];
    }

    for (double& entry : q_matrix) {
        entry /= mu;
    }

    return q_matrix;
}

TreeBuildResult build_tree_from_newick_with_pll(
    const std::vector<std::string>& msa_tip_names,
    const std::string& newick_text)
{
    TreeBuildResult out;
    PllRootedTreePtr rtree = parse_rooted_tree_with_pll(newick_text);

    const unsigned int num_tips   = rtree->tip_count;
    const unsigned int num_inners = rtree->inner_count;
    const unsigned int num_nodes = num_tips + num_inners;

    out.nodes.resize(num_nodes);

    std::unordered_map<std::string,int> msa_idx;
    msa_idx.reserve(msa_tip_names.size() * 2);
    for (size_t i = 0; i < msa_tip_names.size(); ++i) {
        const auto [_, inserted] =
            msa_idx.emplace(msa_tip_names[i], static_cast<int>(i));
        if (!inserted) {
            throw std::runtime_error("Duplicate alignment taxon name: " + msa_tip_names[i]);
        }
    }

    std::vector<pll_rnode_t*> postorder_nodes(num_nodes, nullptr);
    unsigned int count = 0;

    auto cb = [](pll_rnode_t*) -> int {
        return PLL_SUCCESS;
    };

    const int rc = pll_rtree_traverse(
        rtree->root,
        PLL_TREE_TRAVERSE_POSTORDER,
        cb,
        postorder_nodes.data(),
        &count);
    if (rc != PLL_SUCCESS) {
        throw std::runtime_error(
            "pll_rtree_traverse (POSTORDER) failed.");
    }
    if (count != num_nodes) {
        throw std::runtime_error("postorder node count mismatch.");
    }

    std::unordered_map<pll_rnode_t*, int> id_of;
    id_of.reserve(num_nodes * 2);

    // Postorder IDs make children precede their parent, which is also the
    // natural contiguous schedule for the initial upward CLV evaluation.
    for (unsigned int i = 0; i < num_nodes; ++i) {
        const int node_id = static_cast<int>(i);
        id_of[postorder_nodes[i]] = node_id;
        out.nodes[i].id = node_id;
    }

    for (unsigned int i = 0; i < num_nodes; ++i) {
        pll_rnode_t* nd = postorder_nodes[i];
        TreeNode& dst = out.nodes[i];

        const bool is_tip = (nd->left == nullptr && nd->right == nullptr);
        dst.is_tip = is_tip;

        if (is_tip && nd->label) dst.name = std::string(nd->label);
        else dst.name.clear();

        if (nd->parent) {
            dst.parent = id_of[nd->parent];
            dst.branch_length_to_parent = nd->length;
        } else {
            dst.parent = -1;
            dst.branch_length_to_parent = fp_t(0);
            out.root_id = i;
        }

        if (!is_tip) {
            if (!nd->left || !nd->right) {
                throw std::runtime_error(
                    "Internal node missing children (non-binary?).");
            }
            dst.left = id_of[nd->left];
            dst.right = id_of[nd->right];
        }
    }

    // Both directions are checked: matching only tree tips would silently
    // discard alignment taxa that never appear in the topology.
    out.tip_node_by_name.reserve(num_tips * 2);
    std::vector<std::string> missing_tree_tips;
    for (const TreeNode& tn : out.nodes) {
        if (!tn.is_tip) continue;
        if (tn.name.empty()) {
            throw std::runtime_error("Encountered a tip with empty name.");
        }
        const auto [_, inserted] = out.tip_node_by_name.emplace(tn.name, tn.id);
        if (!inserted) {
            throw std::runtime_error("Duplicate tip name in Newick tree: " + tn.name);
        }
        if (msa_idx.find(tn.name) == msa_idx.end()) {
            missing_tree_tips.push_back(tn.name);
        }
    }
    if (!missing_tree_tips.empty()) {
        throw std::runtime_error(
            "Newick tips not found in alignment: " + preview_names(missing_tree_tips));
    }

    std::vector<std::string> missing_alignment_tips;
    for (const auto& entry : msa_idx) {
        if (out.tip_node_by_name.find(entry.first) == out.tip_node_by_name.end()) {
            missing_alignment_tips.push_back(entry.first);
        }
    }
    if (!missing_alignment_tips.empty()) {
        throw std::runtime_error(
            "Alignment taxa not found in Newick tree: " + preview_names(missing_alignment_tips));
    }

    out.postorder.resize(num_nodes);
    for (unsigned int i = 0; i < num_nodes; ++i) out.postorder[i] = out.nodes[i].id;
    out.preorder.resize(num_nodes);
    {
        std::vector<pll_rnode_t*> stack;
        stack.reserve(num_nodes);
        stack.push_back(rtree->root);
        unsigned int idx = 0;
        while (!stack.empty()) {
            pll_rnode_t* nd = stack.back();
            stack.pop_back();
            out.preorder[idx++] = id_of[nd];
            if (nd->right) stack.push_back(nd->right);
            if (nd->left)  stack.push_back(nd->left);
        }
        if (idx != num_nodes) {
            throw std::runtime_error("preorder node count mismatch.");
        }
    }

    return out;
}

std::vector<mlipper::SequenceRecord> build_placement_query(
    const std::vector<std::string>& msa_tip_names,
    const std::vector<std::string>& msa_rows)
{
    if (msa_tip_names.size() != msa_rows.size()) {
        throw std::runtime_error("Alignment names/sequences size mismatch.");
    }

    std::vector<mlipper::SequenceRecord> out_queries;
    out_queries.reserve(msa_tip_names.size());

    for (std::size_t i = 0; i < msa_tip_names.size(); ++i) {
        mlipper::SequenceRecord q;
        q.name = msa_tip_names[i];
        q.sequence = msa_rows[i];
        out_queries.emplace_back(std::move(q));
    }
    return out_queries;
}
