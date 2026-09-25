#include "io/tree_newick.hpp"
#include "tree/tree.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

bool newick_name_requires_quotes(const std::string& name) {
    if (name.empty()) return true;
    for (char ch : name) {
        switch (ch) {
            case '(':
            case ')':
            case '[':
            case ']':
            case ':':
            case ';':
            case ',':
            case '\'':
            case ' ':
            case '\t':
            case '\n':
            case '\r':
                return true;
            default:
                break;
        }
    }
    return false;
}

std::string format_newick_name(const std::string& name) {
    if (!newick_name_requires_quotes(name)) return name;
    std::string quoted;
    quoted.reserve(name.size() + 2);
    quoted.push_back('\'');
    for (char ch : name) {
        quoted.push_back(ch);
        if (ch == '\'') quoted.push_back('\'');
    }
    quoted.push_back('\'');
    return quoted;
}

// Detached serialization view. Keeping it separate lets optional short-edge
// collapse rewrite output multifurcations without mutating the binary resident
// tree or invalidating its node IDs.
struct OutputTreeNode {
    int source_node_id = -1;
    double branch_length_to_parent = 0.0;
    std::string name;
    std::vector<OutputTreeNode> children;
};

OutputTreeNode build_output_subtree(const TreeBuildResult& tree, int node_id) {
    if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) {
        throw std::runtime_error("Invalid node id while preparing Newick output tree.");
    }

    const TreeNode& node = tree.nodes[node_id];
    OutputTreeNode out;
    out.source_node_id = node.id;

    if (node.parent >= 0) {
        out.branch_length_to_parent = static_cast<double>(node.branch_length_to_parent);
        if (!std::isfinite(out.branch_length_to_parent) ||
            out.branch_length_to_parent < 0.0) {
            throw std::runtime_error(
                "Invalid branch length while preparing Newick output tree for node " +
                std::to_string(node.id));
        }
    }

    if (node.is_tip) {
        out.name = node.name.empty() ? ("tip_" + std::to_string(node.id)) : node.name;
        return out;
    }

    if (node.left < 0 || node.right < 0) {
        throw std::runtime_error("Internal node missing child while preparing Newick output tree.");
    }

    out.children.reserve(2);
    out.children.push_back(build_output_subtree(tree, node.left));
    out.children.push_back(build_output_subtree(tree, node.right));
    return out;
}

// Contract an internal edge in the serialization view and transfer its length
// to each promoted child so root-to-tip path lengths remain unchanged.
void collapse_short_internal_output_branches(
    OutputTreeNode& node,
    double epsilon)
{
    if (node.children.empty()) return;

    for (auto& child : node.children) {
        collapse_short_internal_output_branches(child, epsilon);
    }

    std::vector<OutputTreeNode> rewritten_children;
    rewritten_children.reserve(node.children.size());
    for (auto& child : node.children) {
        const bool collapse_child =
            !child.children.empty() &&
            child.branch_length_to_parent <= epsilon;
        if (!collapse_child) {
            rewritten_children.push_back(std::move(child));
            continue;
        }

        const double collapsed_length = child.branch_length_to_parent;
        for (auto& grandchild : child.children) {
            grandchild.branch_length_to_parent += collapsed_length;
            rewritten_children.push_back(std::move(grandchild));
        }
    }

    node.children = std::move(rewritten_children);
}

void write_newick_subtree(
    const TreeBuildResult& tree,
    int node_id,
    std::ostream& os) {
    if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) {
        throw std::runtime_error("Invalid node id while writing Newick tree.");
    }

    const TreeNode& node = tree.nodes[node_id];
    if (node.is_tip) {
        os << format_newick_name(node.name.empty() ? ("tip_" + std::to_string(node.id)) : node.name);
    } else {
        if (node.left < 0 || node.right < 0) {
            throw std::runtime_error("Internal node missing child while writing Newick tree.");
        }
        os << '(';
        write_newick_subtree(tree, node.left, os);
        os << ',';
        write_newick_subtree(tree, node.right, os);
        os << ')';
    }

    if (node.parent >= 0) {
        const double branch_length = static_cast<double>(node.branch_length_to_parent);
        if (!std::isfinite(branch_length) || branch_length < 0.0) {
            throw std::runtime_error(
                "Invalid branch length while writing Newick tree for node " +
                std::to_string(node.id) + " parent " +
                std::to_string(node.parent) + " value " +
                std::to_string(branch_length));
        }
        os << ':' << std::setprecision(17) << branch_length;
    }
}

void write_newick_subtree(
    const OutputTreeNode& node,
    bool is_root,
    std::ostream& os)
{
    if (node.children.empty()) {
        os << format_newick_name(node.name.empty() ? ("tip_" + std::to_string(node.source_node_id)) : node.name);
    } else {
        if (node.children.size() < 2) {
            throw std::runtime_error(
                "Internal output node collapsed below arity 2 while writing Newick tree.");
        }
        os << '(';
        for (size_t i = 0; i < node.children.size(); ++i) {
            if (i > 0) os << ',';
            write_newick_subtree(node.children[i], false, os);
        }
        os << ')';
    }

    if (!is_root) {
        if (!std::isfinite(node.branch_length_to_parent) ||
            node.branch_length_to_parent < 0.0) {
            throw std::runtime_error(
                "Invalid branch length while writing collapsed Newick tree for node " +
                std::to_string(node.source_node_id));
        }
        os << ':' << std::setprecision(17) << node.branch_length_to_parent;
    }
}

} // namespace

namespace mlipper {
namespace treeio {

std::string format_newick_taxon_name(const std::string& name) {
    return format_newick_name(name);
}

void validate_tree_for_output(const TreeBuildResult& tree)
{
    if (tree.root_id < 0 ||
        tree.root_id >= static_cast<int>(tree.nodes.size())) {
        throw std::runtime_error("Cannot serialize tree: invalid root_id.");
    }
    if (tree.nodes[static_cast<size_t>(tree.root_id)].parent != -1) {
        throw std::runtime_error("Cannot serialize tree: root has a parent.");
    }

    struct Visit {
        int node_id;
        bool exiting;
    };
    std::vector<uint8_t> state(tree.nodes.size(), 0);
    std::vector<Visit> stack{{tree.root_id, false}};
    while (!stack.empty()) {
        const Visit visit = stack.back();
        stack.pop_back();
        const size_t index = static_cast<size_t>(visit.node_id);
        if (visit.exiting) {
            state[index] = 2;
            continue;
        }
        if (state[index] != 0) {
            throw std::runtime_error(
                "Cannot serialize tree: cycle or shared child detected.");
        }
        state[index] = 1;
        const TreeNode& node = tree.nodes[index];
        if (node.id != visit.node_id) {
            throw std::runtime_error(
                "Cannot serialize tree: node ID does not match its array slot.");
        }
        if (node.parent >= 0) {
            const double length = node.branch_length_to_parent;
            if (!std::isfinite(length) || length < 0.0) {
                throw std::runtime_error(
                    "Cannot serialize tree: invalid branch length.");
            }
        }
        if (node.is_tip) {
            if (node.left >= 0 || node.right >= 0) {
                throw std::runtime_error(
                    "Cannot serialize tree: tip has children.");
            }
            state[index] = 2;
            continue;
        }
        if (node.left < 0 || node.right < 0 || node.left == node.right ||
            node.left >= static_cast<int>(tree.nodes.size()) ||
            node.right >= static_cast<int>(tree.nodes.size())) {
            throw std::runtime_error(
                "Cannot serialize tree: internal node has invalid children.");
        }
        if (tree.nodes[static_cast<size_t>(node.left)].parent != node.id ||
            tree.nodes[static_cast<size_t>(node.right)].parent != node.id) {
            throw std::runtime_error(
                "Cannot serialize tree: child/parent links disagree.");
        }
        stack.push_back({visit.node_id, true});
        stack.push_back({node.right, false});
        stack.push_back({node.left, false});
    }
    if (std::find(state.begin(), state.end(), uint8_t{0}) != state.end()) {
        throw std::runtime_error(
            "Cannot serialize tree: topology contains disconnected nodes.");
    }
}

std::string write_tree_to_newick_string(const TreeBuildResult& tree)
{
    validate_tree_for_output(tree);
    std::ostringstream oss;
    write_newick_subtree(tree, tree.root_id, oss);
    oss << ';';
    return oss.str();
}

void write_tree_to_newick_file(
    const TreeBuildResult& tree,
    const std::string& path,
    double collapse_internal_epsilon)
{
    validate_tree_for_output(tree);

    std::filesystem::path output_path(path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("Cannot open Newick output path: " + path);
    }
    if (collapse_internal_epsilon < 0.0) {
        ofs << write_tree_to_newick_string(tree) << '\n';
    } else {
        if (!std::isfinite(collapse_internal_epsilon)) {
            throw std::runtime_error(
                "Newick output collapse epsilon must be finite.");
        }
        OutputTreeNode output_root = build_output_subtree(tree, tree.root_id);
        collapse_short_internal_output_branches(
            output_root, collapse_internal_epsilon);
        write_newick_subtree(output_root, true, ofs);
        ofs << ";\n";
    }
    if (!ofs) {
        throw std::runtime_error("Failed while writing Newick output: " + path);
    }
}

} // namespace treeio

} // namespace mlipper
