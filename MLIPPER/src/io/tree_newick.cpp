#include "io/tree_newick.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
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

struct OutputTreeNode {
    int source_node_id = -1;
    bool is_tip = false;
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
    out.is_tip = node.is_tip;

    if (node.parent >= 0) {
        out.branch_length_to_parent = static_cast<double>(node.branch_length_to_parent);
        if (out.branch_length_to_parent < 0.0) {
            throw std::runtime_error(
                "Negative branch length while preparing Newick output tree for node " +
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

size_t collapse_short_internal_output_branches(OutputTreeNode& node, double epsilon) {
    if (node.is_tip) return 0;

    size_t collapsed = 0;
    for (auto& child : node.children) {
        collapsed += collapse_short_internal_output_branches(child, epsilon);
    }

    std::vector<OutputTreeNode> rewritten_children;
    rewritten_children.reserve(node.children.size());
    for (auto& child : node.children) {
        const bool collapse_child =
            !child.is_tip &&
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
        ++collapsed;
    }

    node.children = std::move(rewritten_children);
    return collapsed;
}

void write_newick_subtree(const TreeBuildResult& tree, int node_id, std::ostream& os);
void write_newick_subtree(const OutputTreeNode& node, bool is_root, std::ostream& os);

void write_newick_subtree(const TreeBuildResult& tree, int node_id, std::ostream& os) {
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
        if (branch_length < 0.0) {
            throw std::runtime_error(
                "Negative branch length while writing Newick tree for node " +
                std::to_string(node.id));
        }
        os << ':' << std::setprecision(17) << branch_length;
    }
}

void write_newick_subtree(const OutputTreeNode& node, bool is_root, std::ostream& os) {
    if (node.is_tip) {
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
        if (node.branch_length_to_parent < 0.0) {
            throw std::runtime_error(
                "Negative branch length while writing collapsed Newick tree for node " +
                std::to_string(node.source_node_id));
        }
        os << ':' << std::setprecision(17) << node.branch_length_to_parent;
    }
}

void print_tree_rec(const TreeBuildResult& tree, int node_id, int depth) {
    if (node_id < 0) return;
    const TreeNode& node = tree.nodes[(size_t)node_id];

    for (int i = 0; i < depth; ++i) std::cout << "  ";

    std::cout << "[" << node.id << "]";
    if (node.is_tip) {
        std::cout << " (tip: " << node.name << ")";
    } else {
        std::cout << " (inner)";
    }

    if (node.parent >= 0) {
        std::cout << "  len=" << node.branch_length_to_parent
                  << "  parent=" << node.parent;
    } else {
        std::cout << "  <ROOT>";
    }
    std::cout << "\n";

    if (!node.is_tip) {
        if (node.left >= 0) print_tree_rec(tree, node.left, depth + 1);
        if (node.right >= 0) print_tree_rec(tree, node.right, depth + 1);
    }
}

} // namespace

namespace mlipper {
namespace treeio {

std::string write_tree_to_newick_string(const TreeBuildResult& tree) {
    if (tree.root_id < 0) {
        throw std::runtime_error("Cannot serialize Newick tree: invalid root_id.");
    }
    std::ostringstream oss;
    write_newick_subtree(tree, tree.root_id, oss);
    oss << ';';
    return oss.str();
}

std::string write_tree_to_output_newick_string(
    const TreeBuildResult& tree,
    double collapse_internal_epsilon,
    size_t* collapsed_internal_branches_out) {
    if (collapsed_internal_branches_out) *collapsed_internal_branches_out = 0;
    if (collapse_internal_epsilon < 0.0) {
        return write_tree_to_newick_string(tree);
    }
    if (!std::isfinite(collapse_internal_epsilon)) {
        throw std::runtime_error("Newick output collapse epsilon must be finite.");
    }
    if (tree.root_id < 0) {
        throw std::runtime_error("Cannot serialize Newick tree: invalid root_id.");
    }

    OutputTreeNode output_root = build_output_subtree(tree, tree.root_id);
    const size_t collapsed =
        collapse_short_internal_output_branches(output_root, collapse_internal_epsilon);
    if (collapsed_internal_branches_out) *collapsed_internal_branches_out = collapsed;

    std::ostringstream oss;
    write_newick_subtree(output_root, true, oss);
    oss << ';';
    return oss.str();
}

size_t write_tree_to_newick_file(
    const TreeBuildResult& tree,
    const std::string& path,
    double collapse_internal_epsilon) {
    if (tree.root_id < 0) {
        throw std::runtime_error("Cannot write Newick tree: invalid root_id.");
    }

    std::filesystem::path output_path(path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("Cannot open Newick output path: " + path);
    }
    size_t collapsed_internal_branches = 0;
    const std::string newick = write_tree_to_output_newick_string(
        tree,
        collapse_internal_epsilon,
        &collapsed_internal_branches);
    ofs << newick << '\n';
    return collapsed_internal_branches;
}

void print_tree_structure(const TreeBuildResult& tree) {
    std::cout << "==== Tree structure (indented) ====\n";
    if (tree.root_id < 0) {
        std::cout << "No root_id set!\n";
        return;
    }
    print_tree_rec(tree, tree.root_id, 0);
    std::cout << "===================================\n";
}

} // namespace treeio

} // namespace mlipper
