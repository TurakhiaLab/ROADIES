#include "io/parse_file.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

namespace parse {

namespace detail {


inline std::string trim_copy(const std::string& s) {
    size_t head = 0;
    while (head < s.size() && std::isspace(static_cast<unsigned char>(s[head]))) ++head;
    size_t tail = s.size();
    while (tail > head && std::isspace(static_cast<unsigned char>(s[tail - 1]))) --tail;
    return s.substr(head, tail - head);
}

inline std::string strip_quotes(const std::string& s) {
    if (s.size() >= 2) {
        if ((s.front() == '"' && s.back() == '"') ||
            (s.front() == '\'' && s.back() == '\'')) {
            return s.substr(1, s.size() - 2);
        }
    }
    return s;
}


Alignment read_fasta_alignment(const std::string& path) {
    Alignment aln;
    std::ifstream ifs(path);
    if (!ifs)
        throw std::runtime_error("Cannot open alignment file: " + path);
    std::string line;
    std::string current_name;
    std::string current_seq;
    auto push_current = [&]() {
        if (!current_name.empty()) {
            aln.names.push_back(current_name);
            aln.sequences.push_back(current_seq);
            current_seq.clear();
        }
    };
    while (std::getline(ifs, line)) {
        std::string trimmed = trim_copy(line);
        if (trimmed.empty()) continue;
        if (trimmed[0] == '>') {
            push_current();
            current_name = trim_copy(trimmed.substr(1));
            if (current_name.empty()) {
                throw std::runtime_error("FASTA sequence name is empty.");
            }
        } else if (!current_name.empty()) {
            for (char c : trimmed) {
                if (!std::isspace(static_cast<unsigned char>(c)))
                    current_seq.push_back(c);
            }
        }
    }
    push_current();
    if (aln.names.empty())
        throw std::runtime_error("FASTA alignment has no sequences.");
    aln.sites = aln.sequences[0].size();
    for (const auto& seq : aln.sequences) {
        if (seq.size() != aln.sites)
            throw std::runtime_error("Inconsistent sequence lengths in alignment.");
    }
    return aln;
}

std::vector<std::string> split_on_whitespace(const std::string& text) {
    std::vector<std::string> tokens;
    std::istringstream iss(text);
    std::string token;
    while (iss >> token) tokens.push_back(token);
    return tokens;
}

Alignment read_phylip_alignment(const std::string& path) {
    Alignment aln;
    std::ifstream ifs(path);
    if (!ifs)
        throw std::runtime_error("Cannot open alignment file: " + path);
    std::string header;
    while (std::getline(ifs, header)) {
        if (trim_copy(header).empty()) continue;
        break;
    }
    if (header.empty())
        throw std::runtime_error("PHYLIP file is empty.");
    auto tokens = split_on_whitespace(trim_copy(header));
    if (tokens.size() < 2)
        throw std::runtime_error("PHYLIP header must contain sequence count and site count.");
    int sequence_count = std::stoi(tokens[0]);
    int site_count = std::stoi(tokens[1]);
    std::string line;
    while (aln.names.size() < static_cast<size_t>(sequence_count) && std::getline(ifs, line)) {
        std::string trimmed = trim_copy(line);
        if (trimmed.empty()) continue;
        auto parts = split_on_whitespace(trimmed);
        if (parts.empty()) continue;
        std::string name = parts[0];
        std::string seq;
        for (size_t i = 1; i < parts.size(); ++i)
            seq += parts[i];
        if (seq.empty())
            throw std::runtime_error("Empty sequence encountered in PHYLIP.");
        if (seq.size() != static_cast<size_t>(site_count))
            throw std::runtime_error("PHYLIP sequence length mismatch header.");
        aln.names.push_back(name);
        aln.sequences.push_back(seq);
    }
    if (aln.names.size() != static_cast<size_t>(sequence_count))
        throw std::runtime_error("PHYLIP sequence count mismatch.");
    aln.sites = site_count;
    return aln;
}

Alignment read_alignment(const std::string& path) {
    std::ifstream probe(path);
    if (!probe)
        throw std::runtime_error("Cannot open alignment file: " + path);
    std::string first_line;
    if (!std::getline(probe, first_line))
        throw std::runtime_error("Alignment file is empty.");
    std::string trimmed = trim_copy(first_line);
    probe.close();
    if (trimmed.empty())
        throw std::runtime_error("Alignment file is empty.");
    if (trimmed.front() == '>')
        return read_fasta_alignment(path);
    return read_phylip_alignment(path);
}

namespace newick {

struct SimpleNode {
    std::string label;
    std::string length_str;
    bool has_length = false;
    std::vector<std::unique_ptr<SimpleNode>> children;
};

static void skip_ws(const std::string& text, size_t& pos) {
    while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
}

static std::unique_ptr<SimpleNode> parse_newick_node(const std::string& text, size_t& pos);

static std::unique_ptr<SimpleNode> parse_newick_node(const std::string& text, size_t& pos) {
    skip_ws(text, pos);
    auto node = std::make_unique<SimpleNode>();
    if (pos < text.size() && text[pos] == '(') {
        ++pos;
        while (true) {
            node->children.push_back(parse_newick_node(text, pos));
            skip_ws(text, pos);
            if (pos >= text.size())
                throw std::runtime_error("Unexpected end of Newick while parsing children.");
            if (text[pos] == ',') {
                ++pos;
                continue;
            }
            if (text[pos] == ')') {
                ++pos;
                break;
            }
            throw std::runtime_error("Unexpected character in Newick children list.");
        }
    }

    skip_ws(text, pos);
    size_t label_start = pos;
    while (pos < text.size()) {
        char c = text[pos];
        if (c == ':' || c == ',' || c == ')' || c == ';' || std::isspace(static_cast<unsigned char>(c)))
            break;
        ++pos;
    }
    if (label_start < pos)
        node->label = text.substr(label_start, pos - label_start);

    skip_ws(text, pos);
    if (pos < text.size() && text[pos] == ':') {
        ++pos;
        skip_ws(text, pos);
        size_t len_start = pos;
        while (pos < text.size()) {
            char c = text[pos];
            if (c == ',' || c == ')' || c == ';' || std::isspace(static_cast<unsigned char>(c)))
                break;
            ++pos;
        }
        if (len_start < pos) {
            node->length_str = text.substr(len_start, pos - len_start);
            node->has_length = true;
        }
    }
    return node;
}

static void resolve_polytomies(SimpleNode* node) {
    for (auto& child : node->children) resolve_polytomies(child.get());
    if (node->children.size() <= 2) return;
    std::vector<std::unique_ptr<SimpleNode>> kids = std::move(node->children);
    std::unique_ptr<SimpleNode> combined = std::move(kids[0]);
    for (size_t i = 1; i + 1 < kids.size(); ++i) {
        auto merged = std::make_unique<SimpleNode>();
        merged->children.reserve(2);
        merged->children.push_back(std::move(combined));
        merged->children.push_back(std::move(kids[i]));
        merged->has_length = true;
        merged->length_str = "0";
        combined = std::move(merged);
    }
    node->children.clear();
    node->children.push_back(std::move(combined));
    node->children.push_back(std::move(kids.back()));
}

static void append_newick(const SimpleNode& node, std::string& out) {
    if (!node.children.empty()) {
        out.push_back('(');
        for (size_t i = 0; i < node.children.size(); ++i) {
            append_newick(*node.children[i], out);
            if (i + 1 < node.children.size()) out.push_back(',');
        }
        out.push_back(')');
    }
    if (!node.label.empty()) out += node.label;
    if (node.has_length) {
        out.push_back(':');
        out += node.length_str;
    }
}

static double parse_length_or_zero(const SimpleNode& node) {
    if (!node.has_length || node.length_str.empty()) return 0.0;
    try {
        return std::stod(node.length_str);
    } catch (...) {
        return 0.0;
    }
}

static std::string format_length(double value) {
    std::ostringstream oss;
    oss << std::setprecision(17) << value;
    return oss.str();
}

static std::unique_ptr<SimpleNode> prune_tips_recursive(
    std::unique_ptr<SimpleNode> node,
    const std::unordered_set<std::string>& tip_names_to_remove)
{
    if (!node) return nullptr;

    if (node->children.empty()) {
        const std::string label = strip_quotes(trim_copy(node->label));
        if (!label.empty() && tip_names_to_remove.count(label)) {
            return nullptr;
        }
        return node;
    }

    std::vector<std::unique_ptr<SimpleNode>> kept_children;
    kept_children.reserve(node->children.size());
    for (auto& child : node->children) {
        if (auto kept = prune_tips_recursive(std::move(child), tip_names_to_remove)) {
            kept_children.push_back(std::move(kept));
        }
    }
    node->children = std::move(kept_children);

    if (node->children.empty()) {
        return nullptr;
    }

    if (node->children.size() == 1) {
        std::unique_ptr<SimpleNode> child = std::move(node->children.front());
        const double merged_length = parse_length_or_zero(*node) + parse_length_or_zero(*child);
        child->has_length = true;
        child->length_str = format_length(merged_length);
        return child;
    }

    return node;
}

std::string normalize_newick(const std::string& raw) {
    try {
        size_t pos = 0;
        auto root = parse_newick_node(raw, pos);
        skip_ws(raw, pos);
        if (pos < raw.size() && raw[pos] == ';') ++pos;
        skip_ws(raw, pos);
        if (pos != raw.size())
            throw std::runtime_error("Extra characters after Newick tree.");
        resolve_polytomies(root.get());
        std::string normalized;
        append_newick(*root, normalized);
        normalized.push_back(';');
        return normalized;
    } catch (...) {
        return raw;
    }
}

std::string prune_newick_tips(
    const std::string& raw,
    const std::unordered_set<std::string>& tip_names_to_remove)
{
    size_t pos = 0;
    auto root = parse_newick_node(raw, pos);
    skip_ws(raw, pos);
    if (pos < raw.size() && raw[pos] == ';') ++pos;
    skip_ws(raw, pos);
    if (pos != raw.size()) {
        throw std::runtime_error("Extra characters after Newick tree.");
    }

    root = prune_tips_recursive(std::move(root), tip_names_to_remove);
    if (!root) {
        throw std::runtime_error("All tips were pruned from the Newick tree.");
    }
    resolve_polytomies(root.get());

    std::string pruned;
    append_newick(*root, pruned);
    pruned.push_back(';');
    return pruned;
}

} // namespace newick

} // namespace detail

// }

Alignment read_alignment_file(const std::string& path) {
    return detail::read_alignment(path);
}

std::string normalize_newick(const std::string& raw) {
    return detail::newick::normalize_newick(raw);
}

std::string prune_newick_tips(
    const std::string& raw,
    const std::unordered_set<std::string>& tip_names_to_remove)
{
    return detail::newick::prune_newick_tips(raw, tip_names_to_remove);
}

} // namespace parse
