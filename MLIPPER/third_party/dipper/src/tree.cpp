/**
 * tree.cpp
 *
 * Node and Tree implementation: constructors, DFS, Newick parse/serialization,
 * leaf/internal counts, string split/strip helpers.
 */

#include "tree.hpp"
#include <cassert>
#include <cstdint>
#include <functional>
#include <map>
#include <random>
#include <set>
#include <unordered_set>
Node* lcs(Node* node1, Node* node2);


Node::Node(const std::string& id, double len){
    this->name = id;
    level = 1;
    bl = len;
    parent = nullptr;
}


Node::Node(const std::string& id, Node* par, double len){
    this->name = id;
    bl = len;
    parent = par;
    level = par->level + 1;
    par->children.push_back(this);
}

/** Recursive count of leaves under this node. */
size_t Node::getNumLeaves(){
    size_t num_leaves = 0;
    if (children.size() == 0) return num_leaves;
    for (auto ch: children){
        if (ch->is_leaf()) num_leaves += 1;
        else num_leaves += ch->getNumLeaves();
    }
    return num_leaves;
}

/** Recursive count of nodes in subtree. */
size_t Node::getNumNodes(){
    size_t num_nodes = 1;
    if (children.size() == 0) return num_nodes;
    for (auto ch: children){
        num_nodes += ch->getNumNodes();
    }
    return num_nodes;
}

/** Split on delim; respect quoted segments (odd count of '). */
void stringSplit (std::string const& s, char delim, std::vector<std::string>& words) {
    size_t start_pos = 0, end_pos = 0, temp_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        if (end_pos >= s.length()) {
            break;
        }
        std::string sub;
        if (temp_pos == 0) {
            sub = s.substr(start_pos, end_pos-start_pos);
            if (std::count(sub.begin(), sub.end(), '\'') % 2 == 1) {
                temp_pos = start_pos;
            }
            else {
                words.emplace_back(sub);
            }
        }
        else {
            sub = s.substr(temp_pos, end_pos-temp_pos);
            if (std::count(sub.begin(), sub.end(), '\'') % 2 == 0) {
                temp_pos = 0;
                words.emplace_back(sub);
            }
        }
        // words.emplace_back(s.substr(start_pos, end_pos-start_pos));
        start_pos = end_pos+1;
    }
    auto last = s.substr(start_pos, s.size()-start_pos);
    if (last != "") {
        words.push_back(std::move(last));
    }
}

std::string stripString(std::string s){
    while(s.length() && s[s.length() - 1] == ' '){
        s.pop_back();
    }
    for(size_t i = 0; i < s.length(); i++){
        if(s[i] != ' '){
            return s.substr(i);
        }
    }
    return s;
}

/** Append node and all descendants to vec in DFS order. */
void Tree::dfsExpansion(Node* node,
                                     std::vector< Node* >& vec) {
    vec.push_back(node);
    for(auto child: node->children) {
        dfsExpansion(child, vec);
    }
}

int dfsExpansionSize(Node* node) {
    int c = 0;
    if (node->children.size() == 0) {
        c++;
        return c;
    }
    for (auto &n: node->children){
        c += dfsExpansionSize(n);
    }
    std::cout << node->name << "\t" << c << std::endl;
    return c;
}

/** Serialize subtree rooted at node to Newick (with branch lengths). */
std::string Tree::getNewickString(Node* node) {

    // traversal to print each node subtree size
    // int s = dfsExpansionSize(node);
    // exit(0);
    std::vector< Node* > traversal;
    dfsExpansion(node, traversal);
    std::string newick;

    if (traversal.size() == 1) {
        newick += node->name;
        return newick;
    }

    size_t level_offset = node->level-1;
    size_t curr_level = 0;
    bool prev_open = true;

    std::stack<std::string> node_stack;
    std::stack<float> branch_length_stack;

    for (auto n: traversal) {
        size_t level = n->level-level_offset;
        float branch_length = n->bl;

        if(curr_level < level) {
            if (!prev_open) {
                newick += ',';
            }
            size_t l = level - 1;
            if (curr_level > 1) {
                l = level - curr_level;
            }
            for (size_t i=0; i < l; i++) {
                newick += '(';
                prev_open = true;
            }
            if (n->children.size() == 0) {

                newick += n->name;

                if (branch_length >= 0) {
                    newick += ':';
                    newick += std::to_string(branch_length);
                }
                prev_open = false;
            } else {
                node_stack.push(n->name);
                branch_length_stack.push(branch_length);
            }
        } else if (curr_level > level) {
            prev_open = false;
            for (size_t i = level; i < curr_level; i++) {
                newick += ')';

                newick += node_stack.top();

                if (branch_length_stack.top() >= 0) {
                    newick += ':';
                    newick += std::to_string(branch_length_stack.top());
                }
                node_stack.pop();
                branch_length_stack.pop();
            }
            if (n->children.size() == 0) {
                newick += ',';
                newick += n->name;

                if (branch_length >= 0) {
                    newick += ':';
                    newick += std::to_string(branch_length);
                }
            } else {
                node_stack.push(n->name);
                branch_length_stack.push(branch_length);
            }
        } else {
            prev_open = false;
            if (n->children.size() == 0) {

                newick += ',';
                newick += n->name;

                if (branch_length >= 0) {
                    newick += ':';
                    newick += std::to_string(branch_length);
                }
            } else {
                node_stack.push(n->name);
                branch_length_stack.push(branch_length);
            }
        }
        curr_level = level;
    }
    size_t remaining = node_stack.size();
    for (size_t i = 0; i < remaining; i++) {
        newick += ')';
        newick += node_stack.top();

        if (branch_length_stack.top() >= 0) {
            newick += ':';
            newick += std::to_string(branch_length_stack.top());
        }
        node_stack.pop();
        branch_length_stack.pop();
    }

    newick += ';';
    return newick;
}

void get_node(Node* node, std::string node_id, std::vector<Node*>& foundNode){
    if (node==nullptr) return;
    for (auto &child: node->children)
        get_node(child, node_id, foundNode);

    if (node->name==node_id) foundNode.push_back(node);
}

void rename_tree_ids(Node* node, int& idx, std::string old_name, int& new_idx){
    if (node==nullptr) return;
    for (auto &child: node->children)
        rename_tree_ids(child, idx, old_name, new_idx);
    
    if (node->children.size()!=0){
        if (node->name==old_name)
            new_idx=idx;
        node->name=std::to_string(idx++);
    }
}

void checkNode(Node* node, std::string name){
    if (node==nullptr) return;
    for (auto child: node->children){
        checkNode(child, name);
    }
    if (node->name==name) std::cout << "node still found" << std::endl;
}
std::vector<std::string> Tree::remove_and_rename(std::string node_id){
    std::vector<Node*> foundNode;
    get_node(this->root, node_id, foundNode);

    for (auto& child: foundNode[0]->parent->children) {
        if (child->name != node_id){
            foundNode.push_back(child);
            break;
        }
    }
    std::string old_name = foundNode[1]->name;

    if (foundNode[1]->parent->parent){
        foundNode[1]->parent->parent->children.erase(std::remove(foundNode[1]->parent->parent->children.begin(), foundNode[1]->parent->parent->children.end(), foundNode[1]->parent), foundNode[1]->parent->parent->children.end());
        foundNode[1]->parent->parent->children.push_back(foundNode[1]);
        foundNode[1]->parent = foundNode[1]->parent->parent;
    }

    int new_id;
    int start_id=0;
    rename_tree_ids(this->root, start_id, old_name, new_id);
    
    std::vector<std::string> nodestoroot;
    Node* currNode = foundNode[1];
    while(currNode!=nullptr){
        nodestoroot.push_back(currNode->name);
        currNode=currNode->parent;
    }

    //check
    checkNode(this->root, node_id);

    return nodestoroot;

}

Tree* Tree::fromUnrootedNewick(const std::string& newick, size_t totalLeaves) {
    return new Tree(newick, totalLeaves);
}

namespace {

std::string trimNewickToken(std::string s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) {
        s.erase(s.begin());
    }
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) {
        s.pop_back();
    }
    return s;
}

size_t findMatchingCloseParen(const std::string& s, size_t openPos) {
    if (openPos >= s.size() || s[openPos] != '(') {
        return std::string::npos;
    }
    int depth = 0;
    bool inQuote = false;
    for (size_t i = openPos; i < s.size(); ++i) {
        char c = s[i];
        if (c == '\'') {
            inQuote = !inQuote;
            continue;
        }
        if (inQuote) {
            continue;
        }
        if (c == '(') {
            ++depth;
        } else if (c == ')') {
            --depth;
            if (depth == 0) {
                return i;
            }
        }
    }
    return std::string::npos;
}

void splitTopLevelCommasNewick(const std::string& inner, std::vector<std::string>& parts) {
    parts.clear();
    if (inner.empty()) {
        return;
    }
    int depth = 0;
    bool inQuote = false;
    size_t start = 0;
    for (size_t i = 0; i < inner.size(); ++i) {
        char c = inner[i];
        if (c == '\'') {
            inQuote = !inQuote;
            continue;
        }
        if (inQuote) {
            continue;
        }
        if (c == '(') {
            ++depth;
        } else if (c == ')') {
            --depth;
        } else if (c == ',' && depth == 0) {
            parts.push_back(trimNewickToken(inner.substr(start, i - start)));
            start = i + 1;
        }
    }
    parts.push_back(trimNewickToken(inner.substr(start)));
}

std::pair<std::string, double> parseNameAndBranchLength(const std::string& s) {
    std::string t = trimNewickToken(s);
    bool inQuote = false;
    size_t colonPos = std::string::npos;
    for (size_t i = 0; i < t.size(); ++i) {
        if (t[i] == '\'') {
            inQuote = !inQuote;
        } else if (t[i] == ':' && !inQuote) {
            colonPos = i;
            break;
        }
    }
    std::string name;
    double bl = 1.0;
    if (colonPos == std::string::npos) {
        name = t;
    } else {
        name = trimNewickToken(t.substr(0, colonPos));
        std::string br = trimNewickToken(t.substr(colonPos + 1));
        if (!br.empty()) {
            bl = std::stod(br);
        }
    }
    if (name.size() >= 2 && name.front() == '\'' && name.back() == '\'') {
        name = name.substr(1, name.size() - 2);
    }
    return {name, bl};
}

void treeDepthStats(Node* node, size_t rootLevel, size_t& maxDepth, size_t& sumDepth, size_t& numLeaves) {
    if (node->children.empty()) {
        size_t d = node->level > rootLevel ? node->level - rootLevel : 0;
        maxDepth = std::max(maxDepth, d);
        sumDepth += d;
        ++numLeaves;
        return;
    }
    for (Node* ch : node->children) {
        treeDepthStats(ch, rootLevel, maxDepth, sumDepth, numLeaves);
    }
}

size_t countNewickTips(const std::string& s) {
    std::string t = trimNewickToken(s);
    if (t.empty()) {
        return 0;
    }
    if (t[0] != '(') {
        return 1;
    }
    size_t close = findMatchingCloseParen(t, 0);
    if (close == std::string::npos) {
        return 0;
    }
    std::string inner = t.substr(1, close - 1);
    std::vector<std::string> parts;
    splitTopLevelCommasNewick(inner, parts);
    size_t c = 0;
    for (const auto& p : parts) {
        c += countNewickTips(p);
    }
    return c;
}

Node* parseNewickSubtree(Tree* tree, const std::string& s, Node* parent) {
    std::string t = trimNewickToken(s);
    if (t.empty()) {
        fprintf(stderr, "ERROR: empty Newick subtree\n");
        exit(1);
    }
    if (t[0] == '(') {
        size_t close = findMatchingCloseParen(t, 0);
        if (close == std::string::npos) {
            fprintf(stderr, "ERROR: unmatched '(' in Newick\n");
            exit(1);
        }
        std::string inner = t.substr(1, close - 1);
        std::string rest = trimNewickToken(t.substr(close + 1));
        auto labelLen = parseNameAndBranchLength(rest);
        (void)labelLen.first;
        double edgeBl = labelLen.second;

        std::vector<std::string> childStrs;
        splitTopLevelCommasNewick(inner, childStrs);
        if (childStrs.empty()) {
            fprintf(stderr, "ERROR: internal Newick node with no children\n");
            exit(1);
        }
        int idx = static_cast<int>(tree->m_currInternalNode + 1);
        std::string nid = tree->newInternalNodeId();
        if (labelLen.first.rfind("__mlipper_id_", 0) == 0) {
            nid = labelLen.first;
        }
        Node* inode = nullptr;
        if (parent == nullptr) {
            inode = new Node(nid, edgeBl);
        } else {
            inode = new Node(nid, parent, edgeBl);
        }
        inode->idx = idx;
        tree->allNodes[nid] = inode;
        for (const auto& ch : childStrs) {
            if (ch.empty()) {
                fprintf(stderr, "ERROR: empty child in Newick\n");
                exit(1);
            }
            parseNewickSubtree(tree, ch, inode);
        }
        return inode;
    }

    auto nameLen = parseNameAndBranchLength(t);
    if (nameLen.first.empty()) {
        fprintf(stderr, "ERROR: leaf with empty name in Newick\n");
        exit(1);
    }
    Node* leaf = nullptr;
    if (parent == nullptr) {
        leaf = new Node(nameLen.first, nameLen.second);
    } else {
        leaf = new Node(nameLen.first, parent, nameLen.second);
    }
    leaf->idx = static_cast<int>(tree->m_numLeafID++);
    tree->allNodes[nameLen.first] = leaf;
    return leaf;
}

void unlinkFromParent(Node* child) {
    if (!child->parent) {
        return;
    }
    auto& v = child->parent->children;
    v.erase(std::remove(v.begin(), v.end(), child), v.end());
    child->parent = nullptr;
}

void linkChild(Node* parent, Node* child, double bl) {
    child->parent = parent;
    child->bl = bl;
    parent->children.push_back(child);
}

void fixLevelsFrom(Node* node, size_t level) {
    node->level = level;
    for (Node* ch : node->children) {
        fixLevelsFrom(ch, level + 1);
    }
}

void resolvePolytomyAtNode(Node* p, Tree* tree, std::mt19937& rng) {
    while (p->children.size() > 2) {
        const size_t k = p->children.size();
        std::uniform_int_distribution<size_t> dist(0, k - 1);
        size_t i = dist(rng);
        size_t j = dist(rng);
        while (j == i) {
            j = dist(rng);
        }
        Node* a = p->children[i];
        Node* b = p->children[j];
        const double blA = a->bl;
        const double blB = b->bl;
        unlinkFromParent(a);
        unlinkFromParent(b);
        const int idx = static_cast<int>(tree->m_currInternalNode + 1);
        const std::string nid = tree->newInternalNodeId();
        Node* n = new Node(nid, p, 0.0);
        n->idx = idx;
        tree->allNodes[nid] = n;
        linkChild(n, a, blA);
        linkChild(n, b, blB);
    }
}

void postOrderResolvePolytomies(Node* node, Tree* tree, std::mt19937& rng, bool enabled) {
    if (!node) {
        return;
    }
    for (Node* c : node->children) {
        postOrderResolvePolytomies(c, tree, rng, enabled);
    }
    if (enabled && node->children.size() > 2) {
        resolvePolytomyAtNode(node, tree, rng);
    }
}

void rebuildAllNodesMap(Tree* tree, Node* r) {
    tree->allNodes.clear();
    if (!r) {
        return;
    }
    std::vector<Node*> vec;
    tree->dfsExpansion(r, vec);
    for (Node* n : vec) {
        tree->allNodes[n->name] = n;
    }
}

}  // namespace

Tree::Tree(Node* node) : root(node) {
    if (node) {
        std::vector<Node*> vec;
        Tree::dfsExpansion(node, vec);
        for (auto* n : vec) {
            allNodes[n->name] = n;
        }
    }
}

Tree::Tree(std::string newickString, size_t totalLeaves, bool randomResolvePolytomies,
           uint64_t polytomyRngSeed) {
    newickString = stripString(newickString);
    while (!newickString.empty() && newickString.back() == ';') {
        newickString.pop_back();
    }
    newickString = trimNewickToken(newickString);

    if (newickString.empty()) {
        fprintf(stderr, "WARNING: Tree found empty!\n");
        root = nullptr;
        return;
    }

    m_numLeafID = 0;
    if (totalLeaves == 0) {
        const size_t nTips = countNewickTips(newickString);
        m_currInternalNode = nTips > 0 ? nTips - 1 : 0;
    } else {
        m_currInternalNode = totalLeaves - 1;
    }

    Node* treeRoot = parseNewickSubtree(this, newickString, nullptr);

    std::mt19937 rng;
    if (polytomyRngSeed == 0) {
        std::random_device rd;
        std::seed_seq seq{rd(), rd(), rd(), rd()};
        rng.seed(seq);
    } else {
        rng.seed(static_cast<std::mt19937::result_type>(polytomyRngSeed));
    }
    postOrderResolvePolytomies(treeRoot, this, rng, randomResolvePolytomies);
    if (treeRoot != nullptr) {
        fixLevelsFrom(treeRoot, 1);
    }
    rebuildAllNodesMap(this, treeRoot);

    size_t nLeaf = 0;
    size_t sumDepth = 0;
    size_t maxD = 0;
    if (treeRoot != nullptr) {
        treeDepthStats(treeRoot, treeRoot->level, maxD, sumDepth, nLeaf);
    }
    m_numLeaves = nLeaf;
    m_maxDepth = maxD;
    m_meanDepth = nLeaf > 0 ? static_cast<float>(sumDepth) / static_cast<float>(nLeaf) : 0.f;

    if (treeRoot == nullptr) {
        fprintf(stderr, "WARNING: Tree found empty!\n");
    } else {
        treeRoot->bl = 0;
    }
    root = treeRoot;
}

Tree::~Tree() {
    for (auto n: this->allNodes) {
        delete n.second;
    }
    this->allNodes.clear();
    this->root = nullptr;
}

// =============================================================================
// Unrooted Tree Implementation
// =============================================================================

namespace {
void buildSubtreeFromNode(Node* node, UnrootedNode* unrootedParent, double blToParent, int edgeId,
                          std::unordered_map<std::string, UnrootedNode*>& unrootedNodes,
                          size_t& leafId, size_t& internalId, int& nextEdgeId) {
    UnrootedNode* u = new UnrootedNode();
    u->name = node->name;
    u->idx = node->children.empty() ? static_cast<int>(leafId++) : static_cast<int>(internalId++);
    unrootedNodes[node->name] = u;
    if (unrootedParent) {
        u->neighbors.push_back({unrootedParent, blToParent, edgeId});
        unrootedParent->neighbors.push_back({u, blToParent, edgeId});
    }
    for (auto* ch : node->children) {
        int eid = nextEdgeId++;
        buildSubtreeFromNode(ch, u, ch->bl, eid, unrootedNodes, leafId, internalId, nextEdgeId);
    }
}

void buildAdjacencyFromTreeWithoutRoot(Node* root,
                                       std::unordered_map<std::string, UnrootedNode*>& unrootedNodes,
                                       size_t& leafId, size_t& internalId,
                                       size_t& internalNodeNameCounter, int& nextEdgeId) {
    const auto& children = root->children;
    if (children.empty()) return;
    if (children.size() == 1) {
        buildSubtreeFromNode(children[0], nullptr, 0.0, -1, unrootedNodes, leafId, internalId, nextEdgeId);
        return;
    }
    if (children.size() == 2) {
        Node* c1 = children[0];
        Node* c2 = children[1];
        double combinedBl = c1->bl + c2->bl;
        int eid = nextEdgeId++;
        UnrootedNode* u1 = new UnrootedNode();
        u1->name = c1->name;
        u1->idx = c1->children.empty() ? static_cast<int>(leafId++) : static_cast<int>(internalId++);
        unrootedNodes[c1->name] = u1;
        UnrootedNode* u2 = new UnrootedNode();
        u2->name = c2->name;
        u2->idx = c2->children.empty() ? static_cast<int>(leafId++) : static_cast<int>(internalId++);
        unrootedNodes[c2->name] = u2;
        u1->neighbors.push_back({u2, combinedBl, eid});
        u2->neighbors.push_back({u1, combinedBl, eid});
        for (auto* ch : c1->children) {
            int ceid = nextEdgeId++;
            buildSubtreeFromNode(ch, u1, ch->bl, ceid, unrootedNodes, leafId, internalId, nextEdgeId);
        }
        for (auto* ch : c2->children) {
            int ceid = nextEdgeId++;
            buildSubtreeFromNode(ch, u2, ch->bl, ceid, unrootedNodes, leafId, internalId, nextEdgeId);
        }
    } else {
        std::string internalIdStr = "node_" + std::to_string(++internalNodeNameCounter);
        UnrootedNode* center = new UnrootedNode();
        center->name = internalIdStr;
        center->idx = static_cast<int>(internalId++);
        unrootedNodes[internalIdStr] = center;
        for (auto* ch : children) {
            int eid = nextEdgeId++;
            UnrootedNode* u = new UnrootedNode();
            u->name = ch->name;
            u->idx = ch->children.empty() ? static_cast<int>(leafId++) : static_cast<int>(internalId++);
            unrootedNodes[ch->name] = u;
            u->neighbors.push_back({center, ch->bl, eid});
            center->neighbors.push_back({u, ch->bl, eid});
            for (auto* grandch : ch->children) {
                int geid = nextEdgeId++;
                buildSubtreeFromNode(grandch, u, grandch->bl, geid, unrootedNodes, leafId, internalId, nextEdgeId);
            }
        }
    }
}

double getBranchLength(UnrootedNode* from, UnrootedNode* to) {
    for (const auto& nb : from->neighbors) {
        if (nb.node == to) return nb.length;
    }
    return -1.0;
}

void buildRootedSubtree(UnrootedNode* u, UnrootedNode* exclude, Node* parent,
                       double bl, std::unordered_map<std::string, Node*>& allNodes) {
    Node* n = (parent) ? new Node(u->name, parent, bl) : new Node(u->name, bl);
    n->idx = u->idx;
    allNodes[u->name] = n;
    for (const auto& nb : u->neighbors) {
        if (nb.node != exclude) {
            buildRootedSubtree(nb.node, u, n, nb.length, allNodes);
        }
    }
}

std::string quoteNewickName(const std::string& name) {
    if (name.empty()) return name;
    static const char* special = " \t\r\n()[]':;,_";
    bool needsQuote = (name.find_first_of(special) != std::string::npos);
    if (!needsQuote) return name;
    std::string result = "'";
    for (char c : name) {
        if (c == '\'') result += "''";
        else result += c;
    }
    result += "'";
    return result;
}
}  // namespace

IUnrootedTree* IUnrootedTree::fromNewick(const std::string& newick, size_t totalLeaves) {
    return new UnrootedTree(newick, totalLeaves);
}

UnrootedTree::UnrootedTree(const std::string& newick, size_t totalLeaves) {
    Tree tmp(newick, totalLeaves);
    size_t leafId = 0, internalId = tmp.m_numLeaves;
    size_t internalNodeNameCounter = tmp.m_currInternalNode;
    int nextEdgeId = 0;
    buildAdjacencyFromTreeWithoutRoot(tmp.root, m_nodes, leafId, internalId, internalNodeNameCounter, nextEdgeId);
    m_numLeaves = tmp.m_numLeaves;
    m_numEdges = 0;
    for (const auto& p : m_nodes) {
        m_numEdges += p.second->neighbors.size();
    }
    m_numEdges /= 2;
    m_currInternalNode = internalNodeNameCounter;
    collapseZeroLengthInternalEdges();
}

namespace {
/** Get bipartition signature for node (smallest leaf set when removing node, excluding tip). */
std::vector<std::string> getBipartitionSignature(UnrootedNode* node, UnrootedNode* tipNode,
                                                 const UnrootedTree* tree) {
    std::vector<std::set<std::string>> parts;
    auto getLeaves = [](UnrootedNode* start, UnrootedNode* exclude, const UnrootedTree* t) {
        std::set<std::string> leaves;
        if (!start || !t) return leaves;
        std::stack<std::pair<UnrootedNode*, UnrootedNode*>> stk;
        stk.push({start, exclude});
        while (!stk.empty()) {
            UnrootedNode* u = stk.top().first;
            UnrootedNode* from = stk.top().second;
            stk.pop();
            if (u->is_leaf()) { leaves.insert(u->name); continue; }
            for (const auto& nb : u->neighbors) {
                if (nb.node != from) stk.push({nb.node, u});
            }
        }
        return leaves;
    };
    for (const auto& nb : node->neighbors) {
        if (nb.node == tipNode) continue;
        std::set<std::string> part = getLeaves(nb.node, node, tree);
        part.erase(tipNode->name);
        if (!part.empty()) parts.push_back(std::move(part));
    }
    if (parts.empty()) return {};
    size_t smallestIdx = 0;
    for (size_t i = 1; i < parts.size(); ++i)
        if (parts[i].size() < parts[smallestIdx].size()) smallestIdx = i;
    std::vector<std::string> v(parts[smallestIdx].begin(), parts[smallestIdx].end());
    std::sort(v.begin(), v.end());
    return v;
}

UnrootedNode* findNodeByBipartition(const UnrootedTree* tree, const std::vector<std::string>& sig,
                                   UnrootedNode* excludeTip) {
    for (const auto& p : tree->getNodes()) {
        UnrootedNode* u = p.second;
        if (u == excludeTip || u->is_leaf()) continue;
        if (getBipartitionSignature(u, excludeTip, tree) == sig) return u;
    }
    return nullptr;
}

std::string trimName(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

UnrootedNode* findTipInTree(const UnrootedTree* tree, const std::string& tipName) {
    if (!tree) return nullptr;
    UnrootedNode* n = tree->getNode(tipName);
    if (n && n->is_leaf()) return n;
    std::string trimmed = trimName(tipName);
    if (trimmed.empty()) return nullptr;
    for (const auto& p : tree->getNodes()) {
        UnrootedNode* u = p.second;
        if (!u->is_leaf()) continue;
        if (trimName(u->name) == trimmed) return u;
    }
    return nullptr;
}
}  // namespace

UnrootedTree::UnrootedTree(const std::string& newick, const UnrootedTree* reference,
                          const std::string& tipName, size_t totalLeaves) {
    if (!reference) {
        new (this) UnrootedTree(newick, totalLeaves);
        return;
    }
    UnrootedTree t2Temp(newick, totalLeaves);
    UnrootedNode* tip2 = findTipInTree(&t2Temp, tipName);
    UnrootedNode* tipRef = findTipInTree(reference, tipName);
    if (!tip2 || !tipRef) {
        new (this) UnrootedTree(newick, totalLeaves);
        return;
    }

    m_numLeaves = reference->m_numLeaves;
    m_numEdges = reference->m_numEdges;
    m_currInternalNode = reference->m_currInternalNode;

    std::unordered_map<UnrootedNode*, UnrootedNode*> oldToNew;
    for (const auto& p : reference->getNodes()) {
        UnrootedNode* u = new UnrootedNode();
        u->name = p.second->name;
        u->idx = p.second->idx;
        u->bl = p.second->bl;
        m_nodes[u->name] = u;
        oldToNew[p.second] = u;
    }
    for (const auto& p : reference->getNodes()) {
        UnrootedNode* uOld = p.second;
        UnrootedNode* uNew = oldToNew[uOld];
        for (const auto& nb : uOld->neighbors) {
            uNew->neighbors.push_back({oldToNew[nb.node], nb.length, nb.edge_id});
        }
    }

    UnrootedNode* attach2Temp = tip2->neighbors.empty() ? nullptr : tip2->neighbors[0].node;
    if (!attach2Temp) {
        reindexEdges();
        return;
    }

    std::vector<std::string> sig2 = getBipartitionSignature(attach2Temp, tip2, &t2Temp);
    if (sig2.empty()) {
        reindexEdges();
        return;
    }

    UnrootedNode* attach2Mapped = findNodeByBipartition(reference, sig2, tipRef);
    if (!attach2Mapped) {
        reindexEdges();
        return;
    }

    UnrootedNode* tipNew = oldToNew[tipRef];
    UnrootedNode* attach1New = tipNew->neighbors.empty() ? nullptr : tipNew->neighbors[0].node;
    UnrootedNode* attach2New = oldToNew[attach2Mapped];
    if (attach1New && attach1New != attach2New) {
        auto removeEdge = [](UnrootedNode* a, UnrootedNode* b) {
            auto& nbs = a->neighbors;
            nbs.erase(std::remove_if(nbs.begin(), nbs.end(),
                [b](const UnrootedNeighbor& nb) { return nb.node == b; }), nbs.end());
        };
        removeEdge(tipNew, attach1New);
        removeEdge(attach1New, tipNew);
        double bl = 0.0;
        for (const auto& nb : tip2->neighbors)
            if (nb.node == attach2Temp) { bl = nb.length; break; }
        int eid = static_cast<int>(m_numEdges++);
        tipNew->neighbors.push_back({attach2New, bl, eid});
        attach2New->neighbors.push_back({tipNew, bl, eid});
    }
    reindexEdges();
}

UnrootedTree::~UnrootedTree() {
    for (auto& p : m_nodes) {
        delete p.second;
    }
    m_nodes.clear();
}

UnrootedNode* UnrootedTree::getNode(const std::string& name) const {
    auto it = m_nodes.find(name);
    return (it != m_nodes.end()) ? it->second : nullptr;
}

const std::vector<UnrootedNeighbor>& UnrootedTree::getNeighbors(UnrootedNode* node) const {
    static const std::vector<UnrootedNeighbor> empty;
    if (!node) return empty;
    return node->neighbors;
}

int UnrootedTree::getEdgeId(UnrootedNode* from, UnrootedNode* to) const {
    if (!from || !to) return -1;
    for (const auto& nb : from->neighbors) {
        if (nb.node == to) return nb.edge_id;
    }
    return -1;
}

int UnrootedTree::getEdgeId(UnrootedNode* node, size_t neighborIndex) const {
    if (!node || neighborIndex >= node->neighbors.size()) return -1;
    return node->neighbors[neighborIndex].edge_id;
}

void UnrootedTree::rootBetween(const std::string& nodeName1, const std::string& nodeName2) {
    UnrootedNode* u1 = getNode(nodeName1);
    UnrootedNode* u2 = getNode(nodeName2);
    if (!u1 || !u2) return;
    double bl = getBranchLength(u1, u2);
    if (bl < 0) return;
    double half = bl / 2.0;

    // Remove the existing edge between u1 and u2
    auto removeNeighbor = [](UnrootedNode* from, UnrootedNode* to) {
        auto& nbs = from->neighbors;
        nbs.erase(std::remove_if(nbs.begin(), nbs.end(),
                                 [to](const UnrootedNeighbor& nb) { return nb.node == to; }),
                  nbs.end());
    };
    removeNeighbor(u1, u2);
    removeNeighbor(u2, u1);

    // Add new root node and two new edges
    std::string rootId = "node_" + std::to_string(++m_currInternalNode);
    UnrootedNode* unrootedRoot = new UnrootedNode();
    unrootedRoot->name = rootId;
    unrootedRoot->idx = static_cast<int>(m_numLeaves + m_numEdges);
    m_nodes[rootId] = unrootedRoot;
    m_root = unrootedRoot;

    int edgeId1 = static_cast<int>(m_numEdges);
    int edgeId2 = static_cast<int>(m_numEdges + 1);
    unrootedRoot->neighbors.push_back({u1, half, edgeId1});
    unrootedRoot->neighbors.push_back({u2, half, edgeId2});
    u1->neighbors.push_back({unrootedRoot, half, edgeId1});
    u2->neighbors.push_back({unrootedRoot, half, edgeId2});

    m_numEdges += 1;  // removed 1, added 2 -> net +1
    reindexEdges();
}

void UnrootedTree::rootAt(const std::string& nodeName) {
    UnrootedNode* u = getNode(nodeName);
    if (!u || u->neighbors.empty()) return;
    // Root on first edge (arbitrary choice)
    UnrootedNode* other = u->neighbors[0].node;
    double bl = u->neighbors[0].length;
    double half = bl / 2.0;

    // Remove the existing edge between u and other
    auto removeNeighbor = [](UnrootedNode* from, UnrootedNode* to) {
        auto& nbs = from->neighbors;
        nbs.erase(std::remove_if(nbs.begin(), nbs.end(),
                                 [to](const UnrootedNeighbor& nb) { return nb.node == to; }),
                  nbs.end());
    };
    removeNeighbor(u, other);
    removeNeighbor(other, u);

    // Add new root node and two new edges
    std::string rootId = "node_" + std::to_string(++m_currInternalNode);
    UnrootedNode* unrootedRoot = new UnrootedNode();
    unrootedRoot->name = rootId;
    unrootedRoot->idx = static_cast<int>(m_numLeaves + m_numEdges);
    m_nodes[rootId] = unrootedRoot;

    int edgeId1 = static_cast<int>(m_numEdges);
    int edgeId2 = static_cast<int>(m_numEdges + 1);
    unrootedRoot->neighbors.push_back({u, half, edgeId1});
    unrootedRoot->neighbors.push_back({other, half, edgeId2});
    u->neighbors.push_back({unrootedRoot, half, edgeId1});
    other->neighbors.push_back({unrootedRoot, half, edgeId2});

    m_numEdges += 1;  // removed 1, added 2 -> net +1
    reindexEdges();
}

void UnrootedTree::reindexEdges() {
    std::map<std::pair<uintptr_t, uintptr_t>, int> edgeToNewId;
    int nextId = 0;
    for (const auto& p : m_nodes) {
        UnrootedNode* u = p.second;
        for (const auto& nb : u->neighbors) {
            UnrootedNode* v = nb.node;
            uintptr_t a = reinterpret_cast<uintptr_t>(u);
            uintptr_t b = reinterpret_cast<uintptr_t>(v);
            auto key = std::make_pair(std::min(a, b), std::max(a, b));
            if (edgeToNewId.find(key) == edgeToNewId.end()) {
                edgeToNewId[key] = nextId++;
            }
        }
    }
    for (const auto& p : m_nodes) {
        UnrootedNode* u = p.second;
        for (auto& nb : u->neighbors) {
            UnrootedNode* v = nb.node;
            uintptr_t a = reinterpret_cast<uintptr_t>(u);
            uintptr_t b = reinterpret_cast<uintptr_t>(v);
            auto key = std::make_pair(std::min(a, b), std::max(a, b));
            nb.edge_id = edgeToNewId[key];
        }
    }
    m_numEdges = static_cast<size_t>(nextId);
}

namespace {

void collectLeavesBeyond(UnrootedNode* cur, UnrootedNode* from, std::set<std::string>& out) {
    if (!cur) return;
    if (cur->is_leaf()) {
        out.insert(cur->name);
        return;
    }
    for (const auto& nb : cur->neighbors) {
        if (nb.node == from) continue;
        collectLeavesBeyond(nb.node, cur, out);
    }
}

std::string joinSortedLeafSet(const std::set<std::string>& s) {
    std::string r;
    for (const auto& x : s) {
        if (!r.empty()) r += ',';
        r += x;
    }
    return r;
}

/** Canonical key for a non-trivial split; empty if either side has < 2 leaves. */
std::string nonTrivialSplitKey(const std::set<std::string>& Su, const std::set<std::string>& allLeaves) {
    std::set<std::string> Sv;
    for (const auto& x : allLeaves) {
        if (Su.find(x) == Su.end()) Sv.insert(x);
    }
    if (Su.size() < 2 || Sv.size() < 2) return "";

    const std::set<std::string>* a = &Su;
    const std::set<std::string>* b = &Sv;
    if (Su.size() > Sv.size()) {
        a = &Sv;
        b = &Su;
    } else if (Su.size() == Sv.size() && joinSortedLeafSet(Su) > joinSortedLeafSet(Sv)) {
        a = &Sv;
        b = &Su;
    }
    return joinSortedLeafSet(*a) + std::string("||") + joinSortedLeafSet(*b);
}

}  // namespace

std::set<std::string> collectNonTrivialSplitKeys(const UnrootedTree& tree) {
    std::set<std::string> allLeaves;
    for (const auto& p : tree.getNodes()) {
        if (p.second->is_leaf()) allLeaves.insert(p.second->name);
    }
    std::set<std::string> keys;
    std::set<std::pair<uintptr_t, uintptr_t>> seen;
    for (const auto& p : tree.getNodes()) {
        UnrootedNode* u = p.second;
        for (const auto& nb : u->neighbors) {
            UnrootedNode* v = nb.node;
            uintptr_t a = reinterpret_cast<uintptr_t>(u);
            uintptr_t b = reinterpret_cast<uintptr_t>(v);
            if (a > b) std::swap(a, b);
            if (!seen.insert({a, b}).second) continue;
            std::set<std::string> Su;
            collectLeavesBeyond(u, v, Su);
            std::string key = nonTrivialSplitKey(Su, allLeaves);
            if (!key.empty()) keys.insert(std::move(key));
        }
    }
    return keys;
}

void UnrootedTree::contractUndirectedEdge(UnrootedNode* u, UnrootedNode* v) {
    if (!u || !v) return;
    double len_uv = getBranchLength(u, v);
    if (len_uv < 0.0) return;

    auto removeNeighbor = [](UnrootedNode* from, UnrootedNode* to) {
        auto& nbs = from->neighbors;
        nbs.erase(std::remove_if(nbs.begin(), nbs.end(),
                                 [to](const UnrootedNeighbor& nb) { return nb.node == to; }),
                  nbs.end());
    };

    removeNeighbor(u, v);
    removeNeighbor(v, u);

    std::vector<UnrootedNeighbor> vNbs = v->neighbors;
    for (const auto& nb : vNbs) {
        UnrootedNode* w = nb.node;
        removeNeighbor(w, v);
        double newLen = len_uv + nb.length;
        int eid = static_cast<int>(m_numEdges++);
        u->neighbors.push_back({w, newLen, eid});
        w->neighbors.push_back({u, newLen, eid});
    }

    m_nodes.erase(v->name);
    delete v;
    reindexEdges();
}

void UnrootedTree::collapseZeroLengthInternalEdges(double lengthEps) {
    bool changed = true;
    while (changed) {
        changed = false;
        std::set<std::pair<uintptr_t, uintptr_t>> seen;
        for (const auto& p : m_nodes) {
            UnrootedNode* u = p.second;
            for (const auto& nb : u->neighbors) {
                UnrootedNode* v = nb.node;
                if (nb.length < 0.0 || nb.length > lengthEps) continue;
                if (u->is_leaf() || v->is_leaf()) continue;
                uintptr_t a = reinterpret_cast<uintptr_t>(u);
                uintptr_t b = reinterpret_cast<uintptr_t>(v);
                if (a > b) std::swap(a, b);
                if (!seen.insert({a, b}).second) continue;
                contractUndirectedEdge(u, v);
                changed = true;
                break;
            }
            if (changed) break;
        }
    }
}

void UnrootedTree::collapseToSplitsDisplayedIn(const std::set<std::string>& allowedNonTrivialSplitKeys) {
    bool changed = true;
    while (changed) {
        changed = false;
        std::set<std::string> allLeaves;
        for (const auto& p : m_nodes) {
            if (p.second->is_leaf()) allLeaves.insert(p.second->name);
        }
        std::set<std::pair<uintptr_t, uintptr_t>> seen;
        for (const auto& p : m_nodes) {
            UnrootedNode* u = p.second;
            for (const auto& nb : u->neighbors) {
                UnrootedNode* v = nb.node;
                uintptr_t a = reinterpret_cast<uintptr_t>(u);
                uintptr_t b = reinterpret_cast<uintptr_t>(v);
                if (a > b) std::swap(a, b);
                if (!seen.insert({a, b}).second) continue;
                std::set<std::string> Su;
                collectLeavesBeyond(u, v, Su);
                std::string key = nonTrivialSplitKey(Su, allLeaves);
                if (key.empty()) continue;
                if (allowedNonTrivialSplitKeys.find(key) != allowedNonTrivialSplitKeys.end()) continue;
                contractUndirectedEdge(u, v);
                changed = true;
                break;
            }
            if (changed) break;
        }
    }
}

void UnrootedTree::collapseToSplitsDisplayedIn(const UnrootedTree& other) {
    collapseToSplitsDisplayedIn(collectNonTrivialSplitKeys(other));
}

std::string UnrootedTree::toNewick() const {
    if (m_nodes.empty()) return ";";
    UnrootedNode* start = nullptr;
    for (const auto& p : m_nodes) {
        if (p.second->degree() > 1) {
            start = p.second;
            break;
        }
    }
    if (!start) start = m_nodes.begin()->second;
    std::function<std::string(UnrootedNode*, UnrootedNode*)> recurse = [&](UnrootedNode* n, UnrootedNode* from) {
        if (!n) return std::string();
        if (n->is_leaf()) {
            std::string s = quoteNewickName(n->name);
            if (from) {
                double bl = getBranchLength(n, from);
                if (bl >= 0) s += ":" + std::to_string(bl);
            }
            return s;
        }
        std::string inner;
        for (const auto& nb : n->neighbors) {
            if (nb.node != from) {
                if (!inner.empty()) inner += ",";
                inner += recurse(nb.node, n);
            }
        }
        std::string blStr;
        if (from) {
            double bl = getBranchLength(n, from);
            if (bl >= 0) blStr = ":" + std::to_string(bl);
        }
        std::string namePart;
        if (!n->name.empty() && (n->name.size() < 5 || n->name.substr(0, 5) != "node_")) {
            namePart = quoteNewickName(n->name);
        }
        return "(" + inner + ")" + namePart + blStr;
    };
    return recurse(start, nullptr) + ";";
}

std::vector<UnrootedNode*> UnrootedTree::bfs(UnrootedNode* start) const {
    std::vector<UnrootedNode*> result;
    if (!start) return result;
    std::queue<UnrootedNode*> q;
    std::unordered_set<UnrootedNode*> visited;
    q.push(start);
    visited.insert(start);
    while (!q.empty()) {
        UnrootedNode* u = q.front();
        q.pop();
        result.push_back(u);
        for (const auto& nb : u->neighbors) {
            if (visited.find(nb.node) == visited.end()) {
                visited.insert(nb.node);
                q.push(nb.node);
            }
        }
    }
    return result;
}

std::vector<UnrootedNode*> UnrootedTree::dfs(UnrootedNode* start) const {
    std::vector<UnrootedNode*> result;
    if (!start) return result;
    std::stack<UnrootedNode*> stk;
    std::unordered_set<UnrootedNode*> visited;
    stk.push(start);
    visited.insert(start);
    while (!stk.empty()) {
        UnrootedNode* u = stk.top();
        stk.pop();
        result.push_back(u);
        for (const auto& nb : u->neighbors) {
            if (visited.find(nb.node) == visited.end()) {
                visited.insert(nb.node);
                stk.push(nb.node);
            }
        }
    }
    return result;
}

std::vector<UnrootedNode*> UnrootedTree::nodesBetween(const std::string& nodeName1,
                                                    const std::string& nodeName2) const {
    std::vector<UnrootedNode*> path;
    std::vector<UnrootedEdge> edges = edgesBetween(nodeName1, nodeName2);
    for (auto &e: edges) {
        path.push_back(e.from);
    }
    return path;
}

std::vector<UnrootedEdge> UnrootedTree::edgesBetween(const std::string& nodeName1,
                                                    const std::string& nodeName2) const {
    std::vector<UnrootedEdge> path;
    UnrootedNode* u1 = getNode(nodeName1);
    UnrootedNode* u2 = getNode(nodeName2);
    if (!u1 || !u2 || u1 == u2) return path;

    struct ParentInfo { UnrootedNode* prev; double len; int edge_id; };
    std::unordered_map<UnrootedNode*, ParentInfo> parent;
    std::queue<UnrootedNode*> q;
    std::unordered_set<UnrootedNode*> visited;
    q.push(u1);
    visited.insert(u1);
    parent[u1] = {nullptr, 0.0, -1};
    bool found = false;
    while (!q.empty() && !found) {
        UnrootedNode* u = q.front();
        q.pop();
        for (const auto& nb : u->neighbors) {
            if (visited.find(nb.node) == visited.end()) {
                visited.insert(nb.node);
                parent[nb.node] = {u, nb.length, nb.edge_id};
                if (nb.node == u2) {
                    found = true;
                    break;
                }
                q.push(nb.node);
            }
        }
    }
    if (!found) return path;

    UnrootedNode* curr = u2;
    while (parent[curr].prev != nullptr) {
        const auto& info = parent[curr];
        path.push_back({info.prev, curr, info.len, info.edge_id});
        curr = info.prev;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

namespace {
/** Find a leaf node by name, with flexible matching (exact match, or trim both sides). */
UnrootedNode* findTipNode(const UnrootedTree* tree, const std::string& tipName) {
    if (!tree) return nullptr;
    UnrootedNode* n = tree->getNode(tipName);
    if (n && n->is_leaf()) return n;
    std::string trimmed = trimName(tipName);
    if (trimmed.empty()) return nullptr;
    for (const auto& p : tree->getNodes()) {
        UnrootedNode* u = p.second;
        if (!u->is_leaf()) continue;
        if (trimName(u->name) == trimmed) return u;
    }
    return nullptr;
}

/** BFS distance from start to target, excluding excludeNode from traversal. */
int bfsDistanceExcluding(UnrootedNode* start, UnrootedNode* target, UnrootedNode* excludeNode,
                         const UnrootedTree* tree) {
    if (!start || !target || !tree) return -1;
    if (start == target) return 0;
    std::queue<std::pair<UnrootedNode*, int>> q;
    std::unordered_set<UnrootedNode*> visited;
    q.push({start, 0});
    visited.insert(start);
    if (excludeNode) visited.insert(excludeNode);
    while (!q.empty()) {
        UnrootedNode* u = q.front().first;
        int d = q.front().second;
        q.pop();
        for (const auto& nb : u->neighbors) {
            if (nb.node == excludeNode) continue;
            if (visited.find(nb.node) != visited.end()) continue;
            if (nb.node == target) return d + 1;
            visited.insert(nb.node);
            q.push({nb.node, d + 1});
        }
    }
    return -1;
}
}  // namespace

int computeTipPlacementEdgeDistance(const UnrootedTree* t1, const UnrootedTree* t2,
                                    const std::string& tipName) {
    if (!t1 || !t2) return -1;
    UnrootedNode* tip1 = findTipNode(t1, tipName);
    UnrootedNode* tip2 = findTipNode(t2, tipName);
    if (!tip1 || !tip2) return -1;
    if (!tip1->is_leaf() || !tip2->is_leaf()) return -1;
    if (tip1->neighbors.empty() || tip2->neighbors.empty()) return -1;

    UnrootedNode* attach1 = tip1->neighbors[0].node;
    UnrootedNode* attach2 = tip2->neighbors[0].node;
    UnrootedNode* attach2InT1 = t1->getNode(attach2->name);
    if (!attach2InT1) return -1;

    return bfsDistanceExcluding(attach1, attach2InT1, tip1, t1);
}

void dfsExpansionToUpdateClosestIds(Node* node) {
    if (node->children.size() == 0) return;

    for (auto n: node->children) {
        dfsExpansionToUpdateClosestIds(n);
    }
 
    std::vector<std::pair<int, double>> merged;
    for (auto n: node->children) {
        merged.insert(merged.end(), n->closestIds.begin(), n->closestIds.end());
    }
    // keep top 10 closest ids
    std::sort(merged.begin(), merged.end(),
              [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                  return a.second < b.second;
              });
    if (merged.size() > 10) {
        merged.resize(10);
    }
    node->closestIds = merged;
}

void Tree::updateClosestIds(std::unordered_map<std::string, std::vector<std::pair<int, double>>>& closestIds) {
    std::vector<Node*> traversal;
    dfsExpansion(this->root, traversal);
    for (auto n: traversal) {
        if (n->is_leaf()) {
            n->closestIds = closestIds[n->name];
        }
    }

    // dfs to update internal nodes (keep top 10 closest ids)
    dfsExpansionToUpdateClosestIds(this->root);
    

    // bfs to learn closest ids for parent and merge if needed
    std::queue<Node*> q;
    for (auto child: this->root->children) q.push(child);

    Node* n = this->root;
    std::cout << n->name << ": ";
    for (auto p: n->closestIds) {
        std::cout << "(" << p.first << ", " << p.second << ") ";
    }
    std::cout << std::endl; 
    for (auto n: this->root->children){
        std::cout << n->name << ": ";
        for (auto p: n->closestIds) {
            std::cout << "(" << p.first << ", " << p.second << ") ";
        }
        std::cout << std::endl;
    }

    while (!q.empty()) {
        Node* curr = q.front();
        q.pop();
        for (auto child: curr->children) {
            q.push(child);
        }
        if (curr->parent != nullptr) {
            std::vector<std::pair<int, double>> merged;
            merged.insert(merged.end(), curr->closestIds.begin(), curr->closestIds.end());
            merged.insert(merged.end(), curr->parent->closestIds.begin(), curr->parent->closestIds.end());
            // keep top 10 closest ids
            std::sort(merged.begin(), merged.end(),
                      [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                          return a.second < b.second;
                      });
            if (merged.size() > 10) {
                merged.resize(10);
            }
            curr->closestIds = merged;
        }
    }

    // for (auto n: traversal){
    //     std::cout << n->name << ": ";
    //     for (auto p: n->closestIds) {
    //         std::cout << "(" << p.first << ", " << p.second << ") ";
    //     }
    //     std::cout << std::endl; 
    // }

}

struct placementMetadata_ {
    double distBettwoNodes;
    double longestBl;
    int longestTip;
    std::pair<std::string, std::string> twoNodesitConnects;
    placementMetadata_() : distBettwoNodes(0), longestBl(0), longestTip(-1) {}
    placementMetadata_(double d, double bl, int tip, std::pair<std::string, std::string> conn = std::make_pair("", ""))
        : distBettwoNodes(d), longestBl(bl), longestTip(tip), twoNodesitConnects(conn) {}
};

Node* lcs(Node* node1, Node* node2){
    if (node1 == nullptr || node2 == nullptr) return nullptr;
    std::vector<Node*> node1_parents, node2_parents;
    Node* tmp = node1;
    while(tmp!=nullptr){
        node1_parents.push_back(tmp);
        tmp=tmp->parent;
    }
    tmp=node2;
    while(tmp!=nullptr){
        node2_parents.push_back(tmp);
        tmp=tmp->parent;
    }

    std::reverse(node1_parents.begin(), node1_parents.end());
    std::reverse(node2_parents.begin(), node2_parents.end());

    if (node1_parents.empty() || node2_parents.empty()) return nullptr;

    Node* lcs_result = nullptr;
    size_t itr = 0;
    while (itr < node1_parents.size() && itr < node2_parents.size() &&
           node1_parents[itr]->name == node2_parents[itr]->name) {
        lcs_result = node1_parents[itr];
        itr++;
    }

    if (lcs_result == nullptr || lcs_result->name == "") std::cerr << "Unexecpted behavior: LCS not found" << std::endl;
    return lcs_result;
}
// updateLcsNode
bool checkChildrenStatus(Node* lcs, std::unordered_map<std::string, placementMetadata_>& placementMetadata){
    bool childrenBlDone=true;
    double longestBL=0.0;
    int FarthestTip=-1;
    for (auto child: lcs->children){
        if (placementMetadata.find(child->name) != placementMetadata.end()){
            if (placementMetadata[child->name].distBettwoNodes !=2.0){
                childrenBlDone=false;
                break;
            }
        } else {
            childrenBlDone=false;
            break;
        }
        double currBl = child->bl + placementMetadata[child->name].longestBl;
        if (longestBL<currBl){
            longestBL=currBl;
            FarthestTip = placementMetadata[child->name].longestTip;
        }
    }

    if (childrenBlDone){
        placementMetadata[lcs->name]=placementMetadata_(2.0, longestBL, FarthestTip, placementMetadata[lcs->name].twoNodesitConnects);
        return true;
    } 
    return false;
}

void placementAccuracy(Tree* t, double* h_mashDist, size_t numSequences) {
    std::queue<Node*> placementOrder;
    std::unordered_map<std::string, placementMetadata_> placementMetadata;
    
    auto dfs = [&](auto&& self, Node* node, std::queue<Node*>& placementOrder)->void {
        if (node == nullptr) return;

        for (auto& child: node->children) {
            self(self, child, placementOrder);
        }
        placementOrder.push(node);
        node->bl=0.; //reset branch length to 0. -> very important
        return;
    };

    dfs(dfs, t->root, placementOrder);

    /* Initialize tree */
    Node *a, *b;
    a=placementOrder.front(); placementOrder.pop();
    b=placementOrder.front(); placementOrder.pop();
    assert(a->children.size() == 0);
    assert(b->children.size() == 0);
    
    double dist;
    if (a>b) dist=h_mashDist[a->idx*numSequences+b->idx];
    else dist=h_mashDist[b->idx*numSequences+a->idx];
    a->bl=dist/2;b->bl=dist/2;
    placementMetadata[a->name]=placementMetadata_(2.0,0,a->idx, std::make_pair(a->name, a->parent->name));
    placementMetadata[b->name]=placementMetadata_(2.0,0,b->idx, std::make_pair(b->name, b->parent->name));

    bool specialCase=false;

    // std::cout << a->name << " (" << a->bl << ")"
    //           << '\t' << placementMetadata[a->name].distBettwoNodes  
    //           << '\t' << placementMetadata[a->name].longestBl  
    //           << '\t' << placementMetadata[a->name].longestTip
    //           << '\t' << placementMetadata[a->name].twoNodesitConnects.first
    //           << '\t' << placementMetadata[a->name].twoNodesitConnects.second
    //           << std::endl;

    // std::cout << b->name << " (" << b->bl << ")"
    //           << '\t' << placementMetadata[b->name].distBettwoNodes  
    //           << '\t' << placementMetadata[b->name].longestBl  
    //           << '\t' << placementMetadata[b->name].longestTip
    //           << '\t' << placementMetadata[b->name].twoNodesitConnects.first
    //           << '\t' << placementMetadata[b->name].twoNodesitConnects.second
    //           << std::endl;

    Node *primary=a, *secondary=b, *tertiary=nullptr;
    int count = 0;
    while(placementOrder.size()!=0){
        Node* n=placementOrder.front(); placementOrder.pop();
        // std::cout << "Handling " << n->name << "\t" << primary->name << '\t' << (secondary?secondary->name:"NULL") << std::endl;
        if (n->name == t->root->name){
            double dist = (primary?primary->bl:0)+(secondary?secondary->bl:0);
            primary->bl=dist/2;
            secondary->bl=dist/2;
            continue;
        }

        if (n->children.size()!=0) { //internal node
            // Add internal node to placementMetadata with deepest node id and length
            // if internal node's parent is also done, move forward, otherwise, this internal node becomes right
            checkChildrenStatus(n, placementMetadata);
            if (specialCase){
                tertiary = n;
            }
            else if (primary->parent->name==n->name){
                specialCase=true;
                tertiary=n;
            } else {
                specialCase=false;
                secondary=n;
            }

        } else {
            placementMetadata[n->name]=placementMetadata_(2.0,0,n->idx);
            // primary dist
            int priTip = placementMetadata[primary->name].longestTip;
            double priDist;
            if (priTip>n->idx) priDist=h_mashDist[priTip*numSequences+n->idx];
            else               priDist=h_mashDist[n->idx*numSequences+priTip];
            priDist -= placementMetadata[primary->name].longestBl;

            //primary dist
            double secDist=0.;
            int secTip = placementMetadata[secondary->name].longestTip;
            if (secTip>n->idx) secDist=h_mashDist[secTip*numSequences+n->idx];
            else               secDist=h_mashDist[n->idx*numSequences+secTip];
            secDist -= placementMetadata[secondary->name].longestBl;
            
            double dist = (priDist+secDist-(primary->bl+secondary->bl))/2;
            if (dist<0)dist=0;

            // std::cout << "info \t" 
            //         << priDist << '\t'
            //         << secDist << '\t'
            //         << primary->bl << '\t'
            //         << secondary->bl << '\t'
            //         << dist << '\t'
            //         << std::endl;
            

            priDist-=dist; secDist-=dist;
            if (priDist<0)priDist=0;
            if (secDist<0)secDist=0;
            if (priDist>primary->bl+secondary->bl) dist+=priDist-(primary->bl+secondary->bl); priDist=primary->bl+secondary->bl;
            if (secDist>primary->bl+secondary->bl) dist+=secDist-(primary->bl+secondary->bl); secDist=primary->bl+secondary->bl;
            double rest=primary->bl+secondary->bl-priDist-secDist;
            priDist+=rest/2; secDist+=rest/2;
            
            primary->bl=priDist;
            if (secondary->bl>secDist) secondary->parent->bl = secondary->bl-(secDist);
            secondary->bl=secDist;
            // primary->bl-=dist;
            if (specialCase){
                tertiary->bl=dist;
                primary=tertiary;
            }
            else
                n->bl=dist;
           
            secondary=n;
            specialCase = false;
        }

        // std::cout << n->name << " (" << n->bl << ")"
        //       << '\t' << placementMetadata[n->name].distBettwoNodes  
        //       << '\t' << placementMetadata[n->name].longestBl  
        //       << '\t' << placementMetadata[n->name].longestTip
        //       << '\t' << primary->name << "(" << primary->bl << ")"
        //       << '\t' << (secondary?secondary->name:"NULL") << "(" << (secondary?secondary->bl:-1) << ")"
        //       << std::endl;
        // count++;
        // if (count==13) break;
    }

    // Node* lcs_ = lcs(a,b);
    // placementMetadata[lcs_->name]=placementMetadata_(dist, -1, -1, std::make_pair(a->name, b->name));
    
    // Node* lastNode = lcs_;
    // if (updateLcsNode(lcs_, placementMetadata)){
    //     lastNode = lcs_->parent;
    //     placementMetadata[lastNode->name]=placementMetadata_(-1, placementMetadata[lcs_->name].longestBl, placementMetadata[lcs_->name].longestTip, std::make_pair(lcs_->name, "NULL"));
    // }
    
    // std::cout << lcs_->name << '\t' << placementMetadata[lcs_->name].distBettwoNodes  
    //           << '\t' << placementMetadata[lcs_->name].longestBl  
    //           << '\t' << placementMetadata[lcs_->name].longestTip
    //           << '\t' << placementMetadata[lcs_->name].twoNodesitConnects.first
    //           << '\t' << placementMetadata[lcs_->name].twoNodesitConnects.second
    //           << std::endl;
    
    // int count=0;
    // while(placementOrderLeft.size()>0){
    //     Node* n=placementOrderLeft.front(); placementOrderLeft.pop();
    //     std::cout << "Placing " << n->name << std::endl;

    //     /* find parent node for left deepest node */
    //     Node* leftNode = lastNode;
    //     int leftFarthestNode=placementMetadata[leftNode->name].longestTip;
    //     double leftDist = (leftFarthestNode > n->idx) ? h_mashDist[leftFarthestNode*numSequences+n->idx]:h_mashDist[n->idx*numSequences+leftFarthestNode];
    //     leftDist -= placementMetadata[leftNode->name].longestBl; 
        
    //     /* find child node for right deepest node */
    //     Node* rightNode;
    //     for (auto& child: n->parent->children){
    //         if (child->name != n->name){
    //             rightNode=child;
    //             break;
    //         }
    //     }

    //     if (placementMetadata.find(rightNode->name) == placementMetadata.end()){
    //         placementMetadata[lastNode->name].distBettwoNodes = leftDist;
    //         placementMetadata[lastNode->name].twoNodesitConnects = std::make_pair(placementMetadata[lastNode->name].twoNodesitConnects.first,n->name);
    //     }

    //     std::cout << lastNode->name << '\t' << placementMetadata[lastNode->name].distBettwoNodes  
    //           << '\t' << placementMetadata[lastNode->name].longestBl  
    //           << '\t' << placementMetadata[lastNode->name].longestTip
    //           << '\t' << placementMetadata[lastNode->name].twoNodesitConnects.first
    //           << '\t' << placementMetadata[lastNode->name].twoNodesitConnects.second
    //           << std::endl;
        
    //     if (placementMetadata[lastNode->name].distBettwoNodes == 2)
    //         lastNode=lastNode->parent;

    //     count++;
    //     if (count==2) break;

    // }

    
}
