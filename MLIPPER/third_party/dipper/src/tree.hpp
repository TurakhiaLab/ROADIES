/**
 * tree.hpp
 *
 * In-memory phylogenetic tree: Node (name, branch length, parent, children, level)
 * and Tree (root, allNodes, Newick parse/serialization, DFS, leaf/internal counts).
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>


class Node 
{
public:
    Node(const std::string& id, double len);       /** Root or standalone node. */
    Node(const std::string& id, Node* par, double len);  /** Child of par. */
    size_t getNumLeaves();
    size_t getNumNodes();
    bool is_leaf() {return !(name.substr(0,4) == "node");}
    bool is_leaf_new() {return children.size() == 0;}
    void setNumleaves() {numLeaves = getNumLeaves();};

    std::string name;
    int idx;
    double bl;
    size_t level;
    int lca=0;

    Node* parent;
    std::vector< Node* > children;
    size_t numLeaves = {0};
    float weight = {0};

    std::vector<std::pair<int, double>> closestIds;
};

class Tree
{
    private:
        bool m_isRooted = true;

    public:

        size_t m_currInternalNode{ 0 };
        size_t m_maxDepth{ 0 };
        size_t m_numLeaves{ 0 };
        size_t m_numLeafID{ 0 };
        float m_meanDepth{ 0 };
        std::string newInternalNodeId() { return "node_" + std::to_string(++m_currInternalNode);}

        Node* root;
        std::unordered_map< std::string, Node* > allNodes;
        /** If randomResolvePolytomies, multifurcating nodes are resolved to binary by inserting
         *  internal nodes with branch length 0 to the parent. polytomyRngSeed==0 uses std::random_device. */
        Tree(std::string newick, size_t totalLeaves = 0, bool randomResolvePolytomies = true,
             uint64_t polytomyRngSeed = 0);
        Tree(Node* node);
        Tree() {root = nullptr;}
        ~Tree();
        static Tree* fromUnrootedNewick(const std::string& newick, size_t totalLeaves = 0);
        void reRootTree();
        void setRootBetweenNodes(std::string nodeName1, std::string nodeName2);
        bool isRooted() const { return m_isRooted; }
        std::string getNewickString(Node* node = nullptr);
        std::vector<std::string> remove_and_rename(std::string node_name);
        void dfsExpansion(Node* node, std::vector< Node* >& vec);
        void updateClosestIds(std::unordered_map<std::string, std::vector<std::pair<int, double>>>& closestIds);

};



void placementAccuracy(Tree* t, double* h_mashDist, size_t numSequences);

// =============================================================================
// Unrooted Tree Interface
// =============================================================================

struct UnrootedNode;

/** Neighbor with branch length and edge ID */
struct UnrootedNeighbor {
    UnrootedNode* node;
    double length;
    int edge_id;
};

/** Node for unrooted trees: adjacency-based (neighbors) instead of parent/children */
struct UnrootedNode {
    std::string name;
    int idx;
    double bl;  // branch length to parent (meaningful only after rooting)
    bool is_leaf() const { return neighbors.size() == 1; }
    size_t degree() const { return neighbors.size(); }

    std::vector<UnrootedNeighbor> neighbors;
};

/** Edge in an unrooted tree (from, to, branch length, edge ID) */
struct UnrootedEdge {
    UnrootedNode* from;
    UnrootedNode* to;
    double length;
    int edge_id;
};

/** Abstract interface for unrooted phylogenetic trees */
class IUnrootedTree {
public:
    virtual ~IUnrootedTree() = default;

    virtual size_t numNodes() const = 0;
    virtual size_t numLeaves() const = 0;
    virtual size_t numEdges() const = 0;

    virtual const std::unordered_map<std::string, UnrootedNode*>& getNodes() const = 0;
    virtual UnrootedNode* getNode(const std::string& name) const = 0;

    /** Get neighbors of a node with branch lengths and edge IDs */
    virtual const std::vector<UnrootedNeighbor>& getNeighbors(UnrootedNode* node) const = 0;

    /** Get edge ID for the edge between from and to (-1 if not adjacent) */
    virtual int getEdgeId(UnrootedNode* from, UnrootedNode* to) const = 0;

    /** Get edge ID for the i-th neighbor of node (-1 if index out of range) */
    virtual int getEdgeId(UnrootedNode* node, size_t neighborIndex) const = 0;

    /** Root the tree on the edge between node1 and node2 (modifies in place) */
    virtual void rootBetween(const std::string& nodeName1, const std::string& nodeName2) = 0;

    /** Root the tree at an internal node (creates new root on an edge, modifies in place) */
    virtual void rootAt(const std::string& nodeName) = 0;

    /** Export to Newick string (unrooted) */
    virtual std::string toNewick() const = 0;

    /** BFS from start node, returns nodes in breadth-first order */
    virtual std::vector<UnrootedNode*> bfs(UnrootedNode* start) const = 0;

    /** DFS from start node, returns nodes in depth-first order */
    virtual std::vector<UnrootedNode*> dfs(UnrootedNode* start) const = 0;

    /** Edges along the unique path between two nodes (empty if nodes not found or disconnected) */
    virtual std::vector<UnrootedEdge> edgesBetween(const std::string& nodeName1,
                                                  const std::string& nodeName2) const = 0;

    /** Nodes along the unique path between two nodes (empty if nodes not found or disconnected) */
    virtual std::vector<UnrootedNode*> nodesBetween(const std::string& nodeName1,
                                                  const std::string& nodeName2) const = 0;

    /** Create from Newick string */
    static IUnrootedTree* fromNewick(const std::string& newick, size_t totalLeaves = 0);
};

/** Concrete implementation of unrooted tree using adjacency list */
class UnrootedTree : public IUnrootedTree {
public:
    UnrootedTree() = default;
    explicit UnrootedTree(const std::string& newick, size_t totalLeaves = 0);
    /** Build tree from newick, aligned with reference (same node names, neighbor order).
     *  Used when comparing trees that differ only by one tip placement.
     *  tipName: the tip that was moved between the two trees. */
    UnrootedTree(const std::string& newick, const UnrootedTree* reference,
                 const std::string& tipName, size_t totalLeaves = 0);
    ~UnrootedTree() override;

    size_t numNodes() const override { return m_nodes.size(); }
    size_t numLeaves() const override { return m_numLeaves; }
    size_t numEdges() const override { return m_numEdges; }

    const std::unordered_map<std::string, UnrootedNode*>& getNodes() const override { return m_nodes; }
    UnrootedNode* getNode(const std::string& name) const override;

    const std::vector<UnrootedNeighbor>& getNeighbors(UnrootedNode* node) const override;
    int getEdgeId(UnrootedNode* from, UnrootedNode* to) const override;
    int getEdgeId(UnrootedNode* node, size_t neighborIndex) const override;

    void rootBetween(const std::string& nodeName1, const std::string& nodeName2) override;
    void rootAt(const std::string& nodeName) override;

    std::string toNewick() const override;

    std::vector<UnrootedNode*> bfs(UnrootedNode* start) const override;
    std::vector<UnrootedNode*> dfs(UnrootedNode* start) const override;
    std::vector<UnrootedEdge> edgesBetween(const std::string& nodeName1,
                                          const std::string& nodeName2) const override;
    std::vector<UnrootedNode*> nodesBetween(const std::string& nodeName1,
                                        const std::string& nodeName2) const override;
    std::string getRoot() const { return m_root->name; }

    /** Collapse internal edges whose non-trivial split is not listed in `other` (polytomy vs binary). */
    void collapseToSplitsDisplayedIn(const UnrootedTree& other);
    /** Same, using precomputed keys from the partner tree (before mutating either tree). */
    void collapseToSplitsDisplayedIn(const std::set<std::string>& allowedNonTrivialSplitKeys);

    /** Contract edges with length <= lengthEps between two non-leaf nodes (implicit polytomies from
     *  binary Newick with zero internal branch lengths). */
    void collapseZeroLengthInternalEdges(double lengthEps = 1e-12);

    UnrootedNode* m_root{nullptr};
private:
    std::unordered_map<std::string, UnrootedNode*> m_nodes;
    size_t m_numLeaves{0};
    size_t m_numEdges{0};
    size_t m_currInternalNode{0};

    std::string newInternalNodeId() { return "node_" + std::to_string(++m_currInternalNode); }
    void buildFromNewick(const std::string& newick, size_t totalLeaves);
    void reindexEdges();
    void contractUndirectedEdge(UnrootedNode* u, UnrootedNode* v);
};

std::set<std::string> collectNonTrivialSplitKeys(const UnrootedTree& tree);

/** Compute edge distance between two unrooted trees that differ by one tip placement.
 *  Both trees have the same tips; one tip was moved from its position in t1 to a new
 *  position in t2. Returns the number of edges between the two attachment points,
 *  or -1 on error. */
int computeTipPlacementEdgeDistance(const UnrootedTree* t1, const UnrootedTree* t2,
                                    const std::string& tipName);