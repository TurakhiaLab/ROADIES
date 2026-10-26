/**
 * mash_placement.cuh
 *
 * Declarations for DIPPER MASH/placement pipeline: Param, MashDeviceArrays,
 * MSADeviceArrays, MatrixReader, PlacementDeviceArrays (exact), KPlacementDeviceArrays
 * (k-closest), NJDeviceArrays, and divide-and-conquer variants (DC structs).
 */

#ifndef MASHPL_CUH
#define MASHPL_CUH

#include <stdint.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>
#include <cstdio>
#include <string>
#include "tree.hpp"

namespace MashPlacement
{
    /** Default false: print multifurcating Newick by collapsing ~zero-length internal branches when present.
     *  Set true (e.g. --print-binary-newick) to print the fully resolved binary tree with explicit :0 edges. */
    inline bool g_printBinaryNewick = false;

    struct Param
    {
        uint64_t kmerSize;
        uint64_t sketchSize;
        uint64_t threshold;
        uint64_t distanceType;
        std::string in;
        std::string out;
        uint64_t batchSize;
        uint64_t backboneSize;
        uint64_t totalNumSeqs;
        std::pair<int, int> range;
        bool isProtein=false;
        bool useReducedProtein=false;
        bool useMurphy8=false;
        uint64_t clusterSize=0;

        Param(uint64_t t_kmerSize, uint64_t t_sketchSize, uint64_t t_threshold, uint64_t t_distanceType, std::string t_in, std::string t_out)
        {
            kmerSize = t_kmerSize; sketchSize = t_sketchSize; threshold = t_threshold, distanceType=t_distanceType;
            in = t_in, out = t_out; 
        };
    };

    /** Device arrays for MASH: compressed seqs, prefix/agg lengths, hash list, sketch construction. */
    struct MashDeviceArrays{
        uint64_t * d_compressedSeqs;
        uint64_t * d_prefixCompressed;
        uint64_t * d_aggseqLengths;
        uint64_t * d_seqLengths;
        size_t numSequences;
        uint64_t * d_hashList;

        uint64_t * h_hashList;

        void allocateDeviceArrays (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params);
        void printSketchValues(int numValues);
        void sketchConstructionOnGpu (Param& params);
        void deallocateDeviceArrays ();
        void distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const;

    };
    static MashDeviceArrays mashDeviceArrays;

    /** MASH device arrays for divide-and-conquer (backbone vs full, leaf map, etc.). */
    struct MashDeviceArraysDC{
        uint64_t * d_compressedSeqs;
        uint64_t * d_prefixCompressed;
        uint64_t * d_aggseqLengths;
        uint64_t * d_seqLengths;
        uint64_t * d_hashList;
        uint64_t * d_hashListConst;
        uint64_t * d_hashListBackbone;
        
        uint64_t * h_hashList;

        size_t      totalNumSequences;
        size_t      backboneSize;

        int * d_leafID_map;
        
        void allocateDeviceArraysDC (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params);
        void printSketchValuesDC(int numValues);
        void sketchConstructionOnGpuDC (Param& params, uint64_t** twoBitCompressedSeqs, uint64_t * h_seqLengths, uint64_t numSequences);
        void deallocateDeviceArraysDC ();
        void distConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist) const;
        void distConstructionOnCpu(Param& params, int rowId, double* d_mashDist) const;
        void distConstructionOnGpuForBackboneDC(Param& params, int rowId, double* d_mashDist) const;
        void distRangeConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist, int l, int r, bool clustering = false) const;
        void distRangeConstructionOnCpu(Param& params, int rowId, double* d_mashDist, int l, int r) const;
        void distSpecialIDConstructionOnGpuDC(Param& params, 
                                            int rowId, 
                                            double* d_mashDist, 
                                            int numToConstruct, 
                                            int * d_id,
                                            int * d_leafMap) const;
        void distSpecialIDConstructionOnCpu(Param& params, int rowId, std::vector<double> & h_mashDist, int numToConstruct, std::vector<int> & h_id) const;

    };
    static MashDeviceArraysDC mashDeviceArraysDC;

    /** Device arrays for aligned sequences: 4-bit compressed MSA, seq lengths, distance construction. */
    struct MSADeviceArrays{
        uint64_t * d_compressedSeqs;
        uint64_t * d_seqLengths;
        size_t numSequences;
        int seqLen;

        void allocateDeviceArrays (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params);
        void deallocateDeviceArrays ();
        void distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const;

    };
    static MSADeviceArrays msaDeviceArrays;

    /** MSA device arrays for DC (backbone vs const, seq len). */
    struct MSADeviceArraysDC{
        uint64_t * d_compressedSeqsBackBone;
        uint64_t * d_compressedSeqsConst;
        // uint64_t * d_seqLengths;
        // size_t numSequences;
        int d_seqLen;

        uint64_t * h_compressedSeqs=nullptr;

        size_t      totalNumSequences;
        size_t      backboneSize;

        void allocateDeviceArraysDC (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params, int gpuNum=0);
        void transferToDeviceArraysDC (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, int gpuCluster, Param& params);
        void distConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist) const;
        void distConstructionOnGpuForBackboneDC(Param& params, int rowId, double* d_mashDist) const;
        void distRangeConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist, int l, int r, bool clustering = false) const;
        void distSpecialIDConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist, int numToConstruct, int* d_id, int * d_leafMap) const;
        void deallocateDeviceArraysDC ();

    };
    static MSADeviceArraysDC msaDeviceArraysDC;

    /** PHYLIP distance matrix reader and device copy. */
    struct MatrixReader{
        int numSequences;
        double * h_dist;
        std::vector <std::string> name;
        FILE* filePtr;
        char * buffer;
        void allocateDeviceArrays (int num, FILE* fPtr);
        void distConstructionOnGpu(Param& params, int rowId, double* d_dist);
    };
    static MatrixReader matrixReader;

    /** Exact placement: adjacency, BFS/DFS/depth/levelst/leveled, lim; no k-closest. */
    struct PlacementDeviceArrays{
        int idx, bd;
        int numSequences;
        int * d_head;
        int * d_e;
        int * d_nxt;
        int * d_belong;
        int * d_bfsorder;
        int * d_dfsorder;
        int * d_dep;
        int * d_dfsrk;
        int * d_levelst;
        int * d_leveled;
        int * d_rev;
        double * d_dist;
        double * d_len;
        double * d_lim;

        void allocateDeviceArrays (size_t num);
        void deallocateDeviceArrays ();
        void findPlacementTree(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );
        void printTree(std::vector <std::string> name, std::ofstream& output_);
    };
    static PlacementDeviceArrays placementDeviceArrays;

    /** K-closest placement: closest_id/closest_dis per edge; build tree or add query to backbone. */
    struct KPlacementDeviceArrays{
        int idx, bd;
        int numSequences;
        int backboneSize;
        int * d_head;
        int * d_e;
        int * d_nxt;
        int * d_belong;
        int * d_closest_id;
        double * d_dist;
        double * d_len;
        double * d_closest_dis;

        void allocateDeviceArrays (size_t num, int backboneSize=-1);
        void deallocateDeviceArrays ();
        void initializeDeviceArrays(Tree* t);
        void findPlacementTree(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );

        void addQuery(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );

        // ---- Iterative (one-at-a-time) placement API -------------------------
        /** Temporary device arrays and cursor used by the iterative placement API.
         *  Populated by beginIterativePlacement(); freed by endIterativePlacement(). */
        int*    d_id_iter   = nullptr;
        int*    d_from_iter = nullptr;
        double* d_dis_iter  = nullptr;
        /** Index of the next query sequence to place (backboneSize after begin). */
        int     m_nextQueryIdx = -1;

        /** Allocate iterative-placement temp buffers and set the internal cursor
         *  to the first query sequence (index == backboneSize).
         *  Use when a backbone tree has already been loaded via initializeDeviceArrays(). */
        void beginIterativePlacement();

        /** Build the initial two-sequence tree from scratch (no backbone tree),
         *  allocate temp buffers, and set the cursor to sequence index 2.
         *  Use when no backbone tree is provided — the first two sequences in
         *  the dataset seed the tree, and all remaining sequences are queries. */
        void beginIterativePlacementFromScratch(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );

        /** Describes the edge split performed when one query sequence is placed.
         *
         *  Before placement, edge (splitNodeA ↔ splitNodeB) exists in the tree.
         *  Placement inserts a new internal node and a new tip node, producing:
         *
         *    splitNodeA ──[lenA]── internal ──[lenB]── splitNodeB
         *                              └──[lenTip]── tip (query)
         *
         *  The internal and tip DIPPER node indices are also provided so that
         *  callers can build or update their own node-pointer representation. */
        struct PlacementInfo {
            int    splitNodeA;      ///< DIPPER node index: "from" end of the split edge
            int    splitNodeB;      ///< DIPPER node index: "to"   end of the split edge
            int    internalNodeIdx; ///< DIPPER index of the newly inserted internal node
            int    tipNodeIdx;      ///< DIPPER index of the newly inserted tip (= seqIdx)
            double lenA;            ///< branch length: splitNodeA  → internal
            double lenB;            ///< branch length: splitNodeB  → internal
            double lenTip;          ///< branch length: internal    → tip
        };

        // ---- Two-phase placement (find / commit) --------------------------------
        /** Pending state saved between findBestEdge() and commitQuery(). */
        int    m_pendingSeqIdx  = -1;
        int    m_pendingEid     = -1;
        double m_pendingFracLen = 0.0;
        double m_pendingAddLen  = 0.0;
        /** Host-side per-edge fracLen/addLen computed by calculateBranchLength.
         *  Indexed by edge id; only entries with belong[eid] >= e[eid] are valid.
         *  Used by commitQuery(overrideEid) to look up lengths for an arbitrary edge. */
        std::vector<double> m_hostFracLen;
        std::vector<double> m_hostAddLen;

        /** Phase 1: compute distances and run calculateBranchLength for @p seqIdx,
         *  find the globally optimal edge, save its id/lengths in pending state,
         *  and return a PlacementInfo describing the would-be split.
         *  Does NOT call updateTreeStructure — call commitQuery() to finalise. */
        PlacementInfo findBestEdge(
            int seqIdx,
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );

        /** Return a PlacementInfo describing the split if @p seqIdx were placed
         *  on the given edge @p eid.  Reads edge geometry from GPU and looks up
         *  pre-computed branch lengths from the last findBestEdge() call.
         *  Must be called after findBestEdge() for the same seqIdx. */
        PlacementInfo getEdgeInfo(int eid, int seqIdx) const;

        /** Return the best live edges from the most recent findBestEdge(),
         * ranked by DIPPER's additional-distance objective. */
        std::vector<std::pair<PlacementInfo, double>> getTopKEdgeInfo(
            int seqIdx, int k) const;

        /** Search the current GPU adjacency list for an edge between @p nodeA and
         *  @p nodeB.  Returns the canonical edge id (belong[eid] >= e[eid]) used
         *  by calculateBranchLength and updateTreeStructure, or -1 if not found. */
        int findEdgeBetween(int nodeA, int nodeB) const;

        /** Phase 2: apply the pending placement (updateTreeStructure then
         *  updateClosestNodes) for @p seqIdx.
         *  @param overrideEid  >= 0: place on this edge instead of the DIPPER-chosen
         *                      one.  The edge must be present in the current tree and
         *                      findBestEdge() must have been called for the same seqIdx.
         *                      Pass -1 (default) to use the edge DIPPER selected. */
        void commitQuery(int seqIdx, int overrideEid = -1);
        // -------------------------------------------------------------------------

        /** Place the single query sequence at position @p seqIdx into the tree.
         *  Convenience wrapper: calls findBestEdge() then commitQuery() in one step.
         *
         *  @param out  If non-null, filled with the edge-split details.
         *  @return true if there are still more query sequences to place. */
        bool placeSingleQuery(
            int seqIdx,
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays,
            PlacementInfo* out = nullptr
        );

        /** Free iterative-placement temp buffers and reset the cursor. */
        void endIterativePlacement();
        // ----------------------------------------------------------------------

        void updateBranchLengthsFromTopology(
            Tree* t,
            double* h_mashDist,
            size_t numSeqs,
            uint64_t seed,
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );
        void printTree(std::vector <std::string> name, std::ofstream& output_);
        void printTree(std::vector <std::string> name, std::ofstream& output_, UnrootedTree* t, std::vector<int>& edgeIdsMappingToNewTreeEdges);

        void estimateBranchLengthsFromTopology(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays,
            std::vector<int>& edgeIdsMappingToNewTreeEdges,
            UnrootedTree* t,
            std::vector<std::string>& names
        );
    };
    static KPlacementDeviceArrays kplacementDeviceArrays;

    /** Conventional neighbor-joining: distance matrix and U matrix. */
    struct NJDeviceArrays
    {
        int d_numSequences;
        double * d_mashDist;
        double * d_U;
        void deallocateDeviceArrays ();
        void getDismatrix(
            int numSequences,
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );
        void findNeighbourJoiningTree(std::vector <std::string> &name, std::ofstream& output_);
    };
    static NJDeviceArrays njDeviceArrays;

    /** Host-side DC arrays for cluster trees and k-closest. */
    struct KPlacementDeviceArraysHostDC{
        int idx, bd;
        int numSequences;
        int totalNumSequences;
        int * h_head;
        int * h_e;
        int * h_nxt;
        int * h_belong;
        int * h_closest_id;
        int * h_closest_id_cluster;
        double * h_dist;
        double * h_len;
        double * h_closest_dis;
        double * h_closest_dis_cluster;
        int * clusterID;
        int requestedClusterTopK = 1;
        std::vector<std::vector<int>> clusterTopEdgeIds;
        std::vector<std::vector<double>> clusterTopEdgeAdditionalDistances;

        void allocateHostArraysDC (size_t num, size_t totalNum);
        void deallocateHostArraysDC ();
        void findClusterTreeDC(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArraysDC,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArraysDC
        );
        void printTreeCpuDC(std::vector <std::string> name);
    };
    static KPlacementDeviceArraysHostDC kplacementDeviceArraysHostDC;

    /** Device-side DC k-closest placement: backbone, clusters, per-cluster trees. */
    struct KPlacementDeviceArraysDC{
        int idx, bd;
        int numSequences;
        int totalNumSequences;
        int * d_head;
        int * d_e;
        int * d_nxt;
        int * d_belong;
        int * d_closest_id;
        int * d_closest_id_cluster;
        double * d_dist;
        double * d_len;
        double * d_closest_dis;
        double * d_closest_dis_cluster;

        void allocateDeviceArraysDC (size_t num, size_t totalNum, int gpuNum=0);
        void deallocateDeviceArraysDC ();
        void initializeDeviceArraysFromBackboneTreeDC(Tree* t);
        
        void findBackboneTreeDC(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays,
            const KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC,
            int gpuNum=0
        );

        void findClustersDC(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays,
            KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC
        );
        
        void findClustersDC_batch(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays,
            KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC,
            const int clusteringBatchIdx
        );

        void findClusterTreeDC(
            Param& params,
            MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            MSADeviceArraysDC& msaDeviceArrays,
            KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC,
            std::vector<int>& largeClustersIdx
        );
        void findClusterTreeDC_batch(
            Param& params,
            MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            MSADeviceArraysDC& msaDeviceArrays,
            KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC,
            const std::string dir,
            std::vector<bool>& isCluster
        );

        void findBackboneTreeDCRecursive(
            Param& params,
            MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            MSADeviceArraysDC& msaDeviceArrays,
            const KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC,
            int clusterIdx
        );
        

        void printTreeDC(
            std::vector <std::string> name,
            std::ofstream& output_,
            int rootNode=-1);
    };
    static KPlacementDeviceArraysDC kplacementDeviceArraysDC;

};



#endif
