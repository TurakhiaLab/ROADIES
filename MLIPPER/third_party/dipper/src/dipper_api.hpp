/**
 * dipper_api.hpp
 *
 * Public C++ library API for DIPPER — phylogenetic placement of aligned MSA
 * sequences onto a backbone tree using GPU-accelerated k-closest placement.
 *
 * This header is the single include required to embed DIPPER as a library.
 * Link your target against the same object files that build the `dipper`
 * executable (see CMakeLists.txt) and add src/ to your include path.
 *
 * Two workflows are supported:
 *
 *   WITH backbone tree                 WITHOUT backbone tree
 *   ─────────────────────────────────  ─────────────────────────────────────
 *   1. loadBackboneTree()              (skip — no backbone needed)
 *   2. loadSequences()                 1. loadSequences()
 *   3. initializeGPU()                 2. initializeGPU()
 *   4. placeNextSequence() × N         3. placeNextSequence() × (N-2)
 *   5. writeTree()                     4. writeTree()
 *
 *   When no backbone tree is given, the first two sequences in the input seed
 *   an initial two-node tree; all remaining sequences are placed one at a time.
 *
 * Example
 * -------
 *   See examples/msa_placement_example.cpp.
 */

#pragma once

#include "fourBitCompressor.hpp"
#include "twoBitCompressor.hpp"
#include "mash_placement.cuh"
#include "tree.hpp"

#include <cuda_runtime.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <unistd.h>
#include <algorithm>
#include <functional>
#include <iterator>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

void calculateBranchLengthSpecialIDCpu(
    int num,
    int* head,
    int* nxt,
    std::vector<double>& dis,
    int* e,
    double* len,
    int* belong,
    std::tuple<int, double, double>* minPos,
    int lim,
    double* closest_dis,
    int* closest_id,
    int numToCalculate,
    std::vector<int>& h_edgeMask);

void updateClosestNodesCpu(
    int* head,
    int* nxt,
    int* e,
    double* len,
    double* closest_dis,
    int* closest_id,
    int x,
    int* id,
    int* from,
    double* dis,
    int* h_edgeMaskIndex);

void updateClosestNodesInClusterCpu(
    int* head,
    int* nxt,
    int* e,
    double* len,
    double* closest_dis,
    int* closest_id,
    int x,
    std::vector<int>& id,
    std::vector<int>& from,
    std::vector<double>& dis,
    int cluster_eid,
    int* belong);

void updateTreeStructureInClusterCpu(
    int* head,
    int* nxt,
    int* e,
    double* len,
    double* closest_dis,
    int* closest_id,
    int* belong,
    int eid,
    double fracLen,
    double addLen,
    int placeId,
    int edgeCount,
    int totalNumSequences,
    int placeCount);

void updateClusterInfoCpu(
    int leafID,
    int edgeidx,
    std::vector<int>& h_leafMask,
    std::vector<int>& h_edgeMask,
    int edgeCount,
    int leafCount);

void initializeClusterCpu(
    int eid,
    int* e,
    int* belong,
    int* head,
    int* nxt,
    int* closest_id,
    std::vector<int>& edgeMask,
    std::vector<int>& leafMask);

/**
 * PlacementResult
 *
 * Returned by DipperSession::placeNextSequence().
 *
 * Describes the edge split that occurred when one query sequence was placed
 * onto the tree.  Before placement, edge (splitNodeA ↔ splitNodeB) existed.
 * After placement, that edge is replaced by three new edges:
 *
 *   splitNodeA ──[lenA]── internalNodeIdx ──[lenB]── splitNodeB
 *                                 └──[lenTip]── tipNodeIdx  (query)
 *
 * The node indices are DIPPER's internal numbering.  Map them to sequence
 * names via DipperSession::sortedNames()[idx].
 *
 * For integration with an external node-pointer tree (N nodes before this
 * placement), create:
 *   • external node N+1 corresponding to internalNodeIdx
 *   • external node N+2 corresponding to tipNodeIdx
 * and connect them according to (splitNodeA, splitNodeB, lenA, lenB, lenTip).
 */
struct PlacementResult {
    bool more = false;      ///< true if more queries remain to be placed

    // ---- Split edge (DIPPER internal indices) --------------------------------
    int    splitNodeA = -1; ///< DIPPER idx: "from" end of the edge that was split
    int    splitNodeB = -1; ///< DIPPER idx: "to"   end of the edge that was split

    // ---- Split edge (sequential labels) -------------------------------------
    int    splitLabelA = -1; ///< sequential label of splitNodeA (use session.nodeLabel())
    int    splitLabelB = -1; ///< sequential label of splitNodeB

    // ---- New nodes (DIPPER internal indices) --------------------------------
    int    internalNodeIdx = -1; ///< DIPPER idx of the new internal node
    int    tipNodeIdx      = -1; ///< DIPPER idx of the new tip (query) node

    // ---- New nodes (sequential labels) --------------------------------------
    /// Sequential label of the new internal node.
    /// If N nodes existed before this placement, internalLabel = N.
    int    internalLabel = -1;
    /// Sequential label of the new tip (query) node.
    /// If N nodes existed before this placement, tipLabel = N+1.
    int    tipLabel = -1;

    std::string tipName;    ///< sequence name of the placed query

    // ---- Branch lengths (all >= 0) ------------------------------------------
    double lenA   = 0.0;    ///< splitNodeA  → internalNodeIdx
    double lenB   = 0.0;    ///< splitNodeB  → internalNodeIdx
    double lenTip = 0.0;    ///< internalNodeIdx → tipNodeIdx (query)
};

struct RankedPlacementCandidate {
    int splitLabelA = -1;
    int splitLabelB = -1;
    double additionalDistance = 0.0;
};

/** Name/sequence pair exported from the current DIPPER session inputs. */
struct NamedSequence {
    std::string name;
    std::string sequence;
};

/** Metadata describing the next query that DIPPER would place. */
struct PendingQueryInfo {
    int sequenceIdx = -1;   ///< global sorted sequence index in the current session
    int clusterId   = -1;   ///< divide-and-conquer cluster id (backbone mode only)
    std::string name;
    std::string sequence;
    int splitLabelA = -1;   ///< candidate-edge endpoint label A for D&C scheduling
    int splitLabelB = -1;   ///< candidate-edge endpoint label B for D&C scheduling
    double lenA = 0.0;      ///< endpoint A to the proposed insertion point
    double lenB = 0.0;      ///< endpoint B to the proposed insertion point
    double lenTip = 0.0;    ///< proposed query pendant length
};

struct DivideAndConquerAdvanceResult {
    int splitLabelA = -1;
    int splitLabelB = -1;
    int internalLabel = -1;
    int tipLabel = -1;
    std::string tipName;
};

/** One of DIPPER's ranked immutable-backbone routing candidates. */
struct RankedDipperCandidate {
    int clusterId = -1;
    int splitLabelA = -1;
    int splitLabelB = -1;
    double additionalDistance = std::numeric_limits<double>::infinity();
};

/** Original top-1 assignment plus the requested top-k backbone edges. */
struct QueryTopKBackboneCandidates {
    int sequenceIdx = -1;
    std::string name;
    int originalClusterId = -1;
    std::vector<RankedDipperCandidate> candidates;
};

/**
 * DipperSession
 *
 * Encapsulates one complete DIPPER MSA-placement session.  Each session owns
 * its GPU device arrays independently of the static singletons used by the
 * CLI, so multiple sessions can be constructed in the same process (serially).
 *
 * The class is non-copyable because it owns GPU resources.
 */
class DipperSession {
public:
    DipperSession()  = default;
    ~DipperSession() { _cleanup(); }

    DipperSession(const DipperSession&)            = delete;
    DipperSession& operator=(const DipperSession&) = delete;

    // =========================================================================
    // Step 1 (optional) — Backbone tree
    // =========================================================================

    /**
     * Parse @p newickStr as the backbone tree in Newick format.
     *
     * This step is OPTIONAL.  If omitted, the first two sequences supplied to
     * loadSequences() will seed the initial tree, and all remaining sequences
     * will be treated as queries.
     *
     * @param newickStr   A complete Newick string (single line, terminated by ';').
     *
     * Must be called before loadSequences() when a backbone tree is desired.
     */
    void loadBackboneTree(const std::string& newickStr)
    {
        _requireOneOf({State::Empty}, "loadBackboneTree");
        // Store the raw Newick string.  The Tree object is constructed (with
        // the correct totalLeaves count) later in loadSequences(), once the
        // full sequence set is known.  A temporary parse here is only used to
        // read m_backboneSize for the diagnostic message.
        m_backboneNewick = newickStr;
        delete m_tree;
        m_tree         = new Tree(newickStr);   // temporary; rebuilt in loadSequences
        m_backboneSize = m_tree->m_numLeaves;
        m_state        = State::TreeLoaded;
        std::cerr << "[DIPPER] Backbone tree stored: " << m_backboneSize
                  << " leaves, root = " << m_tree->root->name << "\n";
    }

    /**
     * Read the first line of @p filePath and parse it as a Newick backbone tree.
     *
     * @return true on success, false if the file cannot be opened.
     */
    bool loadBackboneTreeFromFile(const std::string& filePath)
    {
        _requireOneOf({State::Empty}, "loadBackboneTreeFromFile");
        std::ifstream f(filePath);
        if (!f) {
            std::cerr << "[DIPPER] ERROR: cannot open backbone tree file: " << filePath << "\n";
            return false;
        }
        std::string newick;
        std::getline(f, newick);
        loadBackboneTree(newick);
        return true;
    }

    // =========================================================================
    // Step 2 — Aligned sequences
    // =========================================================================

    /**
     * Supply all aligned sequences as pre-loaded in-memory vectors.
     *
     * WITH backbone tree
     * ------------------
     * Each sequence is matched against the backbone tree by name:
     *   • If the name exists in the backbone, the sequence fills its tree-index
     *     slot (used only for distance computation).
     *   • If the name is absent, the sequence is appended as a query and will
     *     be placed during the placeNextSequence() calls.
     *
     * WITHOUT backbone tree
     * ---------------------
     * All sequences are accepted as-is in the order provided.  The first two
     * sequences seed the initial two-node tree; all remaining sequences are
     * placed as queries.  At least 3 sequences are required.
     *
     * @param sequences  Aligned sequence strings (all the same length for MSA).
     * @param names      Sequence names corresponding to each entry in @p sequences.
     */
    void loadSequences(const std::vector<std::string>& sequences,
                       const std::vector<std::string>& names)
    {
        _requireOneOf({State::Empty, State::TreeLoaded}, "loadSequences");
        if (sequences.size() != names.size())
            throw std::invalid_argument(
                "[DIPPER] loadSequences: sequences and names must have the same size.");

        const bool hasBackbone = (m_tree != nullptr);
        if (hasBackbone && sequences.empty())
            throw std::invalid_argument(
                "[DIPPER] loadSequences: at least one sequence is required.");
        if (!hasBackbone && sequences.size() < 3)
            throw std::invalid_argument(
                "[DIPPER] loadSequences: at least 3 sequences are required "
                "when no backbone tree is provided (first two seed the initial tree).");

        m_rawSeqs      = sequences;
        m_rawNames     = names;
        m_numSequences = sequences.size();
        m_rawNameToIndex.clear();
        m_rawNameToIndex.reserve(m_rawNames.size());
        for (size_t idx = 0; idx < m_rawNames.size(); ++idx) {
            const auto [_, inserted] =
                m_rawNameToIndex.emplace(m_rawNames[idx], idx);
            if (!inserted) {
                throw std::invalid_argument(
                    "[DIPPER] loadSequences: duplicate sequence name '" +
                    m_rawNames[idx] + "'.");
            }
        }

        if (hasBackbone) {
            // Rebuild the tree with totalLeaves = m_numSequences so that
            // internal node indices are offset above (m_numSequences - 1) and
            // cannot collide with query sequence slots.  This matches the CLI
            // behaviour (tree_generation.cu: new Tree(newick, namesDump.size())).
            delete m_tree;
            m_tree = new Tree(m_backboneNewick, m_numSequences);
            m_backboneSize = m_tree->m_numLeaves;
            _buildIdMapWithBackbone();
            _assignPreorderLabels();   // label all backbone nodes 0..K-1
            std::cerr << "[DIPPER] Sequences loaded: " << m_numSequences
                      << " total (" << m_backboneSize << " backbone, "
                      << (m_numSequences - m_backboneSize) << " queries).\n" << std::flush;
        } else {
            _buildIdMapNoBackbone();  // seeds get labels 0 and 1
            std::cerr << "[DIPPER] Sequences loaded: " << m_numSequences
                      << " total (no backbone — first 2 will seed the tree, "
                      << (m_numSequences - 2) << " queries)." << std::flush;
        }
        m_state = State::SeqsLoaded;
    }

    // =========================================================================
    // Step 3 — GPU initialization
    // =========================================================================

    /**
     * Compress all sequences with 4-bit encoding, upload them to GPU memory,
     * and prepare the initial tree topology on the GPU.
     *
     * With a backbone tree: loads the provided tree topology onto the GPU.
     * Without a backbone tree: computes the distance between the first two
     *   sequences and builds the initial two-node tree on the GPU.
     *
     * The @p params object is updated in-place:
     *   • params.in  is forced to "m" (MSA / aligned input).
     *   • params.out is forced to "t" (tree output).
     *   • params.totalNumSeqs and params.backboneSize are set automatically.
     *
     * All other fields (range, isProtein, threshold, distanceType, etc.) are
     * respected as supplied by the caller.
     *
     * Must be called after loadSequences().
     */
    void initializeGPU(MashPlacement::Param& params)
    {
        _requireOneOf({State::SeqsLoaded}, "initializeGPU");
        _ensureNoActivePipeline("initializeGPU");
        m_params            = &params;
        params.in           = "m";
        params.out          = "t";
        params.totalNumSeqs = m_numSequences;
        params.backboneSize = m_backboneSize;   // 0 when no backbone

        std::cerr << "[DIPPER] Compressing sequences (4-bit)...\n";
        _compressSequences(params);

        std::cerr << "[DIPPER] Uploading sequences to GPU...\n";
        m_msaArrays.allocateDeviceArrays(
            m_fourBitSeqs, m_seqLengths, m_numSequences, params);

        const bool hasBackbone = (m_tree != nullptr);
        if (hasBackbone) {
            std::cerr << "[DIPPER] Loading backbone tree topology onto GPU...\n";
            m_kpArrays.allocateDeviceArrays(m_numSequences, (int)m_backboneSize);
            m_kpArrays.initializeDeviceArrays(m_tree);
            m_kpArrays.beginIterativePlacement();
            std::cerr << "[DIPPER] GPU initialized. Ready to place "
                      << (m_numSequences - m_backboneSize) << " query sequence(s).\n";
        } else {
            std::cerr << "[DIPPER] Building initial two-sequence tree on GPU...\n";
            // backboneSize == -1 signals scratch mode to allocateDeviceArrays.
            m_kpArrays.allocateDeviceArrays(m_numSequences, /*backboneSize=*/-1);
            m_kpArrays.beginIterativePlacementFromScratch(
                params, m_mashArrays, m_matrixReader, m_msaArrays);
            std::cerr << "[DIPPER] GPU initialized. Ready to place "
                      << (m_numSequences - 2) << " query sequence(s).\n";
        }

        m_activePipeline = ActivePipeline::Iterative;
        m_state = State::GPUReady;
    }

    /**
     * Initialize the legacy divide-and-conquer GPU pipeline.
     *
     * If a backbone tree was loaded, backbone/query partitioning follows the
     * loaded backbone match and @p backboneSize must be 0 or equal to that
     * backbone size. The loaded topology is reused as the fixed DC backbone;
     * only query clustering and cluster-local placement use the DC machinery.
     *
     * If no backbone tree was loaded, backboneSize == 0 uses the legacy CLI
     * heuristic max(2, numSequences / 20).
     */
    void initializeDivideAndConquerGPU(
        MashPlacement::Param& params,
        size_t backboneSize = 0)
    {
        _requireOneOf({State::SeqsLoaded}, "initializeDivideAndConquerGPU");
        _ensureNoActivePipeline("initializeDivideAndConquerGPU");

        const size_t resolved_backbone_size =
            _resolveDivideAndConquerBackboneSize(backboneSize);
        if (resolved_backbone_size < 2) {
            throw std::invalid_argument(
                "[DIPPER] initializeDivideAndConquerGPU: backboneSize must be at least 2.");
        }
        if (resolved_backbone_size > m_numSequences) {
            throw std::invalid_argument(
                "[DIPPER] initializeDivideAndConquerGPU: backboneSize exceeds the number of loaded sequences.");
        }

        m_params            = &params;
        const bool use_mash_routing =
            std::getenv("DIPPER_DC_USE_MASH_ROUTING") != nullptr;
        params.in           = use_mash_routing ? "r" : "m";
        params.out          = "t";
        params.totalNumSeqs = m_numSequences;
        params.backboneSize = resolved_backbone_size;
        params.batchSize    = std::max<uint64_t>(
            1,
            static_cast<uint64_t>(resolved_backbone_size));

        m_dcBackboneSize = resolved_backbone_size;
        m_dcBackboneBuilt = false;
        m_dcTreeBuilt = false;
        m_queryClusterIds.clear();
        m_queryEffectiveClusterIds.clear();
        m_dcScheduledQueryNames.clear();
        m_queryClusterIdsReady = false;
        m_computedQueryClusterTopK = 0;
        m_queryClusterBackboneSize = 0;
        _resetDivideAndConquerSchedulerState();
        if (m_hostDC.clusterID != nullptr) {
            delete[] m_hostDC.clusterID;
            m_hostDC.clusterID = nullptr;
        }

        if (use_mash_routing) {
            std::cerr << "[DIPPER] Compressing sequences (2-bit) and building Mash sketches for divide-and-conquer...\n";
            _compressSequencesTwoBit(params);
            m_mashArraysDC.allocateDeviceArraysDC(
                m_twoBitSeqs, m_seqLengths, m_numSequences, params);
            m_mashArraysDC.sketchConstructionOnGpuDC(
                params, m_twoBitSeqs, m_seqLengths, m_numSequences);
            // Mash owns routing only. Keep an aligned buffer resident for the
            // subsequent cluster-local placement phase in the same session.
            _compressSequences(params);
            m_msaArraysDC.allocateDeviceArraysDC(
                m_fourBitSeqs, m_seqLengths, m_numSequences, params);
        } else {
            std::cerr << "[DIPPER] Compressing sequences (4-bit) for divide-and-conquer...\n";
            _compressSequences(params);
            std::cerr << "[DIPPER] Uploading sequences to divide-and-conquer GPU buffers...\n";
            m_msaArraysDC.allocateDeviceArraysDC(
                m_fourBitSeqs, m_seqLengths, m_numSequences, params);
        }
        m_kpArraysDC.allocateDeviceArraysDC(
            resolved_backbone_size,
            m_numSequences);

        m_hostDC.numSequences = static_cast<int>(resolved_backbone_size);
        m_hostDC.totalNumSequences = static_cast<int>(m_numSequences);
        m_hostDC.bd = 2;
        m_hostDC.idx = 0;
        if (m_tree != nullptr) {
            m_kpArraysDC.initializeDeviceArraysFromBackboneTreeDC(m_tree);
            m_hostDC.bd = static_cast<int>(resolved_backbone_size);
            m_hostDC.idx = static_cast<int>(resolved_backbone_size) * 4 - 4;
            m_dcBackboneBuilt = true;
            m_dcRootNode = m_tree->root != nullptr ? m_tree->root->idx : -1;
            if (m_dcRootNode < 0) {
                throw std::runtime_error(
                    "[DIPPER] initializeDivideAndConquerGPU: fixed backbone root is unavailable.");
            }
            std::cerr << "[DIPPER] Loaded backbone topology reused for divide-and-conquer."
                      << " root=" << m_dcRootNode
                      << "\n";
        }

        std::cerr << "[DIPPER] Divide-and-conquer GPU initialized."
                  << " backbone=" << resolved_backbone_size
                  << " total=" << m_numSequences
                  << " queries=" << (m_numSequences - resolved_backbone_size)
                  << "\n";

        m_activePipeline = ActivePipeline::DivideAndConquer;
        m_state = State::GPUReady;
    }

    void prepareDivideAndConquerAlignedPlacement()
    {
        _ensureDivideAndConquerPipeline(
            "prepareDivideAndConquerAlignedPlacement");
        _ensureDivideAndConquerClustersReady();
        if (m_msaArraysDC.h_compressedSeqs == nullptr) {
            throw std::runtime_error(
                "[DIPPER] prepareDivideAndConquerAlignedPlacement: aligned sequence buffer is unavailable.");
        }
        m_params->in = "m";
        std::cerr
            << "[DIPPER] Divide-and-conquer routing frozen; switched to aligned cluster-local placement.\n";
    }

    /**
     * Rebuild the complete D&C pipeline from an externally updated tree.
     *
     * This is the correctness-first synchronization path for integrations
     * which modify topology outside DIPPER (for example, MLIPPER local SPR).
     * Reconstructing the session avoids relying on incremental directed-edge,
     * closest-leaf, and scheduler-cache repair after a topology rewrite.
     */
    void reloadDivideAndConquerGPU(
        const std::string& newick,
        const std::vector<std::pair<std::string, int>>& tipStableLabels,
        MashPlacement::Param& params)
    {
        _ensureDivideAndConquerPipeline("reloadDivideAndConquerGPU");
        const std::vector<std::string> sequences = m_rawSeqs;
        const std::vector<std::string> names = m_rawNames;

        _cleanup();
        m_state = State::Empty;
        m_backboneNewick.clear();
        m_backboneSize = 0;
        m_dcBackboneTipStableLabels.clear();

        loadBackboneTree(newick);
        setBackboneTipStableLabels(tipStableLabels);
        loadSequences(sequences, names);
        initializeDivideAndConquerGPU(params, 0);
    }

    // =========================================================================
    // Step 4 — Place one query sequence at a time
    // =========================================================================

    /**
     * Place the next query sequence onto the current tree and return the
     * edge-split details needed to update an external node-pointer tree.
     *
     * Each call places exactly one sequence.  Calls must be made in order
     * after initializeGPU().
     *
     * With backbone:    queries start from index backboneSize.
     * Without backbone: queries start from index 2 (first two seeded the tree).
     *
     * @return A PlacementResult whose fields describe:
     *   • more            — whether additional queries remain.
     *   • splitNodeA/B    — DIPPER indices of the edge endpoints that were split.
     *   • internalNodeIdx — DIPPER index of the newly inserted internal node.
     *   • tipNodeIdx      — DIPPER index of the newly inserted tip (query).
     *   • tipName         — sequence name of the placed query.
     *   • lenA/lenB/lenTip — branch lengths for the three new edges.
     *
     * @throws std::logic_error if all queries have already been placed.
     */
    PlacementResult placeNextSequence()
    {
        _requireOneOf({State::GPUReady}, "placeNextSequence");
        _ensureIterativePipeline("placeNextSequence");
        if (m_placementPending)
            throw std::logic_error(
                "[DIPPER] placeNextSequence: a two-phase placement is in progress. "
                "Call commitPlacement() first.");
        if (m_kpArrays.m_nextQueryIdx >= (int)m_numSequences)
            throw std::logic_error(
                "[DIPPER] placeNextSequence: all query sequences have already been placed.");

        const int seqIdx    = m_kpArrays.m_nextQueryIdx;
        const int seedCount = (m_tree != nullptr) ? (int)m_backboneSize : 2;
        const int queryNum  = seqIdx - seedCount + 1;
        const int totalQ    = (int)m_numSequences - seedCount;

        std::cerr << "[DIPPER] Placing query " << queryNum << " / " << totalQ
                  << " (global idx " << seqIdx
                  << ", name: " << m_sortedNames[seqIdx] << ")...\n";

        MashPlacement::KPlacementDeviceArrays::PlacementInfo info{};
        bool more = m_kpArrays.placeSingleQuery(
            seqIdx, *m_params, m_mashArrays, m_matrixReader, m_msaArrays, &info);

        // Assign sequential labels to the two new nodes and register them.
        const int iLabel = m_nextLabel;       // new internal node label (N)
        const int tLabel = m_nextLabel + 1;   // new tip node label      (N+1)
        m_dipperToLabel[info.internalNodeIdx] = iLabel;
        m_dipperToLabel[info.tipNodeIdx]      = tLabel;
        m_labelToDipper[iLabel]               = info.internalNodeIdx;
        m_labelToDipper[tLabel]               = info.tipNodeIdx;
        m_nextLabel += 2;

        PlacementResult result;
        result.more            = more;
        result.splitNodeA      = info.splitNodeA;
        result.splitNodeB      = info.splitNodeB;
        result.splitLabelA     = nodeLabel(info.splitNodeA);
        result.splitLabelB     = nodeLabel(info.splitNodeB);
        result.internalNodeIdx = info.internalNodeIdx;
        result.tipNodeIdx      = info.tipNodeIdx;
        result.internalLabel   = iLabel;
        result.tipLabel        = tLabel;
        result.tipName         = m_sortedNames[seqIdx];
        result.lenA            = info.lenA;
        result.lenB            = info.lenB;
        result.lenTip          = info.lenTip;

        if (!more) {
            m_state = State::Done;
            std::cerr << "[DIPPER] All query sequences placed.\n";
        }
        return result;
    }

    // =========================================================================
    // Step 4 — Two-phase placement: findNextPlacement / commitPlacement
    // =========================================================================

    /**
     * Phase 1 of two-phase placement.
     *
     * Computes DIPPER's optimal placement for the next query sequence and
     * returns the edge-split details, but does NOT insert the tip into the
     * GPU tree yet.  Call commitPlacement() to finalise.
     *
     * The returned PlacementResult reflects the edge that DIPPER recommends.
     * The caller may override this choice in commitPlacement() by supplying
     * different node labels.
     *
     * @throws std::logic_error if commitPlacement() has not been called after
     *         a previous findNextPlacement(), or if all queries are placed.
     */
    PlacementResult findNextPlacement()
    {
        _requireOneOf({State::GPUReady}, "findNextPlacement");
        _ensureIterativePipeline("findNextPlacement");
        if (m_placementPending)
            throw std::logic_error(
                "[DIPPER] findNextPlacement: a placement is already pending. "
                "Call commitPlacement() first.");
        if (m_kpArrays.m_nextQueryIdx >= (int)m_numSequences)
            throw std::logic_error(
                "[DIPPER] findNextPlacement: all query sequences have already been placed.");

        const int seqIdx    = m_kpArrays.m_nextQueryIdx;
        const int seedCount = (m_tree != nullptr) ? (int)m_backboneSize : 2;
        const int queryNum  = seqIdx - seedCount + 1;
        const int totalQ    = (int)m_numSequences - seedCount;

        std::cerr << "[DIPPER] Finding placement for query " << queryNum << " / " << totalQ
                  << " (global idx " << seqIdx
                  << ", name: " << m_sortedNames[seqIdx] << ")...\n";

        MashPlacement::KPlacementDeviceArrays::PlacementInfo info =
            m_kpArrays.findBestEdge(seqIdx, *m_params, m_mashArrays, m_matrixReader, m_msaArrays);

        // Pre-assign sequential labels (fixed regardless of which edge is chosen).
        const int iLabel = m_nextLabel;
        const int tLabel = m_nextLabel + 1;

        PlacementResult result;
        result.more            = (seqIdx + 1 < (int)m_numSequences);
        result.splitNodeA      = info.splitNodeA;
        result.splitNodeB      = info.splitNodeB;
        result.splitLabelA     = nodeLabel(info.splitNodeA);
        result.splitLabelB     = nodeLabel(info.splitNodeB);
        result.internalNodeIdx = info.internalNodeIdx;
        result.tipNodeIdx      = info.tipNodeIdx;
        result.internalLabel   = iLabel;
        result.tipLabel        = tLabel;
        result.tipName         = m_sortedNames[seqIdx];
        result.lenA            = info.lenA;
        result.lenB            = info.lenB;
        result.lenTip          = info.lenTip;

        m_pendingResult    = result;
        m_placementPending = true;
        return result;
    }

    /**
     * Phase 2 of two-phase placement.
     *
     * Inserts the pending query sequence into the GPU tree on the edge whose
     * two endpoints carry the given sequential labels.  The labels may be the
     * same pair returned by findNextPlacement() (confirming DIPPER's choice)
     * or a different pair (overriding it — the chosen edge must exist in the
     * current tree).
     *
     * @param labelA  Sequential label of one endpoint of the target edge.
     * @param labelB  Sequential label of the other endpoint.
     * @return        Final PlacementResult for this query, reflecting the
     *                actual edge used (branch lengths/split info may differ
     *                from findNextPlacement() if an override was requested).
     *
     * @throws std::logic_error if findNextPlacement() has not been called, or
     *         if the edge specified by (labelA, labelB) cannot be found.
     */
    PlacementResult commitPlacement(int labelA, int labelB)
    {
        _requireOneOf({State::GPUReady}, "commitPlacement");
        _ensureIterativePipeline("commitPlacement");
        if (!m_placementPending)
            throw std::logic_error(
                "[DIPPER] commitPlacement: no pending placement. "
                "Call findNextPlacement() first.");

        const int seqIdx = m_kpArrays.m_pendingSeqIdx;

        // Map sequential labels → DIPPER node indices.
        auto itA = m_labelToDipper.find(labelA);
        auto itB = m_labelToDipper.find(labelB);
        if (itA == m_labelToDipper.end() || itB == m_labelToDipper.end())
            throw std::logic_error(
                "[DIPPER] commitPlacement: one or both labels are unknown: "
                + std::to_string(labelA) + ", " + std::to_string(labelB));

        const int dipperA = itA->second;
        const int dipperB = itB->second;

        // Compare against the split edge already stored in m_pendingResult
        // (populated by findBestEdge → getEdgeInfo, no extra GPU reads needed).
        const int defaultDipperA = m_pendingResult.splitNodeA;
        const int defaultDipperB = m_pendingResult.splitNodeB;

        const bool isOverride =
            !((dipperA == defaultDipperA && dipperB == defaultDipperB) ||
              (dipperA == defaultDipperB && dipperB == defaultDipperA));

        int overrideEid = -1;
        MashPlacement::KPlacementDeviceArrays::PlacementInfo finalInfo = {};

        if (isOverride) {
            overrideEid = m_kpArrays.findEdgeBetween(dipperA, dipperB);
            if (overrideEid < 0)
                throw std::logic_error(
                    "[DIPPER] commitPlacement: no edge found between DIPPER nodes "
                    + std::to_string(dipperA) + " and " + std::to_string(dipperB)
                    + " (labels " + std::to_string(labelA) + " / " + std::to_string(labelB) + ").");
            finalInfo = m_kpArrays.getEdgeInfo(overrideEid, seqIdx);
        } else {
            finalInfo.splitNodeA      = m_pendingResult.splitNodeA;
            finalInfo.splitNodeB      = m_pendingResult.splitNodeB;
            finalInfo.internalNodeIdx = m_pendingResult.internalNodeIdx;
            finalInfo.tipNodeIdx      = m_pendingResult.tipNodeIdx;
            finalInfo.lenA            = m_pendingResult.lenA;
            finalInfo.lenB            = m_pendingResult.lenB;
            finalInfo.lenTip          = m_pendingResult.lenTip;
        }

        // Commit the placement on the GPU.
        m_kpArrays.commitQuery(seqIdx, overrideEid);

        // Register the two new nodes in both label maps.
        const int iLabel = m_pendingResult.internalLabel;
        const int tLabel = m_pendingResult.tipLabel;
        m_dipperToLabel[finalInfo.internalNodeIdx] = iLabel;
        m_dipperToLabel[finalInfo.tipNodeIdx]      = tLabel;
        m_labelToDipper[iLabel]                    = finalInfo.internalNodeIdx;
        m_labelToDipper[tLabel]                    = finalInfo.tipNodeIdx;
        m_nextLabel += 2;

        const bool more = (m_kpArrays.m_nextQueryIdx < (int)m_numSequences);

        PlacementResult result;
        result.more            = more;
        result.splitNodeA      = finalInfo.splitNodeA;
        result.splitNodeB      = finalInfo.splitNodeB;
        result.splitLabelA     = nodeLabel(finalInfo.splitNodeA);
        result.splitLabelB     = nodeLabel(finalInfo.splitNodeB);
        result.internalNodeIdx = finalInfo.internalNodeIdx;
        result.tipNodeIdx      = finalInfo.tipNodeIdx;
        result.internalLabel   = iLabel;
        result.tipLabel        = tLabel;
        result.tipName         = m_pendingResult.tipName;
        result.lenA            = finalInfo.lenA;
        result.lenB            = finalInfo.lenB;
        result.lenTip          = finalInfo.lenTip;

        m_placementPending = false;
        m_pendingResult    = {};

        if (!more) {
            m_state = State::Done;
            std::cerr << "[DIPPER] All query sequences placed.\n";
        }
        return result;
    }

    std::vector<RankedPlacementCandidate> exportPendingTopKPlacements(int k) const
    {
        _requireOneOf({State::GPUReady}, "exportPendingTopKPlacements");
        _ensureIterativePipeline("exportPendingTopKPlacements");
        if (!m_placementPending) {
            throw std::logic_error(
                "[DIPPER] exportPendingTopKPlacements: call findNextPlacement first.");
        }
        if (k <= 0) {
            throw std::invalid_argument(
                "[DIPPER] exportPendingTopKPlacements: k must be positive.");
        }
        std::vector<RankedPlacementCandidate> out;
        for (const auto& candidate : m_kpArrays.getTopKEdgeInfo(
                 m_kpArrays.m_pendingSeqIdx, k)) {
            out.push_back({
                nodeLabel(candidate.first.splitNodeA),
                nodeLabel(candidate.first.splitNodeB),
                candidate.second,
            });
        }
        return out;
    }

    // =========================================================================
    // Output
    // =========================================================================

    /**
     * Write the placement tree in Newick format to the file at @p outputPath.
     *
     * Can be called at any point after initializeGPU(), even before all
     * queries are placed, to inspect intermediate results.
     */
    void writeTree(const std::string& outputPath)
    {
        if (m_state < State::GPUReady)
            throw std::logic_error(
                "[DIPPER] writeTree: initializeGPU() must be called first.");
        std::ofstream ofs(outputPath);
        if (!ofs)
            throw std::runtime_error(
                "[DIPPER] writeTree: cannot open output file: " + outputPath);
        if (m_activePipeline == ActivePipeline::DivideAndConquer) {
            if (!m_dcTreeBuilt) {
                throw std::logic_error(
                    "[DIPPER] writeTree: divide-and-conquer tree has not been built yet. "
                    "Call buildDivideAndConquerTree() first.");
            }
            auto& dc_arrays =
                const_cast<MashPlacement::KPlacementDeviceArraysDC&>(
                    m_kpArraysDC);
            dc_arrays.printTreeDC(_leafNamesForCurrentSession(), ofs, m_dcRootNode);
        } else {
            _ensureIterativePipeline("writeTree");
            m_kpArrays.printTree(m_sortedNames, ofs);
        }
        std::cerr << "[DIPPER] Tree written to " << outputPath << "\n";
    }

    /** Write the placement tree to an already-open output stream. */
    void writeTree(std::ofstream& ofs)
    {
        if (m_state < State::GPUReady)
            throw std::logic_error(
                "[DIPPER] writeTree: initializeGPU() must be called first.");
        if (m_activePipeline == ActivePipeline::DivideAndConquer) {
            if (!m_dcTreeBuilt) {
                throw std::logic_error(
                    "[DIPPER] writeTree: divide-and-conquer tree has not been built yet. "
                    "Call buildDivideAndConquerTree() first.");
            }
            auto& dc_arrays =
                const_cast<MashPlacement::KPlacementDeviceArraysDC&>(
                    m_kpArraysDC);
            dc_arrays.printTreeDC(_leafNamesForCurrentSession(), ofs, m_dcRootNode);
        } else {
            _ensureIterativePipeline("writeTree");
            m_kpArrays.printTree(m_sortedNames, ofs);
        }
    }

    // =========================================================================
    // Convenience: place ALL remaining queries in one call
    // =========================================================================

    /**
     * Place every remaining query sequence in insertion order.
     *
     * Equivalent to: while (placeNextSequence()) {}
     */
    void placeAllSequences()
    {
        _ensureIterativePipeline("placeAllSequences");
        while (m_kpArrays.m_nextQueryIdx < (int)m_numSequences)
            (void)placeNextSequence();
    }

    /** Run the full legacy divide-and-conquer tree build. */
    void buildDivideAndConquerTree()
    {
        _requireOneOf({State::GPUReady, State::Done}, "buildDivideAndConquerTree");
        _ensureDivideAndConquerPipeline("buildDivideAndConquerTree");
        if (m_dcTreeBuilt) {
            return;
        }

        _ensureDivideAndConquerClustersReady();

        std::vector<int> large_clusters_idx;
        m_kpArraysDC.findClusterTreeDC(
            *m_params,
            m_mashArraysDC,
            m_matrixReader,
            m_msaArraysDC,
            m_hostDC,
            large_clusters_idx);
        m_dcTreeBuilt = true;
        m_state = State::Done;
        std::cerr << "[DIPPER] Divide-and-conquer tree build completed.\n";
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    /** Number of leaves in the backbone tree (0 when no backbone was given). */
    size_t backboneSize()   const { return _currentBackboneSize(); }

    /** True when no backbone tree was provided (from-scratch mode). */
    bool   hasBackbone()    const { return _currentBackboneSize() > 0; }

    /** Total number of sequences (backbone + query, or all sequences). */
    size_t numSequences()   const { return m_numSequences; }

    /** Total number of nodes in the backbone tree (leaves + internal nodes).
     *
     *  This is the correct starting value for an external N counter:
     *    - Without backbone: returns 2 (the two seed nodes).
     *    - With backbone:    returns allNodes.size() from the parsed Newick.
     *
     *  Each call to placeNextSequence() adds 2 nodes, so after k placements
     *  the external tree has backboneNodeCount() + 2*k nodes. */
    size_t backboneNodeCount() const
    {
        if (m_activePipeline == ActivePipeline::DivideAndConquer) {
            throw std::logic_error(
                "[DIPPER] backboneNodeCount: not defined for divide-and-conquer sessions.");
        }
        if (m_tree) return m_tree->allNodes.size();
        return 2;   // no-backbone mode: two seed nodes
    }

    /** Global index of the next sequence that the active scheduler would place. */
    int nextQueryIndex() const
    {
        if (m_activePipeline == ActivePipeline::DivideAndConquer) {
            if (!hasPendingQuery()) {
                return -1;
            }
            _ensureDivideAndConquerPendingStepReady();
            return m_dcPendingQueryInfo.sequenceIdx;
        }
        return m_kpArrays.m_nextQueryIdx;
    }

    /** True if at least one query sequence still remains to be placed. */
    bool hasPendingQuery() const
    {
        if (m_activePipeline == ActivePipeline::DivideAndConquer) {
            if (m_state != State::GPUReady && m_state != State::Done) {
                return false;
            }
            if (m_dcTreeBuilt) {
                return false;
            }
            _ensureDivideAndConquerSchedulerReady();
            const size_t scheduled_count = m_dcScheduledQueryNames.empty()
                ? m_numSequences - m_dcBackboneSize
                : m_dcScheduledQueryNames.size();
            return m_dcProcessedQueryCount < scheduled_count;
        }
        if (m_state == State::SeqsLoaded) {
            const int initial_query_idx =
                (m_tree != nullptr) ? static_cast<int>(m_backboneSize) : 2;
            return initial_query_idx < static_cast<int>(m_numSequences);
        }

        if (m_state == State::GPUReady || m_state == State::Done) {
            return m_kpArrays.m_nextQueryIdx >= 0 &&
                   m_kpArrays.m_nextQueryIdx < static_cast<int>(m_numSequences);
        }

        return false;
    }

    /** Backward-compatible alias for checking whether any query remains. */
    bool pending_query() const { return hasPendingQuery(); }

    int queryClusterId(int sequence_idx) const
    {
        _ensureQueryClusterAssignments();
        if (sequence_idx < 0 ||
            sequence_idx >= static_cast<int>(m_queryClusterIds.size())) {
            throw std::out_of_range(
                "[DIPPER] queryClusterId: sequence index out of range.");
        }
        if (static_cast<size_t>(sequence_idx) < m_queryClusterBackboneSize) {
            throw std::invalid_argument(
                "[DIPPER] queryClusterId: backbone sequences do not have query cluster IDs.");
        }
        return m_queryClusterIds[static_cast<size_t>(sequence_idx)];
    }

    /** Export the backbone alignment in DIPPER's current backbone-tip order. */
    std::vector<NamedSequence> exportBackboneAlignment() const
    {
        _requireOneOf(
            {State::SeqsLoaded, State::GPUReady, State::Done},
            "exportBackboneAlignment");
        const size_t current_backbone_size = _currentBackboneSize();
        if (current_backbone_size == 0) {
            throw std::logic_error(
                "[DIPPER] exportBackboneAlignment: no backbone partition is available.");
        }
        if (m_sortedNames.size() < current_backbone_size) {
            throw std::runtime_error(
                "[DIPPER] exportBackboneAlignment: sortedNames() is smaller than backboneSize().");
        }

        const auto& raw_name_to_index =
            _buildRawNameToIndexMap("exportBackboneAlignment");
        std::vector<NamedSequence> out;
        out.reserve(current_backbone_size);
        for (size_t idx = 0; idx < current_backbone_size; ++idx) {
            const std::string& name = m_sortedNames[idx];
            const auto it = raw_name_to_index.find(name);
            if (it == raw_name_to_index.end()) {
                throw std::runtime_error(
                    "[DIPPER] exportBackboneAlignment: missing sequence for backbone tip '" +
                    name + "'.");
            }
            out.push_back(NamedSequence{
                name,
                m_rawSeqs[it->second],
            });
        }
        return out;
    }

    std::vector<std::pair<std::string, int>> exportQueryClusterAssignments() const
    {
        _ensureQueryClusterAssignments();
        std::vector<std::pair<std::string, int>> out;
        out.reserve(m_numSequences - m_queryClusterBackboneSize);
        for (size_t idx = m_queryClusterBackboneSize; idx < m_numSequences; ++idx) {
            out.emplace_back(
                m_sortedNames[idx],
                m_queryClusterIds[idx]);
        }
        return out;
    }

    /**
     * Request the query-clustering candidate depth before the first DC
     * clustering pass.  Calling exportQueryTopKBackboneCandidates() with a
     * larger k after clustering is valid, but necessarily recomputes the
     * full query-to-backbone ranking.  Large-backbone clients should set the
     * final depth before initialization so the first pass produces all of
     * the candidates they need.
     */
    void setQueryClusterTopK(int k)
    {
        if (k <= 0) {
            throw std::invalid_argument(
                "[DIPPER] setQueryClusterTopK: k must be positive.");
        }
        if (m_queryClusterIdsReady && k > m_computedQueryClusterTopK) {
            throw std::logic_error(
                "[DIPPER] setQueryClusterTopK: set the requested depth before "
                "the first clustering pass.");
        }
        m_requestedQueryClusterTopK =
            std::max(m_requestedQueryClusterTopK, k);
    }

    std::vector<QueryTopKBackboneCandidates>
    exportQueryTopKBackboneCandidates(int k) const
    {
        if (k <= 0) {
            throw std::invalid_argument(
                "[DIPPER] exportQueryTopKBackboneCandidates: k must be positive.");
        }
        m_requestedQueryClusterTopK = std::max(m_requestedQueryClusterTopK, k);
        _ensureDivideAndConquerClustersReady();
        _copyDivideAndConquerTreeStateToHost();
        if (m_dcRootNode < 0) {
            m_dcRootNode = static_cast<int>(m_numSequences);
        }
        _assignDivideAndConquerShadowLabels();

        std::vector<QueryTopKBackboneCandidates> out;
        out.reserve(m_numSequences - m_dcBackboneSize);
        for (size_t idx = m_dcBackboneSize; idx < m_numSequences; ++idx) {
            QueryTopKBackboneCandidates query;
            query.sequenceIdx = static_cast<int>(idx);
            query.name = m_sortedNames[idx];
            query.originalClusterId = m_queryClusterIds[idx];
            const auto& edge_ids = m_hostDC.clusterTopEdgeIds.at(idx);
            const auto& edge_distances =
                m_hostDC.clusterTopEdgeAdditionalDistances.at(idx);
            if (edge_ids.size() != edge_distances.size()) {
                throw std::runtime_error(
                    "[DIPPER] top-k edge and distance counts differ.");
            }
            const size_t keep = std::min(edge_ids.size(), static_cast<size_t>(k));
            query.candidates.reserve(keep);
            for (size_t rank = 0; rank < keep; ++rank) {
                const int edge_id = edge_ids[rank];
                _requireDivideAndConquerClosestEdgeIndex(
                    edge_id, "exportQueryTopKBackboneCandidates");
                const int node_a = m_hostDC.h_belong[edge_id];
                const int node_b = m_hostDC.h_e[edge_id];
                query.candidates.push_back(RankedDipperCandidate{
                    edge_id,
                    nodeLabel(node_a),
                    nodeLabel(node_b),
                    edge_distances[rank],
                });
            }
            out.push_back(std::move(query));
        }
        return out;
    }

    void setEffectiveQueryClusterAssignments(
        const std::vector<std::pair<std::string, int>>& assignments)
    {
        _ensureDivideAndConquerClustersReady();
        if (m_dcProcessedQueryCount != 0 || m_dcActiveClusterId >= 0 ||
            m_dcPendingStepReady) {
            throw std::logic_error(
                "[DIPPER] effective cluster assignments must be installed before query insertion.");
        }
        m_queryEffectiveClusterIds = m_queryClusterIds;
        std::unordered_map<std::string, size_t> query_index;
        for (size_t idx = m_dcBackboneSize; idx < m_numSequences; ++idx) {
            query_index.emplace(m_sortedNames[idx], idx);
        }
        const int cluster_slots = std::max(0, m_hostDC.numSequences * 4 - 4);
        for (const auto& [name, cluster_id] : assignments) {
            const auto it = query_index.find(name);
            if (it == query_index.end()) {
                throw std::invalid_argument(
                    "[DIPPER] effective assignment names an unknown query: " + name);
            }
            if (cluster_id < 0 || cluster_id >= cluster_slots) {
                throw std::invalid_argument(
                    "[DIPPER] effective assignment has an invalid cluster id.");
            }
            m_queryEffectiveClusterIds[it->second] = cluster_id;
        }
        if (m_dcSchedulerReady) {
            _resetDivideAndConquerSchedulerState();
        }
    }

    /** Restrict the D&C insertion scheduler to a named query shortlist. */
    void setDivideAndConquerScheduledQueries(
        const std::vector<std::string>& query_names)
    {
        _ensureDivideAndConquerClustersReady();
        if (m_dcProcessedQueryCount != 0 || m_dcActiveClusterId >= 0 ||
            m_dcPendingStepReady) {
            throw std::logic_error(
                "[DIPPER] scheduled queries must be installed before query insertion.");
        }
        if (query_names.empty()) {
            throw std::invalid_argument(
                "[DIPPER] scheduled query shortlist must not be empty.");
        }
        std::unordered_set<std::string> available;
        available.reserve(m_numSequences - m_dcBackboneSize);
        for (size_t idx = m_dcBackboneSize; idx < m_numSequences; ++idx) {
            available.insert(m_sortedNames[idx]);
        }
        m_dcScheduledQueryNames.clear();
        m_dcScheduledQueryNames.reserve(query_names.size());
        for (const std::string& name : query_names) {
            if (available.count(name) == 0) {
                throw std::invalid_argument(
                    "[DIPPER] scheduled query shortlist names an unknown query: " +
                    name);
            }
            m_dcScheduledQueryNames.insert(name);
        }
        _resetDivideAndConquerSchedulerState();
    }

    void setBackboneTipStableLabels(
        const std::vector<std::pair<std::string, int>>& labels)
    {
        m_dcBackboneTipStableLabels.clear();
        for (const auto& [name, label] : labels) {
            if (label < 0 ||
                !m_dcBackboneTipStableLabels.emplace(name, label).second) {
                throw std::invalid_argument(
                    "[DIPPER] invalid or duplicate backbone stable label.");
            }
        }
    }

    /** Return the current D&C scheduler order without placing any query. */
    std::vector<std::string> exportDivideAndConquerScheduledQueryNames(
        size_t limit) const
    {
        _ensureDivideAndConquerSchedulerReady();
        std::vector<std::string> out;
        const size_t scheduled_count = m_dcScheduledQueryNames.empty()
            ? m_numSequences - m_dcBackboneSize
            : m_dcScheduledQueryNames.size();
        out.reserve(std::min(limit, scheduled_count));
        for (int cluster_id : m_dcOrderedClusters) {
            for (int query_idx :
                 m_dcClusterQueries[static_cast<size_t>(cluster_id)]) {
                if (out.size() >= limit) {
                    return out;
                }
                out.push_back(m_sortedNames.at(
                    static_cast<size_t>(query_idx)));
            }
        }
        return out;
    }

    std::vector<std::pair<std::string, int>>
    exportEffectiveQueryClusterAssignments() const
    {
        _ensureDivideAndConquerClustersReady();
        const auto& ids = m_queryEffectiveClusterIds.empty()
            ? m_queryClusterIds
            : m_queryEffectiveClusterIds;
        std::vector<std::pair<std::string, int>> out;
        out.reserve(m_numSequences - m_dcBackboneSize);
        for (size_t idx = m_dcBackboneSize; idx < m_numSequences; ++idx) {
            out.emplace_back(m_sortedNames[idx], ids[idx]);
        }
        return out;
    }

    /** Export the next pending query sequence together with its cluster id. */
    PendingQueryInfo exportPendingQuery() const
    {
        _requireOneOf(
            {State::GPUReady, State::Done},
            "exportPendingQuery");
        if (m_activePipeline == ActivePipeline::DivideAndConquer) {
            if (!hasPendingQuery()) {
                throw std::logic_error(
                    "[DIPPER] exportPendingQuery: no pending divide-and-conquer query remains.");
            }
            _ensureDivideAndConquerPendingStepReady();
            return m_dcPendingQueryInfo;
        }
        if (!hasPendingQuery()) {
            throw std::logic_error(
                "[DIPPER] exportPendingQuery: no pending query remains.");
        }

        const int query_idx = nextQueryIndex();
        if (query_idx < 0) {
            throw std::runtime_error(
                "[DIPPER] exportPendingQuery: nextQueryIndex() returned a negative value.");
        }

        const size_t name_idx = static_cast<size_t>(query_idx);
        if (name_idx >= m_sortedNames.size()) {
            throw std::runtime_error(
                "[DIPPER] exportPendingQuery: nextQueryIndex() is out of range for sortedNames().");
        }

        const auto& raw_name_to_index =
            _buildRawNameToIndexMap("exportPendingQuery");
        const std::string& name = m_sortedNames[name_idx];
        const auto raw_it = raw_name_to_index.find(name);
        if (raw_it == raw_name_to_index.end()) {
            throw std::runtime_error(
                "[DIPPER] exportPendingQuery: missing sequence for pending query '" +
                name + "'.");
        }

        PendingQueryInfo out;
        out.sequenceIdx = query_idx;
        out.clusterId = queryClusterId(query_idx);
        out.name = name;
        out.sequence = m_rawSeqs[raw_it->second];
        return out;
    }

    /**
     * Advance the divide-and-conquer scheduler by one query placement step.
     *
     * When the optional committed split labels/lengths are provided, the host
     * shadow tree is updated using the caller's actual committed placement
     * instead of DIPPER's own pending candidate-edge split. This keeps the
     * divide-and-conquer scheduler synchronized with an external placer such as
     * MLIPPER.
     */
    DivideAndConquerAdvanceResult advanceDivideAndConquerStep(
        int committedSplitLabelA = -1,
        int committedSplitLabelB = -1,
        double committedLenA = -1.0,
        double committedLenTip = -1.0)
    {
        _requireOneOf(
            {State::GPUReady, State::Done},
            "advanceDivideAndConquerStep");
        _ensureDivideAndConquerPipeline("advanceDivideAndConquerStep");
        if (!hasPendingQuery()) {
            throw std::logic_error(
                "[DIPPER] advanceDivideAndConquerStep: no pending divide-and-conquer query remains.");
        }
        _ensureDivideAndConquerPendingStepReady();

        DivideAndConquerAdvanceResult advance_result;
        advance_result.tipName = m_dcPendingQueryInfo.name;

        int edge_id = m_dcPendingEdgeId;
        double frac_len = m_dcPendingFracLen;
        double add_len = m_dcPendingAddLen;
        if (committedSplitLabelA >= 0 || committedSplitLabelB >= 0 ||
            committedLenA >= 0.0 || committedLenTip >= 0.0) {
            if (committedSplitLabelA < 0 || committedSplitLabelB < 0) {
                throw std::invalid_argument(
                    "[DIPPER] advanceDivideAndConquerStep: committed split labels must both be provided.");
            }
            if (committedLenA < 0.0 || committedLenTip < 0.0) {
                throw std::invalid_argument(
                    "[DIPPER] advanceDivideAndConquerStep: committed branch lengths must be non-negative.");
            }
            const auto node_a_it =
                m_dcLabelToNode.find(committedSplitLabelA);
            const auto node_b_it =
                m_dcLabelToNode.find(committedSplitLabelB);
            if (node_a_it == m_dcLabelToNode.end() ||
                node_b_it == m_dcLabelToNode.end()) {
                throw std::runtime_error(
                    "[DIPPER] advanceDivideAndConquerStep: committed split labels are unavailable in the current divide-and-conquer shadow tree.");
            }
            _requireDivideAndConquerNodeIndex(
                node_a_it->second,
                "advanceDivideAndConquerStep(committedSplitLabelA)");
            _requireDivideAndConquerNodeIndex(
                node_b_it->second,
                "advanceDivideAndConquerStep(committedSplitLabelB)");
            edge_id =
                _findDivideAndConquerDirectedEdge(
                    node_a_it->second,
                    node_b_it->second);
            if (edge_id < 0) {
                throw std::runtime_error(
                    "[DIPPER] advanceDivideAndConquerStep: committed split edge does not match the current divide-and-conquer shadow tree.");
            }
            frac_len = committedLenA;
            add_len = committedLenTip;
        }

        const int new_internal_idx =
            static_cast<int>(m_numSequences) - 1 + m_dcActiveInsertLeafCount;
        advance_result.splitLabelA =
            nodeLabel(m_hostDC.h_belong[edge_id]);
        advance_result.splitLabelB =
            nodeLabel(m_hostDC.h_e[edge_id]);
        _requireDivideAndConquerPlacementUpdateCapacity(
            edge_id,
            m_dcPendingLeafId,
            new_internal_idx,
            "advanceDivideAndConquerStep");
        updateTreeStructureInClusterCpu(
            m_hostDC.h_head,
            m_hostDC.h_nxt,
            m_hostDC.h_e,
            m_hostDC.h_len,
            m_hostDC.h_closest_dis,
            m_hostDC.h_closest_id,
            m_hostDC.h_belong,
            edge_id,
            frac_len,
            add_len,
            m_dcPendingLeafId,
            m_dcActiveIdx,
            static_cast<int>(m_numSequences),
            m_dcActiveInsertLeafCount);
        _assignDivideAndConquerNodeLabel(new_internal_idx);
        _assignDivideAndConquerNodeLabel(m_dcPendingLeafId);
        advance_result.internalLabel = nodeLabel(new_internal_idx);
        advance_result.tipLabel = nodeLabel(m_dcPendingLeafId);

        m_dcActiveIdx += 4;
        _requireDivideAndConquerActiveClusterCapacity(
            4,
            1,
            "advanceDivideAndConquerStep");
        m_dcActiveLeafMap[static_cast<size_t>(m_dcActiveLeafCount)] =
            m_dcPendingLeafId;
        updateClusterInfoCpu(
            m_dcPendingLeafId,
            m_dcActiveIdx,
            m_dcActiveLeafMask,
            m_dcActiveEdgeMask,
            m_dcActiveEdgeCount,
            m_dcActiveLeafCount);
        if (!m_dcActiveId.empty()) {
            m_dcActiveId[0] = m_dcPendingLeafId;
        }
        if (!m_dcActiveFrom.empty()) {
            m_dcActiveFrom[0] = -1;
        }
        if (!m_dcActiveDis.empty()) {
            m_dcActiveDis[0] = 0.0;
        }
        const int active_cluster_edge_id =
            _resolveDivideAndConquerClusterAnchorEdge(m_dcActiveClusterId);
        _requireDivideAndConquerClosestEdgeIndex(
            active_cluster_edge_id,
            "advanceDivideAndConquerStep(active cluster edge)");
        updateClosestNodesInClusterCpu(
            m_hostDC.h_head,
            m_hostDC.h_nxt,
            m_hostDC.h_e,
            m_hostDC.h_len,
            m_hostDC.h_closest_dis,
            m_hostDC.h_closest_id,
            m_dcPendingLeafId,
            m_dcActiveId,
            m_dcActiveFrom,
            m_dcActiveDis,
            active_cluster_edge_id,
            m_hostDC.h_belong);

        m_dcActiveEdgeCount += 4;
        ++m_dcActiveLeafCount;
        ++m_dcActiveClusterOffset;
        ++m_dcActiveInsertLeafCount;
        ++m_dcProcessedQueryCount;
        if (m_dcPendingLeafId >= 0 &&
            m_dcPendingLeafId <
                static_cast<int>(m_dcInsertedQueryLeafIds.size())) {
            m_dcInsertedQueryLeafIds[static_cast<size_t>(m_dcPendingLeafId)] = 1;
        }
        m_dcPendingStepReady = false;
        m_dcPendingQueryInfo = {};
        m_dcPendingLeafId = -1;
        m_dcPendingEdgeId = -1;
        m_dcPendingFracLen = 0.0;
        m_dcPendingAddLen = 0.0;

        if (m_dcCurrentClusterCursor < m_dcOrderedClusters.size() &&
            m_dcActiveClusterOffset >=
                m_dcClusterQueries[static_cast<size_t>(m_dcActiveClusterId)].size()) {
            ++m_dcCurrentClusterCursor;
            m_dcActiveClusterId = -1;
            m_dcActiveClusterOffset = 0;
        }
        return advance_result;
    }

    void applyDivideAndConquerSprMove(
        int pruneLabel,
        int regraftSplitLabelA,
        int regraftSplitLabelB,
        double regraftLenA,
        double regraftLenTip)
    {
        _requireOneOf(
            {State::GPUReady, State::Done},
            "applyDivideAndConquerSprMove");
        _ensureDivideAndConquerPipeline("applyDivideAndConquerSprMove");
        _ensureDivideAndConquerSchedulerReady();
        if (pruneLabel < 0 || regraftSplitLabelA < 0 ||
            regraftSplitLabelB < 0) {
            throw std::invalid_argument(
                "[DIPPER] applyDivideAndConquerSprMove: labels must be non-negative.");
        }
        if (regraftLenA < 0.0 || regraftLenTip < 0.0) {
            throw std::invalid_argument(
                "[DIPPER] applyDivideAndConquerSprMove: branch lengths must be non-negative.");
        }

        const auto prune_it = m_dcLabelToNode.find(pruneLabel);
        const auto regraft_a_it =
            m_dcLabelToNode.find(regraftSplitLabelA);
        const auto regraft_b_it =
            m_dcLabelToNode.find(regraftSplitLabelB);
        if (prune_it == m_dcLabelToNode.end() ||
            regraft_a_it == m_dcLabelToNode.end() ||
            regraft_b_it == m_dcLabelToNode.end()) {
            throw std::runtime_error(
                "[DIPPER] applyDivideAndConquerSprMove: one or more labels are unavailable in the current shadow tree.");
        }
        _requireDivideAndConquerNodeIndex(
            prune_it->second,
            "applyDivideAndConquerSprMove(pruneLabel)");
        _requireDivideAndConquerNodeIndex(
            regraft_a_it->second,
            "applyDivideAndConquerSprMove(regraftSplitLabelA)");
        _requireDivideAndConquerNodeIndex(
            regraft_b_it->second,
            "applyDivideAndConquerSprMove(regraftSplitLabelB)");

        std::vector<int> parent;
        std::vector<int> parent_edge;
        _buildDivideAndConquerParentView(parent, parent_edge);

        const int prune_root = prune_it->second;
        if (prune_root < 0 ||
            prune_root >= static_cast<int>(parent.size())) {
            throw std::runtime_error(
                "[DIPPER] applyDivideAndConquerSprMove: prune label resolved outside the current shadow tree.");
        }
        if (parent[static_cast<size_t>(prune_root)] < 0 ||
            parent[static_cast<size_t>(prune_root)] == prune_root) {
            throw std::runtime_error(
                "[DIPPER] applyDivideAndConquerSprMove: cannot prune the current shadow root.");
        }

        const int prune_parent =
            parent[static_cast<size_t>(prune_root)];
        _requireDivideAndConquerNodeIndex(
            prune_parent,
            "applyDivideAndConquerSprMove(prune parent)");
        int sibling = -1;
        int sibling_scan_steps = 0;
        const int edge_slots = static_cast<int>(m_numSequences * 8);
        for (int edge_idx = m_hostDC.h_head[prune_parent];
             edge_idx != -1;
             edge_idx = m_hostDC.h_nxt[edge_idx]) {
            if (edge_idx < 0 || edge_idx >= edge_slots) {
                throw std::runtime_error(
                    "[DIPPER] applyDivideAndConquerSprMove: prune-parent adjacency edge is out of range.");
            }
            if (++sibling_scan_steps > edge_slots) {
                throw std::runtime_error(
                    "[DIPPER] applyDivideAndConquerSprMove: prune-parent adjacency list contains a cycle.");
            }
            const int neighbor = m_hostDC.h_e[edge_idx];
            if (neighbor < 0 ||
                neighbor == parent[static_cast<size_t>(prune_parent)] ||
                neighbor == prune_root) {
                continue;
            }
            sibling = neighbor;
            break;
        }
        if (sibling < 0) {
            throw std::runtime_error(
                "[DIPPER] applyDivideAndConquerSprMove: failed to resolve the prune sibling.");
        }
        int grandparent =
            parent[static_cast<size_t>(prune_parent)];
        if (grandparent == prune_parent) {
            grandparent = -1;
        }

        int regraft_parent = regraft_a_it->second;
        int regraft_child = regraft_b_it->second;
        int regraft_edge =
            _findDivideAndConquerDirectedEdge(
                regraft_parent,
                regraft_child);
        double proximal_length = regraftLenA;
        if (regraft_edge < 0 ||
            parent[static_cast<size_t>(regraft_child)] != regraft_parent) {
            regraft_parent = regraft_b_it->second;
            regraft_child = regraft_a_it->second;
            regraft_edge =
                _findDivideAndConquerDirectedEdge(
                    regraft_parent,
                    regraft_child);
            if (regraft_edge < 0 ||
                parent[static_cast<size_t>(regraft_child)] != regraft_parent) {
                throw std::runtime_error(
                    "[DIPPER] applyDivideAndConquerSprMove: regraft split labels do not match the current rooted shadow tree.");
            }
            proximal_length =
                std::max(
                    0.0,
                    m_hostDC.h_len[regraft_edge] - regraftLenA);
        }
        _requireDivideAndConquerClosestEdgeIndex(
            regraft_edge,
            "applyDivideAndConquerSprMove(regraft edge)");

        const double target_edge_length = m_hostDC.h_len[regraft_edge];
        if (proximal_length > target_edge_length) {
            proximal_length = target_edge_length;
        }
        const double distal_length =
            std::max(0.0, target_edge_length - proximal_length);
        const int sibling_edge =
            _findDivideAndConquerDirectedEdge(sibling, prune_parent);
        if (sibling_edge < 0) {
            throw std::runtime_error(
                "[DIPPER] applyDivideAndConquerSprMove: sibling reverse edge is unavailable.");
        }
        _requireDivideAndConquerClosestEdgeIndex(
            sibling_edge,
            "applyDivideAndConquerSprMove(sibling edge)");
        const double sibling_edge_length =
            m_hostDC.h_len[sibling_edge];
        int parent_edge_to_grandparent = -1;
        if (grandparent >= 0) {
            parent_edge_to_grandparent =
                _findDivideAndConquerDirectedEdge(
                    prune_parent,
                    grandparent);
            if (parent_edge_to_grandparent < 0) {
                throw std::runtime_error(
                    "[DIPPER] applyDivideAndConquerSprMove: prune-parent to grandparent edge is unavailable.");
            }
            _requireDivideAndConquerClosestEdgeIndex(
                parent_edge_to_grandparent,
                "applyDivideAndConquerSprMove(parent edge)");
        }
        const double parent_edge_length =
            (grandparent >= 0)
                ? m_hostDC.h_len[parent_edge_to_grandparent]
                : 0.0;

        std::vector<int> reusable_edges;
        const auto require_directed_edge =
            [&](int from_node, int to_node, const char* context) -> int {
                const int edge_idx =
                    _findDivideAndConquerDirectedEdge(from_node, to_node);
                if (edge_idx < 0) {
                    throw std::runtime_error(
                        std::string("[DIPPER] ") + context +
                        ": required directed edge is unavailable before SPR rewiring.");
                }
                _requireDivideAndConquerClosestEdgeIndex(edge_idx, context);
                return edge_idx;
            };
        reusable_edges.push_back(
            require_directed_edge(
                prune_parent,
                prune_root,
                "applyDivideAndConquerSprMove(prune parent -> root)"));
        reusable_edges.push_back(
            require_directed_edge(
                prune_root,
                prune_parent,
                "applyDivideAndConquerSprMove(prune root -> parent)"));
        reusable_edges.push_back(
            require_directed_edge(
                prune_parent,
                sibling,
                "applyDivideAndConquerSprMove(prune parent -> sibling)"));
        reusable_edges.push_back(
            require_directed_edge(
                sibling,
                prune_parent,
                "applyDivideAndConquerSprMove(sibling -> prune parent)"));
        if (grandparent >= 0) {
            reusable_edges.push_back(
                require_directed_edge(
                    grandparent,
                    prune_parent,
                    "applyDivideAndConquerSprMove(grandparent -> prune parent)"));
            reusable_edges.push_back(
                require_directed_edge(
                    prune_parent,
                    grandparent,
                    "applyDivideAndConquerSprMove(prune parent -> grandparent)"));
        }
        reusable_edges.push_back(
            require_directed_edge(
                regraft_parent,
                regraft_child,
                "applyDivideAndConquerSprMove(regraft parent -> child)"));
        reusable_edges.push_back(
            require_directed_edge(
                regraft_child,
                regraft_parent,
                "applyDivideAndConquerSprMove(regraft child -> parent)"));
        {
            std::vector<int> sorted_edges = reusable_edges;
            std::sort(sorted_edges.begin(), sorted_edges.end());
            if (std::adjacent_find(
                    sorted_edges.begin(),
                    sorted_edges.end()) != sorted_edges.end()) {
                throw std::runtime_error(
                    "[DIPPER] applyDivideAndConquerSprMove: SPR rewiring would reuse the same edge slot more than once.");
            }
        }
        // Capture the complete pruned component before changing its boundary.
        // A local SPR repair unit may be a multi-query subtree, not one tip.
        std::vector<int> moved_nodes{prune_root};
        std::vector<int> moved_parents{prune_parent};
        std::vector<int> moved_leaves;
        for (size_t i = 0; i < moved_nodes.size(); ++i) {
            const int node = moved_nodes[i];
            const int parent = moved_parents[i];
            int degree = 0;
            for (int edge_idx = m_hostDC.h_head[node];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                ++degree;
                const int neighbor = m_hostDC.h_e[edge_idx];
                if (neighbor != parent) {
                    moved_nodes.push_back(neighbor);
                    moved_parents.push_back(node);
                }
            }
            if (degree == 1 &&
                node >= 0 && node < static_cast<int>(m_numSequences)) {
                moved_leaves.push_back(node);
            }
        }
        reusable_edges.clear();
        std::vector<int> rewired_edges;
        auto stash_edge = [&](int edge_idx, const char* context) {
            if (edge_idx < 0) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": failed to unlink required edge during SPR rewiring.");
            }
            _requireDivideAndConquerClosestEdgeIndex(edge_idx, context);
            reusable_edges.push_back(edge_idx);
        };
        stash_edge(
            _unlinkDivideAndConquerDirectedEdge(
                prune_parent,
                prune_root),
            "applyDivideAndConquerSprMove(unlink prune parent -> root)");
        stash_edge(
            _unlinkDivideAndConquerDirectedEdge(
                prune_root,
                prune_parent),
            "applyDivideAndConquerSprMove(unlink prune root -> parent)");
        stash_edge(
            _unlinkDivideAndConquerDirectedEdge(
                prune_parent,
                sibling),
            "applyDivideAndConquerSprMove(unlink prune parent -> sibling)");
        stash_edge(
            _unlinkDivideAndConquerDirectedEdge(
                sibling,
                prune_parent),
            "applyDivideAndConquerSprMove(unlink sibling -> prune parent)");
        if (grandparent >= 0) {
            stash_edge(
                _unlinkDivideAndConquerDirectedEdge(
                    grandparent,
                    prune_parent),
                "applyDivideAndConquerSprMove(unlink grandparent -> prune parent)");
            stash_edge(
                _unlinkDivideAndConquerDirectedEdge(
                    prune_parent,
                    grandparent),
                "applyDivideAndConquerSprMove(unlink prune parent -> grandparent)");
        }
        stash_edge(
            _unlinkDivideAndConquerDirectedEdge(
                regraft_parent,
                regraft_child),
            "applyDivideAndConquerSprMove(unlink regraft parent -> child)");
        stash_edge(
            _unlinkDivideAndConquerDirectedEdge(
                regraft_child,
                regraft_parent),
            "applyDivideAndConquerSprMove(unlink regraft child -> parent)");

        const auto take_reusable_edge = [&]() -> int {
            if (reusable_edges.empty()) {
                throw std::runtime_error(
                    "[DIPPER] applyDivideAndConquerSprMove: no reusable shadow-tree edge slot is available.");
            }
            const int edge_idx = reusable_edges.back();
            reusable_edges.pop_back();
            return edge_idx;
        };
        const auto attach_rewired_edge = [&] (
            int from_node, int to_node, double edge_length) {
            const int edge_idx = take_reusable_edge();
            _attachDivideAndConquerDirectedEdge(
                edge_idx, from_node, to_node, edge_length);
            rewired_edges.push_back(edge_idx);
        };

        if (grandparent >= 0) {
            attach_rewired_edge(
                grandparent,
                sibling,
                parent_edge_length + sibling_edge_length);
            attach_rewired_edge(
                sibling,
                grandparent,
                parent_edge_length + sibling_edge_length);
        } else {
            m_dcRootNode = sibling;
        }

        attach_rewired_edge(
            regraft_parent,
            prune_parent,
            proximal_length);
        attach_rewired_edge(
            prune_parent,
            regraft_parent,
            proximal_length);
        attach_rewired_edge(
            prune_parent,
            regraft_child,
            distal_length);
        attach_rewired_edge(
            regraft_child,
            prune_parent,
            distal_length);
        attach_rewired_edge(
            prune_parent,
            prune_root,
            regraftLenTip);
        attach_rewired_edge(
            prune_root,
            prune_parent,
            regraftLenTip);

        _validateDivideAndConquerShadowTree(
            "applyDivideAndConquerSprMove(after_rewire)");
        const bool skip_recompute_closest =
            std::getenv("MLIPPER_DEBUG_SKIP_DIPPER_DC_SPR_CLOSEST") != nullptr;
        if (!skip_recompute_closest) {
            _recomputeDivideAndConquerClosestLeaves();
        }
        m_dcPendingStepReady = false;
        m_dcPendingQueryInfo = {};
        m_dcPendingLeafId = -1;
        m_dcPendingEdgeId = -1;
        m_dcPendingFracLen = 0.0;
        m_dcPendingAddLen = 0.0;
        const bool skip_active_rebuild =
            std::getenv("MLIPPER_DEBUG_SKIP_DIPPER_DC_SPR_ACTIVE_REBUILD") != nullptr;
        if (!skip_active_rebuild) {
            _rebuildDivideAndConquerActiveClusterState();
        }
    }

    /**
     * Convert a DIPPER internal node index to its sequential pre-order label.
     *
     * Valid for any DIPPER index that has appeared in PlacementResult
     * (splitNodeA, splitNodeB, internalNodeIdx, tipNodeIdx) or in the
     * backbone tree nodes.
     *
     * @return  The 0-based sequential label, or -1 if the index is unknown.
     */
    int nodeLabel(int dipperIdx) const
    {
        if (m_activePipeline == ActivePipeline::DivideAndConquer) {
            auto dc_it = m_dcNodeToLabel.find(dipperIdx);
            if (dc_it != m_dcNodeToLabel.end()) {
                return dc_it->second;
            }
        }
        auto it = m_dipperToLabel.find(dipperIdx);
        return (it != m_dipperToLabel.end()) ? it->second : -1;
    }

    /** Names in sorted GPU order.
     *
     *  Layout (with backbone):
     *    [0 .. backboneSize-1]          backbone leaf names
     *    [backboneSize .. numSeqs-1]    query names (not yet placed)
     *    [numSeqs .. ]                  internal backbone node names (e.g. "node_1")
     *
     *  Accessing splitNodeA/splitNodeB from PlacementResult is always safe
     *  because _buildIdMapWithBackbone() ensures internal nodes are included. */
    const std::vector<std::string>& sortedNames() const { return m_sortedNames; }

    /** Export the current divide-and-conquer device tree as a Newick string. */
    std::string exportCurrentDivideAndConquerTreeNewickString(
        bool print_binary_newick = true) const
    {
        _requireOneOf(
            {State::GPUReady, State::Done},
            "exportCurrentDivideAndConquerTreeNewickString");
        _ensureDivideAndConquerPipeline(
            "exportCurrentDivideAndConquerTreeNewickString");
        _ensureDivideAndConquerClustersReady();

        char temp_path[] = "/tmp/dipper_session_dc_tree_XXXXXX";
        const int fd = mkstemp(temp_path);
        if (fd < 0) {
            throw std::runtime_error(
                "[DIPPER] exportCurrentDivideAndConquerTreeNewickString: mkstemp() failed.");
        }
        close(fd);

        const bool previous_binary_mode = MashPlacement::g_printBinaryNewick;
        try {
            MashPlacement::g_printBinaryNewick = print_binary_newick;
            std::ofstream ofs(temp_path);
            if (!ofs) {
                throw std::runtime_error(
                    "[DIPPER] exportCurrentDivideAndConquerTreeNewickString: cannot open temporary file.");
            }
            auto& dc_arrays =
                const_cast<MashPlacement::KPlacementDeviceArraysDC&>(
                    m_kpArraysDC);
            dc_arrays.printTreeDC(_leafNamesForCurrentSession(), ofs, m_dcRootNode);
            ofs.close();
            MashPlacement::g_printBinaryNewick = previous_binary_mode;

            std::ifstream ifs(temp_path);
            if (!ifs) {
                throw std::runtime_error(
                    "[DIPPER] exportCurrentDivideAndConquerTreeNewickString: failed to reopen temporary tree file.");
            }

            std::string newick;
            std::getline(ifs, newick);
            std::remove(temp_path);
            if (newick.empty()) {
                throw std::runtime_error(
                    "[DIPPER] exportCurrentDivideAndConquerTreeNewickString: exported Newick tree is empty.");
            }
            return newick;
        } catch (...) {
            MashPlacement::g_printBinaryNewick = previous_binary_mode;
            std::remove(temp_path);
            throw;
        }
    }

    /** Export the interleaved D&C host shadow tree, including committed queries. */
    std::string exportCurrentDivideAndConquerShadowTreeNewickString() const
    {
        _requireOneOf(
            {State::GPUReady, State::Done},
            "exportCurrentDivideAndConquerShadowTreeNewickString");
        _ensureDivideAndConquerPipeline(
            "exportCurrentDivideAndConquerShadowTreeNewickString");
        _ensureDivideAndConquerSchedulerReady();
        _validateDivideAndConquerShadowTree(
            "exportCurrentDivideAndConquerShadowTreeNewickString");

        std::function<void(int, int, std::ostringstream&)> write_node;
        write_node = [&](int node, int parent, std::ostringstream& out) {
            std::vector<int> child_edges;
            for (int edge_idx = m_hostDC.h_head[node];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                if (m_hostDC.h_e[edge_idx] != parent) {
                    child_edges.push_back(edge_idx);
                }
            }
            if (child_edges.empty()) {
                if (node < 0 ||
                    node >= static_cast<int>(m_numSequences)) {
                    throw std::runtime_error(
                        "[DIPPER] shadow-tree leaf lacks a sequence name.");
                }
                out << m_sortedNames[static_cast<size_t>(node)];
                return;
            }
            out << '(';
            for (size_t idx = 0; idx < child_edges.size(); ++idx) {
                if (idx != 0) {
                    out << ',';
                }
                const int edge_idx = child_edges[idx];
                write_node(m_hostDC.h_e[edge_idx], node, out);
                out << ':' << std::setprecision(17)
                    << m_hostDC.h_len[edge_idx];
            }
            out << ')';
        };

        std::ostringstream out;
        write_node(m_dcRootNode, -1, out);
        out << ';';
        return out.str();
    }

    std::vector<std::pair<int, std::string>>
    exportDivideAndConquerStableLabelCladeSignatures() const
    {
        _requireOneOf(
            {State::GPUReady, State::Done},
            "exportDivideAndConquerStableLabelCladeSignatures");
        _ensureDivideAndConquerPipeline(
            "exportDivideAndConquerStableLabelCladeSignatures");
        _ensureDivideAndConquerSchedulerReady();
        const int edge_slots = static_cast<int>(m_numSequences) * 8;
        std::vector<std::pair<int, std::string>> out;
        std::function<std::vector<std::string>(int, int)> dfs =
            [&](int node, int parent) {
                std::vector<std::string> tips;
                bool has_child = false;
                for (int edge = m_hostDC.h_head[node]; edge != -1;
                     edge = m_hostDC.h_nxt[edge]) {
                    if (edge < 0 || edge >= edge_slots) {
                        throw std::runtime_error(
                            "D&C clade signature adjacency is out of range.");
                    }
                    const int neighbor = m_hostDC.h_e[edge];
                    if (neighbor == parent) continue;
                    has_child = true;
                    auto child_tips = dfs(neighbor, node);
                    tips.insert(tips.end(), child_tips.begin(), child_tips.end());
                }
                if (!has_child) {
                    if (node < 0 ||
                        node >= static_cast<int>(m_sortedNames.size()) ||
                        m_sortedNames[static_cast<size_t>(node)].empty()) {
                        throw std::runtime_error(
                            "D&C clade signature leaf name lookup failed.");
                    }
                    tips.push_back(m_sortedNames[static_cast<size_t>(node)]);
                }
                std::sort(tips.begin(), tips.end());
                const auto label_it = m_dcNodeToLabel.find(node);
                if (label_it == m_dcNodeToLabel.end()) {
                    throw std::runtime_error(
                        "D&C clade signature node has no stable label.");
                }
                std::string signature;
                for (size_t i = 0; i < tips.size(); ++i) {
                    if (i != 0) signature.push_back('\n');
                    signature += tips[i];
                }
                out.emplace_back(label_it->second, std::move(signature));
                return tips;
            };
        dfs(m_dcRootNode, -1);
        std::sort(out.begin(), out.end(),
                  [](const auto& lhs, const auto& rhs) {
                      return lhs.first < rhs.first;
                  });
        return out;
    }

    std::vector<std::tuple<int, int, std::string>>
    exportDivideAndConquerEdgeSplitSignatures() const
    {
        _requireOneOf(
            {State::GPUReady, State::Done},
            "exportDivideAndConquerEdgeSplitSignatures");
        _ensureDivideAndConquerPipeline(
            "exportDivideAndConquerEdgeSplitSignatures");
        _ensureDivideAndConquerSchedulerReady();
        const int edge_slots = static_cast<int>(m_numSequences) * 8;
        std::unordered_map<int, std::vector<std::string>> tips_by_node;
        std::vector<std::pair<int, int>> parent_edges;
        std::function<const std::vector<std::string>&(int, int)> dfs =
            [&](int node, int parent) -> const std::vector<std::string>& {
                std::vector<std::string> tips;
                bool has_child = false;
                for (int edge = m_hostDC.h_head[node]; edge != -1;
                     edge = m_hostDC.h_nxt[edge]) {
                    if (edge < 0 || edge >= edge_slots) {
                        throw std::runtime_error(
                            "D&C edge signature adjacency is out of range.");
                    }
                    const int neighbor = m_hostDC.h_e[edge];
                    if (neighbor == parent) continue;
                    has_child = true;
                    parent_edges.emplace_back(node, neighbor);
                    const auto& child_tips = dfs(neighbor, node);
                    tips.insert(tips.end(), child_tips.begin(), child_tips.end());
                }
                if (!has_child) {
                    tips.push_back(m_sortedNames.at(static_cast<size_t>(node)));
                }
                std::sort(tips.begin(), tips.end());
                return tips_by_node.emplace(node, std::move(tips)).first->second;
            };
        const auto& all_tips = dfs(m_dcRootNode, -1);
        std::vector<std::tuple<int, int, std::string>> out;
        out.reserve(parent_edges.size());
        for (const auto& edge : parent_edges) {
            const auto& one_side = tips_by_node.at(edge.second);
            std::vector<std::string> other;
            std::set_difference(
                all_tips.begin(), all_tips.end(),
                one_side.begin(), one_side.end(),
                std::back_inserter(other));
            const std::vector<std::string>* canonical = &one_side;
            if (other.size() < one_side.size() ||
                (other.size() == one_side.size() && other < one_side)) {
                canonical = &other;
            }
            std::string sig;
            for (size_t i = 0; i < canonical->size(); ++i) {
                if (i != 0) sig.push_back('\n');
                sig += (*canonical)[i];
            }
            out.emplace_back(
                m_dcNodeToLabel.at(edge.first),
                m_dcNodeToLabel.at(edge.second),
                std::move(sig));
        }
        return out;
    }

    /** Export the current GPU tree as a Newick string. */
    std::string exportTreeNewickString(bool print_binary_newick = true)
    {
        _requireOneOf(
            {State::GPUReady, State::Done},
            "exportTreeNewickString");

        char temp_path[] = "/tmp/dipper_session_tree_XXXXXX";
        const int fd = mkstemp(temp_path);
        if (fd < 0) {
            throw std::runtime_error(
                "[DIPPER] exportTreeNewickString: mkstemp() failed.");
        }
        close(fd);

        const bool previous_binary_mode = MashPlacement::g_printBinaryNewick;
        try {
            MashPlacement::g_printBinaryNewick = print_binary_newick;
            writeTree(temp_path);
            MashPlacement::g_printBinaryNewick = previous_binary_mode;

            std::ifstream ifs(temp_path);
            if (!ifs) {
                throw std::runtime_error(
                    "[DIPPER] exportTreeNewickString: failed to reopen temporary tree file.");
            }

            std::string newick;
            std::getline(ifs, newick);
            std::remove(temp_path);
            if (newick.empty()) {
                throw std::runtime_error(
                    "[DIPPER] exportTreeNewickString: exported Newick tree is empty.");
            }
            return newick;
        } catch (...) {
            MashPlacement::g_printBinaryNewick = previous_binary_mode;
            std::remove(temp_path);
            throw;
        }
    }

    /** Export one clade-signature string per live node label in the current tree.
     *
     *  The signature is the sorted set of descendant tip names joined by '\n'.
     *  This gives integration code a topology-stable way to match the same live
     *  node across independently loaded tree representations. */
    std::vector<std::pair<int, std::string>> exportStableLabelCladeSignatures() const
    {
        _requireOneOf({State::GPUReady, State::Done}, "exportStableLabelCladeSignatures");
        _ensureIterativePipeline("exportStableLabelCladeSignatures");

        const int node_slots = static_cast<int>(m_numSequences) * 2;
        const int edge_slots = static_cast<int>(m_numSequences) * 8;
        if (node_slots <= 0 || edge_slots <= 0) {
            throw std::runtime_error(
                "[DIPPER] exportStableLabelCladeSignatures: session has no live nodes.");
        }

        std::vector<int> h_head(static_cast<size_t>(node_slots), -1);
        std::vector<int> h_e(static_cast<size_t>(edge_slots), -1);
        std::vector<int> h_nxt(static_cast<size_t>(edge_slots), -1);
        cudaMemcpy(
            h_head.data(),
            m_kpArrays.d_head,
            static_cast<size_t>(node_slots) * sizeof(int),
            cudaMemcpyDeviceToHost);
        cudaMemcpy(
            h_e.data(),
            m_kpArrays.d_e,
            static_cast<size_t>(edge_slots) * sizeof(int),
            cudaMemcpyDeviceToHost);
        cudaMemcpy(
            h_nxt.data(),
            m_kpArrays.d_nxt,
            static_cast<size_t>(edge_slots) * sizeof(int),
            cudaMemcpyDeviceToHost);

        const int root_idx =
            static_cast<int>(m_numSequences) + m_kpArrays.bd - 2;
        if (root_idx < 0 || root_idx >= node_slots) {
            throw std::runtime_error(
                "[DIPPER] exportStableLabelCladeSignatures: computed root index is out of range.");
        }

        std::unordered_map<int, std::vector<std::string>> tips_by_node;
        std::function<const std::vector<std::string>&(int, int)> dfs =
            [&](int node_idx, int parent_idx) -> const std::vector<std::string>& {
                auto cached = tips_by_node.find(node_idx);
                if (cached != tips_by_node.end()) {
                    return cached->second;
                }

                std::vector<std::string> tips;
                bool has_child = false;
                for (int edge_idx = h_head[static_cast<size_t>(node_idx)];
                     edge_idx != -1;
                     edge_idx = h_nxt[static_cast<size_t>(edge_idx)]) {
                    const int neighbor = h_e[static_cast<size_t>(edge_idx)];
                    if (neighbor == parent_idx) {
                        continue;
                    }
                    has_child = true;
                    const std::vector<std::string>& child_tips =
                        dfs(neighbor, node_idx);
                    tips.insert(
                        tips.end(),
                        child_tips.begin(),
                        child_tips.end());
                }

                if (!has_child) {
                    if (node_idx < 0 ||
                        node_idx >= static_cast<int>(m_sortedNames.size()) ||
                        m_sortedNames[static_cast<size_t>(node_idx)].empty()) {
                        throw std::runtime_error(
                            "[DIPPER] exportStableLabelCladeSignatures: leaf node name lookup failed.");
                    }
                    tips.push_back(m_sortedNames[static_cast<size_t>(node_idx)]);
                }

                std::sort(tips.begin(), tips.end());
                auto [it, inserted] =
                    tips_by_node.emplace(node_idx, std::move(tips));
                (void)inserted;
                return it->second;
            };

        dfs(root_idx, -1);

        std::vector<std::pair<int, std::string>> out;
        out.reserve(tips_by_node.size());
        for (const auto& entry : tips_by_node) {
            const int label = nodeLabel(entry.first);
            if (label < 0) {
                continue;
            }
            out.emplace_back(label, _joinCladeSignature(entry.second));
        }
        std::sort(
            out.begin(),
            out.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.first < rhs.first;
            });
        return out;
    }

private:
    // ---- Session state machine -----------------------------------------------
    enum class State { Empty, TreeLoaded, SeqsLoaded, GPUReady, Done };
    enum class ActivePipeline { None, Iterative, DivideAndConquer };
    State m_state = State::Empty;
    ActivePipeline m_activePipeline = ActivePipeline::None;

    void _requireOneOf(std::initializer_list<State> allowed, const char* caller) const
    {
        for (State s : allowed)
            if (m_state == s) return;
        throw std::logic_error(
            std::string("[DIPPER] ") + caller + "() called out of order "
            "(current state: " + std::to_string(static_cast<int>(m_state)) + ").");
    }

    // ---- Owned data ----------------------------------------------------------
    Tree*                        m_tree           = nullptr;  // null = no backbone
    std::string                  m_backboneNewick;            // stored for deferred rebuild
    std::vector<std::string>     m_rawSeqs;
    std::vector<std::string>     m_rawNames;
    std::unordered_map<std::string, size_t> m_rawNameToIndex;
    std::vector<std::string>     m_sortedNames;
    std::unordered_map<int, int> m_idMap;
    size_t m_numSequences = 0;
    size_t m_backboneSize = 0;   // 0 when no backbone

    // ---- Sequential node labeling (pre-order DFS, 0-based) ------------------
    // Backbone nodes are labeled 0..K-1 at init; each placement adds two more.
    std::unordered_map<int, int> m_dipperToLabel;  // DIPPER idx → sequential label
    std::unordered_map<int, int> m_labelToDipper;  // sequential label → DIPPER idx (reverse)
    int m_nextLabel = 0;                           // next label to assign

    // ---- Two-phase placement state ------------------------------------------
    bool            m_placementPending = false;
    PlacementResult m_pendingResult    {};          // filled by findNextPlacement()

    uint64_t** m_fourBitSeqs = nullptr;
    uint64_t** m_twoBitSeqs  = nullptr;
    uint64_t*  m_seqLengths  = nullptr;

    MashPlacement::Param* m_params = nullptr;

    // GPU device array objects (each session owns its own, not the CLI statics).
    MashPlacement::MSADeviceArrays        m_msaArrays    {};
    MashPlacement::KPlacementDeviceArrays m_kpArrays     {};
    // Unused for MSA input, but required by existing method signatures:
    MashPlacement::MashDeviceArrays       m_mashArrays   {};
    MashPlacement::MatrixReader           m_matrixReader {};
    MashPlacement::MSADeviceArraysDC      m_msaArraysDC  {};
    MashPlacement::MashDeviceArraysDC     m_mashArraysDC {};
    MashPlacement::KPlacementDeviceArraysDC m_kpArraysDC {};
    mutable MashPlacement::KPlacementDeviceArraysHostDC m_hostDC {};
    size_t                                m_dcBackboneSize = 0;
    bool                                  m_dcBackboneBuilt = false;
    bool                                  m_dcTreeBuilt = false;
    mutable std::vector<int>              m_queryClusterIds;
    mutable std::vector<int>              m_queryEffectiveClusterIds;
    mutable std::unordered_set<std::string> m_dcScheduledQueryNames;
    mutable int                           m_requestedQueryClusterTopK = 1;
    mutable int                           m_computedQueryClusterTopK = 0;
    mutable bool                          m_queryClusterIdsReady = false;
    mutable size_t                        m_queryClusterBackboneSize = 0;
    mutable bool                          m_dcSchedulerReady = false;
    mutable size_t                        m_dcProcessedQueryCount = 0;
    mutable std::vector<std::vector<int>> m_dcClusterQueries;
    mutable std::vector<int>              m_dcOrderedClusters;
    mutable std::vector<int>              m_dcClusterPlacementOffsets;
    mutable std::vector<std::pair<int, int>> m_dcClusterAnchorLabels;
    mutable size_t                        m_dcCurrentClusterCursor = 0;
    mutable int                           m_dcActiveClusterId = -1;
    mutable size_t                        m_dcActiveClusterOffset = 0;
    mutable int                           m_dcActiveEdgeCount = 0;
    mutable int                           m_dcActiveLeafCount = 0;
    mutable int                           m_dcActiveIdx = 0;
    mutable int                           m_dcActiveInsertLeafCount = 0;
    mutable std::vector<int>              m_dcActiveEdgeMask;
    mutable std::vector<int>              m_dcActiveLeafMask;
    mutable std::vector<int>              m_dcActiveLeafMap;
    mutable std::vector<int>              m_dcActiveId;
    mutable std::vector<int>              m_dcActiveFrom;
    mutable std::vector<double>           m_dcActiveDis;
    mutable std::vector<double>           m_dcActiveDist;
    mutable std::vector<uint64_t>         m_dcActiveClusterCompressed;
    mutable bool                          m_dcAllSequencesUploaded = false;
    mutable std::vector<unsigned char>    m_dcInsertedQueryLeafIds;
    mutable std::vector<std::tuple<int, double, double>> m_dcMinPos;
    mutable std::unordered_map<int, int>  m_dcNodeToLabel;
    mutable std::unordered_map<int, int>  m_dcLabelToNode;
    std::unordered_map<std::string, int>  m_dcBackboneTipStableLabels;
    mutable int                           m_dcNextLabel = 0;
    mutable int                           m_dcRootNode = -1;
    mutable bool                          m_dcPendingStepReady = false;
    mutable PendingQueryInfo              m_dcPendingQueryInfo {};
    mutable int                           m_dcPendingLeafId = -1;
    mutable int                           m_dcPendingEdgeId = -1;
    mutable double                        m_dcPendingFracLen = 0.0;
    mutable double                        m_dcPendingAddLen = 0.0;
    mutable int*                          m_dcDeviceLeafIds = nullptr;
    mutable int*                          m_dcDeviceLeafMap = nullptr;

    // ---- Helpers -------------------------------------------------------------

    /**
     * Backbone path: backbone sequences fill tree-index slots; novel sequences
     * are appended as queries; internal backbone nodes are appended beyond the
     * query slots so that sortedNames[idx] is always valid for any DIPPER node
     * index that can appear in PlacementResult::splitNodeA / splitNodeB.
     *
     * Final layout of sortedNames:
     *   [0 .. backboneSize-1]   backbone leaf names (at their tree-assigned idx)
     *   [backboneSize .. numSeqs-1]  query names (in input order)
     *   [numSeqs .. ]           internal backbone node names (e.g. "node_1")
     */
    void _buildIdMapWithBackbone()
    {
        m_sortedNames.assign(m_backboneSize, "");
        m_idMap.clear();
        m_queryClusterIds.clear();
        m_queryEffectiveClusterIds.clear();
        m_queryClusterIdsReady = false;
        m_computedQueryClusterTopK = 0;
        m_queryClusterBackboneSize = 0;

        for (int i = 0; i < (int)m_numSequences; ++i) {
            auto it = m_tree->allNodes.find(m_rawNames[i]);
            if (it == m_tree->allNodes.end()) {
                // Novel query: append after backbone slots.
                m_sortedNames.push_back(m_rawNames[i]);
                m_idMap[i] = (int)m_sortedNames.size() - 1;
            } else {
                // Backbone sequence: fill its tree-index slot.
                m_sortedNames[it->second->idx] = m_rawNames[i];
                m_idMap[i] = it->second->idx;
            }
        }

        // Internal backbone nodes have idx >= m_numSequences (because the Tree
        // was constructed with totalLeaves = m_numSequences, pushing internal
        // node indices above the sequence slots).  Extend sortedNames so that
        // sortedNames[node->idx] gives the internal node's Newick name, making
        // PlacementResult::splitNodeA/B always safe to index.
        for (auto& [name, node] : m_tree->allNodes) {
            if (!node->children.empty()) {  // internal node
                const int nidx = node->idx;
                if (nidx >= (int)m_sortedNames.size())
                    m_sortedNames.resize(nidx + 1, "");
                m_sortedNames[nidx] = name;   // e.g. "node_1", "node_2", "node_3"
            }
        }
    }

    /**
     * Pre-order DFS of the backbone tree to assign sequential labels 0, 1, 2, …
     * to every node (root first, then each subtree in child order).
     * After this call, m_nextLabel == backboneNodeCount().
     */
    void _assignPreorderLabels()
    {
        m_dipperToLabel.clear();
        m_labelToDipper.clear();
        m_nextLabel = 0;
        std::function<void(Node*)> dfs = [&](Node* node) {
            if (!node) return;
            const int label = m_nextLabel++;
            m_dipperToLabel[node->idx] = label;
            m_labelToDipper[label]     = node->idx;
            for (Node* child : node->children)
                dfs(child);
        };
        dfs(m_tree->root);
    }

    /**
     * No-backbone path: sequences retain their input order (0, 1, 2, …).
     * Indices 0 and 1 seed the initial tree; the rest are queries.
     * Seed nodes get sequential labels 0 and 1.
     */
    void _buildIdMapNoBackbone()
    {
        m_backboneSize = 0;
        m_sortedNames  = m_rawNames;
        m_idMap.clear();
        m_queryClusterIds.clear();
        m_queryEffectiveClusterIds.clear();
        m_queryClusterIdsReady = false;
        m_computedQueryClusterTopK = 0;
        m_queryClusterBackboneSize = 0;
        for (int i = 0; i < (int)m_numSequences; ++i)
            m_idMap[i] = i;

        // Seed nodes 0 and 1 get the first two labels.
        m_dipperToLabel.clear();
        m_labelToDipper.clear();
        m_dipperToLabel[0] = 0;  m_labelToDipper[0] = 0;
        m_dipperToLabel[1] = 1;  m_labelToDipper[1] = 1;
        m_nextLabel = 2;
    }

    /** 4-bit compress all sequences in parallel using TBB. */
    void _compressSequences(const MashPlacement::Param& params)
    {
        m_fourBitSeqs = new uint64_t*[m_numSequences];
        m_seqLengths  = new uint64_t [m_numSequences];

        const bool rangeModify =
            (params.range.first > 0 || params.range.second > -1);

        tbb::parallel_for(
            tbb::blocked_range<int>(0, (int)m_numSequences),
            [&](tbb::blocked_range<int> r) {
                for (int i = r.begin(); i < r.end(); ++i) {
                    int seqLen = (int)m_rawSeqs[i].size();
                    if (rangeModify) {
                        if (params.range.second > -1) seqLen = params.range.second + 1;
                        if (params.range.first  > 0)  seqLen -= params.range.first;
                    }
                    const uint64_t compSize = params.isProtein
                        ? (uint64_t)(seqLen + 7)  / 8
                        : (uint64_t)(seqLen + 15) / 16;
                    uint64_t* buf = new uint64_t[compSize];
                    fourBitCompressor(
                        m_rawSeqs[i], m_rawSeqs[i].size(), buf,
                        params.range.first, params.range.second, params.isProtein);
                    const int newId      = m_idMap[i];
                    m_seqLengths[newId]  = (uint64_t)seqLen;
                    m_fourBitSeqs[newId] = buf;
                }
            });
    }

    void _compressSequencesTwoBit(const MashPlacement::Param& params)
    {
        m_twoBitSeqs = new uint64_t*[m_numSequences];
        m_seqLengths = new uint64_t[m_numSequences];
        tbb::parallel_for(
            tbb::blocked_range<int>(0, static_cast<int>(m_numSequences)),
            [&](tbb::blocked_range<int> r) {
                for (int i = r.begin(); i < r.end(); ++i) {
                    const size_t seq_len = m_rawSeqs[i].size();
                    const size_t words = params.isProtein
                        ? (seq_len + 7) / 8
                        : (seq_len + 31) / 32;
                    uint64_t* buffer = new uint64_t[words];
                    twoBitCompressor(
                        m_rawSeqs[i], seq_len, buffer, params.isProtein,
                        params.useReducedProtein, params.useMurphy8);
                    const int new_id = m_idMap[i];
                    m_seqLengths[new_id] = seq_len;
                    m_twoBitSeqs[new_id] = buffer;
                }
            });
    }

    void _cleanup()
    {
        if (m_state >= State::GPUReady) {
            if (m_activePipeline == ActivePipeline::Iterative) {
                m_kpArrays.endIterativePlacement();
                m_msaArrays.deallocateDeviceArrays();
                m_kpArrays.deallocateDeviceArrays();
            } else if (m_activePipeline == ActivePipeline::DivideAndConquer) {
                if (m_params != nullptr && m_params->in == "r") {
                    m_mashArraysDC.deallocateDeviceArraysDC();
                } else {
                    m_msaArraysDC.deallocateDeviceArraysDC();
                }
                m_kpArraysDC.deallocateDeviceArraysDC();
                m_hostDC.deallocateHostArraysDC();
                if (m_hostDC.clusterID != nullptr) {
                    delete[] m_hostDC.clusterID;
                    m_hostDC.clusterID = nullptr;
                }
            }
        }
        if (m_dcDeviceLeafIds != nullptr) {
            cudaFree(m_dcDeviceLeafIds);
            m_dcDeviceLeafIds = nullptr;
        }
        if (m_dcDeviceLeafMap != nullptr) {
            cudaFree(m_dcDeviceLeafMap);
            m_dcDeviceLeafMap = nullptr;
        }
        if (m_fourBitSeqs) {
            for (size_t i = 0; i < m_numSequences; ++i)
                delete[] m_fourBitSeqs[i];
            delete[] m_fourBitSeqs;
            m_fourBitSeqs = nullptr;
        }
        if (m_twoBitSeqs) {
            for (size_t i = 0; i < m_numSequences; ++i)
                delete[] m_twoBitSeqs[i];
            delete[] m_twoBitSeqs;
            m_twoBitSeqs = nullptr;
        }
        delete[] m_seqLengths;
        m_seqLengths = nullptr;
        delete m_tree;
        m_tree = nullptr;
        m_queryClusterIds.clear();
        m_queryEffectiveClusterIds.clear();
        m_queryClusterIdsReady = false;
        m_computedQueryClusterTopK = 0;
        m_queryClusterBackboneSize = 0;
        m_dcBackboneSize = 0;
        m_dcBackboneBuilt = false;
        m_dcTreeBuilt = false;
        _resetDivideAndConquerSchedulerState();
        m_activePipeline = ActivePipeline::None;
    }

    void _ensureNoActivePipeline(const char* caller) const
    {
        if (m_activePipeline == ActivePipeline::None) {
            return;
        }
        throw std::logic_error(
            std::string("[DIPPER] ") + caller +
            ": a GPU pipeline is already active in this session.");
    }

    void _ensureIterativePipeline(const char* caller) const
    {
        if (m_activePipeline == ActivePipeline::Iterative) {
            return;
        }
        throw std::logic_error(
            std::string("[DIPPER] ") + caller +
            ": this API is only available for the standard iterative placement pipeline.");
    }

    void _ensureDivideAndConquerPipeline(const char* caller) const
    {
        if (m_activePipeline == ActivePipeline::DivideAndConquer) {
            return;
        }
        throw std::logic_error(
            std::string("[DIPPER] ") + caller +
            ": this API is only available for the divide-and-conquer pipeline.");
    }

    size_t _currentBackboneSize() const
    {
        if (m_activePipeline == ActivePipeline::DivideAndConquer &&
            m_dcBackboneSize > 0) {
            return m_dcBackboneSize;
        }
        return m_backboneSize;
    }

    size_t _resolveDivideAndConquerBackboneSize(size_t requested_backbone_size) const
    {
        if (m_tree != nullptr && m_backboneSize > 0) {
            if (requested_backbone_size != 0 &&
                requested_backbone_size != m_backboneSize) {
                throw std::invalid_argument(
                    "[DIPPER] initializeDivideAndConquerGPU: when a backbone tree is loaded, backboneSize must be 0 or equal to the loaded backbone size.");
            }
            return m_backboneSize;
        }

        if (requested_backbone_size != 0) {
            return requested_backbone_size;
        }
        return std::max<size_t>(2, m_numSequences / 20);
    }

    std::vector<std::string> _leafNamesForCurrentSession() const
    {
        if (m_sortedNames.size() < m_numSequences) {
            throw std::runtime_error(
                "[DIPPER] _leafNamesForCurrentSession: sortedNames() is smaller than numSequences().");
        }
        return std::vector<std::string>(
            m_sortedNames.begin(),
            m_sortedNames.begin() + static_cast<long>(m_numSequences));
    }

    static void _checkCudaStatus(cudaError_t status, const char* context)
    {
        if (status == cudaSuccess) {
            return;
        }
        throw std::runtime_error(
            std::string("[DIPPER] ") + context +
            ": CUDA failure: " + cudaGetErrorString(status));
    }

    void _resetDivideAndConquerSchedulerState()
    {
        m_dcSchedulerReady = false;
        m_dcProcessedQueryCount = 0;
        m_dcClusterQueries.clear();
        m_dcOrderedClusters.clear();
        m_dcClusterPlacementOffsets.clear();
        m_dcClusterAnchorLabels.clear();
        m_dcCurrentClusterCursor = 0;
        m_dcActiveClusterId = -1;
        m_dcActiveClusterOffset = 0;
        m_dcActiveEdgeCount = 0;
        m_dcActiveLeafCount = 0;
        m_dcActiveIdx = 0;
        m_dcActiveInsertLeafCount = 0;
        m_dcActiveEdgeMask.clear();
        m_dcActiveLeafMask.clear();
        m_dcActiveLeafMap.clear();
        m_dcActiveId.clear();
        m_dcActiveFrom.clear();
        m_dcActiveDis.clear();
        m_dcActiveDist.clear();
        m_dcActiveClusterCompressed.clear();
        m_dcAllSequencesUploaded = false;
        m_dcInsertedQueryLeafIds.clear();
        m_dcMinPos.clear();
        m_dcNodeToLabel.clear();
        m_dcLabelToNode.clear();
        m_dcNextLabel = 0;
        m_dcRootNode = -1;
        m_dcPendingStepReady = false;
        m_dcPendingQueryInfo = {};
        m_dcPendingLeafId = -1;
        m_dcPendingEdgeId = -1;
        m_dcPendingFracLen = 0.0;
        m_dcPendingAddLen = 0.0;
    }

    void _ensureDivideAndConquerHostArraysAllocated() const
    {
        if (m_hostDC.h_head != nullptr) {
            return;
        }
        const_cast<MashPlacement::KPlacementDeviceArraysHostDC&>(m_hostDC)
            .allocateHostArraysDC(m_dcBackboneSize, m_numSequences);
    }

    void _copyDivideAndConquerTreeStateToHost() const
    {
        _ensureDivideAndConquerHostArraysAllocated();
        const size_t node_slots = m_numSequences * 2;
        const size_t edge_slots = m_numSequences * 8;
        const size_t closest_slots = edge_slots * 5;
        _checkCudaStatus(
            cudaMemcpy(
                m_hostDC.h_head,
                m_kpArraysDC.d_head,
                node_slots * sizeof(int),
                cudaMemcpyDeviceToHost),
            "_copyDivideAndConquerTreeStateToHost(head)");
        _checkCudaStatus(
            cudaMemcpy(
                m_hostDC.h_e,
                m_kpArraysDC.d_e,
                edge_slots * sizeof(int),
                cudaMemcpyDeviceToHost),
            "_copyDivideAndConquerTreeStateToHost(e)");
        _checkCudaStatus(
            cudaMemcpy(
                m_hostDC.h_nxt,
                m_kpArraysDC.d_nxt,
                edge_slots * sizeof(int),
                cudaMemcpyDeviceToHost),
            "_copyDivideAndConquerTreeStateToHost(nxt)");
        _checkCudaStatus(
            cudaMemcpy(
                m_hostDC.h_belong,
                m_kpArraysDC.d_belong,
                edge_slots * sizeof(int),
                cudaMemcpyDeviceToHost),
            "_copyDivideAndConquerTreeStateToHost(belong)");
        _checkCudaStatus(
            cudaMemcpy(
                m_hostDC.h_len,
                m_kpArraysDC.d_len,
                edge_slots * sizeof(double),
                cudaMemcpyDeviceToHost),
            "_copyDivideAndConquerTreeStateToHost(len)");
        _checkCudaStatus(
            cudaMemcpy(
                m_hostDC.h_closest_dis,
                m_kpArraysDC.d_closest_dis,
                closest_slots * sizeof(double),
                cudaMemcpyDeviceToHost),
            "_copyDivideAndConquerTreeStateToHost(closest_dis)");
        _checkCudaStatus(
            cudaMemcpy(
                m_hostDC.h_closest_id,
                m_kpArraysDC.d_closest_id,
                closest_slots * sizeof(int),
                cudaMemcpyDeviceToHost),
            "_copyDivideAndConquerTreeStateToHost(closest_id)");
    }

    int _assignDivideAndConquerNodeLabel(int node_idx) const
    {
        const auto existing = m_dcNodeToLabel.find(node_idx);
        if (existing != m_dcNodeToLabel.end()) {
            return existing->second;
        }
        int label = -1;
        if (node_idx >= 0 &&
            node_idx < static_cast<int>(m_sortedNames.size())) {
            const std::string& name =
                m_sortedNames[static_cast<size_t>(node_idx)];
            int degree = 0;
            for (int edge_idx = m_hostDC.h_head[node_idx];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                if (m_hostDC.h_e[edge_idx] >= 0) {
                    ++degree;
                }
            }
            const bool is_shadow_leaf = degree <= 1;
            const auto tip_it = m_dcBackboneTipStableLabels.find(name);
            if (is_shadow_leaf &&
                tip_it != m_dcBackboneTipStableLabels.end()) {
                label = tip_it->second;
            } else if (name.rfind("__mlipper_id_", 0) == 0) {
                label = std::stoi(name.substr(13));
            }
        }
        if (label < 0) {
            const auto label_reserved_for_tip = [&](int candidate) {
                return std::any_of(
                    m_dcBackboneTipStableLabels.begin(),
                    m_dcBackboneTipStableLabels.end(),
                    [&](const auto& entry) {
                        return entry.second == candidate;
                    });
            };
            while (m_dcLabelToNode.count(m_dcNextLabel) != 0 ||
                   label_reserved_for_tip(m_dcNextLabel)) {
                ++m_dcNextLabel;
            }
            label = m_dcNextLabel;
        }
        if (m_dcLabelToNode.count(label) != 0) {
            throw std::runtime_error(
                "[DIPPER] duplicate stable node label in backbone tree.");
        }
        m_dcNextLabel = std::max(m_dcNextLabel, label + 1);
        m_dcNodeToLabel[node_idx] = label;
        m_dcLabelToNode[label] = node_idx;
        return label;
    }

    void _assignDivideAndConquerShadowLabels() const
    {
        m_dcNodeToLabel.clear();
        m_dcLabelToNode.clear();
        m_dcNextLabel = 0;
        const int root_idx =
            (m_dcRootNode >= 0)
                ? m_dcRootNode
                : static_cast<int>(m_numSequences);
        if (root_idx < 0 ||
            root_idx >= static_cast<int>(m_numSequences * 2) ||
            m_hostDC.h_head[root_idx] == -1) {
            throw std::runtime_error(
                "[DIPPER] _assignDivideAndConquerShadowLabels: divide-and-conquer root is unavailable.");
        }
        std::function<void(int, int)> dfs = [&](int node_idx, int parent_idx) {
            _assignDivideAndConquerNodeLabel(node_idx);
            for (int edge_idx = m_hostDC.h_head[node_idx];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                const int neighbor = m_hostDC.h_e[edge_idx];
                if (neighbor == parent_idx || neighbor < 0) {
                    continue;
                }
                if (m_dcNodeToLabel.find(neighbor) != m_dcNodeToLabel.end()) {
                    continue;
                }
                dfs(neighbor, node_idx);
            }
        };
        dfs(root_idx, -1);
    }

    void _ensureDivideAndConquerDeviceScratch() const
    {
        const size_t leaf_slots = std::max<size_t>(m_numSequences * 2, 1);
        if (m_dcDeviceLeafIds == nullptr) {
            _checkCudaStatus(
                cudaMalloc(&m_dcDeviceLeafIds, leaf_slots * sizeof(int)),
                "_ensureDivideAndConquerDeviceScratch(leafIds)");
        }
        if (m_dcDeviceLeafMap == nullptr) {
            _checkCudaStatus(
                cudaMalloc(&m_dcDeviceLeafMap, leaf_slots * sizeof(int)),
                "_ensureDivideAndConquerDeviceScratch(leafMap)");
        }
    }

    void _buildDivideAndConquerClusterSchedule() const
    {
        const int cluster_slots =
            std::max(0, m_hostDC.numSequences * 4 - 4);
        m_dcClusterQueries.assign(
            static_cast<size_t>(cluster_slots),
            std::vector<int>{});
        m_dcClusterAnchorLabels.assign(
            static_cast<size_t>(cluster_slots),
            std::pair<int, int>{-1, -1});
        for (size_t idx = m_dcBackboneSize; idx < m_numSequences; ++idx) {
            if (!m_dcScheduledQueryNames.empty() &&
                m_dcScheduledQueryNames.count(m_sortedNames[idx]) == 0) {
                continue;
            }
            const std::vector<int>& scheduled_cluster_ids =
                m_queryEffectiveClusterIds.empty()
                    ? m_queryClusterIds
                    : m_queryEffectiveClusterIds;
            const int cluster_id = scheduled_cluster_ids[idx];
            if (cluster_id < 0 || cluster_id >= cluster_slots) {
                throw std::runtime_error(
                    "[DIPPER] _buildDivideAndConquerClusterSchedule: cluster id is out of range.");
            }
            m_dcClusterQueries[static_cast<size_t>(cluster_id)].push_back(
                static_cast<int>(idx));
        }

        m_dcOrderedClusters.clear();
        m_dcClusterPlacementOffsets.clear();
        int prefix = 0;
        for (int cluster_id = 0; cluster_id < cluster_slots; ++cluster_id) {
            const std::vector<int>& queries =
                m_dcClusterQueries[static_cast<size_t>(cluster_id)];
            if (m_hostDC.h_belong[cluster_id] >= 0 &&
                m_hostDC.h_e[cluster_id] >= 0) {
                m_dcClusterAnchorLabels[static_cast<size_t>(cluster_id)] =
                    std::pair<int, int>{
                        nodeLabel(m_hostDC.h_belong[cluster_id]),
                        nodeLabel(m_hostDC.h_e[cluster_id])};
            }
            if (queries.empty()) {
                continue;
            }
            m_dcOrderedClusters.push_back(cluster_id);
            m_dcClusterPlacementOffsets.push_back(prefix);
            prefix += static_cast<int>(queries.size());
        }
    }

    void _uploadDivideAndConquerActiveClusterSequences() const
    {
        if (m_activePipeline != ActivePipeline::DivideAndConquer) {
            return;
        }
        if (m_dcAllSequencesUploaded) {
            return;
        }
        if (m_params == nullptr || m_params->in != "m") {
            throw std::runtime_error(
                "[DIPPER] _uploadDivideAndConquerActiveClusterSequences: only aligned-sequence D&C scheduling is currently supported.");
        }
        if (m_msaArraysDC.h_compressedSeqs == nullptr) {
            throw std::runtime_error(
                "[DIPPER] _uploadDivideAndConquerActiveClusterSequences: missing compressed sequence buffer.");
        }

        const size_t max_length_compressed =
            static_cast<size_t>(m_msaArraysDC.d_seqLen + 15) / 16;
        m_dcActiveClusterCompressed.assign(
            m_numSequences * max_length_compressed,
            0);
        std::copy_n(
            m_msaArraysDC.h_compressedSeqs,
            m_dcActiveClusterCompressed.size(),
            m_dcActiveClusterCompressed.data());
        _checkCudaStatus(
            cudaMemcpy(
                m_msaArraysDC.d_compressedSeqsConst,
                m_dcActiveClusterCompressed.data(),
                m_dcActiveClusterCompressed.size() * sizeof(uint64_t),
                cudaMemcpyHostToDevice),
            "_uploadDivideAndConquerActiveClusterSequences");
        m_dcAllSequencesUploaded = true;
    }

    void _requireDivideAndConquerHostArraysReady(const char* context) const
    {
        if (m_hostDC.h_head == nullptr || m_hostDC.h_e == nullptr ||
            m_hostDC.h_nxt == nullptr || m_hostDC.h_belong == nullptr ||
            m_hostDC.h_len == nullptr || m_hostDC.h_closest_dis == nullptr ||
            m_hostDC.h_closest_id == nullptr) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": divide-and-conquer host arrays are not allocated.");
        }
    }

    int _divideAndConquerNodeSlots() const
    {
        return static_cast<int>(m_numSequences * 2);
    }

    int _divideAndConquerEdgeSlots() const
    {
        return static_cast<int>(m_numSequences * 8);
    }

    int _divideAndConquerClosestEdgeSlots() const
    {
        return static_cast<int>(m_numSequences * 8);
    }

    void _requireDivideAndConquerNodeIndex(
        int node_idx,
        const char* context) const
    {
        const int node_slots = _divideAndConquerNodeSlots();
        if (node_idx < 0 || node_idx >= node_slots) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": node index " + std::to_string(node_idx) +
                " is outside divide-and-conquer node capacity " +
                std::to_string(node_slots) + ".");
        }
    }

    void _requireDivideAndConquerClosestEdgeIndex(
        int edge_idx,
        const char* context) const
    {
        const int closest_edge_slots =
            _divideAndConquerClosestEdgeSlots();
        if (edge_idx < 0 || edge_idx >= closest_edge_slots) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": edge index " + std::to_string(edge_idx) +
                " is outside divide-and-conquer closest-array capacity " +
                std::to_string(closest_edge_slots) + ".");
        }
    }

    void _requireDivideAndConquerPlacementUpdateCapacity(
        int edge_id,
        int pending_leaf_id,
        int new_internal_idx,
        const char* context) const
    {
        _requireDivideAndConquerHostArraysReady(context);
        _requireDivideAndConquerClosestEdgeIndex(edge_id, context);
        _requireDivideAndConquerNodeIndex(pending_leaf_id, context);
        _requireDivideAndConquerNodeIndex(new_internal_idx, context);
        const int closest_edge_slots =
            _divideAndConquerClosestEdgeSlots();
        if (m_dcActiveIdx < 0 || m_dcActiveIdx + 3 >= closest_edge_slots) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": active edge insertion index would exceed closest-array capacity.");
        }
    }

    void _requireDivideAndConquerActiveClusterCapacity(
        int added_edges,
        int added_leaves,
        const char* context) const
    {
        if (added_edges < 0 || added_leaves < 0) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": active-cluster capacity increments must be non-negative.");
        }
        if (m_dcActiveEdgeCount < 0 || m_dcActiveLeafCount < 0) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": active-cluster counts are negative.");
        }
        if (m_dcActiveEdgeCount + added_edges >
                static_cast<int>(m_dcActiveEdgeMask.size()) ||
            m_dcActiveLeafCount + added_leaves >
                static_cast<int>(m_dcActiveLeafMask.size()) ||
            m_dcActiveLeafCount + added_leaves >
                static_cast<int>(m_dcActiveLeafMap.size())) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": active-cluster edge/leaf mask capacity would be exceeded.");
        }
        const int node_slots = _divideAndConquerNodeSlots();
        if (static_cast<int>(m_dcActiveId.size()) < node_slots ||
            static_cast<int>(m_dcActiveFrom.size()) < node_slots ||
            static_cast<int>(m_dcActiveDis.size()) < node_slots) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": active-cluster traversal buffers are smaller than node capacity.");
        }
    }

    void _requireDivideAndConquerClosestLeafIds(
        int edge_idx,
        const char* context) const
    {
        _requireDivideAndConquerClosestEdgeIndex(edge_idx, context);
        for (int rank = 0; rank < 5; ++rank) {
            const int leaf_id = m_hostDC.h_closest_id[edge_idx * 5 + rank];
            if (leaf_id == -1) {
                continue;
            }
            if (leaf_id < 0 || leaf_id >= static_cast<int>(m_numSequences)) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": closest leaf id " + std::to_string(leaf_id) +
                    " is outside loaded sequence capacity.");
            }
        }
    }

    void _requireDivideAndConquerClusterAnchorReady(
        int anchor_edge_id,
        const char* context) const
    {
        _requireDivideAndConquerHostArraysReady(context);
        _requireDivideAndConquerClosestEdgeIndex(anchor_edge_id, context);
        const int from_node = m_hostDC.h_belong[anchor_edge_id];
        const int to_node = m_hostDC.h_e[anchor_edge_id];
        _requireDivideAndConquerNodeIndex(from_node, context);
        _requireDivideAndConquerNodeIndex(to_node, context);
        const int reverse_edge =
            _findDivideAndConquerDirectedEdge(to_node, from_node);
        if (reverse_edge < 0) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": cluster anchor edge is missing its reverse direction.");
        }
        _requireDivideAndConquerClosestLeafIds(anchor_edge_id, context);
        _requireDivideAndConquerClosestLeafIds(reverse_edge, context);
    }

    void _initializeDivideAndConquerClusterMasks(
        int anchor_edge_id,
        std::vector<int>& edge_mask,
        std::vector<int>& leaf_mask,
        const char* context) const
    {
        if (edge_mask.size() < 2 || leaf_mask.size() < 10) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": active cluster masks are smaller than the legacy initializer requires.");
        }
        _requireDivideAndConquerClusterAnchorReady(anchor_edge_id, context);
        initializeClusterCpu(
            anchor_edge_id,
            m_hostDC.h_e,
            m_hostDC.h_belong,
            m_hostDC.h_head,
            m_hostDC.h_nxt,
            m_hostDC.h_closest_id,
            edge_mask,
            leaf_mask);
        for (int idx = 0; idx < 2; ++idx) {
            _requireDivideAndConquerClosestEdgeIndex(
                edge_mask[static_cast<size_t>(idx)],
                context);
        }
        for (int idx = 0; idx < 10; ++idx) {
            const int leaf_id = leaf_mask[static_cast<size_t>(idx)];
            if (leaf_id == -1) {
                continue;
            }
            if (leaf_id < 0 || leaf_id >= static_cast<int>(m_numSequences)) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": initialized active leaf id is outside loaded sequence capacity.");
            }
        }
    }

    void _requireDivideAndConquerActiveLeafMaskReady(
        const char* context) const
    {
        if (m_dcActiveLeafCount < 0 ||
            m_dcActiveLeafCount >
                static_cast<int>(m_dcActiveLeafMask.size()) ||
            m_dcActiveLeafCount >
                static_cast<int>(m_dcActiveLeafMap.size())) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": active leaf-mask count exceeds buffer capacity.");
        }
        if (m_dcActiveClusterId < 0 ||
            m_dcActiveClusterId >=
                static_cast<int>(m_dcClusterQueries.size())) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": active cluster id is out of range.");
        }
        for (int idx = 0; idx < m_dcActiveLeafCount; ++idx) {
            const int leaf_id =
                m_dcActiveLeafMask[static_cast<size_t>(idx)];
            const int leaf_map =
                m_dcActiveLeafMap[static_cast<size_t>(idx)];
            if (leaf_id == -1) {
                continue;
            }
            if (leaf_id < 0 || leaf_id >= static_cast<int>(m_numSequences)) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": active leaf id is outside loaded sequence capacity.");
            }
            if (leaf_id < static_cast<int>(m_dcBackboneSize)) {
                if (leaf_map != leaf_id) {
                    throw std::runtime_error(
                        std::string("[DIPPER] ") + context +
                        ": active backbone leaf does not map to its global sequence id.");
                }
                continue;
            }
            if (leaf_map != leaf_id) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": active query leaf does not map to its global sequence id.");
            }
            if (m_dcInsertedQueryLeafIds.empty() ||
                leaf_id >=
                    static_cast<int>(m_dcInsertedQueryLeafIds.size()) ||
                !m_dcInsertedQueryLeafIds[static_cast<size_t>(leaf_id)]) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": active query leaf has not been inserted into the shadow tree.");
            }
        }
    }

    void _refreshDivideAndConquerSeedLeafMaps(
        int seed_leaf_count,
        const char* context) const
    {
        if (seed_leaf_count < 0 ||
            seed_leaf_count > static_cast<int>(m_dcActiveLeafMask.size()) ||
            seed_leaf_count > static_cast<int>(m_dcActiveLeafMap.size())) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": seed leaf count exceeds active mask capacity.");
        }
        if (m_dcActiveClusterId < 0 ||
            m_dcActiveClusterId >=
                static_cast<int>(m_dcClusterQueries.size())) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": active cluster id is out of range.");
        }
        for (int idx = 0; idx < seed_leaf_count; ++idx) {
            const int leaf_id =
                m_dcActiveLeafMask[static_cast<size_t>(idx)];
            if (leaf_id == -1) {
                m_dcActiveLeafMap[static_cast<size_t>(idx)] = -1;
                continue;
            }
            if (leaf_id < 0 || leaf_id >= static_cast<int>(m_numSequences)) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": seed leaf id is outside loaded sequence capacity.");
            }
            if (leaf_id < static_cast<int>(m_dcBackboneSize)) {
                m_dcActiveLeafMap[static_cast<size_t>(idx)] = leaf_id;
                continue;
            }
            if (m_dcInsertedQueryLeafIds.empty() ||
                leaf_id >=
                    static_cast<int>(m_dcInsertedQueryLeafIds.size()) ||
                !m_dcInsertedQueryLeafIds[static_cast<size_t>(leaf_id)]) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": seed closest leaf refers to a query that has not been inserted.");
            }
            m_dcActiveLeafMap[static_cast<size_t>(idx)] = leaf_id;
        }
    }

    void _requireDivideAndConquerActiveEdgeMaskReady(
        int edge_count,
        const char* context) const
    {
        if (edge_count <= 0 ||
            edge_count > static_cast<int>(m_dcActiveEdgeMask.size()) ||
            edge_count > static_cast<int>(m_dcMinPos.size())) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": active edge-mask count is outside buffer capacity.");
        }
        for (int idx = 0; idx < edge_count; ++idx) {
            const int edge_id =
                m_dcActiveEdgeMask[static_cast<size_t>(idx)];
            _requireDivideAndConquerClosestEdgeIndex(edge_id, context);
            const int from_node = m_hostDC.h_belong[edge_id];
            const int to_node = m_hostDC.h_e[edge_id];
            _requireDivideAndConquerNodeIndex(from_node, context);
            _requireDivideAndConquerNodeIndex(to_node, context);
            const int reverse_edge =
                _findDivideAndConquerDirectedEdge(to_node, from_node);
            if (reverse_edge < 0) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": active edge is missing its reverse direction.");
            }
            _requireDivideAndConquerClosestLeafIds(edge_id, context);
            _requireDivideAndConquerClosestLeafIds(reverse_edge, context);
        }
    }

    int _findDivideAndConquerDirectedEdge(int from_node, int to_node) const
    {
        _requireDivideAndConquerHostArraysReady(
            "_findDivideAndConquerDirectedEdge");
        if (from_node < 0 || to_node < 0) {
            return -1;
        }
        const int node_slots = _divideAndConquerNodeSlots();
        if (from_node >= node_slots || to_node >= node_slots) {
            return -1;
        }
        const int edge_slots = _divideAndConquerEdgeSlots();
        int chain_steps = 0;
        for (int edge_idx = m_hostDC.h_head[from_node];
             edge_idx != -1;
             edge_idx = m_hostDC.h_nxt[edge_idx]) {
            if (edge_idx < 0 || edge_idx >= edge_slots) {
                throw std::runtime_error(
                    "[DIPPER] _findDivideAndConquerDirectedEdge: adjacency edge index is out of range.");
            }
            if (++chain_steps > edge_slots) {
                throw std::runtime_error(
                    "[DIPPER] _findDivideAndConquerDirectedEdge: adjacency list contains a cycle.");
            }
            if (m_hostDC.h_e[edge_idx] == to_node) {
                return edge_idx;
            }
        }
        return -1;
    }

    int _unlinkDivideAndConquerDirectedEdge(int from_node, int to_node) const
    {
        _requireDivideAndConquerHostArraysReady(
            "_unlinkDivideAndConquerDirectedEdge");
        if (from_node < 0 || to_node < 0) {
            return -1;
        }
        const int node_slots = _divideAndConquerNodeSlots();
        if (from_node >= node_slots || to_node >= node_slots) {
            return -1;
        }
        const int edge_slots = _divideAndConquerEdgeSlots();
        int prev_edge_idx = -1;
        int chain_steps = 0;
        for (int edge_idx = m_hostDC.h_head[from_node];
             edge_idx != -1;
             edge_idx = m_hostDC.h_nxt[edge_idx]) {
            if (edge_idx < 0 || edge_idx >= edge_slots) {
                throw std::runtime_error(
                    "[DIPPER] _unlinkDivideAndConquerDirectedEdge: adjacency edge index is out of range.");
            }
            if (++chain_steps > edge_slots) {
                throw std::runtime_error(
                    "[DIPPER] _unlinkDivideAndConquerDirectedEdge: adjacency list contains a cycle.");
            }
            if (m_hostDC.h_e[edge_idx] != to_node) {
                prev_edge_idx = edge_idx;
                continue;
            }
            const int next_edge_idx = m_hostDC.h_nxt[edge_idx];
            if (prev_edge_idx < 0) {
                m_hostDC.h_head[from_node] = next_edge_idx;
            } else {
                m_hostDC.h_nxt[prev_edge_idx] = next_edge_idx;
            }
            m_hostDC.h_nxt[edge_idx] = -1;
            return edge_idx;
        }
        return -1;
    }

    void _attachDivideAndConquerDirectedEdge(
        int edge_idx,
        int from_node,
        int to_node,
        double edge_length) const
    {
        _requireDivideAndConquerHostArraysReady(
            "_attachDivideAndConquerDirectedEdge");
        _requireDivideAndConquerClosestEdgeIndex(
            edge_idx,
            "_attachDivideAndConquerDirectedEdge");
        _requireDivideAndConquerNodeIndex(
            from_node,
            "_attachDivideAndConquerDirectedEdge(from)");
        _requireDivideAndConquerNodeIndex(
            to_node,
            "_attachDivideAndConquerDirectedEdge(to)");
        m_hostDC.h_e[edge_idx] = to_node;
        m_hostDC.h_len[edge_idx] = edge_length;
        m_hostDC.h_belong[edge_idx] = from_node;
        m_hostDC.h_nxt[edge_idx] = m_hostDC.h_head[from_node];
        m_hostDC.h_head[from_node] = edge_idx;
    }

    void _validateDivideAndConquerShadowTree(const char* context) const
    {
        const int node_slots = static_cast<int>(m_numSequences * 2);
        const int edge_slots = static_cast<int>(m_numSequences * 8);
        const int closest_edge_slots = static_cast<int>(m_numSequences * 8);
        if (node_slots <= 0 || edge_slots <= 0) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": shadow-tree capacity is empty.");
        }
        if (m_dcRootNode < 0 || m_dcRootNode >= node_slots) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": shadow-tree root is out of range.");
        }

        std::vector<int> edge_owner(static_cast<size_t>(edge_slots), -1);
        std::vector<int> degree(static_cast<size_t>(node_slots), 0);
        int active_directed_edges = 0;
        for (int node = 0; node < node_slots; ++node) {
            int chain_steps = 0;
            for (int edge_idx = m_hostDC.h_head[node];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                if (edge_idx < 0 || edge_idx >= edge_slots) {
                    throw std::runtime_error(
                        std::string("[DIPPER] ") + context +
                        ": adjacency edge index is out of range for node " +
                        std::to_string(node) + ".");
                }
                if (++chain_steps > edge_slots) {
                    throw std::runtime_error(
                        std::string("[DIPPER] ") + context +
                        ": adjacency list contains a cycle for node " +
                        std::to_string(node) + ".");
                }
                if (edge_owner[static_cast<size_t>(edge_idx)] >= 0) {
                    throw std::runtime_error(
                        std::string("[DIPPER] ") + context +
                        ": edge slot " + std::to_string(edge_idx) +
                        " is linked from multiple adjacency lists.");
                }
                if (m_hostDC.h_belong[edge_idx] != node) {
                    throw std::runtime_error(
                        std::string("[DIPPER] ") + context +
                        ": edge slot " + std::to_string(edge_idx) +
                        " has mismatched belong/head ownership.");
                }
                if (edge_idx >= closest_edge_slots) {
                    throw std::runtime_error(
                        std::string("[DIPPER] ") + context +
                        ": edge slot " + std::to_string(edge_idx) +
                        " exceeds the divide-and-conquer closest-array capacity.");
                }
                const int neighbor = m_hostDC.h_e[edge_idx];
                if (neighbor < 0 || neighbor >= node_slots) {
                    throw std::runtime_error(
                        std::string("[DIPPER] ") + context +
                        ": edge slot " + std::to_string(edge_idx) +
                        " points to an out-of-range node.");
                }
                const int next_edge = m_hostDC.h_nxt[edge_idx];
                if (next_edge < -1 || next_edge >= edge_slots) {
                    throw std::runtime_error(
                        std::string("[DIPPER] ") + context +
                        ": edge slot " + std::to_string(edge_idx) +
                        " has an out-of-range nxt pointer.");
                }
                edge_owner[static_cast<size_t>(edge_idx)] = node;
                ++degree[static_cast<size_t>(node)];
                ++active_directed_edges;
            }
        }

        if (active_directed_edges == 0) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": shadow tree has no active edges.");
        }

        std::vector<char> visited(static_cast<size_t>(node_slots), 0);
        std::vector<int> queue;
        queue.reserve(static_cast<size_t>(node_slots));
        queue.push_back(m_dcRootNode);
        visited[static_cast<size_t>(m_dcRootNode)] = 1;
        for (size_t queue_idx = 0; queue_idx < queue.size(); ++queue_idx) {
            const int node = queue[queue_idx];
            for (int edge_idx = m_hostDC.h_head[node];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                const int neighbor = m_hostDC.h_e[edge_idx];
                const int reverse_edge =
                    _findDivideAndConquerDirectedEdge(neighbor, node);
                if (reverse_edge < 0) {
                    throw std::runtime_error(
                        std::string("[DIPPER] ") + context +
                        ": edge " + std::to_string(node) + " -> " +
                        std::to_string(neighbor) +
                        " is missing its reverse direction.");
                }
                if (!visited[static_cast<size_t>(neighbor)]) {
                    visited[static_cast<size_t>(neighbor)] = 1;
                    queue.push_back(neighbor);
                }
            }
        }

        int reachable_nodes = 0;
        for (int node = 0; node < node_slots; ++node) {
            const int node_degree = degree[static_cast<size_t>(node)];
            if (node_degree == 0) {
                continue;
            }
            if (!visited[static_cast<size_t>(node)]) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": shadow tree became disconnected at node " +
                    std::to_string(node) + ".");
            }
            ++reachable_nodes;
            if (node != m_dcRootNode && node_degree != 1 && node_degree != 3) {
                throw std::runtime_error(
                    std::string("[DIPPER] ") + context +
                    ": node " + std::to_string(node) +
                    " has invalid degree " + std::to_string(node_degree) +
                    " for a binary shadow tree.");
            }
        }

        if (active_directed_edges != 2 * std::max(0, reachable_nodes - 1)) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": directed-edge count " +
                std::to_string(active_directed_edges) +
                " is inconsistent with reachable node count " +
                std::to_string(reachable_nodes) + ".");
        }
    }

    std::vector<int> _buildDivideAndConquerNodePath(
        int start_node,
        int end_node) const
    {
        const int node_slots = static_cast<int>(m_numSequences * 2);
        if (start_node < 0 || start_node >= node_slots ||
            end_node < 0 || end_node >= node_slots) {
            throw std::runtime_error(
                "[DIPPER] _buildDivideAndConquerNodePath: node id is out of range.");
        }
        std::vector<int> parent(node_slots, -1);
        std::queue<int> pending;
        pending.push(start_node);
        parent[static_cast<size_t>(start_node)] = start_node;
        while (!pending.empty() &&
               parent[static_cast<size_t>(end_node)] < 0) {
            const int node = pending.front();
            pending.pop();
            for (int edge_idx = m_hostDC.h_head[node];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                const int neighbor = m_hostDC.h_e[edge_idx];
                if (neighbor < 0 ||
                    neighbor >= node_slots ||
                    parent[static_cast<size_t>(neighbor)] >= 0) {
                    continue;
                }
                parent[static_cast<size_t>(neighbor)] = node;
                pending.push(neighbor);
            }
        }
        if (parent[static_cast<size_t>(end_node)] < 0) {
            throw std::runtime_error(
                "[DIPPER] _buildDivideAndConquerNodePath: failed to connect anchor labels in the current shadow tree.");
        }

        std::vector<int> path;
        for (int node = end_node; node != start_node;
             node = parent[static_cast<size_t>(node)]) {
            path.push_back(node);
        }
        path.push_back(start_node);
        std::reverse(path.begin(), path.end());
        return path;
    }

    void _buildDivideAndConquerParentView(
        std::vector<int>& parent,
        std::vector<int>& parent_edge) const
    {
        const int node_slots = static_cast<int>(m_numSequences * 2);
        parent.assign(static_cast<size_t>(node_slots), -1);
        parent_edge.assign(static_cast<size_t>(node_slots), -1);
        const int root_idx =
            (m_dcRootNode >= 0)
                ? m_dcRootNode
                : static_cast<int>(m_numSequences);
        if (root_idx < 0 || root_idx >= node_slots) {
            throw std::runtime_error(
                "[DIPPER] _buildDivideAndConquerParentView: root is out of range.");
        }

        std::queue<int> pending;
        pending.push(root_idx);
        parent[static_cast<size_t>(root_idx)] = root_idx;
        while (!pending.empty()) {
            const int node = pending.front();
            pending.pop();
            for (int edge_idx = m_hostDC.h_head[node];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                const int neighbor = m_hostDC.h_e[edge_idx];
                if (neighbor < 0 ||
                    neighbor >= node_slots ||
                    parent[static_cast<size_t>(neighbor)] >= 0) {
                    continue;
                }
                parent[static_cast<size_t>(neighbor)] = node;
                parent_edge[static_cast<size_t>(neighbor)] =
                    _findDivideAndConquerDirectedEdge(node, neighbor);
                pending.push(neighbor);
            }
        }
    }

    int _resolveDivideAndConquerClusterAnchorEdge(int cluster_id) const
    {
        const int cluster_slots =
            static_cast<int>(m_dcClusterAnchorLabels.size());
        if (cluster_id < 0 || cluster_id >= cluster_slots) {
            throw std::runtime_error(
                "[DIPPER] _resolveDivideAndConquerClusterAnchorEdge: cluster id is out of range.");
        }
        const std::pair<int, int> anchor_labels =
            m_dcClusterAnchorLabels[static_cast<size_t>(cluster_id)];
        if (anchor_labels.first < 0 || anchor_labels.second < 0) {
            return cluster_id;
        }
        const auto node_a_it = m_dcLabelToNode.find(anchor_labels.first);
        const auto node_b_it = m_dcLabelToNode.find(anchor_labels.second);
        if (node_a_it == m_dcLabelToNode.end() ||
            node_b_it == m_dcLabelToNode.end()) {
            throw std::runtime_error(
                "[DIPPER] _resolveDivideAndConquerClusterAnchorEdge: cluster anchor labels are unavailable.");
        }

        const std::vector<int> path =
            _buildDivideAndConquerNodePath(
                node_a_it->second,
                node_b_it->second);
        if (path.size() < 2) {
            throw std::runtime_error(
                "[DIPPER] _resolveDivideAndConquerClusterAnchorEdge: anchor path is empty.");
        }

        double total_length = 0.0;
        for (size_t path_idx = 1; path_idx < path.size(); ++path_idx) {
            const int edge_idx =
                _findDivideAndConquerDirectedEdge(
                    path[path_idx - 1],
                    path[path_idx]);
            if (edge_idx < 0) {
                throw std::runtime_error(
                    "[DIPPER] _resolveDivideAndConquerClusterAnchorEdge: anchor path edge is unavailable.");
            }
            total_length += m_hostDC.h_len[edge_idx];
        }

        double remaining = total_length * 0.5;
        for (size_t path_idx = 1; path_idx < path.size(); ++path_idx) {
            const int edge_idx =
                _findDivideAndConquerDirectedEdge(
                    path[path_idx - 1],
                    path[path_idx]);
            const double edge_length = m_hostDC.h_len[edge_idx];
            if (remaining <= edge_length ||
                path_idx + 1 == path.size()) {
                return edge_idx;
            }
            remaining -= edge_length;
        }
        throw std::runtime_error(
            "[DIPPER] _resolveDivideAndConquerClusterAnchorEdge: failed to resolve anchor edge.");
    }

    void _updateDivideAndConquerClosestLeavesAfterSubtreeSpr(
        int moved_root,
        const std::vector<int>& moved_nodes,
        const std::vector<int>& moved_leaves,
        const std::vector<int>& rewired_edges) const
    {
        const int node_slots = static_cast<int>(m_numSequences * 2);
        const int edge_slots = static_cast<int>(m_numSequences * 8);
        _requireDivideAndConquerNodeIndex(
            moved_root,
            "_updateDivideAndConquerClosestLeavesAfterSubtreeSpr(moved root)");
        std::vector<unsigned char> is_moved_leaf(
            static_cast<size_t>(node_slots), 0);
        std::vector<unsigned char> is_moved_node(
            static_cast<size_t>(node_slots), 0);
        for (int moved_node : moved_nodes) {
            _requireDivideAndConquerNodeIndex(
                moved_node,
                "_updateDivideAndConquerClosestLeavesAfterSubtreeSpr(moved node)");
            is_moved_node[static_cast<size_t>(moved_node)] = 1;
        }
        for (int moved_leaf : moved_leaves) {
            _requireDivideAndConquerNodeIndex(
                moved_leaf,
                "_updateDivideAndConquerClosestLeavesAfterSubtreeSpr(moved leaf)");
            is_moved_leaf[static_cast<size_t>(moved_leaf)] = 1;
        }

        std::vector<unsigned char> needs_recompute(
            static_cast<size_t>(edge_slots), 0);
        for (int edge_idx : rewired_edges) {
            _requireDivideAndConquerClosestEdgeIndex(
                edge_idx,
                "_updateDivideAndConquerClosestLeavesAfterSubtreeSpr(rewired edge)");
            needs_recompute[static_cast<size_t>(edge_idx)] = 1;
        }

        // Only the pruned query tip changes its distance to unchanged edges.
        // Remove its stale contribution everywhere while preserving the
        // sorted top-5 entries of every other leaf.
        for (int edge_idx = 0; edge_idx < edge_slots; ++edge_idx) {
            int write_rank = 0;
            for (int rank = 0; rank < 5; ++rank) {
                const size_t slot =
                    static_cast<size_t>(edge_idx) * 5 + rank;
                const int leaf_id = m_hostDC.h_closest_id[slot];
                if (leaf_id >= 0 &&
                    is_moved_leaf[static_cast<size_t>(leaf_id)]) {
                    // The cache stores only five entries.  Once one of them
                    // disappears, the previous sixth-nearest leaf is not
                    // available, so this directed edge must be rebuilt.
                    needs_recompute[static_cast<size_t>(edge_idx)] = 1;
                }
                if (leaf_id < 0 ||
                    is_moved_leaf[static_cast<size_t>(leaf_id)]) {
                    continue;
                }
                const size_t write_slot =
                    static_cast<size_t>(edge_idx) * 5 + write_rank;
                m_hostDC.h_closest_id[write_slot] = leaf_id;
                m_hostDC.h_closest_dis[write_slot] =
                    m_hostDC.h_closest_dis[slot];
                ++write_rank;
            }
            while (write_rank < 5) {
                const size_t slot =
                    static_cast<size_t>(edge_idx) * 5 + write_rank;
                m_hostDC.h_closest_id[slot] = -1;
                m_hostDC.h_closest_dis[slot] = 2.0;
                ++write_rank;
            }
        }

        auto insert_closest = [&](int edge_idx, int leaf_id, double distance) {
            for (int rank = 0; rank < 5; ++rank) {
                const size_t slot =
                    static_cast<size_t>(edge_idx) * 5 + rank;
                if (m_hostDC.h_closest_id[slot] == leaf_id) {
                    if (distance < m_hostDC.h_closest_dis[slot]) {
                        m_hostDC.h_closest_dis[slot] = distance;
                    }
                    return;
                }
                if (distance < m_hostDC.h_closest_dis[slot]) {
                    for (int shift = 4; shift > rank; --shift) {
                        const size_t dst =
                            static_cast<size_t>(edge_idx) * 5 + shift;
                        const size_t src = dst - 1;
                        m_hostDC.h_closest_id[dst] =
                            m_hostDC.h_closest_id[src];
                        m_hostDC.h_closest_dis[dst] =
                            m_hostDC.h_closest_dis[src];
                    }
                    m_hostDC.h_closest_id[slot] = leaf_id;
                    m_hostDC.h_closest_dis[slot] = distance;
                    return;
                }
            }
        };

        // For edges inside the moved component, the direction facing back
        // toward its root contains the outside tree and therefore changes.
        std::vector<int> component_nodes{moved_root};
        std::vector<int> component_parents{-1};
        for (size_t i = 0; i < component_nodes.size(); ++i) {
            const int node = component_nodes[i];
            const int parent = component_parents[i];
            for (int edge_idx = m_hostDC.h_head[node];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                const int neighbor = m_hostDC.h_e[edge_idx];
                if (neighbor == parent) {
                    continue;
                }
                if (is_moved_node[static_cast<size_t>(neighbor)]) {
                    const int reverse = _findDivideAndConquerDirectedEdge(
                        neighbor, node);
                    if (reverse >= 0) {
                        needs_recompute[static_cast<size_t>(reverse)] = 1;
                    }
                    component_nodes.push_back(neighbor);
                    component_parents.push_back(node);
                }
            }
        }

        // Reinsert every moved leaf into unchanged outside-facing caches.
        for (int moved_leaf : moved_leaves) {
            std::vector<int> traversal_nodes{moved_leaf};
            std::vector<int> traversal_parents{-1};
            std::vector<double> traversal_distances{0.0};
            for (size_t index = 0; index < traversal_nodes.size(); ++index) {
                const int node = traversal_nodes[index];
                const int parent = traversal_parents[index];
                const double distance = traversal_distances[index];
                for (int edge_idx = m_hostDC.h_head[node];
                     edge_idx != -1;
                     edge_idx = m_hostDC.h_nxt[edge_idx]) {
                    const int neighbor = m_hostDC.h_e[edge_idx];
                    if (neighbor == parent) {
                        continue;
                    }
                    if (!needs_recompute[static_cast<size_t>(edge_idx)]) {
                        insert_closest(edge_idx, moved_leaf, distance);
                    }
                    traversal_nodes.push_back(neighbor);
                    traversal_parents.push_back(node);
                    traversal_distances.push_back(
                        distance + m_hostDC.h_len[edge_idx]);
                }
            }
        }

        // Reused adjacency slots now represent different directed edges, so
        // all five entries—not just the moved tip—must be recomputed.  Search
        // only the source-side component and stop after the nearest five tips.
        std::vector<int> visit_stamp(static_cast<size_t>(node_slots), 0);
        int stamp = 0;
        using QueueEntry = std::pair<double, int>;
        for (int edge_idx = 0; edge_idx < edge_slots; ++edge_idx) {
            if (!needs_recompute[static_cast<size_t>(edge_idx)]) {
                continue;
            }
            for (int rank = 0; rank < 5; ++rank) {
                const size_t slot =
                    static_cast<size_t>(edge_idx) * 5 + rank;
                m_hostDC.h_closest_id[slot] = -1;
                m_hostDC.h_closest_dis[slot] = 2.0;
            }
            const int source = m_hostDC.h_belong[edge_idx];
            const int excluded_neighbor = m_hostDC.h_e[edge_idx];
            _requireDivideAndConquerNodeIndex(
                source,
                "_updateDivideAndConquerClosestLeavesAfterSubtreeSpr(source)");
            _requireDivideAndConquerNodeIndex(
                excluded_neighbor,
                "_updateDivideAndConquerClosestLeavesAfterSubtreeSpr(neighbor)");
            ++stamp;
            std::priority_queue<
                QueueEntry,
                std::vector<QueueEntry>,
                std::greater<QueueEntry>> pending;
            pending.emplace(0.0, source);
            int found = 0;
            while (!pending.empty() && found < 5) {
                const auto [distance, node] = pending.top();
                pending.pop();
                if (visit_stamp[static_cast<size_t>(node)] == stamp) {
                    continue;
                }
                visit_stamp[static_cast<size_t>(node)] = stamp;
                int degree = 0;
                for (int adjacent = m_hostDC.h_head[node];
                     adjacent != -1;
                     adjacent = m_hostDC.h_nxt[adjacent]) {
                    ++degree;
                }
                if (degree == 1 &&
                    node >= 0 &&
                    node < static_cast<int>(m_numSequences)) {
                    if (distance >= 2.0) {
                        break;
                    }
                    m_hostDC.h_closest_id[
                        static_cast<size_t>(edge_idx) * 5 + found] = node;
                    m_hostDC.h_closest_dis[
                        static_cast<size_t>(edge_idx) * 5 + found] = distance;
                    ++found;
                }
                for (int adjacent = m_hostDC.h_head[node];
                     adjacent != -1;
                     adjacent = m_hostDC.h_nxt[adjacent]) {
                    const int neighbor = m_hostDC.h_e[adjacent];
                    if ((node == source && neighbor == excluded_neighbor) ||
                        visit_stamp[static_cast<size_t>(neighbor)] == stamp) {
                        continue;
                    }
                    pending.emplace(
                        distance + m_hostDC.h_len[adjacent],
                        neighbor);
                }
            }
        }

        if (std::getenv(
                "MLIPPER_VALIDATE_DIPPER_DC_SPR_CLOSEST_INCREMENTAL") !=
            nullptr) {
            std::vector<double> incremental_dis(
                m_hostDC.h_closest_dis,
                m_hostDC.h_closest_dis +
                    static_cast<size_t>(edge_slots) * 5);
            std::vector<int> incremental_id(
                m_hostDC.h_closest_id,
                m_hostDC.h_closest_id +
                    static_cast<size_t>(edge_slots) * 5);
            _recomputeDivideAndConquerClosestLeaves();
            size_t mismatches = 0;
            size_t reported_mismatches = 0;
            for (int node = 0; node < node_slots; ++node) {
                for (int edge_idx = m_hostDC.h_head[node];
                     edge_idx != -1;
                     edge_idx = m_hostDC.h_nxt[edge_idx]) {
                    for (int rank = 0; rank < 5; ++rank) {
                        const size_t slot =
                            static_cast<size_t>(edge_idx) * 5 + rank;
                        if (incremental_id[slot] !=
                                m_hostDC.h_closest_id[slot] ||
                            std::abs(
                                incremental_dis[slot] -
                                m_hostDC.h_closest_dis[slot]) > 1e-10) {
                            if (reported_mismatches < 12) {
                                std::cerr
                                    << "[DIPPER] closest mismatch edge="
                                    << edge_idx << " source=" << node
                                    << " neighbor=" << m_hostDC.h_e[edge_idx]
                                    << " rank=" << rank
                                    << " rebuilt="
                                    << static_cast<int>(needs_recompute[
                                           static_cast<size_t>(edge_idx)])
                                    << " incremental=("
                                    << incremental_id[slot] << ","
                                    << incremental_dis[slot] << ") full=("
                                    << m_hostDC.h_closest_id[slot] << ","
                                    << m_hostDC.h_closest_dis[slot] << ")"
                                    << " incremental_is_moved="
                                    << (incremental_id[slot] >= 0
                                        ? static_cast<int>(is_moved_leaf[
                                              static_cast<size_t>(incremental_id[slot])])
                                        : 0)
                                    << " full_is_moved="
                                    << (m_hostDC.h_closest_id[slot] >= 0
                                        ? static_cast<int>(is_moved_leaf[
                                              static_cast<size_t>(m_hostDC.h_closest_id[slot])])
                                        : 0)
                                    << " moved_count=" << moved_leaves.size()
                                    << "\n";
                                ++reported_mismatches;
                            }
                            ++mismatches;
                        }
                    }
                }
            }
            if (mismatches != 0) {
                throw std::runtime_error(
                    "[DIPPER] incremental SPR closest-leaf update disagrees with full recomputation in " +
                    std::to_string(mismatches) + " entries.");
            }
            std::cerr
                << "[DIPPER] Incremental SPR closest-leaf validation passed.\n";
        }
    }

    void _recomputeDivideAndConquerClosestLeaves() const
    {
        const size_t closest_edge_slots = m_numSequences * 8;
        for (size_t edge_idx = 0; edge_idx < closest_edge_slots; ++edge_idx) {
            for (int rank = 0; rank < 5; ++rank) {
                m_hostDC.h_closest_dis[edge_idx * 5 + rank] = 2.0;
                m_hostDC.h_closest_id[edge_idx * 5 + rank] = -1;
            }
        }

        const int node_slots = static_cast<int>(m_numSequences * 2);
        int root = m_dcRootNode;
        if (root < 0 || root >= node_slots || m_hostDC.h_head[root] == -1) {
            root = -1;
            for (int node = 0; node < node_slots; ++node) {
                if (m_hostDC.h_head[node] != -1) {
                    root = node;
                    break;
                }
            }
        }
        if (root < 0) {
            return;
        }

        std::vector<int> degree(static_cast<size_t>(node_slots), 0);
        for (int node_idx = 0; node_idx < node_slots; ++node_idx) {
            for (int edge_idx = m_hostDC.h_head[node_idx];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                if (edge_idx < 0 ||
                    edge_idx >= static_cast<int>(closest_edge_slots)) {
                    throw std::runtime_error(
                        "[DIPPER] _recomputeDivideAndConquerClosestLeaves: active edge id exceeds closest-array capacity.");
                }
                ++degree[static_cast<size_t>(node_idx)];
            }
        }

        std::vector<int> parent(static_cast<size_t>(node_slots), -2);
        std::vector<int> parent_out_edge(static_cast<size_t>(node_slots), -1);
        std::vector<int> parent_in_edge(static_cast<size_t>(node_slots), -1);
        std::vector<int> order;
        order.reserve(static_cast<size_t>(node_slots));
        parent[static_cast<size_t>(root)] = -1;
        order.push_back(root);
        for (size_t i = 0; i < order.size(); ++i) {
            const int node = order[i];
            for (int edge_idx = m_hostDC.h_head[node];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                const int neighbor = m_hostDC.h_e[edge_idx];
                if (neighbor == parent[static_cast<size_t>(node)]) {
                    continue;
                }
                if (parent[static_cast<size_t>(neighbor)] != -2) {
                    continue;
                }
                parent[static_cast<size_t>(neighbor)] = node;
                parent_out_edge[static_cast<size_t>(neighbor)] = edge_idx;
                parent_in_edge[static_cast<size_t>(neighbor)] =
                    _findDivideAndConquerDirectedEdge(neighbor, node);
                order.push_back(neighbor);
            }
        }

        using ClosestEntry = std::pair<double, int>;
        auto write_message = [&](int edge_idx,
                                 std::vector<ClosestEntry> candidates) {
            if (edge_idx < 0) {
                return;
            }
            std::sort(candidates.begin(), candidates.end(),
                [](const ClosestEntry& lhs, const ClosestEntry& rhs) {
                    if (lhs.first != rhs.first) return lhs.first < rhs.first;
                    return lhs.second < rhs.second;
                });
            int rank = 0;
            int previous_id = -1;
            for (const auto& candidate : candidates) {
                if (candidate.first >= 2.0 || candidate.second == previous_id) {
                    continue;
                }
                const size_t slot = static_cast<size_t>(edge_idx) * 5 + rank;
                m_hostDC.h_closest_dis[slot] = candidate.first;
                m_hostDC.h_closest_id[slot] = candidate.second;
                previous_id = candidate.second;
                if (++rank == 5) break;
            }
        };
        auto append_message = [&](std::vector<ClosestEntry>& candidates,
                                  int edge_idx,
                                  double offset) {
            if (edge_idx < 0) return;
            for (int rank = 0; rank < 5; ++rank) {
                const size_t slot = static_cast<size_t>(edge_idx) * 5 + rank;
                if (m_hostDC.h_closest_id[slot] < 0) break;
                candidates.emplace_back(
                    m_hostDC.h_closest_dis[slot] + offset,
                    m_hostDC.h_closest_id[slot]);
            }
        };

        // Postorder: source-side messages from each child toward its parent.
        for (auto it = order.rbegin(); it != order.rend(); ++it) {
            const int node = *it;
            if (node == root) continue;
            std::vector<ClosestEntry> candidates;
            if (degree[static_cast<size_t>(node)] == 1) {
                candidates.emplace_back(0.0, node);
            }
            for (int edge_idx = m_hostDC.h_head[node];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                const int neighbor = m_hostDC.h_e[edge_idx];
                if (neighbor == parent[static_cast<size_t>(node)]) continue;
                append_message(candidates,
                    parent_in_edge[static_cast<size_t>(neighbor)],
                    m_hostDC.h_len[edge_idx]);
            }
            write_message(parent_in_edge[static_cast<size_t>(node)],
                          std::move(candidates));
        }

        // Preorder: complement messages from each parent toward its child.
        for (int node : order) {
            for (int edge_idx = m_hostDC.h_head[node];
                 edge_idx != -1;
                 edge_idx = m_hostDC.h_nxt[edge_idx]) {
                const int child = m_hostDC.h_e[edge_idx];
                if (parent[static_cast<size_t>(child)] != node) continue;
                std::vector<ClosestEntry> candidates;
                if (degree[static_cast<size_t>(node)] == 1) {
                    candidates.emplace_back(0.0, node);
                }
                if (parent[static_cast<size_t>(node)] >= 0) {
                    append_message(candidates,
                        parent_out_edge[static_cast<size_t>(node)],
                        m_hostDC.h_len[parent_in_edge[static_cast<size_t>(node)]]);
                }
                for (int sibling_edge = m_hostDC.h_head[node];
                     sibling_edge != -1;
                     sibling_edge = m_hostDC.h_nxt[sibling_edge]) {
                    const int sibling = m_hostDC.h_e[sibling_edge];
                    if (sibling == child ||
                        sibling == parent[static_cast<size_t>(node)]) continue;
                    append_message(candidates,
                        parent_in_edge[static_cast<size_t>(sibling)],
                        m_hostDC.h_len[sibling_edge]);
                }
                write_message(edge_idx, std::move(candidates));
            }
        }

        if (std::getenv("MLIPPER_VALIDATE_DIPPER_DC_CLOSEST_LINEAR") != nullptr) {
            std::vector<double> linear_dis(
                m_hostDC.h_closest_dis,
                m_hostDC.h_closest_dis + closest_edge_slots * 5);
            std::vector<int> linear_id(
                m_hostDC.h_closest_id,
                m_hostDC.h_closest_id + closest_edge_slots * 5);
            for (size_t edge_idx = 0; edge_idx < closest_edge_slots; ++edge_idx) {
                for (int rank = 0; rank < 5; ++rank) {
                    m_hostDC.h_closest_dis[edge_idx * 5 + rank] = 2.0;
                    m_hostDC.h_closest_id[edge_idx * 5 + rank] = -1;
                }
            }
            std::vector<int> edge_mask_index(closest_edge_slots, -1);
            std::vector<int> leaf_nodes;
            for (int node = 0; node < node_slots; ++node) {
                for (int edge_idx = m_hostDC.h_head[node];
                     edge_idx != -1;
                     edge_idx = m_hostDC.h_nxt[edge_idx]) {
                    edge_mask_index[static_cast<size_t>(edge_idx)] = edge_idx;
                }
                if (degree[static_cast<size_t>(node)] == 1) {
                    leaf_nodes.push_back(node);
                }
            }
            std::vector<int> bfs_id(static_cast<size_t>(node_slots), -1);
            std::vector<int> bfs_from(static_cast<size_t>(node_slots), -1);
            std::vector<double> bfs_dis(static_cast<size_t>(node_slots), 0.0);
            for (int leaf : leaf_nodes) {
                updateClosestNodesCpu(
                    m_hostDC.h_head, m_hostDC.h_nxt, m_hostDC.h_e,
                    m_hostDC.h_len, m_hostDC.h_closest_dis,
                    m_hostDC.h_closest_id, leaf, bfs_id.data(),
                    bfs_from.data(), bfs_dis.data(), edge_mask_index.data());
            }
            size_t mismatches = 0;
            for (int node : order) {
                for (int edge_idx = m_hostDC.h_head[node];
                     edge_idx != -1;
                     edge_idx = m_hostDC.h_nxt[edge_idx]) {
                    for (int rank = 0; rank < 5; ++rank) {
                        const size_t slot = static_cast<size_t>(edge_idx) * 5 + rank;
                        if (linear_id[slot] != m_hostDC.h_closest_id[slot] ||
                            std::abs(linear_dis[slot] - m_hostDC.h_closest_dis[slot]) > 1e-10) {
                            if (mismatches < 12) {
                                std::cerr << "[DIPPER] linear closest mismatch edge="
                                          << edge_idx << " rank=" << rank
                                          << " linear=(" << linear_id[slot] << ","
                                          << linear_dis[slot] << ") legacy=("
                                          << m_hostDC.h_closest_id[slot] << ","
                                          << m_hostDC.h_closest_dis[slot] << ")\n";
                            }
                            ++mismatches;
                        }
                    }
                }
            }
            if (mismatches != 0) {
                throw std::runtime_error(
                    "[DIPPER] linear closest-leaf recomputation disagrees with legacy in " +
                    std::to_string(mismatches) + " entries.");
            }
            std::copy(linear_dis.begin(), linear_dis.end(), m_hostDC.h_closest_dis);
            std::copy(linear_id.begin(), linear_id.end(), m_hostDC.h_closest_id);
            std::cerr << "[DIPPER] Linear closest-leaf validation passed.\n";
        }
    }

    void _rebuildDivideAndConquerActiveClusterState() const
    {
        if (m_dcActiveClusterId < 0) {
            return;
        }

        const int anchor_edge_id =
            _resolveDivideAndConquerClusterAnchorEdge(m_dcActiveClusterId);
        const int current_cluster_cursor =
            static_cast<int>(m_dcCurrentClusterCursor);
        if (current_cluster_cursor < 0 ||
            current_cluster_cursor >=
                static_cast<int>(m_dcClusterPlacementOffsets.size())) {
            throw std::runtime_error(
                "[DIPPER] _rebuildDivideAndConquerActiveClusterState: cluster cursor is out of range.");
        }

        const int placement_offset =
            m_dcClusterPlacementOffsets[static_cast<size_t>(
                current_cluster_cursor)];
        const std::vector<int>& cluster_queries =
            m_dcClusterQueries[static_cast<size_t>(m_dcActiveClusterId)];

        m_dcActiveEdgeMask.assign(m_numSequences * 4, -1);
        m_dcActiveLeafMask.assign(m_numSequences * 2, -1);
        m_dcActiveLeafMap.assign(m_numSequences * 2, -1);
        m_dcActiveId.assign(m_numSequences * 2, -1);
        m_dcActiveFrom.assign(m_numSequences * 2, -1);
        m_dcActiveDis.assign(m_numSequences * 2, 0.0);
        m_dcActiveDist.assign(m_numSequences, 0.0);
        m_dcMinPos.assign(
            m_numSequences * 4,
            std::make_tuple(0, 0.0, 2.0));

        _initializeDivideAndConquerClusterMasks(
            anchor_edge_id,
            m_dcActiveEdgeMask,
            m_dcActiveLeafMask,
            "_rebuildDivideAndConquerActiveClusterState");
        _refreshDivideAndConquerSeedLeafMaps(
            10,
            "_rebuildDivideAndConquerActiveClusterState(seed leaves)");

        const std::pair<int, int> anchor_labels =
            m_dcClusterAnchorLabels[static_cast<size_t>(m_dcActiveClusterId)];
        const auto anchor_a_it =
            m_dcLabelToNode.find(anchor_labels.first);
        const auto anchor_b_it =
            m_dcLabelToNode.find(anchor_labels.second);
        if (anchor_a_it == m_dcLabelToNode.end() ||
            anchor_b_it == m_dcLabelToNode.end()) {
            throw std::runtime_error(
                "[DIPPER] _rebuildDivideAndConquerActiveClusterState: cluster anchor labels are unavailable.");
        }
        const std::vector<int> anchor_path =
            _buildDivideAndConquerNodePath(
                anchor_a_it->second,
                anchor_b_it->second);
        std::vector<char> anchor_node_mask(
            static_cast<size_t>(m_numSequences * 2),
            0);
        std::vector<char> edge_selected(
            static_cast<size_t>(m_numSequences * 8),
            0);
        auto append_edge = [&](int edge_idx) {
            if (edge_idx < 0 ||
                edge_idx >= static_cast<int>(edge_selected.size()) ||
                edge_selected[static_cast<size_t>(edge_idx)]) {
                return;
            }
            if (m_dcActiveEdgeCount >=
                static_cast<int>(m_dcActiveEdgeMask.size())) {
                throw std::runtime_error(
                    "[DIPPER] _rebuildDivideAndConquerActiveClusterState: active edge mask capacity exceeded.");
            }
            edge_selected[static_cast<size_t>(edge_idx)] = 1;
            m_dcActiveEdgeMask[static_cast<size_t>(m_dcActiveEdgeCount++)] =
                edge_idx;
        };
        m_dcActiveEdgeCount = 0;
        for (size_t path_idx = 0; path_idx < anchor_path.size(); ++path_idx) {
            anchor_node_mask[static_cast<size_t>(anchor_path[path_idx])] = 1;
            if (path_idx == 0) {
                continue;
            }
            const int parent_node = anchor_path[path_idx - 1];
            const int child_node = anchor_path[path_idx];
            append_edge(
                _findDivideAndConquerDirectedEdge(parent_node, child_node));
            append_edge(
                _findDivideAndConquerDirectedEdge(child_node, parent_node));
        }
        if (!anchor_path.empty() &&
            m_dcRootNode >= 0 &&
            m_dcRootNode < static_cast<int>(anchor_node_mask.size()) &&
            !anchor_node_mask[static_cast<size_t>(m_dcRootNode)]) {
            const std::vector<int> root_connector_path =
                _buildDivideAndConquerNodePath(
                    m_dcRootNode,
                    anchor_path.front());
            for (size_t path_idx = 0;
                 path_idx < root_connector_path.size();
                 ++path_idx) {
                anchor_node_mask[
                    static_cast<size_t>(root_connector_path[path_idx])] = 1;
                if (path_idx == 0) {
                    continue;
                }
                const int parent_node = root_connector_path[path_idx - 1];
                const int child_node = root_connector_path[path_idx];
                append_edge(
                    _findDivideAndConquerDirectedEdge(parent_node, child_node));
                append_edge(
                    _findDivideAndConquerDirectedEdge(child_node, parent_node));
            }
        }

        std::vector<int> parent;
        std::vector<int> parent_edge;
        _buildDivideAndConquerParentView(parent, parent_edge);

        m_dcActiveLeafCount = 10;
        const size_t processed_queries =
            std::min(m_dcActiveClusterOffset, cluster_queries.size());
        for (size_t query_offset = 0;
             query_offset < processed_queries;
             ++query_offset) {
            const int query_idx = cluster_queries[query_offset];
            if (m_dcActiveLeafCount >=
                static_cast<int>(m_dcActiveLeafMask.size()) ||
                m_dcActiveLeafCount >=
                    static_cast<int>(m_dcActiveLeafMap.size())) {
                throw std::runtime_error(
                    "[DIPPER] _rebuildDivideAndConquerActiveClusterState: active leaf mask capacity exceeded.");
            }
            m_dcActiveLeafMask[static_cast<size_t>(m_dcActiveLeafCount)] =
                query_idx;
            m_dcActiveLeafMap[static_cast<size_t>(m_dcActiveLeafCount)] =
                query_idx;
            ++m_dcActiveLeafCount;

            int node = query_idx;
            while (node >= 0 &&
                   node < static_cast<int>(parent.size()) &&
                   !anchor_node_mask[static_cast<size_t>(node)]) {
                const int edge_down =
                    parent_edge[static_cast<size_t>(node)];
                const int parent_node =
                    parent[static_cast<size_t>(node)];
                if (edge_down < 0 || parent_node < 0) {
                    throw std::runtime_error(
                        "[DIPPER] _rebuildDivideAndConquerActiveClusterState: failed to reconnect a processed cluster query to the anchor path.");
                }
                append_edge(edge_down);
                append_edge(
                    _findDivideAndConquerDirectedEdge(node, parent_node));
                node = parent_node;
            }
        }

        m_dcActiveIdx =
            static_cast<int>(m_dcBackboneSize) * 4 +
            (placement_offset + static_cast<int>(processed_queries)) * 4;
        m_dcActiveInsertLeafCount =
            static_cast<int>(m_dcBackboneSize) +
            placement_offset + static_cast<int>(processed_queries);
        _uploadDivideAndConquerActiveClusterSequences();
    }

    bool _startNextDivideAndConquerCluster() const
    {
        if (m_dcCurrentClusterCursor >= m_dcOrderedClusters.size()) {
            m_dcActiveClusterId = -1;
            m_dcActiveClusterOffset = 0;
            return false;
        }

        m_dcActiveClusterId =
            m_dcOrderedClusters[m_dcCurrentClusterCursor];
        m_dcActiveClusterOffset = 0;
        m_dcActiveEdgeCount = 2;
        m_dcActiveLeafCount = 10;
        m_dcActiveEdgeMask.assign(m_numSequences * 4, -1);
        m_dcActiveLeafMask.assign(m_numSequences * 2, -1);
        m_dcActiveLeafMap.assign(m_numSequences * 2, -1);
        m_dcActiveId.assign(m_numSequences * 2, -1);
        m_dcActiveFrom.assign(m_numSequences * 2, -1);
        m_dcActiveDis.assign(m_numSequences * 2, 0.0);
        m_dcActiveDist.assign(m_numSequences, 0.0);
        m_dcMinPos.assign(
            m_numSequences * 4,
            std::make_tuple(0, 0.0, 2.0));

        const int anchor_edge_id =
            _resolveDivideAndConquerClusterAnchorEdge(m_dcActiveClusterId);
        _initializeDivideAndConquerClusterMasks(
            anchor_edge_id,
            m_dcActiveEdgeMask,
            m_dcActiveLeafMask,
            "_startNextDivideAndConquerCluster");
        _refreshDivideAndConquerSeedLeafMaps(
            m_dcActiveLeafCount,
            "_startNextDivideAndConquerCluster(seed leaves)");

        const int placement_offset =
            m_dcClusterPlacementOffsets[m_dcCurrentClusterCursor];
        m_dcActiveIdx = static_cast<int>(m_dcBackboneSize) * 4 +
                        placement_offset * 4;
        m_dcActiveInsertLeafCount =
            static_cast<int>(m_dcBackboneSize) + placement_offset;
        _uploadDivideAndConquerActiveClusterSequences();
        return true;
    }

    void _ensureDivideAndConquerSchedulerReady() const
    {
        _ensureDivideAndConquerClustersReady();
        if (m_dcSchedulerReady) {
            return;
        }
        _copyDivideAndConquerTreeStateToHost();
        if (m_dcRootNode < 0) {
            m_dcRootNode = static_cast<int>(m_numSequences);
        }
        _assignDivideAndConquerShadowLabels();
        _validateDivideAndConquerShadowTree(
            "_ensureDivideAndConquerSchedulerReady(initial_shadow_tree)");
        _buildDivideAndConquerClusterSchedule();
        _ensureDivideAndConquerDeviceScratch();
        m_dcProcessedQueryCount = 0;
        m_dcInsertedQueryLeafIds.assign(m_numSequences, 0);
        m_dcCurrentClusterCursor = 0;
        m_dcActiveClusterId = -1;
        m_dcActiveClusterOffset = 0;
        m_dcPendingStepReady = false;
        m_dcPendingQueryInfo = {};
        m_dcPendingLeafId = -1;
        m_dcPendingEdgeId = -1;
        m_dcPendingFracLen = 0.0;
        m_dcPendingAddLen = 0.0;
        m_dcSchedulerReady = true;
    }

    void _copyActiveDivideAndConquerLeafStateToDevice() const
    {
        const size_t bytes =
            static_cast<size_t>(m_dcActiveLeafCount) * sizeof(int);
        _checkCudaStatus(
            cudaMemcpy(
                m_dcDeviceLeafIds,
                m_dcActiveLeafMask.data(),
                bytes,
                cudaMemcpyHostToDevice),
            "_copyActiveDivideAndConquerLeafStateToDevice(leafIds)");
        _checkCudaStatus(
            cudaMemcpy(
                m_dcDeviceLeafMap,
                m_dcActiveLeafMap.data(),
                bytes,
                cudaMemcpyHostToDevice),
            "_copyActiveDivideAndConquerLeafStateToDevice(leafMap)");
    }

    void _ensureDivideAndConquerPendingStepReady() const
    {
        _ensureDivideAndConquerSchedulerReady();
        if (m_dcPendingStepReady) {
            return;
        }
        const size_t scheduled_count = m_dcScheduledQueryNames.empty()
            ? m_numSequences - m_dcBackboneSize
            : m_dcScheduledQueryNames.size();
        if (m_dcProcessedQueryCount >= scheduled_count) {
            return;
        }
        if (m_dcActiveClusterId < 0 && !_startNextDivideAndConquerCluster()) {
            throw std::runtime_error(
                "[DIPPER] _ensureDivideAndConquerPendingStepReady: no active cluster is available.");
        }

        const std::vector<int>& cluster_queries =
            m_dcClusterQueries[static_cast<size_t>(m_dcActiveClusterId)];
        if (m_dcActiveClusterOffset >= cluster_queries.size()) {
            throw std::runtime_error(
                "[DIPPER] _ensureDivideAndConquerPendingStepReady: active cluster cursor is out of range.");
        }
        const int query_idx =
            cluster_queries[m_dcActiveClusterOffset];
        _requireDivideAndConquerActiveLeafMaskReady(
            "_ensureDivideAndConquerPendingStepReady(active leaf mask)");
        _copyActiveDivideAndConquerLeafStateToDevice();

        m_msaArraysDC.distSpecialIDConstructionOnGpuDC(
            *m_params,
            query_idx,
            m_kpArraysDC.d_dist,
            m_dcActiveLeafCount,
            m_dcDeviceLeafIds,
            m_dcDeviceLeafMap);
        _checkCudaStatus(
            cudaDeviceSynchronize(),
            "_ensureDivideAndConquerPendingStepReady(distSpecialID)");
        _checkCudaStatus(
            cudaMemcpy(
                m_dcActiveDist.data(),
                m_kpArraysDC.d_dist,
                m_numSequences * sizeof(double),
                cudaMemcpyDeviceToHost),
            "_ensureDivideAndConquerPendingStepReady(distCopy)");

        _requireDivideAndConquerActiveEdgeMaskReady(
            m_dcActiveEdgeCount,
            "_ensureDivideAndConquerPendingStepReady(active edge mask)");
        calculateBranchLengthSpecialIDCpu(
            m_dcActiveClusterId,
            m_hostDC.h_head,
            m_hostDC.h_nxt,
            m_dcActiveDist,
            m_hostDC.h_e,
            m_hostDC.h_len,
            m_hostDC.h_belong,
            m_dcMinPos.data(),
            static_cast<int>(m_numSequences) * 4 - 4,
            m_hostDC.h_closest_dis,
            m_hostDC.h_closest_id,
            m_dcActiveEdgeCount,
            m_dcActiveEdgeMask);
        const auto iter = std::min_element(
            m_dcMinPos.begin(),
            m_dcMinPos.begin() + m_dcActiveEdgeCount,
            [](const auto& lhs, const auto& rhs) {
                return std::get<2>(lhs) < std::get<2>(rhs);
            });
        if (iter == m_dcMinPos.begin() + m_dcActiveEdgeCount) {
            throw std::runtime_error(
                "[DIPPER] _ensureDivideAndConquerPendingStepReady: failed to find a candidate edge.");
        }

        const int candidate_edge_id = std::get<0>(*iter);
        _requireDivideAndConquerClosestEdgeIndex(
            candidate_edge_id,
            "_ensureDivideAndConquerPendingStepReady(candidate edge)");
        const int node_a = m_hostDC.h_belong[candidate_edge_id];
        const int node_b = m_hostDC.h_e[candidate_edge_id];
        _requireDivideAndConquerNodeIndex(
            node_a,
            "_ensureDivideAndConquerPendingStepReady(candidate node A)");
        _requireDivideAndConquerNodeIndex(
            node_b,
            "_ensureDivideAndConquerPendingStepReady(candidate node B)");
        const int label_a = nodeLabel(node_a);
        const int label_b = nodeLabel(node_b);
        if (label_a < 0 || label_b < 0) {
            throw std::runtime_error(
                "[DIPPER] _ensureDivideAndConquerPendingStepReady: candidate-edge labels are unavailable.");
        }

        const auto& raw_name_to_index =
            _buildRawNameToIndexMap("exportPendingQuery");
        const std::string& name =
            m_sortedNames[static_cast<size_t>(query_idx)];
        const auto raw_it = raw_name_to_index.find(name);
        if (raw_it == raw_name_to_index.end()) {
            throw std::runtime_error(
                "[DIPPER] _ensureDivideAndConquerPendingStepReady: missing sequence for pending query '" +
                name + "'.");
        }

        m_dcPendingQueryInfo = {};
        m_dcPendingQueryInfo.sequenceIdx = query_idx;
        m_dcPendingQueryInfo.clusterId = m_dcActiveClusterId;
        m_dcPendingQueryInfo.name = name;
        m_dcPendingQueryInfo.sequence = m_rawSeqs[raw_it->second];
        m_dcPendingQueryInfo.splitLabelA = label_a;
        m_dcPendingQueryInfo.splitLabelB = label_b;
        m_dcPendingQueryInfo.lenA = std::get<1>(*iter);
        m_dcPendingQueryInfo.lenB = std::max(
            0.0,
            m_hostDC.h_len[candidate_edge_id] - std::get<1>(*iter));
        m_dcPendingQueryInfo.lenTip = std::get<2>(*iter);
        m_dcPendingLeafId = query_idx;
        m_dcPendingEdgeId = candidate_edge_id;
        m_dcPendingFracLen = std::get<1>(*iter);
        m_dcPendingAddLen = std::get<2>(*iter);
        m_dcPendingStepReady = true;
    }

    void _ensureDivideAndConquerClustersReady() const
    {
        _requireOneOf(
            {State::GPUReady, State::Done},
            "buildDivideAndConquerClusters");
        if (m_activePipeline != ActivePipeline::DivideAndConquer) {
            throw std::logic_error(
                "[DIPPER] buildDivideAndConquerClusters: divide-and-conquer pipeline is not active.");
        }
        if (!m_dcBackboneBuilt) {
            const_cast<MashPlacement::KPlacementDeviceArraysDC&>(m_kpArraysDC)
                .findBackboneTreeDC(
                    *m_params,
                    m_mashArraysDC,
                    const_cast<MashPlacement::MatrixReader&>(m_matrixReader),
                    m_msaArraysDC,
                    m_hostDC);
            const_cast<DipperSession*>(this)->m_dcBackboneBuilt = true;
        }
        if (m_queryClusterIdsReady &&
            m_computedQueryClusterTopK >= m_requestedQueryClusterTopK) {
            return;
        }

        if (m_hostDC.clusterID != nullptr) {
            delete[] m_hostDC.clusterID;
            m_hostDC.clusterID = nullptr;
        }
        m_hostDC.requestedClusterTopK = m_requestedQueryClusterTopK;

        const_cast<MashPlacement::KPlacementDeviceArraysDC&>(m_kpArraysDC)
            .findClustersDC(
                *m_params,
                m_mashArraysDC,
                const_cast<MashPlacement::MatrixReader&>(m_matrixReader),
                m_msaArraysDC,
                const_cast<MashPlacement::KPlacementDeviceArraysHostDC&>(m_hostDC));

        std::vector<int> cluster_ids(m_numSequences, -1);
        for (size_t idx = m_dcBackboneSize; idx < m_numSequences; ++idx) {
            cluster_ids[idx] = m_hostDC.clusterID[idx];
        }
        const_cast<DipperSession*>(this)->m_queryClusterIds = std::move(cluster_ids);
        const_cast<DipperSession*>(this)->m_queryEffectiveClusterIds =
            m_queryClusterIds;
        const_cast<DipperSession*>(this)->m_queryClusterIdsReady = true;
        const_cast<DipperSession*>(this)->m_computedQueryClusterTopK =
            m_requestedQueryClusterTopK;
        const_cast<DipperSession*>(this)->m_queryClusterBackboneSize =
            m_dcBackboneSize;
    }

    const std::unordered_map<std::string, size_t>& _buildRawNameToIndexMap(
        const char* context) const
    {
        if (m_rawNameToIndex.size() != m_rawNames.size()) {
            throw std::runtime_error(
                std::string("[DIPPER] ") + context +
                ": raw sequence-name index is inconsistent.");
        }
        return m_rawNameToIndex;
    }

    void _ensureQueryClusterAssignments() const
    {
        _requireOneOf(
            {State::GPUReady, State::Done},
            "exportQueryClusterAssignments");
        if (m_activePipeline == ActivePipeline::DivideAndConquer) {
            _ensureDivideAndConquerClustersReady();
            return;
        }
        if (m_queryClusterIdsReady) {
            return;
        }
        if (m_tree == nullptr || m_backboneSize == 0) {
            throw std::logic_error(
                "[DIPPER] exportQueryClusterAssignments: a backbone tree is required.");
        }
        if (m_params == nullptr || m_fourBitSeqs == nullptr || m_seqLengths == nullptr) {
            throw std::logic_error(
                "[DIPPER] exportQueryClusterAssignments: session is missing initialized sequence buffers.");
        }

        MashPlacement::Param dc_params = *m_params;
        dc_params.in = "m";
        dc_params.out = "t";
        dc_params.totalNumSeqs = m_numSequences;
        dc_params.backboneSize = m_backboneSize;
        dc_params.batchSize = std::max<uint64_t>(
            1,
            static_cast<uint64_t>(m_backboneSize));

        MashPlacement::MashDeviceArraysDC mash_arrays_dc{};
        MashPlacement::MSADeviceArraysDC msa_arrays_dc{};
        MashPlacement::MatrixReader matrix_reader{};
        MashPlacement::KPlacementDeviceArraysDC kp_arrays_dc{};
        MashPlacement::KPlacementDeviceArraysHostDC host_dc{};

        const auto cleanup = [&]() {
            if (host_dc.clusterID != nullptr) {
                delete[] host_dc.clusterID;
                host_dc.clusterID = nullptr;
            }
            msa_arrays_dc.deallocateDeviceArraysDC();
            kp_arrays_dc.deallocateDeviceArraysDC();
        };

        try {
            msa_arrays_dc.allocateDeviceArraysDC(
                m_fourBitSeqs,
                m_seqLengths,
                m_numSequences,
                dc_params);
            kp_arrays_dc.allocateDeviceArraysDC(
                m_backboneSize,
                m_numSequences);
            kp_arrays_dc.initializeDeviceArraysFromBackboneTreeDC(m_tree);
            kp_arrays_dc.findClustersDC(
                dc_params,
                mash_arrays_dc,
                matrix_reader,
                msa_arrays_dc,
                host_dc);

            std::vector<int> cluster_ids(m_numSequences, -1);
            for (size_t idx = m_backboneSize; idx < m_numSequences; ++idx) {
                cluster_ids[idx] = host_dc.clusterID[idx];
            }
            m_queryClusterIds = std::move(cluster_ids);
            m_queryClusterIdsReady = true;
            m_queryClusterBackboneSize = m_backboneSize;
        } catch (...) {
            cleanup();
            throw;
        }

        cleanup();
    }

    static std::string _joinCladeSignature(
        const std::vector<std::string>& tip_names)
    {
        std::string signature;
        for (size_t idx = 0; idx < tip_names.size(); ++idx) {
            if (idx > 0) {
                signature.push_back('\n');
            }
            signature += tip_names[idx];
        }
        return signature;
    }
};
