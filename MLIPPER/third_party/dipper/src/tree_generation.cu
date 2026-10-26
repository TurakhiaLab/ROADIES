/**
 * tree_generation.cu
 *
 * DIPPER main entry point and phylogenetic tree generation pipeline.
 * Supports multiple input formats (distance matrix, unaligned/aligned FASTA),
 * algorithms (placement, conventional NJ, divide-and-conquer), and output modes.
 * Orchestrates sequence I/O, compression, GPU sketch construction, and tree building.
 */

#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <boost/program_options.hpp> 
#include "../src/kseq.h"
#include "zlib.h"
#include <cuda_runtime.h>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include "version.hpp"

#ifndef TWOBITCOMPRESSOR_HPP
#include "../src/twoBitCompressor.hpp"
#endif

#ifndef FOURBITCOMPRESSOR_HPP
#include "../src/fourBitCompressor.hpp"
#endif

#ifndef MASHPL_CUH
#include "../src/mash_placement.cuh"
#endif

namespace po = boost::program_options;

KSEQ_INIT2(, gzFile, gzread)

po::options_description mainDesc("DIPPER Command Line Arguments");

/** Define and register all command-line options (required, optional, help). */
void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    po::options_description requiredDesc("Required Options");
    requiredDesc.add_options()
        ("input-format,i",     po::value<std::string>()->required(),
        "Input format:\n"
        "  d - distance matrix in PHYLIP format\n"
        "  r - unaligned sequences in FASTA format\n"
        "  m - aligned sequences in FASTA format")

        ("input-file,I",       po::value<std::string>()->required(),
        "Input file path:\n"
        "  PHYLIP format for distance matrix\n"
        "  FASTA format for aligned or unaligned sequences")

        ("output-file,O",      po::value<std::string>()->required(),
        "Output file path");


    po::options_description optionalDesc("Optional Options");
    optionalDesc.add_options()
        ("output-format,o",    po::value<std::string>(),
        "Output format:\n"
        "  t - phylogenetic tree in Newick format (default)\n"
        "  d - distance matrix in PHYLIP format (coming soon)")

        ("algorithm,m",        po::value<std::string>(),
        "Algorithm selection:\n"
        "  0 - default mode\n"
        "  1 - force placement\n"
        "  2 - force conventional NJ\n"
        "  3 - force divide-and-conquer")

        ("K-closest,K",   po::value<std::string>(),
        "Placement mode:\n"
        "  -1 - exact mode\n"
        "  10 - default")

        ("kmer-size,k",        po::value<std::string>(),
        "K-mer size:\n"
        "  Valid range: 2-15 (default: 15)")

        ("sketch-size,s",      po::value<std::string>(),
        "Sketch size (default: 1000)")

        ("distance-type,d",    po::value<std::string>(),
        "Distance type to calculate:\n"
        "  1 - uncorrected\n"
        "  2 - JC (default)\n"
        "  3 - Tajima-Nei\n"
        "  4 - K2P\n"
        "  5 - Tamura\n"
        "  6 - Jinnei")

        ("add,a",
        "Add query to backbone using k-closest placement")

        ("protein,p",
        "Input sequences are protein sequences (default: DNA/RNA sequences)")
        
        ("protein-reduced",
        "Use reduced 4-group 2-bit protein alphabet for raw input (-i r): "
        "L=LVIMC, A=AGSTP, F=FYW, E=EDNQKRH")
        
        ("protein-murphy8",
        "Use Murphy8 3-bit protein alphabet for raw input (-i r)")

        ("input-tree,t",       po::value<std::string>(),
        "Input backbone tree (Newick format), required with --add option")

        ("input-tree-2",       po::value<std::string>(),
        "Input backbone tree (Newick format), required with --node-dist option")

        ("range,r",        po::value<std::string>(),
        "Restrict processing to a subset of alignment coordinates. Provide as start,end (e.g., --range 0,100)")

        ("update-bl", "Update branch lengths of a tree with DIPPER-computed distances without changing topology")

        ("seed", po::value<std::string>(),
        "Random seed for tip selection in --update-bl (default: 0)")

        ("node-dist",           po::value<std::string>(),
        "find distance between a tip into trees, differ by one placement")

        ("max-cluster-size, C",           po::value<std::string>(),
        "cluster size for divide-and-conquer")
        
        ("help,h",
        "Print this help message")
        
        ("version,v", "Print DIPPER version")

        ("force-binary",
         "Force output of fully resolved binary Newick with explicit :0 branches; default is to collapse near-zero "
         "internal edges to multifurcations when writing trees");

    mainDesc.add(requiredDesc).add(optionalDesc);

}

/** Read FASTA sequences into seqs/names, using nameToIdx to match backbone order (e.g. for --add). */
void readAllSequences(po::variables_map& vm, std::vector<std::string>& seqs, std::vector<std::string>& names, std::unordered_map<std::string, int>& nameToIdx)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["input-file"].as<std::string>();

    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    seqs.resize(names.size());

    while (kseq_read(kseq_rd) >= 0) {
        size_t seqLen = kseq_rd->seq.l;
        if (nameToIdx.find(std::string(kseq_rd->name.s, kseq_rd->name.l)) == nameToIdx.end()) {
            seqs.push_back(std::string(kseq_rd->seq.s, seqLen));
            names.push_back(std::string(kseq_rd->name.s, kseq_rd->name.l));
        } else {
            int id = nameToIdx[std::string(kseq_rd->name.s, kseq_rd->name.l)];
            seqs[id] = std::string(kseq_rd->seq.s, seqLen);
        }
    }

    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    // std::cout << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
}

void readSequences(po::variables_map& vm, std::vector<std::string>& seqs, std::vector<std::string>& names, std::unordered_map<std::string, bool>namesInTree={})
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["input-file"].as<std::string>();

    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    while (kseq_read(kseq_rd) >= 0) {
        if (namesInTree.empty() || namesInTree.find(std::string(kseq_rd->name.s, kseq_rd->name.l)) != namesInTree.end()) {
            size_t seqLen = kseq_rd->seq.l;
            seqs.push_back(std::string(kseq_rd->seq.s, seqLen));
            names.push_back(std::string(kseq_rd->name.s, kseq_rd->name.l));
        }
    }
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    // std::cout << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
}




int main(int argc, char** argv) {
    auto inputStart = std::chrono::high_resolution_clock::now();

    parseArguments(argc, argv);

    po::variables_map vm;


    try{
        po::store(po::command_line_parser(argc, argv).options(mainDesc).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        if(vm.count("help")) {
            std::cerr << mainDesc << std::endl;
            return 0;
        }
        else if (vm.count("version")) {
            std::cout << "DIPPER Version " << PROJECT_VERSION << std::endl;
            return 0;
        } 
        std::cerr << "\033[31m" << e.what() << "\033[0m"  << std::endl;
        std::cerr << mainDesc << std::endl;
        return 1;
    }

    MashPlacement::g_printBinaryNewick = vm.count("force-binary") > 0;

    

    
    if (vm.count("add")) {
        if (!vm.count("input-tree")) {
            std::cerr << "\033[31m" << "Backbone tree (--input-tree/-t) is required with --add option" << "\033[0m" << std::endl;
            std::cerr << mainDesc << std::endl;
            return 1;
        }
    }
    

    // Kmer Size
    uint64_t k = 15;
    try {k= (uint64_t)std::stoi(vm["kmer-size"].as<std::string>());}
    catch(std::exception &e){}

    // Sketch Size
    uint64_t sketchSize = 1000;
    // try {sketchSize= (uint64_t)std::stoi(vm["sketch-size"].as<std::string>());}
    // catch(std::exception &e){}

    // Erroneous k-mer thresold
    uint64_t threshold = 1;
    // try {threshold= (uint64_t)std::stoi(vm["threshold"].as<std::string>());}
    // catch(std::exception &e){}

    uint64_t distanceType = 2;
    try {distanceType= (uint64_t)std::stoi(vm["distance-type"].as<std::string>());}
    catch(std::exception &e){}

    std::string in = "r";
    try {in = vm["input-format"].as<std::string>();}
    catch(std::exception &e){}

    std::string out = "t";
    try {out = vm["output-format"].as<std::string>();}
    catch(std::exception &e){}

    std::string algo = "0";
    try {algo = vm["algorithm"].as<std::string>();}
    catch(std::exception &e){}

    std::string placemode = "10";
    try {placemode = vm["K-closest"].as<std::string>();}
    catch(std::exception &e){}

    bool add = false;
    if (vm.count("add")) add = true;

    bool isProtein = false;
    if (vm.count("protein")) isProtein = true;
    bool useReducedProtein = vm.count("protein-reduced") > 0;
    bool useMurphy8 = vm.count("protein-murphy8") > 0;
    if (useReducedProtein && useMurphy8) {
        std::cerr << "ERROR: --protein-reduced and --protein-murphy8 are mutually exclusive.\n";
        return 1;
    }
    if (useReducedProtein) isProtein = true;
    if (useMurphy8) isProtein = true;

    uint64_t clusterSize = 0;
    try {clusterSize = std::stoull(vm["max-cluster-size"].as<std::string>());}
    catch(std::exception &e){}

    std::pair<int, int> range({-1,-1});
    if (vm.count("range")) {
        std::string rangeStr;
        try {
            rangeStr = vm["range"].as<std::string>();
        } catch (std::exception &e) {
            std::cerr << "ERROR: Unable to read --range option: " << e.what() << "\n";
            return 1;
        }
        auto pos = rangeStr.find(',');
        if (pos == std::string::npos) {
            std::cerr << "ERROR: Unable to parse --range option. Expect \"start,end\".\n";
            return 1;
        }
        try {
            range.first = std::stoi(rangeStr.substr(0, pos));
            range.second = std::stoi(rangeStr.substr(pos+1));
        } catch (std::exception &e) {
            std::cerr << "ERROR: Unable to parse --range option. Expect integers: " << e.what() << "\n";
            return 1;
        }
        if (range.first < 0 || range.second < range.first) {
            std::cerr << "ERROR: Invalid --range values. Ensure start>=0 and end>=start.\n";
            return 1;
        }
    }

    std::string treeFile = "";
    try {treeFile = vm["input-tree"].as<std::string>();}
    catch(std::exception &e){}
    if (add && treeFile == "") {
        std::cerr << "ERROR: Input tree file is required for adding query to a backbone tree.\n";
        return 1;
    }
    if (vm.count("update-bl") && treeFile == "") {
        std::cerr << "ERROR: Input tree (--input-tree/-t) is required for --update-bl.\n";
        return 1;
    }

    std::string treeFile2 = "";
    try {treeFile2 = vm["input-tree-2"].as<std::string>();}
    catch(std::exception &e){}

    std::string outputFile = vm["output-file"].as<std::string>();
    std::ofstream output_(outputFile.c_str());

    /* Algorithm selection thresholds: use placement for 30k--1M seqs, divide-and-conquer above. */
    int placement_thr = 30000; 
    int dc_thr = 1000000; 

    MashPlacement::Param params(k, sketchSize, threshold, distanceType, in, out);
    params.range.first = range.first;
    params.range.second = range.second;
    params.isProtein = isProtein;
    params.useReducedProtein = useReducedProtein;
    params.useMurphy8 = useMurphy8;

    if (in == "r" && isProtein && useReducedProtein && k > 32) {
        std::cerr << "ERROR: --protein-reduced supports k-mer size up to 32.\n";
        return 1;
    }
    if (in == "r" && isProtein && useMurphy8 && k > 21) {
        std::cerr << "ERROR: --protein-murphy8 supports k-mer size up to 21.\n";
        return 1;
    }

    bool updateBranchLengths = false;
    bool updateBranchLengths_old = false;
    if (vm.count("update-bl")) updateBranchLengths = true;

    uint64_t seed = 0;
    if (vm.count("seed")) {
        try { seed = std::stoull(vm["seed"].as<std::string>()); }
        catch (std::exception& e) {
            std::cerr << "ERROR: Invalid --seed value.\n";
            return 1;
        }
    }

    std::string node_dist="";
    try {node_dist = vm["node-dist"].as<std::string>();}
    catch(std::exception &e){}

    if (node_dist != "") {
        while (!node_dist.empty() && (node_dist.back() == '\n' || node_dist.back() == '\r' || node_dist.back() == ' ')) {
            node_dist.pop_back();
        }
        while (!node_dist.empty() && (node_dist.front() == ' ' || node_dist.front() == '\t')) {
            node_dist.erase(0, 1);
        }
        if (node_dist.empty()) {
            std::cerr << "ERROR: --node-dist requires a non-empty tip name.\n";
            return 1;
        }
        if (treeFile == "") {
            std::cerr << "ERROR: Input tree (--input-tree/-t) is required for --node-dist option.\n";
            std::cerr << mainDesc << std::endl;
            return 1;
        }
        if (treeFile2 == "") {
            std::cerr << "ERROR: Input tree 2 (--input-tree-2) is required for --node-dist option.\n";
            std::cerr << mainDesc << std::endl;
            return 1;
        }
        std::ifstream treeFileStream(treeFile);
        if (!treeFileStream) {
            std::cerr << "ERROR: Unable to open input tree file: " << treeFile << "\n";
            return 1;
        }
        std::string newickTree;
        std::getline(treeFileStream, newickTree);
        treeFileStream.close();
        UnrootedTree* t1 = new UnrootedTree(newickTree, 0);

        std::ifstream treeFileStream2(treeFile2);
        if (!treeFileStream2) {
            std::cerr << "ERROR: Unable to open input tree file 2: " << treeFile2 << "\n";
            delete t1;
            return 1;
        }
        std::string newickTree2;
        std::getline(treeFileStream2, newickTree2);
        treeFileStream2.close();

        // Align resolution when topologies match but one tree has polytomies: collapse
        // edges in each tree whose non-trivial splits are not displayed in the partner.
        UnrootedTree* t2work = new UnrootedTree(newickTree2, 0);
        std::set<std::string> splitsTree1 = collectNonTrivialSplitKeys(*t1);
        std::set<std::string> splitsTree2 = collectNonTrivialSplitKeys(*t2work);
        t1->collapseToSplitsDisplayedIn(splitsTree2);
        t2work->collapseToSplitsDisplayedIn(splitsTree1);
        std::string newickTree2Aligned = t2work->toNewick();
        delete t2work;

        UnrootedTree* t2 = new UnrootedTree(newickTree2Aligned, t1, node_dist, 0);

        int edgeDist = computeTipPlacementEdgeDistance(t1, t2, node_dist);
        delete t1;
        delete t2;

        if (edgeDist < 0) {
            std::cerr << "ERROR: Could not compute edge distance (tip '" << node_dist << "' not found or invalid trees).\n";
            return 1;
        }
        std::cout << edgeDist << std::endl;
        return 0;
    }

    if (updateBranchLengths_old) {
        std::ifstream treeFileStream(treeFile);
        if (!treeFileStream) {
            std::cerr << "ERROR: Unable to open input tree file: " << treeFile << "\n";
            return 1;
        }
        std::vector<std::string> seqs, names, namesDump;
        readSequences(vm, seqs, namesDump);
        std::cerr << "Read " << seqs.size() << " sequences from input file.\n";
        assert(seqs.size() > 0 && "No sequences found in the input file.");
        std::string newickTree;
        std::getline(treeFileStream, newickTree);
        Tree *t = new Tree(newickTree, namesDump.size());
        std::cerr << "Tree loaded successfully with "<< t->allNodes.size()<<" nodes and root " << t->root->name << ".\n";
        size_t backboneSize = t->m_numLeaves;
        size_t numSequences = seqs.size();
        std::unordered_map<int, int> idMap;
        names.resize(backboneSize);
        for (int i=0; i<(int)numSequences;i++){
            if (t->allNodes.find(namesDump[i]) == t->allNodes.end()) {
                names.push_back(namesDump[i]);
                idMap[i] = names.size()-1;
            } else {
                names[t->allNodes[namesDump[i]]->idx] = namesDump[i];
                idMap[i]=t->allNodes[namesDump[i]]->idx;
            }
        }
        if (in == "m" && out == "t") {
            uint64_t ** fourBitCompressedSeqs = new uint64_t*[numSequences];
            uint64_t * seqLengths = new uint64_t[numSequences];
            bool alignmentLengthModify = false;
            if (params.range.first > 0 || params.range.second > -1) alignmentLengthModify=true;
            tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
            for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
                uint64_t i = static_cast<uint64_t>(idx_);
                int localSeqLength = seqs[i].size();
                if (alignmentLengthModify) {
                    if (params.range.second > -1) localSeqLength=params.range.second+1;
                    if (params.range.first > 0) localSeqLength-=params.range.first;
                }
                uint64_t fourBitCompressedSize = (localSeqLength+15)/16;
                uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
                fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed, params.range.first, params.range.second);
                int newId = idMap[i];
                seqLengths[newId] = localSeqLength;
                fourBitCompressedSeqs[newId] = fourBitCompressed;
            }});
            MashPlacement::msaDeviceArrays.allocateDeviceArrays(fourBitCompressedSeqs, seqLengths, numSequences, params);
            MashPlacement::njDeviceArrays.getDismatrix(
                numSequences,params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays
            );
            std::cout << "Distance matrix construction completed on GPU" << std::endl;
            auto h_mashdist = new double[1ll*numSequences*numSequences];
            cudaMemcpy(h_mashdist, MashPlacement::njDeviceArrays.d_mashDist, sizeof(double)*1ll*numSequences*numSequences, cudaMemcpyDeviceToHost);
            placementAccuracy(t, h_mashdist,numSequences);
            std::string newick= t->getNewickString(t->root);
            output_ << newick;
            output_.close();
            return 0;
        } else {
            std::cerr << "Adding new sequnces only supported with input aligned and unaligned sequences\n";
            exit(1);
        }
        return 0;
    }

    if (updateBranchLengths) {
        std::ifstream treeFileStream(treeFile);
        if (!treeFileStream) {
            std::cerr << "ERROR: Unable to open input tree file: " << treeFile << "\n";
            return 1;
        }
        std::string newickTree;
        std::getline(treeFileStream, newickTree);
        UnrootedTree *t = new UnrootedTree(newickTree);
        std::cerr << "Tree loaded successfully with " << t->numNodes() << " nodes.\n";
        std::unordered_map<std::string, bool> namesInTree;
        for (auto &a: t->getNodes()){
            if (a.second->is_leaf()) {
                namesInTree[a.second->name] = true;
            }
        }

        std::vector<std::string> seqs, namesDump, names;
        readSequences(vm, seqs, namesDump, namesInTree);

        size_t numSequences = seqs.size();
        names.resize(numSequences);
        std::vector<int> ids(numSequences);
        for(int i=0;i<numSequences;i++) ids[i]=i;
        std::mt19937 rnd(time(NULL));
        std::shuffle(ids.begin(),ids.end(),rnd);


        /* print all sequences */
        std::cerr << "Read " << seqs.size() << " sequences from input file.\n";
        assert(seqs.size() > 0 && "No sequences found in the input file.");

        
    

        // if (in == "r" && out == "t") {
        //     uint64_t ** twoBitCompressedSeqs = new uint64_t*[numSequences];
        //     uint64_t * seqLengths = new uint64_t[numSequences];
        //     tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range) {
        //         for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
        //             size_t i = static_cast<size_t>(idx_);
        //             int newId = idMap[i];
        //             uint64_t twoBitCompressedSize = (seqs[i].size() + 31) / 32;
        //             uint64_t * twoBitCompressed = new uint64_t[twoBitCompressedSize];
        //             twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed);
        //             seqLengths[newId] = seqs[i].size();
        //             twoBitCompressedSeqs[newId] = twoBitCompressed;
        //         }
        //     });
        //     MashPlacement::mashDeviceArrays.allocateDeviceArrays(twoBitCompressedSeqs, seqLengths, numSequences, params);
        //     MashPlacement::mashDeviceArrays.sketchConstructionOnGpu(params);
        //     MashPlacement::njDeviceArrays.getDismatrix(
        //         numSequences, params, MashPlacement::mashDeviceArrays,
        //         MashPlacement::matrixReader, MashPlacement::msaDeviceArrays
        //     );
        //     std::cerr << "Distance matrix construction completed on GPU\n";
        // } 
        if (in == "m" && out == "t") {
            uint64_t ** fourBitCompressedSeqs = new uint64_t*[numSequences];
            uint64_t * seqLengths = new uint64_t[numSequences];
            bool alignmentLengthModify = false;
            if (params.range.first > 0 || params.range.second > -1) alignmentLengthModify=true;
            tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
            for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
                uint64_t i = static_cast<uint64_t>(idx_);
                int localSeqLength = seqs[i].size();
                if (alignmentLengthModify) {
                    if (params.range.second > -1) localSeqLength=params.range.second+1;
                    if (params.range.first > 0) localSeqLength-=params.range.first;
                }
                uint64_t fourBitCompressedSize = (localSeqLength+15)/16;
                uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
                fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed, params.range.first, params.range.second);
                
                seqLengths[ids[i]] = localSeqLength;
                fourBitCompressedSeqs[ids[i]] = fourBitCompressed;
                names[ids[i]] = namesDump[i];
            }});
            
            /* Validation 
            // print all edges in the tree
            std::set<int> edgeIds;
            for (auto &a: t->getNodes()){
                for (auto &b: a.second->neighbors) {
                    edgeIds.insert(b.edge_id);
                }
            }

            std::cout << "Number of edges in the tree: " << edgeIds.size() << std::endl;
            for (auto &e: edgeIds) {
                std::cout << e << std::endl;
            }
            */
            /* Add root between first two tips */
            // std::set<UnrootedEdge> edgeIds;
            std::vector<UnrootedNode*> path = t->nodesBetween(names[0], names[1]);
            t->rootBetween(path[path.size()/2-1]->name, path[path.size()/2]->name);

            // Print tree after rooting
            // edgeIds.clear();
            // for (auto &a: t->getNodes()){
            //     for (auto &b: a.second->neighbors) {
            //         std::cout << b.edge_id << " from " << a.second->name << " to " << b.node->name << std::endl;
            //     }
            // }

            // // print root 
            // std::cout << "Root: " << t->m_root->name << std::endl;
            // for (auto &a: t->m_root->neighbors){
            //     std::cout << a.edge_id  << std::endl;
            // }

            std::vector<int> edgeIdsMappingToNewTreeEdges(t->numEdges(),-1);
            MashPlacement::msaDeviceArrays.allocateDeviceArrays(fourBitCompressedSeqs, seqLengths, numSequences, params);
            MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences);
            MashPlacement::kplacementDeviceArrays.estimateBranchLengthsFromTopology(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays, edgeIdsMappingToNewTreeEdges, t, names);
            
            /* update edge lengths in the tree */
            MashPlacement::kplacementDeviceArrays.printTree(names, output_, t, edgeIdsMappingToNewTreeEdges);
        } 

        


        
        
        // else {
        //     std::cerr << "ERROR: --update-bl requires input format 'r' (raw) or 'm' (aligned) with output format 't' (tree).\n";
        //     return 1;
        // }

        // double * h_mashdist = new double[1ll * numSequences * numSequences];
        // cudaMemcpy(h_mashdist, MashPlacement::njDeviceArrays.d_mashDist,
        //            sizeof(double) * 1ll * numSequences * numSequences, cudaMemcpyDeviceToHost);

        // MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences, numSequences);
        // MashPlacement::kplacementDeviceArrays.updateBranchLengthsFromTopology(
        //     t, h_mashdist, numSequences, seed, params,
        //     MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays
        // );

        // std::vector<std::string> placementNames(numSequences);
        // std::vector<int> tip_order = t->getPlacementOrder(seed);
        // for (size_t i = 0; i < numSequences; i++) {
        //     Node* n = t->getNodeByIdx(tip_order[i]);
        //     placementNames[i] = n ? n->name : "";
        // }
        // MashPlacement::kplacementDeviceArrays.printTree(placementNames, output_);
        // output_.close();

        // delete[] h_mashdist;
        // MashPlacement::kplacementDeviceArrays.deallocateDeviceArrays();
        // MashPlacement::njDeviceArrays.deallocateDeviceArrays();
        // if (in == "r") MashPlacement::mashDeviceArrays.deallocateDeviceArrays();
        // else MashPlacement::msaDeviceArrays.deallocateDeviceArrays();

        return 0;
    }

    if (add) {
        /* Load backbone tree from file and build name/index mapping. */
        std::ifstream treeFileStream(treeFile);
        if (!treeFileStream) {
            std::cerr << "ERROR: Unable to open input tree file: " << treeFile << "\n";
            return 1;
        }
        std::vector<std::string> seqs, names, namesDump;
        readSequences(vm, seqs, namesDump);
        std::cerr << "Read " << seqs.size() << " sequences from input file.\n";
        assert(seqs.size() > 0 && "No sequences found in the input file.");

        std::string newickTree;
        std::getline(treeFileStream, newickTree);
        Tree *t = new Tree(newickTree, namesDump.size());
        std::cerr << "Tree loaded successfully with "<< t->allNodes.size()<<" nodes and root " << t->root->name << ".\n";
        size_t backboneSize = t->m_numLeaves;
        size_t numSequences = seqs.size();

        std::unordered_map<int, int> idMap;

        names.resize(backboneSize);
        for (int i=0; i<numSequences;i++){
            if (t->allNodes.find(namesDump[i]) == t->allNodes.end()) {
                names.push_back(namesDump[i]);
                idMap[i] = names.size()-1;
            } else {
                names[t->allNodes[namesDump[i]]->idx] = namesDump[i];
                idMap[i]=t->allNodes[namesDump[i]]->idx;
            }
        }

        /* Unaligned input: 2-bit compress, build MASH sketches, then k-closest placement. */
        if (in == "r" && out == "t") {
            uint64_t ** twoBitCompressedSeqs = new uint64_t*[numSequences];
            uint64_t * seqLengths = new uint64_t[numSequences];
            tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
            for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
                uint64_t i = static_cast<uint64_t>(idx_);
                uint64_t twoBitCompressedSize = (isProtein && !useReducedProtein && !useMurphy8)
                    ? (seqs[i].size()+7)/8
                    : (isProtein && useMurphy8) ? (seqs[i].size()+20)/21
                    : (seqs[i].size()+31)/32;
                uint64_t * twoBitCompressed = new uint64_t[twoBitCompressedSize];
                twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed, isProtein, useReducedProtein, useMurphy8);

                int newId = idMap[i];
                seqLengths[newId] = seqs[i].size();
                twoBitCompressedSeqs[newId] = twoBitCompressed;
            }});
            std::cerr << "Allocating Mash Device Arrays" << std::endl;
            MashPlacement::mashDeviceArrays.allocateDeviceArrays(twoBitCompressedSeqs, seqLengths, numSequences, params);
            
            std::cerr << "Sketch Construction in Progress" << std::endl;
            MashPlacement::mashDeviceArrays.sketchConstructionOnGpu(params);

            MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences, backboneSize);
            MashPlacement::kplacementDeviceArrays.initializeDeviceArrays(t);
            MashPlacement::kplacementDeviceArrays.addQuery(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
            MashPlacement::kplacementDeviceArrays.printTree(names, output_);
        } else if (in == "m" && out == "t") {
            /* Aligned input: 4-bit compress, MSA distance, k-closest placement onto backbone. */
            uint64_t ** fourBitCompressedSeqs = new uint64_t*[numSequences];
            uint64_t * seqLengths = new uint64_t[numSequences];
            bool alignmentLengthModify = false;
            if (params.range.first > 0 || params.range.second > -1) alignmentLengthModify=true;
            tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
            for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
                uint64_t i = static_cast<uint64_t>(idx_);
                int localSeqLength = seqs[i].size();
                if (alignmentLengthModify) {
                    if (params.range.second > -1) localSeqLength=params.range.second+1;
                    if (params.range.first > 0) localSeqLength-=params.range.first;
                }
                uint64_t fourBitCompressedSize = (localSeqLength+15)/16;
                uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
                fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed, params.range.first, params.range.second);
                
                
                int newId = idMap[i];
                seqLengths[newId] = localSeqLength;
                fourBitCompressedSeqs[newId] = fourBitCompressed;
            }});
            MashPlacement::msaDeviceArrays.allocateDeviceArrays(fourBitCompressedSeqs, seqLengths, numSequences, params);
            MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences, backboneSize);
            MashPlacement::kplacementDeviceArrays.initializeDeviceArrays(t);
            MashPlacement::kplacementDeviceArrays.addQuery(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
            MashPlacement::kplacementDeviceArrays.printTree(names, output_);
        } else {
            std::cerr << "Adding new sequnces only supported with input aligned and unaligned sequences\n";
            exit(1);
        }
        return;
    }

    /* Aligned FASTA -> Newick tree. */
    if (in == "m" && out == "t"){
        std::vector<std::string> seqs,names_, names;

        readSequences(vm, seqs, names_);
        size_t numSequences = seqs.size();
        names.resize(numSequences);
        std::vector<int> ids(numSequences);
        for(int i=0;i<numSequences;i++) ids[i]=i;
        std::mt19937 rnd(time(NULL));
        std::shuffle(ids.begin(),ids.end(),rnd);

        /* 4-bit compress aligned sequences (DNA/AA); optionally restrict to --range. */
        auto compressStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "Compressing input sequence using two-bit encoding.\n");
        uint64_t ** fourBitCompressedSeqs = new uint64_t*[numSequences];
        uint64_t * seqLengths = new uint64_t[numSequences];
        bool alignmentLengthModify=false;
        if (params.range.first > 0 || params.range.second > -1) alignmentLengthModify=true;
        tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
        for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
            uint64_t i = static_cast<uint64_t>(idx_);
            
            int localSeqLength = seqs[i].size();
            if (alignmentLengthModify) {
                if (params.range.second > -1) localSeqLength=params.range.second+1;
                if (params.range.first > 0) localSeqLength-=params.range.first;
            }
            uint64_t fourBitCompressedSize = isProtein ? (localSeqLength+7)/8 : (localSeqLength+15)/16;

            uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
            fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed, params.range.first, params.range.second, isProtein);

            seqLengths[ids[i]]=localSeqLength;
            fourBitCompressedSeqs[ids[i]] = fourBitCompressed;
            names[ids[i]] = names_[i];
        }});
        
        auto compressEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds compressTime = compressEnd - compressStart;
        // std::cout << "Compressed in: " <<  compressTime.count() << " ns\n";
        auto inputEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds inputTime = inputEnd - inputStart; 
        std::cerr << "Input in: " <<  inputTime.count()/1000000 << " ms\n";

        /* Allocate MSA device arrays; then run placement, k-closest, DC, or conventional NJ. */
        auto createArrayStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "\nAllocating Gpu device arrays.\n");
        // std::cerr<<"########\n";
        // std::cerr<<"########\n";
        MashPlacement::msaDeviceArrays.allocateDeviceArrays(fourBitCompressedSeqs, seqLengths, numSequences, params);
        if(algo=="1"||algo=="0"&&numSequences>=placement_thr&&numSequences<dc_thr){
            std::cerr<<"Using ";
            if(placemode=="-1"){
                /* Exact placement: consider all branches for each new sequence. */
                std::cerr<<" exact placement mode\n";
                MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
                auto createArrayEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
                std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";


                //Build Tree on Gpu
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::placementDeviceArrays.printTree(names, output_);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

                // Print first 10 hash values corresponding to each sequence
                // MashPlacement::mashDeviceArrays.printSketchValues(10);
                MashPlacement::msaDeviceArrays.deallocateDeviceArrays();
                MashPlacement::placementDeviceArrays.deallocateDeviceArrays();
            }
            else{
                /* K-closest placement: restrict placement to k nearest backbone nodes. */
                std::cerr<<"k-closest placement mode\n";
                MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences);
                auto createArrayEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
                std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";


                //Build Tree on Gpu
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::kplacementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::kplacementDeviceArrays.printTree(names, output_);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

                // Print first 10 hash values corresponding to each sequence
                // MashPlacement::mashDeviceArrays.printSketchValues(10);
                MashPlacement::msaDeviceArrays.deallocateDeviceArrays();
                MashPlacement::kplacementDeviceArrays.deallocateDeviceArrays();
            }
        }
        else if (algo=="3"|| algo=="0"&&numSequences>=dc_thr){
            /* Divide-and-conquer: backbone tree, cluster non-backbone, then place per cluster. */
            std::cerr<<"Using divide-and-conquer mode\n";
            int totalNumSequences = numSequences;
	        int backboneSize = numSequences/100;
            
            if (clusterSize != 0) {
                if (clusterSize < numSequences/1000) {
                    std::cerr << "Warning: cluster size set to too small, setting to default " << numSequences/100 << std::endl;
                } else {
                    backboneSize = clusterSize;
                    std::cerr << "Cluster size set to " << backboneSize << std::endl;
                }
            } else {
                std::cerr << "Cluster size not set, setting to default " << numSequences/100 << std::endl;
            }
            
            params.batchSize = backboneSize;
            params.backboneSize = backboneSize;

            std::vector<int> largeClustersIdx; // cluster idx mapped to closest id's in the backbone tree

            MashPlacement::msaDeviceArraysDC.allocateDeviceArraysDC(fourBitCompressedSeqs, seqLengths, numSequences, params);
            MashPlacement::kplacementDeviceArraysDC.allocateDeviceArraysDC(backboneSize, totalNumSequences);
            auto createArrayEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
            std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";

            //Build Tree on Gpu
            auto createTreeStart = std::chrono::high_resolution_clock::now();
            MashPlacement::kplacementDeviceArraysDC.findBackboneTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            MashPlacement::kplacementDeviceArraysDC.findClustersDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            MashPlacement::kplacementDeviceArraysDC.findClusterTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC, largeClustersIdx);

            auto createTreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
            MashPlacement::kplacementDeviceArraysDC.printTreeDC(names, output_);
            std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

            // Print first 10 hash values corresponding to each sequence
            // MashPlacement::mashDeviceArrays.printSketchValues(10);
            MashPlacement::msaDeviceArraysDC.deallocateDeviceArraysDC();
            MashPlacement::kplacementDeviceArraysDC.deallocateDeviceArraysDC();

            // for (int i=0; i<largeClustersIdx.size(); i++){
            //     std::cout << "Handling cluster " << largeClustersIdx[i] << std::endl;
            //     MashPlacement::kplacementDeviceArraysDC.findBackboneTreeDCRecursive(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC, largeClustersIdx[i]);
            //     break;
            // }

            return 0;
        }
        else{
            /* Conventional neighbor-joining on full distance matrix. */
            std::cerr<<"Using conventional NJ\n";
            if(numSequences>=40000){
                std::cerr<<"Warning: forcing conventional NJ on large datasets might result in unexpected behavior\n";
            }
            MashPlacement::njDeviceArrays.getDismatrix(
                numSequences,params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays
            );
            MashPlacement::njDeviceArrays.findNeighbourJoiningTree(names, output_);
            MashPlacement::msaDeviceArrays.deallocateDeviceArrays();
            MashPlacement::njDeviceArrays.deallocateDeviceArrays();
        }
    }
    /* Unaligned FASTA -> Newick tree. */
    else if (in == "r" && out == "t"){
        std::vector<std::string> seqs,names_, names;

        readSequences(vm, seqs, names_);
        size_t numSequences = seqs.size();
        names.resize(numSequences);
        std::vector<int> ids(numSequences);
        for(int i=0;i<numSequences;i++) ids[i]=i;
        std::mt19937 rnd(time(NULL));
        std::shuffle(ids.begin(),ids.end(),rnd);

        /* 2-bit compress unaligned DNA/RNA for MASH sketching. */
        auto compressStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "Compressing input sequence using two-bit encoding.\n");
        uint64_t ** twoBitCompressedSeqs = new uint64_t*[numSequences];
        uint64_t * seqLengths = new uint64_t[numSequences];
        tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
        for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
            uint64_t i = static_cast<uint64_t>(idx_);
            uint64_t twoBitCompressedSize = (isProtein && !useReducedProtein && !useMurphy8)
                ? (seqs[i].size()+7)/8
                : (isProtein && useMurphy8) ? (seqs[i].size()+20)/21
                : (seqs[i].size()+31)/32;
            uint64_t * twoBitCompressed = new uint64_t[twoBitCompressedSize];
            twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed, isProtein, useReducedProtein, useMurphy8);

            seqLengths[ids[i]] = seqs[i].size();
            twoBitCompressedSeqs[ids[i]] = twoBitCompressed;
            names[ids[i]] = names_[i];
        }});
        
        auto compressEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds compressTime = compressEnd - compressStart;
        // std::cout << "Compressed in: " <<  compressTime.count() << " ns\n";
        auto inputEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds inputTime = inputEnd - inputStart; 
        std::cerr << "Input in: " <<  inputTime.count()/1000000 << " ms\n";

        /* Build tree on GPU: placement (exact or k-closest), DC, or conventional NJ. */
        if(algo=="1"||algo=="0"&&numSequences>=placement_thr&&numSequences<dc_thr){
            auto createArrayStart = std::chrono::high_resolution_clock::now();
            // fprintf(stdout, "\nAllocating Gpu device arrays.\n");
            MashPlacement::mashDeviceArrays.allocateDeviceArrays(twoBitCompressedSeqs, seqLengths, numSequences, params);
            auto createArrayEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
            std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";

            // Build sketch on Gpu
            auto createSketchStart = std::chrono::high_resolution_clock::now();
            MashPlacement::mashDeviceArrays.sketchConstructionOnGpu(params);
            auto createSketchEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createSketchTime = createSketchEnd - createSketchStart; 
            std::cerr << "Sketch Created in: " <<  createSketchTime.count()/1000000 << " ms\n";
            if(placemode=="-1"){
                std::cerr<<"Using exact placement mode\n";
                MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::placementDeviceArrays.printTree(names, output_);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";
                MashPlacement::mashDeviceArrays.deallocateDeviceArrays();
                MashPlacement::placementDeviceArrays.deallocateDeviceArrays();
            }
            else{
                std::cerr<<"Using k-closest placement mode\n";
                MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences);
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::kplacementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::kplacementDeviceArrays.printTree(names, output_);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";
                MashPlacement::mashDeviceArrays.deallocateDeviceArrays();
                MashPlacement::kplacementDeviceArrays.deallocateDeviceArrays();
            }
        }
        else if (algo=="3"||algo=="0"&&numSequences>=dc_thr){
            std::cerr<<"Using divide-and-conquer mode\n";
            
            int totalNumSequences = numSequences;

            int backboneSize = numSequences/100;
            
            if (clusterSize != 0) {
                if (clusterSize < numSequences/1000) {
                    std::cerr << "Warning: cluster size set to too small, setting to default " << numSequences/100 << std::endl;
                } else {
                    backboneSize = clusterSize;
                    std::cerr << "Cluster size set to " << backboneSize << std::endl;
                }
            } else {
                std::cerr << "Cluster size not set, setting to default " << numSequences/100 << std::endl;
            }
            
            if (totalNumSequences < 30000) backboneSize = numSequences/4;

	        params.batchSize = backboneSize;
            params.backboneSize = backboneSize;
            std::vector<int> largeClustersIdx;

            auto createArrayStart = std::chrono::high_resolution_clock::now();
            MashPlacement::mashDeviceArraysDC.allocateDeviceArraysDC(twoBitCompressedSeqs, seqLengths, numSequences, params);
            auto createArrayEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
            std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";

            auto createSketchStart = std::chrono::high_resolution_clock::now();
            MashPlacement::mashDeviceArraysDC.sketchConstructionOnGpuDC(params, twoBitCompressedSeqs, seqLengths, numSequences);
            auto createSketchEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createSketchTime = createSketchEnd - createSketchStart; 
            std::cerr << "Sketch Created in: " <<  createSketchTime.count()/1000000 << " ms\n";
            
            MashPlacement::kplacementDeviceArraysDC.allocateDeviceArraysDC(backboneSize, totalNumSequences);
            MashPlacement::kplacementDeviceArraysHostDC.allocateHostArraysDC(backboneSize, totalNumSequences);
            auto createTreeStart = std::chrono::high_resolution_clock::now();
            
            MashPlacement::kplacementDeviceArraysDC.findBackboneTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            MashPlacement::kplacementDeviceArraysDC.findClustersDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            auto createTreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart;

            MashPlacement::kplacementDeviceArraysDC.findClusterTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC, largeClustersIdx);
            MashPlacement::kplacementDeviceArraysDC.printTreeDC(names, output_);
            MashPlacement::kplacementDeviceArraysDC.deallocateDeviceArraysDC();
            MashPlacement::mashDeviceArraysDC.deallocateDeviceArraysDC();
            std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";
        }
        else{
            std::cerr<<"Using conventional NJ\n";
            if(numSequences>=40000){
                std::cerr<<"Warning: forcing conventional NJ on large datasets might result in unexpected behavior\n";
            }
            // Create arrays
            auto createArrayStart = std::chrono::high_resolution_clock::now();
            MashPlacement::mashDeviceArrays.allocateDeviceArrays(twoBitCompressedSeqs, seqLengths, numSequences, params);
            auto createArrayEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
            std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";

            // Build sketch on Gpu
            auto createSketchStart = std::chrono::high_resolution_clock::now();
            MashPlacement::mashDeviceArrays.sketchConstructionOnGpu(params);
            auto createSketchEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createSketchTime = createSketchEnd - createSketchStart; 
            std::cerr << "Sketch Created in: " <<  createSketchTime.count()/1000000 << " ms\n";

            MashPlacement::njDeviceArrays.getDismatrix(
                numSequences,params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays
            );
            MashPlacement::njDeviceArrays.findNeighbourJoiningTree(names, output_);
            MashPlacement::mashDeviceArrays.deallocateDeviceArrays();
            MashPlacement::njDeviceArrays.deallocateDeviceArrays();
        }

        // Print first 10 hash values corresponding to each sequence
        // MashPlacement::mashDeviceArrays.printSketchValues(10);

    }
    else if(in == "d" && out == "t") {
        /* PHYLIP distance matrix -> Newick tree (placement or NJ; no DC). */
        std::string fileName = vm["input-file"].as<std::string>();
        FILE* filePtr = fopen(fileName.c_str(), "r");
        if (filePtr == nullptr){
            std::cerr << "Cannot open file: " << fileName << std::endl;
            return 1;
        }
        const size_t bufferSize = 64 * 1024 * 1024; 
        char* buffer = new char[bufferSize];
        if (setvbuf(filePtr, buffer, _IOFBF, bufferSize) != 0) {
            std::cerr << "Failed in setting buffer" << std::endl;
            delete[] buffer;
            fclose(filePtr);
            return 1;
        }
        char *temp = new char[20];
        int numSequences;
        fscanf(filePtr, "%d", &numSequences);
        fgets(temp, 20, filePtr);
        MashPlacement::matrixReader.allocateDeviceArrays(numSequences, filePtr);
        if(algo=="1"||algo=="0"&&numSequences>=placement_thr&&numSequences<dc_thr){
            if(placemode=="-1"){
                std::cerr<<"Using exact placement mode\n";
                MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
                MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                MashPlacement::placementDeviceArrays.printTree(MashPlacement::matrixReader.name, output_);
            }
            else{
                std::cerr<<"Using k-closest placement mode\n";
                MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences);
                MashPlacement::kplacementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                MashPlacement::kplacementDeviceArrays.printTree(MashPlacement::matrixReader.name, output_);
            }
        } else if (algo=="3"|| algo=="0"&&numSequences>=dc_thr){
            std::cerr<<"Divide-and-conquer mode not supported with input matrix\n";
            exit(1);
        }
        else{
            std::cerr<<"Using conventional NJ\n";
            if(numSequences>=40000){
                std::cerr<<"Warning: forcing conventional NJ on large datasets might result in unexpected behavior\n";
            }
            MashPlacement::njDeviceArrays.getDismatrix(
                numSequences,params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays
            );
            MashPlacement::njDeviceArrays.findNeighbourJoiningTree(MashPlacement::matrixReader.name, output_);
            MashPlacement::njDeviceArrays.deallocateDeviceArrays();
        }
        fclose(filePtr);
    }
    else{
        printf("Invalid input-output combinations!!!!!\n");
        exit(1);
    }
    return 0;
}
