#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
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
#include <filesystem>

#include <cstring>   
#include <unistd.h>  
#include <sys/stat.h> 


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
po::options_description mergeDesc("DIPPER Merge Trees Command Line Arguments");

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

        ("batch-size,B",        po::value<std::string>(),
        "Batch size:\n (default: 100000)")

        ("add,a",
        "Add query to backbone using k-closest placement")

        ("input-tree,t",       po::value<std::string>(),
        "Input backbone tree (Newick format), required with --add option")

        ("shuffle,S",       po::value<std::string>(),
        "Shuffle input sequences")

        ("cluster-archieve,C",       po::value<std::string>(),
        "Cluster archieve directory path")

        ("rename,R",       
        "Cluster archieve directory path")

        ("protein,p",
        "Input sequences are protein sequences (default: DNA/RNA sequences)")

        ("help,h",
        "Print this help message")
        
        ("version,v", "Print DIPPER version")

        ("print-binary-newick",
         "Print fully resolved binary Newick with explicit :0 branches; default is to collapse near-zero "
         "internal edges to multifurcations when writing trees");

    mainDesc.add(requiredDesc).add(optionalDesc);

}

void parseArgumentsMerge(int argc, char** argv)
{
    // Setup boost::program_options
    po::options_description requiredDesc("Required Options");
    requiredDesc.add_options()
        ("trees,t",     po::value<std::vector<std::string>>()->multitoken()->required(),
        //  po::value<std::vector<std::string>>()->required(),
        "Input trees in Newick format")

        ("output-file,O",      po::value<std::string>()->required(),
        "Output file path");

    po::options_description optionalDesc("Optional Options");
    optionalDesc.add_options()
        ("help,h",
        "Print this help message")
        
        ("version,v", "Print DIPPER version");
    mergeDesc.add(requiredDesc).add(optionalDesc);

}

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
    // std::cerr << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
}

void readSequences(po::variables_map& vm, std::vector<std::string>& seqs, std::vector<std::string>& names)
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
        size_t seqLen = kseq_rd->seq.l;
        seqs.push_back(std::string(kseq_rd->seq.s, seqLen));
        names.push_back(std::string(kseq_rd->name.s, kseq_rd->name.l));
    }
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    // std::cerr << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
}

struct GzFastaReader {
    gzFile gz = nullptr;
    kseq_t* ks = nullptr;
    bool eof = false;
};

static GzFastaReader* open_gz_fasta_reader(const std::string &path) {
    gzFile f = gzopen(path.c_str(), "r");
    if (!f) {
        std::cerr << "ERROR: cant open file: " << path << "\n";
        return nullptr;
    }
    kseq_t* k = kseq_init(f);
    if (!k) {
        gzclose(f);
        std::cerr << "ERROR: kseq_init failed for: " << path << "\n";
        return nullptr;
    }
    auto *r = new GzFastaReader();
    r->gz = f;
    r->ks = k;
    r->eof = false;
    return r;
}

static int read_next_batch(GzFastaReader* reader, int count, std::vector<std::string>& seqs, std::vector<std::string>& names) {
    if (!reader || !reader->ks || reader->eof) return 0;
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    int read = 0;
    while (read < count && kseq_read(reader->ks) >= 0) {
        size_t seqLen = reader->ks->seq.l;
        seqs.emplace_back(reader->ks->seq.s, seqLen);
        names.emplace_back(reader->ks->name.s, reader->ks->name.l);
        ++read;
    }
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cerr << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
    return read;
}

static void close_gz_fasta_reader(GzFastaReader* reader) {
    if (!reader) return;
    if (reader->ks) {
        kseq_destroy(reader->ks);
        reader->ks = nullptr;
    }
    if (reader->gz) {
        gzclose(reader->gz);
        reader->gz = nullptr;
    }
    delete reader;
}

void readSequences_batch(po::variables_map& vm, int start, int count, std::vector<std::string>& seqs, std::vector<std::string>& names)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["input-file"].as<std::string>();

    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    int idx = 0;
    while (idx < start && kseq_read(kseq_rd) >= 0) {
        idx++;
    }

    while (kseq_read(kseq_rd) >= 0 && count > 0) {
        size_t seqLen = kseq_rd->seq.l;
        seqs.push_back(std::string(kseq_rd->seq.s, seqLen));
        names.push_back(std::string(kseq_rd->name.s, kseq_rd->name.l));
        count--;
    }
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    
}

std::string make_temp_dir_mkdtemp() {
    char tmpl[] = "tmp_XXXXXX";
    char* ret = mkdtemp(tmpl);
    if (!ret) {
        std::cerr << "mkdtemp failed: " << strerror(errno) << "\n";
        return {};
    }
    return std::string(ret);
}

bool append_to_gzip(const std::string &gzip_path, const char *buf, size_t len) {
    gzFile gz = gzopen(gzip_path.c_str(), "ab");
    if (!gz) {
        std::cerr << "gzopen failed: " << strerror(errno) << "\n";
        return false;
    }

    int written = gzwrite(gz, buf, static_cast<unsigned int>(len));
    if (written == 0) {
        int errnum = 0;
        const char *errstr = gzerror(gz, &errnum);
        std::cerr << "gzwrite failed: " << (errstr ? errstr : "unknown")
                  << " (errnum=" << errnum << ")\n";
        gzclose(gz);
        return false;
    }

    if (gzclose(gz) != Z_OK) {
        std::cerr << "gzclose returned error\n";
        return false;
    }
    return true;
}

bool append_to_gzip_without_open_close(gzFile gz, const std::string &gzip_path, const char *buf, size_t len) {

    int written = gzwrite(gz, buf, static_cast<unsigned int>(len));
    if (written == 0) {
        int errnum = 0;
        const char *errstr = gzerror(gz, &errnum);
        std::cerr << "gzwrite failed: " << (errstr ? errstr : "unknown")
                  << " (errnum=" << errnum << ")\n";
        // gzclose(gz); // DO NOT close the handle here.
        return false;
    }
    return true;

}

bool read_line_local(gzFile file, std::string& line, int maxLength) {
    char buffer[maxLength];
    if (gzgets(file, buffer, maxLength) == Z_NULL) {
        return false;
    }
    line = buffer;
    if (!line.empty() && line.back() == '\n') line.pop_back();
    return true;
}

int main(int argc, char** argv) {
    auto inputStart = std::chrono::high_resolution_clock::now();

    if (std::string(argv[1]) == "merge"){
        parseArgumentsMerge(argc, argv);
        po::variables_map vmMerge;

        try{
            po::store(po::command_line_parser(argc, argv).options(mergeDesc).run(), vmMerge);
            po::notify(vmMerge);
        }
        catch(std::exception &e){
            if(vmMerge.count("help")) {
                std::cerr << mergeDesc << std::endl;
                return 0;
            }
            else if (vmMerge.count("version")) {
                std::cerr << "DIPPER Version " << PROJECT_VERSION << std::endl;
                return 0;
            }
            std::cerr << "\033[31m" << e.what() << "\033[0m"  << std::endl;
            std::cerr << mergeDesc << std::endl;
            return 1;
        }

        std::vector<std::string> treeFiles = vmMerge["trees"].as<std::vector<std::string>>();
        std::string outputFile = vmMerge["output-file"].as<std::string>();

        
        return 0;
    }

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
            std::cerr << "DIPPER Version " << PROJECT_VERSION << std::endl;
            return;
        } 
        std::cerr << "\033[31m" << e.what() << "\033[0m"  << std::endl;
        std::cerr << mainDesc << std::endl;
        return 1;
    }

    MashPlacement::g_printBinaryNewick = vm.count("print-binary-newick") > 0;

    
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
    try {sketchSize= (uint64_t)std::stoi(vm["sketch-size"].as<std::string>());}
    catch(std::exception &e){}

    // Erroneous k-mer thresold
    uint64_t threshold = 1;
    // try {threshold= (uint64_t)std::stoi(vm["threshold"].as<std::string>());}
    // catch(std::exception &e){}

    uint64_t distanceType = 1;
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

    int batchsize = 100000;
    try {batchsize = std::stoi(vm["batch-size"].as<std::string>());}
    catch(std::exception &e){}

    bool add = false;
    if (vm.count("add")) add = true;

    bool shuffle = false;
    if (vm.count("shuffle")) shuffle = true;

    bool rename = false;
    if (vm.count("rename")) rename = true;

    bool isProtein = false;
    if (vm.count("protein")) isProtein = true;

    std::string inClusterArchieve = "";
    try {inClusterArchieve = vm["cluster-archieve"].as<std::string>();}
    catch(std::exception &e){}

    std::string treeFile = "";
    try {treeFile = vm["input-tree"].as<std::string>();}
    catch(std::exception &e){}
    if (add && treeFile == "") {
        std::cerr << "ERROR: Input tree file is required for adding query to a backbone tree.\n";
        return 1;
    }

    std::string outputFile = vm["output-file"].as<std::string>();
    std::ofstream output_(outputFile.c_str());

    int placement_thr = 30000; 
    int dc_thr = 1000000; 
    bool optimized = false;
    bool merge_cluster = true;

    MashPlacement::Param params(k, sketchSize, threshold, distanceType, in, out);
    params.isProtein = isProtein;

    if (rename) {
        // Load the tree from the file
        std::ifstream treeFileStream(treeFile);
        if (!treeFileStream) {
            std::cerr << "ERROR: Unable to open input tree file: " << treeFile << "\n";
            return 1;
        }

        std::string newickTree;
        std::getline(treeFileStream, newickTree);
        Tree *t = new Tree(newickTree);
        std::cerr << "Tree loaded successfully with "<< t->allNodes.size()<<" nodes and root " << t->root->name << ".\n";

        t->updateNewick();
        std::string out=t->getNewickString(t->root);
        return;

    }

    if (shuffle) {
        std::cerr << "Shuffling input sequences" << std::endl;
        std::vector<std::string> seqs,names;

        readSequences(vm, seqs, names);
        size_t numSequences = seqs.size();
        std::vector<int> ids(numSequences);
        for(int i=0;i<numSequences;i++) ids[i]=i;
        std::mt19937 rnd(time(NULL));
        std::shuffle(ids.begin(),ids.end(),rnd);
        
        std::ofstream outFile(vm["output-file"].as<std::string>());
        for (int idx=0; idx<seqs.size(); ++idx) {
            int printIdx = ids[idx];
            outFile << ">" << names[printIdx] << "\n" << seqs[printIdx] << "\n";
        };
        outFile.close();

        std::cerr << "Sequences shuffled" << std::endl;
        return 0;
    }

    if (add) {
        // Load the tree from the file
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

        if (in == "r" && out == "t") {
            uint64_t ** twoBitCompressedSeqs = new uint64_t*[numSequences];
            uint64_t * seqLengths = new uint64_t[numSequences];
            tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
            for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
                uint64_t i = static_cast<uint64_t>(idx_);
                uint64_t twoBitCompressedSize = isProtein ? (seqs[i].size()+7)/8 : (seqs[i].size()+31)/32;
                uint64_t * twoBitCompressed = new uint64_t[twoBitCompressedSize];
                twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed, isProtein);

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
            uint64_t ** fourBitCompressedSeqs = new uint64_t*[numSequences];
            uint64_t * seqLengths = new uint64_t[numSequences];
            tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
            for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
                uint64_t i = static_cast<uint64_t>(idx_);
                uint64_t fourBitCompressedSize = (seqs[i].size()+15)/16;
                uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
                fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed);

                int newId = idMap[i];
                seqLengths[newId] = seqs[i].size();
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

    if (in == "m" && out == "t" && algo == "3" && optimized){
        std::cerr<<"Using divide-and-conquer Optmized mode\n";

        /* Muli-GPU support available here only */
        int maxGpuNum;
        cudaError_t err = cudaGetDeviceCount(&maxGpuNum);
        if (err != cudaSuccess) {
            std::cerr << "No GPU detected." << std::endl;
            maxGpuNum = 0; // or handle accordingly
            return 0;
        }
        // int gpuNum = (vm.count("gpu")) ? vm["gpu"].as<int>() : maxGpuNum;
        int gpuNum = maxGpuNum;
        std::cerr << "Using " << gpuNum << " GPUs" << std::endl;
        tbb::global_control init(tbb::global_control::max_allowed_parallelism, gpuNum);
        
        std::string seqFileName = vm["input-file"].as<std::string>();
        gzFile f_rd = gzopen(seqFileName.c_str(), "r");
        if (!f_rd) {
            fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
            exit(1);
        }

        kseq_t* kseq_rd = kseq_init(f_rd);

        int totalSequences = 0;
        std::unordered_map<std::string, int> globalNameToIdx;
        while (kseq_read(kseq_rd) >= 0) {
            globalNameToIdx[std::string(kseq_rd->name.s, kseq_rd->name.l)] = totalSequences;
            totalSequences++;
        }
        kseq_destroy(kseq_rd);
        gzclose(f_rd);

        GzFastaReader* batchReader = open_gz_fasta_reader(seqFileName);
        if (!batchReader) exit(1);

        // Read Batched input for backbone tree construction
        std::vector<std::string> seqs, names, globalNames;
        globalNames.resize(totalSequences);
        int got = read_next_batch(batchReader, batchsize, seqs, names);
        size_t numSequences = seqs.size();

        
        
        // Compress Sequences (2-bit compressor)
        auto compressStart = std::chrono::high_resolution_clock::now();
        
        uint64_t ** fourBitCompressedSeqs = new uint64_t*[numSequences];
        uint64_t * seqLengths = new uint64_t[numSequences];
        size_t maxLengthCompressed = 0;

        tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
        for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
            uint64_t i = static_cast<uint64_t>(idx_);
            uint64_t fourBitCompressedSize = (seqs[i].size()+15)/16;
            uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
            fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed);

            seqLengths[i] = seqs[i].size();
            fourBitCompressedSeqs[i] = fourBitCompressed;
            globalNames[i] = names[i];
        }});

        auto compressEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds compressTime = compressEnd - compressStart;
        auto inputEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds inputTime = inputEnd - inputStart; 

        // Create arrays
        auto createArrayStart = std::chrono::high_resolution_clock::now();
  
        int totalNumSequences = totalSequences;
        int backboneSize = batchsize;
        params.totalNumSeqs = totalNumSequences;
        params.batchSize = backboneSize;
        params.backboneSize = backboneSize;
        
        int gn;
        std::vector<MashPlacement::MSADeviceArraysDC> msaDeviceArraysDCs(gpuNum);
        std::vector<MashPlacement::KPlacementDeviceArraysDC> kplacementDeviceArraysDCs(gpuNum);
        tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range) { 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            cudaSetDevice(gn);
            MashPlacement::MSADeviceArraysDC msaDeviceArraysDC;
            MashPlacement::KPlacementDeviceArraysDC kplacementDeviceArraysDC;
            msaDeviceArraysDCs[gn] = msaDeviceArraysDC;
            kplacementDeviceArraysDCs[gn] = kplacementDeviceArraysDC;
        }});

        tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range) { 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            cudaSetDevice(gn);
            msaDeviceArraysDCs[gn].allocateDeviceArraysDC(fourBitCompressedSeqs, seqLengths, numSequences, params);
            kplacementDeviceArraysDCs[gn].allocateDeviceArraysDC(backboneSize, totalNumSequences);
            auto createArrayEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 

            auto createTreeStart = std::chrono::high_resolution_clock::now();
            kplacementDeviceArraysDCs[gn].findBackboneTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, msaDeviceArraysDCs[gn], MashPlacement::kplacementDeviceArraysHostDC);
            auto createTreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
            std::cerr << "Backbone tree created in: " <<  createTreeTime.count()/1000000 << " ms" << std::endl;
        }});

        maxLengthCompressed=(seqLengths[0]+15)/16;

        
        /* Delete pointers */
        for (uint64_t i=0; i<numSequences; i++){
            delete[] fourBitCompressedSeqs[i];
        }
        delete[] fourBitCompressedSeqs;
        delete[] seqLengths;

        // Clustering
        int clusteringSeqIdx = 0; // seq idx within the batch
        int numClusteres = totalNumSequences/backboneSize + (totalNumSequences%backboneSize!=0);
        MashPlacement::kplacementDeviceArraysHostDC.clusterID = new int[totalNumSequences];
        
        if(merge_cluster){
            cudaSetDevice(1);
            std::vector<bool> isCluster(backboneSize*4-4, false);

            /* iterate over all files in the archieve */
            int test_var = 0;
            int actual_placement = 0;
            int clustersPerBatchFile = 1000;
            int clusterFiles = (numSequences*4-4)/clustersPerBatchFile + ((numSequences*4-4)%clustersPerBatchFile != 0);
            for (int bc=0; bc<clusterFiles;bc++){
                std::string file = inClusterArchieve + "/" + std::to_string(bc) + ".gz";
                /* if file exist */
                gzFile gzfile = gzopen(file.c_str(), "rb");
                if (gzfile != Z_NULL){
                    uint64_t counter=0;
                    std::string line;
                    std::string name, id, cluster;
                    while (read_line_local(gzfile, line, maxLengthCompressed)) {
                        if (line.empty() || line[0] != '>') continue;

                        std::istringstream ss(line.substr(1));
                        std::getline(ss, name, '\t');
                        std::getline(ss, cluster, '\t');
                        std::getline(ss, id, '\t');
                        MashPlacement::kplacementDeviceArraysHostDC.clusterID[globalNameToIdx[name]] = std::stoi(cluster);
                        if (std::stoi(cluster) < 0 or std::stoi(cluster) > backboneSize*4-4){
                            std::cerr << "Error in cluster archieve file: " << file << " for " << name << "\t" << cluster << std::endl;
                        }
                        isCluster[std::stoi(cluster)] = true;
                        counter++;
                    }
                }
                gzclose(gzfile);
            }
            kplacementDeviceArraysDCs[1].findClusterTreeDC_batch(params, 
                                                                MashPlacement::mashDeviceArraysDC, 
                                                                MashPlacement::matrixReader, 
                                                                msaDeviceArraysDCs[1], 
                                                                MashPlacement::kplacementDeviceArraysHostDC, 
                                                                inClusterArchieve, 
                                                                isCluster);
            return;
        }
        
        std::cerr << "Total clusters to be processed: " << numClusteres << std::endl;
        

        // prepare directory and files
        int clusterPerBatchFile = 1000;
        int clusterFiles = (batchsize*4-4)/clusterPerBatchFile + ((batchsize*4-4)%clusterPerBatchFile!=0);
        auto dir = make_temp_dir_mkdtemp();
        std::cerr << "Created temporary directory: " << dir << std::endl;
        for (int i=0; i<clusterFiles; i++){
            std::string filename = dir + "/" + std::to_string(i) + ".gz";
            std::ofstream ofs(filename, std::ios::binary);
            ofs.close();
        }

        int cpuBatchCount = 6;

        uint64_t ** fourBitCompressedSeqsCluster = new uint64_t*[batchsize * cpuBatchCount];
        uint64_t * seqLengthsCluster = new uint64_t[batchsize * cpuBatchCount];
        std::vector<bool> isCluster(backboneSize*4-4, false);
        
        int tmp_totalSeq_write = 0;
        for (int clusteringBatchIdx=1; clusteringBatchIdx < numClusteres; clusteringBatchIdx+=cpuBatchCount){
            std::cerr << "Processing clusters " << clusteringBatchIdx << " to " << (cpuBatchCount+clusteringBatchIdx>=numClusteres ? numClusteres-1 : clusteringBatchIdx+cpuBatchCount-1) << std::endl;
            int globalSeqIdx = clusteringBatchIdx*backboneSize;
            // std::vector<std::string> seqs, names;
            seqs.clear();
            names.clear();
            // readSequences_batch(vm, globalSeqIdx, backboneSize, seqs, names);
            int got = read_next_batch(batchReader, batchsize*cpuBatchCount, seqs, names);
            size_t cpuBatchSize_local = seqs.size();
            std::cerr << "Working on " << cpuBatchSize_local << " seqs" << std::endl;
            assert(cpuBatchSize_local > 0);
            
            auto seqReadStartG = std::chrono::high_resolution_clock::now();
            auto seqReadStart = std::chrono::high_resolution_clock::now();
            tbb::parallel_for(tbb::blocked_range<int>(0, cpuBatchSize_local), [&](tbb::blocked_range<int> range){
            for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
                uint64_t i = static_cast<uint64_t>(idx_);
                uint64_t fourBitCompressedSize = (seqs[i].size()+15)/16;
                uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
                fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed);

                seqLengthsCluster[i] = seqs[i].size();
                fourBitCompressedSeqsCluster[i] = fourBitCompressed;
                globalNames[globalSeqIdx+i] = names[i];
            }});
            auto seqReadEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
            std::cerr << "Compressed in: " <<  seqReadTime.count() << " ns" << std::endl;

            seqReadStart = std::chrono::high_resolution_clock::now();
            // for (int gpuCluster=0; gpuCluster<cpuBatchCount; gpuCluster++){
            std::atomic<int> gpuCluster;
            gpuCluster.store(0);
            
            tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){
            for (int gn = range.begin(); gn < range.end(); ++gn) {
                cudaSetDevice(gn);
                while (gpuCluster < cpuBatchCount) {
                    int gC = gpuCluster.fetch_add(1);
                    if (clusteringBatchIdx+gC >= numClusteres) break;
                    msaDeviceArraysDCs[gn].transferToDeviceArraysDC(fourBitCompressedSeqsCluster, seqLengthsCluster, cpuBatchSize_local, gC, params);
                    // MashPlacement::msaDeviceArraysDC.transferToDeviceArraysDC(fourBitCompressedSeqsCluster, seqLengthsCluster, cpuBatchSize_local, gpuCluster, params);
                    int gIdx = clusteringBatchIdx+gC;
                    kplacementDeviceArraysDCs[gn].findClustersDC_batch(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, msaDeviceArraysDCs[gn], MashPlacement::kplacementDeviceArraysHostDC, gIdx);
                    // MashPlacement::kplacementDeviceArraysDC.findClustersDC_batch(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC, gIdx);
                }            
            }});    

        
            seqReadEnd = std::chrono::high_resolution_clock::now();
            seqReadTime = seqReadEnd - seqReadStart;
            std::cerr << "Clustered in: " <<  seqReadTime.count() << " ns" << std::endl;

            std::string cudaError = cudaGetErrorString(cudaGetLastError());
            if (cudaError != "no error")
                std::cerr << "CUDA Error after clustering: " << cudaError << std::endl;
            
            // Writing seqs to cluster files
            seqReadStart = std::chrono::high_resolution_clock::now();
            uint64_t fourBitCompressedSize = (seqs[0].size()+15)/16;
            for (int bc=0; bc<clusterFiles; bc++){
                
                std::string filename = dir + "/" + std::to_string(bc) + ".gz";

                gzFile gz = gzopen(filename.c_str(), "ab");
                if (!gz) {
                    std::cerr << "gzopen failed: " << strerror(errno) << "\n";
                    return false;
                }

                for (int f=0;f<cpuBatchSize_local;f++){
                    int cluster = MashPlacement::kplacementDeviceArraysHostDC.clusterID[globalSeqIdx + f];
                    int batchFile = (cluster/clusterPerBatchFile); 
                    if (batchFile == bc){
                        isCluster[cluster] = true;
                        tmp_totalSeq_write++;
                        std::string header = ">"+names[f]+"\t"
                                                +std::to_string(cluster)+"\t"
                                                +std::to_string(globalSeqIdx+f)+ "\n";
                        append_to_gzip_without_open_close(gz, filename, header.c_str(), header.size());
                        size_t byteCount = static_cast<size_t>(fourBitCompressedSize) * sizeof(uint64_t);
                        if (byteCount > 0 && fourBitCompressedSeqsCluster[f] != nullptr) {
                            const char* data = reinterpret_cast<const char*>(fourBitCompressedSeqsCluster[f]);
                            if (!append_to_gzip_without_open_close(gz, filename, data, byteCount)) {
                                std::cerr << "Failed to append compressed bytes to " << filename << "\n";
                            }
                        }

                        // terminate entry
                        const char nl = '\n';
                        append_to_gzip_without_open_close(gz, filename, &nl, 1);
                    }
                }
                if (gzclose(gz) != Z_OK) {
                    std::cerr << "gzclose returned error\n";
                    return false;
                }
            }
            seqReadEnd = std::chrono::high_resolution_clock::now();
            seqReadTime = seqReadEnd - seqReadStart;
            std::cerr << "Written to files in: " <<  seqReadTime.count() << " ns" << std::endl;

            auto seqReadEndG = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds seqReadTimeG = seqReadEndG - seqReadStartG;
            std::cerr << "Clusters processed in: " <<  seqReadTimeG.count() << " ns" << std::endl;

            for (uint64_t i=0; i<cpuBatchSize_local; i++){
                delete[] fourBitCompressedSeqsCluster[i];
            }
        
            seqReadEnd = std::chrono::high_resolution_clock::now();
            seqReadTime = seqReadEnd - seqReadStart;
            
        }
        std::cerr << "Total sequences written to cluster files: " << tmp_totalSeq_write << std::endl;
        delete[] fourBitCompressedSeqsCluster;
        delete[] seqLengthsCluster;
        close_gz_fasta_reader(batchReader);
    

        auto seqReadStartG = std::chrono::high_resolution_clock::now();
        cudaSetDevice(1);
        kplacementDeviceArraysDCs[1].findClusterTreeDC_batch(params, 
                                                            MashPlacement::mashDeviceArraysDC, 
                                                            MashPlacement::matrixReader, 
                                                            msaDeviceArraysDCs[1], 
                                                            MashPlacement::kplacementDeviceArraysHostDC, 
                                                            dir, 
                                                            isCluster);
        
        auto seqReadEndG = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds seqReadTimeG = seqReadEndG - seqReadStartG;
        std::cerr << "Restricted placedment done in: " <<  seqReadTimeG.count() << " ns" << std::endl;

        kplacementDeviceArraysDCs[1].printTreeDC(globalNames, output_);
        msaDeviceArraysDCs[1].deallocateDeviceArraysDC();
        kplacementDeviceArraysDCs[1].deallocateDeviceArraysDC();

        // std::cerr << "Removing files" << std::endl;
        // std::error_code ec;
        // bool removed = std::filesystem::remove_all(dir, ec);
        // if (!ec && removed) {
        //     std::cerr << "Removed: " << dir << "\n";
        // } else if (ec) {
        //     std::cerr << "Failed to remove " << dir << ": " << ec.message() << "\n";
        //     return 1;
        // } else {
        //     std::cerr << dir << " did not exist or was not removable (maybe not empty)\n";
        // }
        return 0;
    } else if (in == "m" && out == "t"){
        std::vector<std::string> seqs,names_, names;

        // Read Input Sequences (Fasta format)
        readSequences(vm, seqs, names_);
        size_t numSequences = seqs.size();
        names.resize(numSequences);
        std::vector<int> ids(numSequences);
        for(int i=0;i<numSequences;i++) ids[i]=i;
        std::mt19937 rnd(time(NULL));
        std::shuffle(ids.begin(),ids.end(),rnd);

    
        // Compress Sequences (2-bit compressor)
        auto compressStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "Compressing input sequence using two-bit encoding.\n");
        uint64_t ** fourBitCompressedSeqs = new uint64_t*[numSequences];
        uint64_t * seqLengths = new uint64_t[numSequences];
        tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
        for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
            uint64_t i = static_cast<uint64_t>(idx_);
            uint64_t fourBitCompressedSize = (seqs[i].size()+15)/16;
            uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
            fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed);

            seqLengths[ids[i]] = seqs[i].size();
            fourBitCompressedSeqs[ids[i]] = fourBitCompressed;
            names[ids[i]] = names_[i];
        }});
        
        auto compressEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds compressTime = compressEnd - compressStart;
        // std::cerr << "Compressed in: " <<  compressTime.count() << " ns\n";
        auto inputEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds inputTime = inputEnd - inputStart; 
        std::cerr << "Input in: " <<  inputTime.count()/1000000 << " ms\n";

        // Create arrays
        auto createArrayStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "\nAllocating Gpu device arrays.\n");
        // std::cerr<<"########\n";
        // std::cerr<<"########\n";
        
        if(algo=="1"||algo=="0"&&numSequences>=placement_thr&&numSequences<dc_thr){
            MashPlacement::msaDeviceArrays.allocateDeviceArrays(fourBitCompressedSeqs, seqLengths, numSequences, params);
            std::cerr<<"Using ";
            if(placemode=="-1"){
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
            std::cerr<<"Using divide-and-conquer mode\n";
            int totalNumSequences = numSequences;
            int backboneSize = numSequences/20;
            params.batchSize = backboneSize;
            params.backboneSize = backboneSize;
            MashPlacement::msaDeviceArraysDC.allocateDeviceArraysDC(fourBitCompressedSeqs, seqLengths, numSequences, params);
            MashPlacement::kplacementDeviceArraysDC.allocateDeviceArraysDC(backboneSize, totalNumSequences);
            auto createArrayEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
            std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";

            //Build Tree on Gpu
            auto createTreeStart = std::chrono::high_resolution_clock::now();
            MashPlacement::kplacementDeviceArraysDC.findBackboneTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            MashPlacement::kplacementDeviceArraysDC.findClustersDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            MashPlacement::kplacementDeviceArraysDC.findClusterTreeDCRecursive(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);

            auto createTreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
            MashPlacement::kplacementDeviceArraysDC.printTreeDC(names, output_);
            std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

            // Print first 10 hash values corresponding to each sequence
            // MashPlacement::mashDeviceArrays.printSketchValues(10);
            MashPlacement::msaDeviceArraysDC.deallocateDeviceArraysDC();
            MashPlacement::kplacementDeviceArraysDC.deallocateDeviceArraysDC();
        }
        else{
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
    else if (in == "r" && out == "t"){
        std::vector<std::string> seqs,names_, names;

        // Read Input Sequences (Fasta format)
        readSequences(vm, seqs, names_);
        size_t numSequences = seqs.size();
        names.resize(numSequences);
        std::vector<int> ids(numSequences);
        for(int i=0;i<numSequences;i++) ids[i]=i;
        std::mt19937 rnd(time(NULL));
        std::shuffle(ids.begin(),ids.end(),rnd);
        
        // Compress Sequences (2-bit compressor)
        auto compressStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "Compressing input sequence using two-bit encoding.\n");
        uint64_t ** twoBitCompressedSeqs = new uint64_t*[numSequences];
        uint64_t * seqLengths = new uint64_t[numSequences];
        tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
        for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
            uint64_t i = static_cast<uint64_t>(idx_);
            uint64_t twoBitCompressedSize = isProtein ? (seqs[i].size()+7)/8 : (seqs[i].size()+31)/32;
            uint64_t * twoBitCompressed = new uint64_t[twoBitCompressedSize];
            twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed, isProtein);

            seqLengths[ids[i]] = seqs[i].size();
            twoBitCompressedSeqs[ids[i]] = twoBitCompressed;
            names[ids[i]] = names_[i];
        }});
        
        auto compressEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds compressTime = compressEnd - compressStart;
        // std::cerr << "Compressed in: " <<  compressTime.count() << " ns\n";
        auto inputEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds inputTime = inputEnd - inputStart; 
        std::cerr << "Input in: " <<  inputTime.count()/1000000 << " ms\n";


        //Build Tree on Gpu
        if(algo=="1"||algo=="0"&&numSequences>=placement_thr&&numSequences<dc_thr){
            // Create arrays
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
            params.batchSize = backboneSize;
            params.backboneSize = backboneSize;

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

            MashPlacement::kplacementDeviceArraysDC.findClusterTreeDCRecursive(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
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
    else {
        printf("Invalid input-output combinations!!!!!\n");
        exit(1);
    }
    return 0;
}
     
