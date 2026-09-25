#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <boost/program_options.hpp> 
#include "../../src/kseq.h"
#include "zlib.h"



#ifndef TWOBITCOMPRESSOR_HPP
#include "../../src/twoBitCompressor.hpp"
#endif

#ifndef FOURBITCOMPRESSOR_HPP
#include "../../src/fourBitCompressor.hpp"
#endif

#ifndef MASHPL_CUH
#include "../../src/divide_and_conquer/mash_placement.cuh"
#endif

namespace po = boost::program_options;

KSEQ_INIT2(, gzFile, gzread)

po::options_description mainDesc("Placement Command Line Arguments");


void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        // ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
        ("input-file,f", po::value<std::string>()->required(), "Input format, phylip format for distance matrix, fasta format for MSA or raw sequences, required")
        ("kmer-size,k", po::value<std::string>(), "Kmer-size (Valid: 2-15, Default = 15)")
        ("sketch-size,s", po::value<std::string>(), "Sketch-size (Default = 10000)")
        ("threshold,r", po::value<std::string>(), "Erroneous k-mer threshold (Default = 1)")
        ("distance-type,t", po::value<std::string>(), "Distance type to be calculated, range from 1 to 6: 1 - uncorrected, 2 - JC, 3 - tajima-nei, 4 - K2P, 5 - Tamura, 6 - Jinnei")
        ("input-format,i", po::value<std::string>()->required(), "Input format (d - distance matrix, r - raw sequences, m - MSA), required")
        ("output-format,o", po::value<std::string>()->required(), "Output format (d - distance matrix, t - phylogenetic tree), required")
        ("algorithm,a", po::value<std::string>(), "Algorithm selection (0 - default mode, 1 - force placement, 2 - force conventional NJ)")
        ("placement-mode,p", po::value<std::string>(), "Placement mode selection (0 - exact mode, 1 - k-closest mode), default is k-closest")
        // ("batch-size,b", po::value<std::string>(), "Batch size for GPU processing (Default = 100000)")
        ("print-binary-newick",
         "Print fully resolved binary Newick with explicit :0 branches; default collapses near-zero internal edges")
        ("protein", "Input sequences are protein (for -i r raw sequences; default DNA/RNA)")
        ("help,h", "Print help messages");

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
    // std::cout << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
}


void readSequencesBatch(po::variables_map& vm, 
                        std::vector<std::string>& seqs, 
                        std::vector<std::string>& names)
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
    // std::cout << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
}

void readSequencesPhylip(po::variables_map& vm, 
                        std::vector<std::string>& seqs, 
                        std::vector<std::string>& names,
                        uint64_t ** fourBitCompressedSeqs,
                        uint64_t * seqLengths)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["input-file"].as<std::string>();

    std::ifstream file(seqFileName);
    if (!file) throw std::runtime_error("Could not open file: " + seqFileName);

    int nseq, seqlen;
    file >> nseq >> seqlen;
    std::cerr << "Number of sequences: " << nseq << ", Length of sequences: " << seqlen << std::endl;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip rest of line

    fourBitCompressedSeqs = new uint64_t*[nseq];
    seqLengths = new uint64_t[nseq];

    std::string line;

    for (int i = 0; i < nseq; ++i) {
        if (!std::getline(file, line)) {
            throw std::runtime_error("Unexpected end of file while reading sequences.");
        }

        if ((int)line.size() < 10) {
            throw std::runtime_error("Line too short for name and sequence: " + line);
        }

        std::string name = line.substr(0, 10);
        std::string seq = line.substr(10);

        // Trim whitespace
        name.erase(name.find_last_not_of(" \t\r\n") + 1);
        seq.erase(std::remove_if(seq.begin(), seq.end(), ::isspace), seq.end());

        if ((int)seq.size() != seqlen) {
            throw std::runtime_error("Sequence length mismatch for " + name + ": found " + std::to_string(seq.size()) + ", expected " + std::to_string(seqlen));
        }


        // seqs.push_back(seq);
        names.push_back(name);

        uint64_t fourBitCompressedSize = (seq.size()+15)/16;
        uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
        fourBitCompressor(seq, seq.size(), fourBitCompressed);

        seqLengths[i] = seq.size();
        fourBitCompressedSeqs[i] = fourBitCompressed;
    }

    // gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    // if (!f_rd) {
    //     fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
    //     exit(1);
    // }

    // kseq_t* kseq_rd = kseq_init(f_rd);

    // while (kseq_read(kseq_rd) >= 0) {
    //     size_t seqLen = kseq_rd->seq.l;
    //     seqs.push_back(std::string(kseq_rd->seq.s, seqLen));
    //     names.push_back(std::string(kseq_rd->name.s, kseq_rd->name.l));
    // }
    
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
        std::cerr << mainDesc << std::endl;
        // Return with error code 1 unless the user specifies help
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    MashPlacement::g_printBinaryNewick = vm.count("print-binary-newick") > 0;

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
    try {threshold= (uint64_t)std::stoi(vm["threshold"].as<std::string>());}
    catch(std::exception &e){}

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

    std::string placemode = "1";
    try {placemode = vm["algorithm"].as<std::string>();}
    catch(std::exception &e){}

    // uint64_t batchSize = 50000;
    // try {batchSize = (uint64_t)std::stoi(vm["batch-size"].as<std::string>());}
    // catch(std::exception &e){}

    int defau_thre = 10000; // Default threshold between conventional NJ and placement

    bool isProtein = vm.count("protein") > 0;

    MashPlacement::Param params(k, sketchSize, threshold, distanceType, in, out);
    params.isProtein = isProtein;

    size_t totalNumSequences = 0;
    size_t backboneSize = 0;

    if (in == "m" && out == "t"){
        std::vector<std::string> seqs,names;
        
        uint64_t ** fourBitCompressedSeqs;
        uint64_t * seqLengths;
        
        // readSequencesPhylip(vm, seqs, names, fourBitCompressedSeqs, seqLengths);
        std::string seqFileName = vm["input-file"].as<std::string>();

        std::ifstream file(seqFileName);
        if (!file) throw std::runtime_error("Could not open file: " + seqFileName);

        int nseq, seqlen;
        file >> nseq >> seqlen;
        std::cerr << "Number of sequences: " << nseq << ", Length of sequences: " << seqlen << std::endl;
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip rest of line

        std::vector<int> ids(nseq);
        for(int i=0;i<nseq;i++) ids[i]=i;
        std::mt19937 rnd(time(NULL));
        std::shuffle(ids.begin(),ids.end(),rnd);

        names.resize(nseq);

        fourBitCompressedSeqs = new uint64_t*[nseq];
        seqLengths = new uint64_t[nseq];

        std::string line;

        for (int i = 0; i < nseq; ++i) {
            if (!std::getline(file, line)) {
                throw std::runtime_error("Unexpected end of file while reading sequences.");
            }

            if ((int)line.size() < 10) {
                throw std::runtime_error("Line too short for name and sequence: " + line);
            }

            std::string name = line.substr(0, 10);
            std::string seq = line.substr(10);

            // Trim whitespace
            name.erase(name.find_last_not_of(" \t\r\n") + 1);
            seq.erase(std::remove_if(seq.begin(), seq.end(), ::isspace), seq.end());

            if ((int)seq.size() != seqlen) {
                throw std::runtime_error("Sequence length mismatch for " + name + ": found " + std::to_string(seq.size()) + ", expected " + std::to_string(seqlen));
            }


            // seqs.push_back(seq);
            names[ids[i]] = name;

            uint64_t fourBitCompressedSize = (seq.size()+15)/16;
            uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
            fourBitCompressor(seq, seq.size(), fourBitCompressed);

            seqLengths[ids[i]] = seq.size();
            fourBitCompressedSeqs[ids[i]] = fourBitCompressed;
        }
        
        std::cerr << "Sequences read successfully\n";

        size_t numSequences = names.size();
        // std::vector<int> ids(numSequences);
        // std::vector<std::string> temp2(numSequences);
        // // ,temp1(numSequences);
        // for(int i=0;i<numSequences;i++) ids[i]=i;
        // std::mt19937 rnd(time(NULL));
        // std::shuffle(ids.begin(),ids.end(),rnd);

        // for(int i=0;i<numSequences;i++){
        //     // temp1[i]=seqs[ids[i]];
        //     temp2[i]=names[ids[i]];
        // }
        // // seqs=temp1;
        // names=temp2;

        // // write the shuffled sequences to a file
        // std::ofstream outfile("shuffled_sequences.fasta");
        // for (size_t i = 0; i < seqs.size(); ++i) {
        //     outfile << ">" << names[i] << "\n" << seqs[i] << "\n";
        // }
        // outfile.close();


        // Compress Sequences (2-bit compressor)
        auto compressStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "Compressing input sequence using two-bit encoding.\n");
        
        auto compressEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds compressTime = compressEnd - compressStart;
        // std::cout << "Compressed in: " <<  compressTime.count() << " ns\n";
        auto inputEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds inputTime = inputEnd - inputStart; 
        std::cerr << "Input in: " <<  inputTime.count()/1000000 << " ms\n";

        // Create arrays
        auto createArrayStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "\nAllocating Gpu device arrays.\n");
        // std::cerr<<"########\n";
        
        totalNumSequences = numSequences;
        backboneSize = numSequences/20;
        params.batchSize = backboneSize;
        params.backboneSize = backboneSize;

        MashPlacement::msaDeviceArraysDC.allocateDeviceArraysDC(fourBitCompressedSeqs, seqLengths, numSequences, params);
        std::cerr << "Allocated Device Arrays" << std::endl;
        if(algo=="1"||algo=="0"&&numSequences>=defau_thre){
            if(placemode=="0"){
                std::cerr<<"Using exact placement mode\n";
                MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
                auto createArrayEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
                std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";


                //Build Tree on Gpu
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::placementDeviceArrays.printTree(names);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

                // Print first 10 hash values corresponding to each sequence
                // MashPlacement::mashDeviceArrays.printSketchValues(10);
                MashPlacement::msaDeviceArraysDC.deallocateDeviceArraysDC();
                MashPlacement::placementDeviceArrays.deallocateDeviceArrays();
            }
            else{
                std::cerr<<"Using k-closest placement mode\n";
                MashPlacement::kplacementDeviceArraysDC.allocateDeviceArraysDC(backboneSize, totalNumSequences);
                auto createArrayEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
                std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";

                //Build Tree on Gpu
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::kplacementDeviceArraysDC.findBackboneTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHost);
                MashPlacement::kplacementDeviceArraysDC.findClustersDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHost);
                MashPlacement::kplacementDeviceArraysDC.findClusterTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHost);

                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::kplacementDeviceArraysDC.printTreeDC(names);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

                // Print first 10 hash values corresponding to each sequence
                // MashPlacement::mashDeviceArrays.printSketchValues(10);
                MashPlacement::msaDeviceArraysDC.deallocateDeviceArraysDC();
                MashPlacement::kplacementDeviceArraysDC.deallocateDeviceArraysDC();
            }
        }
        else{
            std::cerr<<"Using conventional NJ\n";
            if(numSequences>=40000){
                std::cerr<<"Warning: forcing conventional NJ on large datasets might result in unexpected behavior\n";
            }
            MashPlacement::njDeviceArrays.getDismatrix(
                numSequences,params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC
            );
            MashPlacement::njDeviceArrays.findNeighbourJoiningTree(names);
            MashPlacement::msaDeviceArraysDC.deallocateDeviceArraysDC();
            MashPlacement::njDeviceArrays.deallocateDeviceArrays();
        }
    }
    else if (in == "r" && out == "t"){
        std::vector<std::string> seqs,names_, names;

        // Read Input Sequences (Fasta format)
        readSequences(vm, seqs, names_);
        size_t numSequences = seqs.size();
        std::vector<int> ids(numSequences);
        std::vector<std::string> temp1(numSequences),temp2(numSequences);
        for(int i=0;i<numSequences;i++) ids[i]=i;
        std::mt19937 rnd(time(NULL));
        std::shuffle(ids.begin(),ids.end(),rnd);
        names.resize(numSequences);

        // Compress Sequences (2-bit compressor)
        auto compressStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "Compressing input sequence using two-bit encoding.\n");
        uint64_t ** twoBitCompressedSeqs = new uint64_t*[numSequences];
        uint64_t * seqLengths = new uint64_t[numSequences];
        tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
        for (int idx_= range.begin(); idx_ < range.end(); ++idx_) 
        {   
            int i = idx_;
            uint64_t twoBitCompressedSize = isProtein ? (seqs[i].size()+7)/8 : (seqs[i].size()+31)/32;
            uint64_t * twoBitCompressed = new uint64_t[twoBitCompressedSize];
            twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed, isProtein);

            seqLengths[ids[i]] = seqs[i].size();
            twoBitCompressedSeqs[ids[i]] = twoBitCompressed;
            names[ids[i]] = names_[i];

        }});
        // delete seqs
        std::vector<std::string>().swap(seqs);

        auto compressEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds compressTime = compressEnd - compressStart;
        // std::cout << "Compressed in: " <<  compressTime.count() << " ns\n";
        auto inputEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds inputTime = inputEnd - inputStart; 
        std::cerr << "Input in: " <<  inputTime.count()/1000000 << " ms\n";

        totalNumSequences = numSequences;
        backboneSize = numSequences/20;
        params.batchSize = backboneSize;
        params.backboneSize = backboneSize;

        // Create arrays
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
        
        // MashPlacement::mashDeviceArrays.deallocateDeviceArrays();
        //Build Tree on Gpu
        if(algo=="1"||algo=="0"&&numSequences>=defau_thre){
            if(placemode=="0"){
                std::cerr<<"Using exact placement mode\n";
                MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::placementDeviceArrays.printTree(names);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";
                MashPlacement::mashDeviceArraysDC.deallocateDeviceArraysDC();
                MashPlacement::placementDeviceArrays.deallocateDeviceArrays();
            }
            else{
                std::cerr<<"Using k-closest placement mode\n";
                MashPlacement::kplacementDeviceArraysDC.allocateDeviceArraysDC(backboneSize, totalNumSequences);
                MashPlacement::kplacementDeviceArraysHost.allocateHostArrays(backboneSize, totalNumSequences);
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                
                MashPlacement::kplacementDeviceArraysDC.findBackboneTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHost);
                MashPlacement::kplacementDeviceArraysDC.findClustersDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHost);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart;

                MashPlacement::kplacementDeviceArraysDC.findClusterTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHost);
                MashPlacement::kplacementDeviceArraysDC.printTreeDC(names);
                MashPlacement::kplacementDeviceArraysDC.deallocateDeviceArraysDC();
                MashPlacement::mashDeviceArraysDC.deallocateDeviceArraysDC();
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";
                
            }
        }
        else{
            std::cerr<<"Using conventional NJ\n";
            if(numSequences>=40000){
                std::cerr<<"Warning: forcing conventional NJ on large datasets might result in unexpected behavior\n";
            }
            MashPlacement::njDeviceArrays.getDismatrix(
                numSequences,params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC
            );
            MashPlacement::njDeviceArrays.findNeighbourJoiningTree(names);
            MashPlacement::mashDeviceArraysDC.deallocateDeviceArraysDC();
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

        totalNumSequences = numSequences;
        backboneSize = numSequences/20;
        params.batchSize = backboneSize;
        params.backboneSize = backboneSize;

        if(algo=="1"||algo=="0"&&numSequences>=defau_thre){
            if(placemode=="0"){
                std::cerr<<"Using exact placement mode\n";
                MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
                MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC);
                MashPlacement::placementDeviceArrays.printTree(MashPlacement::matrixReader.name);
            }
            else{
                std::cerr<<"Using k-closest placement mode\n";
                MashPlacement::kplacementDeviceArraysDC.allocateDeviceArraysDC(backboneSize, totalNumSequences);
                MashPlacement::kplacementDeviceArraysHost.allocateHostArrays(backboneSize, totalNumSequences);
                MashPlacement::kplacementDeviceArraysDC.findBackboneTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHost);
                MashPlacement::kplacementDeviceArraysDC.printTreeDC(MashPlacement::matrixReader.name);
            }
        }
        else{
            std::cerr<<"Using conventional NJ\n";
            if(numSequences>=40000){
                std::cerr<<"Warning: forcing conventional NJ on large datasets might result in unexpected behavior\n";
            }
            MashPlacement::njDeviceArrays.getDismatrix(
                numSequences,params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC
            );
            MashPlacement::njDeviceArrays.findNeighbourJoiningTree(MashPlacement::matrixReader.name);
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
