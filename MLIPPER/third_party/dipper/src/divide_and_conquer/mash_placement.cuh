/**
 * divide_and_conquer/mash_placement.cuh
 *
 * Declarations for divide-and-conquer mode: Param, MashDeviceArraysDC,
 * MSADeviceArraysDC, MatrixReader, PlacementDeviceArrays, KPlacementDeviceArraysHostDC,
 * NJDeviceArrays, KPlacementDeviceArraysDC. Used by placement_close_k.cu for
 * backbone construction, clustering, and per-cluster tree building.
 */

#ifndef MASHDC_CUH
#define MASHDC_CUH

#include <stdint.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>
#include <cstdio>
#include <string>
#include <utility>

namespace MashPlacement
{
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
        std::pair<int, int> range = {-1, -1};
        bool isProtein = false;

        Param(uint64_t t_kmerSize, uint64_t t_sketchSize, uint64_t t_threshold, uint64_t t_distanceType, std::string t_in, std::string t_out)
        {
            kmerSize = t_kmerSize; sketchSize = t_sketchSize; threshold = t_threshold, distanceType=t_distanceType;
            in = t_in, out = t_out; 
        };
    };

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

    struct MSADeviceArraysDC{
        uint64_t * d_compressedSeqsBackBone;
        uint64_t * d_compressedSeqsConst;
        // uint64_t * d_seqLengths;
        // size_t numSequences;
        int d_seqLen;

        uint64_t * h_compressedSeqs;

        size_t      totalNumSequences;
        size_t      backboneSize;

        void allocateDeviceArraysDC (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params);
        void distConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist) const;
        void distConstructionOnGpuForBackboneDC(Param& params, int rowId, double* d_mashDist) const;
        void distRangeConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist, int l, int r, bool clustering = false) const;
        void distSpecialIDConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist, int numToConstruct, int* d_id, int * d_leafMap) const;
        void deallocateDeviceArraysDC ();

    };
    static MSADeviceArraysDC msaDeviceArraysDC;

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
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays
        );
        void printTree(std::vector <std::string> name);
    };
    static PlacementDeviceArrays placementDeviceArrays;

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

        void allocateHostArraysDC (size_t num, size_t totalNum);
        void deallocateHostArraysDC ();
        void findClusterTreeDC(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays
        );
        void printTreeCpuDC(std::vector <std::string> name);
    };
    static KPlacementDeviceArraysHostDC kplacementDeviceArraysHostDC;

    struct NJDeviceArrays
    {
        int d_numSequences;
        double * d_mashDist;
        double * d_U;
        void deallocateDeviceArrays ();
        void getDismatrix(
            int numSequences,
            Param& params,
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays
        );
        void findNeighbourJoiningTree(std::vector <std::string> &name);
    };
    static NJDeviceArrays njDeviceArrays;

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

        void allocateDeviceArraysDC (size_t num, size_t totalNum);
        void deallocateDeviceArraysDC ();
        
        void findBackboneTreeDC(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays,
            const KPlacementDeviceArraysHost& kplacementDeviceArraysHost
        );

        void findClustersDC(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays,
            KPlacementDeviceArraysHost& kplacementDeviceArraysHost
        );

        void findClusterTreeDC(
            Param& params,
            MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            MSADeviceArraysDC& msaDeviceArrays,
            KPlacementDeviceArraysHost& kplacementDeviceArraysHost
        );
        void printTreeDC(std::vector <std::string> name);
    };
    static KPlacementDeviceArraysDC kplacementDeviceArraysDC;

};

#endif