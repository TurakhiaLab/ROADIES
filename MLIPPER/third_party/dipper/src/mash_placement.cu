/**
 * mash_placement.cu
 *
 * MASH sketch construction on GPU and k-closest placement. Allocates device arrays for
 * 2-bit compressed sequences, hash lists, prefix sums; builds MinHash sketches via
 * MurmurHash3; computes MASH distances; maintains k-closest backbone nodes per edge
 * and runs placement (build tree or add query onto backbone). Also holds 
 * DeviceArrays / findPlacementTree k-closest logic.
 */

#ifndef MASHPL_CUH
#include "mash_placement.cuh"
#endif

#include <stdio.h>
#include <queue>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <chrono>
#include <iostream>
#include <cub/cub.cuh>

/* d_aggseqLengths = prefix sum of compressed lengths, e.g. [3,5,3] -> [1,2,3]. */
void MashPlacement::DeviceArrays::allocateDeviceArrays(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params)
{
    cudaError_t err;

    d_numSequences = numSequences = int(num);

    uint64_t kmerSize = params.kmerSize;
    size_t hashListLength = 0;   

    // Allocate memory
    err = cudaMalloc(&d_aggseqLengths, numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_seqLengths, numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    uint64_t * h_aggseqLengths = new uint64_t[numSequences];
    uint64_t flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) flatStringLength+= (h_seqLengths[i]+31)/32;
    uint64_t * h_flattenCompressSeqs = new uint64_t[flatStringLength];
    flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) 
    {
        uint64_t flatStringLengthLocal = (h_seqLengths[i]+31)/32;
        hashListLength += h_seqLengths[i] - kmerSize + 1;
        flatStringLength+=flatStringLengthLocal;
        for (size_t j=0; j<flatStringLengthLocal;j++)  
        {
            h_flattenCompressSeqs[j] = h_compressedSeqs[i][j];
            // if (i==9) printf("%u\n",h_flattenCompressSeqs[j]); 
        }
        h_flattenCompressSeqs += flatStringLengthLocal;
        h_aggseqLengths[i] = flatStringLength;
    }

    h_flattenCompressSeqs -= flatStringLength;
    //printf("%d", flatStringLength);

    // err = cudaMalloc(&d_hashList, hashListLength*sizeof(uint64_t));
    err = cudaMalloc(&d_hashList, 1000*numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    // printf("%p", d_hashList);


    err = cudaMalloc(&d_compressedSeqs, flatStringLength*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_prefixCompressed, d_numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }


    // Transfer data
    err = cudaMemcpy(d_aggseqLengths, h_aggseqLengths, numSequences*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_seqLengths, h_seqLengths, numSequences*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_compressedSeqs, h_flattenCompressSeqs, flatStringLength*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    thrust::device_ptr<uint64_t> dev_seqLengths(d_seqLengths);
    thrust::device_ptr<uint64_t> dev_prefixComp(d_prefixCompressed);

    thrust::transform(thrust::device,
        dev_seqLengths, dev_seqLengths + d_numSequences, dev_prefixComp, 
        [] __device__ (const uint64_t& x) -> uint64_t { 
            return (x + 31) / 32;
        }
    );

    thrust::exclusive_scan(dev_prefixComp, dev_prefixComp + d_numSequences, dev_prefixComp);

    cudaDeviceSynchronize();

    bd = 2, idx = 0;
    err = cudaMalloc(&d_dist, numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_head, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_e, numSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_len, numSequences*8*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_nxt, numSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_belong, numSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_dis, numSequences*20*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_id, numSequences*20*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
}

void MashPlacement::DeviceArrays::deallocateDeviceArrays(){
    cudaFree(d_compressedSeqs);
    cudaFree(d_aggseqLengths);
    cudaFree(d_seqLengths);
    cudaFree(d_hashList);
    cudaFree(d_prefixCompressed);
    // cudaFree(d_mashDist);
}

#define BIG_CONSTANT(x) (x##LLU)

__device__ inline uint64_t rotl64 ( uint64_t x, int8_t r )
{
    return (x << r) | (x >> (64 - r));
}

#define ROTL64(x,y)    rotl64(x,y)

__device__ uint64_t getblock64 ( const uint64_t * p, int i )
{
    return p[i];
}

/** MurmurHash3 64-bit mix (fmix). */
__device__ uint64_t fmix64 ( uint64_t k )
{
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;

    return k;
}

/** MurmurHash3 x64 128-bit on raw k-mer; used for MASH MinHash. */
__device__ void MurmurHash3_x64_128_MASH ( void * key, const int len,
                           const uint32_t seed, void * out)
{
    const uint8_t * data = (const uint8_t*)key;
    const int nblocks = len / 16;

    uint64_t h1 = seed;
    uint64_t h2 = seed;

    const uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
    const uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

    //----------
    // body

    const uint64_t * blocks = (const uint64_t *)(data);

    for(int i = 0; i < nblocks; i++)
    {
        uint64_t k1 = getblock64(blocks,i*2+0);
        uint64_t k2 = getblock64(blocks,i*2+1);

        k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;

        h1 = ROTL64(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

        k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

        h2 = ROTL64(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
    }

    //----------
    // tail

    const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

    uint64_t k1 = 0;
    uint64_t k2 = 0;

    switch(len & 15)
    {
        case 15: k2 ^= ((uint64_t)tail[14]) << 48;
        case 14: k2 ^= ((uint64_t)tail[13]) << 40;
        case 13: k2 ^= ((uint64_t)tail[12]) << 32;
        case 12: k2 ^= ((uint64_t)tail[11]) << 24;
        case 11: k2 ^= ((uint64_t)tail[10]) << 16;
        case 10: k2 ^= ((uint64_t)tail[ 9]) << 8;
        case  9: k2 ^= ((uint64_t)tail[ 8]) << 0;
               k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

        case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
        case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
        case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
        case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
        case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
        case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
        case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
        case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
               k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
    };

    //----------
    // finalization

    h1 ^= len; h2 ^= len;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    ((uint64_t*)out)[0] = h1;
    ((uint64_t*)out)[1] = h2;
}


/** Decode 2-bit compressed k-mer into forward and reverse complement char arrays. */
__device__ void decompress(uint64_t compressedSeq, uint64_t kmerSize, char * decompressedSeq_fwd, char * decompressedSeq_rev) {
    static const char lookupTable[4] = {'A', 'C', 'G', 'T'};
    for (int i = kmerSize - 1; i >= 0; i--) {
        uint64_t twoBitVal = (compressedSeq >> (2 * i)) & 0x3;
        decompressedSeq_fwd[i] = lookupTable[twoBitVal];
        decompressedSeq_rev[kmerSize - 1 - i] = lookupTable[3 - twoBitVal];
    }
}

/** Lexicographic compare of forward vs reverse k-mer; used to pick canonical form. */
__device__ int memcmp_device(const char* kmer_fwd, const char* kmer_rev, int kmerSize) {
    for (int i = 0; i < kmerSize; i++) {
        if (kmer_fwd[i] < kmer_rev[i]) {
            return -1;
        }
        if (kmer_fwd[i] > kmer_rev[i]) {
            return 1;
        }
    }
    return 0;
}

/** Build MinHash sketch per sequence: slide over 2-bit compressed seq, canonical k-mer ->
 * MurmurHash3, block radix sort, keep top sketchSize hashes per sequence. */
__global__ void sketchConstruction(
    uint64_t * d_compressedSeqs,
    uint64_t * d_seqLengths,
    uint64_t * d_prefixCompressed,
    size_t d_numSequences,
    uint64_t * d_hashList,
    uint64_t kmerSize
) {
    extern __shared__ uint64_t stored[];

    typedef cub::BlockRadixSort<uint64_t, 512, 3> BlockRadixSort;
    __shared__ typename BlockRadixSort::TempStorage temp_storage;

    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int stride = gridDim.x;

    uint64_t kmer = 0;
    char kmer_fwd[32];
    char kmer_rev[32];
    uint8_t out[16];
    
    for (size_t i = bx; i < d_numSequences; i += stride) {
        //if (tx == 0) printf("Started block %d\n", i);
        stored[tx] = 0xFFFFFFFFFFFFFFFF; // reset block's stored values
        stored[tx + 500] = 0xFFFFFFFFFFFFFFFF;

        uint64_t  * hashList = d_hashList + 1000*i;

        uint64_t seqLength = d_seqLengths[i];
        uint64_t * compressedSeqs = d_compressedSeqs + d_prefixCompressed[i];

        size_t j = tx;
        size_t iteration = 0;
        while (iteration <= seqLength - kmerSize) {

            uint64_t keys[3];

            if (j <= seqLength - kmerSize) {

                uint64_t index = j/32;
                uint64_t shift1 = 2*(j%32);

                if (shift1>0) {
                    uint64_t shift2 = 64-shift1;
                    kmer = ((compressedSeqs[index] >> shift1) | (compressedSeqs[index+1] << shift2)); //& mask;
                }
                else {   
                    kmer = compressedSeqs[index];// & mask;
                }

                decompress(kmer, kmerSize, kmer_fwd, kmer_rev);

                // if ((i == 0) && (tx == 0)) {
                //     for (char c : kmer_fwd) printf("%c", c);   
                
                //     printf("\n");
                // }

                // convert to char representation and call w/ original
                MurmurHash3_x64_128_MASH( (memcmp_device(kmer_fwd, kmer_rev, kmerSize) <= 0) 
                    ? kmer_fwd : kmer_rev, kmerSize, 42, out);
                
                uint64_t hash = *((uint64_t *)out);

                // Combine stored and computed to sort and rank
                keys[0] = (tx < 500) ? stored[tx] : 0xFFFFFFFFFFFFFFFF;
                keys[1] = (tx < 500) ? stored[tx + 500] : 0xFFFFFFFFFFFFFFFF;
                keys[2] = hash;
            } else {
                keys[0] = (tx < 500) ? stored[tx] : 0xFFFFFFFFFFFFFFFF;
                keys[1] = (tx < 500) ? stored[tx + 500] : 0xFFFFFFFFFFFFFFFF;
                keys[2] = 0xFFFFFFFFFFFFFFFF;
            }

            BlockRadixSort(temp_storage).Sort(keys);
  
            __syncthreads();
            // // Move top 1000 hashes back to stored
            if (tx < 333) {
                stored[3*tx] = keys[0];
                stored[3*tx + 1] = keys[1];
                stored[3*tx + 2] = keys[2];
            } else if (tx == 333) {
                stored[999] = keys[0];
            }

            __syncthreads();

            iteration += blockDim.x;
            j += blockDim.x;

        }
    
        // Result writing back to global memory.
        if (tx < 500) {
            // d_hashList[i*1000 + tx] = stored[tx];
            // d_hashList[i*1000 + tx + 500] = stored[tx + 500];
            hashList[tx] = stored[tx];
            hashList[tx + 500] = stored[tx + 500];
        }

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));
        }
       
    }

}

/** Launch sketchConstruction for all sequences (single block per sequence in current impl). */
void MashPlacement::sketchConstructionOnGpu
(
    uint64_t * d_compressedSeqs,
    uint64_t * d_prefixCompressed,
    uint64_t * d_aggseqLengths,
    uint64_t * d_seqLengths,
    size_t d_numSequences,
    uint64_t * d_hashList,
    uint64_t * h_seqLengths,
    Param& params
){

    const uint64_t kmerSize = params.kmerSize; // Extract kmerSize

    auto timerStart = std::chrono::high_resolution_clock::now();

    int threadsPerBlock = 512;
    int blocksPerGrid = (d_numSequences + threadsPerBlock - 1) / threadsPerBlock;
    size_t sharedMemorySize = sizeof(uint64_t) * (2000);
    sketchConstruction<<<1, threadsPerBlock, sharedMemorySize>>>(
        d_compressedSeqs, d_seqLengths, d_prefixCompressed, d_numSequences, d_hashList, kmerSize
    );

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();

    auto timerEnd = std::chrono::high_resolution_clock::now();
    auto time = timerEnd - timerStart;
    // std::cout << "Time to generate hashes: " << time.count() << "ns\n";

}

// __device__ double mashDistance
// (
//     uint64_t * A,
//     uint64_t * B,
//     uint64_t kmerSize,
//     uint64_t sketchSize
// ){
//     uint32_t APtr = 0, BPtr = 0, inter = 0, uni = 0;

//     for(;uni<sketchSize;uni++)
//     {
//         if (A[APtr]==B[BPtr]) inter++, APtr++, BPtr++;
//         else if (A[APtr]>B[BPtr]) BPtr++;
//         else APtr++;
//     }

//     double jaccardEstimate = double(inter)/uni;

//     double mashDist = abs((log(2.0*jaccardEstimate/(1.0+jaccardEstimate)))/kmerSize);

//     return mashDist;

// }

// __global__ void mashDistConstruction
// (
//     uint64_t rowId,
//     uint64_t * d_hashList,
//     uint64_t * d_seqLengths,
//     double * d_mashDist,
//     uint64_t kmerSize,
//     uint64_t sketchSize
// ){
//     int tx = threadIdx.x;
//     int bx = blockIdx.x;
//     int bs = blockDim.x;
//     int idx = tx+bx*bs;
//     if(idx>=rowId) return;
//     // printf("%d %d %d\n",tx,bx,rowId);
//     uint64_t * hashList = d_hashList;
//     double mashDist = mashDistance(hashList+rowId*sketchSize, hashList+idx*sketchSize, kmerSize, sketchSize);
//     d_mashDist[idx] = mashDist;
// }

/** Compute MASH distances from row rowId to rows 0..rowId-1 using Jaccard via sketch overlap. */
__global__ void mashDistConstruction(
    int rowId,
    uint64_t * d_hashList,
    double * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences
) {
    int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x;
    int idx = tx+bx*bs;
    if(idx>=rowId) return;
    int uni = 0, bPos = rowId, inter = 0;
    uint64_t aval, bval;
    for(int i=idx; uni < sketchSize; i+=numSequences, uni++){
        aval = d_hashList[i];
        while(uni < sketchSize && bPos < numSequences * sketchSize){
            bval = d_hashList[bPos];
            if(bval > aval) break;
            if(bval < aval) uni++;
            else inter++;
            bPos += numSequences;
        }
        if(uni >= sketchSize) break;
    }
    double jaccardEstimate = double(inter)/uni;
    d_mashDist[idx] = abs((log(2.0*jaccardEstimate/(1.0+jaccardEstimate)))/kmerSize);
}


void MashPlacement::DeviceArrays::printSketchValues(int numValues) 
{
    uint64_t * h_hashList = new uint64_t[1000*d_numSequences];


    uint64_t * hashList = d_hashList;

    cudaError_t err;

    //printf("Total Hashes: %d", d_numSequences*1000);

    err = cudaMemcpy(h_hashList, hashList, d_numSequences*1000*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    // printf("i\thashList[i] (%zu)\n");
    for (int j = 0; j < d_numSequences; j++) {
        printf("Sequence (%d)\n", j);
        for (int i=0; i<numValues; i++) {
            printf("%i\t%lu\n", i, h_hashList[j*1000 + i]);
        }
    }
    

}




void MashPlacement::DeviceArrays::printMashDist(uint64_t h_numSequences, std::vector<std::string> seqs) 
{
    
    double * h_mashDist = new double[h_numSequences];

    double * mashDist = d_dist;

    cudaError_t err;


    err = cudaMemcpy(h_mashDist, mashDist, (h_numSequences)*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }
    printf("%d\n",h_numSequences);
    for (int i=0; i<h_numSequences; i++) 
    {
        std::cout<<seqs[i]<<'\t';
        for (int j=0; j<h_numSequences; j++) 
        {
            printf("%f\t", *h_mashDist);
            h_mashDist++;
        }
        printf("\n");        
    }
}


/** Reset adjacency and closest_id/closest_dis to sentinels. */
__global__ void initialize(
    int lim,
    int nodes,
    double * d_closest_dis,
    int * d_closest_id,
    int * head,
    int * nxt,
    int * belong,
    int * e
){
    int tx=threadIdx.x,bs=blockDim.x;
    int bx=blockIdx.x,gs=gridDim.x;
    int idx=tx+bs*bx;
    if(idx<lim){
        for(int i=0;i<5;i++){
            d_closest_dis[idx*5+i]=2;
            d_closest_id[idx*5+i]=-1;
        }
        nxt[idx] = -1;
        e[idx] = -1;
        belong[idx] = -1;
    }
    if(idx<nodes) head[idx] = -1;
}

/** Compare placement tuples by addLen (third element). */
struct compare_tuple
{
  __host__ __device__
  bool operator()(thrust::tuple<int,double,double> lhs, thrust::tuple<int,double,double> rhs)
  {
    return thrust::get<2>(lhs) < thrust::get<2>(rhs);
  }
};
/* Tuple: (edge_id, fracLen, addLen). */

__global__ void calculateBranchLength(
    int num, // should be bd, not numSequences 
    int * head,
    int * nxt,
    double * dis, 
    int * e, 
    double * len, 
    int * belong,
    thrust::tuple<int,double,double> * minPos,
    int lim,
    double * closest_dis,
    int * closest_id,
    int totSeqNum
){
    int tx=threadIdx.x,bs=blockDim.x,bx=blockIdx.x,gs=gridDim.x;
    int idx=tx+bs*bx;
    if(idx>=lim) return;
    if(idx>=num*4-4||belong[idx]<e[idx]){
        thrust::tuple <int,double,double> minTuple(0,0,2);
        minPos[bx*bs+tx]=minTuple;
        return;
    }
    int x=belong[idx],oth=e[idx];
    int eid=idx,otheid;
    double dis1=0, dis2=0, val;
    for(int i=0;i<5;i++)
        if(closest_id[eid*5+i]!=-1){
            val = dis[closest_id[eid*5+i]]-closest_dis[eid*5+i];
            if(val>dis1) dis1=val;
        }
    otheid=head[oth];
    while(e[otheid]!=x) assert(otheid!=-1),otheid=nxt[otheid];
    for(int i=0;i<5;i++)
        if(closest_id[otheid*5+i]!=-1){
            val = dis[closest_id[otheid*5+i]]-closest_dis[otheid*5+i];
            if(val>dis2) dis2=val;
        }
    double additional_dis=(dis1+dis2-len[eid])/2;
    if(additional_dis<0) additional_dis=0;
    dis1-=additional_dis,dis2-=additional_dis;
    if(dis1<0) dis1=0;
    if(dis2<0) dis2=0;
    if(dis1>len[eid]) additional_dis+=dis1-len[eid],dis1=len[eid];
    if(dis2>len[eid]) additional_dis+=dis2-len[eid],dis2=len[eid];
    // assert(dis1+dis2-1e-6<=len[eid]);
    double rest=len[eid]-dis1-dis2;
    dis1+=rest/2,dis2+=rest/2;
    thrust::tuple <int,double,double> minTuple(eid,dis1,additional_dis);
    minPos[bx*bs+tx]=minTuple;
}

/** BFS from x: update closest_id/closest_dis per edge when improving. */
__global__ void updateClosestNodes(
    int * head,
    int * nxt,
    int * e,
    double * len,
    double * closest_dis,
    int * closest_id,
    int x,
    int * id,
    int * from,
    double * dis
){
    int l=0,r=-1;
    id[++r]=x,dis[x]=0,from[x]=-1;
    while(l<=r){
        int node=id[l],fb=from[l];
        double d=dis[l];
        l++;
        for(int i=head[node];i!=-1;i=nxt[i]){
            if(e[i]==fb) continue;
            for(int j=0;j<5;j++){
                double nowd=closest_dis[i*5+j];
                if(nowd>d){
                    for(int k=4;k>j;k--){
                        closest_dis[i*5+k]=closest_dis[i*5+k-1];
                        closest_id[i*5+k]=closest_id[i*5+k-1];
                    }
                    closest_dis[i*5+j]=d;
                    closest_id[i*5+j]=x;
                    id[++r]=e[i],dis[r]=d+len[i],from[r]=node;
                    break;
                }
            }
        }
    }
}

/** Insert new leaf on edge eid; update closest_id/closest_dis for new edges. */
__global__ void updateTreeStructure(
    int * head,
    int * nxt,
    int * e,
    double * len,
    double * closest_dis,
    int * closest_id,
    int * belong,
    int eid,
    double fracLen,
    double addLen,
    int placeId,
    int edgeCount,
    int numSequences
){
    int middle=placeId+numSequences-1, outside=placeId;
    int x=belong[eid],y=e[eid];
    double originalDis=len[eid];
    int xe,ye;
    for(int i=head[x];i!=-1;i=nxt[i])
        if(e[i]==y){
            e[i]=middle,len[i]=fracLen,xe=i;
            break;
        }
    for(int i=head[y];i!=-1;i=nxt[i])
        if(e[i]==x){
            e[i]=middle,len[i]-=fracLen,ye=i;
            break;
        }
    /*
    Need to update:
    e, len, nxt, head, belong, closest_dis, closest_id
    */
    //middle -> x
    e[edgeCount]=x,len[edgeCount]=fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    for(int i=0;i<5;i++)
        if(closest_id[ye*5+i]!=-1){
            closest_id[edgeCount*5+i]=closest_id[ye*5+i];
            closest_dis[edgeCount*5+i]=closest_dis[ye*5+i]+originalDis-fracLen;
        }
    edgeCount++;
    //middle -> y
    e[edgeCount]=y,len[edgeCount]=originalDis-fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    for(int i=0;i<5;i++)
        if(closest_id[xe*5+i]!=-1){
            closest_id[edgeCount*5+i]=closest_id[xe*5+i];
            closest_dis[edgeCount*5+i]=closest_dis[xe*5+i]+fracLen;
        }
    edgeCount++;
    //outside -> middle
    e[edgeCount]=middle,len[edgeCount]=addLen,nxt[edgeCount]=head[outside],head[outside]=edgeCount,belong[edgeCount]=outside;
    edgeCount++;
    //middle -> outside
    e[edgeCount]=outside,len[edgeCount]=addLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    int e1=edgeCount-2, e2=edgeCount-3;
    for(int i=0;i<5;i++){
        if(closest_id[e1*5+i]==-1) break;
        for(int j=0;j<5;j++)
            if(closest_dis[edgeCount*5+j]>closest_dis[e1*5+i]){
                for(int k=4;k>j;k--){
                    closest_dis[edgeCount*5+k]=closest_dis[edgeCount*5+k-1];
                    closest_id[edgeCount*5+k]=closest_id[edgeCount*5+k-1];
                }
                closest_dis[edgeCount*5+j]=closest_dis[e1*5+i];
                closest_id[edgeCount*5+j]=closest_id[e1*5+i];
                break;
            }
    }
    for(int i=0;i<5;i++){
        if(closest_id[e2*5+i]==-1) break;
        for(int j=0;j<5;j++)
            if(closest_dis[edgeCount*5+j]>closest_dis[e2*5+i]){
                for(int k=4;k>j;k--){
                    closest_dis[edgeCount*5+k]=closest_dis[edgeCount*5+k-1];
                    closest_id[edgeCount*5+k]=closest_id[edgeCount*5+k-1];
                }
                closest_dis[edgeCount*5+j]=closest_dis[e2*5+i];
                closest_id[edgeCount*5+j]=closest_id[e2*5+i];
                break;
            }
    }
    edgeCount++;
}

/** Build 3-node tree (0, 1, nv) with two edges of length dis[0]/2. */
__global__ void buildInitialTree(
    int numSequences,
    int * head,
    int * e,
    double * len,
    int * nxt,
    int * belong,
    double * dis,
    int edgeCount
){
    int nv = numSequences;
    double d = dis[0];
    // 0 -> nv
    e[edgeCount]=nv,len[edgeCount]=d/2,nxt[edgeCount]=head[0],head[0]=edgeCount,belong[edgeCount]=0;
    edgeCount++;
    // 1 -> nv
    e[edgeCount]=nv,len[edgeCount]=d/2,nxt[edgeCount]=head[1],head[1]=edgeCount,belong[edgeCount]=1;
    edgeCount++;
    // nv -> 0
    e[edgeCount]=0,len[edgeCount]=d/2,nxt[edgeCount]=head[nv],head[nv]=edgeCount,belong[edgeCount]=nv;
    edgeCount++;
    // nv -> 0
    e[edgeCount]=1,len[edgeCount]=d/2,nxt[edgeCount]=head[nv],head[nv]=edgeCount,belong[edgeCount]=nv;
    edgeCount++;
}

/** Transpose hash list from [seq][sketch] to [sketch][seq] for distance kernel. */
__global__ void rearrangeHashList(
    int numSequences,
    int sketchSize,
    uint64_t * original,
    uint64_t * target
){
    int tx = threadIdx.x, bx = blockIdx.x;
    int bs = blockDim.x;
    int idx = tx+bs*bx;
    if(idx>=numSequences) return;
    for(int i=0;i<sketchSize;i++){
        target[i*numSequences+idx] = original[idx*sketchSize + i];
    }
}

/** K-closest placement: rearrange hashes, init, build initial tree, updateClosestNodes;
 * then for each new taxon: MASH dist -> calculateBranchLength -> updateTree -> updateClosestNodes. */
void MashPlacement::findPlacementTree(
    int numSequences,
    int bd,
    int idx,
    double * d_dist,
    int * d_head,
    int * d_e,
    double * d_len,
    int * d_nxt,
    int * d_belong,
    double * d_closest_dis,
    int * d_closest_id,
    std::vector<std::string> name,
    uint64_t * d_hashList,
    uint64_t * d_seqLengths,
    Param& params
){
    int * h_head = new int[numSequences*2];
    int * h_e = new int[numSequences*8];
    int * h_nxt = new int[numSequences*8];
    double * h_len = new double[numSequences*8];
    double * h_closest_dis = new double[numSequences*20];
    int * h_closest_id = new int[numSequences*20];
    int * d_id;
    uint64_t * temp_hashList;
    auto err = cudaMalloc(&d_id, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int * d_from;
    err = cudaMalloc(&d_from, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    double * d_dis;
    err = cudaMalloc(&d_dis, numSequences*2*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&temp_hashList, numSequences*1000*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    std::function<void(int,int)>  print=[&](int node, int from){
        if(h_nxt[h_head[node]]!=-1){
            printf("(");
            std::vector <int> pos;
            for(int i=h_head[node];i!=-1;i=h_nxt[i])
                if(h_e[i]!=from)
                    pos.push_back(i);
            for(size_t i=0;i<pos.size();i++){
                print(h_e[pos[i]],node);
                printf(":");
                printf("%.5g%c",h_len[pos[i]],i+1==pos.size()?')':',');
            }
        }
        else std::cout<<name[node];
    };
    /*
    Initialize closest nodes by inifinite
    */
    int threadNum = 256, blockNum = (numSequences+255)/threadNum;
    rearrangeHashList <<<blockNum, threadNum >>>(
        numSequences,
        int(params.sketchSize),
        d_hashList,
        temp_hashList
    );
    std::swap(d_hashList, temp_hashList);
    cudaFree(temp_hashList);
    threadNum = 256, blockNum = (numSequences*4-4+threadNum-1)/threadNum;
    initialize <<<blockNum, threadNum>>> (
        numSequences*4-4,
        numSequences*2,
        d_closest_dis,
        d_closest_id,
        d_head,
        d_nxt,
        d_belong,
        d_e
    );
    /*
    Build Initial Tree
    */
    mashDistConstruction <<<1, 1>>> (
        1,
        d_hashList,
        d_dist,
        params.kmerSize, 
        params.sketchSize,
        numSequences
    );

    buildInitialTree <<<1,1>>> (
        numSequences,
        d_head,
        d_e,
        d_len,
        d_nxt,
        d_belong,
        d_dist,
        idx
    );
    idx += 4;
    /*
    Initialize closest nodes by inital tree
    */
    for(int i=0;i<bd;i++){
        updateClosestNodes <<<1,1>>> (
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis
        );
    }
 
    thrust::device_vector <thrust::tuple<int,double,double>> minPos(numSequences*4-4);
    // std::cout<<"FFF\n";
    std::chrono::nanoseconds disTime(0), treeTime(0);
    for(int i=bd;i<numSequences;i++){
        auto disStart = std::chrono::high_resolution_clock::now();
        blockNum = (i + 255) / 256;
        mashDistConstruction <<<blockNum, threadNum>>> (
            i,
            d_hashList,
            d_dist,
            params.kmerSize, 
            params.sketchSize,
            numSequences
        );
        cudaDeviceSynchronize();
        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        blockNum = (numSequences*4-4 + 255) / 256;
        calculateBranchLength <<<blockNum,threadNum>>> (
            i,
            d_head,
            d_nxt,
            d_dist,
            d_e,
            d_len,
            d_belong,
            thrust::raw_pointer_cast(minPos.data()),
            numSequences*4-4,
            d_closest_dis,
            d_closest_id,
            numSequences
        );
        auto iter=thrust::min_element(minPos.begin(),minPos.end(),compare_tuple());
        thrust::tuple<int,double,double> smallest=*iter;
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid=thrust::get<0>(smallest);
        double fracLen=thrust::get<1>(smallest),addLen=thrust::get<2>(smallest);
        updateTreeStructure <<<1,1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            d_belong,
            eid,
            fracLen,
            addLen,
            i,
            idx,
            numSequences
        );
        idx+=4;
        /*
        Update closest nodes
        */
        updateClosestNodes <<<1,1>>> (
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis
        );
        cudaDeviceSynchronize();
        auto treeEnd = std::chrono::high_resolution_clock::now();
        disTime += disEnd - disStart;
        treeTime += treeEnd - treeStart;
        /*
        ----------------------Testing-------------------------
        */
        // std::cout<<i<<" "<<eid<<" "<<fracLen<<" "<<addLen<<'\n';
        // cudaError_t err;
        // err = cudaMemcpy(h_head, d_head, numSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // err = cudaMemcpy(h_e, d_e, numSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // err = cudaMemcpy(h_nxt, d_nxt, numSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // err = cudaMemcpy(h_len, d_len, numSequences*8*sizeof(double),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // err = cudaMemcpy(h_closest_dis, d_closest_dis, numSequences*20*sizeof(double),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // err = cudaMemcpy(h_closest_id, d_closest_id, numSequences*20*sizeof(int),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // print(numSequences+bd-2,-1);
        // std::cout<<";\n";
        // for(int j=0;j<i*4;j++){
        //     std::cout<<"Edge "<<j<<" to "<<h_e[j]<<":\n";
        //     for(int k=0;k<5;k++) printf("%.5g %d\n",h_closest_dis[j*5+k],h_closest_id[j*5+k]);
        // }
        /*
        ----------------Testing Ends-----------------------------
        */
    }
    /*
    Copy the tree back to CPU and output
    */
    // cudaError_t err;
    err = cudaMemcpy(h_head, d_head, numSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_e, d_e, numSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_nxt, d_nxt, numSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_len, d_len, numSequences*8*sizeof(double),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    print(numSequences+bd-2,-1);
    std::cout<<";\n";
    std::cerr << "Distance Operation Time " <<  disTime.count()/1000000 << " ms\n";
    std::cerr << "Tree Operation Time " <<  treeTime.count()/1000000 << " ms\n";
}