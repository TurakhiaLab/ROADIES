
#include <hip/hip_runtime.h>
#include "mash_placement.cuh.hip"

#include <stdio.h>
#include <queue>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <chrono>
#include <iostream>
#include <hipcub/hipcub.hpp>

void MashPlacement::MashDeviceArrays::allocateDeviceArrays(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params)
{
    hipError_t err;

    numSequences = int(num);

    uint64_t kmerSize = params.kmerSize;
    size_t hashListLength = 0;   

    // Allocate memory
    err = hipMalloc(&d_aggseqLengths, numSequences*sizeof(uint64_t));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }
    err = hipMalloc(&d_seqLengths, numSequences*sizeof(uint64_t));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }

    /* Flatten data: word count per sequence is DNA (seq+31)/32 or protein (seq+7)/8. */
    uint64_t * h_aggseqLengths = new uint64_t[numSequences];
    uint64_t flatStringLength=0;
    auto wordsPerSeq = [&](uint64_t len) -> uint64_t {
        return params.isProtein ? (len + 7) / 8 : (len + 31) / 32;
    };
    for (size_t i =0; i<numSequences; i++) flatStringLength+= wordsPerSeq(h_seqLengths[i]);
    uint64_t * h_flattenCompressSeqs = new uint64_t[flatStringLength];
    flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) 
    {
        uint64_t flatStringLengthLocal = wordsPerSeq(h_seqLengths[i]);
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
    err = hipMalloc(&d_hashList, 1000*numSequences*sizeof(uint64_t));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }

    // printf("%p", d_hashList);


    err = hipMalloc(&d_compressedSeqs, flatStringLength*sizeof(uint64_t));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }

    err = hipMalloc(&d_prefixCompressed, numSequences*sizeof(uint64_t));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }


    // Transfer data
    err = hipMemcpy(d_aggseqLengths, h_aggseqLengths, numSequences*sizeof(uint64_t), hipMemcpyHostToDevice);
    if (err != hipSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    err = hipMemcpy(d_seqLengths, h_seqLengths, numSequences*sizeof(uint64_t), hipMemcpyHostToDevice);
    if (err != hipSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    err = hipMemcpy(d_compressedSeqs, h_flattenCompressSeqs, flatStringLength*sizeof(uint64_t), hipMemcpyHostToDevice);
    if (err != hipSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    // generate prefix array (per-sequence compressed word count)
    thrust::device_ptr<uint64_t> dev_seqLengths(d_seqLengths);
    thrust::device_ptr<uint64_t> dev_prefixComp(d_prefixCompressed);

    if (params.isProtein) {
        thrust::transform(thrust::device,
            dev_seqLengths, dev_seqLengths + numSequences, dev_prefixComp,
            [] __device__ (const uint64_t& x) -> uint64_t { return (x + 7) / 8; });
    } else {
        thrust::transform(thrust::device,
            dev_seqLengths, dev_seqLengths + numSequences, dev_prefixComp,
            [] __device__ (const uint64_t& x) -> uint64_t { return (x + 31) / 32; });
    }

    thrust::exclusive_scan(dev_prefixComp, dev_prefixComp + numSequences, dev_prefixComp);

    hipDeviceSynchronize();
}

void MashPlacement::MashDeviceArrays::deallocateDeviceArrays(){
    hipFree(d_compressedSeqs);
    hipFree(d_aggseqLengths);
    hipFree(d_seqLengths);
    hipFree(d_prefixCompressed);
    // cudaFree(d_hashList);
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

__device__ uint64_t fmix64 ( uint64_t k )
{
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;

    return k;
}

// First hashing function using raw sequence
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


__device__ void decompress(uint64_t compressedSeq, uint64_t kmerSize, char * decompressedSeq_fwd, char * decompressedSeq_rev) {
    static const char lookupTable[4] = {'A', 'C', 'G', 'T'};
    for (int i = kmerSize - 1; i >= 0; i--) {
        uint64_t twoBitVal = (compressedSeq >> (2 * i)) & 0x3;
        decompressedSeq_fwd[i] = lookupTable[twoBitVal];
        decompressedSeq_rev[kmerSize - 1 - i] = lookupTable[3 - twoBitVal];
    }
}

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

__device__ void decompress_protein_kmer(uint64_t * compressedSeqs, uint64_t j, uint64_t kmerSize, char * kmer_fwd) {
    static const char dec[21] = {
        'A','C','G','T','D','E','F','H','I','K','L','M','N','P','Q','R','S','V','W','Y','X'
    };
    for (uint64_t p = 0; p < kmerSize; p++) {
        uint64_t pos = j + p;
        uint64_t wi = pos / 8;
        int sh = (int)((pos % 8) * 8);
        uint64_t val = (compressedSeqs[wi] >> sh) & 0xFFULL;
        if (val > 20) val = 20;
        kmer_fwd[p] = dec[val];
    }
}

__global__ void sketchConstruction(
    uint64_t * d_compressedSeqs,
    uint64_t * d_seqLengths,
    uint64_t * d_prefixCompressed,
    size_t numSequences,
    uint64_t * d_hashList,
    uint64_t kmerSize,
    int isProtein
) {
    extern __shared__ uint64_t stored[];

    typedef hipcub::BlockRadixSort<uint64_t, 512, 3> BlockRadixSort;
    __shared__ typename BlockRadixSort::TempStorage temp_storage;

    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int stride = gridDim.x;

    uint64_t kmer = 0;
    char kmer_fwd[32];
    char kmer_rev[32];
    uint8_t out[16];
    
    for (size_t i = bx; i < numSequences; i += stride) {
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

                if (isProtein) {
                    decompress_protein_kmer(compressedSeqs, j, kmerSize, kmer_fwd);
                    MurmurHash3_x64_128_MASH(kmer_fwd, (int)kmerSize, 42, out);
                } else {
                    uint64_t index = j/32;
                    uint64_t shift1 = 2*(j%32);

                    if (shift1>0) {
                        uint64_t shift2 = 64-shift1;
                        kmer = ((compressedSeqs[index] >> shift1) | (compressedSeqs[index+1] << shift2));
                    }
                    else {
                        kmer = compressedSeqs[index];
                    }

                    decompress(kmer, kmerSize, kmer_fwd, kmer_rev);
                    MurmurHash3_x64_128_MASH( (memcmp_device(kmer_fwd, kmer_rev, kmerSize) <= 0)
                        ? kmer_fwd : kmer_rev, kmerSize, 42, out);
                }

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
       
    }

}

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

void MashPlacement::MashDeviceArrays::sketchConstructionOnGpu(Param& params){
    const uint64_t kmerSize = params.kmerSize; // Extract kmerSize
    auto timerStart = std::chrono::high_resolution_clock::now();

    int threadsPerBlock = 512;
    int blocksPerGrid = (numSequences + threadsPerBlock - 1) / threadsPerBlock;
    size_t sharedMemorySize = sizeof(uint64_t) * (2000);
    sketchConstruction<<<blocksPerGrid, threadsPerBlock, sharedMemorySize>>>(
        d_compressedSeqs, d_seqLengths, d_prefixCompressed, numSequences, d_hashList, kmerSize,
        params.isProtein ? 1 : 0
    );

    uint64_t * temp_hashList;
    auto err = hipMalloc(&temp_hashList, numSequences*int(params.sketchSize)*sizeof(uint64_t));
    if (err != hipSuccess){
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }
    rearrangeHashList <<<blocksPerGrid, threadsPerBlock >>>(
        numSequences,
        int(params.sketchSize),
        d_hashList,
        temp_hashList
    );
    std::swap(d_hashList, temp_hashList);
    hipFree(temp_hashList);

    
    err = hipGetLastError();
    if (err != hipSuccess) {
        printf("CUDA Error: %s\n", hipGetErrorString(err));
    }
    hipDeviceSynchronize();
    h_hashList = new uint64_t[1000*numSequences];
    hipMemcpy(h_hashList,d_hashList,1000*numSequences*sizeof(uint64_t),hipMemcpyDeviceToHost);
    auto timerEnd = std::chrono::high_resolution_clock::now();
    auto time = timerEnd - timerStart;
    // std::cout << "Time to generate hashes: " << time.count() << "ns\n";

}

__global__ void mashDistConstruction(
    int rowId,
    uint64_t * d_hashList,
    double * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences
) {
    int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x, gs = gridDim.x;
    int idx = tx+bx*bs;
    for (int idx_=idx; idx_<rowId;idx_+=bs*gs){
        if(idx_>=rowId) return;
        int uni = 0, bPos = rowId, inter = 0;
        uint64_t aval, bval;
        for(int i=idx_; uni < sketchSize; i+=numSequences, uni++){
            aval = d_hashList[i];
            while(uni < sketchSize && bPos < numSequences * sketchSize){
                bval = d_hashList[bPos];
                // printf("%ull %ull\n",aval,bval);
                if(bval > aval) break;
                if(bval < aval) uni++;
                else inter++;
                bPos += numSequences;
            }
            if(uni >= sketchSize) break;
        }
        assert(uni==1000);
        // assert(inter!=0);
        double jaccardEstimate = max(double(inter),1.0)/uni;
        d_mashDist[idx_] = min(1.0, abs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
    }
}

void MashPlacement::MashDeviceArrays::distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const{
    int threadNum = 1024, blockNum = 1024;
    mashDistConstruction <<<blockNum, threadNum>>> (
        rowId, 
        d_hashList, 
        d_mashDist, 
        params.kmerSize, 
        params.sketchSize, 
        numSequences
    );

    hipError_t err = hipDeviceSynchronize();
    if (err != hipSuccess) {
        printf("CUDA Error: %s\n", hipGetErrorString(err));
    }

    // double * h_dist = new double[numSequences];
    // auto err = cudaMemcpy(h_dist, d_mashDist, numSequences*sizeof(double), cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess) {
    //     fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
    //     exit(1);
    // }

    // fprintf(stderr, "Row (%d)\n", rowId);
    // for(int i=0;i<10;i++) fprintf(stderr, "%.5lf\n", h_dist[i]);
    // fprintf(stderr, "\n");
}


void MashPlacement::MashDeviceArrays::printSketchValues(int numValues) 
{
    uint64_t * h_hashList = new uint64_t[1000*numSequences];


    uint64_t * hashList = d_hashList;

    hipError_t err;

    //printf("Total Hashes: %d", numSequences*1000);

    err = hipMemcpy(h_hashList, hashList, numSequences*1000*sizeof(uint64_t), hipMemcpyDeviceToHost);
    if (err != hipSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    // printf("i\thashList[i] (%zu)\n");
    for (int j = 0; j < numSequences; j++) {
        fprintf(stderr, "Sequence (%d)\n", j);
        for (int i=0; i<numValues; i++) {
            fprintf(stderr, "%i\t%lu\n", i, h_hashList[i*numSequences+j]);
        }
    }
    

}