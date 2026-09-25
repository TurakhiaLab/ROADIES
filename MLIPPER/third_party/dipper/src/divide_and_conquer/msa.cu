#include "../mash_placement.cuh"

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

void MashPlacement::MSADeviceArraysDC::transferToDeviceArraysDC(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, int gpuCluster, Param& params)
{
    cudaError_t err;
    
    size_t maxLengthCompressed = (this->d_seqLen + 15) / 16;

    /* Flatten complete data */
    
    int start = gpuCluster*params.backboneSize;
    int end = (start + params.backboneSize <= num) ? (start + params.backboneSize) : num;
    uint64_t flatStringLength= maxLengthCompressed*(end-start);

    uint64_t* compressedSeqs = new uint64_t[flatStringLength];
    // if (gpuCluster == 0 || this->h_compressedSeqs == nullptr)
    //     this->h_compressedSeqs = new uint64_t[flatStringLength];
    std::cout << "Initializing for " << gpuCluster << std::endl;
    
    for (size_t i = start; i < end; i++) 
    {
        for (size_t j=0; j<maxLengthCompressed;j++)  
        {
            // this->h_compressedSeqs[j] = h_compressedSeqs[i][j];
            compressedSeqs[j] = h_compressedSeqs[i][j];
        }
        // this->h_compressedSeqs += maxLengthCompressed;
        compressedSeqs += maxLengthCompressed;
    }
    // this->h_compressedSeqs -= flatStringLength;
    compressedSeqs -= flatStringLength;


    /* Transfer only the backbone data */
    err = cudaMemcpy(d_compressedSeqsConst, compressedSeqs, 1ll*(flatStringLength)*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }
    delete[] compressedSeqs;

    cudaDeviceSynchronize();
}


void MashPlacement::MSADeviceArraysDC::allocateDeviceArraysDC(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params, int gpuNum)
{
    cudaError_t err;
    this->totalNumSequences = int(num);
    this->backboneSize = params.backboneSize;
    
    this->d_seqLen = h_seqLengths[0];
    
    size_t maxLengthCompressed = (this->d_seqLen + 15) / 16;

    err = cudaMalloc(&d_compressedSeqsBackBone, 2*maxLengthCompressed*this->backboneSize*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    // The interleaved host-side D&C scheduler may upload an entire cluster of
    // query sequences here, which can be larger than the initial backbone.
    // Keep backboneSize as the semantic threshold used by the kernels, but
    // size this scratch buffer for the full sequence set.
    err = cudaMalloc(&d_compressedSeqsConst, maxLengthCompressed*this->totalNumSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    /* Flatten complete data */
    uint64_t flatStringLength= maxLengthCompressed*this->totalNumSequences;
    this->h_compressedSeqs = new uint64_t[flatStringLength];
    for (size_t i =0; i<this->totalNumSequences; i++) 
    {
        for (size_t j=0; j<maxLengthCompressed;j++)  
        {
            this->h_compressedSeqs[j] = h_compressedSeqs[i][j];
        }
        this->h_compressedSeqs += maxLengthCompressed;
    }
    this->h_compressedSeqs -= flatStringLength;


    /* Transfer only the backbone data */
    err = cudaMemcpy(d_compressedSeqsBackBone, this->h_compressedSeqs, 1ll*(maxLengthCompressed*this->backboneSize)*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    cudaDeviceSynchronize();
}

void MashPlacement::MSADeviceArraysDC::deallocateDeviceArraysDC(){
    cudaFree(d_compressedSeqsBackBone);
    cudaFree(d_compressedSeqsConst);
    d_compressedSeqsBackBone = nullptr;
    d_compressedSeqsConst = nullptr;
    delete[] h_compressedSeqs;
    h_compressedSeqs = nullptr;
    // cudaFree(d_seqLengths);
    // cudaFree(d_hashList);
    // cudaFree(d_mashDist);
}

#define DIST_UNCORRECTED 1
#define DIST_JUKESCANTOR 2
#define DIST_TAJIMANEI 3
#define DIST_KIMURA2P 4
#define DIST_TAMURA 5
#define DIST_JINNEI 6


__device__ void calculateParamsDC(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int & useful, int & match){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et<4&&ec<4) useful++;
            if(et<4&&et==ec) match++;
        }
    }
}

__device__ void calculateParamsDCParallel(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int & useful, int & match){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;

    __shared__ int sharedUseful[1024];
    __shared__ int sharedMatch[1024];
    sharedUseful[tx]=0;
    sharedMatch[tx]=0;

    for(int i=tx;i<compLen;i+=1024){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et<4&&ec<4) sharedUseful[tx]++;
            if(et<4&&et==ec) sharedMatch[tx]++;
        }
    }

    __syncthreads();

    // reduction
    for(int stride=1024/2; stride>0; stride/=2){
        if(tx<stride){
            sharedUseful[tx] += sharedUseful[tx + stride];
            sharedMatch[tx] += sharedMatch[tx + stride];
        }
        __syncthreads();
    }

    if(tx==0){
        useful=sharedUseful[0];
        match=sharedMatch[0];
    }
    __syncthreads();
}

__device__ void calculateParamsBatchDC(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, uint64_t * compressedSeqsConst, int & useful, int & match){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;

    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqsConst[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et<4&&ec<4) useful++;
            if(et<4&&et==ec) match++;
        }
    }
}

__device__ void calculateParamsBatchDCParallel(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, uint64_t * compressedSeqsConst, int & useful, int & match){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;

    __shared__ int sharedUseful[1024];
    __shared__ int sharedMatch[1024];
    sharedUseful[tx]=0;
    sharedMatch[tx]=0;

    for(int i=tx;i<compLen;i+=1024){
        long long vt=compressedSeqs[px+i], vc=compressedSeqsConst[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et<4&&ec<4) sharedUseful[tx]++;
            if(et<4&&et==ec) sharedMatch[tx]++;
        }
    }
    __syncthreads();

    // reduction
    for(int stride=1024/2; stride>0; stride/=2){
        if(tx<stride){
            sharedUseful[tx] += sharedUseful[tx + stride];
            sharedMatch[tx] += sharedMatch[tx + stride];
        }
        __syncthreads();
    }

    if(tx==0){
        useful=sharedUseful[0];
        match=sharedMatch[0];
    }
    __syncthreads();
}

__device__ void calculateParamsDC_TJ(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int * frac, int &tot, int &match, int * pr){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            frac[ec]++, frac[et]++, tot++;
            if(ec>et){
                int temp=ec;
                ec=et,et=temp;
            }
            if(ec==et) match++;
            if(ec==0&&et==2) pr[0]++;
            else if(ec==0&&et==3) pr[1]++;
            else if(ec==1&&et==2) pr[2]++;
            else if(ec==1&&et==3) pr[3]++;
        }
    }
}

__device__ void calculateParamsDCParallel_TJ(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int * frac, int &tot, int &match, int * pr){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;

    __shared__ int sharedP[4][128];
    __shared__ int sharedFrac[4][128];
    __shared__ int sharedTot[128];
    __shared__ int sharedMatch[128];

    const bool active = tx < 128;
    if (active) {
        for(int i=0;i<4;i++) sharedP[i][tx]=0;
        for(int i=0;i<4;i++) sharedFrac[i][tx]=0;
        sharedTot[tx]=0;
        sharedMatch[tx]=0;
    }

    __syncthreads();

    if (active) {
        for(int i=tx;i<compLen;i+=128){
            long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
            for(int j=0;j<16&&i*16+j<seqLen;j++){
                int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
                if(et>=4||ec>=4) continue;
                sharedFrac[ec][tx]++, sharedFrac[et][tx]++, sharedTot[tx]++;
                if(ec>et){
                    int temp=ec;
                    ec=et,et=temp;
                }
                if(ec==et) sharedMatch[tx]++;
                if(ec==0&&et==2) sharedP[0][tx]++;
                else if(ec==0&&et==3) sharedP[1][tx]++;
                else if(ec==1&&et==2) sharedP[2][tx]++;
                else if(ec==1&&et==3) sharedP[3][tx]++;
            }
        }
    }

    __syncthreads();
    // reduction
    
    for(int stride=128/2; stride>0; stride/=2){
        if(active && tx<stride){
            for(int i=0;i<4;i++) sharedP[i][tx] += sharedP[i][tx + stride];
            for(int i=0;i<4;i++) sharedFrac[i][tx] += sharedFrac[i][tx + stride];
            sharedTot[tx] += sharedTot[tx + stride];
            sharedMatch[tx] += sharedMatch[tx + stride];
        }
        __syncthreads();
    }

    if(tx==0){
        for(int i=0;i<4;i++) pr[i]=sharedP[i][0];
        for(int i=0;i<4;i++) frac[i]=sharedFrac[i][0];
        tot=sharedTot[0];
        match=sharedMatch[0];
    }
    __syncthreads();
}

__device__ void calculateParamsBatchDC_TJ(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, uint64_t * compressedSeqsConst, int * frac, int &tot, int &match, int * pr){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqsConst[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            frac[ec]++, frac[et]++, tot++;
            if(ec>et){
                int temp=ec;
                ec=et,et=temp;
            }
            if(ec==et) match++;
            if(ec==0&&et==2) pr[0]++;
            else if(ec==0&&et==3) pr[1]++;
            else if(ec==1&&et==2) pr[2]++;
            else if(ec==1&&et==3) pr[3]++;
        }
    }
}

__device__ void calculateParamsBatchDCParallel_TJ(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, uint64_t * compressedSeqsConst, int * frac, int &tot, int &match, int * pr){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;

    __shared__ int sharedP[4][128];
    __shared__ int sharedFrac[4][128];
    __shared__ int sharedTot[128];
    __shared__ int sharedMatch[128];
    const bool active = tx < 128;
    if (active) {
        for(int i=0;i<4;i++) sharedP[i][tx]=0;
        for(int i=0;i<4;i++) sharedFrac[i][tx]=0;
        sharedTot[tx]=0;
        sharedMatch[tx]=0;
    }

    __syncthreads();

    if (active) {
        for(int i=tx;i<compLen;i+=128){
            long long vt=compressedSeqs[px+i], vc=compressedSeqsConst[py+i];
            for(int j=0;j<16&&i*16+j<seqLen;j++){
                int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
                if(et>=4||ec>=4) continue;
                sharedFrac[ec][tx]++, sharedFrac[et][tx]++, sharedTot[tx]++;
                if(ec>et){
                    int temp=ec;
                    ec=et,et=temp;
                }
                if(ec==et) sharedMatch[tx]++;
                if(ec==0&&et==2) sharedP[0][tx]++;
                else if(ec==0&&et==3) sharedP[1][tx]++;
                else if(ec==1&&et==2) sharedP[2][tx]++;
                else if(ec==1&&et==3) sharedP[3][tx]++;
            }
        }
    }
    __syncthreads();

    // reduction
    for(int stride=128/2; stride>0; stride/=2){
        if(active && tx<stride){
            for(int i=0;i<4;i++) sharedP[i][tx] += sharedP[i][tx + stride];
            for(int i=0;i<4;i++) sharedFrac[i][tx] += sharedFrac[i][tx + stride];
            sharedTot[tx] += sharedTot[tx + stride];
            sharedMatch[tx] += sharedMatch[tx + stride];
        }
        __syncthreads();
    }
    if(tx==0){
        for(int i=0;i<4;i++) pr[i]=sharedP[i][0];
        for(int i=0;i<4;i++) frac[i]=sharedFrac[i][0];
        tot=sharedTot[0];
        match=sharedMatch[0];
    }
    __syncthreads();

}

__device__ void calculateParamsDC_K2P(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int &p, int &q, int &tot){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            tot++;
            if(et==ec) continue;
            if(et%2==ec%2) p++;
            else q++;
        }
    }
}

__device__ void calculateParamsDCParallel_K2P(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int &p, int &q, int &tot){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;

    __shared__ int sharedP[1024];
    __shared__ int sharedQ[1024];
    __shared__ int sharedTot[1024];

    sharedP[tx]=0;
    sharedQ[tx]=0;
    sharedTot[tx]=0;

    for(int i=tx;i<compLen;i+=1024){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            sharedTot[tx]++;
            if(et==ec) continue;
            if(et%2==ec%2) sharedP[tx]++;
            else sharedQ[tx]++;
        }
    }

    __syncthreads();

    //reduce
    for (int stride = 1024 / 2; stride > 0; stride /= 2) {
        if (tx < stride) {
            sharedP[tx] += sharedP[tx + stride];
            sharedQ[tx] += sharedQ[tx + stride];
            sharedTot[tx] += sharedTot[tx + stride];
        }
        __syncthreads();
    }
    // write the final results to the first thread
    if (tx == 0) {
        p = sharedP[0];
        q = sharedQ[0];
        tot = sharedTot[0];
    }
    __syncthreads();

}

__device__ void calculateParamsBatchDC_K2P(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, uint64_t * compressedSeqsConst, int &p, int &q, int &tot){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqsConst[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            tot++;
            if(et==ec) continue;
            if(et%2==ec%2) p++;
            else q++;
        }
    }
}

__device__ void calculateParamsBatchDCParallel_K2P(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, uint64_t * compressedSeqsConst, int &p, int &q, int &tot){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x; 
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    
    __shared__ int sharedP[1024];
    __shared__ int sharedQ[1024];
    __shared__ int sharedTot[1024];

    sharedP[tx]=0;
    sharedQ[tx]=0;
    sharedTot[tx]=0;
    
    for(int i=tx;i<compLen;i+=1024){
        long long vt=compressedSeqs[px+i], vc=compressedSeqsConst[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            sharedTot[tx]++;
            if(et==ec) continue;
            if(et%2==ec%2) sharedP[tx]++;
            else sharedQ[tx]++;
        }
    }

    __syncthreads();

    //reduce
    for (int stride = 1024 / 2; stride > 0; stride /= 2) {
        if (tx < stride) {
            sharedP[tx] += sharedP[tx + stride];
            sharedQ[tx] += sharedQ[tx + stride];
            sharedTot[tx] += sharedTot[tx + stride];
        }
        __syncthreads();
    }
    // write the final results to the first thread
    if (tx == 0) {
        p = sharedP[0];
        q = sharedQ[0];
        tot = sharedTot[0];
    }
    __syncthreads();

}

__device__ void calculateParamsDC_TAMURA(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int &p, int &q, int &tot, int &gc1, int &gc2){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            tot++;
            if(et==ec) continue;
            if(et%2==ec%2) p++;
            else q++;
            if(ec==1||ec==2) gc1++;
            if(et==1||et==2) gc2++;
        }
    }
}

__device__ void calculateParamsDCParallel_TAMURA(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int &p, int &q, int &tot, int &gc1, int &gc2){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;

    // create a shared memory array to store results
    __shared__ int sharedP[512];
    __shared__ int sharedQ[512];
    __shared__ int sharedTot[512];
    __shared__ int sharedGC1[512];
    __shared__ int sharedGC2[512];
    const bool active = tx < 512;
    if (active) {
        sharedP[tx] = 0;
        sharedQ[tx] = 0;
        sharedTot[tx] = 0;
        sharedGC1[tx] = 0;
        sharedGC2[tx] = 0;
    }

    __syncthreads();

    if (active) {
        for(int i=tx;i<compLen;i+=512){
            long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
            for(int j=0;j<16&&i*16+j<seqLen;j++){
                int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
                if(et>=4||ec>=4) continue;
                sharedTot[tx]++;
                if(et==ec) continue;
                if(et%2==ec%2) sharedP[tx]++;
                else sharedQ[tx]++;
                if(ec==1||ec==2) sharedGC1[tx]++;
                if(et==1||et==2) sharedGC2[tx]++;
            }
        }
    }

    __syncthreads();
    // reduce the results in shared memory
    for (int stride = 512 / 2; stride > 0; stride /= 2) {
        if (active && tx < stride) {
            sharedP[tx] += sharedP[tx + stride];
            sharedQ[tx] += sharedQ[tx + stride];
            sharedTot[tx] += sharedTot[tx + stride];
            sharedGC1[tx] += sharedGC1[tx + stride];
            sharedGC2[tx] += sharedGC2[tx + stride];
        }
        __syncthreads();
    }
    // write the final results to the first thread
    if (tx == 0) {
        p = sharedP[0];
        q = sharedQ[0];
        tot = sharedTot[0];
        gc1 = sharedGC1[0];
        gc2 = sharedGC2[0];
    }
    __syncthreads();
}

__device__ void calculateParamsBatchDC_TAMURA(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, uint64_t * compressedSeqsConst, int &p, int &q, int &tot, int &gc1, int &gc2){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqsConst[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            tot++;
            if(et==ec) continue;
            if(et%2==ec%2) p++;
            else q++;
            if(ec==1||ec==2) gc1++;
            if(et==1||et==2) gc2++;
        }
    }
}

__device__ void calculateParamsBatchDCParallel_TAMURA(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, uint64_t * compressedSeqsConst, int &p, int &q, int &tot, int &gc1, int &gc2){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;

    // create a shared memory array to store results
    __shared__ int sharedP[512];
    __shared__ int sharedQ[512];
    __shared__ int sharedTot[512];
    __shared__ int sharedGC1[512];
    __shared__ int sharedGC2[512];
    const bool active = tx < 512;
    if (active) {
        sharedP[tx] = 0;
        sharedQ[tx] = 0;
        sharedTot[tx] = 0;
        sharedGC1[tx] = 0;
        sharedGC2[tx] = 0;
    }

    __syncthreads();

    if (active) {
        for(int i=tx;i<compLen;i+=512){
            long long vt=compressedSeqs[px+i], vc=compressedSeqsConst[py+i];
            for(int j=0;j<16&&i*16+j<seqLen;j++){
                int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
                if(et>=4||ec>=4) continue;
                sharedTot[tx]++;
                if(et==ec) continue;
                if(et%2==ec%2) sharedP[tx]++;
                else sharedQ[tx]++;
                if(ec==1||ec==2) sharedGC1[tx]++;
                if(et==1||et==2) sharedGC2[tx]++;
            }
        }
    }

    __syncthreads();
    // reduce the results in shared memory
    for (int stride = 512 / 2; stride > 0; stride /= 2) {
        if (active && tx < stride) {
            sharedP[tx] += sharedP[tx + stride];
            sharedQ[tx] += sharedQ[tx + stride];
            sharedTot[tx] += sharedTot[tx + stride];
            sharedGC1[tx] += sharedGC1[tx + stride];
            sharedGC2[tx] += sharedGC2[tx + stride];
        }
        __syncthreads();
    }
    // write the final results to the first thread
    if (tx == 0) {
        p = sharedP[0];
        q = sharedQ[0];
        tot = sharedTot[0];
        gc1 = sharedGC1[0];
        gc2 = sharedGC2[0];
    }
    __syncthreads();
}

__global__ void MSADistConstructionDC(
    int rowId,
    uint64_t * compressedSeqs,
    double * dist,
    int seqLen,
    int numSequences,
    int distanceType
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int idx=tx+bs*bx;
    if(idx>=rowId) return;
    if(distanceType==DIST_UNCORRECTED||distanceType==DIST_JUKESCANTOR){
        int useful=0, match=0;
        calculateParamsDC(rowId, idx, seqLen, compressedSeqs, useful, match);
        double uncor=1-double(match)/useful;
        if(distanceType==DIST_UNCORRECTED) dist[idx]=uncor;
        else dist[idx]=-0.75*log(1.0-uncor/0.75);
        // printf("%d %d %d %d\n",rowId, idx, match, useful);
    }
    else if(distanceType==DIST_TAJIMANEI){
        int frac[4]={},pr[4]={},tot=0,match=0;
        double fr[4]={};
        calculateParamsDC_TJ(rowId, idx, seqLen, compressedSeqs, frac, tot, match, pr);
        for(int i=0;i<4;i++) fr[i]=double(frac[i])/tot/2.0;
        double h=0;
        h+=0.5*pr[0]*fr[0]*fr[2];
        h+=0.5*pr[1]*fr[0]*fr[3];
        h+=0.5*pr[2]*fr[1]*fr[2];
        h+=0.5*pr[3]*fr[1]*fr[3];
        double D=double(tot-match)/tot;
        double b=0.5*(1.0-fr[0]*fr[0]-fr[2]*fr[2]+D*D/h);
        dist[idx]=-b*log(1.0-D/b);
    }
    else if(distanceType==DIST_KIMURA2P||distanceType==DIST_JINNEI){
        int p=0,q=0,tot=0;
        calculateParamsDC_K2P(rowId, idx, seqLen, compressedSeqs, p, q, tot);
        double pp=double(p)/tot,qq=double(q)/tot;
        if(distanceType==DIST_KIMURA2P) dist[idx]=-0.5*log((1-2*pp-qq)*sqrt(1-2*qq));
        else dist[idx]=0.5*(1.0/(1-2*pp-qq)+0.5/(1-qq*2)-1.5);
    }
    else if(distanceType==DIST_TAMURA){
        int p=0,q=0,tot=0,gc1=0,gc2=0;
        calculateParamsDC_TAMURA(rowId, idx, seqLen, compressedSeqs, p, q, tot, gc1, gc2);
        double pp=double(p)/tot,qq=double(q)/tot, c=double(gc1)/tot+double(gc2)/tot-2*double(gc1)*double(gc2)/tot/tot;
        dist[idx]=-c*log(1-pp/c-qq)-0.5*(1-c)*log(1-2*qq);
    }
    else dist[idx]=0.0;
}


__global__ void MSADistConstructionRangeDC(
    int rowId,
    uint64_t * compressedSeqs,
    double * dist,
    int seqLen,
    int numSequences,
    int distanceType,
    int st,
    int ed
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    // int idx=tx+bs*bx;
    // if(idx>ed-st) return;
    // idx+=st;
    int idx;
    for(int blockID = bx; blockID <= ed-st; blockID += gridDim.x) {
        if (blockID > ed-st) return;
        idx = blockID+st;
    
        if(distanceType==DIST_UNCORRECTED||distanceType==DIST_JUKESCANTOR){
            int useful=0, match=0;
            calculateParamsDCParallel(rowId, idx, seqLen, compressedSeqs, useful, match);
            if (tx == 0) {
                double uncor=1-double(match)/useful;
                if(distanceType==DIST_UNCORRECTED) dist[idx]=uncor;
                else dist[idx]=-0.75*log(1.0-uncor/0.75);
            }
            // printf("%d %d %d %d\n",rowId, idx, match, useful);
        }
        else if(distanceType==DIST_TAJIMANEI){
            int frac[4]={},pr[4]={},tot=0,match=0;
            double fr[4]={};
            calculateParamsDCParallel_TJ(rowId, idx, seqLen, compressedSeqs, frac, tot, match, pr);
            if (tx == 0){
                for(int i=0;i<4;i++) fr[i]=double(frac[i])/tot/2.0;
                double h=0;
                h+=0.5*pr[0]*fr[0]*fr[2];
                h+=0.5*pr[1]*fr[0]*fr[3];
                h+=0.5*pr[2]*fr[1]*fr[2];
                h+=0.5*pr[3]*fr[1]*fr[3];
                double D=double(tot-match)/tot;
                double b=0.5*(1.0-fr[0]*fr[0]-fr[2]*fr[2]+D*D/h);
                dist[idx]=-b*log(1.0-D/b);
            }
        }
        else if(distanceType==DIST_KIMURA2P||distanceType==DIST_JINNEI){
            int p=0,q=0,tot=0;
            calculateParamsDCParallel_K2P(rowId, idx, seqLen, compressedSeqs, p, q, tot);
            if (tx == 0){
                double pp=double(p)/tot,qq=double(q)/tot;
                if(distanceType==DIST_KIMURA2P) dist[idx]=-0.5*log((1-2*pp-qq)*sqrt(1-2*qq));
                else dist[idx]=0.5*(1.0/(1-2*pp-qq)+0.5/(1-qq*2)-1.5);
            }
        }
        else if(distanceType==DIST_TAMURA){
            int p=0,q=0,tot=0,gc1=0,gc2=0;
            calculateParamsDCParallel_TAMURA(rowId, idx, seqLen, compressedSeqs, p, q, tot, gc1, gc2);
            if (tx == 0) {
                double pp=double(p)/tot,qq=double(q)/tot, c=double(gc1)/tot+double(gc2)/tot-2*double(gc1)*double(gc2)/tot/tot;
                dist[idx]=-c*log(1-pp/c-qq)-0.5*(1-c)*log(1-2*qq);
            }
        }
        else dist[idx]=0.0;
    }
}

__global__ void MSADistConstructionRangeForClusteringDC(
    int rowId,
    uint64_t * compressedSeqs,
    uint64_t * compressedSeqsConst,
    double * dist,
    int seqLen,
    int numSequences,
    int distanceType,
    int st,
    int ed
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    // int idx=tx+bs*bx;
    // if(idx>=ed-st) return;
    // idx+=st;
    int idx;
    for(int blockID = bx; blockID <= ed-st; blockID += gridDim.x) {
        if (blockID > ed-st) return;
        idx = blockID+st;
        if(distanceType==DIST_UNCORRECTED||distanceType==DIST_JUKESCANTOR){
            int useful=0, match=0;
            calculateParamsBatchDCParallel(rowId, idx, seqLen, compressedSeqs, compressedSeqsConst,useful, match);
            if (tx == 0) {
                double uncor=1-double(match)/useful;
                if(distanceType==DIST_UNCORRECTED) dist[idx]=uncor;
                else dist[idx]=-0.75*log(1.0-uncor/0.75);
            }// printf("%d %d %d %d\n",rowId, idx, match, useful);
        }
        else if(distanceType==DIST_TAJIMANEI){
            int frac[4]={},pr[4]={},tot=0,match=0;
            double fr[4]={};
            calculateParamsBatchDCParallel_TJ(rowId, idx, seqLen, compressedSeqs, compressedSeqsConst, frac, tot, match, pr);
            if (tx == 0) {
                for(int i=0;i<4;i++) fr[i]=double(frac[i])/tot/2.0;
                double h=0;
                h+=0.5*pr[0]*fr[0]*fr[2];
                h+=0.5*pr[1]*fr[0]*fr[3];
                h+=0.5*pr[2]*fr[1]*fr[2];
                h+=0.5*pr[3]*fr[1]*fr[3];
                double D=double(tot-match)/tot;
                double b=0.5*(1.0-fr[0]*fr[0]-fr[2]*fr[2]+D*D/h);
                dist[idx]=-b*log(1.0-D/b);
            }
        }
        else if(distanceType==DIST_KIMURA2P||distanceType==DIST_JINNEI){
            int p=0,q=0,tot=0;
            calculateParamsBatchDCParallel_K2P(rowId, idx, seqLen, compressedSeqs, compressedSeqsConst, p, q, tot);
            if (tx == 0) {
                double pp=double(p)/tot,qq=double(q)/tot;
                if(distanceType==DIST_KIMURA2P) dist[idx]=-0.5*log((1-2*pp-qq)*sqrt(1-2*qq));
                else dist[idx]=0.5*(1.0/(1-2*pp-qq)+0.5/(1-qq*2)-1.5);
        
            }
        }
        else if(distanceType==DIST_TAMURA){
            int p=0,q=0,tot=0,gc1=0,gc2=0;
            calculateParamsBatchDCParallel_TAMURA(rowId, idx, seqLen, compressedSeqs, compressedSeqsConst, p, q, tot, gc1, gc2);
            if (tx == 0) {
                double pp=double(p)/tot,qq=double(q)/tot, c=double(gc1)/tot+double(gc2)/tot-2*double(gc1)*double(gc2)/tot/tot;
                dist[idx]=-c*log(1-pp/c-qq)-0.5*(1-c)*log(1-2*qq);
        
                }
        }
        else dist[idx]=0.0;
    }
}

__global__ void MSADistConstructionSpecialIDDC(
    int rowId,
    uint64_t * compressedSeqsBackbone,
    uint64_t * compressedSeqsConst,
    double * dist,
    int seqLen,
    int backboneSize,
    int distanceType,
    int numToConstruct,
    int * d_id,
    int * d_leafMap
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int idx=tx+bs*bx;
    if(idx>=numToConstruct) return;
    int mapIdx=d_leafMap[idx];
    idx = d_id[idx];
    if(idx==-1) return;

    int idx_const = idx;
    if (idx >= backboneSize) idx_const = mapIdx;
    uint64_t * compressedSeqs = compressedSeqsBackbone;
    if (idx >= backboneSize) compressedSeqs = compressedSeqsConst;

    if(distanceType==DIST_UNCORRECTED||distanceType==DIST_JUKESCANTOR){
        int useful=0, match=0;
        calculateParamsBatchDC(rowId, idx_const, seqLen, compressedSeqs, compressedSeqsConst,useful, match);
        double uncor=1-double(match)/useful;
        if(distanceType==DIST_UNCORRECTED) dist[idx]=uncor;
        else dist[idx]=-0.75*log(1.0-uncor/0.75);
        // printf("%d %d %d %d\n",rowId, idx, match, useful);
    }
    else if(distanceType==DIST_TAJIMANEI){
        int frac[4]={},pr[4]={},tot=0,match=0;
        double fr[4]={};
        calculateParamsBatchDC_TJ(rowId, idx_const, seqLen, compressedSeqs, compressedSeqsConst, frac, tot, match, pr);
        for(int i=0;i<4;i++) fr[i]=double(frac[i])/tot/2.0;
        double h=0;
        h+=0.5*pr[0]*fr[0]*fr[2];
        h+=0.5*pr[1]*fr[0]*fr[3];
        h+=0.5*pr[2]*fr[1]*fr[2];
        h+=0.5*pr[3]*fr[1]*fr[3];
        double D=double(tot-match)/tot;
        double b=0.5*(1.0-fr[0]*fr[0]-fr[2]*fr[2]+D*D/h);
        dist[idx]=-b*log(1.0-D/b);
    }
    else if(distanceType==DIST_KIMURA2P||distanceType==DIST_JINNEI){
        int p=0,q=0,tot=0;
        calculateParamsBatchDC_K2P(rowId, idx_const, seqLen, compressedSeqs, compressedSeqsConst, p, q, tot);
        double pp=double(p)/tot,qq=double(q)/tot;
        if(distanceType==DIST_KIMURA2P) dist[idx]=-0.5*log((1-2*pp-qq)*sqrt(1-2*qq));
        else dist[idx]=0.5*(1.0/(1-2*pp-qq)+0.5/(1-qq*2)-1.5);
    }
    else if(distanceType==DIST_TAMURA){
        int p=0,q=0,tot=0,gc1=0,gc2=0;
        calculateParamsBatchDC_TAMURA(rowId, idx_const, seqLen, compressedSeqs, compressedSeqsConst, p, q, tot, gc1, gc2);
        double pp=double(p)/tot,qq=double(q)/tot, c=double(gc1)/tot+double(gc2)/tot-2*double(gc1)*double(gc2)/tot/tot;
        dist[idx]=-c*log(1-pp/c-qq)-0.5*(1-c)*log(1-2*qq);
    }
    else dist[idx]=0.0;
}

void MashPlacement::MSADeviceArraysDC::distRangeConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist, int l, int r, bool clustering) const{
    int threadNum = 1024, blockNum = 1024;
    if (!clustering) {
        MSADistConstructionRangeDC <<<1024, 1024>>>  (
            rowId, 
            d_compressedSeqsBackBone, 
            d_mashDist, 
            d_seqLen,
            backboneSize,
            params.distanceType,
            l,
            r
        );
    } else {
        MSADistConstructionRangeForClusteringDC <<<1024, 1024>>> (
            rowId, 
            d_compressedSeqsBackBone, 
            d_compressedSeqsConst,
            d_mashDist, 
            d_seqLen,
            backboneSize,
            params.distanceType,
            l,
            r
        );
    }
}

void MashPlacement::MSADeviceArraysDC::distConstructionOnGpuForBackboneDC(Param& params, int rowId, double* d_mashDist) const{
    int threadNum = 1024, blockNum = 1024;
    MSADistConstructionDC<<<1024, 1024>>>  (
        rowId, 
        d_compressedSeqsBackBone, 
        d_mashDist, 
        d_seqLen,
        backboneSize,
        params.distanceType
    );
}

void MashPlacement::MSADeviceArraysDC::distConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist) const{
    int threadNum = 1024, blockNum = 1024;
    MSADistConstructionDC<<<1024, 1024>>> (
        rowId, 
        d_compressedSeqsBackBone, 
        d_mashDist, 
        d_seqLen,
        backboneSize,
        params.distanceType
    );
}

void MashPlacement::MSADeviceArraysDC::distSpecialIDConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist, int numToConstruct, int* d_id, int * d_leafMap) const{
    int threadNum = 1024, blockNum = 1024;
    // std::cerr << "rowId: " << rowId << ", params.distanceType: " << params.distanceType << ", numToConstruct: " << numToConstruct << std::endl;
    MSADistConstructionSpecialIDDC <<<1024, 1024>>> (
        rowId, 
        d_compressedSeqsBackBone,
        d_compressedSeqsConst, 
        d_mashDist, 
        d_seqLen,
        backboneSize,
        params.distanceType,
        numToConstruct,
        d_id,
        d_leafMap
    );
}
