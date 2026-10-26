
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

void MashPlacement::MSADeviceArrays::allocateDeviceArrays(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params)
{
    hipError_t err;

    numSequences = int(num);
    seqLen = h_seqLengths[0];
    // Allocate memory
    err = hipMalloc(&d_seqLengths, 1ll*numSequences*sizeof(uint64_t));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }
    // std::cerr<<"????????\n";
    /* Flatten data */
    uint64_t flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) flatStringLength+= (h_seqLengths[i]+15)/16;
    uint64_t * h_flattenCompressSeqs = new uint64_t[flatStringLength];
    flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) 
    {
        uint64_t flatStringLengthLocal = (h_seqLengths[i]+15)/16;
        flatStringLength+=flatStringLengthLocal;
        for (size_t j=0; j<flatStringLengthLocal;j++)  
        {
            h_flattenCompressSeqs[j] = h_compressedSeqs[i][j];
            // if (i==9) printf("%u\n",h_flattenCompressSeqs[j]); 
        }
        h_flattenCompressSeqs += flatStringLengthLocal;
    }
    h_flattenCompressSeqs -= flatStringLength;
    // std::cerr<<"?????????\n";

    err = hipMalloc(&d_compressedSeqs, 1ll*flatStringLength*sizeof(uint64_t));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }


    // Transfer data

    err = hipMemcpy(d_seqLengths, h_seqLengths, 1ll*numSequences*sizeof(uint64_t), hipMemcpyHostToDevice);
    if (err != hipSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    err = hipMemcpy(d_compressedSeqs, h_flattenCompressSeqs, 1ll*flatStringLength*sizeof(uint64_t), hipMemcpyHostToDevice);
    if (err != hipSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    hipDeviceSynchronize();
}

void MashPlacement::MSADeviceArrays::deallocateDeviceArrays(){
    hipFree(d_compressedSeqs);
    hipFree(d_seqLengths);
    // cudaFree(d_hashList);
    // cudaFree(d_mashDist);
}

#define DIST_UNCORRECTED 1
#define DIST_JUKESCANTOR 2
#define DIST_TAJIMANEI 3
#define DIST_KIMURA2P 4
#define DIST_TAMURA 5
#define DIST_JINNEI 6


__device__ void calculateParams(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int & useful, int & match){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et<4||ec<4) useful++;
            if(et<4&&et==ec) match++;
        }
    }
}


__device__ void calculateParamsParallel(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int & useful, int & match){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;

    // create a shared memory array to store results
    __shared__ int sharedUseful[1024];
    __shared__ int sharedMatch[1024];
    sharedUseful[tx] = 0;
    sharedMatch[tx] = 0;

    if (tx >= compLen) {
        return; // If thread index is out of bounds, exit early
    }
    for (int i=tx; i<compLen; i+=1024) {
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et<4||ec<4) sharedUseful[tx]++;
            if(et<4&&et==ec) sharedMatch[tx]++;
        }
    }
    __syncthreads();

    // reduce the results in shared memory
    for (int stride = bs / 2; stride > 0; stride /= 2) {
        if (tx < stride) {
            sharedUseful[tx] += sharedUseful[tx + stride];
            sharedMatch[tx] += sharedMatch[tx + stride];
        }
        __syncthreads();
    }

    // write the final results to the first thread
    if (tx == 0) {
        // for (int i=0; i<1024; i+=1) {
        //     if (i > compLen) break;
        //     useful += sharedUseful[i];
        //     match += sharedMatch[i];
        // }
        useful = sharedUseful[0];
        match = sharedMatch[0];
    }
    __syncthreads();

    // for(int i=0;i<compLen;i++){
    //     long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
    //     for(int j=0;j<16&&i*16+j<seqLen;j++){
    //         int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
    //         if(et<4||ec<4) useful++;
    //         if(et<4&&et==ec) match++;
    //     }
    // }
}

__device__ void calculateParams_TJ(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int * frac, int &tot, int &match, int * pr){
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

__device__ void calculateParams_K2P(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int &p, int &q, int &tot){
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

__device__ void calculateParams_TAMURA(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int &p, int &q, int &tot, int &gc1, int &gc2){
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

__global__ void MSADistConstruction(
    int rowId,
    uint64_t * compressedSeqs,
    double * dist,
    int seqLen,
    int numSequences,
    int distanceType
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int idx=tx+bs*bx;
    // if(idx>=rowId) return;
    for (int blockID = bx; blockID < rowId; blockID += gridDim.x) {
        if (blockID >= rowId) return; // Ensure we don't access out of bounds
        // printf("bx: %d, rowId: %d\n", blockID, rowId);
        if(distanceType==DIST_UNCORRECTED||distanceType==DIST_JUKESCANTOR){
            int useful=0, match=0;
            calculateParamsParallel(rowId, blockID, seqLen, compressedSeqs, useful, match);
            // calculateParams(rowId, idx, seqLen, compressedSeqs, useful, match);
            if (tx == 0) {   
                double uncor=1-double(match)/useful;
                if(distanceType==DIST_UNCORRECTED) dist[blockID]=uncor;
                else dist[blockID]=-0.75*log(1.0-uncor/0.75);
                // printf("%d %d %d %d\n",rowId, blockID, match, useful);
            }
        }
        else if(distanceType==DIST_TAJIMANEI){
            int frac[4]={},pr[4]={},tot=0,match=0;
            double fr[4]={};
            calculateParams_TJ(rowId, idx, seqLen, compressedSeqs, frac, tot, match, pr);
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
            calculateParams_K2P(rowId, idx, seqLen, compressedSeqs, p, q, tot);
            double pp=double(p)/tot,qq=double(q)/tot;
            if(distanceType==DIST_KIMURA2P) dist[idx]=-0.5*log((1-2*pp-qq)*sqrt(1-2*qq));
            else dist[idx]=0.5*(1.0/(1-2*pp-qq)+0.5/(1-qq*2)-1.5);
        }
        else if(distanceType==DIST_TAMURA){
            int p=0,q=0,tot=0,gc1=0,gc2=0;
            calculateParams_TAMURA(rowId, idx, seqLen, compressedSeqs, p, q, tot, gc1, gc2);
            double pp=double(p)/tot,qq=double(q)/tot, c=double(gc1)/tot+double(gc2)/tot-2*double(gc1)*double(gc2)/tot/tot;
            dist[idx]=-c*log(1-pp/c-qq)-0.5*(1-c)*log(1-2*qq);
        }
        else dist[idx]=0.0;
    }
}


void MashPlacement::MSADeviceArrays::distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const{
    int threadNum = 1024, blockNum = 1024; // dont change threadNUM, interally it is used to calculate the distance
    // printf("rowId: %d params.distanceType %d \n", rowId, params.distanceType);
    MSADistConstruction <<<blockNum, threadNum>>> (
        rowId, 
        d_compressedSeqs, 
        d_mashDist, 
        seqLen,
        numSequences,
        params.distanceType
    );
}


