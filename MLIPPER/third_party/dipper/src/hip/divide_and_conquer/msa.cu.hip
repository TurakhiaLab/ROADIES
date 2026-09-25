
#include <hip/hip_runtime.h>
#include "../mash_placement.cuh.hip"

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

void MashPlacement::MSADeviceArraysDC::transferToDeviceArraysDC(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, int gpuCluster, Param& params)
{
    hipError_t err;
    
    size_t maxLengthCompressed = (this->d_seqLen + 15) / 16;

    /* Flatten complete data */
    
    int start = gpuCluster*params.backboneSize;
    int end = (start + params.backboneSize <= num) ? (start + params.backboneSize) : num;
    uint64_t flatStringLength= maxLengthCompressed*(end-start);
    if (gpuCluster == 0)
        this->h_compressedSeqs = new uint64_t[flatStringLength];
    
    for (size_t i = start; i < end; i++) 
    {
        for (size_t j=0; j<maxLengthCompressed;j++)  
        {
            this->h_compressedSeqs[j] = h_compressedSeqs[i][j];
        }
        this->h_compressedSeqs += maxLengthCompressed;
    }
    this->h_compressedSeqs -= flatStringLength;


    /* Transfer only the backbone data */
    err = hipMemcpy(d_compressedSeqsConst, this->h_compressedSeqs, 1ll*(flatStringLength)*sizeof(uint64_t), hipMemcpyHostToDevice);
    if (err != hipSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    hipDeviceSynchronize();
}


void MashPlacement::MSADeviceArraysDC::allocateDeviceArraysDC(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params)
{
    hipError_t err;
    this->totalNumSequences = int(num);
    this->backboneSize = params.backboneSize;
    
    this->d_seqLen = h_seqLengths[0];
    
    size_t maxLengthCompressed = (this->d_seqLen + 15) / 16;

    err = hipMalloc(&d_compressedSeqsBackBone, maxLengthCompressed*this->backboneSize*sizeof(uint64_t));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }

    err = hipMalloc(&d_compressedSeqsConst, maxLengthCompressed*this->backboneSize*sizeof(uint64_t));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
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
    err = hipMemcpy(d_compressedSeqsBackBone, this->h_compressedSeqs, 1ll*(maxLengthCompressed*this->backboneSize)*sizeof(uint64_t), hipMemcpyHostToDevice);
    if (err != hipSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    hipDeviceSynchronize();
}

void MashPlacement::MSADeviceArraysDC::deallocateDeviceArraysDC(){
    hipFree(d_compressedSeqsBackBone);
    hipFree(d_compressedSeqsConst);
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
            if(et<4||ec<4) useful++;
            if(et<4&&et==ec) match++;
        }
    }
}

__device__ void calculateParamsBatchDC(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, uint64_t * compressedSeqsConst, int & useful, int & match){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    // printf("px: %lld, py: %lld\n", px, py);
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqsConst[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et<4||ec<4) useful++;
            if(et<4&&et==ec) match++;
        }
    }
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
    int idx=tx+bs*bx;
    if(idx>ed-st) return;
    idx+=st;
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
    int idx=tx+bs*bx;
    if(idx>=ed-st) return;
    idx+=st;
    if(distanceType==DIST_UNCORRECTED||distanceType==DIST_JUKESCANTOR){
        int useful=0, match=0;
        calculateParamsBatchDC(rowId, idx, seqLen, compressedSeqs, compressedSeqsConst,useful, match);
        double uncor=1-double(match)/useful;
        if(distanceType==DIST_UNCORRECTED) dist[idx]=uncor;
        else dist[idx]=-0.75*log(1.0-uncor/0.75);
        // printf("%d %d %d %d\n",rowId, idx, match, useful);
    }
    else if(distanceType==DIST_TAJIMANEI){
        int frac[4]={},pr[4]={},tot=0,match=0;
        double fr[4]={};
        calculateParamsBatchDC_TJ(rowId, idx, seqLen, compressedSeqs, compressedSeqsConst, frac, tot, match, pr);
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
        calculateParamsBatchDC_K2P(rowId, idx, seqLen, compressedSeqs, compressedSeqsConst, p, q, tot);
        double pp=double(p)/tot,qq=double(q)/tot;
        if(distanceType==DIST_KIMURA2P) dist[idx]=-0.5*log((1-2*pp-qq)*sqrt(1-2*qq));
        else dist[idx]=0.5*(1.0/(1-2*pp-qq)+0.5/(1-qq*2)-1.5);
    }
    else if(distanceType==DIST_TAMURA){
        int p=0,q=0,tot=0,gc1=0,gc2=0;
        calculateParamsBatchDC_TAMURA(rowId, idx, seqLen, compressedSeqs, compressedSeqsConst, p, q, tot, gc1, gc2);
        double pp=double(p)/tot,qq=double(q)/tot, c=double(gc1)/tot+double(gc2)/tot-2*double(gc1)*double(gc2)/tot/tot;
        dist[idx]=-c*log(1-pp/c-qq)-0.5*(1-c)*log(1-2*qq);
    }
    else dist[idx]=0.0;
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
    // printf("rowId %d idx %d mapIdx %d idx_new %d numSequences %d\n", rowId, idx, d_leafMap[idx], d_id[idx], backboneSize);
    int mapIdx=d_leafMap[idx];
    idx = d_id[idx];
    if(idx==-1) return;

    int idx_const = idx;
    if (idx > backboneSize) idx_const = mapIdx;
    uint64_t * compressedSeqs = compressedSeqsBackbone;
    if (idx > backboneSize) compressedSeqs = compressedSeqsConst;
    // printf("idx: %d, idx_const: %d, mapIdx: %d compressedSeqs: %p compressedSeqsConst: %p\n", idx, idx_const, mapIdx, compressedSeqsBackbone, compressedSeqsConst);
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
    int threadNum = 256, blockNum = (rowId+threadNum-1)/threadNum;
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
    int threadNum = 256, blockNum = (rowId+threadNum-1)/threadNum;
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
    int threadNum = 256, blockNum = (rowId+threadNum-1)/threadNum;
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
    int threadNum = 256, blockNum = (rowId+threadNum-1)/threadNum;
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


