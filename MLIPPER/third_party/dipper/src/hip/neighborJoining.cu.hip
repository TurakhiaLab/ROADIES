#ifndef NJ_CUH
#include "mash_placement.cuh.hip"
#endif

#include <stdio.h>
#include <ostream>
#include <fstream>
#include <queue>
#include <vector>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/tuple.h>
#include <hip/hip_runtime.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <cmath>

__global__ void fillDismatrix(int numSequences, double* d_mashDist){
    int tx=threadIdx.x, bx=blockIdx.x;
    int bs=blockDim.x, gs=gridDim.x;
    int st=bx, step=gs;
    for(int i=st;i<numSequences;i+=step){
        for(int j=tx;j<numSequences;j+=bs){
            if(i==j) d_mashDist[i*numSequences+j]=0; // Distances on diagnal is set to 0
            if(i<=j) continue; // Only do the following process when we are in lower triangle
            double val=d_mashDist[i*numSequences+j];
            d_mashDist[j*numSequences+i]=val;
        }
    }
}


void MashPlacement::NJDeviceArrays::getDismatrix(
    int numSequences,
    Param& params,
    const MashDeviceArrays& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArrays& msaDeviceArrays
){
    //Allocate memories on GPU and copy data to GPU
    d_numSequences = numSequences;
    auto err = hipMalloc(&d_mashDist, 1ll*numSequences*numSequences*sizeof(double));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }
    err = hipMalloc(&d_U, numSequences*sizeof(double));
    if (err != hipSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: hipMalloc failed!\n");
        exit(1);
    }
    hipDeviceSynchronize();
    if(params.in == "d"){
        matrixReader.distConstructionOnGpu(params, 0, d_mashDist);
    }
    for(int i=1;i<numSequences;i++){
        if(params.in == "r"){
            mashDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_mashDist+i*numSequences
            );
        }
        else if(params.in == "d"){
            matrixReader.distConstructionOnGpu(
                params,
                i,
                d_mashDist+i*numSequences
            );
        }
        else if(params.in == "m"){
            msaDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_mashDist+i*numSequences
            );
        }
    }
    fillDismatrix <<<1024,1024>>> (d_numSequences, d_mashDist);
}

void MashPlacement::NJDeviceArrays::deallocateDeviceArrays(){
    //Free space on GPU at the end
    hipFree(d_mashDist);
    hipFree(d_U);
    hipDeviceSynchronize();
}

__global__ void calculateU(int d_numSequences, double *d_mashDist, double *d_U){
    //Only executes once, calculate the sum of distances in each row
    int numSequences=d_numSequences;
    __shared__ double temp_U[128]; // Make sure to set this value larger than numSequences/gridSize
    for(auto i=0;i<128;i++)
        if(threadIdx.x==0)
            temp_U[i]=0;
    __syncthreads();
    for(auto i=blockIdx.x;i<numSequences;i+=gridDim.x){
        double temp=0;
        for(auto j=threadIdx.x;j<numSequences;j+=blockDim.x)
            if(i!=j) temp += d_mashDist[i*numSequences+j];
        atomicAdd(&temp_U[i/gridDim.x],temp);
        // For each block, calculate the sum of distances in some rows
        // For each thread, use atomicAdd to update the sum to shared memory
    }
    __syncthreads();
    if(threadIdx.x==0){
        for(auto i=blockIdx.x;i<numSequences;i+=gridDim.x) d_U[i]=temp_U[i/gridDim.x];
        // For each block, update the sum to global memory
    }
}

__global__ void findMinDist(int d_numSequences, double *d_mashDist, double *d_U, thrust::tuple<int,int,double> *minPos, int len){
    __shared__ double U[128];
    //Required length of each block equal to d_numSequences/gridDim.x+1;
    //Make sure to set the size of shared memory large enough
    int tx=threadIdx.x, bx=blockIdx.x;
    int bs=blockDim.x, gs=gridDim.x;
    int numSequences=d_numSequences;
    int sz=numSequences/gs,st=sz*bx;
    if(numSequences%gs>bx) sz++;
    st+=min(bx,numSequences%gs);
    int ed=st+sz;
    // Calculate the start and end row for the block
    for(auto i=st+tx;i<ed;i+=bs)
        U[i-st]=d_U[i]/(numSequences-2);
    // Store the U value corresponding to the block in shared memory
    __syncthreads();
    double colU;
    double minD=10000;
    int x=0,y=0;
    for(auto j=tx;j<numSequences;j+=bs){
        colU=d_U[j]/(numSequences-2); // Store the U value corresponding to j-th column in register
        for(auto i=st;i<ed;i++){
            double temp=d_mashDist[i*len+j]-U[i-st]-colU;
            if(i!=j&&temp<minD)
                minD=temp,x=i,y=j; // Update if we find a smaller new value
        }
    }

    thrust::tuple <int,int,double> minTuple(x,y,minD);
    minPos[bx*bs+tx]=minTuple;
    //Store the minimum value found by this thread to global memory for further processing
}

struct compare_tuple
{
  __host__ __device__
  bool operator()(thrust::tuple<int,int,double> lhs, thrust::tuple<int,int,double> rhs)
  {
    return thrust::get<2>(lhs) < thrust::get<2>(rhs);
    //Always find the tuple whose third value (the criteria we want to minimize) is minimized
  }
};


__global__ void updateDisMatrix(int d_numSequences, double *d_mashDist, double disxy, int x,int y, double* d_U, int total){
    // Erase node x and y from the matrix (x<y)
    // Insert the distances to the new node into the matrix, 
    // at the original column/row of x for coalasced memory access.
    // Always move the last row/column (distances related to last node) to row/column of y to avoid gap in matrix
    int tx=threadIdx.x, bx=blockIdx.x;
    int bs=blockDim.x, gs=gridDim.x;
    int start=bx*bs+tx,step=bs*gs,numSequences=total;
    for(int i=start;i<d_numSequences-1;i+=step)
        if(i!=x&&i!=y){
            double val=(d_mashDist[x*numSequences+i]+d_mashDist[y*numSequences+i]-disxy)*0.5;
            // Calculate distance from i to the new node
            double fardis=d_mashDist[(d_numSequences-1)*numSequences+i];
            // Calcualte distance from i to the last node
            d_U[i]+=-d_mashDist[x*numSequences+i]-d_mashDist[y*numSequences+i]+val;
            atomicAdd(&d_U[x],val);
            // Update U values (sum of distances in each row)
            d_mashDist[x*numSequences+i]=val;
            d_mashDist[i*numSequences+x]=val;
            d_mashDist[y*numSequences+i]=fardis;
            d_mashDist[i*numSequences+y]=fardis;
            // Update distance matrix
        }
    if(tx==0&&bx==0){
        // Deal with the influence of distance between new node and last node
        d_U[y]=d_U[d_numSequences-1];
        int i=d_numSequences-1;
        double val=(d_mashDist[x*numSequences+i]+d_mashDist[y*numSequences+i]-disxy)*0.5;
        d_U[y]+=-d_mashDist[x*numSequences+i]-d_mashDist[y*numSequences+i]+val;
        atomicAdd(&d_U[x],val);
        d_mashDist[x*numSequences+y]=val;
        d_mashDist[y*numSequences+x]=val;
    }
}


void MashPlacement::NJDeviceArrays::findNeighbourJoiningTree(std::vector <std::string> &name, std::ofstream& output_){
    std::vector <std::vector<std::pair<int,double>>> tree(d_numSequences*2);
    // Store the tree in a vector of vectors, while each small vector consists of several pairs
    // First element in the pair is the index of its child
    // Second element in the pair is the distance to its child
    int cntThread=256,cntBlock=256; // Could be adjusted based on experiment
    int realID[d_numSequences]; // Actual sequence id corresponding to i-th row/column in the matrix
    for(int i=0;i<d_numSequences;i++) realID[i]=i;
    // Build the full distance matrix
    thrust::device_vector <thrust::tuple<int,int,double>> minPos(cntThread*cntBlock);
    // Device_vector that stores minimum value found by each thread
    calculateU <<<cntBlock,cntThread>>> (d_numSequences, d_mashDist, d_U);
    // Calculate initial sum of distances in each row
    int ID=d_numSequences;
    for(int i=0;i<d_numSequences-2;i++){
        findMinDist <<<cntBlock,cntThread>>> (d_numSequences-i,d_mashDist,d_U,thrust::raw_pointer_cast(minPos.data()),d_numSequences);
        hipDeviceSynchronize();
        auto iter=thrust::min_element(minPos.begin(),minPos.end(),compare_tuple());
        // Find the global minimum
        thrust::tuple<int,int,double> smallest=*iter;
        int x=thrust::get<0>(smallest),y=thrust::get<1>(smallest);
        double modified_dis=thrust::get<2>(smallest);
        if(x>y) std::swap(x,y);
        // Ensure x is smaller than y
        double* dis=new double;
        double *ux=new double;
        double *uy=new double;
        hipMemcpy(dis,d_mashDist+x*d_numSequences+y,sizeof(double),hipMemcpyDeviceToHost);
        hipMemcpy(ux,d_U+x,sizeof(double),hipMemcpyDeviceToHost);
        hipMemcpy(uy,d_U+y,sizeof(double),hipMemcpyDeviceToHost);
        double blX = (*dis+(*ux)/(d_numSequences-i-2)-(*uy)/(d_numSequences-i-2))*0.5;
        double blY = *dis - blX;
        // Calculate branch lengths
        if(blX<0) blY+=blX, blX=0;
        if(blY<0) blX+=blY, blY=0;
        // Avoid negative branches
        tree[ID].push_back({realID[x],blX});
        tree[ID].push_back({realID[y],blY});
        // Insert the two edges into the tree
        double temp=0;
        realID[x]=ID++, realID[y]=realID[d_numSequences-i-1]; 
        // Update the actual sequence indices corresponding to rows/columns in distance matrix
        hipMemcpy(d_U+x,&temp,sizeof(double),hipMemcpyHostToDevice);
        updateDisMatrix <<<cntBlock,cntThread>>> (d_numSequences-i,d_mashDist,*dis,x,y,d_U, d_numSequences);
        // Update the distance matrix and U array
        hipDeviceSynchronize();
    }
    // If only two nodes are left, their distance is fixed
    int x=0,y=1;
    double dis=0;
    hipMemcpy(&dis,d_mashDist+x*d_numSequences+y,sizeof(double),hipMemcpyDeviceToHost);
    tree[d_numSequences*2-2].push_back({realID[x],dis*0.5});
    tree[d_numSequences*2-2].push_back({realID[y],dis*0.5});
    hipDeviceSynchronize();
    const double kNjLenEps = 1e-12;
    auto collectFlatKids = [&](int ch, double L, std::vector<std::pair<int, double>>& out) {
        std::function<void(int, double)> go;
        go = [&](int c, double len) {
            if (std::fabs(len) <= kNjLenEps && !tree[c].empty()) {
                for (const auto& p : tree[c]) {
                    go(p.first, p.second);
                }
            } else {
                out.push_back({c, len});
            }
        };
        go(ch, L);
    };
    // Output the tree recursively
    std::function<void(int)> print = [&](int node) {
        if (tree[node].size()) {
            output_ << "(";
            std::vector<std::pair<int, double>> kids;
            if (g_printBinaryNewick) {
                kids = tree[node];
            } else {
                for (const auto& p : tree[node]) {
                    collectFlatKids(p.first, p.second, kids);
                }
            }
            for (size_t i = 0; i < kids.size(); ++i) {
                print(kids[i].first);
                output_ << ":";
                output_ << kids[i].second << (i + 1 == kids.size() ? ')' : ',');
            }
        } else {
            output_ << name[node];
        }
    };
    // Root of the tree has an index of d_numSequneces*2-2
    print(d_numSequences*2-2);
    // std::cout<<";\n";
    output_<<";\n";
}

