#include "mash_placement.cuh"

#include <stdio.h>
#include <cstdio>
#include <queue>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <chrono>
#include <iostream>
#include <cub/cub.cuh>

void MashPlacement::MatrixReader::allocateDeviceArrays(int num, FILE* fPtr){
    numSequences = num;
    h_dist = new double[num];
    name.resize(num);
    buffer = new char[num*20];
    filePtr = fPtr;
}

void MashPlacement::MatrixReader::distConstructionOnGpu(Param& params, int rowId, double* d_mashDist){
    fgets(buffer, numSequences*20, filePtr);
    char *p=buffer;
    name[rowId]="";
    while(1){
        char c=*p;
        p++;
        if(c=='\t'||c=='\n'||c==' ') break;
        name[rowId]+=c;
    }
    if(rowId == 0) return;
    for(int i=0;i<rowId;i++){
        std::string num="";
        while(1){
            char c=*p;
            p++;
            if(c=='\t'||c=='\n'||c==' ') break;
            num+=c;
        }
        h_dist[i]=stof(num);
    }
    cudaMemcpy(d_mashDist, h_dist, rowId*sizeof(double), cudaMemcpyHostToDevice);
}