/**
 * divide_and_conquer/placement_close_k.cu
 *
 * Divide-and-conquer placement: build backbone tree on subset, cluster non-backbone
 * by k-closest edge, then build per-cluster trees (edge/leaf masks, special-ID
 * distance construction). Handles MASH and MSA input, batch clustering, and
 * cluster-tree output.
 */

#include "../mash_placement.cuh"

#include <stdio.h>
#include <zlib.h>
#include <queue>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <chrono>
#include <cmath>
#include <iostream>
#include <cub/cub.cuh>
#include <cuda_runtime.h>

/** Allocate DC placement arrays: adjacency, closest_id/dis (and cluster copies), dist, len. */
void MashPlacement::KPlacementDeviceArraysDC::allocateDeviceArraysDC(size_t num, size_t totalNum, int gpuNum){
    cudaError_t err;

    numSequences = int(num);
    totalNumSequences = int(totalNum);
    bd = 2, idx = 0;

    err = cudaMalloc(&d_dist, totalNumSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_head, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_e, totalNumSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_len, totalNumSequences*8*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_nxt, totalNumSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_belong, totalNumSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_dis, totalNumSequences*40*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_id, totalNumSequences*40*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_closest_dis_cluster, totalNumSequences*40*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_id_cluster, totalNumSequences*40*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    cudaDeviceSynchronize();
}

/** Reset adjacency and closest_id/closest_dis to sentinels; used in DC backbone init. */
__global__ void initializeDC(
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
    for (int idx_=idx; idx_<lim; idx_+=gs*bs) {
        if(idx_>=lim) return;
        for(int i=0;i<5;i++){
            d_closest_dis[idx_*5+i]=2;
            d_closest_id[idx_*5+i]=-1;
        }
        nxt[idx_] = -1;
        e[idx_] = -1;
        belong[idx_] = -1;
        if(idx_<nodes) head[idx_] = -1;
    }
}

/** Compare placement tuples by addLen (third element). */
struct compare_tupleDC {
  __host__ __device__
  bool operator()(thrust::tuple<int,double,double> lhs, thrust::tuple<int,double,double> rhs)
  {
    return thrust::get<2>(lhs) < thrust::get<2>(rhs);
  }
};

/** For each candidate edge, compute (eid, fracLen, addLen) from closest_id/closest_dis. */
__global__ void calculateBranchLengthDC(
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
    bool penalize_negative_pendant
){
    int tx=threadIdx.x,bs=blockDim.x,bx=blockIdx.x,gs=gridDim.x;
    int idx=tx+bs*bx;
    // if(idx>=lim) return;
    for (int idx_=idx; idx_<lim; idx_+=gs*bs) {
        if(idx_>=num*4-4||belong[idx_]<e[idx_]){
            thrust::tuple <int,double,double> minTuple(0,0,2);
            minPos[idx_]=minTuple;
            continue;
        }
        int x=belong[idx_],oth=e[idx_];
        int eid=idx_,otheid;
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
        const double routing_deficit = additional_dis < 0
            ? -additional_dis
            : additional_dis;
        if(additional_dis<0) additional_dis=0;
        dis1-=additional_dis,dis2-=additional_dis;
        if(dis1<0) dis1=0;
        if(dis2<0) dis2=0;
        if(dis1>len[eid]) additional_dis+=dis1-len[eid],dis1=len[eid];
        if(dis2>len[eid]) additional_dis+=dis2-len[eid],dis2=len[eid];
        // assert(dis1+dis2-1e-6<=len[eid]);
        double rest=len[eid]-dis1-dis2;
        dis1+=rest/2,dis2+=rest/2;
        thrust::tuple <int,double,double> minTuple(
            eid, dis1,
            penalize_negative_pendant ? routing_deficit : additional_dis);
        minPos[idx_]=minTuple;
    }
}

/** Same as calculateBranchLengthDC but over edge mask d_edgeMask (subset of edges). */
__global__ void calculateBranchLengthSpecialIDDC(
    int num, // useless here
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
    int numToCalculate,
    int * d_edgeMask 
){
    int tx=threadIdx.x,bs=blockDim.x,bx=blockIdx.x,gs=gridDim.x;
    int idx=tx+bs*bx;
    // if(idx>=numToCalculate) return;
    for (int slot=idx; slot<numToCalculate; slot+=gs*bs) {
        const int edge_idx = d_edgeMask[slot];
        if(belong[edge_idx]<e[edge_idx]){
            thrust::tuple <int,double,double> minTuple(0,0,2);
            minPos[slot]=minTuple;
            continue;
        }
        int x=belong[edge_idx],oth=e[edge_idx];
        int eid=edge_idx,otheid;
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
        minPos[slot]=minTuple;
    }
}

/** BFS from x; update closest_id/closest_dis per edge when improving. */
__global__ void updateClosestNodesDC(
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
    id[++r]=x,dis[r]=0,from[r]=-1;
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

/** BFS from x restricted by d_edgeMaskIndex (cluster edges only). */
__global__ void updateClosestNodesClusterDC(
    int * head,
    int * nxt,
    int * e,
    double * len,
    double * closest_dis,
    int * closest_id,
    int x,
    int * id,
    int * from,
    double * dis,
    int * d_edgeMaskIndex
){
    int l=0,r=-1;
    id[++r]=x,dis[r]=0,from[r]=-1;
    while(l<=r){
        int node=id[l],fb=from[l];
        double d=dis[l];
        l++;
        for(int i=head[node];i!=-1;i=nxt[i]){
            if(d_edgeMaskIndex[i]!=i) continue;
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

/** BFS from x within cluster: skip cluster edge endpoints, use d_edgeMaskIndex. */
__global__ void updateClosestNodesInClusterDC(
    int * head,
    int * nxt,
    int * e,
    double * len,
    double * closest_dis,
    int * closest_id,
    int x,
    int * id,
    int * from,
    double * dis,
    int cluster_eid,
    int * belong,
    int * d_edgeMaskIndex
){
    int l=0,r=-1;
    id[++r]=x,dis[r]=0,from[r]=-1;
    int ed1=e[cluster_eid], ed2=belong[cluster_eid];
    while(l<=r){
        int node=id[l],fb=from[l];
        double d=dis[l];
        l++;
        if(node==ed1||node==ed2) continue;
        
        for(int i=head[node];i!=-1;i=nxt[i]){
            // printf("node: %d, head[node]: %d, belong[node]: %d i : %d \n", node, head[node], belong[node], i);
            if (d_edgeMaskIndex[i]!=i) continue;
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
__global__ void updateTreeStructureDC(
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
    int placeId, // Id of the newly placed node
    int edgeCount, // Position to insert a new edge in linked list
    int totalNumSequences
){
    int middle=placeId+totalNumSequences-1, outside=placeId;
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

/** Same as updateTreeStructureDC but uses placeCount for middle node indexing in cluster. */
__global__ void updateTreeStructureInClusterDC(
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
    int placeId, // Id of the newly placed node
    int edgeCount, // Position to insert a new edge in linked list
    int totalNumSequences,
    int placeCount // this is the placeCount-th leave
){
    int middle=placeCount+totalNumSequences-1, outside=placeId;
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

/** Build 3-node backbone (0, 1, nv) with two edges of length dis[0]/2. */
__global__ void buildInitialTreeDC(
    int totalNumSequences,
    int * head,
    int * e,
    double * len,
    int * nxt,
    int * belong,
    double * dis,
    int edgeCount
){
    int nv = totalNumSequences;
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
    // nv -> 1
    e[edgeCount]=1,len[edgeCount]=d/2,nxt[edgeCount]=head[nv],head[nv]=edgeCount,belong[edgeCount]=nv;
    edgeCount++;
}

/** Record new leaf and its edges in leaf/edge masks and leafMap for cluster subtree. */
__global__ void updateClusterInfoDC (
    int leafID,
    int edgeidx,
    int * d_leafMask,
    int * d_edgeMask,
    int * d_edgeMaskIndex,
    int edgeCount,
    int leafCount,
    int * d_leafMap,
    int leaf_idx_in_cluster
){

    d_leafMap[leafCount]=leaf_idx_in_cluster;
    d_leafMask[leafCount++]=leafID;
    for(int i=1;i<=4;i++) {
        d_edgeMask[edgeCount++]=edgeidx-i;
        d_edgeMaskIndex[edgeidx-i]=edgeidx-i;
    }

}

/** Copy closest_id/closest_dis to cluster copies before cluster placement. */
__global__ void copyClosestIdsDC(
    int * d_closest_id,
    int * d_closest_id_cluster,
    double * d_closest_dis,
    double * d_closest_dis_cluster,
    int edgeCount
){
    for(int i=0;i<edgeCount;i++){
        for(int j=0;j<5;j++){
            d_closest_id_cluster[i*5+j]=d_closest_id[i*5+j];
            d_closest_dis_cluster[i*5+j]=d_closest_dis[i*5+j];
        }
    }
}

/** Copy cluster closest_id/dis back to main arrays after cluster placement. */
__global__ void copyBackClosestIdsDC(
    int * d_closest_id,
    int * d_closest_id_cluster,
    double * d_closest_dis,
    double * d_closest_dis_cluster,
    int edgeCount
){
    for(int i=0;i<edgeCount;i++){
        for(int j=0;j<5;j++){
            d_closest_id[i*5+j]=d_closest_id_cluster[i*5+j];
            d_closest_dis[i*5+j]=d_closest_dis_cluster[i*5+j];
        }
    }
}

/** Set up cluster edge/leaf masks and leafMap from closest_id of edge eid. */
__global__ void initializeClusterDC (
    int eid,
    int * e,
    int * belong,
    int * head,
    int * nxt,
    int * closest_id,
    int * edgeMask,
    int * leafMask,
    int * edgeMaskIndex,
    int * d_leafMap
){
    int x=belong[eid],y=e[eid];
    int otheid=head[y];
    while(e[otheid]!=x) assert(e[otheid]!=-1),otheid=nxt[otheid];
    
    int leafCount=0;
    for(int i=0;i<5;i++) leafMask[leafCount++]=closest_id[eid*5+i];
    for(int i=0;i<5;i++) d_leafMap[i]=leafMask[i];
    for(int i=0;i<5;i++) leafMask[leafCount++]=closest_id[otheid*5+i];
    for(int i=0;i<5;i++) d_leafMap[i+5]=leafMask[i+5];

    int edgeCount=0;
    edgeMask[edgeCount++]=eid, edgeMask[edgeCount++]=otheid;
    edgeMaskIndex[eid]=eid, edgeMaskIndex[otheid]=otheid;

    // printf("closest ids:\n");
    // for(int i=0;i<5;i++) printf("%d\t", closest_id[624*5+i]);
    // for(int i=0;i<5;i++) printf("%d\t", closest_id[541*5+i]);
    // printf("\n");
    
}



void MashPlacement::KPlacementDeviceArraysDC::deallocateDeviceArraysDC(){
    cudaFree(d_head);
    cudaFree(d_e);
    cudaFree(d_nxt);
    cudaFree(d_belong);
    cudaFree(d_closest_id);
    cudaFree(d_closest_id_cluster);
    cudaFree(d_dist);
    cudaFree(d_len);
    cudaFree(d_closest_dis);
    cudaFree(d_closest_dis_cluster);

    d_head = nullptr;
    d_e = nullptr;
    d_nxt = nullptr;
    d_belong = nullptr;
    d_closest_id = nullptr;
    d_closest_id_cluster = nullptr;
    d_dist = nullptr;
    d_len = nullptr;
    d_closest_dis = nullptr;
    d_closest_dis_cluster = nullptr;
}

void MashPlacement::KPlacementDeviceArraysDC::initializeDeviceArraysFromBackboneTreeDC(Tree *t)
{
    if (t == nullptr || t->root == nullptr) {
        fprintf(stderr, "Error: Backbone tree is null in initializeDeviceArraysFromBackboneTreeDC.\n");
        exit(1);
    }

    const size_t node_slots = static_cast<size_t>(totalNumSequences) * 2;
    const size_t edge_slots = static_cast<size_t>(totalNumSequences) * 8;

    std::vector<int> h_head(node_slots, -1);
    std::vector<int> h_e(edge_slots, -1);
    std::vector<int> h_nxt(edge_slots, -1);
    std::vector<int> h_belong(edge_slots, -1);
    std::vector<double> h_len(edge_slots, 2.0);

    size_t edgeCount = 0;
    std::function<void(Node*)> dfs = [&](Node* node) {
        if (node == nullptr) {
            return;
        }
        for (Node* child : node->children) {
            dfs(child);
        }
        if (node->parent == nullptr) {
            return;
        }
        if (edgeCount + 2 > edge_slots) {
            fprintf(stderr, "Error: Backbone tree exceeds DC edge capacity.\n");
            exit(1);
        }
        const int child_idx = node->idx;
        const int parent_idx = node->parent->idx;
        if (child_idx < 0 || parent_idx < 0 ||
            child_idx >= static_cast<int>(node_slots) ||
            parent_idx >= static_cast<int>(node_slots)) {
            fprintf(stderr, "Error: Backbone tree node index exceeds DC node capacity.\n");
            exit(1);
        }

        h_e[edgeCount] = parent_idx;
        h_len[edgeCount] = node->bl;
        h_belong[edgeCount] = child_idx;
        h_nxt[edgeCount] = h_head[static_cast<size_t>(child_idx)];
        h_head[static_cast<size_t>(child_idx)] = static_cast<int>(edgeCount);
        ++edgeCount;

        h_e[edgeCount] = child_idx;
        h_len[edgeCount] = node->bl;
        h_belong[edgeCount] = parent_idx;
        h_nxt[edgeCount] = h_head[static_cast<size_t>(parent_idx)];
        h_head[static_cast<size_t>(parent_idx)] = static_cast<int>(edgeCount);
        ++edgeCount;
    };
    dfs(t->root);

    const size_t expected_edges =
        static_cast<size_t>(std::max(0, numSequences * 4 - 4));
    if (edgeCount != expected_edges) {
        fprintf(stderr,
                "Error: Fixed-backbone DC expected %zu directed edges, found %zu.\n",
                expected_edges,
                edgeCount);
        exit(1);
    }

    auto checkCopy = [](cudaError_t err, const char* context) {
        if (err != cudaSuccess) {
            fprintf(stderr, "Gpu_ERROR: %s failed: %s\n", context, cudaGetErrorString(err));
            exit(1);
        }
    };

    checkCopy(cudaMemcpy(d_head, h_head.data(), node_slots * sizeof(int), cudaMemcpyHostToDevice),
              "cudaMemcpy d_head");
    checkCopy(cudaMemcpy(d_e, h_e.data(), edge_slots * sizeof(int), cudaMemcpyHostToDevice),
              "cudaMemcpy d_e");
    checkCopy(cudaMemcpy(d_nxt, h_nxt.data(), edge_slots * sizeof(int), cudaMemcpyHostToDevice),
              "cudaMemcpy d_nxt");
    checkCopy(cudaMemcpy(d_belong, h_belong.data(), edge_slots * sizeof(int), cudaMemcpyHostToDevice),
              "cudaMemcpy d_belong");
    checkCopy(cudaMemcpy(d_len, h_len.data(), edge_slots * sizeof(double), cudaMemcpyHostToDevice),
              "cudaMemcpy d_len");

    int *d_id = nullptr;
    int *d_from = nullptr;
    double *d_dis = nullptr;
    checkCopy(cudaMalloc(&d_id, node_slots * sizeof(int)), "cudaMalloc d_id");
    checkCopy(cudaMalloc(&d_from, node_slots * sizeof(int)), "cudaMalloc d_from");
    checkCopy(cudaMalloc(&d_dis, node_slots * sizeof(double)), "cudaMalloc d_dis");

    initializeDC<<<1024, 1024>>>(
        static_cast<int>(edge_slots),
        static_cast<int>(node_slots),
        d_closest_dis,
        d_closest_id,
        d_head,
        d_nxt,
        d_belong,
        d_e);
    checkCopy(cudaGetLastError(), "initializeDC");
    checkCopy(cudaDeviceSynchronize(), "initializeDC synchronize");

    // Restore the loaded topology after initializeDC reset the adjacency arrays.
    checkCopy(cudaMemcpy(d_head, h_head.data(), node_slots * sizeof(int), cudaMemcpyHostToDevice),
              "cudaMemcpy restore d_head");
    checkCopy(cudaMemcpy(d_e, h_e.data(), edge_slots * sizeof(int), cudaMemcpyHostToDevice),
              "cudaMemcpy restore d_e");
    checkCopy(cudaMemcpy(d_nxt, h_nxt.data(), edge_slots * sizeof(int), cudaMemcpyHostToDevice),
              "cudaMemcpy restore d_nxt");
    checkCopy(cudaMemcpy(d_belong, h_belong.data(), edge_slots * sizeof(int), cudaMemcpyHostToDevice),
              "cudaMemcpy restore d_belong");
    checkCopy(cudaMemcpy(d_len, h_len.data(), edge_slots * sizeof(double), cudaMemcpyHostToDevice),
              "cudaMemcpy restore d_len");

    for (int tip_id = 0; tip_id < numSequences; ++tip_id) {
        updateClosestNodesDC<<<1, 1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            tip_id,
            d_id,
            d_from,
            d_dis);
        checkCopy(cudaGetLastError(), "updateClosestNodesDC");
    }
    checkCopy(cudaDeviceSynchronize(), "updateClosestNodesDC synchronize");
    cudaFree(d_id);
    cudaFree(d_from);
    cudaFree(d_dis);

    idx = static_cast<int>(edgeCount);
    bd = numSequences;
}


/** Copy tree arrays to host and emit Newick via recursive DFS. */
void MashPlacement::KPlacementDeviceArraysDC::printTreeDC(
    std::vector <std::string> name,
    std::ofstream& output_,
    int rootNode){
    int * h_head = new int[totalNumSequences*2];
    int * h_e = new int[totalNumSequences*8];
    int * h_nxt = new int[totalNumSequences*8];
    double * h_len = new double[totalNumSequences*8];
    double * h_closest_dis = new double[totalNumSequences*40];
    int * h_closest_id = new int[totalNumSequences*40];
    auto err = cudaMemcpy(h_head, d_head, totalNumSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_e, d_e, totalNumSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_nxt, d_nxt, totalNumSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_len, d_len, totalNumSequences*8*sizeof(double),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }

    
    std::function<void(int,int)> print = [&](int node, int from) {
        // Build pos first, before deciding if this is a leaf or internal node
        std::vector<std::pair<int,int>> pos; // {edge_index, parent_node}
    
        if (h_head[node] != -1) {
            // Use a queue to recursively collapse chains of zero-length edges
            std::queue<std::pair<int,int>> toExpand; // {collapsed_node, came_from}
            toExpand.push({node, from});
    
            while (!toExpand.empty()) {
                auto [cur, curFrom] = toExpand.front();
                toExpand.pop();
    
                for (int i = h_head[cur]; i != -1; i = h_nxt[i]) {
                    if (h_e[i] == curFrom) continue; // skip parent direction
    
                    if (h_len[i] == 0) {
                        int collapsed = h_e[i];
                        // Check if collapsed is a leaf (only edge is back to cur)
                        bool isLeaf = (h_head[collapsed] == -1 ||
                                       h_nxt[h_head[collapsed]] == -1);
                        if (isLeaf) {
                            // Zero-length edge to a leaf: still include it
                            pos.push_back({i, cur});
                        } else {
                            // Internal collapsed node: expand it further
                            toExpand.push({collapsed, cur});
                        }
                    } else {
                        pos.push_back({i, cur});
                    }
                }
            }
        }
    
        if (pos.empty()) {
            // Leaf node (or degenerate case): just print name
            output_ << name[node];
        } else {
            output_ << "(";
            for (size_t i = 0; i < pos.size(); i++) {
                auto [edgeIdx, parent] = pos[i];
                print(h_e[edgeIdx], parent);
                output_ << ":" << h_len[edgeIdx];
                output_ << (i + 1 == pos.size() ? ')' : ',');
            }
        }
    };

    std::function<void(int,int)>  print_binary=[&](int node, int from){
        if(h_nxt[h_head[node]]!=-1){
            // printf("(");
            output_ << "(";
            std::vector <int> pos;
            for(int i=h_head[node];i!=-1;i=h_nxt[i])
                if(h_e[i]!=from)
                    pos.push_back(i);
            for(size_t i=0;i<pos.size();i++){
                print_binary(h_e[pos[i]],node);
                // printf(":");
                // printf("%.5g%c",h_len[pos[i]],i+1==pos.size()?')':',');
                output_ << ":";
                // output_ << "%.5g%c",h_len[pos[i]],i+1==pos.size()?')':',';
                output_ << h_len[pos[i]] << (i+1==pos.size()?')':',');
            }
        }
        // else std::cout<<name[node];
        else output_<<name[node];
    };

    const int print_root =
        (rootNode >= 0) ? rootNode : (totalNumSequences + bd - 2);
    if (!g_printBinaryNewick) {
        print(print_root, -1);
    } else {
        print_binary(print_root, -1);
    }
    output_ << ";";
}

__global__
void printSeqsDC(
    uint64_t * d_compressedSeqs,
    int num,
    int size
){
    for (int i=490;i<490+num;i++){
        for (int j=0;j<size;j++){
            printf("%ld\n",(d_compressedSeqs[i*size+j]));
        }
    }
}

/** Build backbone tree on first numSequences taxa: init, distances, initial tree,
 * updateClosestNodes, then sequential placement with distRangeConstruction and
 * calculateBranchLengthDC / updateTreeStructureDC / updateClosestNodesDC. */
void MashPlacement::KPlacementDeviceArraysDC::findBackboneTreeDC(
    Param& params,
    const MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArraysDC& msaDeviceArrays,
    const KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost,
    int gpuNum
){ 
    cudaError_t err;
    if(params.in == "d"){
        matrixReader.distConstructionOnGpu(params, 0, d_dist);
    }
    int * d_id;
    err = cudaMalloc(&d_id, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed 1!\n");
        exit(1);
    }
    int * d_from;
    err = cudaMalloc(&d_from, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed 2 !\n");
        exit(1);
    }
    double * d_dis;
    err = cudaMalloc(&d_dis, totalNumSequences*2*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed 3!\n");
        exit(1);
    }
    
    

    /*
    Initialize closest nodes by inifinite
    */
    // int threadNum = 1024, blockNum = (totalNumSequences*4-4+threadNum-1)/threadNum;
    int threadNum = 1024, blockNum = 1024;
    initializeDC <<<blockNum, threadNum>>> (
        totalNumSequences*8,
        totalNumSequences*2,
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
    if(params.in == "r"){
        mashDeviceArrays.distConstructionOnGpuForBackboneDC(
            params,
            1,
            d_dist
        );
    }
    else if(params.in == "d"){
        matrixReader.distConstructionOnGpu(
            params,
            1,
            d_dist
        );
    }
    else if(params.in == "m"){
        msaDeviceArrays.distConstructionOnGpuForBackboneDC(
            params,
            1,
            d_dist
        );
    }
    
    buildInitialTreeDC <<<1,1>>> (
        totalNumSequences,
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
        updateClosestNodesDC <<<1,1>>> (
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
 
    std::cerr << "Finished initial tree construction" << std::endl;
    auto backboneStart = std::chrono::high_resolution_clock::now();
    thrust::device_vector <thrust::tuple<int,double,double>> minPos(numSequences*4-4);


    std::chrono::nanoseconds disTime(0), treeTime(0);
    for(int i=bd;i<numSequences;i++){
        auto disStart = std::chrono::high_resolution_clock::now();
        if(params.in == "r"){
            mashDeviceArrays.distRangeConstructionOnGpuDC(
                params,
                i,
                d_dist,
                0,
                i-1
            );
        }
        else if(params.in == "d"){
            matrixReader.distConstructionOnGpu(
                params,
                i,
                d_dist
            );
        }
        else if(params.in == "m"){
            msaDeviceArrays.distRangeConstructionOnGpuDC(
                params,
                i,
                d_dist,
                0,
                i-1
            );
        }

        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        calculateBranchLengthDC <<<blockNum,threadNum>>> (
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
            false
        );

        auto iter=thrust::min_element(minPos.begin(),minPos.begin()+numSequences*4-4,compare_tupleDC());
        thrust::tuple<int,double,double> smallest=*iter;
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid=thrust::get<0>(smallest);
        double fracLen=thrust::get<1>(smallest),addLen=thrust::get<2>(smallest);
        updateTreeStructureDC <<<1,1>>>(
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
            totalNumSequences
        );
        idx+=4;

        /*
        Update closest nodes
        */
        updateClosestNodesDC <<<1,1>>> (
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
        auto treeEnd = std::chrono::high_resolution_clock::now();
        disTime += disEnd - disStart;
        treeTime += treeEnd - treeStart;
    }

    cudaDeviceSynchronize();
    auto backboneEnd = std::chrono::high_resolution_clock::now();
    auto backboneTime = backboneEnd - backboneStart;
    std::cerr << "Finished backbone construction in: "<< backboneTime.count()/1000000 << " ms\n";
    
    return;
}

/** Assign each non-backbone sequence to k-closest edge (clusterID); batch distance
 * construction and calculateBranchLengthDC per sequence. */
void MashPlacement::KPlacementDeviceArraysDC::findClustersDC(
    Param& params,
    const MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArraysDC& msaDeviceArrays,
    KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost
){ 
    cudaError_t err;
    int idx=params.backboneSize*4-4;
    
    kplacementDeviceArraysHost.clusterID = new int[totalNumSequences];
    kplacementDeviceArraysHost.clusterTopEdgeIds.assign(
        totalNumSequences, std::vector<int>{});
    kplacementDeviceArraysHost.clusterTopEdgeAdditionalDistances.assign(
        totalNumSequences, std::vector<double>{});
    thrust::device_vector <thrust::tuple<int,double,double>> minPos(totalNumSequences*4-4);
    thrust::host_vector <thrust::tuple<int,double,double>> hostMinPos(
        totalNumSequences*4-4);
    int threadNum = 1024, blockNum = 1024;

    auto clusterStart = std::chrono::high_resolution_clock::now();
    const bool audit_distances =
        std::getenv("DIPPER_DC_DISTANCE_AUDIT") != nullptr;
    const bool nearest_anchor_routing =
        std::getenv("DIPPER_DC_NEAREST_ANCHOR_ROUTING") != nullptr;
    std::vector<int> routing_head;
    std::vector<int> routing_nxt;
    std::vector<int> routing_e;
    std::vector<int> routing_belong;
    if (nearest_anchor_routing) {
        routing_head.resize(totalNumSequences * 2);
        routing_nxt.resize(totalNumSequences * 8);
        routing_e.resize(totalNumSequences * 8);
        routing_belong.resize(totalNumSequences * 8);
        cudaMemcpy(routing_head.data(), d_head,
                   routing_head.size() * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(routing_nxt.data(), d_nxt,
                   routing_nxt.size() * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(routing_e.data(), d_e,
                   routing_e.size() * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(routing_belong.data(), d_belong,
                   routing_belong.size() * sizeof(int), cudaMemcpyDeviceToHost);
    }
    std::vector<size_t> nearest_backbone_counts(
        audit_distances ? numSequences : 0, 0);

    uint64_t localBatchSize = params.batchSize;
    for(int i=numSequences;i<totalNumSequences;i+=localBatchSize){
        if (i+localBatchSize>totalNumSequences) localBatchSize=totalNumSequences-i;
        // std::cerr << "Processing batch: "<< i << " to " << i+localBatchSize << "\n";
        
        /* copy data to d_hashListConst or d_compressedSeqsConst */
        if (params.in == "r") 
            err = cudaMemcpy(mashDeviceArrays.d_hashListConst, mashDeviceArrays.h_hashList+i*params.sketchSize, localBatchSize*params.sketchSize*sizeof(uint64_t),cudaMemcpyHostToDevice);
        else if (params.in == "m"){
            size_t maxLengthCompressed = (msaDeviceArrays.d_seqLen + 15) / 16;
            err = cudaMemcpy(msaDeviceArrays.d_compressedSeqsConst, msaDeviceArrays.h_compressedSeqs+i*maxLengthCompressed, 1ll*localBatchSize*maxLengthCompressed*sizeof(uint64_t),cudaMemcpyHostToDevice);
        }
        else {
            std::cerr << "Error: Input type must be unaligned or aligned for clustering based approach\n";
            exit(1);
        }
        if (err != cudaSuccess) {
            fprintf(stderr, "Gpu_ERROR: d_hashListConst cudaMemcpy failed!\n");
            exit(1);
        }

        for (int j=i;j<i+localBatchSize;j++){
            
            auto disStart = std::chrono::high_resolution_clock::now();
            if(params.in == "r"){
                mashDeviceArrays.distRangeConstructionOnGpuDC(
                    params,
                    j-i,
                    d_dist,
                    0,
                    numSequences-1,
                    true
                );
            }
            else if(params.in == "d"){
                matrixReader.distConstructionOnGpu(
                    params,
                    j,
                    d_dist
                );
            }
            else if(params.in == "m"){
                msaDeviceArrays.distRangeConstructionOnGpuDC(
                    params,
                    j-i,
                    d_dist,
                    0,
                    numSequences-1,
                    true
                );
            }

            std::vector<double> host_dist;
            if (audit_distances || nearest_anchor_routing) {
                host_dist.resize(numSequences);
                cudaMemcpy(
                    host_dist.data(), d_dist,
                    numSequences * sizeof(double), cudaMemcpyDeviceToHost);
            }
            if (audit_distances) {
                size_t finite_count = 0;
                size_t nan_count = 0;
                size_t inf_count = 0;
                double min_distance = std::numeric_limits<double>::infinity();
                int min_index = -1;
                for (int ref = 0; ref < numSequences; ++ref) {
                    const double value = host_dist[ref];
                    if (std::isnan(value)) {
                        ++nan_count;
                    } else if (!std::isfinite(value)) {
                        ++inf_count;
                    } else {
                        ++finite_count;
                        if (value < min_distance) {
                            min_distance = value;
                            min_index = ref;
                        }
                    }
                }
                if (min_index >= 0) {
                    ++nearest_backbone_counts[static_cast<size_t>(min_index)];
                }
                if (j == i || j == i + 1 || j == i + localBatchSize - 1) {
                    std::cerr
                        << "[DIPPER] D&C distance audit: query_index=" << j
                        << " finite=" << finite_count
                        << " nan=" << nan_count
                        << " inf=" << inf_count
                        << " min=" << min_distance
                        << " min_backbone_index=" << min_index << "\n";
                }
            }

            // double * h_dis = new double[numSequences];
            // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
            // fprintf(stderr, "%d\n",i);
            // for(int j=0;j<i;j++) std::cerr<<h_dis[j]<<" ";std::cerr<<'\n';


            auto disEnd = std::chrono::high_resolution_clock::now();
            auto treeStart = std::chrono::high_resolution_clock::now();
            
            calculateBranchLengthDC <<<blockNum,threadNum>>> (
                j,
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
                std::getenv("DIPPER_DC_PENALIZE_NEGATIVE_PENDANT") != nullptr
            );
            auto iter=thrust::min_element(minPos.begin(),minPos.begin()+numSequences*4-4,compare_tupleDC());
            thrust::tuple<int,double,double> smallest=*iter;
            kplacementDeviceArraysHost.clusterID[j] = thrust::get<0>(smallest);
            int requested_topk = std::max(
                1, kplacementDeviceArraysHost.requestedClusterTopK);
            if (const char* value = std::getenv("DIPPER_DC_PROVENANCE_TOPK")) {
                requested_topk = std::max(requested_topk, std::max(1, std::atoi(value)));
            }
            hostMinPos = minPos;
            std::vector<std::pair<double, int>> ranked_edges;
            ranked_edges.reserve(numSequences * 2 - 3);
            for (int edge_idx = 0;
                 edge_idx < numSequences * 4 - 4;
                 ++edge_idx) {
                const int candidate_id = thrust::get<0>(hostMinPos[edge_idx]);
                const double additional_distance =
                    thrust::get<2>(hostMinPos[edge_idx]);
                if (candidate_id != edge_idx ||
                    !std::isfinite(additional_distance) ||
                    additional_distance >= 2.0) {
                    continue;
                }
                ranked_edges.emplace_back(additional_distance, candidate_id);
            }
            const int keep = std::min(
                requested_topk,
                static_cast<int>(ranked_edges.size()));
            std::partial_sort(
                ranked_edges.begin(),
                ranked_edges.begin() + keep,
                ranked_edges.end(),
                [](const auto& lhs, const auto& rhs) {
                    if (lhs.first != rhs.first) {
                        return lhs.first < rhs.first;
                    }
                    return lhs.second < rhs.second;
                });
            auto& top_edges =
                kplacementDeviceArraysHost.clusterTopEdgeIds[j];
            auto& top_distances =
                kplacementDeviceArraysHost.clusterTopEdgeAdditionalDistances[j];
            top_edges.reserve(requested_topk);
            top_distances.reserve(requested_topk);
            const int original_top1 =
                kplacementDeviceArraysHost.clusterID[j];
            top_edges.push_back(original_top1);
            top_distances.push_back(thrust::get<2>(smallest));
            for (int rank = 0; rank < keep; ++rank) {
                if (ranked_edges[rank].second == original_top1) {
                    continue;
                }
                top_edges.push_back(ranked_edges[rank].second);
                top_distances.push_back(ranked_edges[rank].first);
                if (static_cast<int>(top_edges.size()) == requested_topk) {
                    break;
                }
            }
            if (nearest_anchor_routing) {
                std::vector<int> anchors(numSequences);
                std::iota(anchors.begin(), anchors.end(), 0);
                std::partial_sort(
                    anchors.begin(),
                    anchors.begin() + std::min(
                        static_cast<int>(anchors.size()), requested_topk),
                    anchors.end(),
                    [&](int lhs, int rhs) {
                        if (host_dist[lhs] != host_dist[rhs]) {
                            return host_dist[lhs] < host_dist[rhs];
                        }
                        return lhs < rhs;
                    });
                top_edges.clear();
                top_distances.clear();
                std::unordered_set<int> seen_edges;
                for (const int anchor : anchors) {
                    if (!std::isfinite(host_dist[anchor])) continue;
                    std::vector<int> incident_edges;
                    for (int edge_idx = routing_head[anchor];
                         edge_idx != -1;
                         edge_idx = routing_nxt[edge_idx]) {
                        int canonical_edge = edge_idx;
                        if (routing_belong[canonical_edge] <
                            routing_e[canonical_edge]) {
                            const int neighbor =
                                routing_e[canonical_edge];
                            for (int reverse =
                                     routing_head[neighbor];
                                 reverse != -1;
                                 reverse = routing_nxt[reverse]) {
                                if (routing_e[reverse] == anchor) {
                                    canonical_edge = reverse;
                                    break;
                                }
                            }
                        }
                        if (seen_edges.insert(canonical_edge).second) {
                            incident_edges.push_back(canonical_edge);
                        }
                    }
                    std::sort(incident_edges.begin(), incident_edges.end());
                    for (const int edge_idx : incident_edges) {
                        top_edges.push_back(edge_idx);
                        top_distances.push_back(host_dist[anchor]);
                        if (static_cast<int>(top_edges.size()) == requested_topk) {
                            break;
                        }
                    }
                    if (static_cast<int>(top_edges.size()) == requested_topk) {
                        break;
                    }
                }
                if (top_edges.empty()) {
                    throw std::runtime_error(
                        "nearest-anchor routing produced no candidate edge");
                }
                kplacementDeviceArraysHost.clusterID[j] = top_edges.front();
            }
            // std::cerr << "Sequence " << j << " assigned to cluster " << kplacementDeviceArraysHost.clusterID[j] << "\n";
        }
    }
    auto clusterEnd = std::chrono::high_resolution_clock::now();
    auto clusterTime = clusterEnd - clusterStart;


    std::cerr << "Finished clustering in: "<< clusterTime.count()/1000000 << " ms\n";
    if (const char* provenance_path = std::getenv("DIPPER_DC_PROVENANCE")) {
        std::ofstream provenance(provenance_path);
        if (!provenance) {
            throw std::runtime_error(
                std::string("failed to open D&C provenance output: ") +
                provenance_path);
        }
        provenance
            << "query_index\tquery_name\tsector_edge\ttop_edge_ids"
               "\ttop_edge_scores\ttop1_top2_margin\n";
        for (int query = numSequences; query < totalNumSequences; ++query) {
            const auto& edges =
                kplacementDeviceArraysHost.clusterTopEdgeIds[query];
            const auto& scores =
                kplacementDeviceArraysHost.clusterTopEdgeAdditionalDistances[query];
            provenance << query << '\t';
            if (query < static_cast<int>(matrixReader.name.size())) {
                provenance << matrixReader.name[query];
            }
            provenance << '\t' << kplacementDeviceArraysHost.clusterID[query]
                       << '\t';
            for (size_t rank = 0; rank < edges.size(); ++rank) {
                if (rank) provenance << ',';
                provenance << edges[rank];
            }
            provenance << '\t';
            for (size_t rank = 0; rank < scores.size(); ++rank) {
                if (rank) provenance << ',';
                provenance << scores[rank];
            }
            provenance << '\t';
            if (scores.size() >= 2) {
                provenance << scores[1] - scores[0];
            }
            provenance << '\n';
        }
        std::cerr << "[DIPPER] D&C provenance written to "
                  << provenance_path << "\n";
    }
    if (audit_distances) {
        size_t distinct_nearest = 0;
        size_t largest_group = 0;
        for (const size_t count : nearest_backbone_counts) {
            if (count > 0) {
                ++distinct_nearest;
                largest_group = std::max(largest_group, count);
            }
        }
        std::cerr
            << "[DIPPER] D&C distance audit summary: distinct_nearest_backbone="
            << distinct_nearest << " largest_group=" << largest_group << "\n";
    }

    /* Copy data from device to host */
    // err = cudaMemcpy(kplacementDeviceArraysHost.h_dist, d_dist, totalNumSequences*sizeof(double),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_dist cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_head, d_head, totalNumSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_head cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_e, d_e, totalNumSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_e cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_len, d_len, totalNumSequences*8*sizeof(double),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_len cudaMemcpy failed!\n");
    //     exit(1);
    // }
    // err = cudaMemcpy(kplacementDeviceArraysHost.h_nxt, d_nxt, totalNumSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_nxt cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_belong, d_belong, totalNumSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_belong cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_closest_dis, d_closest_dis, totalNumSequences*40*sizeof(double),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_closest_dis cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_closest_id, d_closest_id, totalNumSequences*40*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_closest_id cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_closest_dis_cluster, d_closest_dis_cluster, totalNumSequences*40*sizeof(double),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_closest_dis_cluster cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_closest_id_cluster, d_closest_id_cluster, totalNumSequences*40*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_closest_id_cluster cudaMemcpy failed!\n");
    //     exit(1);
    // }

    cudaDeviceSynchronize();

    std::cerr << "Finished data transfer\n";
    return;
}

/** Batch variant of findClustersDC for a single clustering batch index. */
void MashPlacement::KPlacementDeviceArraysDC::findClustersDC_batch(
    Param& params,
    const MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArraysDC& msaDeviceArrays,
    KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost,
    const int clusteringBatchIdx
){ 
    cudaError_t err;
    thrust::device_vector <thrust::tuple<int,double,double>> minPos(totalNumSequences*4-4);
    int threadNum = 1024, blockNum = 1024;

    auto clusterStart = std::chrono::high_resolution_clock::now();

    for(int i=0;i<params.batchSize && i<params.totalNumSeqs-clusteringBatchIdx*params.batchSize;i++){
        
        if(params.in == "r"){
            mashDeviceArrays.distRangeConstructionOnGpuDC(
                params,
                i,
                d_dist,
                0,
                numSequences-1,
                true
            );
        }
        else if(params.in == "d"){
            matrixReader.distConstructionOnGpu(
                params,
                i,
                d_dist
            );
        }
        else if(params.in == "m"){
            msaDeviceArrays.distRangeConstructionOnGpuDC(
                params,
                i,
                d_dist,
                0,
                numSequences-1,
                true
            );
        }


        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        
        calculateBranchLengthDC <<<blockNum,threadNum>>> (
            clusteringBatchIdx*params.batchSize + i,
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
            std::getenv("DIPPER_DC_PENALIZE_NEGATIVE_PENDANT") != nullptr
        );
        auto iter=thrust::min_element(minPos.begin(),minPos.begin()+numSequences*4-4,compare_tupleDC());
        thrust::tuple<int,double,double> smallest=*iter;
        kplacementDeviceArraysHost.clusterID[clusteringBatchIdx*params.batchSize + i] = thrust::get<0>(smallest);

    }
    cudaDeviceSynchronize();

    return;
}

__global__ void rearrangeHashListInClusterDC(
    int numSequences,
    int sketchSize,
    uint64_t * original,
    uint64_t * target
){
    int tx = threadIdx.x, bx = blockIdx.x;
    int bs = blockDim.x;
    int idx = tx+bs*bx;
    // if(idx>=numSequences) return;
    for (int idx_=idx; idx_<numSequences; idx_+=bs*gridDim.x){
        if (idx_ >= numSequences) return;
        for(int i=0;i<sketchSize;i++){
            target[i*numSequences+idx_] = original[idx_*sketchSize + i];
        }
    }
    // for(int i=0;i<sketchSize;i++){
    //     target[i*numSequences+idx] = original[idx*sketchSize + i];
    // }
}

void transferMashClusterInfoDC(
    MashPlacement::MashDeviceArraysDC& mashDeviceArrays,
    std::vector<int> leafList,
    MashPlacement::Param& params
){
    uint64_t * hashListLocal = new uint64_t[params.backboneSize*params.sketchSize];
    uint64_t * hashList = mashDeviceArrays.h_hashList;
    int l=0;
    for (auto &leaf: leafList){
        memcpy(hashListLocal+l*params.sketchSize, hashList+leaf*params.sketchSize, params.sketchSize*sizeof(uint64_t));
        l++;
    }

    auto err = cudaMemcpy(mashDeviceArrays.d_hashListConst, hashListLocal, params.backboneSize*params.sketchSize*sizeof(uint64_t),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: d_hashListConst cudaMemcpy failed!\n");
        exit(1);
    }

    /* Rearrange only for backbone tree */
    uint64_t * temp_hashList;
    err = cudaMalloc(&temp_hashList, params.sketchSize*params.backboneSize*sizeof(uint64_t));
    if (err != cudaSuccess){
        fprintf(stderr, "Gpu_ERROR: temp_hashList cudaMalloc failed!\n");
        exit(1);
    }
    int threadsPerBlock = 1024;
    int blocksPerGrid = 1024;
    rearrangeHashListInClusterDC <<<blocksPerGrid, threadsPerBlock >>>(
        params.backboneSize,
        int(params.sketchSize),
        mashDeviceArrays.d_hashListConst,
        temp_hashList
    );
    std::swap(mashDeviceArrays.d_hashListConst, temp_hashList);
    cudaFree(temp_hashList);
    
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
    delete[] hashListLocal;
    return;
}

void transferMashClusterInfoDCRecursiveBackbone(
    MashPlacement::MashDeviceArraysDC& mashDeviceArrays,
    std::vector<int> leafList,
    MashPlacement::Param& params
){
    uint64_t * hashListLocal = new uint64_t[params.backboneSize*params.sketchSize];
    uint64_t * hashList = mashDeviceArrays.h_hashList;
    int l=0;
    for (auto &leaf: leafList){
        memcpy(hashListLocal+l*params.sketchSize, hashList+leaf*params.sketchSize, params.sketchSize*sizeof(uint64_t));
        l++;
    }

    auto gpuMemLoc=mashDeviceArrays.d_hashListBackbone+params.backboneSize*params.sketchSize;
    auto err = cudaMemcpy(gpuMemLoc, hashListLocal, params.backboneSize*params.sketchSize*sizeof(uint64_t),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: d_hashListBackbone cudaMemcpy failed!\n");
        exit(1);
    }

    uint64_t * temp_hashList;
    err = cudaMalloc(&temp_hashList, params.sketchSize*params.backboneSize*sizeof(uint64_t));
    if (err != cudaSuccess){
        fprintf(stderr, "Gpu_ERROR: temp_hashList cudaMalloc failed!\n");
        exit(1);
    }
    int threadsPerBlock = 1024;
    int blocksPerGrid = 1024;
    rearrangeHashListInClusterDC <<<blocksPerGrid, threadsPerBlock >>>(
        params.backboneSize,
        int(params.sketchSize),
        gpuMemLoc,
        temp_hashList
    );
    std::swap(gpuMemLoc, temp_hashList);
    cudaFree(temp_hashList);
    
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
    delete[] hashListLocal;
    return;
}


void transferMsaClusterInfoDC(
    MashPlacement::MSADeviceArraysDC& mashDeviceArrays,
    std::vector<int> leafList,
    MashPlacement::Param& params
){
    size_t maxLengthCompressed = (mashDeviceArrays.d_seqLen + 15) / 16;
    uint64_t * compressedSeqs_local = new uint64_t[params.backboneSize*maxLengthCompressed];
    uint64_t * compressedSeqs = mashDeviceArrays.h_compressedSeqs;

    int l=0;
    for (auto &leaf: leafList){
        memcpy(compressedSeqs_local+1ll*l*maxLengthCompressed, compressedSeqs+1ll*leaf*maxLengthCompressed, 1ll*maxLengthCompressed*sizeof(uint64_t));
        l++;
    }
    auto err = cudaMemcpy(mashDeviceArrays.d_compressedSeqsConst, compressedSeqs_local, 1ll*params.backboneSize*maxLengthCompressed*sizeof(uint64_t),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: d_hashListConst cudaMemcpy failed!\n");
        exit(1);
    }
    cudaDeviceSynchronize();
    delete[] compressedSeqs_local;
    return;
}

void transferMsaClusterInfoDCRecursiveBackbone(
    MashPlacement::MSADeviceArraysDC& mashDeviceArrays,
    std::vector<int> leafList,
    MashPlacement::Param& params
){
    size_t maxLengthCompressed = (mashDeviceArrays.d_seqLen + 15) / 16;
    uint64_t * compressedSeqs_local = new uint64_t[params.backboneSize*maxLengthCompressed];
    uint64_t * compressedSeqs = mashDeviceArrays.h_compressedSeqs;

    int l=0;
    for (auto &leaf: leafList){
        memcpy(compressedSeqs_local+1ll*l*maxLengthCompressed, compressedSeqs+1ll*leaf*maxLengthCompressed, 1ll*maxLengthCompressed*sizeof(uint64_t));
        l++;
    }
    auto err = cudaMemcpy(mashDeviceArrays.d_compressedSeqsConst, compressedSeqs_local, 1ll*params.backboneSize*maxLengthCompressed*sizeof(uint64_t),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: d_hashListConst cudaMemcpy failed!\n");
        exit(1);
    }
    cudaDeviceSynchronize();
    delete[] compressedSeqs_local;
    return;
}




__global__ 
void print_d_hashListConstDC(
    uint64_t * d_hashListConst,
    int sketchSize,
    int batchSize
){
    for (int i=0;i<batchSize;i++){
        for (int j=0;j<sketchSize;j++){
            printf("%llu ", d_hashListConst[i*sketchSize+j]);
        }
        printf("\n");
    }
}

__global__ 
void resetEdgeMaskIndexDC(int * d_edgeMaskIndex, int size){
    int tx = threadIdx.x, bx = blockIdx.x;
    int bs = blockDim.x, gs = gridDim.x;
    int idx = tx+bs*bx;
    // if(idx>=size) return;
    for (int idx_=idx; idx_<size; idx_+=bs*gs){
        if (idx_ >= size) return;
        d_edgeMaskIndex[idx_] = -1;
    }
    // d_edgeMaskIndex[idx] = -1;
    
}

__global__ 
void resetIdFromDisDC(int * id, int * from, double * dis, int size){
    int tx = threadIdx.x, bx = blockIdx.x;
    int bs = blockDim.x;
    int idx = tx+bs*bx;
    if(idx>=size) return;
    id[idx] = -1;
    from[idx] = -1;
    dis[idx] = 2.0;   
}

/** For each cluster: transfer MASH/MSA data, initializeClusterDC, then place each leaf
 * via distSpecialIDConstruction, calculateBranchLengthSpecialIDDC, updateTreeStructureInClusterDC,
 * updateClusterInfoDC, updateClosestNodesInClusterDC. */
void MashPlacement::KPlacementDeviceArraysDC::findClusterTreeDC(
    Param& params,
    MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    MSADeviceArraysDC& msaDeviceArrays,
    KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost,
    std::vector<int>& largeClustersIdx
){ 
    int idx=params.backboneSize*4-4;
    int threadNum = 1024, blockNum = 1024;
    int * d_id;
    auto err = cudaMalloc(&d_id, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int * d_from;
    err = cudaMalloc(&d_from, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    double * d_dis;
    err = cudaMalloc(&d_dis, totalNumSequences*2*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    auto * cluster_id = kplacementDeviceArraysHost.clusterID;
    std::vector<std::vector <int>> contains(numSequences*4-4);

    for(int i=numSequences;i<totalNumSequences;i++) contains[cluster_id[i]].push_back(i);

    
    int * d_edgeMask;
    err = cudaMalloc(&d_edgeMask, totalNumSequences*4*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int * d_edgeMaskIndex;
    err = cudaMalloc(&d_edgeMaskIndex, totalNumSequences*4*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int * d_leafMask;
    err = cudaMalloc(&d_leafMask, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int * d_leafMap;
    err = cudaMalloc(&d_leafMap, params.backboneSize*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int insertLeafCount=numSequences;
    thrust::device_vector <thrust::tuple<int,double,double>> minPos(totalNumSequences*4-4);
    
    
    std::vector<int> leafList (params.backboneSize);
    std::unordered_map <int,int> h_leafMap;
    
    // for(int i=0;i<numSequences*4-4;i++){ 
    int actual_placement=0;
    int i=0;
    while (i<numSequences*4-4) {
        std::cerr << "Processing batch "<< i << " out of " << numSequences*4-4 <<  std::endl;
        h_leafMap.clear();
        int startClusterID = i;
        int batchSize=0;
        if (params.backboneSize<contains[i].size()) {
            std::cerr << "Cluster " << i << " size (" << contains[i].size() <<") is larger than backbone size\n";
            exit(1);
        }
        while (batchSize<params.backboneSize && i<numSequences*4-4) {
            if (contains[i].size() + batchSize >= params.backboneSize) break;
            for (auto &leaf: contains[i]) { 
                h_leafMap[leaf] = batchSize;
                leafList[batchSize++] = leaf; 
            }
            i++;
        }

        if(params.in == "r")
            transferMashClusterInfoDC(mashDeviceArrays, leafList, params);
        else if (params.in == "m")
            transferMsaClusterInfoDC(msaDeviceArrays, leafList, params);
        else {
            std::cerr << "Error: Input type must be unaligned or aligned for clustering based approach\n";
            exit(1);
        }

        int localCount_=0;
        for (int j=startClusterID;j<i;j++) {

            if (contains[j].size() == 0) continue;
            resetEdgeMaskIndexDC<<<blockNum,threadNum>>>(d_edgeMaskIndex, totalNumSequences*4);
            cudaDeviceSynchronize();
            initializeClusterDC <<<1,1>>>(
                j,
                d_e,
                d_belong,
                d_head,
                d_nxt,
                d_closest_id,
                d_edgeMask,
                d_leafMask,
                d_edgeMaskIndex,
                d_leafMap
            );
            int edgeCount=2, leafCount=10;
            for(auto &leaf:contains[j]) {
                actual_placement++;

                int leaf_idx_in_cluster = h_leafMap[leaf];
                if(params.in == "r"){
                    
                    mashDeviceArrays.distSpecialIDConstructionOnGpuDC(
                        params,
                        localCount_,
                        d_dist,
                        leafCount,
                        d_leafMask,
                        d_leafMap
                    );

                } else if(params.in == "m"){
                    msaDeviceArrays.distSpecialIDConstructionOnGpuDC(
                        params,
                        localCount_,
                        d_dist,
                        leafCount,
                        d_leafMask,
                        d_leafMap
                    );
                } 
                localCount_++;    

                calculateBranchLengthSpecialIDDC <<<blockNum,threadNum>>> (
                    j,
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
                    edgeCount,
                    d_edgeMask
                );
                auto iter=thrust::min_element(minPos.begin(),minPos.begin()+edgeCount,compare_tupleDC());
                thrust::tuple<int,double,double> smallest=*iter;

                /*
                Update Tree Structure
                */

                int eid=thrust::get<0>(smallest);
                double fracLen=thrust::get<1>(smallest),addLen=thrust::get<2>(smallest);
                // std::cerr << "eid: " << eid << " Cluster ID: " << j << " dist " << addLen << " dist2 " << fracLen << "\n";
                updateTreeStructureInClusterDC <<<1,1>>>(
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
                    leaf,
                    idx,
                    totalNumSequences,
                    insertLeafCount
                );
                idx+=4, insertLeafCount++;

                /*
                Update edgeMask and leafMask
                */

                updateClusterInfoDC<<<1,1>>> (
                    leaf,
                    idx,
                    d_leafMask,
                    d_edgeMask,
                    d_edgeMaskIndex,
                    edgeCount,
                    leafCount,
                    d_leafMap,
                    leaf_idx_in_cluster
                );

                /*
                Update closest nodes
                */

                updateClosestNodesInClusterDC <<<1,1>>> (
                    d_head,
                    d_nxt,
                    d_e,
                    d_len,
                    d_closest_dis,
                    d_closest_id,
                    leaf,
                    d_id,
                    d_from,
                    d_dis,
                    j,
                    d_belong,
                    d_edgeMaskIndex
                );

                edgeCount+=4, leafCount++;
                // std::cerr << "leaf: " << leaf << " Cluster ID: " << j << " eid " << eid << " dist " << addLen << " dist2 " << fracLen << "\n";
                // if (leaf == 879) exit(0);
                // exit(0);
            }
        }
        // exit(0);
        assert(localCount_==batchSize);
        
    }

    cudaDeviceSynchronize();

}

bool read_line(gzFile file, std::string& line, int maxLength) {
    char buffer[maxLength];
    if (gzgets(file, buffer, maxLength) == Z_NULL) {
        return false;
    }
    line = buffer;
    if (!line.empty() && line.back() == '\n') line.pop_back();
    return true;
}

bool read_binary_blob(gzFile file, size_t n, uint64_t* out_data) {
    size_t total_bytes = n * sizeof(uint64_t);
    char* raw_ptr = reinterpret_cast<char*>(out_data);

    size_t bytes_read = 0;
    while (bytes_read < total_bytes) {
        int chunk = gzread(file, raw_ptr + bytes_read, total_bytes - bytes_read);
        if (chunk <= 0) {
            std::cerr << "Error: Unexpected end of file while reading binary data." << std::endl;
            delete[] out_data;
            out_data = nullptr;
            return false;
        }
        bytes_read += chunk;
    }
    return true;
}

void MashPlacement::KPlacementDeviceArraysDC::findClusterTreeDC_batch(
    Param& params,
    MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    MSADeviceArraysDC& msaDeviceArrays,
    KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost,
    const std::string dir,
    std::vector<bool>& isCluster
){ 
    cudaError_t err;
    int idx=params.backboneSize*4-4;
    int threadNum = 1024, blockNum = 1024;
    int * d_id;
    err = cudaMalloc(&d_id, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int * d_from;
    err = cudaMalloc(&d_from, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    double * d_dis;
    err = cudaMalloc(&d_dis, totalNumSequences*2*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    auto * cluster_id = kplacementDeviceArraysHost.clusterID;
    std::vector<std::vector <int>> contains(numSequences*4-4);

    for(int i=numSequences;i<totalNumSequences;i++) contains[cluster_id[i]].push_back(i);
    
   
    int * d_edgeMask;
    err = cudaMalloc(&d_edgeMask, totalNumSequences*4*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int * d_edgeMaskIndex;
    err = cudaMalloc(&d_edgeMaskIndex, totalNumSequences*4*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int * d_leafMask;
    err = cudaMalloc(&d_leafMask, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int * d_leafMap;
    err = cudaMalloc(&d_leafMap, params.backboneSize*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int insertLeafCount=numSequences;
    thrust::device_vector <thrust::tuple<int,double,double>> minPos(totalNumSequences*4-4);
    
    
    std::vector<int> leafList (params.backboneSize);
    std::unordered_map <int,int> h_leafMap;
    
    int test_var = 0;
    int actual_placement = 0;
    int clustersPerBatchFile = 1000;
    int clusterFiles = (numSequences*4-4)/clustersPerBatchFile + ((numSequences*4-4)%clustersPerBatchFile != 0);
    
    for (int bc=0; bc<clusterFiles;bc++){
        /* Find number of sequences in the file*/
        int seqInFile=0;
        std::unordered_map<int, std::vector<int>> clusterToCompressSeqIdxMap;
        std::unordered_map<int, int> localIdxToOriginalIdxMap;
        for (int i=bc*clustersPerBatchFile;i<(bc+1)*clustersPerBatchFile && i<numSequences*4-4;i++){
            if (!isCluster[i]) continue;
            clusterToCompressSeqIdxMap[i] = std::vector<int>();
            seqInFile+=contains[i].size();
        }

        if (seqInFile == 0) continue;


        /* Read all sequences in the file*/
        std::string path = dir + "/" + std::to_string(bc) + ".gz";
        std::cerr << "Reading file: " << path << " with " << seqInFile << " sequences\n";
        size_t maxLengthCompressed = (msaDeviceArrays.d_seqLen + 15) / 16;
        uint64_t * compressedSeqs_local = new uint64_t[seqInFile*maxLengthCompressed];

        gzFile file = gzopen(path.c_str(), "rb");
        if (!file) {
            std::cerr << "gzopen failed: " << path << " " << strerror(errno) << "\n";
            exit(0);
        }
        std::string line;
        uint64_t counter=0;
        while (read_line(file, line, maxLengthCompressed)) {
            if (line.empty() || line[0] != '>') continue;

            std::string name, id;
            std::istringstream ss(line.substr(1));
            std::getline(ss, name, '\t');
            std::getline(ss, id, '\t');
            localIdxToOriginalIdxMap[counter] = std::stoi(id);

            if (!read_binary_blob(file, maxLengthCompressed, compressedSeqs_local+1ll*counter*maxLengthCompressed)) {
                std::cerr << "Failed to read binary data for sequence: " << name << std::endl;
                break;
            }

            if (clusterToCompressSeqIdxMap.find(cluster_id[std::stoi(id)]) != clusterToCompressSeqIdxMap.end()) {
                clusterToCompressSeqIdxMap[cluster_id[std::stoi(id)]].push_back(counter);
            }
            counter++;
        }
        gzclose(file);


        for(int i=bc*clustersPerBatchFile;i<(bc+1)*clustersPerBatchFile && i<numSequences*4-4;i++){
            if (!isCluster[i]) continue;
            std::cerr << "Handling cluster " << i << " with size " << contains[i].size() << "\n";
            
            // open file
            uint64_t * compressedSeqs_local_per_cluster = new uint64_t[params.backboneSize*maxLengthCompressed];

            int counter=0;
            for (auto &compressIdx: clusterToCompressSeqIdxMap[i]){
                if (counter + 1> numSequences) {
                    std::cerr << "Cluster " << i << " size is larger than backbone size\n";
                    exit(0);
                }
                
                memcpy(compressedSeqs_local_per_cluster+1ll*counter*maxLengthCompressed, compressedSeqs_local+1ll*compressIdx*maxLengthCompressed, 1ll*maxLengthCompressed*sizeof(uint64_t));
                
                h_leafMap[localIdxToOriginalIdxMap[compressIdx]] = counter;
                leafList[counter] = localIdxToOriginalIdxMap[compressIdx];
                counter++;
                
            }

            if(params.in == "m") {
                auto err = cudaMemcpy(msaDeviceArrays.d_compressedSeqsConst, compressedSeqs_local_per_cluster, 1ll*params.backboneSize*maxLengthCompressed*sizeof(uint64_t),cudaMemcpyHostToDevice);
                if (err != cudaSuccess)
                {
                    fprintf(stderr, "Gpu_ERROR: d_hashListConst cudaMemcpy failed!\n");
                    exit(1);
                }
                cudaDeviceSynchronize();
            }
            

            if (contains[i].size() == 0) continue;
            int localCount_=0;

            resetEdgeMaskIndexDC<<<blockNum,threadNum>>>(d_edgeMaskIndex, totalNumSequences*4);
            cudaDeviceSynchronize();
            initializeClusterDC <<<1,1>>>(
                i,
                d_e,
                d_belong,
                d_head,
                d_nxt,
                d_closest_id,
                d_edgeMask,
                d_leafMask,
                d_edgeMaskIndex,
                d_leafMap
            );
            int edgeCount=2, leafCount=10;
            for(auto &leaf:contains[i]) {
                actual_placement++;
                int leaf_idx_in_cluster = h_leafMap[leaf];
                if(params.in == "m"){
                    msaDeviceArrays.distSpecialIDConstructionOnGpuDC(
                        params,
                        localCount_,
                        d_dist,
                        leafCount,
                        d_leafMask,
                        d_leafMap
                    );
                } 
                localCount_++;

                calculateBranchLengthSpecialIDDC <<<blockNum,threadNum>>> (
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
                    edgeCount,
                    d_edgeMask
                );
                auto iter=thrust::min_element(minPos.begin(),minPos.begin()+edgeCount,compare_tupleDC());
                thrust::tuple<int,double,double> smallest=*iter;

                /*
                Update Tree Structure
                */

                int eid=thrust::get<0>(smallest);
                double fracLen=thrust::get<1>(smallest),addLen=thrust::get<2>(smallest);
                // std::cerr << "eid: " << eid << " Cluster ID: " << j << " dist " << addLen << " dist2 " << fracLen << "\n";
                updateTreeStructureInClusterDC <<<1,1>>>(
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
                    leaf,
                    idx,
                    totalNumSequences,
                    insertLeafCount
                );
                idx+=4, insertLeafCount++;

                /*
                Update edgeMask and leafMask
                */

                updateClusterInfoDC<<<1,1>>> (
                    leaf,
                    idx,
                    d_leafMask,
                    d_edgeMask,
                    d_edgeMaskIndex,
                    edgeCount,
                    leafCount,
                    d_leafMap,
                    leaf_idx_in_cluster
                );

                /*
                Update closest nodes
                */

                updateClosestNodesInClusterDC <<<1,1>>> (
                    d_head,
                    d_nxt,
                    d_e,
                    d_len,
                    d_closest_dis,
                    d_closest_id,
                    leaf,
                    d_id,
                    d_from,
                    d_dis,
                    i,
                    d_belong,
                    d_edgeMaskIndex
                );

                edgeCount+=4, leafCount++;
                // std::cerr << "leaf: " << leaf << " Cluster ID: " << j << " eid " << eid << " dist " << addLen << " dist2 " << fracLen << "\n";
                // if (leaf == 879) exit(0);
                // exit(0);
            }
            delete[] compressedSeqs_local_per_cluster;

        }

        delete[] compressedSeqs_local;
    }
    cudaDeviceSynchronize();


    return;
}


void MashPlacement::KPlacementDeviceArraysDC::findBackboneTreeDCRecursive(
    Param& params,
    MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    MSADeviceArraysDC& msaDeviceArrays,
    const KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost,
    int clusterIdx
){ 
    cudaError_t err;
    
    int * d_id;
    err = cudaMalloc(&d_id, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed 1!\n");
        exit(1);
    }
    int * d_from;
    err = cudaMalloc(&d_from, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed 2 !\n");
        exit(1);
    }
    double * d_dis;
    err = cudaMalloc(&d_dis, totalNumSequences*2*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed 3!\n");
        exit(1);
    }

    int * d_edgeMask;
    err = cudaMalloc(&d_edgeMask, totalNumSequences*4*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int * d_edgeMaskIndex;
    err = cudaMalloc(&d_edgeMaskIndex, totalNumSequences*4*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int * d_leafMask;
    err = cudaMalloc(&d_leafMask, totalNumSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    int * d_leafMap;
    err = cudaMalloc(&d_leafMap, params.backboneSize*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    auto * cluster_id = kplacementDeviceArraysHost.clusterID;
    std::vector<std::vector <int>> contains(numSequences*4-4);

    for(int i=numSequences;i<totalNumSequences;i++) contains[cluster_id[i]].push_back(i);

    std::vector<int> leafList (params.backboneSize);
    std::unordered_map <int,int> h_leafMap;

    int newBackboneSize=params.backboneSize;

    for (auto i: contains[clusterIdx]){
        if (newBackboneSize==2*params.backboneSize) break;
        h_leafMap[contains[clusterIdx][i]]=newBackboneSize;
        leafList[newBackboneSize-params.backboneSize]=contains[clusterIdx][i];
        newBackboneSize++;
    }

    if(params.in == "r")
        transferMashClusterInfoDC(mashDeviceArrays, leafList, params);
    else if (params.in == "m")
        transferMsaClusterInfoDC(msaDeviceArrays, leafList, params);
    else {
        std::cerr << "Error: Input type must be unaligned or aligned for clustering based approach\n";
        exit(1);
    }



    // int threadNum = 1024, blockNum = 1024;
    // initializeClusterDC <<<1,1>>>(
    //     clusterIdx,
    //     d_e,
    //     d_belong,
    //     d_head,
    //     d_nxt,
    //     d_closest_id,
    //     d_edgeMask,
    //     d_leafMask,
    //     d_edgeMaskIndex,
    //     d_leafMap
    // );
    
    // std::chrono::nanoseconds disTime(0), treeTime(0);
    // for(int i=bd;i<numSequences;i++){
    //     auto disStart = std::chrono::high_resolution_clock::now();
    //     if(params.in == "r"){
    //         mashDeviceArrays.distRangeConstructionOnGpuDC(
    //             params,
    //             i,
    //             d_dist,
    //             0,
    //             i-1
    //         );
    //     }
    //     else if(params.in == "d"){
    //         matrixReader.distConstructionOnGpu(
    //             params,
    //             i,
    //             d_dist
    //         );
    //     }
    //     else if(params.in == "m"){
    //         msaDeviceArrays.distRangeConstructionOnGpuDC(
    //             params,
    //             i,
    //             d_dist,
    //             0,
    //             i-1
    //         );
    //     }

    //     auto disEnd = std::chrono::high_resolution_clock::now();
    //     auto treeStart = std::chrono::high_resolution_clock::now();
    //     calculateBranchLengthDC <<<blockNum,threadNum>>> (
    //         i,
    //         d_head,
    //         d_nxt,
    //         d_dist,
    //         d_e,
    //         d_len,
    //         d_belong,
    //         thrust::raw_pointer_cast(minPos.data()),
    //         numSequences*4-4,
    //         d_closest_dis,
    //         d_closest_id
    //     );

    //     auto iter=thrust::min_element(minPos.begin(),minPos.begin()+numSequences*4-4,compare_tupleDC());
    //     thrust::tuple<int,double,double> smallest=*iter;
    //     /*
    //     Update Tree (and assign closest nodes to newly added nodes)
    //     */
    //     int eid=thrust::get<0>(smallest);
    //     double fracLen=thrust::get<1>(smallest),addLen=thrust::get<2>(smallest);
    //     updateTreeStructureDC <<<1,1>>>(
    //         d_head,
    //         d_nxt,
    //         d_e,
    //         d_len,
    //         d_closest_dis,
    //         d_closest_id,
    //         d_belong,
    //         eid,
    //         fracLen,
    //         addLen,
    //         i,
    //         idx,
    //         totalNumSequences
    //     );
    //     idx+=4;

    //     /*
    //     Update closest nodes
    //     */
    //     updateClosestNodesDC <<<1,1>>> (
    //         d_head,
    //         d_nxt,
    //         d_e,
    //         d_len,
    //         d_closest_dis,
    //         d_closest_id,
    //         i,
    //         d_id,
    //         d_from,
    //         d_dis
    //     );
    //     auto treeEnd = std::chrono::high_resolution_clock::now();
    //     disTime += disEnd - disStart;
    //     treeTime += treeEnd - treeStart;
    // }

    // cudaDeviceSynchronize();
    // auto backboneEnd = std::chrono::high_resolution_clock::now();
    // auto backboneTime = backboneEnd - backboneStart;
    // std::cerr << "Finished backbone construction in: "<< backboneTime.count()/1000000 << " ms\n";
    
    return;
}
