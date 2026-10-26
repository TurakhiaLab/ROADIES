#include "mash_placement.cuh"

#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cub/cub.cuh>

void MashPlacement::KPlacementDeviceArrays::allocateDeviceArrays(size_t num, int backboneSize)
{
    cudaError_t err;
    numSequences = int(num);
    bd = 2, idx = 0;
    this->backboneSize = backboneSize;
    err = cudaMalloc(&d_dist, numSequences * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_head, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_e, numSequences * 8 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_len, numSequences * 8 * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_nxt, numSequences * 8 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_belong, numSequences * 8 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_dis, numSequences * 20 * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_id, numSequences * 20 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
}

__global__ void initializeID(
    int lim,
    double *d_closest_dis,
    int *d_closest_id)
{
    int tx = threadIdx.x, bs = blockDim.x;
    int bx = blockIdx.x, gs = gridDim.x;
    int idx = tx + bs * bx;
    if (idx < lim)
    {
        for (int i = 0; i < 5; i++)
        {
            d_closest_dis[idx * 5 + i] = 2;
            d_closest_id[idx * 5 + i] = -1;
        }
    }
}

__global__ void updateClosestNodes(
    int *head,
    int *nxt,
    int *e,
    double *len,
    double *closest_dis,
    int *closest_id,
    int x,
    int *id,
    int *from,
    double *dis)
{
    int l = 0, r = -1;
    id[++r] = x, dis[x] = 0, from[x] = -1;
    while (l <= r)
    {
        int node = id[l], fb = from[l];
        double d = dis[l];
        l++;
        for (int i = head[node]; i != -1; i = nxt[i])
        {
            // printf("%d %d: \n", node, head[node]);
            if (e[i] == fb)
                continue;
            for (int j = 0; j < 5; j++)
            {
                double nowd = closest_dis[i * 5 + j];
                if (nowd > d)
                {
                    for (int k = 4; k > j; k--)
                    {
                        closest_dis[i * 5 + k] = closest_dis[i * 5 + k - 1];
                        closest_id[i * 5 + k] = closest_id[i * 5 + k - 1];
                    }
                    // printf("%d: (%d %lf)\t", i*5+j, x, d);
                    closest_dis[i * 5 + j] = d;
                    closest_id[i * 5 + j] = x;
                    id[++r] = e[i], dis[r] = d + len[i], from[r] = node;
                    break;
                }
            }
            // printf("\n");
        }
    }
}

void MashPlacement::KPlacementDeviceArrays::initializeDeviceArrays(Tree *t)
{
    size_t numSequences = this->numSequences;
    size_t totalNodes = t->allNodes.size();

    Node *root = t->root;
    if (root == nullptr)
    {
        fprintf(stderr, "Error: Tree root is null.\n");
        return;
    }
    if (root->children.empty())
    {
        fprintf(stderr, "Error: Tree has no children.\n");
        return;
    }
    if (totalNodes < 2)
    {
        fprintf(stderr, "Error: Tree has less than two nodes.\n");
        return;
    }

    int *h_head = new int[numSequences * 2];
    int *h_e = new int[numSequences * 8];
    int *h_nxt = new int[numSequences * 8];
    int *h_belong = new int[numSequences * 8];
    double *h_len = new double[numSequences * 8];

    // initialize arrays
    for (int i = 0; i < numSequences * 2; i++)
        h_head[i] = -1;
    for (int i = 0; i < numSequences * 8; i++)
    {
        h_e[i] = -1;
        h_nxt[i] = -1;
        h_len[i] = 2;
        h_belong[i] = -1;
    }

    // dfs postorder traversal to initialize the tree structure
    std::function<void(Node *, size_t &)> dfs = [&](Node *node, size_t &edgeCount)
    {
        if (node == nullptr)
            return;
        for (Node *child : node->children)
        {
            dfs(child, edgeCount);
        }
        if (node->parent == nullptr)
            return;
        // if (node->parent->parent == nullptr) return; // skip root's children
        int x = node->idx;
        int y = node->parent->idx;
        // child->parent edge (x->y)
        h_e[edgeCount] = y;
        h_len[edgeCount] = node->bl;
        h_belong[edgeCount] = x;
        h_nxt[edgeCount] = h_head[x];
        h_head[x] = edgeCount;
        edgeCount++;

        // parent->child edge (y->x)
        h_e[edgeCount] = x;
        h_len[edgeCount] = node->bl;
        h_belong[edgeCount] = y;
        h_nxt[edgeCount] = h_head[y];
        h_head[y] = edgeCount;
        edgeCount++;
        // }
    };
    size_t edgeCount = 0;
    dfs(root, edgeCount);

    // transfer data to device
    cudaError_t err;
    err = cudaMemcpy(d_head, h_head, numSequences * 2 * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed for d_head!\n");
        exit(1);
    }
    err = cudaMemcpy(d_e, h_e, numSequences * 8 * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed for d_e!\n");
        exit(1);
    }
    err = cudaMemcpy(d_len, h_len, numSequences * 8 * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed for d_len!\n");
        exit(1);
    }
    err = cudaMemcpy(d_nxt, h_nxt, numSequences * 8 * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed for d_nxt!\n");
        exit(1);
    }
    err = cudaMemcpy(d_belong, h_belong, numSequences * 8 * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed for d_belong!\n");
        exit(1);
    }

    // initialize closest_dis and closest_id
    int *d_id;
    err = cudaMalloc(&d_id, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int *d_from;
    err = cudaMalloc(&d_from, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    double *d_dis;
    err = cudaMalloc(&d_dis, numSequences * 2 * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    initializeID<<<1024, 1024>>>(
        numSequences * 4 - 4,
        d_closest_dis,
        d_closest_id);

    for (int i = 0; i < this->backboneSize; i++)
    {
        updateClosestNodes<<<1, 1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis);
    }

    return;
}

__global__ void initialize(
    int lim,
    int nodes,
    double *d_closest_dis,
    int *d_closest_id,
    int *head,
    int *nxt,
    int *belong,
    int *e)
{
    int tx = threadIdx.x, bs = blockDim.x;
    int bx = blockIdx.x, gs = gridDim.x;
    int idx = tx + bs * bx;
    for (int t = idx; t < lim; t += bs * gs)
    {
        if (t < lim)
        {
            for (int i = 0; i < 5; i++)
            {
                d_closest_dis[t * 5 + i] = 2;
                d_closest_id[t * 5 + i] = -1;
            }
            nxt[t] = -1;
            e[t] = -1;
            belong[t] = -1;
        }
        if (t < nodes)
            head[t] = -1;
    }
}




struct compare_tuple {
  __host__ __device__
  bool operator()(thrust::tuple<int,double,double> lhs, thrust::tuple<int,double,double> rhs)
  {
    return thrust::get<2>(lhs) < thrust::get<2>(rhs);
    //Always find the tuple whose third value (the criteria we want to minimize) is minimized
  }
};


struct find_if_tuple_bidir {
    int bidirNewedgeIdToAttach;
    __host__ __device__
    find_if_tuple_bidir(int id) : bidirNewedgeIdToAttach(id) {}
    __host__ __device__
    bool operator()(thrust::tuple<int,double,double> tup) const
    {
      return thrust::get<0>(tup) == bidirNewedgeIdToAttach;
    }
};
/*
Three variables in tuple:
ID of branch in linked list,
distance to new node inserted on branch from starting vertex (belong[id]),
distance from new node inserted on branch to new node inserted outside branch
*/

__global__ void calculateBranchLength(
    int num, // should be bd, not numSequences
    int *head,
    int *nxt,
    double *dis,
    int *e,
    double *len,
    int *belong,
    thrust::tuple<int, double, double> *minPos,
    int lim,
    double * closest_dis,
    int * closest_id,
    int totSeqNum
){
    int tx=threadIdx.x,bs=blockDim.x,bx=blockIdx.x,gs=gridDim.x;
    int idx_=tx+bs*bx;
    for (int idx=idx_; idx<lim; idx+=bs*gs){
        if(idx>=lim) return;
        if(idx>=num*4-4||belong[idx]<e[idx]){
            thrust::tuple <int,double,double> minTuple(0,0,2);
            minPos[idx]=minTuple;
            continue;
        }
        int x = belong[idx], oth = e[idx];
        int eid = idx, otheid;
        double dis1 = 0, dis2 = 0, val;
        for (int i = 0; i < 5; i++)
            if (closest_id[eid * 5 + i] != -1)
            {
                val = dis[closest_id[eid * 5 + i]] - closest_dis[eid * 5 + i];
                if (val > dis1)
                    dis1 = val;
            }
        otheid=head[oth];
        /* modified for fixed topology */
        while(otheid!=-1 && e[otheid]!=x) otheid=nxt[otheid];
        if(otheid==-1){
            thrust::tuple <int,double,double> minTuple(0,0,2);
            minPos[idx]=minTuple;
            continue;
        }
        for(int i=0;i<5;i++)
            if(closest_id[otheid*5+i]!=-1){
                val = dis[closest_id[otheid*5+i]]-closest_dis[otheid*5+i];
                if(val>dis2) dis2=val;
            }
        double additional_dis = (dis1 + dis2 - len[eid]) / 2;
        if (additional_dis < 0)
            additional_dis = 0;
        dis1 -= additional_dis, dis2 -= additional_dis;
        if (dis1 < 0)
            dis1 = 0;
        if (dis2 < 0)
            dis2 = 0;
        if (dis1 > len[eid])
            additional_dis += dis1 - len[eid], dis1 = len[eid];
        if (dis2 > len[eid])
            additional_dis += dis2 - len[eid], dis2 = len[eid];
        // assert(dis1+dis2-1e-6<=len[eid]);
        double rest=len[eid]-dis1-dis2;
        dis1+=rest/2,dis2+=rest/2;
        // print eif, dis1, and additional_dis
        // printf("eid %d (nodes %d-%d): dis1=%.8lf, dis2=%.8lf, len=%.8lf, add=%.8lf\n", eid, x, oth, dis1, dis2, len[eid], additional_dis);
        thrust::tuple <int,double,double> minTuple(eid,dis1,additional_dis);
        minPos[idx]=minTuple;
    }
}

__global__ void updateTreeStructuretoAddQuery(
    int *head,
    int *nxt,
    int *e,
    double *len,
    double *closest_dis,
    int *closest_id,
    int *belong,
    int eid,
    double fracLen,
    double addLen,
    int placeId,   // Id of the newly placed node
    int edgeCount, // Position to insert a new edge in linked list
    int numSequences)
{
    int middle = placeId + numSequences - 1, outside = placeId;
    int x = belong[eid], y = e[eid];
    double originalDis = len[eid];
    int xe, ye;
    for (int i = head[x]; i != -1; i = nxt[i])
        if (e[i] == y)
        {
            e[i] = middle, len[i] = fracLen, xe = i;
            break;
        }
    for (int i = head[y]; i != -1; i = nxt[i])
        if (e[i] == x)
        {
            e[i] = middle, len[i] -= fracLen, ye = i;
            break;
        }
    /*
    Need to update:
    e, len, nxt, head, belong, closest_dis, closest_id
    */
    // middle -> x
    e[edgeCount] = x, len[edgeCount] = fracLen, nxt[edgeCount] = head[middle], head[middle] = edgeCount, belong[edgeCount] = middle;
    for (int i = 0; i < 5; i++)
        if (closest_id[ye * 5 + i] != -1)
        {
            closest_id[edgeCount * 5 + i] = closest_id[ye * 5 + i];
            closest_dis[edgeCount * 5 + i] = closest_dis[ye * 5 + i] + originalDis - fracLen;
        }
    edgeCount++;
    // middle -> y
    e[edgeCount] = y, len[edgeCount] = originalDis - fracLen, nxt[edgeCount] = head[middle], head[middle] = edgeCount, belong[edgeCount] = middle;
    for (int i = 0; i < 5; i++)
        if (closest_id[xe * 5 + i] != -1)
        {
            closest_id[edgeCount * 5 + i] = closest_id[xe * 5 + i];
            closest_dis[edgeCount * 5 + i] = closest_dis[xe * 5 + i] + fracLen;
        }
    edgeCount++;
    // outside -> middle
    e[edgeCount] = middle, len[edgeCount] = addLen, nxt[edgeCount] = head[outside], head[outside] = edgeCount, belong[edgeCount] = outside;
    edgeCount++;
    // middle -> outside
    e[edgeCount] = outside, len[edgeCount] = addLen, nxt[edgeCount] = head[middle], head[middle] = edgeCount, belong[edgeCount] = middle;
    int e1 = edgeCount - 2, e2 = edgeCount - 3;
    for (int i = 0; i < 5; i++)
    {
        if (closest_id[e1 * 5 + i] == -1)
            break;
        for (int j = 0; j < 5; j++)
            if (closest_dis[edgeCount * 5 + j] > closest_dis[e1 * 5 + i])
            {
                for (int k = 4; k > j; k--)
                {
                    closest_dis[edgeCount * 5 + k] = closest_dis[edgeCount * 5 + k - 1];
                    closest_id[edgeCount * 5 + k] = closest_id[edgeCount * 5 + k - 1];
                }
                closest_dis[edgeCount * 5 + j] = closest_dis[e1 * 5 + i];
                closest_id[edgeCount * 5 + j] = closest_id[e1 * 5 + i];
                break;
            }
    }
    for (int i = 0; i < 5; i++)
    {
        if (closest_id[e2 * 5 + i] == -1)
            break;
        for (int j = 0; j < 5; j++)
            if (closest_dis[edgeCount * 5 + j] > closest_dis[e2 * 5 + i])
            {
                for (int k = 4; k > j; k--)
                {
                    closest_dis[edgeCount * 5 + k] = closest_dis[edgeCount * 5 + k - 1];
                    closest_id[edgeCount * 5 + k] = closest_id[edgeCount * 5 + k - 1];
                }
                closest_dis[edgeCount * 5 + j] = closest_dis[e2 * 5 + i];
                closest_id[edgeCount * 5 + j] = closest_id[e2 * 5 + i];
                break;
            }
    }
    edgeCount++;
}

__global__ void updateTreeStructure(
    int *head,
    int *nxt,
    int *e,
    double *len,
    double *closest_dis,
    int *closest_id,
    int *belong,
    int eid,
    double fracLen,
    double addLen,
    int placeId,   // Id of the newly placed node
    int edgeCount, // Position to insert a new edge in linked list
    int numSequences
){
    int middle=placeId+numSequences-1, outside=placeId;
    int x=belong[eid],y=e[eid];
    double originalDis=len[eid];
    int xe,ye;
    for(int i=head[x];i!=-1;i=nxt[i])
        if(e[i]==y){
            e[i]=middle,len[i]=fracLen,xe=i;
            // printf("x -> middle: eid: %d\n", i);
            break;
        }
    for(int i=head[y];i!=-1;i=nxt[i])
        if(e[i]==x){
            e[i]=middle,len[i]-=fracLen,ye=i;
            // printf("y -> middle: eid: %d\n", i);
            break;
        }
    /*
    Need to update:
    e, len, nxt, head, belong, closest_dis, closest_id
    */
    //middle -> x
    e[edgeCount]=x,len[edgeCount]=fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    // printf("middle -> x: eid: %d\n", edgeCount);
    for(int i=0;i<5;i++)
        if(closest_id[ye*5+i]!=-1){
            closest_id[edgeCount*5+i]=closest_id[ye*5+i];
            closest_dis[edgeCount*5+i]=closest_dis[ye*5+i]+originalDis-fracLen;
        }
    edgeCount++;
    //middle -> y
    e[edgeCount]=y,len[edgeCount]=originalDis-fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    // printf("middle -> y: eid: %d\n", edgeCount);
    for(int i=0;i<5;i++)
        if(closest_id[xe*5+i]!=-1){
            closest_id[edgeCount*5+i]=closest_id[xe*5+i];
            closest_dis[edgeCount*5+i]=closest_dis[xe*5+i]+fracLen;
        }
    edgeCount++;
    //outside -> middle
    e[edgeCount]=middle,len[edgeCount]=addLen,nxt[edgeCount]=head[outside],head[outside]=edgeCount,belong[edgeCount]=outside;
    // printf("outside -> middle: eid: %d\n", edgeCount);
    edgeCount++;
    //middle -> outside
    e[edgeCount]=outside,len[edgeCount]=addLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    // printf("middle -> outside: eid: %d\n", edgeCount);
    int e1=edgeCount-2, e2=edgeCount-3;
    for(int i=0;i<5;i++){
        if(closest_id[e1*5+i]==-1) break;
        for(int j=0;j<5;j++)
            if(closest_dis[edgeCount*5+j]>closest_dis[e1*5+i]){
                for(int k=4;k>j;k--){
                    closest_dis[edgeCount*5+k]=closest_dis[edgeCount*5+k-1];
                    closest_id[edgeCount*5+k]=closest_id[edgeCount*5+k-1];
                }
                closest_dis[edgeCount * 5 + j] = closest_dis[e1 * 5 + i];
                closest_id[edgeCount * 5 + j] = closest_id[e1 * 5 + i];
                break;
            }
    }
    for (int i = 0; i < 5; i++)
    {
        if (closest_id[e2 * 5 + i] == -1)
            break;
        for (int j = 0; j < 5; j++)
            if (closest_dis[edgeCount * 5 + j] > closest_dis[e2 * 5 + i])
            {
                for (int k = 4; k > j; k--)
                {
                    closest_dis[edgeCount * 5 + k] = closest_dis[edgeCount * 5 + k - 1];
                    closest_id[edgeCount * 5 + k] = closest_id[edgeCount * 5 + k - 1];
                }
                closest_dis[edgeCount * 5 + j] = closest_dis[e2 * 5 + i];
                closest_id[edgeCount * 5 + j] = closest_id[e2 * 5 + i];
                break;
            }
    }
    edgeCount++;
}

__global__ void buildInitialTree(
    int numSequences,
    int *head,
    int *e,
    double *len,
    int *nxt,
    int *belong,
    double *dis,
    int edgeCount)
{
    int nv = numSequences;
    double d = dis[0];
    // 0 -> nv
    e[edgeCount]=nv,len[edgeCount]=d/2,nxt[edgeCount]=head[0],head[0]=edgeCount,belong[edgeCount]=0;
    // printf("0 -> nv: eid: %d, from: %d, to: %d\n", edgeCount);
    edgeCount++;
    // 1 -> nv
    e[edgeCount]=nv,len[edgeCount]=d/2,nxt[edgeCount]=head[1],head[1]=edgeCount,belong[edgeCount]=1;
    // printf("1 -> nv: eid: %d\n", edgeCount);
    edgeCount++;
    // nv -> 0
    e[edgeCount] = 0, len[edgeCount] = d / 2, nxt[edgeCount] = head[nv], head[nv] = edgeCount, belong[edgeCount] = nv;
    edgeCount++;
    // nv -> 0
    e[edgeCount] = 1, len[edgeCount] = d / 2, nxt[edgeCount] = head[nv], head[nv] = edgeCount, belong[edgeCount] = nv;
    edgeCount++;
}

void MashPlacement::KPlacementDeviceArrays::deallocateDeviceArrays()
{
    cudaFree(d_head);
    cudaFree(d_e);
    cudaFree(d_nxt);
    cudaFree(d_belong);
    cudaFree(d_closest_id);
    cudaFree(d_dist);
    cudaFree(d_len);
    cudaFree(d_closest_dis);
}


void MashPlacement::KPlacementDeviceArrays::printTree(std::vector <std::string> name, std::ofstream& output_){
    int    * h_head = new int[numSequences*2];
    int    * h_e    = new int[numSequences*8];
    int    * h_nxt  = new int[numSequences*8];
    double * h_len  = new double[numSequences*8];

    auto checkCopy = [](cudaError_t err) {
        if (err != cudaSuccess) {
            fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
            exit(1);
        }
    };
    checkCopy(cudaMemcpy(h_head, d_head, numSequences*2*sizeof(int),    cudaMemcpyDeviceToHost));
    checkCopy(cudaMemcpy(h_e,    d_e,    numSequences*8*sizeof(int),    cudaMemcpyDeviceToHost));
    checkCopy(cudaMemcpy(h_nxt,  d_nxt,  numSequences*8*sizeof(int),    cudaMemcpyDeviceToHost));
    checkCopy(cudaMemcpy(h_len,  d_len,  numSequences*8*sizeof(double), cudaMemcpyDeviceToHost));

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

    if (!g_printBinaryNewick) {
        print(numSequences + bd - 2, -1);
    } else {
        print_binary(numSequences + bd - 2, -1);
    }
    output_ << ";";

    delete[] h_head;
    delete[] h_e;
    delete[] h_nxt;
    delete[] h_len;
}

void MashPlacement::KPlacementDeviceArrays::printTree(std::vector <std::string> name, std::ofstream& output_, UnrootedTree* t, std::vector<int>& edgeIdsMappingToNewTreeEdges){
    double * h_len = new double[numSequences*8];
    auto err = cudaMemcpy(h_len, d_len, numSequences*8*sizeof(double),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }

    // Iterate over all edges in the tree and update the edge length based on edgeIdsMappingToNewTreeEdges
    for (auto &nodePair : t->getNodes()) {
        UnrootedNode* node = nodePair.second;
        for (auto &edge : node->neighbors) {
            int edge_id = edge.edge_id;
            assert(edge_id >= 0 && edge_id < numSequences * 8 && "Edge ID out of bounds");
            edge.length = h_len[edge_id];
        }
    }

    if (!g_printBinaryNewick) {
        t->collapseZeroLengthInternalEdges();
    }
    output_ << t->toNewick() << std::endl;
    return;

}

void MashPlacement::KPlacementDeviceArrays::findPlacementTree(
    Param &params,
    const MashDeviceArrays &mashDeviceArrays,
    MatrixReader &matrixReader,
    const MSADeviceArrays &msaDeviceArrays)
{
    if (params.in == "d")
    {
        matrixReader.distConstructionOnGpu(params, 0, d_dist);
    }
    int *d_id;
    auto err = cudaMalloc(&d_id, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int *d_from;
    err = cudaMalloc(&d_from, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    double *d_dis;
    err = cudaMalloc(&d_dis, numSequences * 2 * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    /*
    Initialize closest nodes by inifinite
    */
    // int threadNum = 256, blockNum = (numSequences*4-4+threadNum-1)/threadNum;
    int threadNum = 1024, blockNum = 1024;
    initialize<<<blockNum, threadNum>>>(
        numSequences * 4 - 4,
        numSequences * 2,
        d_closest_dis,
        d_closest_id,
        d_head,
        d_nxt,
        d_belong,
        d_e);
    /*
    Build Initial Tree
    */
    if (params.in == "r")
    {
        mashDeviceArrays.distConstructionOnGpu(
            params,
            1,
            d_dist);
    }
    else if (params.in == "d")
    {
        matrixReader.distConstructionOnGpu(
            params,
            1,
            d_dist);
    }
    else if (params.in == "m")
    {
        msaDeviceArrays.distConstructionOnGpu(
            params,
            1,
            d_dist);
    }
    // cudaDeviceSynchronize();

    // return;
    // double * h_dis = new double[numSequences];
    // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
    // for(int j=0;j<1;j++) fprintf(stderr,"%.8lf ",h_dis[j]);std::cerr<<'\n';

    buildInitialTree <<<1,1>>> (
        numSequences,
        d_head,
        d_e,
        d_len,
        d_nxt,
        d_belong,
        d_dist,
        idx);
    idx += 4;
    /*
    Initialize closest nodes by inital tree
    */
    for (int i = 0; i < bd; i++)
    {
        updateClosestNodes<<<1, 1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis);
    }

    thrust::device_vector<thrust::tuple<int, double, double>> minPos(numSequences * 4 - 4);
    // std::cout<<"FFF\n";
    std::chrono::nanoseconds disTime(0), treeTime(0);
    for (int i = bd; i < numSequences; i++)
    {
        auto disStart = std::chrono::high_resolution_clock::now();
        // blockNum = (i + 255) / 256;
        // blockNum = 1024;
        if (params.in == "r")
        {
            mashDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        else if (params.in == "d")
        {
            matrixReader.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        else if (params.in == "m")
        {
            msaDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        cudaDeviceSynchronize();

        // double * h_dis = new double[numSequences];
        // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
        // fprintf(stderr, "%d\n",i);
        // for(int j=0;j<i;j++) std::cerr<<h_dis[j]<<" ";std::cerr<<'\n';

        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        
        calculateBranchLength <<<blockNum,threadNum>>> (
            i,
            d_head,
            d_nxt,
            d_dist,
            d_e,
            d_len,
            d_belong,
            thrust::raw_pointer_cast(minPos.data()),
            numSequences * 4 - 4,
            d_closest_dis,
            d_closest_id,
            numSequences
        );
        
        auto iter=thrust::min_element(minPos.begin(),minPos.end(),compare_tuple());
        thrust::tuple<int,double,double> smallest=*iter;
        // /* print top 5 sorted elements */
        // for(int j=0;j<5;j++){
        //     thrust::tuple<int,double,double> s=*iter;
        //     std::cerr<<thrust::get<0>(s)<<" "<<thrust::get<1>(s)<<" "<<thrust::get<2>(s)<<'\n';
        //     iter++;
        // } 
        
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid = thrust::get<0>(smallest);
        double fracLen = thrust::get<1>(smallest), addLen = thrust::get<2>(smallest);
        updateTreeStructure<<<1, 1>>>(
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
            numSequences);
        idx += 4;
        /*
        Update closest nodes
        */
        updateClosestNodes<<<1, 1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis);
        // cudaDeviceSynchronize();
        auto treeEnd = std::chrono::high_resolution_clock::now();
        disTime += disEnd - disStart;
        treeTime += treeEnd - treeStart;
        // std::cerr << "Seq " << i << " at " << eid << " with fracLen " << fracLen 
        //           << " and addLen " << addLen << std::endl;
    }
    std::cerr << "Distance Operation Time " << disTime.count() / 1000000 << " ms\n";
    std::cerr << "Tree Operation Time " << treeTime.count() / 1000000 << " ms\n";
}





void MashPlacement::KPlacementDeviceArrays::addQuery(
    Param &params,
    const MashDeviceArrays &mashDeviceArrays,
    MatrixReader &matrixReader,
    const MSADeviceArrays &msaDeviceArrays)
{

    int *d_id;
    auto err = cudaMalloc(&d_id, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int *d_from;
    err = cudaMalloc(&d_from, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    double *d_dis;
    err = cudaMalloc(&d_dis, numSequences * 2 * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int threadNum = 256, blockNum = (numSequences * 4 - 4 + threadNum - 1) / threadNum;
    idx += 4 * this->backboneSize - 4; // Adjust idx to account for the backbone tree size

    thrust::device_vector<thrust::tuple<int, double, double>> minPos(numSequences * 4 - 4);
    std::chrono::nanoseconds disTime(0), treeTime(0);
    for (int i = this->backboneSize; i < numSequences; i++)
    {
        auto disStart = std::chrono::high_resolution_clock::now();
        blockNum = (i + 255) / 256;
        // blockNum = 1024;
        if (params.in == "r")
        {
            mashDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        else if (params.in == "d")
        {
            matrixReader.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        else if (params.in == "m")
        {
            msaDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        cudaDeviceSynchronize();

        // double * h_dis = new double[numSequences];
        // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
        // fprintf(stderr, "%d\n",i);
        // for(int j=0;j<i;j++) std::cerr<<h_dis[j]<<" ";std::cerr<<'\n';

        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        blockNum = (numSequences * 4 - 4 + 255) / 256;
        // blockNum = 1024;
        calculateBranchLength<<<blockNum, threadNum>>>(
            i,
            d_head,
            d_nxt,
            d_dist,
            d_e,
            d_len,
            d_belong,
            thrust::raw_pointer_cast(minPos.data()),
            numSequences * 4 - 4,
            d_closest_dis,
            d_closest_id,
            numSequences);
        auto iter = thrust::min_element(minPos.begin(), minPos.end(), compare_tuple());
        thrust::tuple<int, double, double> smallest = *iter;
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid = thrust::get<0>(smallest);
        double fracLen = thrust::get<1>(smallest), addLen = thrust::get<2>(smallest);
        updateTreeStructure<<<1, 1>>>(
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
            numSequences);
        idx += 4;
        /*
        Update closest nodes
        */
        updateClosestNodes<<<1, 1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis);
        // cudaDeviceSynchronize();
        auto treeEnd = std::chrono::high_resolution_clock::now();
        disTime += disEnd - disStart;
        treeTime += treeEnd - treeStart;
    }
    std::cerr << "Distance Operation Time " <<  disTime.count()/1000000 << " ms\n";
    std::cerr << "Tree Operation Time " <<  treeTime.count()/1000000 << " ms\n";

}

// ---- Iterative (one-at-a-time) placement API --------------------------------

void MashPlacement::KPlacementDeviceArrays::beginIterativePlacement()
{
    cudaError_t err;
    err = cudaMalloc(&d_id_iter, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess) { fprintf(stderr, "Gpu_ERROR: cudaMalloc failed (d_id_iter)!\n"); exit(1); }
    err = cudaMalloc(&d_from_iter, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess) { fprintf(stderr, "Gpu_ERROR: cudaMalloc failed (d_from_iter)!\n"); exit(1); }
    err = cudaMalloc(&d_dis_iter, numSequences * 2 * sizeof(double));
    if (err != cudaSuccess) { fprintf(stderr, "Gpu_ERROR: cudaMalloc failed (d_dis_iter)!\n"); exit(1); }

    /* Mirror the idx adjustment from addQuery: skip the backbone's edge slots. */
    idx += 4 * this->backboneSize - 4;
    m_nextQueryIdx = this->backboneSize;
}

void MashPlacement::KPlacementDeviceArrays::beginIterativePlacementFromScratch(
    Param& params,
    const MashDeviceArrays& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArrays& msaDeviceArrays)
{
    cudaError_t err;
    err = cudaMalloc(&d_id_iter, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess) { fprintf(stderr, "Gpu_ERROR: cudaMalloc failed (d_id_iter)!\n"); exit(1); }
    err = cudaMalloc(&d_from_iter, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess) { fprintf(stderr, "Gpu_ERROR: cudaMalloc failed (d_from_iter)!\n"); exit(1); }
    err = cudaMalloc(&d_dis_iter, numSequences * 2 * sizeof(double));
    if (err != cudaSuccess) { fprintf(stderr, "Gpu_ERROR: cudaMalloc failed (d_dis_iter)!\n"); exit(1); }

    /* Initialize k-closest arrays and adjacency lists to sentinel values. */
    int threadNum = 1024, blockNum = 1024;
    initialize<<<blockNum, threadNum>>>(
        numSequences * 4 - 4,
        numSequences * 2,
        d_closest_dis,
        d_closest_id,
        d_head,
        d_nxt,
        d_belong,
        d_e);

    /* Compute distances from seq 1 to seq 0 to seed the initial two-node tree. */
    if (params.in == "r") {
        mashDeviceArrays.distConstructionOnGpu(params, 1, d_dist);
    } else if (params.in == "d") {
        matrixReader.distConstructionOnGpu(params, 1, d_dist);
    } else if (params.in == "m") {
        msaDeviceArrays.distConstructionOnGpu(params, 1, d_dist);
    }

    /* Build the initial two-sequence tree (seq 0 and seq 1). */
    buildInitialTree<<<1, 1>>>(
        numSequences,
        d_head,
        d_e,
        d_len,
        d_nxt,
        d_belong,
        d_dist,
        idx);
    idx += 4;   /* bd == 2: two nodes × 2 directed edges each. */

    /* Seed the k-closest lists from the initial two-node tree. */
    for (int i = 0; i < bd; i++) {
        updateClosestNodes<<<1, 1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id_iter,
            d_from_iter,
            d_dis_iter);
    }

    /* All sequences from index bd onward are queries. */
    m_nextQueryIdx = bd;
}

// -----------------------------------------------------------------------------
// Two-phase placement: findBestEdge / getEdgeInfo / findEdgeBetween / commitQuery
// -----------------------------------------------------------------------------

MashPlacement::KPlacementDeviceArrays::PlacementInfo
MashPlacement::KPlacementDeviceArrays::findBestEdge(
    int seqIdx,
    Param& params,
    const MashDeviceArrays& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArrays& msaDeviceArrays)
{
    /* Compute pairwise distances from seqIdx to all currently-placed nodes. */
    if (params.in == "r")
        mashDeviceArrays.distConstructionOnGpu(params, seqIdx, d_dist);
    else if (params.in == "d")
        matrixReader.distConstructionOnGpu(params, seqIdx, d_dist);
    else if (params.in == "m")
        msaDeviceArrays.distConstructionOnGpu(params, seqIdx, d_dist);
    cudaDeviceSynchronize();

    /* Run calculateBranchLength over every edge slot. */
    const int edgeCount = numSequences * 4 - 4;
    int threadNum = 256, blockNum = 1024;
    thrust::device_vector<thrust::tuple<int, double, double>> minPos(edgeCount);

    calculateBranchLength<<<blockNum, threadNum>>>(
        seqIdx,
        d_head, d_nxt, d_dist, d_e, d_len, d_belong,
        thrust::raw_pointer_cast(minPos.data()),
        edgeCount,
        d_closest_dis, d_closest_id,
        numSequences);

    /* Find globally optimal edge. */
    auto iter = thrust::min_element(minPos.begin(), minPos.end(), compare_tuple());
    thrust::tuple<int, double, double> smallest = *iter;
    m_pendingEid     = thrust::get<0>(smallest);
    m_pendingFracLen = thrust::get<1>(smallest);
    m_pendingAddLen  = thrust::get<2>(smallest);
    m_pendingSeqIdx  = seqIdx;

    /* Save host-side copy of fracLen/addLen for every edge slot so that
     * commitQuery(overrideEid) can serve arbitrary edge overrides. */
    {
        std::vector<thrust::tuple<int, double, double>> h_minPos(edgeCount);
        cudaMemcpy(h_minPos.data(),
                   thrust::raw_pointer_cast(minPos.data()),
                   edgeCount * sizeof(thrust::tuple<int, double, double>),
                   cudaMemcpyDeviceToHost);
        m_hostFracLen.resize(edgeCount);
        m_hostAddLen.resize(edgeCount);
        for (int i = 0; i < edgeCount; ++i) {
            m_hostFracLen[i] = thrust::get<1>(h_minPos[i]);
            m_hostAddLen[i]  = thrust::get<2>(h_minPos[i]);
        }
    }

    return getEdgeInfo(m_pendingEid, seqIdx);
}

MashPlacement::KPlacementDeviceArrays::PlacementInfo
MashPlacement::KPlacementDeviceArrays::getEdgeInfo(int eid, int seqIdx) const
{
    int    h_nodeA, h_nodeB;
    double h_origLen;
    cudaMemcpy(&h_nodeA,   d_belong + eid, sizeof(int),    cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_nodeB,   d_e      + eid, sizeof(int),    cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_origLen, d_len    + eid, sizeof(double), cudaMemcpyDeviceToHost);

    const double frac = m_hostFracLen[eid];
    const double add  = m_hostAddLen[eid];

    PlacementInfo info{};
    info.splitNodeA      = h_nodeA;
    info.splitNodeB      = h_nodeB;
    info.internalNodeIdx = seqIdx + numSequences - 1;
    info.tipNodeIdx      = seqIdx;
    info.lenA            = frac;
    info.lenB            = h_origLen - frac;
    info.lenTip          = add;
    return info;
}

std::vector<std::pair<MashPlacement::KPlacementDeviceArrays::PlacementInfo, double>>
MashPlacement::KPlacementDeviceArrays::getTopKEdgeInfo(int seqIdx, int k) const
{
    if (k <= 0 || seqIdx != m_pendingSeqIdx || m_hostAddLen.empty()) {
        return {};
    }
    const int active_edge_slots = seqIdx * 4 - 4;
    std::vector<int> h_belong(static_cast<size_t>(active_edge_slots));
    std::vector<int> h_e(static_cast<size_t>(active_edge_slots));
    cudaMemcpy(
        h_belong.data(), d_belong,
        static_cast<size_t>(active_edge_slots) * sizeof(int),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(
        h_e.data(), d_e,
        static_cast<size_t>(active_edge_slots) * sizeof(int),
        cudaMemcpyDeviceToHost);
    std::vector<int> ranked_eids;
    ranked_eids.reserve(static_cast<size_t>(active_edge_slots / 2));
    for (int eid = 0; eid < active_edge_slots; ++eid) {
        if (h_belong[static_cast<size_t>(eid)] >= 0 &&
            h_e[static_cast<size_t>(eid)] >= 0 &&
            h_belong[static_cast<size_t>(eid)] >= h_e[static_cast<size_t>(eid)] &&
            std::isfinite(m_hostAddLen[static_cast<size_t>(eid)])) {
            ranked_eids.push_back(eid);
        }
    }
    std::stable_sort(
        ranked_eids.begin(), ranked_eids.end(),
        [&](int lhs, int rhs) {
            return m_hostAddLen[static_cast<size_t>(lhs)] <
                m_hostAddLen[static_cast<size_t>(rhs)];
        });
    if (ranked_eids.size() > static_cast<size_t>(k)) {
        ranked_eids.resize(static_cast<size_t>(k));
    }
    std::vector<std::pair<PlacementInfo, double>> out;
    out.reserve(ranked_eids.size());
    for (int eid : ranked_eids) {
        out.emplace_back(
            getEdgeInfo(eid, seqIdx),
            m_hostAddLen[static_cast<size_t>(eid)]);
    }
    return out;
}

int MashPlacement::KPlacementDeviceArrays::findEdgeBetween(int nodeA, int nodeB) const
{
    /* Copy the full d_belong and d_e arrays to host and search linearly.
     * These arrays are sized numSequences*8, typically a few thousand ints. */
    const int edgeSlots = numSequences * 8;
    const int edgeCount = numSequences * 4 - 4;
    std::vector<int> h_belong(edgeSlots), h_e(edgeSlots);
    cudaMemcpy(h_belong.data(), d_belong, edgeSlots * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_e.data(),      d_e,      edgeSlots * sizeof(int), cudaMemcpyDeviceToHost);

    /* Return the canonical direction eid where belong[eid] >= e[eid].
     * This is the direction that calculateBranchLength writes valid results for. */
    for (int i = 0; i < edgeCount; ++i) {
        if (((h_belong[i] == nodeA && h_e[i] == nodeB) ||
             (h_belong[i] == nodeB && h_e[i] == nodeA)) &&
            h_belong[i] >= h_e[i])
        {
            return i;
        }
    }
    return -1;
}

void MashPlacement::KPlacementDeviceArrays::commitQuery(int seqIdx, int overrideEid)
{
    const int    eid     = (overrideEid >= 0) ? overrideEid : m_pendingEid;
    const double fracLen = (overrideEid >= 0) ? m_hostFracLen[eid] : m_pendingFracLen;
    const double addLen  = (overrideEid >= 0) ? m_hostAddLen[eid]  : m_pendingAddLen;

    updateTreeStructure<<<1, 1>>>(
        d_head, d_nxt, d_e, d_len,
        d_closest_dis, d_closest_id, d_belong,
        eid, fracLen, addLen,
        seqIdx, idx, numSequences);
    idx += 4;

    updateClosestNodes<<<1, 1>>>(
        d_head, d_nxt, d_e, d_len,
        d_closest_dis, d_closest_id,
        seqIdx,
        d_id_iter, d_from_iter, d_dis_iter);

    m_nextQueryIdx = seqIdx + 1;
}

// -----------------------------------------------------------------------------

bool MashPlacement::KPlacementDeviceArrays::placeSingleQuery(
    int seqIdx,
    Param& params,
    const MashDeviceArrays& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArrays& msaDeviceArrays,
    PlacementInfo* out)
{
    if (seqIdx >= numSequences) return false;
    PlacementInfo info = findBestEdge(seqIdx, params, mashDeviceArrays, matrixReader, msaDeviceArrays);
    if (out) *out = info;
    commitQuery(seqIdx);
    return m_nextQueryIdx < numSequences;
}

void MashPlacement::KPlacementDeviceArrays::endIterativePlacement()
{
    if (d_id_iter)   { cudaFree(d_id_iter);   d_id_iter   = nullptr; }
    if (d_from_iter) { cudaFree(d_from_iter); d_from_iter = nullptr; }
    if (d_dis_iter)  { cudaFree(d_dis_iter);  d_dis_iter  = nullptr; }
    m_nextQueryIdx = -1;
}

// -----------------------------------------------------------------------------

void MashPlacement::KPlacementDeviceArrays::estimateBranchLengthsFromTopology(
    Param& params,
    const MashDeviceArrays& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArrays& msaDeviceArrays,
    std::vector<int>& edgeMapOldToNew,
    UnrootedTree* t,
    std::vector<std::string>& names
){

    std::unordered_map<int, int> bidirEdgeMap;
    std::unordered_map<int, std::vector<int>> edgeMapNewToOld;
    std::unordered_map<int, bool> visistedEdges;
    std::unordered_map<int, bool> visistedNodes;

    
    int * d_id;
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

    
    int threadNum = 1024, blockNum = 1024;
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
    if(params.in == "r"){
        mashDeviceArrays.distConstructionOnGpu(
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
        msaDeviceArrays.distConstructionOnGpu(
            params,
            1,
            d_dist
        );
    }
    // cudaDeviceSynchronize();

    // return;
    // double * h_dis = new double[numSequences];
    // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
    // for(int j=0;j<1;j++) fprintf(stderr,"%.8lf ",h_dis[j]);std::cerr<<'\n';

    // std::cout << "==================== Building initial tree ====================" << std::endl;
    // std::cout << names[0] << " " << names[1] << std::endl;

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
    cudaDeviceSynchronize();
    idx += 4;

    /* edges 0 and 2 are same edge in different directions */
    std::vector<UnrootedEdge> edges = t->edgesBetween(names[0], t->getRoot());
    edgeMapNewToOld[0]=std::vector<int>();
    for(auto &e: edges){
        visistedEdges[e.edge_id]=true;
        edgeMapOldToNew[e.edge_id]=0;
        edgeMapNewToOld[0].push_back(e.edge_id);
    }
    edges = t->edgesBetween(names[1], t->getRoot());
    edgeMapNewToOld[1]=std::vector<int>();
    for(auto &e: edges){
        visistedEdges[e.edge_id]=true;
        edgeMapOldToNew[e.edge_id]=1;
        edgeMapNewToOld[1].push_back(e.edge_id);
    }
    
    /* edges 0 and 2 are same edge in different directions (based on buildInitialTree) */
    bidirEdgeMap[0]=2;
    bidirEdgeMap[2]=0;
    /* edges 1 and 3 are same edge in different directions (based on buildInitialTree) */
    bidirEdgeMap[1]=3;
    bidirEdgeMap[3]=1;

    // std::cout << "bidirEdgeMap: " << std::endl;
    // for(auto &e: bidirEdgeMap){
    //     std::cout << "edge " << e.first << " -> " << e.second << std::endl;
    // }

    // std::cout << "edgeMapNewToOld: " << std::endl;
    // for(auto &e: edgeMapNewToOld){
    //     std::cout << "edge " << e.first << " -> ";
    //     for(auto &f: e.second){
    //         std::cout << f << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "edgeMapOldToNew: " << std::endl;
    // for (int i=0;i<edgeMapOldToNew.size();i++){
    //     std::cout << "edge " << i << " -> " << edgeMapOldToNew[i] << std::endl;
    // }
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

    for(int i=bd;i<numSequences;i++){
        std::string queryName = names[i];
        std::vector<UnrootedEdge> edges = t->edgesBetween(names[i], names[i-1]);

        int oldedgeIdToAttach = -1;
        int newedgeIdToAttach = -1;
        for (auto &e: edges){
            if (visistedEdges.find(e.edge_id) != visistedEdges.end() && visistedEdges[e.edge_id]) {
                oldedgeIdToAttach = e.edge_id;
                break;
            }
        }
        if(oldedgeIdToAttach == -1){
            std::cerr << "Error: Edge not found" << std::endl;
            exit(1);
        }
        // std::cout << "oldedgeIdToAttach: " << oldedgeIdToAttach << std::endl;

        /* update edgeMapNewToOld */
        int bidirNewedgeIdToAttach=-1;
        assert(edgeMapOldToNew[oldedgeIdToAttach] != -1 && "Edge not found in edgeMapOldToNew");
        newedgeIdToAttach=edgeMapOldToNew[oldedgeIdToAttach];
        if (edgeMapNewToOld.find(newedgeIdToAttach) == edgeMapNewToOld.end()){
            newedgeIdToAttach = bidirEdgeMap[newedgeIdToAttach];
        }
        bidirNewedgeIdToAttach = bidirEdgeMap[newedgeIdToAttach];
        
        edgeMapNewToOld[bidirNewedgeIdToAttach]=std::vector<int>();
        // std::cout << "moving edges: ";
        for (auto &e: edges){
            if (visistedEdges[e.edge_id] && edgeMapOldToNew[e.edge_id] == newedgeIdToAttach) {
                // std::cout << e.edge_id << " ";
                edgeMapNewToOld[bidirNewedgeIdToAttach].push_back(e.edge_id);
                // Remove e.edge_id from edgeMapNewToOld[newedgeIdToAttach] vector
                auto& vec = edgeMapNewToOld[newedgeIdToAttach];
                vec.erase(std::remove(vec.begin(), vec.end(), e.edge_id), vec.end());
                edgeMapOldToNew[e.edge_id]=bidirNewedgeIdToAttach;
            }
        }
        // std::cout << std::endl;

        /* update bidirEdgeMap */
        bidirEdgeMap[newedgeIdToAttach]=idx;
        bidirEdgeMap[idx]=newedgeIdToAttach;
        bidirEdgeMap[bidirNewedgeIdToAttach]=idx+1;
        bidirEdgeMap[idx+1]=bidirNewedgeIdToAttach;

        /* add new tip information */
        edgeMapNewToOld[idx+2]=std::vector<int>();
        for (auto &e: edges){
            if (!visistedEdges[e.edge_id]) {
                edgeMapNewToOld[idx+2].push_back(e.edge_id);
                edgeMapOldToNew[e.edge_id]=idx+2;
                visistedEdges[e.edge_id]=true;
            }
        }
        bidirEdgeMap[idx+2]=idx+3;
        bidirEdgeMap[idx+3]=idx+2;

        // std::cout << "bidirEdgeMap: " << std::endl;
        // for(auto &e: bidirEdgeMap){
        //     std::cout << "edge " << e.first << " -> " << e.second << std::endl;
        // }

        // std::cout << "edgeMapNewToOld: " << std::endl;
        // for(auto &e: edgeMapNewToOld){
        //     std::cout << "edge " << e.first << " -> ";
        //     for(auto &f: e.second){
        //         std::cout << f << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // std::cout << "edgeMapOldToNew: " << std::endl;
        // for (int i=0;i<edgeMapOldToNew.size();i++){
        //     std::cout << "edge " << i << " -> " << edgeMapOldToNew[i] << std::endl;
        // }
        
        if(params.in == "r"){
            mashDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist
            );
        } else if(params.in == "m"){
            msaDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist
            );
        }
        cudaDeviceSynchronize();

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
        
        auto iter = thrust::find_if(minPos.begin(), minPos.end(), find_if_tuple_bidir(newedgeIdToAttach));
        if (iter == minPos.end()) {
            iter = thrust::find_if(minPos.begin(), minPos.end(), find_if_tuple_bidir(bidirNewedgeIdToAttach));
        }
        if (iter == minPos.end()) {
            std::cerr << "Error: Edge not found in minPos" << std::endl;
            exit(1);
        }
        
        thrust::tuple<int,double,double> smallest=*iter;
        
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid=newedgeIdToAttach;
        double fracLen=thrust::get<1>(smallest),addLen=thrust::get<2>(smallest);
        // std::cout << "eid: " << eid << " fracLen: " << fracLen << " addLen: " << addLen << std::endl;
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
        cudaDeviceSynchronize();
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
    }
}
