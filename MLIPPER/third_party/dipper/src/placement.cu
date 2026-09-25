/* placement.cu: Exact placement mode. Sequential GPU placement: distances -> lim updates ->
 * best-branch choice -> tree structure update (adjacency, depth, BFS/DFS, levelst/leveled). */

#include "mash_placement.cuh"

#include <stdio.h>
#include <queue>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cub/cub.cuh>

/** Check CUDA errors and sync; exit on failure. */
void checkCudaErrorsHere(const char* location) {
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA error before %s: %s\n", location, cudaGetErrorString(err));
        exit(-1);
    }
    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        printf("CUDA error after sync at %s: %s\n", location, cudaGetErrorString(err));
        exit(-1);
    }
}

/* Allocate GPU arrays for exact placement: adjacency, BFS/DFS/depth/levels, dist, len, lim. */
void MashPlacement::PlacementDeviceArrays::allocateDeviceArrays(size_t num){
    cudaError_t err;
    numSequences = int(num);
    bd = 2, idx = 0;
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
    err = cudaMalloc(&d_rev, numSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_bfsorder, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_dfsorder, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_dep, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_dfsrk, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_levelst, numSequences*sizeof(int)*2);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_leveled, numSequences*sizeof(int)*2);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_dist, numSequences*sizeof(double));
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
    err = cudaMalloc(&d_lim, numSequences*8*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
}

/** Reset adjacency slots and node metadata (head, dep, dfsrk, levelst, leveled) to sentinels. */
__global__ void initialize(
    int lim,
    int nodes,
    int * head,
    int * nxt,
    int * belong,
    int * e,
    int * dep,
    int * dfsrk,
    int * levelst,
    int * leveled
){
    int tx=threadIdx.x,bs=blockDim.x;
    int bx=blockIdx.x,gs=gridDim.x;
    int idx=tx+bs*bx;
    if(idx<lim){
        nxt[idx] = -1;
        e[idx] = -1;
        belong[idx] = -1;
    }
    if(idx<nodes) head[idx] = -1, dep[idx]=nodes*10, dfsrk[idx]=levelst[idx]=leveled[idx]=-1;
}

struct compare_tuple {
  __host__ __device__
  bool operator()(thrust::tuple<int,double,double> lhs, thrust::tuple<int,double,double> rhs)
  {
    return thrust::get<2>(lhs) < thrust::get<2>(rhs);
    //Always find the tuple whose third value (the criteria we want to minimize) is minimized
  }
};
/*
Three variables in tuple:
ID of branch in linked list,
distance to new node inserted on branch from starting vertex (belong[id]),
distance from new node inserted on branch to new node inserted outside branch
*/

/** For each candidate edge, compute fracLen and addLen from lim; output (eid, fracLen, addLen). */
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
    int totSeqNum,
    double * d_lim,
    int * d_dep
){
    int tx=threadIdx.x,bs=blockDim.x,bx=blockIdx.x,gs=gridDim.x;
    int idx=tx+bs*bx;
    if(idx>=lim) return;
    if(idx>=num*4-4||d_dep[belong[idx]]>d_dep[e[idx]]){
        thrust::tuple <int,double,double> minTuple(0,0,2);
        minPos[bx*bs+tx]=minTuple;
        return;
    }
    int x=belong[idx],oth=e[idx];
    int eid=idx,otheid;
    double dis1=0, dis2=0, val;
    dis1 = d_lim[eid];
    otheid=head[oth];
    while(e[otheid]!=x) assert(otheid!=-1),otheid=nxt[otheid];
    dis2 = d_lim[otheid];
    double additional_dis=(dis1+dis2-len[eid])/2;
    // printf("%d %d %d %lf %lf %lf\n",eid, x,oth,dis1,dis2,additional_dis);
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

/** Insert new leaf on edge eid: split edge, add middle node, attach outside (placeId). Update rev. */
__global__ void updateTreeStructure(
    int * head,
    int * nxt,
    int * e,
    double * len,
    int * belong,
    int * rev,
    int eid,
    double fracLen,
    double addLen,
    int placeId, // Id of the newly placed node
    int edgeCount, // Position to insert a new edge in linked list
    int numSequences,
    int * dfsrk,
    int * dep
){
    int middle=placeId+numSequences-1, outside=placeId;
    int x=belong[eid],y=e[eid];
    double originalDis=len[eid];
    int xe,ye;
    for(int i=head[x];i!=-1;i=nxt[i])
        if(e[i]==y){
            e[i]=middle,len[i]=fracLen,xe=i,rev[xe]=edgeCount;
            break;
        }
    for(int i=head[y];i!=-1;i=nxt[i])
        if(e[i]==x){
            e[i]=middle,len[i]-=fracLen,ye=i,rev[ye]=edgeCount+1;
            break;
        }
    /*
    Need to update:
    e, len, nxt, head, belong, closest_dis, closest_id
    */
    //middle -> x
    e[edgeCount]=x,len[edgeCount]=fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle,rev[edgeCount]=xe;
    edgeCount++;
    //middle -> y
    e[edgeCount]=y,len[edgeCount]=originalDis-fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle,rev[edgeCount]=ye;
    edgeCount++;
    //outside -> middle
    e[edgeCount]=middle,len[edgeCount]=addLen,nxt[edgeCount]=head[outside],head[outside]=edgeCount,belong[edgeCount]=outside,rev[edgeCount]=edgeCount+1;
    edgeCount++;
    //middle -> outside
    e[edgeCount]=outside,len[edgeCount]=addLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle,rev[edgeCount]=edgeCount-1;
    edgeCount++;
    if(dfsrk[x]>dfsrk[y]){
        int temp=x;
        y=x, x=temp;
    }
    dfsrk[middle]=dfsrk[y];
    dfsrk[outside]=dfsrk[middle]+1;
    dep[middle]=dep[x],     dep[outside]=dep[middle]+1;
}

/** Build 3-node tree: nodes 0, 1, and nv=numSequences; two edges of length dis[0]/2. */
__global__ void buildInitialTree(
    int numSequences,
    int * head,
    int * e,
    double * len,
    int * nxt,
    int * belong,
    double * dis,
    int edgeCount,
    int * bfsorder,
    int * dfsorder,
    int * dep,
    int * dfsrk,
    int * levelst,
    int * leveled,
    int * rev
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
    // nv -> 1
    e[edgeCount]=1,len[edgeCount]=d/2,nxt[edgeCount]=head[nv],head[nv]=edgeCount,belong[edgeCount]=nv;
    edgeCount++;
    bfsorder[0]=nv,bfsorder[1]=0,bfsorder[2]=1;
    dfsorder[0]=nv,dfsorder[1]=0,dfsorder[2]=1;
    dep[nv]=0, dep[0]=dep[1]=1;
    dfsrk[nv]=0,dfsrk[0]=1,dfsrk[1]=2;
    levelst[0]=leveled[0]=0, levelst[1]=1,leveled[1]=2;
    rev[0]=2, rev[2]=0, rev[1]=3,rev[3]=1;
    dep[0]=1, dep[1]=1, dep[nv]=0;
    dfsrk[0]=1, dfsrk[1]=2, dfsrk[nv]=0;
}


__global__ void updateFromBottomToTop(
    int tot,
    int level,
    int * head,
    int * e,
    int * nxt,
    double * len,
    int * levelst,
    int * leveled,
    int * bfsorder,
    double * dist,
    double * lim,
    int * rev,
    int * dep
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int idx=tx+bs*bx;
    idx=levelst[level]+idx;
    if(idx>leveled[level]) return;
    idx=bfsorder[idx];
    double mx=0;
    if(idx<tot) mx=dist[idx];
    for(int i=head[idx];i!=-1;i=nxt[i]){
        if(dep[e[i]]>dep[idx]){
            double req=lim[rev[i]]-len[i];
            if(req>mx) mx=req;
        }
    }
    for(int i=head[idx];i!=-1;i=nxt[i])
        if(dep[e[i]]<dep[idx])
            lim[i]=mx;
}

/** Propagate lim from root toward children: sibling max of (lim - len). */
__global__ void updateFromTopToBottom(
    int tot,
    int level,
    int * head,
    int * e,
    int * nxt,
    double * len,
    int * levelst,
    int * leveled,
    int * bfsorder,
    double * dist,
    double * lim,
    int * rev,
    int * dep
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int idx=tx+bs*bx;
    idx=levelst[level]+idx;
    if(idx>leveled[level]) return;
    idx=bfsorder[idx];
    for(int i=head[idx];i!=-1;i=nxt[i]){
        if(dep[e[i]]>dep[idx]){
            double mx=0;
            for(int j=head[idx];j!=-1;j=nxt[j]){
                if(e[j]!=e[i]){
                    double req=lim[rev[j]]-len[j];
                    if(req>mx) mx=req;
                }
            }
            lim[i]=mx;
            // printf("%d %lf\n",i,lim[i]);
        }
    }
}

/** Shift DFS ranks >= dfsrk[ex1] by 2 for nodes in [id, num), excluding ex1/ex2. */
__global__ void updateDfsRk(
    int tot,
    int * dfsrk,
    int ex1,
    int ex2,
    int id,
    int num
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int idx=tx+bs*bx;
    if(idx>id&&idx<num) return;
    if(idx>=tot||idx==ex1||idx==ex2) return;
    if(dfsrk[idx]>=dfsrk[ex1]) dfsrk[idx]+=2;
}

/** Write temp[idx] = dfsrk[idx]-1 for valid insertion-range nodes, else sentinel. */
__global__ void findEndRk(
    int tot,
    int * dfsrk,
    int * dep,
    int * temp,
    int ref,
    int id,
    int num
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int idx=tx+bs*bx;
    if(idx>=tot) return;
    // printf("%d %d %d %d\n",idx,dfsrk[idx],dep[idx],dep[ref]);
    if(idx>id&&idx<num) temp[idx]=1000000000;
    else if(dfsrk[idx]<=dfsrk[ref]+2||dep[idx]>dep[ref]+1) temp[idx]=1000000000;
    else temp[idx]=dfsrk[idx]-1;
}

__global__ void updateDepth(
    int tot,
    int * dep,
    int * dfsrk,
    int ref,
    int edrk,
    int * bfsorder,
    int id,
    int num
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int idx=tx+bs*bx;
    if(idx>id&&idx<num) return;
    if(idx>=tot) return;
    bfsorder[idx]=idx;
    if(dfsrk[idx]<=edrk&&dfsrk[idx]>=dfsrk[ref]) dep[idx]++;
}

/** Update levelst[level] and leveled[level] from BFS order and depth. */
__global__ void updateLevelStEd(
    int tot,
    int * bfsorder,
    int * dep,
    int * levelst,
    int * leveled,
    int id,
    int num
){
    int tx=threadIdx.x, bs=blockDim.x, bx=blockIdx.x;
    int idx=tx+bs*bx;
    if(idx>=id*2+1) return;
    // printf("%d %d %d\n",idx,bfsorder[idx],dep[bfsorder[idx]]);
    if(idx==0||dep[bfsorder[idx-1]]!=dep[bfsorder[idx]]) levelst[dep[bfsorder[idx]]]=idx;
    if(idx+1==id*2+1||dep[bfsorder[idx+1]]!=dep[bfsorder[idx]]) leveled[dep[bfsorder[idx]]]=idx;
}

void MashPlacement::PlacementDeviceArrays::deallocateDeviceArrays(){
    cudaFree(d_head);
    cudaFree(d_e);
    cudaFree(d_nxt);
    cudaFree(d_belong);
    cudaFree(d_bfsorder);
    cudaFree(d_dfsorder);
    cudaFree(d_dep);
    cudaFree(d_dist);
    cudaFree(d_len);
    cudaFree(d_lim);
    cudaFree(d_dfsrk);
    cudaFree(d_levelst);
    cudaFree(d_leveled);
}

/** Copy tree arrays to host and emit Newick via recursive DFS over adjacency. */
void MashPlacement::PlacementDeviceArrays::printTree(std::vector <std::string> name, std::ofstream& output_){
    int * h_head = new int[numSequences*2];
    int * h_e = new int[numSequences*8];
    int * h_nxt = new int[numSequences*8];
    double * h_len = new double[numSequences*8];
    auto err = cudaMemcpy(h_head, d_head, numSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
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
}

/** Exact placement loop: init, build initial tree, then for each taxon compute distances,
 * update lim (bottom-up then top-down), choose best branch, update structure and DFS/depth/levels. */
void MashPlacement::PlacementDeviceArrays::findPlacementTree(
    Param& params,
    const MashDeviceArrays& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArrays& msaDeviceArrays
){ 
    if(params.in == "d"){
        matrixReader.distConstructionOnGpu(params, 0, d_dist);
    }

    int threadNum = 1024, blockNum = 1024;
    initialize <<<blockNum, threadNum>>> (
        numSequences*4-4,
        numSequences*2-1,
        d_head,
        d_nxt,
        d_belong,
        d_e,
        d_dep,
        d_dfsrk,
        d_levelst,
        d_leveled
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
    // double * h_dis = new double[numSequences];
    // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
    // fprintf(stderr, "%d\n",1);
    // for(int j=0;j<1;j++) fprintf(stderr,"%.8lf ",h_dis[j]);std::cerr<<'\n';
    buildInitialTree <<<1,1>>> (
        numSequences,
        d_head,
        d_e,
        d_len,
        d_nxt,
        d_belong,
        d_dist,
        idx,
        d_bfsorder,
        d_dfsorder,
        d_dep,
        d_dfsrk,
        d_levelst,
        d_leveled,
        d_rev
    );
    idx += 4;
    // cudaDeviceSynchronize();
    // std::cerr<<"FFFF\n";
    thrust::device_vector <thrust::tuple<int,double,double>> minPos(numSequences*4-4);
    std::chrono::nanoseconds disTime(0), treeTime(0);
    int *id_maxdep=new int[1], *maxdep=new int[1];
    int *levelst=new int[numSequences], *leveled=new int[numSequences];
    int * d_temp;
    auto err = cudaMalloc(&d_temp, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    for(int i=bd;i<numSequences;i++){
        auto disStart = std::chrono::high_resolution_clock::now();

        if(params.in == "r"){
            mashDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist
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
            msaDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist
            );
        }
        // cudaDeviceSynchronize();
        // double * h_dis = new double[numSequences];
        // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
        // fprintf(stderr, "%d\n",i);
        // for(int j=0;j<i;j++) std::cerr<<h_dis[j]<<" ";std::cerr<<'\n';


        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();

        cudaMemcpy(id_maxdep, d_bfsorder+i*2-2, sizeof(int), cudaMemcpyDeviceToHost);
        int id=*id_maxdep;
        cudaMemcpy(maxdep, d_dep+id,sizeof(int),cudaMemcpyDeviceToHost);
        int mx=*maxdep;
        // std::cerr<<"-------"<<mx<<'\n';
        cudaMemcpy(levelst, d_levelst, sizeof(int)*(mx+1), cudaMemcpyDeviceToHost);
        cudaMemcpy(leveled, d_leveled, sizeof(int)*(mx+1), cudaMemcpyDeviceToHost);
        for(int j=mx;j>=0;j--){
            // std::cerr<<j<<" "<<levelst[j]<<" "<<leveled[j]<<'\n';
            blockNum = (leveled[j]-levelst[j]+1+255)/256;
            updateFromBottomToTop <<<blockNum, threadNum>>> (
                numSequences,
                j,
                d_head,
                d_e,
                d_nxt,
                d_len,
                d_levelst,
                d_leveled,
                d_bfsorder,
                d_dist,
                d_lim,
                d_rev,
                d_dep
            );
        }
        for(int j=0;j<=mx;j++){
            blockNum = (leveled[j]-levelst[j]+1+255)/256;
            updateFromTopToBottom <<<blockNum, threadNum>>> (
                numSequences,
                j,
                d_head,
                d_e,
                d_nxt,
                d_len,
                d_levelst,
                d_leveled,
                d_bfsorder,
                d_dist,
                d_lim,
                d_rev,
                d_dep
            );
        }
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
            numSequences,
            d_lim,
            d_dep
        );
        auto iter=thrust::min_element(minPos.begin(),minPos.end(),compare_tuple());
        thrust::tuple<int,double,double> smallest=*iter;

        int eid=thrust::get<0>(smallest);
        double fracLen=thrust::get<1>(smallest),addLen=thrust::get<2>(smallest);
        // std::cerr<<eid<<" "<<fracLen<<" "<<addLen<<'\n';
        updateTreeStructure <<<1,1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_belong,
            d_rev,
            eid,
            fracLen,
            addLen,
            i,
            idx,
            numSequences,
            d_dfsrk,
            d_dep
        );
        idx+=4;

        /* Update DFS rank (exclude new nodes), then depth, BFS order, levelst/leveled. */
        // cudaMemcpy(levelst, d_dfsrk, (numSequences+i)*sizeof(int), cudaMemcpyDeviceToHost);
        // for(int j=0;j<numSequences+i;j++) std::cerr<<levelst[j]<<' ';std::cerr<<'\n';
        blockNum = (numSequences+i+255)/256;
        updateDfsRk <<<blockNum,threadNum>>>(
            numSequences+i,
            d_dfsrk,
            numSequences+i-1,
            i,
            i,
            numSequences
        );
        // cudaMemcpy(levelst, d_dfsrk, (numSequences+i)*sizeof(int), cudaMemcpyDeviceToHost);
        // for(int j=0;j<numSequences+i;j++) std::cerr<<levelst[j]<<' ';std::cerr<<'\n';
        // cudaMemcpy(levelst, d_dep, (numSequences+i)*sizeof(int), cudaMemcpyDeviceToHost);
        // for(int j=0;j<numSequences+i;j++) std::cerr<<levelst[j]<<' ';std::cerr<<'\n';
        /*
        Update depth based on DFS order, new nodes included
        */
        findEndRk <<<blockNum, threadNum>>>(
            numSequences+i,
            d_dfsrk,
            d_dep,
            d_temp,
            numSequences+i-1,
            i,
            numSequences
        );
        // cudaMemcpy(levelst, d_temp, (numSequences+i)*sizeof(int), cudaMemcpyDeviceToHost);
        // for(int j=0;j<numSequences+i;j++) std::cerr<<levelst[j]<<' ';std::cerr<<'\n';
        int small = thrust::reduce(thrust::device, d_temp,d_temp+numSequences+i,numSequences+i-1, thrust::minimum<int>());
        // std::cerr<<small<<'\n';
        // cudaMemcpy(levelst, d_dep, (numSequences+i)*sizeof(int), cudaMemcpyDeviceToHost);
        // for(int j=0;j<numSequences+i;j++) std::cerr<<levelst[j]<<' ';std::cerr<<'\n';
        /* Update depth from DFS range, then BFS order and levelst/leveled. */
        updateDepth <<<blockNum, threadNum>>>(
            numSequences+i,
            d_dep,
            d_dfsrk,
            numSequences+i-1,
            small,
            d_bfsorder,
            i,
            numSequences
        );
        // cudaMemcpy(levelst, d_dep, (numSequences+i)*sizeof(int), cudaMemcpyDeviceToHost);
        // for(int j=0;j<numSequences+i;j++) std::cerr<<levelst[j]<<' ';std::cerr<<'\n';
        thrust::copy(thrust::device, d_dep, d_dep+numSequences+i, d_temp);
        thrust::stable_sort_by_key(thrust::device, d_temp, d_temp+numSequences+i, d_bfsorder);
        /*
        Update level st/ed
        */
        // std::cerr<<"#########\n";
        // cudaMemcpy(levelst, d_bfsorder, (numSequences+i)*sizeof(int), cudaMemcpyDeviceToHost);
        // for(int j=0;j<numSequences+i;j++) std::cerr<<levelst[j]<<' ';std::cerr<<'\n';
        updateLevelStEd <<<blockNum, threadNum>>>(
            numSequences+i,
            d_bfsorder,
            d_dep,
            d_levelst,
            d_leveled,
            i,
            numSequences
        );
        // cudaDeviceSynchronize();
        auto treeEnd = std::chrono::high_resolution_clock::now();
        disTime += disEnd - disStart;
        treeTime += treeEnd - treeStart;
    }
    std::cerr << "Distance Operation Time " <<  disTime.count()/1000000 << " ms\n";
    std::cerr << "Tree Operation Time " <<  treeTime.count()/1000000 << " ms\n";
}
