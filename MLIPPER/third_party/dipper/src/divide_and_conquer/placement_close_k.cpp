#include "../mash_placement.cuh"

#include <stdio.h>
#include <queue>
#include <chrono>
#include <iostream>
#include <sstream>
#include <tuple>
#include <cassert>
#include <algorithm>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>

void MashPlacement::KPlacementDeviceArraysHostDC::allocateHostArraysDC(size_t num, size_t totalNum){
    numSequences = int(num);
    totalNumSequences = int(totalNum);
    bd = 2, idx = 0;

    h_dist = (double *)malloc(totalNumSequences*sizeof(double));
    h_head = (int *)malloc(totalNumSequences*2*sizeof(int));    
    h_e = (int *)malloc(totalNumSequences*8*sizeof(int));
    h_len = (double *)malloc(totalNumSequences*8*sizeof(double));
    h_nxt = (int *)malloc(totalNumSequences*8*sizeof(int));
    h_belong = (int *)malloc(totalNumSequences*8*sizeof(int));
    h_closest_dis = (double *)malloc(totalNumSequences*40*sizeof(double));
    h_closest_id = (int *)malloc(totalNumSequences*40*sizeof(int));
    h_closest_dis_cluster = (double *)malloc(totalNumSequences*40*sizeof(double));
    h_closest_id_cluster = (int *)malloc(totalNumSequences*40*sizeof(int));

    return;
}

// void initialize(
//     int lim,
//     int nodes,
//     double * d_closest_dis,
//     int * d_closest_id,
//     int * head,
//     int * nxt,
//     int * belong,
//     int * e
// ){
//     for (int idx=0;idx<lim;idx++){
//         for(int i=0;i<5;i++){
//             d_closest_dis[idx*5+i]=2;
//             d_closest_id[idx*5+i]=-1;
//         }
//         nxt[idx] = -1;
//         e[idx] = -1;
//         belong[idx] = -1;
//     }
    
//     for (int idx=0;idx<nodes;idx++){
//         head[idx] = -1;
//     }
//     return;
// }

/*
Three variables in tuple:
ID of branch in linked list,
distance to new node inserted on branch from starting vertex (belong[id]),
distance from new node inserted on branch to new node inserted outside branch
*/
struct compare_tuple {
  bool operator()(std::tuple<int,double,double> lhs, std::tuple<int,double,double> rhs)
  {
    return std::get<2>(lhs) < std::get<2>(rhs);
    //Always find the tuple whose third value (the criteria we want to minimize) is minimized
  }
};


void calculateBranchLengthCpu(
    int num, // should be bd, not numSequences 
    int * head,
    int * nxt,
    double * dis, 
    int * e, 
    double * len, 
    int * belong,
    std::tuple<int,double,double> * minPos,
    int lim,
    double * closest_dis,
    int * closest_id
){
    for (int idx=0; idx<lim; idx++) {
        if(idx>=num*4-4||belong[idx]<e[idx]){
            std::tuple <int,double,double> minTuple(0,0,2);
            minPos[idx]=minTuple;
            continue;
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

        double rest=len[eid]-dis1-dis2;
        dis1+=rest/2,dis2+=rest/2;
        std::tuple <int,double,double> minTuple(eid,dis1,additional_dis);
        minPos[idx]=minTuple;
    }
    return;

}


void calculateBranchLengthSpecialIDCpu(
    int num, // useless here
    int * head,
    int * nxt,
    // double * dis, 
    std::vector<double> & dis,
    int * e, 
    double * len, 
    int * belong,
    std::tuple<int,double,double> * minPos,
    int lim,
    double * closest_dis,
    int * closest_id,
    int numToCalculate,
    std::vector<int> & h_edgeMask
    // int * h_edgeMask 
){
    // printf("numToCalculate: %d\n", numToCalculate);                                      
    for (int idx_=0; idx_<numToCalculate; idx_++) {
    // tbb::parallel_for(tbb::blocked_range<int>(0, numToCalculate), [&](tbb::blocked_range<int> range){
    // for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
        int idx = h_edgeMask[idx_];
        if(belong[idx]<e[idx]){
            std::tuple <int,double,double> minTuple(0,0,2);
            minPos[idx_]=minTuple;
            continue;
        }
        // printf("idx: %d, belong: %d, e: %d\n", idx, belong[idx], e[idx]);
        int x=belong[idx],oth=e[idx];
        int eid=idx,otheid;
        double dis1=0, dis2=0, val;
        for(int i=0;i<5;i++)
            if(closest_id[eid*5+i]!=-1){
                val = dis[closest_id[eid*5+i]]-closest_dis[eid*5+i];
                if(val>dis1) dis1=val;
            }
        otheid=head[oth];
        while(e[otheid]!=x) {
            // printf("othL %d, otheid: %d, e[otheid]: %d, x: %d\n", oth, otheid, e[otheid], x);
            assert(otheid!=-1);
            otheid=nxt[otheid];
        }
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
        std::tuple <int,double,double> minTuple(eid,dis1,additional_dis);
        minPos[idx_]=minTuple;
        // printf("idx: %d, 1: %d, 2: %lf, 3: %lf\n", idx_, eid, dis1, additional_dis);
    }
    // }});
    
}


void updateClosestNodesCpu(
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
    int * h_edgeMaskIndex
){
    int l=0,r=-1;
    id[++r]=x,dis[r]=0,from[r]=-1;
    while(l<=r){
        int node=id[l],fb=from[l];
        double d=dis[l];
        l++;
        // 
        for(int i=head[node];i!=-1;i=nxt[i]){
            if(h_edgeMaskIndex[i]==-1) continue;
            // std::cerr << "updating " << node << " " << head[node]<< " " << e[head[node]] << " " << h_edgeMask[head[node]] << std::endl;
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

void updateClosestNodesInClusterCpu(
    int * head,
    int * nxt,
    int * e,
    double * len,
    double * closest_dis,
    int * closest_id,
    int x,
    std::vector <int> & id,
    std::vector <int> & from,
    std::vector<double> & dis,
    // int * id,
    // int * from,
    // double * dis,
    int cluster_eid,
    int * belong
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


void updateTreeStructureCpu(
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

void updateTreeStructureInClusterCpu(
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
    int placeId, // Id of the newly placed node (node id = leave_count*4-4)
    int edgeCount, // Position to insert a new edge in linked list (edge index)
    int totalNumSequences,
    int placeCount // this is the placeCount-th leave (nth-leave)
){
    int middle=placeCount+totalNumSequences-1, outside=placeId;
    int x=belong[eid],y=e[eid];
    double originalDis=len[eid];
    int xe,ye;
    for(int i=head[x];i!=-1;i=nxt[i]) {
        if(e[i]==y){
            e[i]=middle,len[i]=fracLen,xe=i;
            break;
        }
    }
    for(int i=head[y];i!=-1;i=nxt[i]){
        if(e[i]==x){
            e[i]=middle,len[i]-=fracLen,ye=i;
            break;
        }
    }

    // printf("updateTreeStructureInClusterCpu: eid: %d, x: %d, y: %d, xe: %d, ye: %d\n", eid, x, y, xe, ye);
    // for(int i=0;i<5;i++){
    //     printf("closest_id[%d]: %d, closest_dis[%d]: %lf\n", ye*5+i, closest_id[ye*5+i], ye*5+i, closest_dis[ye*5+i]);
    // }
    // for(int i=0;i<5;i++){
    //     printf("closest_id[%d]: %d, closest_dis[%d]: %lf\n", xe*5+i, closest_id[xe*5+i], xe*5+i, closest_dis[xe*5+i]);
    // }

    /*
    Need to update:
    e, len, nxt, head, belong, closest_dis, closest_id
    */
    //middle -> x
    e[edgeCount]=x,len[edgeCount]=fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    for(int i=0;i<5;i++){
        if(closest_id[ye*5+i]!=-1){
            closest_id[edgeCount*5+i]=closest_id[ye*5+i];
            closest_dis[edgeCount*5+i]=closest_dis[ye*5+i]+originalDis-fracLen;
        }
    }
        
    edgeCount++;
    //middle -> y
    e[edgeCount]=y,len[edgeCount]=originalDis-fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    for(int i=0;i<5;i++) {
        if(closest_id[xe*5+i]!=-1){
            closest_id[edgeCount*5+i]=closest_id[xe*5+i];
            closest_dis[edgeCount*5+i]=closest_dis[xe*5+i]+fracLen;
        }
    }
        
    edgeCount++;
    //outside -> middle
    e[edgeCount]=middle,len[edgeCount]=addLen,nxt[edgeCount]=head[outside],head[outside]=edgeCount,belong[edgeCount]=outside;
    edgeCount++;
    //middle -> outside
    e[edgeCount]=outside,len[edgeCount]=addLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    
    // printf("edgecount and e1 e2: %d, %d, %d\n", edgeCount, edgeCount-2, edgeCount-3);
    // for(int i=0;i<5;i++){
    //     printf("closest_id[%d]: %d, closest_dis[%d]: %lf\n", edgeCount*5+i, closest_id[edgeCount*5+i], edgeCount*5+i, closest_dis[edgeCount*5+i]);
    // }
    int e1=edgeCount-2, e2=edgeCount-3;
    for(int i=0;i<5;i++){
        if(closest_id[e1*5+i]==-1) break;
        for(int j=0;j<5;j++) {
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
            
    }
    // for(int i=0;i<5;i++){
    //     printf("closest_id[%d]: %d, closest_dis[%d]: %lf\n", edgeCount*5+i, closest_id[edgeCount*5+i], edgeCount*5+i, closest_dis[edgeCount*5+i]);
    // }
    for(int i=0;i<5;i++){
        if(closest_id[e2*5+i]==-1) break;
        for(int j=0;j<5;j++){
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
            
    }
    // for(int i=0;i<5;i++){
    //     printf("closest_id[%d]: %d, closest_dis[%d]: %lf\n", edgeCount*5+i, closest_id[edgeCount*5+i], edgeCount*5+i, closest_dis[edgeCount*5+i]);
    // }
    edgeCount++;

}

void buildInitialTreeCpu(
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
    // nv -> 0
    e[edgeCount]=1,len[edgeCount]=d/2,nxt[edgeCount]=head[nv],head[nv]=edgeCount,belong[edgeCount]=nv;
    edgeCount++;
}

void updateClusterInfoCpu (
    int leafID,
    int edgeidx,
    std::vector <int> & h_leafMask,
    std::vector <int> & h_edgeMask,
    // int * h_leafMask,
    // int * h_edgeMask,
    int edgeCount,
    int leafCount
    // int * h_edgeMaskIndex
){
    h_leafMask[leafCount++]=leafID;
    for(int i=1;i<=4;i++) {
        h_edgeMask[edgeCount++]=edgeidx-i;
        // h_edgeMaskIndex[edgeidx-i]=edgeidx-i;
    }
}

void copyClosestIdsCpu(
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

void copyBackClosestIdsCpu(
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

void initializeClusterCpu (
    int eid,
    int * e,
    int * belong,
    int * head,
    int * nxt,
    int * closest_id,
    std::vector <int> & edgeMask,
    std::vector <int> & leafMask
    // int * edgeMask,
    // int * leafMask
){
    int x=belong[eid],y=e[eid];
    int otheid=head[y];
    while(e[otheid]!=x) assert(e[otheid]!=-1),otheid=nxt[otheid];
    
    int leafCount=0;
    for(int i=0;i<5;i++) leafMask[leafCount++]=closest_id[eid*5+i];
    for(int i=0;i<5;i++) leafMask[leafCount++]=closest_id[otheid*5+i];

    int edgeCount=0;
    edgeMask[edgeCount++]=eid;
    edgeMask[edgeCount++]=otheid;
}



void MashPlacement::KPlacementDeviceArraysHostDC::deallocateHostArraysDC(){
    
    free(h_head);
    free(h_e);
    free(h_nxt);
    free(h_belong);
    free(h_closest_id);
    free(h_closest_id_cluster);
    free(h_closest_dis);
    free(h_closest_dis_cluster);
    free(h_dist);
    free(h_len);

    h_head = nullptr;
    h_e = nullptr;
    h_nxt = nullptr;
    h_belong = nullptr;
    h_closest_id = nullptr;
    h_closest_id_cluster = nullptr;
    h_closest_dis = nullptr;
    h_closest_dis_cluster = nullptr;
    h_dist = nullptr;
    h_len = nullptr;

}


void MashPlacement::KPlacementDeviceArraysHostDC::printTreeCpuDC(std::vector <std::string> name){
    auto print=[&](int node, int from, auto&& print)->void {
        if(h_nxt[h_head[node]]!=-1){
            printf("(");
            std::vector <int> pos;
            for(int i=h_head[node];i!=-1;i=h_nxt[i])
                if(h_e[i]!=from)
                    pos.push_back(i);
            for(size_t i=0;i<pos.size();i++){
                print(h_e[pos[i]],node, print);
                printf(":");
                printf("%.5g%c",h_len[pos[i]],i+1==pos.size()?')':',');
            }
        }
        else std::cout<<name[node];
    };
    
    print(totalNumSequences+bd-2,-1, print);
    std::cout<<";";
}

/* Clusterting function on CPU - > might need modification
void MashPlacement::KPlacementDeviceArraysHost::findCluster(
    Param& params,
    const MashDeviceArrays& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArrays& msaDeviceArrays
)
 auto clusterStart = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds disTime(0), treeTime(0);
    
    tbb::parallel_for(tbb::blocked_range<int>(numSequences, totalNumSequences), [&](tbb::blocked_range<int> range){
    for (int i = range.begin(); i < range.end(); ++i) {
        std::vector <std::tuple<int,double,double>> minPos(totalNumSequences*4-4);
        double * h_dist_local = (double *)malloc(totalNumSequences*sizeof(double));
    // for(int i=numSequences;i<totalNumSequences;i++){
        auto disStart = std::chrono::high_resolution_clock::now();
        if(params.in == "r"){
            mashDeviceArrays.distRangeConstructionOnCpu(
                params,
                i,
                h_dist_local,
                0,
                numSequences-1
            );
        }
        // else if(params.in == "d"){
        //     matrixReader.distConstructionOnGpu(
        //         params,
        //         i,
        //         h_dist
        //     );
        // }
        // else if(params.in == "m"){
        //     msaDeviceArrays.distConstructionOnGpu(
        //         params,
        //         i,
        //         h_dist
        //     );
        // }

        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        calculateBranchLengthCpu(
            i,
            h_head,
            h_nxt,
            h_dist_local,
            h_e,
            h_len,
            h_belong,
            minPos.data(),
            numSequences*4-4,
            h_closest_dis,
            h_closest_id
        );
        auto iter=std::min_element(minPos.begin(),minPos.begin()+numSequences*4-4,compare_tuple());
        std::tuple<int,double,double> smallest=*iter;
        cluster_id[i] = std::get<0>(smallest);
        // std::cerr << "Cluster ID of " << i << " is " << cluster_id[i] << "\n";
    }});
    // }

    auto clusterEnd = std::chrono::high_resolution_clock::now();
    auto clusterTime = clusterEnd - clusterStart;


    std::cerr << "Finished clustering in: "<< clusterTime.count()/1000000 << " ms\n";
} 
*/

void MashPlacement::KPlacementDeviceArraysHostDC::findClusterTreeDC(
    Param& params,
    const MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArraysDC& msaDeviceArrays
){ 
    // tbb::global_control c(tbb::global_control::max_allowed_parallelism, 32);
    int idx=params.backboneSize*4-4;
    
    int *cluster_id = clusterID;
    
    auto inClusterStart = std::chrono::high_resolution_clock::now();
    
    std::vector<std::vector <int>> contains(numSequences*4-4);

    for(int i=numSequences;i<totalNumSequences;i++) contains[cluster_id[i]].push_back(i);

    int insertLeafCount=numSequences;
    
    std::vector <int> containsPrefixSum(numSequences*4-4);
    tbb::parallel_scan( tbb::blocked_range<size_t>(0, contains.size()),0,
        [&](const tbb::blocked_range<size_t>& r, int sum, bool is_final_scan) -> int {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                int sz = contains[i].size();
                if (is_final_scan) {
                    containsPrefixSum[i] = sum + sz;
                }
                sum += sz;
            }
            return sum;
        },
        [](int left, int right) -> int {
            return left + right;
        }
    );
    
    // for(int i=0;i<numSequences*4-4;i++){ 
    
    
    
    tbb::enumerable_thread_specific<std::vector<int>> h_id_global([&]() {
        return std::vector<int>(totalNumSequences*2);  
    });

    tbb::enumerable_thread_specific<std::vector<int>> h_from_global([&]() {
        return std::vector<int>(totalNumSequences*2);  
    });

    tbb::enumerable_thread_specific<std::vector<double>> h_dis_global([&]() {
        return std::vector<double>(totalNumSequences*2);  
    });

    tbb::enumerable_thread_specific<std::vector<double>> h_dist_global([&]() {
        return std::vector<double>(totalNumSequences);  
    });

    tbb::enumerable_thread_specific<std::vector<int>> h_edgeMask_global([&]() {
        return std::vector<int>(totalNumSequences*4);  
    });

    tbb::enumerable_thread_specific<std::vector<int>> h_leafMask_global([&]() {
        return std::vector<int>(totalNumSequences*2);  
    });

    // i is the cluster it belongs to (we are working on)
    tbb::parallel_for(tbb::blocked_range<int>(0,numSequences*4-4), [&](tbb::blocked_range<int> range){
    for (int i = range.begin(); i < range.end(); ++i) {
        if (contains[i].size() == 0) continue;
        
        /////// VVV IMPORTANT VVV //////////
        int idx_ = params.backboneSize*4;
        int insertLeafCount_ = params.backboneSize;
        if (i>0) idx_ += containsPrefixSum[i-1]*4;
        if (i>0) insertLeafCount_ += containsPrefixSum[i-1];
        std::vector <std::tuple<int,double,double>> minPos(totalNumSequences*4-4); // move this out of for loop if not using tbb
        ////////////////////////////////////
        
        // std::vector<int>& h_edgeMask = edgeMaskGlobal.local();
        // std::vector<int>& h_leafMask = leafMaskGlobal.local();
        // h_edgeMask.assign(totalNumSequences*4, 0);
        // h_leafMask.assign(totalNumSequences*2, 0);
        // int tid = tbb::this_task_arena::current_thread_index();
        // std::vector<int>& h_edgeMask = edgeMaskBuffers[tid];
        // std::vector<int>& h_leafMask = leafMaskBuffers[tid];
        // h_edgeMask.assign(totalNumSequences*4, -1);
        // h_leafMask.assign(totalNumSequences*2, -1);
        
        
        // std::vector<int> h_edgeMask(totalNumSequences*4);
        // std::vector<int> h_leafMask(totalNumSequences*2);
        // double * h_dist_local = (double *)malloc(totalNumSequences*sizeof(double));
        std::vector<int> & h_edgeMask = h_edgeMask_global.local();
        std::vector<int> & h_leafMask = h_leafMask_global.local();
        std::vector<double> & h_dist_local = h_dist_global.local();
        std::vector<int> & h_id = h_id_global.local();
        std::vector<int> & h_from = h_from_global.local();
        std::vector<double> & h_dis = h_dis_global.local();
        
        // h_id.assign(totalNumSequences*2, -1);
        // h_from.assign(totalNumSequences*2, -1);
        // h_dis.assign(totalNumSequences*2, 2.0);

        // int * h_id = (int *)malloc(totalNumSequences*2*sizeof(int));
        // int * h_from = (int *)malloc(totalNumSequences*2*sizeof(int));
        // double * h_dis = (double *)malloc(totalNumSequences*2*sizeof(double));
        
        

        initializeClusterCpu (
            i,
            h_e,
            h_belong,
            h_head,
            h_nxt,
            h_closest_id,
            h_edgeMask,
            h_leafMask
        );
        int edgeCount=2, leafCount=10;


        for(auto &leaf:contains[i]){
            if(params.in == "r"){
                mashDeviceArrays.distSpecialIDConstructionOnCpu(
                    params,
                    leaf,
                    h_dist_local,
                    leafCount,
                    h_leafMask
                );
            }

            calculateBranchLengthSpecialIDCpu(
                i,
                h_head,
                h_nxt,
                h_dist_local,
                h_e,
                h_len,
                h_belong,
                minPos.data(),
                numSequences*4-4,
                h_closest_dis,
                h_closest_id,
                edgeCount,
                h_edgeMask
            );
            auto iter=std::min_element(minPos.begin(),minPos.begin()+edgeCount,compare_tuple());
            std::tuple<int,double,double> smallest=*iter;
            /*
            Update Tree Structure
            */

            int eid=std::get<0>(smallest);
            double fracLen=std::get<1>(smallest),addLen=std::get<2>(smallest);

            updateTreeStructureInClusterCpu (
                h_head,
                h_nxt,
                h_e,
                h_len,
                h_closest_dis,
                h_closest_id,
                h_belong,
                eid,
                fracLen,
                addLen,
                leaf,
                idx_,
                totalNumSequences,
                insertLeafCount_
            );
            idx_+=4, insertLeafCount_++;
            
            /*
            Update edgeMask and leafMask
            */

            updateClusterInfoCpu (
                leaf,
                idx_,
                h_leafMask,
                h_edgeMask,
                edgeCount,
                leafCount
                // h_edgeMaskIndex
            );
            
            /*
            Update closest nodes
            */
            
            updateClosestNodesInClusterCpu (
                h_head,
                h_nxt,
                h_e,
                h_len,
                h_closest_dis,
                h_closest_id,
                leaf,
                h_id,
                h_from,
                h_dis,
                i,
                h_belong
            );

            edgeCount+=4, leafCount++;
            
        }
        
        // free(h_edgeMask);
        // free(h_leafMask);
        // free(h_dist_local);
        // free(h_id);
        // free(h_from);
        // free(h_dis);
        
        
    // }
    }});


    auto inClusterEnd = std::chrono::high_resolution_clock::now();
    auto inClusterTime = inClusterEnd - inClusterStart;


    std::cerr << "Finished in cluster placement in: "<< inClusterTime.count()/1000000 << " ms\n";

    std::cerr<<"Finished complete tree construction\n";
}
