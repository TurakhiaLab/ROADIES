#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace std;
const int N=16000;
double dis[N+5][N+5];
int n;
string name[N*2+5];
vector <pair<int,pair<double,int>>> v[N*2+5];
int fa[N*2+5];
double dn[N*2+5],up[N*2+5],falen[N*2+5],mn,ublen,dblen;
double subd[N*2+5],upd[N*2+5],m;
int vis[N*2+5],sz[N*2+5],bd,nowmn,sel,viscnt;
int bel[N*2+5];
double dep[N*2+5];
int global_cnt;
int att;

vector <int> to_be_placed[N*2+5];

mt19937 rnd(time(NULL));

void dfs1(int x,int from,int nid){
    double mx=0;
    if(x<n){
        dn[x]=dis[nid][x];
        return;
    }
    for(auto t:v[x])
        if(t.first!=from)
            falen[t.first]=t.second.first,fa[t.first]=x,dfs1(t.first,x,nid),mx=max(mx,dn[t.first]-t.second.first),bel[t.first]=t.second.second;
    dn[x]=mx;
    // cout<<x<<" "<<dn[x]<<'\n';
}
void dfs2(int x,int from,int nid){
    for(auto t:v[x])
        if(t.first!=from){
            double mx=0;
            for(auto tt:v[x])
                if(tt.first!=t.first&&tt.first!=from)
                    mx=max(mx,dn[tt.first]-tt.second.first);
            mx=max(mx,up[x]-falen[x]);
            up[t.first]=mx;
            dfs2(t.first,x,nid);
        }
}
void place(int x,int from,int nid,int cid){
    // cout<<"----------  "<<x<<'\n';
    if(fa[x]!=-1){
        // cout<<x<<" ";
        double total=up[x]+dn[x];
        double addi=max(0.0,(total-falen[x])/2);
        double ub=up[x]-addi,db=dn[x]-addi;
        if(ub<0) addi-=ub,ub=0;
        if(db<0) addi-=db,db=0;
        if(addi<mn&&(cid==-1||cid==bel[x])){
            mn=addi,sel=x,ublen=ub,dblen=db;
            // if(sel<=1) cerr<<x<<" "<<from<<'\n';
        }
        if(cid==-1||cid==bel[x]) att++;
        // cout<<min(x,fa[x])<<" "<<max(x,fa[x])<<" "<<addi<<" "<<mn<<'\n';
    }
    if(!vis[x]){
        for(auto t:v[x])
            if(t.first!=from)
                place(t.first,x,nid,cid);
    }
}

void place(int x,int rt,int cid){
    // memset(fa,0,sizeof(fa)),up[rt]=0;
    fa[rt]=-1,up[rt]=0;
    dfs1(rt,-1,x),dfs2(rt,-1,x);
    att=0,mn=1e6,place(rt,-1,x,cid);
    // cerr<<x<<" "<<att<<" "<<cid<<'\n';
    // cerr<<sel<<'\n';
    // cerr<<fa[sel]<<'\n';
    // if(rt==1399) cerr<<sel<<" "<<x<<'\n';
    int mid=++global_cnt,out=x;
    int did=sel,uid=fa[sel];
    // cerr<<mid<<' '<<out<<'\n';
    if(cid==-1){
        for(auto &t:v[sel])
            if(t.first==fa[sel]){
                t.first=mid;
                t.second.first=dblen;
                t.second.second=sel;
                break;
            }
        for(auto &t:v[uid])
            if(t.first==sel){
                t.first=mid;
                t.second.first=ublen;
                t.second.second=mid;
                break;
            }
        // cerr<<"Reached here\n";
        assert(ublen>=0);
        assert(dblen>=0);
        assert(mn>=0);
        v[mid].push_back({did,{dblen,did}}),v[mid].push_back({uid,{ublen,mid}});
        v[mid].push_back({out,{mn,out}}),v[out].push_back({mid,{mn,out}});
    }
    else{
        assert(bel[sel]==cid);
        for(auto &t:v[sel])
            if(t.first==fa[sel]){
                t.first=mid;
                t.second.first=dblen;
                t.second.second=cid;
                break;
            }
        for(auto &t:v[uid])
            if(t.first==sel){
                t.first=mid;
                t.second.first=ublen;
                t.second.second=cid;
                break;
            }
        // cerr<<"Reached here\n";
        assert(ublen>=0);
        assert(dblen>=0);
        assert(mn>=0);
        v[mid].push_back({did,{dblen,cid}}),v[mid].push_back({uid,{ublen,cid}});
        v[mid].push_back({out,{mn,cid}}),v[out].push_back({mid,{mn,cid}});
    }
}

void find_place_branch(int x,int rt){
    // memset(fa,0,sizeof(fa)),up[rt]=0;
    fa[rt]=-1,up[rt]=0;
    dfs1(rt,-1,x),dfs2(rt,-1,x);
    mn=1e6,place(rt,-1,x,-1);
    // cout<<sel<<" "<<x<<'\n';
}

void print(int x,int from){
    if(x<n) cout<<name[x];
    else{
        printf("(");
        int fst=0,f=0;
        for(auto t:v[x])
            if(t.first!=from){
                if(!f) f=1;
                else printf(",");
                print(t.first,x),printf(":%.6lf",t.second.first),assert(t.second.first>=0);
            }
        printf(")");
    }
}

int build_placement_tree(vector<int> nodes){
    // cerr<<nodes.size()<<'\n';
    if(nodes.size()==0) return -1;
    if(nodes.size()==1) return nodes[0];
    double st=dis[nodes[0]][nodes[1]];
    int rt=++global_cnt;
    v[nodes[0]].push_back({rt,{st/2,nodes[0]}}),v[nodes[1]].push_back({rt,{st/2,nodes[1]}});
    v[rt].push_back({nodes[0],{st/2,nodes[0]}}),v[rt].push_back({nodes[1],{st/2,nodes[1]}});
    for(int i=2;i<nodes.size();i++){
        place(nodes[i],rt,-1);
        // cerr<<"nodes count "<<i<<" out of "<<nodes.size()<<'\n';
        // if(rt==1399) print(rt,-1),puts("");
    }
    return rt;
}

vector <pair<int,int>> udpair;

void final_place(int x,int from){
    // cerr<<x<<" "<<from<<'\n';
    for(auto t:v[x])
        if(t.first!=from){
            final_place(t.first,x);
            udpair.push_back({t.first,x});
        }
}

map <string,vector<pair<string,double>>> reftree;
vector <string> path;
string reference_tree;
vector <int> oth;
int ccc;

// string create_ref(int l,int r){
//     if(reference_tree[l]=='('){
//         string myname="useless"+to_string(ccc);
//         ccc++;
//         int now=l+1;
//         while(now<r){
//             if(reference_tree[now]=='('){
//                 string temp=create_ref(now,oth[now]);
//                 now=oth[now]+2;
//                 string num="";
//                 while(reference_tree[now]!=','&&reference_tree[now]!=')'&&reference_tree[now]!=';'){
//                     num+=reference_tree[now++];
//                 }
//                 now++;
//                 reftree[myname].push_back({temp,stod(num)});
//             }
//             else{
//                 // if(reference_tree[now]!='i'){
//                 //     cerr<<reference_tree[now]<<'\n';
//                 //     for(int i=now-2;i<=now+2;i++) cerr<<reference_tree[i];
//                 //     cerr<<'\n';
//                 //     cerr<<l<<" "<<r<<'\n';
//                 //     cerr<<oth[now-1]<<" "<<now-1<<'\n';
//                 // }
//                 int end=now;
//                 while(reference_tree[end]!=':') end++;
//                 end--;
//                 string temp=create_ref(now,end);
//                 now=end+2;
//                 string num="";
//                 while(reference_tree[now]!=','&&reference_tree[now]!=')'&&reference_tree[now]!=';'){
//                     num+=reference_tree[now++];
//                 }
//                 now++;
//                 reftree[myname].push_back({temp,stod(num)});
//                 // cerr<<temp<<" "<<num<<'\n';
//                 // cerr<<num<<" "<<stod(num)<<'\n';
//             }
//         }
//         return myname;
//     }
//     else{
//         int now=l;
//         string myname="";
//         while(now<=r) myname+=reference_tree[now],now++;
//         return myname;
//     }
// }

map <string,int> clustered;
std::ofstream* outFile = nullptr;
std::ofstream* outFile2 = nullptr;

// int special_print_ref(string x){
//     if(x[0]=='i'){
//         *outFile<<x;
//         if(clustered.count(x)) return 0;
//         else return 1;
//     }
//     else{
//         int tot=0;
//         *outFile<<'(';
//         int fst=0,f=0;
//         for(auto t:reftree[x]){
//             if(!f) f=1;
//             else *outFile<<',';
//             int temp=special_print_ref(t.first);
//             tot+=temp;
//             // cerr<<t.second<<'\n';
//             if(temp){
//                 *outFile<<":"<<fixed<<setprecision(8)<<t.second;
//             }
//             else *outFile<<":0.0";
//         }
//         *outFile<<')';
//         return tot;
//     }
// }

// int special_print(int x,int from){
//     if(x<n){
//         *outFile2<<name[x];
//         if(clustered.count(name[x])) return 0;
//         else return 1;
//     }
//     else{
//         int tot=0;
//         *outFile2<<'(';
//         int fst=0,f=0;
//         for(auto t:v[x]) if(t.first!=from){
//             if(!f) f=1;
//             else *outFile2<<',';
//             int temp=special_print(t.first,x);
//             tot+=temp;
//             // cerr<<t.second<<'\n';
//             if(temp){
//                 *outFile2<<":"<<fixed<<setprecision(8)<<t.second;
//             }
//             else *outFile2<<":0.0";
//         }
//         *outFile2<<')';
//         return tot;
//     }
// }

int main(){
    cin>>n;
    for(int i=0;i<n;i++){
        cin>>name[i];
        for(int j=0;j<n;j++) scanf("%lf",&dis[i][j]);
    }

    global_cnt=n-1;
    vector <int> all;
    for(int i=0;i<n;i++) all.push_back(i);
    shuffle(all.begin(),all.end(),rnd);

    int backbone_size=10;
    vector <int> backbone;
    for(int i=0;i<backbone_size;i++){
        backbone.push_back(all[i]);
    }
    assert(build_placement_tree(backbone)==n);

    // cerr<<"Backbone is done\n";

    for(int i=backbone_size;i<n;i++){
        find_place_branch(all[i],n);
        to_be_placed[sel].push_back(all[i]);
        // if(name[all[i]]=="id_522") cerr<<sel<<" "<<all[i]<<'\n';
        // cerr<<sel<<" "<<all[i]<<'\n';
    }

    // cerr<<"Clustering is done\n";

    final_place(n,-1);

    // cerr<<"UDpair found\n";
    // cerr<<global_cnt<<'\n';
    int sum=0;
    for(int i=0;i<global_cnt;i++) sum+=to_be_placed[i].size();
    // cerr<<sum<<" "<<global_cnt<<'\n';
    for(auto temp:udpair){
        int did=temp.first,uid=temp.second;
        for(auto x:to_be_placed[did]) place(x,n,did);
    }

    // ////////////////////////////////////
    // // Start creating clustered list
    // clustered.clear();
    // for(int i=backbone_size;i<n;i++) clustered[name[all[i]]]=1;
    // cerr<<clustered.size()<<'\n';

    // ////////////////////////////////////
    // // Start processing reference tree
    // // Need to modify this line
    // ifstream file("/data/zec022/phastsim_datasets/dataset_8000/sars-cov-2_simulation_output.tree"); 

    // getline(file,reference_tree);
    // // cerr<<reference_tree<<'\n';
    // file.close(); 

    // string clean_reference_tree="";
    // int flag=0;
    // for(int i=0;i<reference_tree.size();i++){
    //     if(reference_tree[i]=='[') flag=1;
    //     if(flag==0&&reference_tree[i]!=';'){
    //         clean_reference_tree+=reference_tree[i];
    //     }
    //     if(reference_tree[i]==']') flag=0;
    // }
    // reference_tree=clean_reference_tree;
    // reference_tree.pop_back();
    // reference_tree="("+reference_tree+");";
    // oth.resize(reference_tree.size());
    // vector <int> pos;
    // for(int i=0;i<reference_tree.size();i++){
    //     if(reference_tree[i]=='(') pos.push_back(i);
    //     else if(reference_tree[i]==')'){
    //         oth[pos.back()]=i;
    //         oth[i]=pos.back();
    //         pos.pop_back();
    //     }
    // }
    // ccc=0;
    // string rt=create_ref(0,reference_tree.size()-1);
    // outFile = new std::ofstream("/data/zec022/phastsim_datasets/dataset_8000/sars-cov-2_simulation_output_polytomy.tree");
    // special_print_ref(rt);
    // *outFile<<';';
    // delete outFile;

    // // Finished processing reference tree

    // /////////////////////////////////////
    // // Start processing generated tree

    // outFile2 = new std::ofstream("/data/zec022/phastsim_datasets/dataset_8000/generated_polytomy.tree");
    // special_print(n,-1);
    // *outFile2<<';';
    // delete outFile2;

    // // Finished processing generated tree


    print(n,-1),printf(";");

    return 0;
}