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
const int N=8000;
double dis[N*2+5][N*2+5];
int n;
string name[N*2+5];
vector <pair<int,double>> v[N*2+5];
int fa[N*2+5];
double dn[N*2+5],up[N*2+5],falen[N*2+5],mn,ublen,dblen;
double subd[N*2+5],upd[N*2+5],m;
int vis[N*2+5],sz[N*2+5],bd,nowmn,sel,viscnt;
double dep[N*2+5];
mt19937 rnd(time(NULL));
void dfs1(int x,int from,int nid){
    double mx=0;
    if(x<n){
        dn[x]=dis[nid][x];
        return;
    }
    for(auto t:v[x])
        if(t.first!=from)
            falen[t.first]=t.second,fa[t.first]=x,dfs1(t.first,x,nid),mx=max(mx,dn[t.first]-t.second);
    dn[x]=mx;
    // cout<<x<<" "<<dn[x]<<'\n';
}
void dfs2(int x,int from,int nid){
    for(auto t:v[x])
        if(t.first!=from){
            double mx=0;
            for(auto tt:v[x])
                if(tt.first!=t.first&&tt.first!=from)
                    mx=max(mx,dn[tt.first]-tt.second);
            mx=max(mx,up[x]-falen[x]);
            up[t.first]=mx;
            dfs2(t.first,x,nid);
        }
}
void place(int x,int from,int nid){
    // cout<<"----------  "<<x<<'\n';
    if(fa[x]!=-1){
        // cout<<x<<" ";
        double total=up[x]+dn[x];
        double addi=max(0.0,(total-falen[x])/2);
        double ub=up[x]-addi,db=dn[x]-addi;
        if(ub<0) addi-=ub,ub=0;
        if(db<0) addi-=db,db=0;
        if(addi<mn) mn=addi,sel=x,ublen=ub,dblen=db;
        // cout<<min(x,fa[x])<<" "<<max(x,fa[x])<<" "<<addi<<" "<<mn<<'\n';
    }
    if(!vis[x]){
        for(auto t:v[x])
            if(t.first!=from)
                place(t.first,x,nid);
    }
}
void place(int x){
    memset(fa,0,sizeof(fa)),up[n]=0;
    dfs1(n,-1,x),dfs2(n,-1,x);
    mn=1e6,place(n,-1,x);
    int mid=n+x-1,out=x;
    int did=sel,uid=fa[sel];
    for(auto &t:v[sel])
        if(t.first==fa[sel]){
            t.first=mid;
            t.second=dblen;
            break;
        }
    for(auto &t:v[uid])
        if(t.first==sel){
            t.first=mid;
            t.second=ublen;
            break;
        }
    assert(ublen>=0);
    assert(dblen>=0);
    assert(mn>=0);
    v[mid].push_back({did,dblen}),v[mid].push_back({uid,ublen});
    v[mid].push_back({out,mn}),v[out].push_back({mid,mn});
}
pair<pair<int,int>,double> check_place(int x,int st){
    memset(fa,-1,sizeof(fa)),up[st]=0;
    dfs1(st,-1,x),dfs2(st,-1,x);
    mn=1e6,place(st,-1,x);
    // puts("");
    // cout<<"!!!"<<mn<<'\n';
    int did=sel,uid=fa[sel];
    return make_pair(make_pair(did,uid),mn);
}
void dfssz(int x,int from,int tot){
    int mnsz=INT_MAX,used=0;
    sz[x]=x<n?1:0;
    subd[x]=0;
    for(auto t:v[x])
        if(t.first!=from&&!vis[t.first])
            dfssz(t.first,x,tot),sz[x]+=sz[t.first],mnsz=min(mnsz,sz[t.first]),used++;
    mnsz=min(mnsz,tot-1-sz[x]);
    if(from!=-1) used++;
    // cout<<x<<" "<<mnsz<<'\n';
    if(mnsz>nowmn&&used==3) nowmn=mnsz,sel=x;
}
vector <int> ct[3],allct[3];
void pb(int x,int from,int group){
    if(x<n) ct[group].push_back(x);
    allct[group].push_back(x);
    for(auto t:v[x])
        if(t.first!=from&&!vis[t.first])
            dep[t.first]=dep[x]+t.second,pb(t.first,x,group);
}
int center, samp=7, threshold=5;
int vote(int x,std::vector<int> c1,std::vector<int> c2,std::vector<int> c3){
    std::vector<int> a[3]={c1,c2,c3};
    int best=-1;
    double bl[3]={dep[allct[0][0]],dep[allct[1][0]],dep[allct[2][0]]};
    double mn=1e9;
    for(int i=0;i<3;i++){
        double l1=0;
        for(auto t:a[i]) l1=std::max(l1,dis[t][x]-(dep[t]-bl[i]));
        double l2=0;
        for(int j=0;j<3;j++) if(j!=i) for(auto t:a[j]) l2=std::max(l2,dis[t][x]-dep[t]);
        double total=l1+l2;
        double addi=max(0.0,(total-bl[i])/2);
        double ub=l1-addi,db=l2-addi;
        // if(ub<0) addi-=ub,ub=0;
        // if(db<0) addi-=db,db=0;
        if(addi<mn) mn=addi,best=i,ublen=ub,dblen=db;
        // cout<<x<<" "<<i<<" "<<addi<<" "<<dis[a[i]][x]<<" "<<dep[a[i]]<<'\n'; 
    }
    return best;
}

int vote(int x,int c1,int c2,int c3){
    int a[3]={c1,c2,c3},best=-1;
    double mn=1e9;
    for(int i=0;i<3;i++){
        double l1=dis[a[i]][x];
        double l2=0;
        for(int j=0;j<3;j++) if(j!=i) l2=max(l2,dis[a[j]][x]-dep[a[j]]);
        double total=l1+l2;
        double addi=max(0.0,(total-dep[a[i]])/2);
        double ub=l1-addi,db=l2-addi;
        if(ub<0) addi-=ub,ub=0;
        if(db<0) addi-=db,db=0;
        if(addi<mn) mn=addi,best=i,ublen=ub,dblen=db;
        // cout<<x<<" "<<i<<" "<<addi<<" "<<dis[a[i]][x]<<" "<<dep[a[i]]<<'\n'; 
    }
    return best;
}

int divide_and_conquer(int x,int st,int tot){
    // cout<<"Start checking a size of "<<tot<<" from "<<st<<'\n';
    nowmn=0;
    dfssz(st,-1,tot);
    assert(v[sel].size()==3);
    int cen=sel;
    for(int i=0;i<3;i++)
        allct[i].clear(),ct[i].clear(),dep[v[sel][i].first]=v[sel][i].second,pb(v[sel][i].first,sel,i),shuffle(ct[i].begin(),ct[i].end(),rnd);
    auto temp=check_place(x,sel);
    // cout<<"Correct is ------"<<temp.first.first<<" "<<temp.first.second<<'\n';
    // cout<<"Total size is -------"<<tot<<'\n';
    // cout<<"Three subtrees are "<<allct[0][0]<<" "<<allct[1][0]<<" "<<allct[2][0]<<'\n';
    // cout<<"Leave sizes are "<<ct[0].size()<<" "<<ct[1].size()<<" "<<ct[2].size()<<'\n';
    // int correct=-1;
    // for(int i=0;i<3;i++)
    //     for(auto t:allct[i])
    //         if(t==temp.first.first||t==temp.first.second)
    //             correct=i;
    int cnt=0,mx=0,subtreeid=-1,c[3]={};
    std::vector<int> sp[3];
    int sz=min(samp,int(min(min(ct[0].size(),ct[1].size()),ct[2].size())));
    for(int i=0;i<sz;i++)
        for(int j=0;j<3;j++)
            sp[j].push_back(ct[j][i]);
    // c[vote(x,sp[0],sp[1],sp[2])]++;
    for(int i=0;i<min(samp,int(min(min(ct[0].size(),ct[1].size()),ct[2].size())));i++){
        c[vote(x,ct[0][i],ct[1][i],ct[2][i])]++;
    }
    for(int i=0;i<3;i++) if(c[i]>mx) mx=c[i],subtreeid=i;
    if(ct[subtreeid].size()<=threshold) return cen;
    else{
        vis[cen]=1;
        return divide_and_conquer(x,allct[subtreeid][0],ct[subtreeid].size());
    }
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
                print(t.first,x),printf(":%.6lf",t.second),assert(t.second>=0);
            }
        printf(")");
    }
}
int main(){
    cin>>n;
    for(int i=0;i<n;i++){
        cin>>name[i];
        for(int j=0;j<n;j++) scanf("%lf",&dis[i][j]);
    }
    double st=dis[0][1];
    v[0].push_back({n,st/2}),v[1].push_back({n,st/2});
    v[n].push_back({0,st/2}),v[n].push_back({1,st/2});
    bd=n/10;
    threshold=max(5,n/50);
    samp=1000;
    for(int i=2;i<bd;i++){
        place(i);
        // cout<<"#######"<<i<<'\n';
        // for(int j=0;j<=i;j++) for(auto t:v[j]) cout<<j<<" "<<t.first<<'\n';
        // for(int j=n;j<=n+i-1;j++) for(auto t:v[j]) cout<<j<<" "<<t.first<<'\n';
    }
    int correct=0,choice;
    for(int i=bd;i<n;i++){
        memset(vis,0,sizeof(vis));
        choice=divide_and_conquer(i,n,i);

        assert(v[choice].size()>1);
        // cout<<"Start checking..............."<<'\n';
        auto divans=check_place(i,choice);
        int mid=n+i-1,out=i;
        int did=divans.first.first,uid=divans.first.second;
        for(auto &t:v[did])
            if(t.first==uid){
                t.first=mid;
                t.second=dblen;
                break;
            }
        for(auto &t:v[uid])
            if(t.first==did){
                t.first=mid;
                t.second=ublen;
                break;
            }
        assert(ublen>=0);
        assert(dblen>=0);
        assert(mn>=0);
        v[mid].push_back({did,dblen}),v[mid].push_back({uid,ublen});
        v[mid].push_back({out,mn}),v[out].push_back({mid,mn});

        // memset(vis,0,sizeof(vis));
        // auto temp=check_place(i,n);
        // cout<<"End checking..................."<<'\n';
        // cout<<divans.first.first<<" "<<divans.first.second<<" "<<divans.second<<'\n';
        // cout<<temp.first.first<<" "<<temp.first.second<<" "<<temp.second<<'\n'<<'\n';
        // if(divans==temp) correct++;
    }
    // printf("%d out of %d is correct\n", correct, n-bd);
    print(n,-1),printf(";");
    return 0;
}