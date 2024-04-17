#ifndef KPLEX_GRAPH_H
#define KPLEX_GRAPH_H
#include "utils.h"


enum : uint8_t
{ // 奇数表示连接
    UNLINK2LESS = 0,
    LINK2LESS = 1,
    UNLINK2EQUAL = 2,
    LINK2EQUAL = 3,
    UNLINK2MORE = 4,
    LINK2MORE = 5
};
enum MODE
{
    INIT,
    LBCHANGE,
    NONE
};


class KPlexGraph 
{
public:
    vector<vecui> adjList;
    ui V, E;

    KPlexGraph(ui n)
    {
        V = 0;
        E = 0;
        adjList.resize(n);
        cnMat.resize(n * n);
        adjMat.resize(n*n);
        cn.resize(n);
        initVectors(n);
    }

    int k, tau_e, tau_v;
    vector<vecui> cn;
    vector<pair<ui, ui>> Qe;
    vecui Qv;
    vecui degree, minv;
    vecui pruned;
    vecui lookup;
    vecui kplex;
    vecui peelSeq;
    vecui rec;
    MBitSet vis;
    vecui cnMat, adjMat;
    ui rv = 0, re = 0, n;
    ui lb = 0, minindex = 0;
    // ListLinearHeap heap;

    KPlexGraph(ui _k, ui n) : k(_k)
    {
        // this constructor is called for two-hop graph gi, hence allocating such large memory here.
        adjList.resize(n);
        for(auto& adj: adjList)
            adj.reserve(n);
        cnMat.resize(n * n);
        cn.resize(n);
        initVectors(n);
    }
    void load(const std::vector<std::pair<int,int> > &vp, ui _n) {
        n = _n;
        V = n;
        for(ui i = 0;i < vp.size();i ++) {
            assert(vp[i].first >= 0&&vp[i].first < n&&vp[i].second >= 0&&vp[i].second < n);
        	ui a = vp[i].first, b = vp[i].second;
            adjMat[a*n + b] = adjMat[b*n + a] = 1;
        }
        for(ui i=0;i<n;i++)
            for(ui j=0;j<n;j++){
                if(adjMat[i*n+j])
                    adjList[i].push_back(j);
            }
    }

    void unload(const std::vector<std::pair<int,int> > &vp, ui _n) {
        n = _n;
        for(ui i = 0;i < vp.size();i ++) {
            assert(vp[i].first >= 0&&vp[i].first < n&&vp[i].second >= 0&&vp[i].second < n);
        	ui a = vp[i].first, b = vp[i].second;
            adjMat[a*n + b] = adjMat[b*n + a] = 0;
        }
        for(ui i=0;i<n;i++)
            adjList[i].clear();
        V = 0;
    }


    void initMin()
    {
        minv.reserve(V);
        for (ui i = 0; i < V; i++)
        {
            // cout<<degree[i]<<" ";
            if (degree[i] != 0)
            {
                minv.push_back(i);
            }
        }
    }
    void initVectors(ui sz)
    {
        peelSeq.resize(sz);
        degree.resize(sz);
        lookup.resize(sz);
        pruned.resize(sz);
        Qv.reserve(sz);
        Qe.reserve(sz); // todo need to comeup with better estimate for Qe.
        vis.init(sz);
    }
        void addEdge(ui u, ui v)
    {
        adjList[u].push_back(v);
        adjList[v].push_back(u);

    }
    void resize(ui n)
    {
        V = n;
        adjList.resize(n);
    }
    ui addVertex(){
        return V++;
    }

    void sort()
    {
        for (auto &adj : adjList)
            std::sort(adj.begin(), adj.end());
    }
};
#endif // KPLEX_GRAPH_H
