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

    KPlexGraph(ui n) : heap(n, n - 1)
    {
        V = 0;
        E = 0;
        adjList.resize(n);
        cnMat.resize(n * n);
        adjMat.resize(n * n);
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
    vecui peelSeq;
    vecui core;
    vecui rec;
    MBitSet vis;
    vecui cnMat, adjMat;
    ui rv = 0, re = 0, n;
    ui lb = 0, minindex = 0;
    ListLinearHeap heap;

    KPlexGraph(ui _k, ui n) : k(_k), heap(n, n - 1)
    {
        // this constructor is called for two-hop graph gi, hence allocating such large memory here.
        adjList.resize(n);
        for (auto &adj : adjList)
            adj.reserve(n);
        cnMat.resize(n * n);
        cn.resize(n);
        initVectors(n);
    }
    void print(){
        for(ui i=0;i<V;i++){
            cout<<i<<"["<<adjList[i].size()<<"]: ";
            // for(ui u: adjList[i])
            //     cout<<u<<" ";
            cout<<endl;
        }
    }
    ui degeneracyKPlex(vecui& kplex)
    {
        Timer t;
        vis.reset();
        ui ub = 0;
        heap.init(V, adjList);
        ui max_core = 0;
        ui idx = V;
        for (ui i = 0; i < V; i++)
        {
            ui u, key;
            heap.pop_min(u, key);
            if (key > max_core)
                max_core = key;
            peelSeq[i] = u;
            core[u]=max_core;
            // ui t_UB = min(max_core + k, V - i);
            // if (V - i < t_UB)
            //     t_UB = V - i;
            // if (t_UB > ub)
            //     ub = t_UB;
            ub = max(ub, min(max_core + k, V - i));
            if (idx == V && key + k >= V - i)
                idx = i;
            vis.set(u);

            for (ui v : adjList[u])
                if (!vis.test(v))
                    heap.decrement(v, 1);
        }

        printf("*** Degeneracy k-plex size: %u, max_core: %u, ub: %u, Time: %lu (microseconds)\n", V , max_core, ub, t.elapsed());

        if (V - idx > kplex.size())
        {
            kplex.clear();
            for (ui i = idx; i < V; i++)
                kplex.push_back(peelSeq[i]);
        }
        return ub;
    }

    void load(const std::vector<std::pair<int, int>> &vp, ui _n)
    {
        V = n = _n;
        for (ui i = 0; i < vp.size(); i++)
        {
            assert(vp[i].first >= 0 && vp[i].first < n && vp[i].second >= 0 && vp[i].second < n);
            ui a = vp[i].first, b = vp[i].second;
            adjMat[a * n + b] = adjMat[b * n + a] = 1;
        }
        for (ui i = 0; i < n; i++)
            for (ui j = 0; j < n; j++)
            {
                if (adjMat[i * n + j])
                    adjList[i].push_back(j);
            }
    }

    void unload(const std::vector<std::pair<int, int>> &vp, ui _n)
    {
        n = _n;
        for (ui i = 0; i < vp.size(); i++)
        {
            assert(vp[i].first >= 0 && vp[i].first < n && vp[i].second >= 0 && vp[i].second < n);
            ui a = vp[i].first, b = vp[i].second;
            adjMat[a * n + b] = adjMat[b * n + a] = 0;
        }
        for (ui i = 0; i < n; i++)
            adjList[i].clear();
    }

    ui &cnMatrix(ui i, ui j)
    {
        return cnMat[i * V + j];
    }
    void buildCommonMatrix(ui sz1h)
    {

        auto intersect = [](const vecui &lst1, const vecui &lst2, const ui sz)
        {
            int i = 0, j = 0;
            int common = 0;
            while (i < lst1.size() and j < lst2.size() and lst1[i] < sz and lst2[j] < sz)
            {
                if (lst1[i] < lst2[j])
                    ++i;
                else if (lst1[i] > lst2[j])
                    ++j;
                else
                {
                    if (lst1[i] != 0)
                        common++;
                    i++, j++;
                }
            }
            return common;
        };

        // 0 is seed vertex, 1...sz1h is 1-hop vertices, sz1h+1...end is 2-hop
        // const ui sz1h = adjList[0].size() + 1;
        const int thresPP1 = lb - k - 2 * max(k - 2, 0), thresPP2 = lb - k - 2 * max(k - 3, 0);
        const int thresPC1 = lb - 2 * k - max(k - 2, 0), thresPC2 = lb - k - 1 - max(k - 2, 0) - max((int)k - 3, 0);
        const int thresCC1 = lb - 2 * k - (k - 1), thresCC2 = lb - 2 * k + 2 - (k - 1);
        // 两个1-hop或本身
        for (int i = 1; i < sz1h; ++i)
        {
            for (int j = 1; j < i; ++j)
            {
                const int common = intersect(adjList[i], adjList[j], sz1h);

                if (isNeighbor(adjList[i], j))
                {
                    if (common > thresCC1)
                        cnMatrix(i, j) = LINK2MORE;
                    else if (common == thresCC1)
                        cnMatrix(i, j) = LINK2EQUAL;
                    else
                        cnMatrix(i, j) = LINK2LESS;
                }
                else
                {
                    if (common > thresCC2)
                        cnMatrix(i, j) = UNLINK2MORE;
                    else if (common == thresCC2)
                        cnMatrix(i, j) = UNLINK2EQUAL;
                    else
                        cnMatrix(i, j) = UNLINK2LESS;
                }
                cnMatrix(j, i) = cnMatrix(i, j);
            }
        }
        // 一个1-hop 一个2-hop
        for (int i = sz1h; i < V; ++i)
        {
            for (int j = 1; j < sz1h; ++j)
            {
                const int common = intersect(adjList[i], adjList[j], sz1h);

                if (isNeighbor(adjList[i], j))
                {
                    if (common > thresPC1)
                        cnMatrix(i, j) = LINK2MORE;
                    else if (common == thresPC1)
                        cnMatrix(i, j) = LINK2EQUAL;
                    else
                        cnMatrix(i, j) = LINK2LESS;
                }
                else
                {
                    if (common > thresPC2)
                        cnMatrix(i, j) = UNLINK2MORE;
                    else if (common == thresPC2)
                        cnMatrix(i, j) = UNLINK2EQUAL;
                    else
                        cnMatrix(i, j) = UNLINK2LESS;
                }
                cnMatrix(j, i) = cnMatrix(i, j);
            }
            if (k == 2)
                continue; // why
            // 两个2-hop点在lead顶点的1-hop集合中的交集数量
            for (int j = sz1h; j < i; ++j)
            {
                const int common = intersect(adjList[i], adjList[j], sz1h);

                if (isNeighbor(adjList[i], j))
                {
                    if (common > thresPP1)
                        cnMatrix(i, j) = LINK2MORE;
                    else if (common == thresPP1)
                        cnMatrix(i, j) = LINK2EQUAL;
                    else
                        cnMatrix(i, j) = LINK2LESS;
                }
                else
                {
                    if (common > thresPP2)
                        cnMatrix(i, j) = UNLINK2MORE;
                    else if (common == thresPP2)
                        cnMatrix(i, j) = UNLINK2EQUAL;
                    else
                        cnMatrix(i, j) = UNLINK2LESS;
                }
                cnMatrix(j, i) = cnMatrix(i, j);
            }
        }
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
        core.resize(sz);
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
    ui addVertex()
    {
        return V++;
    }

    void sort()
    {
        for (auto &adj : adjList)
            std::sort(adj.begin(), adj.end());
    }
    void clear()
    {
        for (ui u = 0; u < V; u++)
        {
            cn[u].clear();
            pruned[u] = 0;
        }
        peelSeq.clear();
        // for(ui u=0;u<V; u++)
        //     adjList[u].clear();
        V = 0;
    }
};
#endif // KPLEX_GRAPH_H
