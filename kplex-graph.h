#ifndef KPLEX_GRAPH_H
#define KPLEX_GRAPH_H
#include <stdc++.h>
#include "utils.h"
using namespace std;
using namespace chrono;
typedef unsigned int ui;
typedef vector<ui> vecui;
typedef pair<ui, ui> Edge;

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
    void degreePrune()
    {
        while (!Qv.empty())
            removeVertex();
    }
    void maxDegreeKplex(ui thresh);
    void corePrune(ui);
    void applyCTCP(MODE);
    void initCTCP();
    void trussPrune();
    void removeVertex();
    void decDegree(ui);
    void decCN(Edge);
    void populate();

    Edge getEdge(ui u, ui j)
    {
        ui v = adjList[u][j];
        if (u > v)
            return {u, j};
        else
        {
            ui k = getLowerBound(adjList[v], u);
            return {v, k};
        }
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
    bool exists(Edge e)
    {
        return cn[e.first][e.second] != V + 1;
    }

    bool exists(ui u)
    {
        return degree[u];
        // return not pruned[u];
    }

    void removeEdge(Edge e)
    {
        cn[e.first][e.second] = V + 1;
    }

    void scheduleRemoval(Edge e)
    {
        if (cn[e.first][e.second] != V)
        {
            Qe.push_back(e);
            cn[e.first][e.second] = V;
        }
    }
    void scheduleRemoval(ui v)
    {
        if (pruned[v])
            return;
        pruned[v] = 1;
        Qv.push_back(v);
    }
    ui &getCN(Edge e)
    {
        return cn[e.first][e.second];
    }

    void shrinkEdges()
    {
        for (ui u = 0; u < V; u++)
        {
            ui j = 0;
            for (ui i = 0; i < adjList[u].size(); i++)
            {
                ui v = adjList[u][i];
                if (exists(v) and exists(getEdge(u, i)))
                    adjList[u][j++] = v;
            }
            adjList[u].resize(j);
            pruned[u] = 0;
        }
    }

    void shrink(bool checkEdge = false)
    {
        rec.resize(V);
        ui rid = 0;
        E = 0;
        for (ui i = 0; i < V; i++)
        {
            if (exists(i))
                rec[i] = rid++;
            // recoded ids are also in increasing order in any adjacency list. hence no need to sort after shrinking
        }

        for (ui u = 0; u < V; u++)
        {
            ui j = 0;
            if (not exists(u))
                continue;
            for (ui i = 0; i < adjList[u].size(); i++)
            {
                ui v = adjList[u][i];
                if (exists(v))
                    if (!checkEdge or exists(getEdge(u, i)))
                        adjList[u][j++] = rec[v];
            }

            adjList[u].resize(j);
            // if(j>100) cout<<j<<" ";
            E += j;
        }
        // removing zero degree vertices.
        ui nu = 0;
        // check2.tick();
        for (ui u = 0; u < V; u++)
        {
            if (checkEdge)
                cn[u].clear();
            if (exists(u))
                swap(adjList[u], adjList[nu++]);
            // adjList[u].swap(adjList[nu++]);
            pruned[u] = 0;
        }
        // check2.tock();
        V = nu;
        E = 0;
        for (ui u = 0; u < V; u++)
            E += adjList[u].size();
        cout << nu << " rcoded " << E / 2 << endl;
        E /= 2;
        rec.clear();
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
};
void KPlexGraph::corePrune(ui q)
{
    degree.clear();
    for (ui i = 0; i < V; i++)
        degree.push_back(adjList[i].size());

    for (ui i = 0; i < V; i++)
        if (degree[i] < q)
        {
            Qv.push_back(i);
            degree[i] = 0;
        }
    for (ui i = 0; i < Qv.size(); i++)
    {
        ui u = Qv[i];
        for (ui v : adjList[u])
            if (exists(v) and degree[v]-- == q)
            {
                Qv.push_back(v);
                degree[v] = 0;
            }
        adjList[u].clear();
        degree[u] = 0;
    }
    Qv.clear();
    shrink();
}
void KPlexGraph::initCTCP()
{
    // number of common neighbors are stored in one direction only i.e. for an edge (u, v)
    // where u > v, the cn[u][ind(v)] stores number of common neighbors of (u, v)
    // tau_e = lb + 1 - 2 * k;
    // tau_v = lb + 1 - k;
    if (cn.empty())
        cn.resize(V); // this condition is true only for G graph, the gi graph will allocate cn in its constructor
    for (ui u = 0; u < V; u++)
    {
        degree[u] = adjList[u].size();
        ui sz = 0;
        for (ui v : adjList[u])
            if (u > v)
                sz++;
            else
                break;
        cn[u].clear();
        cn[u].resize(sz);
        // fill(cn[u].begin(), cn[u].begin()+sz, 0);
    }
    cnMat.resize(V*V);
    vecui edges, degrees(V, 0), stp;
    edges.reserve(E);
    stp.reserve(V + 1);
    stp.push_back(0);
    ui sum = 0;
    degenerate();
    for (ui i = 0; i < V; i++)
    {
        ui j = 0;
        for(ui v:adjList[i])
            if(peelSeq[v]>peelSeq[i])
                edges.push_back(v), j++;
        degrees[i] = j;
        sum += j;
        stp.push_back(edges.size());
    }
    cout << "sum is: " << edges.size()<<" "<<sum << endl;
    // printvec("deg", degrees);
    vector<ui> cnm(V, 0);
    // counting common neighbors
    check.tick();

    // for (ui u = 0; u < V; u++)
    // {
    //     Lookup neigh(&lookup, &adjList[u]);
    //     for (ui i = stp[u]; i < stp[u+1]; i++)
    //     {
    //         ui v = edges[i];
    //         for (ui j = stp[v]; j < stp[v+1]; j++)
    //         {
    //             ui w = edges[j];
    //             cnm[u]++;
    //         }
    //     }
    // }

    // ui sum=0;
    for (ui u = 0; u < V; u++)
    {
        for (ui i = 0; i < degrees[u]; i++)
        {
            // ui v = adjList[u][i];
            ui v = edges[stp[u] + i];
            // if(v>u) break;
            lookup[v] = i + 1;
        }

        for (ui i = 0; i < degrees[u]; i++)
        {
            // ui v = adjList[u][i];
            ui v = edges[stp[u] + i];
            // if (v > u)
            //     break;
            for (ui j = 0; j < degrees[v]; j++)
            {
                // ui w = adjList[v][j];
                ui w = edges[stp[v] + j];
                // if (w > v)
                //     break;
                if (lookup[w])
                {
                    cnMat[u*V+i]++;
                    cnMat[i*V+i]++;
                    cnMat[j*V+i]++;
                    // cn[u][i]++;
                    // cn[u][j]++;
                    // cn[u][lookup[w] - 1]++;
                }
            }
        }
        for (ui i = 0; i < degrees[u]; i++)
        {
            ui v = edges[stp[u] + i];
            // ui v = adjList[u][i];
            // if(v>u) break;
            lookup[v] = 0;
        }
    }
    cout << sum;
    check.tock();
    cout << "initialized in " << check.ticktock() << endl;
}

void KPlexGraph::populate()
{
    tau_e = lb + 1 - 2 * k;
    tau_v = lb + 1 - k;
    // cout << "tau_v: " << tau_v << " tau_e: " << tau_e << endl;
    // populate Qe based on cn[u][v] < lb+1-2k
    for (ui u = 0; u < V; u++)
    {
        if (degree[u] < tau_v and exists(u))
        {
            scheduleRemoval(u);
            continue;
        }
        for (ui j = 0; j < adjList[u].size(); j++)
        {
            if (adjList[u][j] > u)
                break;
            if (cn[u][j] < tau_e)
            {
                scheduleRemoval({u, j});
            }
        }
    }
    // cout << "intialization Qe: " << Qe.size() << endl;
}

void KPlexGraph::applyCTCP(MODE md)
{

    if (md == INIT)
    {
        initCTCP();
        populate();
    }
    else if (md == LBCHANGE)
        populate();
    // printvec("degree: ", degree);
    while (true)
    {
        // if (!Qe.empty())
        //     trussPrune();
        // else
        if (!Qv.empty())
            removeVertex();
        else if (!Qe.empty())
            trussPrune();
        else
            break;
    }
}

void KPlexGraph::decDegree(ui u)
{
    if (not pruned[u])
    {
        degree[u]--;
        // heap.decrement(u, 1);
        if (degree[u] == tau_v)
            scheduleRemoval(u);
    }
}
void KPlexGraph::decCN(Edge e)
{
    if (exists(e) and getCN(e)-- == tau_e)
    {
        scheduleRemoval(e);
    }
}

void KPlexGraph::trussPrune()
{
    for (ui i = 0; i < Qe.size(); i++)
    {
        auto e = Qe[i];
        ui u = e.first;
        ui j = e.second;
        ui v = adjList[u][j];
        // if ((not exists(u))  or (not exists(v)))
        if (pruned[u] or pruned[v] or not exists(e))
            continue;
        removeEdge(e);
        decDegree(u);
        decDegree(v);

        Lookup neigh(&lookup, &adjList[u]);
        for (ui j = 0; j < adjList[v].size(); j++)
        {
            ui w = adjList[v][j];
            // if (!exists(w))
            if (pruned[w])
                continue;
            if (neigh[w])
            {
                Edge e1 = getEdge(v, j);
                Edge e2 = getEdge(u, neigh[w] - 1);
                if (exists(e1) and exists(e2))
                {
                    // w is a common neighbor of u, v
                    decCN(e1);
                    decCN(e2);
                }
            }
        }
    }
    re += Qe.size();
    // cout<<"edges removed: "<<Qe.size()<<endl;
    Qe.clear();
}
void KPlexGraph::removeVertex()
{
    ui u = Qv.back();
    Qv.pop_back();

    rv++;
    Lookup neigh(&lookup, &(adjList[u]));
    for (ui i = 0; i < adjList[u].size(); i++)
    {
        Edge e = getEdge(u, i);
        if (!exists(e))
            continue;
        removeEdge(e);
        ui v = adjList[u][i];
        decDegree(v);
        for (ui j = 0; j < adjList[v].size(); j++)
        {
            ui w = adjList[v][j];
            if (!exists(w))
                continue;
            auto e = getEdge(v, j);
            // if(pruned[w] or !exists(e))
            if (neigh[w] and exists(getEdge(u, neigh[w] - 1)))
                // w is a common neighbor of u, v
                decCN(e);
        }
    }
    neigh.erase();
    adjList[u].clear();
    degree[u] = 0;
}



void KPlexGraph::maxDegreeKplex(ui thresh)
{
    Timer t;
    assert(kplex.empty());
    ui *head = new ui[V];
    ui *next = new ui[V];
    ui *degree = new ui[V];

    ui *vis = new ui[V];
    memset(vis, 0, sizeof(ui) * V);

    int max_degree = 0;
    for (ui i = 0; i < V; i++)
        head[i] = V;
    for (ui i = 0; i < V; i++)
    {
        degree[i] = adjList[i].size();
        if (degree[i] > max_degree)
            max_degree = degree[i];
        next[i] = head[degree[i]];
        head[degree[i]] = i;
    }

    for (ui processed_vertices = 0; max_degree + k >= kplex.size() && processed_vertices < thresh; processed_vertices++)
    {
        ui u = V;
        while (max_degree >= 0 && max_degree + k >= kplex.size() && u == V)
        {
            for (ui v = head[max_degree]; v != V;)
            {
                ui tmp = next[v];
                if (degree[v] == max_degree)
                {
                    u = v;
                    head[max_degree] = tmp;
                    break;
                }
                else if (degree[v] + k >= kplex.size())
                {
                    next[v] = head[degree[v]];
                    head[degree[v]] = v;
                }
                v = tmp;
            }
            if (u == V)
            {
                head[max_degree] = V;
                --max_degree;
            }
        }
        if (u == V)
            break;

        vis[u] = 1;
        for (ui v : adjList[u])
            if (!vis[v])
                --degree[v];

        vector<ui> vs;
        for (ui v : adjList[u])
            if (!vis[v])
                vs.push_back(v);

        vector<ui> vs_deg(vs.size());
        for (ui j = 0; j < vs.size(); j++)
            vis[vs[j]] = 2;
        for (ui j = 0; j < vs.size(); j++)
        {
            ui v = vs[j], d = 0;
            for (ui w : adjList[v])
            {
                if (vis[w] == 2)
                    ++d;
            }
            vs_deg[j] = d;
        }
        for (ui j = 0; j < vs.size(); j++)
            vis[vs[j]] = 0;

        vector<ui> res;
        res.push_back(u);
        ui vs_size = vs.size();
        while (vs_size > 0 && res.size() + vs_size + k - 1 > kplex.size())
        {
            // while(vs_size > 0) {
            ui idx = 0;
            for (ui j = 1; j < vs_size; j++)
            {
                if (vs_deg[j] > vs_deg[idx])
                    idx = j;
                else if (vs_deg[j] == vs_deg[idx] && degree[vs[j]] > degree[vs[idx]])
                    idx = j;
            }
            u = vs[idx];

            ui V = 0;
            for (ui v : adjList[u])
                if (!vis[v])
                    vis[v] = 2;
            for (ui j = 0; j < vs_size; j++)
                if (vis[vs[j]])
                {
                    if (j != V)
                        swap(vs[V], vs[j]);
                    vs_deg[V] = vs_deg[j];
                    ++V;
                }
            for (ui v : adjList[u])
                if (vis[v] == 2)
                    vis[v] = 0;

            res.push_back(u);
            for (ui k = 0; k < V; k++)
                vis[vs[k]] = k + 2;
            for (ui j = V; j < vs_size; j++)
            {
                ui v = vs[j];
                for (ui w : adjList[v])
                {
                    if (vis[w] >= 2)
                        --vs_deg[vis[w] - 2];
                }
            }
            for (ui k = 0; k < V; k++)
                vis[vs[k]] = 0;

            vs_size = V;
        }

        // TO DO: extend res to be a maximal k-plex

        if (res.size() > kplex.size())
            kplex = res;
    }

    delete[] vis;
    delete[] head;
    delete[] next;
    delete[] degree;

    printf("*** Heuristic kplex size: %lu, time: %lu (microseconds)\n", kplex.size(), t.elapsed());
}
#endif // KPLEX_GRAPH_H
