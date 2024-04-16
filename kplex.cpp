#include "utils.h"
Timer gtime, check, sr;
#include "kplex-graph.h"
#define PuCSize (P.size() + C.size())
#define PuCuMSize (P.size() + C.size() + M.size())
ui k;
bool flag = false;
#define TIMEOVER (gtime.elapsed() / 1000000 > 60 * 60)
// #define NAIVE
// #define INIT_SEESAW

#define RULE2
#define SEESAW
#define CTCP

#ifdef NAIVE
#define RECSEARCH naiveSearch
#else
#define RECSEARCH recSearch
#endif
class MaxKPlex
{
    KPlexGraph &G;
    KPlexGraph g;
    vecui kplex, bestKPlex;
    ListLinearHeap heap;
    vecui removed;
    vecui dG, dP;
    RandList P, C, M;
    // Partition upper bound vectors...
    MBitSet bmp;
    vector<vecui> PI;
    vecui lookup, ISc, ISp;
    RandList block;
    vecui temp;

public:
    MaxKPlex(KPlexGraph &_G)
        : G(_G), kplex(_G.kplex), g(k, _G.V)
    {
        P.init(G.V);
        C.init(G.V);
        block.init(G.V);
        M.init(G.V);
        dP.resize(G.V);
        dG.resize(G.V);
        removed.resize(G.V);
        heap.init(G);

        PI.resize(G.V);
        ISc.reserve(G.V);
        ISp.reserve(G.V);
        lookup.resize(G.V);
        bmp.init(G.V);
    }
    ui updateC_K(ui u)
    {
        ui sz = C.size();
        for (int i = 0; i < C.size();)
        {
            int ele = C[i];
            if (UNLINK2EQUAL > g.cnMatrix(u, ele))
                removeFromC(ele, true);
            else
                ++i;
        }
        return sz - C.size();
    }
    ui updateSecNeigh_K(ui u)
    {
        ui sz = M.size();
        for (int i = 0; i < M.size();)
        {
            int ele = M[i];
            if (UNLINK2EQUAL > g.cnMatrix(u, ele))
                M.fakeRemove(ele);
            else
                ++i;
        }
        return sz - M.size();
    }

    void search()
    {
        ui v, key;
        G.initCTCP(); // initializes cn and degrees...
        g.lb = G.lb;
        // G.initMin();
        // while (heap.pop_min(v, key)) // returns false when all vertices processed
        // for(ui v=0;v<G.V;v++)
        while (true)
        {
            v = G.popMin();
            // cout<<v<<" ";

            if (v == G.V)
                break;
            if (!G.exists(v))
            {
                continue;
            }
            // g is materialized in createblock with a grpah spanned by v, 1hop and 2hop neighbors
            ui sz1h = createBlk(v);
            //  createBlock(v);
            ui ub = g.degeneracyKPlex();
            if (g.lb > G.lb)
            {
                G.lb = g.lb;
                cout << "Degen found a larger kplex of size: " << g.lb << endl;
            }

            if (ub <= g.lb)
            {
                // cout<<"skipping "<<v<<endl;
                continue;
            }
#ifdef CTCP
            g.applyCTCP(INIT);
            // if seed vertex was removed by CTCP, skip this iteration
            if (g.exists(0) and g.V > g.lb)
                g.shrinkEdges();
#endif
            initContainers(sz1h);
#ifdef INIT_SEESAW
            for (ui i = 0; i < M.size(); i++)
                addToC(M[i]);
            ui ub = seesawUB();
            for (ui i = 0; i < M.size(); i++)
                removeFromC(M[i]);
            if (ub > g.lb)
#endif
            {
                g.buildCommonMatrix(sz1h);
                // g.degenerate();
                sr.tick();
                kSearch(k - 1);
                sr.tock();
            }
            // naiveSearch();

            G.scheduleRemoval(v); // set v to delete...
            if (g.lb > G.lb)
            {
                G.lb = g.lb;
                G.applyCTCP(LBCHANGE); // calculate cn according to new lb...
            }
            else
                G.applyCTCP(NONE);
        }
        cout << "check time: " << check.ticktock() << endl;
        cout << "search time: " << sr.ticktock() << endl;
    }

    void kSearch(ui m)
    {
        if (PuCuMSize <= g.lb or TIMEOVER)
            return;

        if (M.size() == 0)
        {
            flag = false;
            RECSEARCH(FIRST);
            return;
        }
        // first branch
        // i think there is no need to put sec neigh to x, only removing from
        // secneigh sould be sufficient

        // ui u = secNeighToX();
        // kSearch(m);
        // XToSecNeigh(u);

        M.fakePop();
        kSearch(m);
        M.fakeRecPop();
        // vecui rsn(m);
        ui br = 1, rc = 0, rsn = 0;
        while (br < m)
        {
            ui u = M.fakePop();
            rsn++;                      // so that the fakepop here is recoved in the end of this
                                        // function
            addToP_K(u);                // u is added to P, that causes C, M and X to
                                        // shrink... as per theorem 9, 10, 11
            rc += updateC_K(u);         // Applying Theorem 10 to update C
            rsn += updateSecNeigh_K(u); // applying Theorem 9 to update second hop neighbors

            if (M.empty())
            {
                kSearch(m - br);
                break;
            }
            else
            {
                M.fakePop();
                kSearch(m - br);
                M.fakeRecPop();
            }
            br++;
        }
        if (br == m)
        {
            ui u = M.fakePop();
            addToP_K(u);
            rsn++;
            // todo prune 1hop neighbors, and then call recSearch if we can find a
            // larger kplex than ub size
            rc += updateC_K(u);
            // recSearch();
            flag = false;
            RECSEARCH(FIRST);
        }
        while (br >= 1)
        {
            removeFromP_K();
            br--;
        }
        M.fakeRecover(rsn);
        recoverC(rc);
    }
    void recSearch(RecLevel level)
    {
        if ((level == OTHER and flag) or PuCSize <= g.lb or TIMEOVER)
            return;
        ui cp = moveDirectlyToP();
        ui rc = updateC(), ub = 0;
        if (C.empty())
        {
            if (P.size() > g.lb)
            {
                cout << "RecSearch found a larger kplex of size: " << P.size() << endl;
                // g.print();
                g.lb = P.size();
                flag = true;
                // checking validity of kplex
                for (ui i = 0; i < P.size(); i++)
                    if (dP[P[i]] < P.size() - k)
                        cout << " Invalid " << P[i];
            }
            goto RECOVER;
        }
        // ui ub = relaxGCB();
#ifdef SEESAW
        ub = seesawUB();
        if (ub > g.lb)
#endif
        {
            auto B = getBranchings();
            while (B.first < B.second)
            {
                if (level == FIRST)
                    flag = false;
                else if (flag)
                {
                    // return to root level of branchings, but before that recover all branching vertices to C
                    while (B.first++ < B.second)
                        addToC(0, true);
                    break;
                }
                ui bn = maxDegenVertex(B.first++, B.second);
                addToP_K(bn);
                ui rc = pruneC(bn); // apply theorem 11 to remove such vertices in C that can't co-exist with bn
                recSearch(OTHER);
                recoverC(rc);
                removeFromP(bn);
                C.fakeRecPop();
            }
        }
    RECOVER:
        recoverC(rc);
        // recover cp number of vertices directly moved to P
        for (ui i = 0; i < cp; i++)
        {
            PToC(P.top());
        }
        // updateC have done fakeRemove rc vertices, now recover
    }
    void recoverC(ui rc)
    {
        for (ui i = 0; i < rc; i++)
        {
            ui u = C.fakeRecPop();
            for (ui v : g.adjList[u])
                dG[v]++;
        }
    }
    pair<ui, ui> getBranchings()
    {
        for (ui i = 0; i < P.size(); i++)
        {
            ui u = P[i];
            if (support(u) == 0)
                continue;
            // skipping it, because this is a boundary vertex, and it can't have any non-neighbor candidate
            // Lookup neig(&lookup, &g.adjList[u]);
            bmp.setup(g.adjList[u], g.V);
            for (ui j = 0; j < C.size(); j++)
            {
                ui v = C[j];
                if (!bmp.test(v))
                    PI[u].push_back(v);
            }
        }
        ui beta = g.lb - P.size();
        ui cend = 0;
        while (true)
        {
            ui maxpi = -1;
            double maxdise = 0;
            for (ui i = 0; i < P.size(); i++)
            {
                ui u = P[i];
                if (PI[u].empty())
                    continue;
                // cost(pi) = P.size()-dP[pi]
                double cost = min(support(u), (ui)PI[u].size());
                double dise = PI[u].size() / cost;
                if (cost <= beta and dise > maxdise)
                    maxpi = u, maxdise = dise;
            }
            if (maxpi != -1)
            {

                bmp.setup(PI[maxpi], g.V);
                // remove pi* from C
                for (ui i = cend; i < C.size(); i++)
                {
                    ui v = C[i];
                    if (bmp.test(v))
                    {
                        // rather than removing from C, we are changing the positions within C.
                        // When function completes
                        // [0...cend) holds all vertices in C, and [cend, sz) holds the branching vertices.
                        C.swapElements(cend++, i);
                    }
                    // else
                    //     i++;
                }
                // beta-=cost(pi*)
                beta -= min(support(maxpi), (ui)PI[maxpi].size());
                // remove maxpi from every pi
                for (ui i = 0; i < P.size(); i++)
                {
                    // Removing pi* from all pi in PI
                    ui u = P[i];
                    ui j = 0;
                    for (ui k = 0; k < PI[u].size(); k++)
                        if (!bmp.test(PI[u][k]))
                            PI[u][j++] = PI[u][k];
                    PI[u].resize(j);
                }
                // remove maxpi...
                PI[maxpi].clear();
            }
            else
                break;
            if (beta == 0)
                break;
        }

        if (beta > 0)
        {
            cend += min(beta, C.size() - cend);
        }

        ui sz = C.size();
        for (ui i = cend; i < sz; i++)
        {
            // vertices in [cend, sz) range are Branching vertices, just fakeremove them from C
            removeFromC(C.top(), true);
        }

        // clear PI
        for (ui i = 0; i < P.size(); i++)
            PI[P[i]].clear();
        return {cend, sz};
    }
    ui tryPartition()
    {
        for (ui i = 0; i < P.size(); i++)
        {
            ui u = P[i];
            if (support(u) == 0)
                continue;
            // Lookup neig(&lookup, &g.adjList[u]);
            bmp.setup(g.adjList[u], g.V);
            for (ui j = 0; j < C.size(); j++)
            {
                // PI[u] = non-neighbors of u in C
                ui v = C[j];
                if (!bmp.test(v))
                    PI[u].push_back(v);
            }
        }
        ui maxpi = P[0];
        double maxdise = 0;
        ui ub = 0;
        for (ui i = 1; i < P.size(); i++)
        {
            ui u = P[i];
            if (PI[u].empty())
                continue;
            // cost(pi) = P.size()-dP[pi]
            double cost = min(support(u), (ui)PI[u].size());
            double dise = PI[u].size() / cost;
            if (dise > maxdise or (dise == maxdise and PI[u].size() > PI[maxpi].size()))
                maxpi = u, maxdise = dise, ub = cost;
        }

        for (ui u : PI[maxpi])
            ISp.push_back(u);

        // clear PI
        for (ui i = 0; i < P.size(); i++)
            PI[P[i]].clear();
        return ub;
    }
    ui relaxGCB()
    {
        ui UB = P.size();
        ui sz = C.size();
        while (C.size())
        {
            ISc.clear();
            ui ub = tryColor();
            for (ui v : ISc)
                C.fakeRemove(v);
            UB += ub;
        }
        C.fakeRecover(sz);
        return UB;
    }

    ui seesawUB()
    {
        ui UB = P.size();
        ui sz = C.size();
        while (C.size())
        {
            ISp.clear();
            ISc.clear();

            double ubp = tryPartition();
            // ubp = 0;
            // gtime.tick();
            double ubc = tryColor();
            // gtime.tock();
            if (ubp == 0 or
                ISc.size() / ubc > ISp.size() / ubp or
                (ISc.size() / ubc == ISp.size() / ubp and ISc.size() > ISp.size()))

            {
                for (ui v : ISc)
                    C.fakeRemove(v);
                UB += ubc;
            }
            else
            {

                for (ui v : ISp)
                    C.fakeRemove(v);
                UB += ubp;
            }
            // cout<<C.size()<<" "<<ISp.size();
        }
        C.fakeRecover(sz);
        return UB;
    }
    void createIS()
    {
        ISc.push_back(C[0]);

        // check.tick();
        for (ui i = 1; i < C.size(); i++)
        {
            ui u = C[i];
            bool flag = true;
            // for(ui v:g.adjList[u])
            //     bm.set(v);
            // Lookup l1(&lookup, &g.adjList[u], true);
            // Lookup neigh(&lookup, &g.adjList[u], true);
            bmp.setup(g.adjList[u], g.V);
            for (ui v : ISc)
                // if (bmp.test(v))
                // if(l1[v])
                // if(isNeighbor(g.adjList[u], v))
                if (bmp.test(v))
                {
                    flag = false;
                    break;
                }
            if (flag)
                ISc.push_back(u);
        }
        // check.tock();
    }
    ui support(ui u)
    {
        return k - (P.size() - dP[u]);
    }

    ui TISUB()
    {
        ui maxsup = 0;
        for (ui i = 0; i < ISc.size(); i++)
        {
            for (ui j = i + 1; j < ISc.size(); j++)
            {
                if (support(ISc[j]) > support(ISc[i]))
                    swap(ISc[i], ISc[j]);
            }
            // if (i+1 > support(ISc[i])){
            //     return i;
            // }
            // not using <= condition because i is starting from 0...
            if (i < support(ISc[i]))
                maxsup++;
            else
                break;
        }
        return maxsup;
    }
    ui tryColor()
    {
        createIS();
        ui ub = TISUB();
        ui vlc = 0;
        // collect loose vertices i.e. v \in ISc | support(v) > ub
        for (ui i = 0; i < ISc.size(); i++)
        {
            if (support(ISc[i]) > ub)
            {
                if (i != vlc)
                    swap(ISc[i], ISc[vlc]);
                vlc++;
            }
        }
        // ISc[0... vlc) we have loose vertices

        // Lookup inIS(&lookup, &ISc, true);
        bmp.setup(ISc, g.V);
        for (ui i = 0; vlc < ub and i < C.size(); i++)
        {
            ui u = C[i];
            if (bmp.test(u)) // this loop running for C\ISc
                continue;
            ui vc = 0;
            for (ui j = vlc; j < ISc.size(); j++) // this loop runs in ISc\LC
            {
                ui v = ISc[j];
                if (isNeighbor(g.adjList[v], u))
                {
                    if (j != vlc + vc)
                        swap(ISc[vlc + vc], ISc[j]);
                    vc++;
                }
            }
            if (vlc + vc + 1 <= ub)
            {
                vlc += vc;
                ISc.push_back(u);
                swap(ISc.back(), ISc[vlc++]);
                bmp.set(u);
            }
        }
        for (ui i = 0; i < C.size(); i++)
        {
            ui u = C[i];
            if (bmp.test(u) or support(u) >= ub) // this loop running for C\ISc
                continue;
            ui nv = 0;
            for (ui v : g.adjList[u])
            {
                if (bmp.test(v))
                    nv++;
            }
            if (nv <= ub - support(u))
            {
                ISc.push_back(u);
                bmp.set(u);
            }
        }
        return ub;
    }
    // bool upperboundK()
    // {
    //     int cnt = P.size();
    //     for (int i = 0; i < P.size(); i++)
    //     {
    //         const int v = P[i];
    //         nonadjInP[i] = P.size() - dP[v];
    //     }
    //     for (int i = 0; i < C.size(); i++)
    //     {
    //         int max_noadj = -1;
    //         int max_index = 0;
    //         const int ele = C[i];
    //         for (int j = 1; j < P.size(); j++)
    //         {
    //             const int v = P[j];
    //             if (!isAdjMtx(ele, v) && nonadjInP[j] > max_noadj)
    //             {
    //                 max_noadj = nonadjInP[j];
    //                 max_index = j;
    //             }
    //         }
    //         if (max_noadj < k)
    //         {
    //             cnt++;
    //             nonadjInP[max_index]++;
    //         }
    //         if (cnt >= G.lb)
    //             return true;
    //     }
    //     if (cnt >= G.lb)
    //         return true;
    //     return false;
    // }

    ui moveDirectlyToP()
    {
        ui sz = 0;
        vecui temp = P.getData();
        for (ui i = 0; i < C.size();)
        {
            ui u = C[i];
            if (!canMoveToP(u))
            {
                i++;
                continue;
            }
            // vldb BR1 rule
            if (dG[u] + 2 >= PuCSize)
            {
                CToP(u);
                sz++;
            }
            else
#ifdef RULE2
            {
                // xiao2017 Lemma4 rule. u and all non-nieghbors vertices of u shoulud be k-satisfied
                //  u is already k-satisfied. we check here non-neighbors be k-satisfied.
                bool flag = true;
                bmp.setup(g.adjList[u], g.V);
                for (ui i = 0; flag and i < P.size(); i++)
                {
                    ui v = P[i];
                    if (!bmp.test(v) and dG[v] + k < PuCSize)
                        flag = false;
                }

                for (ui i = 0; flag and i < C.size(); i++)
                {
                    ui v = C[i];
                    if (!bmp.test(v) and dG[v] + k < PuCSize)
                        flag = false;
                }
                if (flag)
                {
                    sz++;
                    CToP(u);
                }
                else
                    i++;
            }
#else
                i++;
#endif
        }
        return sz;
    }
    ui pruneC(ui u)
    {
        // apply theorem 11 to remove those vertices in C that can't exist with u.
        ui sz = C.size();
        for (int i = 0; i < C.size();)
        {
            int ele = C[i];
            if (UNLINK2EQUAL > g.cnMatrix(u, ele) or !canMoveToP(ele))
                removeFromC(ele, true); // fake remove when flag is true
            else
                ++i;
        }
        return sz - C.size();
    }
    ui maxDegenVertex(ui st, ui en)
    {
        // using C.get because this function can return even those vertices which are fake-removed from C
        ui bn = C.get(st), ind = st;
        for (ui i = st + 1; i < en; i++)
        {
            ui v = C.get(i);
            if (g.peelSeq[v] > g.peelSeq[bn])
                // if(g.adjList[v].size()>g.adjList[bn].size())
                bn = v, ind = i;
        }
        C.swapElements(ind, C.size());
        return bn;
    }
    ui createBlk(ui v)
    {
        reset();
        block.add(v);
        ui vid = 0;
        for (ui i = 0; i < G.adjList[v].size(); i++)
        {
            ui u = G.adjList[v][i];
            // cout<<u<<" "<<G.V<<endl;
            if (G.exists(u) and G.exists(G.getEdge(v, i)))
            {
                block.add(u);
                g.adjList[0].push_back(++vid);
            }
        }
        ui sz1h = vid;
        temp.reserve(sz1h);
        temp.clear();
        // 0 is seed vertex ID, [1...sz1h] is 1hop neighbors, then 2hop neighbors
        for (ui i = 1; i <= sz1h; i++)
        {
            ui u = block[i];
            g.adjList[i].push_back(0);
            for (ui j = 0; j < G.adjList[u].size(); j++)
            {
                ui w = G.adjList[u][j];
                if (w == v or !G.exists(w) or !G.exists(G.getEdge(u, j)))
                    continue;

                if (not block.contains(w))
                    // it's a newly discovered two hope neighbor, so add it to block
                    block.add(w);
                vid = block.getIndex(w);
                // it's a two-hop neighbor, put them in temp and insert in the end after sorting
                if (vid > sz1h)
                    temp.push_back(vid);
                else
                    g.adjList[i].push_back(vid);
            }
            sort(temp.begin(), temp.end());
            g.adjList[i].insert(g.adjList[i].end(), temp.begin(), temp.end());
            temp.clear();
        }

        for (ui i = sz1h + 1; i < block.size(); i++)
        {
            ui w = block[i];
            for (ui j = 0; j < G.adjList[w].size(); j++)
            {
                ui z = G.adjList[w][j];
                // z and w both are two-hop neighbors
                if (block.contains(z) and G.exists(G.getEdge(w, j)))
                {
                    vid = block.getIndex(z);
                    // it's a two-hop neighbor, put them in temp and insert in the end after sorting
                    if (vid > sz1h)
                        temp.push_back(vid);
                    else
                        g.adjList[i].push_back(vid);
                }
            }
            sort(temp.begin(), temp.end());
            g.adjList[i].insert(g.adjList[i].end(), temp.begin(), temp.end());
            temp.clear();
        }
        g.V = block.size();
        // g.sort();
        // block.clear();
        // g.print();

        return sz1h;
    }

    ui createBlock(ui v)
    {
        reset();
        block.add(v);
        ui vid = g.addVertex();

        for (ui i = 0; i < G.adjList[v].size(); i++)
        {
            ui u = G.adjList[v][i];
            // cout<<u<<" "<<G.V<<endl;
            if (G.exists(u) and G.exists(G.getEdge(v, i)))
            {
                block.add(u);
                vid = g.addVertex();
                g.addEdge(0, vid);
            }
        }
        ui sz1h = vid;
        // 0 is seed vertex ID, [1...sz1h] is 1hop neighbors, then 2hop neighbors
        for (ui i = 1; i <= sz1h; i++)
        {
            ui u = block[i];
            for (ui j = 0; j < G.adjList[u].size(); j++)
            {
                ui w = G.adjList[u][j];
                if (w == v or !G.exists(w) or !G.exists(G.getEdge(u, j)))
                    continue;

                if (block.contains(w))
                {
                    ui ind = block.getIndex(w);
                    // it's a one hope neighbor
                    if (ind > i and ind <= sz1h)
                        g.addEdge(i, ind);
                    // else this edge is already added by i
                }
                else
                {
                    // it's a 2 hope neighbor
                    block.add(w);
                    vid = g.addVertex();
                    // g.addEdge(i, vid);
                }
            }
        }

        for (ui i = sz1h + 1; i < block.size(); i++)
        {
            ui w = block[i];
            for (ui j = 0; j < G.adjList[w].size(); j++)
            {
                ui z = G.adjList[w][j];
                // z and w both are two-hop neighbors
                if (G.exists(G.getEdge(w, j)))
                    if (block.contains(z) and block.getIndex(z) < i)
                        g.addEdge(i, block.getIndex(z));
            }
        }
        g.sort();
        block.clear();
        // g.print();

        return sz1h;
    }
    void initContainers(ui sz1h)
    {
        addToP_K(0);
        for (ui u = 1; u < g.V; u++)
        {
#ifdef CTCP
            if (g.exists(u))
#endif
                if (u <= sz1h)
                    addToC(u);
                else
                    M.add(u);
        }
    }

    void addToC(ui u, bool fake = false)
    {
        if (fake)
            u = C.fakeRecPop();
        else
            C.add(u);
        for (ui v : g.adjList[u])
            dG[v]++;
    }

    void removeFromC(ui u, bool fake = false)
    {
        if (fake)
            C.fakeRemove(u);
        else
            C.remove(u);
        for (ui v : g.adjList[u])
            dG[v]--;
    }

    void addToP(ui u)
    {
        P.add(u);
        for (ui v : g.adjList[u])
            dP[v]++;
    }
    void removeFromP(ui u)
    {
        P.remove(u);
        for (ui v : g.adjList[u])
            dP[v]--;
    }

    ui addToP_K(ui u)
    {
        addToP(u);
        for (ui v : g.adjList[u])
            dG[v]++;
        return u;
    }
    ui removeFromP_K()
    {
        ui u = P.top();
        removeFromP(u);
        // M.add(u);
        for (ui v : g.adjList[u])
            dG[v]--;
        return u;
    }
    ui updateC()
    {
        ui sz = C.size();
        for (ui i = 0; i < C.size();)
        {
            ui u = C[i];
            if (!canMoveToP(u))
                removeFromC(u, true);
            else
                i++;
        }

        return sz - C.size();
    }
    bool canMoveToP(ui u)
    {
        // u is not yet in P, hence checking <=
        if (dP[u] + k <= P.size())
            return false;
        for (ui i = 0; i < P.size(); i++)
        {
            ui v = P.get(i);
            if (dP[v] + k == P.size() && !isNeighbor(g.adjList[u], v))
                return false;
        }
        return true;
    }
    void CToP(ui u)
    {
        C.remove(u);
        addToP(u);
    }
    void PToC(ui u)
    {
        removeFromP(u);
        C.add(u);
    }
    void reset()
    {
        fill(dP.begin(), dP.begin() + g.V, 0);
        fill(dG.begin(), dG.begin() + g.V, 0);
        M.clear();
        block.clear();
        P.clear();
        C.clear();
        g.clear();
    }
    void naiveSearch()
    {
        if (PuCSize < g.lb)
            return;
        ui rc = updateC();
        if (C.empty())
        {
            if (P.size() > g.lb)
            {
                cout << P.size() << " kp:";
                P.print();
                g.lb = P.size();
                for (ui i = 0; i < P.size(); i++)
                    if (dP[P[i]] < P.size() - k)
                    {
                        cout << " Invalid " << P[i];
                    }
                cout << endl;
            }
            C.fakeRecover(rc);
            return;
        }
        ui u = C.top();
        // cout<<u<<" c ";C.print();
        CToP(u);
        naiveSearch();
        PToC(u);
        C.fakePop();
        naiveSearch();
        C.fakeRecPop();
        C.fakeRecover(rc);
    }
};

int main(int argc, char *argv[])
{
    CommandLine cmd(argc, argv);
    std::string file = cmd.GetOptionValue("-g", "");
    k = cmd.GetOptionIntValue("-k", 2);

    if (file == "")
    {
        cout << "Please provide data file" << endl;
        exit(-1);
    }
    KPlexGraph Gr(k, file);
    // Gr.maxDegreeKP(10);
    ui ub = Gr.degeneracyKPlex();
    ui n = Gr.V, m = Gr.E / 2, initkpx = Gr.lb;
    cout << "Init: n = " << n << " m = " << m << endl;
    cout << "InitKPX: " << initkpx << " UB: " << ub << endl;
    Timer init, search;
    if (Gr.lb < ub)
    {
        init.tick();
        Gr.corePrune(Gr.lb + 1 - k);
        cout << "After core-shrink n = " << Gr.V << " , m = " << Gr.E << endl;
        Gr.applyCTCP(INIT);
        Gr.shrink(true);
        cout << "After core-truss pruning n = " << Gr.V << " , m = " << Gr.E << endl;
        // when true passed, it will also check existency of edges. because some edges are removed by ctcp
        init.tock();

        if (Gr.V > Gr.lb)
        {
            MaxKPlex kps(Gr);
            search.tick();
            kps.search();
            search.tock();
        }
    }

    cout << ">>" << file;
    cout << " n " << n << " , m " << m;
    cout << " K " << k;
    cout << " MaxKPX " << Gr.lb;
    cout << " InitKPX " << initkpx;
    cout << " Init_Time " << init.ticktock();
    cout << " Search_Time " << search.ticktock();
    cout << " Total_Time " << init.ticktock() + search.ticktock();
    cout << endl;
    return 0;
}