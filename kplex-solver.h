#include "kplex-graph.h"
Timer gtime, check, sr;
#define PuCSize (P.size() + C.size())
#define PuCuMSize (P.size() + C.size() + M.size())
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
    KPlexGraph g;
    vecui &kplex;
    vecui dG, dP;
    RandList P, C, M;
    // Partition upper bound vectors...
    MBitSet bmp;
    vector<vecui> PI;
    vecui lookup, ISc, ISp;
    RandList block;
    vecui temp;
    vecui matrix;
    ui n, sz1h, k;

public:

    MaxKPlex(ui m, ui _k, vecui& kp)
        :  g(m), kplex(kp), k(_k)
    {
        P.init(m);
        C.init(m);
        block.init(m);
        M.init(m);
        bmp.init(m);

        dP.resize(m);
        dG.resize(m);
        PI.resize(m);
        lookup.resize(m);
        
        kplex.reserve(m);
        ISc.reserve(m);
        ISp.reserve(m);
    }
    void solve_instance(ui n, ui sz1h, const auto& vp){
        g.load(vp, n);
        g.buildCommonMatrix(sz1h);
        initContainers(sz1h);

        kSearch(k-1);

        g.unload(vp, n);
        reset();
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

    void kSearch(ui m)
    {
        if (PuCuMSize <= kplex.size() or TIMEOVER)
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
        if ((level == OTHER and flag) or PuCSize <= kplex.size() or TIMEOVER)
            return;
        ui cp = moveDirectlyToP();
        ui rc = updateC(), ub = 0;
        if (C.empty())
        {
            if (P.size() > kplex.size())
            {
                kplex.clear();
                cout << "RecSearch found a larger kplex of size: " << P.size() << endl;
                flag = true;
                // checking validity of kplex
                for (ui i = 0; i < P.size(); i++){
                    ui u = P[i];
                    if (dP[u] < P.size() - k)
                        cout << " Invalid " << u;
                    kplex.push_back(u);
                }
            }
            goto RECOVER;
        }
        // ui ub = relaxGCB();
#ifdef SEESAW
        ub = seesawUB();
        if (ub > kplex.size())
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
        ui beta = kplex.size() - P.size();
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
    //         if (cnt >= kplex.size())
    //             return true;
    //     }
    //     if (cnt >= kplex.size())
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
    

    void initContainers(ui sz1h)
    {
        addToP_K(0);
        for (ui u = 1; u < g.V; u++)
        {
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
        if (PuCSize < kplex.size())
            return;
        ui rc = updateC();
        if (C.empty())
        {
            if (P.size() > kplex.size())
            {
                kplex.clear();
                cout << P.size() << " kp:";
                P.print();
                for (ui i = 0; i < P.size(); i++)
                    if (dP[P[i]] < P.size() - k)
                        cout << " Invalid " << P[i];
                    else
                        kplex.push_back(P[i]);
                    
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
