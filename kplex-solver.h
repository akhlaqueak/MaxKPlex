#include "kplex-graph.h"
Timer gtime, check, sr;
#define PuCSize (P.size() + C.size())
#define PuCuMSize (P.size() + C.size() + M.size())
bool flag = false;
#define TIMEOVER (gtime.elapsed() / 1000000 > 60 * 60)
// #define NAIVE
// #define INIT_SEESAW

// #define RULE2
// #define SEESAW
// #define CTCP
// #defeine _SECOND_ORDER_PRUNING_

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
    ui sz1h;

    // BBMATRIX Stuff

    ui n;
    char *matrix;
    long long matrix_size;
#ifdef _SECOND_ORDER_PRUNING_
    ui *cn;
    std::queue<std::pair<ui, ui>> Qe;
    std::vector<std::pair<ui, ui>> removed_edges;
    long long removed_edges_n;
#endif

    ui *degree;
    ui *degree_in_S;

    ui K;
    vecui best_solution;
    ui best_solution_size;

    ui *neighbors;
    ui *nonneighbors;

    ui *SR;     // union of S and R, where S is at the front
    ui *SR_rid; // reverse ID for SR
    std::queue<ui> Qv;
    ui *level_id;

    std::vector<ui> must_include_vertices;
    ui *peelOrder;
    ui R_end;

public:
    MaxKPlex(ui m, ui _k, vecui &kp)
        : g(_k, m), kplex(kp), K(_k)
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

        matrix = new char[m * m];
        neighbors = new ui[m];
        nonneighbors = new ui[m];
        degree = new ui[m];
        SR = new ui[m];
        SR_rid = new ui[m];
        level_id = new ui[m];
        peelOrder = new ui[m];
        best_solution_size = 0;
        R_end=0;
    }

    void initialization(const auto &vp, bool must_include_0)
    {
#ifdef _SECOND_ORDER_PRUNING_
        delete[] cn;
        cn = new ui[matrix_size];
#endif

#ifdef _SECOND_ORDER_PRUNING_
        memset(cn, 0, sizeof(ui) * ((long long)n) * n);
#endif
        // memset(matrix, 0, sizeof(char) * ((long long)n) * n);
        // fill(matrix, matrix+(n*n), 0);

        fill(degree, degree+n, 0);
        for (auto & e: vp)
        {
            assert(e.first >= 0 && e.first < n && e.second >= 0 && e.second < n);
            ui a = e.first, b = e.second;
            degree[a]++;
            degree[b]++;
            matrix[a * n + b] = matrix[b * n + a] = 1;
        }
        // the following computes a degeneracy ordering and a heuristic solution
        ui *peel_sequence = neighbors;
        ui *core = nonneighbors;
        ui *vis = SR;
        memset(vis, 0, sizeof(ui) * n);
        ui max_core = 0, UB = 0, idx = n;
        for (ui i = 0; i < n; i++)
        {
            // todo check if heap works better here.
            ui u, min_degree = n;
            for (ui j = 0; j < n; j++)
                if (!vis[j] && degree[j] < min_degree)
                {
                    u = j;
                    min_degree = degree[j];
                }
            if (min_degree > max_core)
                max_core = min_degree;
            core[u] = max_core;
            peel_sequence[i] = u;
            peelOrder[u] = i;
            vis[u] = 1;

            ui t_UB = core[u] + K;
            if (n - i < t_UB)
                t_UB = n - i;
            if (t_UB > UB)
                UB = t_UB;

            if (idx == n && min_degree + K >= n - i)
                idx = i;

            for (ui j = 0; j < n; j++)
                if (!vis[j] && matrix[u * n + j])
                    --degree[j];
        }
        if (n - idx > best_solution_size)
        {
            best_solution_size = n - idx;
            best_solution.clear();
            for (ui i = idx; i < n; i++)
                // best_solution[i-idx] = peel_sequence[i];
                best_solution.push_back(peel_sequence[i]);
            printf("Degen find a solution of size %u\n", best_solution_size);
        }

        // memset(degree_in_S, 0, sizeof(ui)*n);
        R_end = 0;
        for (ui i = 0; i < n; i++)
            SR_rid[i] = n;
        for (ui i = 0; i < n; i++)
            if (core[i] + K > best_solution_size)
            {
                SR[R_end] = i;
                SR_rid[i] = R_end;
                ++R_end;
            }

        if (must_include_0 && SR_rid[0] == n)
        {
            R_end = 0;
            return;
        }

        for (ui i = 0; i < R_end; i++)
        {
            ui u = SR[i];
            degree[u] = 0;
            for (ui j = 0; j < R_end; j++)
                if (matrix[u * n + SR[j]])
                    ++degree[u];
        }

        memset(level_id, 0, sizeof(ui) * n);
        for (ui i = 0; i < R_end; i++)
            level_id[SR[i]] = n;

        assert(Qv.empty());

#ifdef _SECOND_ORDER_PRUNING_
        for (ui i = 0; i < R_end; i++)
        {
            ui neighbors_n = 0;
            char *t_matrix = matrix + SR[i] * n;
            for (ui j = 0; j < R_end; j++)
                if (t_matrix[SR[j]])
                    neighbors[neighbors_n++] = SR[j];
            for (ui j = 0; j < neighbors_n; j++)
                for (ui k = j + 1; k < neighbors_n; k++)
                {
                    ++cn[neighbors[j] * n + neighbors[k]];
                    ++cn[neighbors[k] * n + neighbors[j]];
                }
        }

        while (!Qe.empty())
            Qe.pop();
        for (ui i = 0; i < R_end; i++)
            for (ui j = i + 1; j < R_end; j++)
            {
                if (matrix[SR[i] * n + SR[j]] && upper_bound_based_prune(0, SR[i], SR[j]))
                {
                    Qe.push(std::make_pair(SR[i], SR[j]));
                }
            }
        removed_edges_n = 0;
#endif

        must_include_vertices.resize(n);

        if (remove_vertices_and_edges_with_prune(0, R_end, 0))
            R_end = 0;
    }

    bool
    remove_vertices_and_edges_with_prune(ui S_end, ui &R_end, ui level)
    {
#ifdef _SECOND_ORDER_PRUNING_
        while (!Qv.empty() || !Qe.empty())
        {
#else
        while (!Qv.empty())
        {
#endif
            while (!Qv.empty())
            {
                ui u = Qv.front();
                Qv.pop(); // remove u
                assert(SR[SR_rid[u]] == u);
                assert(SR_rid[u] >= S_end && SR_rid[u] < R_end);
                --R_end;
                swap_pos(SR_rid[u], R_end);

                bool terminate = false;
                ui neighbors_n = 0;
                char *t_matrix = matrix + u * n;
                for (ui i = 0; i < R_end; i++)
                    if (t_matrix[SR[i]])
                    {
                        ui w = SR[i];
                        neighbors[neighbors_n++] = w;
                        --degree[w];
                        if (degree[w] + K <= best_solution_size)
                        {
                            if (i < S_end)
                                terminate = true; // UB1
                            else if (level_id[w] > level)
                            { // RR3
                                level_id[w] = level;
                                Qv.push(w);
                            }
                        }
                    }
                // UB1
                if (terminate)
                {
                    for (ui i = 0; i < neighbors_n; i++)
                        ++degree[neighbors[i]];
                    level_id[u] = n;
                    ++R_end;
                    return true;
                }

#ifdef _SECOND_ORDER_PRUNING_
                for (ui i = 1; i < neighbors_n; i++)
                {
                    ui w = neighbors[i];
                    for (ui j = 0; j < i; j++)
                    {
                        ui v = neighbors[j];
                        assert(cn[v * n + w]);
#ifndef NDEBUG
                        ui common_neighbors = 0;
                        for (ui k = S_end; k <= R_end; k++)
                            if (matrix[SR[k] * n + v] && matrix[SR[k] * n + w])
                                ++common_neighbors;
                        assert(cn[v * n + w] == common_neighbors);
                        assert(cn[w * n + v] == common_neighbors);
#endif
                        --cn[v * n + w];
                        --cn[w * n + v];
#ifndef NDEBUG
                        common_neighbors = 0;
                        for (ui k = S_end; k < R_end; k++)
                            if (matrix[SR[k] * n + v] && matrix[SR[k] * n + w])
                                ++common_neighbors;
                        assert(cn[v * n + w] == common_neighbors);
                        assert(cn[w * n + v] == common_neighbors);
#endif

                        if (!upper_bound_based_prune(S_end, v, w))
                            continue;

                        if (SR_rid[w] < S_end)
                            terminate = true; // v, w \in S --- UB2
                        else if (SR_rid[v] >= S_end)
                        { // v, w, \in R --- RR5
                            if (matrix[v * n + w])
                                Qe.push(std::make_pair(v, w));
                        }
                        else if (level_id[w] > level)
                        { // RR4
                            level_id[w] = level;
                            Qv.push(w);
                        }
                    }
                }
                if (terminate)
                    return true;
#endif
            }

#ifdef _SECOND_ORDER_PRUNING_
            if (Qe.empty())
                break;

            ui v = Qe.front().first, w = Qe.front().second;
            Qe.pop();
            if (level_id[v] <= level || level_id[w] <= level || !matrix[v * n + w])
                continue;
            assert(SR_rid[v] >= S_end && SR_rid[v] < R_end && SR_rid[w] >= S_end && SR_rid[w] < R_end);

            if (degree[v] + K <= best_solution_size + 1)
            {
                level_id[v] = level;
                Qv.push(v);
            }
            if (degree[w] + K <= best_solution_size + 1)
            {
                level_id[w] = level;
                Qv.push(w);
            }
            if (!Qv.empty())
                continue;

#ifndef NDEBUG
                // printf("remove edge between %u and %u\n", v, w);
#endif

            assert(matrix[v * n + w]);
            matrix[v * n + w] = matrix[w * n + v] = 0;
            --degree[v];
            --degree[w];

            if (removed_edges.size() == removed_edges_n)
            {
                removed_edges.push_back(std::make_pair(v, w));
                ++removed_edges_n;
            }
            else
                removed_edges[removed_edges_n++] = std::make_pair(v, w);

            char *t_matrix = matrix + v * n;
            for (ui i = 0; i < R_end; i++)
                if (t_matrix[SR[i]])
                {
                    --cn[w * n + SR[i]];
                    --cn[SR[i] * n + w];
                    if (!upper_bound_based_prune(S_end, w, SR[i]))
                        continue;
                    if (i < S_end)
                    {
                        if (level_id[w] > level)
                        {
                            level_id[w] = level;
                            Qv.push(w);
                        }
                    }
                    else if (matrix[w * n + SR[i]])
                        Qe.push(std::make_pair(w, SR[i]));
                }
            t_matrix = matrix + w * n;
            for (ui i = 0; i < R_end; i++)
                if (t_matrix[SR[i]])
                {
                    --cn[v * n + SR[i]];
                    --cn[SR[i] * n + v];
                    if (!upper_bound_based_prune(S_end, v, SR[i]))
                        continue;
                    if (i < S_end)
                    {
                        if (level_id[v] > level)
                        {
                            level_id[v] = level;
                            Qv.push(v);
                        }
                    }
                    else if (matrix[v * n + SR[i]])
                        Qe.push(std::make_pair(v, SR[i]));
                }
#endif
        }

        return false;
    }
    void swap_pos(ui i, ui j)
    {
        std::swap(SR[i], SR[j]);
        SR_rid[SR[i]] = i;
        SR_rid[SR[j]] = j;
    }

    void solve_instance(ui _n, ui sz1h, const auto &vp)
    {
        n = _n;
        cout<<n<<" "<<_n;
        initialization(vp, true);
        if (R_end)
        {
            initContainers(sz1h);
            kSearch(K - 1);
        }
        reset(vp);
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
        rc += updateC_SecondOrder();
        if (C.empty())
        {
            if (P.size() > kplex.size())
            {
                kplex.clear();
                cout << "RecSearch found a larger kplex of size: " << P.size() << endl;
                flag = true;
                // checking validity of kplex
                for (ui i = 0; i < P.size(); i++)
                {
                    ui u = P[i];
                    if (dP[u] < P.size() - K)
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
            if (secondOrderUB())
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
    ui updateC_SecondOrder()
    {
        ui sz = C.size();
        for (ui j = 0; j < P.size(); j++)
            for (ui i = 0; i < C.size();)
            {
                ui u = C[i], v = P[j];
                if (P.size() + 1 + support(u) + support(v) + g.soMatrix(u, v) <= kplex.size())
                    removeFromC(u, true);
                else
                    i++;
            }

        return sz - C.size();
    }
    ui secondOrderUB()
    {
        ui ub = PuCSize;
        for (ui i = 0; i < P.size(); i++)
        {
            for (ui j = 0; j < P.size(); j++)
                if (i != j)
                {
                    ui u = P[i], v = P[j];
                    ui b = support(u) + support(v) + g.soMatrix(u, v) + P.size();
                    if (b < ub)
                        ub = b;
                }
        }
        return ub > kplex.size();
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
        return K - (P.size() - dP[u]);
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
    //         if (max_noadj < K)
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
                    if (!bmp.test(v) and dG[v] + K < PuCSize)
                        flag = false;
                }

                for (ui i = 0; flag and i < C.size(); i++)
                {
                    ui v = C[i];
                    if (!bmp.test(v) and dG[v] + K < PuCSize)
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
            if (g.core[v] > g.core[bn])
                // if(g.adjList[v].size()>g.adjList[bn].size())
                bn = v, ind = i;
        }
        C.swapElements(ind, C.size());
        return bn;
    }

    void initContainers(ui sz1h)
    {
        for (ui i = 0; i < R_end; i++)
        {
            for (ui j = 0; j < R_end; j++)
            {
                ui u = SR[i], v = SR[j];
                if (matrix[u * n + v])
                    g.adjList[i].push_back(j);
            }
        }
        g.V = R_end;

        addToP_K(0);
        for (ui i = 1; i < R_end; i++)
        {
            ui u = SR[i];
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
        {
            dP[v]++;
            g.soMatrix(u, v)--;
        }
    }
    void removeFromP(ui u)
    {
        P.remove(u);
        for (ui v : g.adjList[u])
        {
            dP[v]--;
            g.soMatrix(u, v)++;
        }
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
        if (dP[u] + K <= P.size())
            return false;
        for (ui i = 0; i < P.size(); i++)
        {
            ui v = P.get(i);
            if (dP[v] + K == P.size() && !isNeighbor(g.adjList[u], v))
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
    void reset(const auto& vp)
    {
        for (auto &e : vp)
            matrix[e.first * n + e.second] = matrix[e.second * n + e.first] = 0;
        fill(dP.begin(), dP.begin() + n, 0);
        fill(dG.begin(), dG.begin() + n, 0);
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
                    if (dP[P[i]] < P.size() - K)
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
