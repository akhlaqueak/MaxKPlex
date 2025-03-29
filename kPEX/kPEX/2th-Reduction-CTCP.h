#ifndef REDUCTION_H
#define REDUCTION_H

#include "Graph.h"

/**
 * @brief second order reduction
 */
class Reduction
{
private:
    ui *triangles;
    Graph *g;
    ui n, m;
    vector<bool> vis; // a cache when counting triangles
    ui *edge_to, *pstart, *deg;
    /**
     * @brief count the number of triangles for each edge and save the results in triangls[]
     * T(n)=O(m + sigma( min(d[from], d[to]) | (from,to) in E) )
     */
    void countTriangles()
    {
        if (triangles == nullptr)
        {
            triangles = new ui[m];
        }
        double start_list = get_system_time_microsecond();
        for (ui i = 0; i < n; i++)
        {
            for (ui j = pstart[i]; j < pstart[i + 1]; j++)
            {
                assert(j < m);
                assert(edge_to[j] < n);
                vis[edge_to[j]] = 1;
            }
            for (ui j = pstart[i]; j < pstart[i + 1]; j++)
            {
                ui v = edge_to[j];
                if (g->d[i] < g->d[v])
                    continue;
                int cnt = 0;
                for (ui a = pstart[v]; a < pstart[v + 1]; a++)
                {
                    ui w = edge_to[a];
                    if (vis[w])
                        cnt++;
                }
                triangles[j] = cnt;
            }
            for (ui j = pstart[i]; j < pstart[i + 1]; j++)
                vis[edge_to[j]] = 0;
        }
        for (ui i = 0; i < n; i++)
        {
            for (ui j = pstart[i]; j < pstart[i + 1]; j++)
            {
                ui v = edge_to[j];
                if (g->d[i] > g->d[v])
                    continue;
                ui ano_edge = lower_bound(edge_to + pstart[v], edge_to + pstart[v + 1], i) - edge_to;
                triangles[j] = triangles[ano_edge];
            }
        }
        list_triangle_time += get_system_time_microsecond() - start_list;
    }
    /**
     * @brief weak reduce
     *
     * @param q containg the vertices that we need to remove; Initially, q is not empty
     * @param vertex_removed a mask to record whether a vertex is in queue
     * @param edge_removed a mask recording whether an edge is removed by strong reduction
     *
     * @return the number of vertices that we remove
     */
    int remove_vertex(queue<ui> &q, vector<bool> &vertex_removed, vector<bool> &edge_removed)
    {
        int cnt = 0;
        while (q.size())
        {
            ui u = q.front();
            q.pop();
            cnt++;
            for (ui i = pstart[u]; i < pstart[u + 1]; i++)
            {
                ui v = edge_to[i];
                if (vertex_removed[v] || edge_removed[i])
                    continue;
                if (--deg[v] + paramK <= lb)
                    q.push(v), vertex_removed[v] = 1;
            }
        }
        return cnt;
    }

public:
    Reduction(Graph *_g) : g(_g), triangles(nullptr)
    {
        n = g->n;
        m = g->m;
        vis.resize(n, 0);
    }
    ~Reduction()
    {
        if (triangles != nullptr)
        {
            delete[] triangles;
            triangles = nullptr;
        }
    }
    /**
     * @brief we first remove edges with updating triangles[]; also, we put core-reduce on highest priority.
     * assume we already conduct weak reduce
     */
    void strong_reduce(int lb, int start_u = 0)
    {
        CTCP(lb);
    }
    void rebuild_graph(vector<bool> &vertex_removed, vector<bool> &edge_removed, ui *que = nullptr)
    {
        ui *q = que;
        if (q == nullptr)
        {
            q = new ui[n];
        }
        ui new_n = 0;
        assert(n == g->map_refresh_id.size());
        vector<ui> new_map(n);
        for (ui i = 0; i < n; i++)
            if (!vertex_removed[i])
            {
                new_map[new_n] = g->map_refresh_id[i];
                q[i] = new_n++;
            }
        new_map.resize(new_n);
        g->map_refresh_id = new_map;
        ui *new_pstart = new ui[new_n + 1];
        ui *new_d = new ui[new_n];
        ui j = 0;
        for (ui i = 0; i < n; i++)
        {
            if (vertex_removed[i])
                continue;
            ui u = q[i];
            new_pstart[u] = j;
            for (ui p = pstart[i]; p < pstart[i + 1]; p++)
            {
                ui v = edge_to[p];
                if (!vertex_removed[v] && !edge_removed[p])
                    edge_to[j++] = q[edge_to[p]];
            }
            new_d[u] = j - new_pstart[u];
        }
        new_pstart[new_n] = j;
        delete[] g->d;
        delete[] g->pstart;
        g->d = new_d;
        g->pstart = new_pstart;
        deg = new_d;
        pstart = new_pstart;
        if (j * 2 < m)
        {
            ui *new_edge_to = new ui[j];
            memcpy(new_edge_to, edge_to, sizeof(ui) * j);
            delete[] g->edge_to;
            g->edge_to = new_edge_to;
            edge_to = new_edge_to;

            if (triangles != nullptr)
                delete[] triangles;
            triangles = new ui[j];
        }
        g->m = m = j;
        g->n = n = new_n;
        if (que == nullptr)
        {
            delete[] q;
        }
    }
    /**
     * @brief core-truss co-pruning
     */
    void CTCP(int lb)
    {
        edge_to = g->edge_to;
        pstart = g->pstart;
        deg = g->d;
        queue<pii> q_edge; //(edge_id, from) edge_to[edge_id]=to
        queue<ui> q_vertex;
        countTriangles();
        vector<bool> edge_removed(m, 0); // in_queue
        vector<bool> vertex_removed(n, 0);
        for (ui u = 0; u < n; u++)
        {
            if (deg[u] + paramK <= lb)
            {
                q_vertex.push(u);
                vertex_removed[u] = 1;
                continue;
            }
            for (ui i = pstart[u]; i < pstart[u + 1]; i++)
            {
                ui v = edge_to[i];
                if (triangles[i] + 2 * paramK <= lb)
                {
                    edge_removed[i] = 1;
                    if (u < v)
                        q_edge.push({i, u});
                }
            }
        }
        vector<bool> actually_rm(m, 0); // pop from queue
        while (q_edge.size() || q_vertex.size())
        {
            while (q_edge.size())
            {
                auto edge = q_edge.front();
                q_edge.pop();
                ui edge_id = edge.x, u = edge.y, v = edge_to[edge_id];
                actually_rm[edge_id] = 1;
                actually_rm[find(edge_to + pstart[v], edge_to + pstart[v + 1], u) + pstart[v]] = 1;
                if (--deg[u] + paramK <= lb && !vertex_removed[u])
                {
                    q_vertex.push(u);
                    vertex_removed[u] = 1;
                }
                if (--deg[v] + paramK <= lb && !vertex_removed[v])
                {
                    q_vertex.push(v);
                    vertex_removed[v] = 1;
                }
                ui *a = edge_to + pstart[u], *b = edge_to + pstart[u + 1];
                ui *l = edge_to + pstart[v], *r = edge_to + pstart[v + 1];
                // update the info of other edge triangles
                while (a < b && l < r)
                {
                    // note that we need to consider a situation where an edge is in the q_edge but not removed yet
                    if (actually_rm[a - edge_to] || vertex_removed[*a])
                    {
                        a++;
                        continue;
                    }
                    if (actually_rm[l - edge_to] || vertex_removed[*l])
                    {
                        l++;
                        continue;
                    }
                    if (*a == *l)
                    {
                        ui w = *a; // w is the common neighbor of u,v
                        ui id_uw = a - edge_to;
                        ui id_wu = find(edge_to + pstart[w], edge_to + pstart[w + 1], u) + pstart[w];
                        if (!edge_removed[id_uw])
                        {
                            --triangles[id_uw];
                            if (triangles[id_uw] + paramK * 2 <= lb)
                            {
                                edge_removed[id_uw] = edge_removed[id_wu] = 1;
                                q_edge.push({id_uw, u});
                            }
                            triangles[id_wu]--;
                        }
                        ui id_vw = l - edge_to;
                        ui id_wv = find(edge_to + pstart[w], edge_to + pstart[w + 1], v) + pstart[w];
                        if (!edge_removed[id_vw])
                        {
                            --triangles[id_vw];
                            if (triangles[id_vw] + paramK * 2 <= lb)
                            {
                                edge_removed[id_vw] = edge_removed[id_wv] = 1;
                                q_edge.push({id_vw, v});
                            }
                            triangles[id_wv]--;
                        }
                        a++;
                        l++;
                    }
                    else if (*a < *l)
                        a++;
                    else
                        l++;
                }
            }
            if (q_vertex.size())
            {
                ui u = q_vertex.front();
                q_vertex.pop();
                for (ui i = pstart[u]; i < pstart[u + 1]; i++)
                {
                    if (edge_removed[i])
                        continue;
                    ui v = edge_to[i];
                    if (vertex_removed[v])
                        continue;
                    if (--deg[v] + paramK <= lb)
                    {
                        vertex_removed[v] = 1;
                        q_vertex.push(v);
                    }
                }
                // update the triangles containing u
                for (ui i = pstart[u]; i < pstart[u + 1]; i++)
                {
                    if (edge_removed[i])
                        continue;
                    ui v = edge_to[i];
                    if (vertex_removed[v])
                        continue;
                    ui *st = edge_to + pstart[v], *ed = edge_to + pstart[v + 1];
                    for (ui j = i + 1; j < pstart[u + 1]; j++)
                    {
                        if (edge_removed[j])
                            continue;
                        ui w = edge_to[j]; // v w are neighbors of u
                        if (vertex_removed[w])
                            continue;
                        if (has(st, ed, w)) // v is connected to w
                        {
                            ui id_vw = find(st, ed, w) + pstart[v];
                            ui id_wv = find(edge_to + pstart[w], edge_to + pstart[w + 1], v) + pstart[w];
                            assert(edge_to[id_vw] == w);
                            assert(edge_to[id_wv] == v);
                            if (edge_removed[id_vw])
                                continue;
                            assert(triangles[id_vw] == triangles[id_wv]);
                            if (--triangles[id_vw] + 2 * paramK <= lb)
                            {
                                edge_removed[id_vw] = edge_removed[id_wv] = 1;
                                q_edge.push({id_vw, v});
                            }
                            else
                                --triangles[id_wv];
                        }
                    }
                }
            }
        }
        rebuild_graph(vertex_removed, edge_removed);
    }
};

#endif