#include "Graph.h"
double threshold = 1e9;
Timer thresh, branchings, bounding;
#include <numeric>
// #include "MSearcher.h"
#include "KPlex_BB_matrix.h"
#include "CTPrune.h"
#define K_THRESH 10
using namespace std;

Graph::Graph(const char *_dir, const int _K)
{
	dir = string(_dir);
	K = _K;

	n = m = 0;

	pstart = nullptr;
	pend = pend_buf = nullptr;
	edges = nullptr;

	kplex.clear();

	s_degree = s_edges = NULL;
	s_pstart = s_pend = NULL;
	s_peel_sequence = s_core = NULL;
	s_vis = NULL;
	s_heap = NULL;

	s_edgelist_pointer = NULL;
	s_tri_cnt = s_edge_list = NULL;
	s_active_edgelist = NULL;
	s_deleted = NULL;
}

Graph::~Graph()
{
	if (pstart != nullptr)
	{
		delete[] pstart;
		pstart = nullptr;
	}
	if (pend != nullptr)
	{
		delete[] pend;
		pend = nullptr;
	}
	if (pend_buf != NULL)
	{
		delete[] pend_buf;
		pend_buf = NULL;
	}
	if (edges != nullptr)
	{
		delete[] edges;
		edges = nullptr;
	}
	if (s_degree != NULL)
	{
		delete[] s_degree;
		s_degree = NULL;
	}
	if (s_pstart != NULL)
	{
		delete[] s_pstart;
		s_pstart = NULL;
	}
	if (s_pend != NULL)
	{
		delete[] s_pend;
		s_pend = NULL;
	}
	if (s_edges != NULL)
	{
		delete[] s_edges;
		s_edges = NULL;
	}
	if (s_peel_sequence != NULL)
	{
		delete[] s_peel_sequence;
		s_peel_sequence = NULL;
	}
	if (s_core != NULL)
	{
		delete[] s_core;
		s_core = NULL;
	}
	if (s_vis != NULL)
	{
		delete[] s_vis;
		s_vis = NULL;
	}
	if (s_heap != NULL)
	{
		delete s_heap;
		s_heap = NULL;
	}
	if (s_edgelist_pointer != NULL)
	{
		delete[] s_edgelist_pointer;
		s_edgelist_pointer = NULL;
	}
	if (s_active_edgelist != NULL)
	{
		delete[] s_active_edgelist;
		s_active_edgelist = NULL;
	}
	if (s_deleted != NULL)
	{
		delete[] s_deleted;
		s_deleted = NULL;
	}
}

void Graph::read()
{
	FILE *f = Utility::open_file(dir.c_str(), "rb");

	ui tt;
	fread(&tt, sizeof(int), 1, f);
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	// printf("\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, f);
	if (pstart == nullptr)
		pstart = new ept[n + 1];
	if (edges == nullptr)
		edges = new ui[m];
	fread(edges, sizeof(int), m, f);

	pstart[0] = 0;
	std::partial_sum(degree, degree + n, pstart + 1);
	// for(ui i = 0;i < n;i ++) {
	// if(degree[i] > 0) {
	// 	fread(edges+pstart[i], sizeof(int), degree[i], f);

	// 	// remove self loops and parallel edges
	// 	ui *buff = edges+pstart[i];
	// 	sort(buff, buff+degree[i]);
	// 	ui idx = 0;
	// 	for(ui j = 0;j < degree[i];j ++) {
	// 		if(buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);
	// 		if(buff[j] == i||(j > 0&&buff[j] == buff[j-1])) continue;
	// 		buff[idx ++] = buff[j];
	// 	}
	// 	degree[i] = idx;
	// }

	// pstart[i+1] = pstart[i] + degree[i];
	// }

	fclose(f);
	delete[] degree;
}

void Graph::write()
{
	FILE *fout = Utility::open_file("kplexes.txt", "w");
	fprintf(fout, "%lu\n", kplex.size());
	sort(kplex.begin(), kplex.end());
	for (ui i = 0; i < kplex.size(); i++)
		fprintf(fout, " %u", kplex[i]);
	fprintf(fout, "\n");
	fclose(fout);
}

void Graph::verify_kplex()
{
	char *vis = new char[n];
	memset(vis, 0, sizeof(char) * n);

	FILE *fin = Utility::open_file("kplexes.txt", "r");

	ui kplex_size = n, kplex_n, idx = 0;
	char ok = 1;
	while (fscanf(fin, "%u", &kplex_n) == 1)
	{
		++idx;
		if (kplex_size == n)
		{
			kplex_size = kplex_n;
			printf("k-plex sizes: %u\n", kplex_size);
		}
		if (kplex_n != kplex_size)
			printf("!!! WA k-plex size: %u!\n", kplex_n);
		vector<ui> kplex;
		for (ui i = 0; i < kplex_n; i++)
		{
			ui tmp;
			fscanf(fin, "%u", &tmp);
			kplex.pb(tmp);
		}

		for (ui i = 0; i < kplex.size(); i++)
		{
			if (vis[kplex[i]])
			{
				printf("WA k-plex! Duplicate vertex: %u\n", idx);
				ok = 0;
				break;
			}
			vis[kplex[i]] = 1;
		}
		for (ui i = 0; i < kplex.size(); i++)
		{
			ui d = 0;
			for (ui j = pstart[kplex[i]]; j < pstart[kplex[i] + 1]; j++)
				if (vis[edges[j]])
					++d;
			if (d + K < kplex.size())
			{
				ok = 0;
				printf("WA k-plex! Not enough neighbors!\n");
			}
		}
		for (ui i = 0; i < kplex.size(); i++)
			vis[kplex[i]] = 0;
	}
	if (ok)
		printf("Correct k-plexes!\n");
	fclose(fin);

	delete[] vis;
}

void Graph::search()
{
	cout << "graph_size: " << n << " " << m << endl;
	Timer t;
	kplex.resize(1); // screen out trivial cases
	ui max_degree = 0;
	for (ui i = 0; i < n; i++)
	{
		if (pstart[i + 1] - pstart[i] > max_degree)
			max_degree = pstart[i + 1] - pstart[i];
	}

	cout << "initial_max_degree: " << max_degree << endl;
	// heuristic_kplex_max_degree(10);
	ui oldn = n;
	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *degree = new ui[n];
	char *vis = new char[n];

	ListLinearHeap *heap = new ListLinearHeap(n, n - 1);

	const ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, false);
	if (kplex.size() < UB)
	{
		ui old_size = kplex.size();
		ui *out_mapping = new ui[n];
		ui *rid = new ui[n];
		ui *edgelist_pointer = new ui[m];

		shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);
		ego_degen(n, m, peel_sequence, pstart, edges, degree, rid, vis, heap, edgelist_pointer, true);

		if (kplex.size() > old_size)
		{
			old_size = kplex.size();
			for (ui i = 0; i < kplex.size(); i++)
			{
				assert(kplex[i] < n);
				kplex[i] = out_mapping[kplex[i]];
			}
		}

		if (kplex.size() + 1 > 2 * K)
			CTPrune::core_truss_copruning(n, m, kplex.size() + 1 - K, kplex.size() + 1 - 2 * K, peel_sequence, out_mapping, rid, pstart, edges, degree, true);
		// else
		// 	shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);
		ui max_degree = 0;
		cout << "after shrinking: " << n << endl;
		degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, true);

		for (ui i = 0; i < n; i++)
			max_degree = max(max_degree, pstart[i + 1] - pstart[i]);
		cout << "residual_max_degree: " << max_degree << endl;
	}
}

void Graph::search_dense()
{

	Timer t;
	dense_search = true;
	read(); // read the graph again...
	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *degree = new ui[n];
	char *vis = new char[n];
	vector<ui> dense_kplex = kplex;
	kplex.pop_back(); // removing one item from kplex, so that degen can get |P|-k core
	ui max_n_edges = dense_kplex.size() * (dense_kplex.size() - 1);

	ListLinearHeap *heap = new ListLinearHeap(n, n - 1);
	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, false);

	// delete heap;
	// delete[] vis;
	// delete[] degree;
	ui init_edges = 0;

	if (n > dense_kplex.size() && best_n_edges < max_n_edges)
	{

		ui *out_mapping = new ui[n];
		ui *rid = new ui[n];
		ui *edgelist_pointer = new ui[m];

		shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);

		delete[] core;
		core = NULL;

		ui *degree = new ui[n];
		for (ui i = 0; i < n; i++)
			degree[i] = pstart[i + 1] - pstart[i];

		ListLinearHeap *linear_heap = new ListLinearHeap(n, n - 1);
		linear_heap->init(n, n - 1, peel_sequence, degree);

		// assert(pend == nullptr);
		pend = new ept[n];

		oriented_triangle_counting(n, m, peel_sequence, pstart, pend, edges, edgelist_pointer, rid); // edgelist_pointer currently stores triangle_counts

		// delete[] peel_sequence; peel_sequence = NULL;

		pend_buf = new ept[n];
		ui *edge_list = new ui[m];
		ui *tri_cnt = new ui[m / 2];
		reorganize_oriented_graph(n, tri_cnt, edge_list, pstart, pend, pend_buf, edges, edgelist_pointer, rid);

		for (ui i = 0; i < n; i++)
			pend[i] = pstart[i + 1];

		ui *active_edgelist = new ui[m >> 1];
		ui active_edgelist_n = m >> 1;
		for (ui i = 0; i < (m >> 1); i++)
			active_edgelist[i] = i;

		ui *Qe = new ui[m >> 1];
		char *deleted = new char[m >> 1];
		memset(deleted, 0, sizeof(char) * (m >> 1));
		char *exists = new char[n];
		memset(exists, 0, sizeof(char) * n);

		ui *Qv = new ui[n];
		ui Qv_n = 0;
		KPLEX_BB_MATRIX *kplex_solver = new KPLEX_BB_MATRIX(true); // true means solving for dense search
		kplex_solver->allocateMemory(n);
		kplex_solver->best_n_edges = best_n_edges;

		if (kplex.size() > 2 * K - 2)
		{
			m -= 2 * peeling(n, linear_heap, Qv, Qv_n, kplex.size() + 1 - K, Qe, true, kplex.size() + 1 - 2 * K, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
			printf("*** After core-truss co-pruning: n = %s, m = %s, density = %.4lf\n", Utility::integer_to_string(n - Qv_n).c_str(), Utility::integer_to_string(m / 2).c_str(), double(m) / (n - Qv_n) / (n - Qv_n - 1));
		}

		Timer tt;

		int max_n = n - Qv_n;
		s_degree = new ui[max_n];
		s_pstart = new ept[max_n + 1];
		s_pend = new ept[max_n];
		s_edges = new ui[m];
		s_peel_sequence = new ui[max_n];
		s_core = new ui[max_n];
		s_vis = new char[max_n];
		s_heap = new ListLinearHeap(max_n, max_n - 1);
		s_edgelist_pointer = new ui[m];
		s_tri_cnt = new ui[m / 2];
		s_edge_list = new ui[m];
		s_active_edgelist = new ui[m / 2];
		s_deleted = new char[m / 2];

		vector<pair<ui, ui>> vp;
		vp.reserve(m / 2);
		ui *t_degree = new ui[n];

		ui max_n_prune = 0, max_n_search = 0, prune_cnt = 0, search_cnt = 0;
		double min_density_prune = 1, min_density_search = 1, total_density_prune = 0, total_density_search = 0;
		ui last_m = 0;

		for (ui i = 0; i < n && m && best_n_edges < max_n_edges; i++)
		{
			ui u, key;
			linear_heap->pop_min(u, key);
			// if(key != 0) printf("u = %u, key = %u\n", u, key);
			if (key < kplex.size() + 1 - K)
			{
				if (degree[u] != 0)
				{ // degree[u] == 0 means u is deleted. it could be the case that degree[u] == 0, but key[u] > 0, as key[u] is not fully updated in linear_heap
					Qv[0] = u;
					Qv_n = 1;
					if (kplex.size() + 1 > 2 * K)
						m -= 2 * peeling(n, linear_heap, Qv, Qv_n, kplex.size() + 1 - K, Qe, false, kplex.size() + 1 - 2 * K, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
					else
						m -= 2 * peeling(n, linear_heap, Qv, Qv_n, kplex.size() + 1 - K, Qe, false, 0, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
				}
				continue;
			}
			if (m == 0)
				break;
			assert(degree[u] == key);
			if (thresh.elapsed() / 1e6 >= threshold)
				break;

			ui *ids = Qv;
			ui ids_n = 0;
			bool mflag = false;
			// cout<<u<<endl;

			bool check = false;

			extract_subgraph_and_prune(u, ids, ids_n, rid, vp, Qe, t_degree, exists, pend, deleted, edgelist_pointer);
			if (ids_n)
			{
				double density = (double(vp.size() * 2)) / ids_n / (ids_n - 1);
				total_density_prune += density;
				++prune_cnt;
				if (density < min_density_prune)
					min_density_prune = density;
				if (ids_n > max_n_prune)
					max_n_prune = ids_n;
				// cout<<"Density"<<density<<" ";
			}
			if (K < K_THRESH && ids_n > kplex.size() && vp.size() * 2 < m)
				subgraph_prune(ids, ids_n, vp, rid, Qv, Qe, exists);

			if (ids_n > kplex.size())
			{
				double density = (double(vp.size() * 2)) / ids_n / (ids_n - 1);
				total_density_search += density;
				++search_cnt;
				if (density < min_density_search)
					min_density_search = density;
				if (ids_n > max_n_search)
					max_n_search = ids_n;
				ui presize = kplex.size();
				// cout<<"searching: "<<u<<" -> ids_n "<<ids_n<<" density: "<<density<<endl;
				kplex_solver->load_graph(ids_n, vp);
				kplex_solver->kPlex(K, UB, kplex, true);
				if (init_edges == 0 && best_n_edges > 0)
					init_edges = best_n_edges;

				if (kplex.size() > presize)
				{
					if (kplex_solver->best_n_edges > best_n_edges)
					{
						best_n_edges = kplex_solver->best_n_edges;
						cout << "A denser kplex found with #edges: " << best_n_edges << endl;
						dense_kplex.clear();
						for (ui j = 0; j < kplex.size(); j++)
							dense_kplex.push_back(ids[kplex[j]]);
					}
					kplex.pop_back();
				}
			}
			Qv[0] = u;
			Qv_n = 1;
			m -= 2 * peeling(n, linear_heap, Qv, Qv_n, kplex.size() + 1 - K, Qe, false, kplex.size() + 1 - 2 * K, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
		}

		if (prune_cnt == 0)
			++prune_cnt;
		if (search_cnt == 0)
			++search_cnt;
		printf("prune_cnt: %u, max_n: %u, min_density: %.4lf, avg_density: %.4lf\n", prune_cnt, max_n_prune, min_density_prune, total_density_prune / prune_cnt);
		printf("search_cnt: %u, max_n: %u, min_density: %.4lf, avg_density: %.4lf\n", search_cnt, max_n_search, min_density_search, total_density_search / search_cnt);
		// printf("*** Search time: %s \n", Utility::integer_to_string(tt.elapsed()).c_str());
		// printf(">>%s t_Search: %f", dir.substr(dir.find_last_of("/")).c_str(), tt.elapsed()/1000000.0);

		printf(">>%s-dense \tMaxKPlex_Size: %lu t_Total: %f n_mkp: %d initial_edges: %d densest_kplex_edges: %d\n", dir.substr(dir.find_last_of("/") + 1).c_str(), dense_kplex.size(), t.elapsed() / 1e6, kplex_solver->nmkp, init_edges, best_n_edges);

		delete kplex_solver;
		delete linear_heap;
		delete[] t_degree;
		delete[] exists;
		delete[] out_mapping;
		delete[] rid;
		delete[] degree;
		delete[] edgelist_pointer;
		delete[] tri_cnt;
		delete[] active_edgelist;
		delete[] Qe;
		delete[] Qv;
		delete[] deleted;
	}
	delete[] core;
	delete[] peel_sequence;

	// printf("\tMaxKPlex_Size: %lu t_Total: %f t_Seesaw: %f\n", kplex.size(), t.elapsed()/1000000.0, 0);
	// printf("\tMaximum kPlex Size: %lu, Total Time: %s (microseconds)\n", kplex.size(), Utility::integer_to_string(t.elapsed()).c_str());
	// printf("*** Node count: %lld\n",nodeCnt);
	// printf("*** Gamma avg: %.4lf\n",(edgeCnt+0.1)/(nodeCnt+0.1));
	// printf("*** Sub avg: %.4lf\n",(subSum+0.1)/(subCnt+0.1));
	// printf("*** Sub max: %lld\n",subMax);
	// printf("*** Max degree: %d\n",max_degree);
	if (kplex.size() == 2 * K - 2)
		printf("!!! Warning: Trivial Case !!!\n");
}

void Graph::write_subgraph(ui n, const vector<pair<ui, ui>> &edge_list)
{
	FILE *fout = Utility::open_file("edges.txt", "w");

	fprintf(fout, "%u %lu\n", n, edge_list.size());
	for (ui i = 0; i < edge_list.size(); i++)
		fprintf(fout, "%d %d\n", edge_list[i].first, edge_list[i].second);

	fclose(fout);
}

void Graph::subgraph_prune(ui *ids, ui &_n, vector<pair<ui, ui>> &edge_list, ui *rid, ui *Qv, ui *Qe, char *exists)
{
	ui s_n;
	ept s_m;
	load_graph_from_edgelist(_n, edge_list, s_n, s_m, s_degree, s_pstart, s_edges);
	degen(s_n, s_peel_sequence, s_core, s_pstart, s_edges, s_degree, s_vis, s_heap, false);
	shrink_graph(s_n, s_m, s_peel_sequence, s_core, ids, ids, rid, s_pstart, s_edges, false);

	if (s_n > 0 && kplex.size() + 1 > 2 * K)
	{
		// if(false) {
		//  printf("before n = %u, m = %lu\n", s_n, s_m);
		oriented_triangle_counting(s_n, s_m, s_peel_sequence, s_pstart, s_pend, s_edges, s_edgelist_pointer, rid);
		reorganize_oriented_graph(s_n, s_tri_cnt, s_edge_list, s_pstart, s_pend, pend_buf, s_edges, s_edgelist_pointer, rid);
		for (ui i = 0; i < s_n; i++)
		{
			s_pend[i] = s_pstart[i + 1];
			s_degree[i] = s_pend[i] - s_pstart[i];
		}
		for (ept i = 0; i < (s_m >> 1); i++)
			s_active_edgelist[i] = i;
		memset(s_deleted, 0, sizeof(char) * (s_m >> 1));

		ui s_active_edgelist_n = s_m >> 1;
		ui Qv_n = 0;
		s_m -= 2 * peeling(0, NULL, Qv, Qv_n, kplex.size() + 1 - K, Qe, true, kplex.size() + 1 - 2 * K, s_tri_cnt, s_active_edgelist, s_active_edgelist_n, s_edge_list, s_edgelist_pointer, s_deleted, s_degree, s_pstart, s_pend, s_edges, exists);
		for (ui i = 0; i < Qv_n; i++)
			if (Qv[i] == 0)
			{
				_n = 0;
				return;
			}
		extract_subgraph(0, Qv, s_n, rid, edge_list, exists, s_pstart, s_pend, s_edges, s_deleted, s_edgelist_pointer);
		for (ui i = 0; i < s_n; i++)
			Qv[i] = ids[Qv[i]];
		for (ui i = 0; i < s_n; i++)
			ids[i] = Qv[i];
		_n = s_n;
		assert(edge_list.size() * 2 == s_m);
		// printf("*after n = %u, m = %lu\n", s_n, s_m);
	}
	else
	{
		_n = s_n;
		edge_list.clear();
		for (ui i = 0; i < s_n; i++)
			for (ept j = s_pstart[i]; j < s_pstart[i + 1]; j++)
				if (s_edges[j] > i)
				{
					edge_list.push_back(make_pair(s_edges[j], i));
				}
		assert(edge_list.size() * 2 == s_m);
	}
}

void Graph::load_graph_from_edgelist(ui _n, const vector<pair<ui, ui>> &edge_list, ui &n, ept &m, ui *degree, ept *pstart, ui *edges)
{
	n = _n;
	m = (ui)edge_list.size() * 2;
	for (ui i = 0; i < n; i++)
		degree[i] = 0;
	for (ept i = 0; i < m / 2; i++)
	{
		assert(edge_list[i].first >= 0 && edge_list[i].first < n && edge_list[i].second >= 0 && edge_list[i].second < n);
		degree[edge_list[i].first]++;
		degree[edge_list[i].second]++;
	}

	pstart[0] = 0;
	for (ui i = 0; i < n; i++)
		pstart[i + 1] = pstart[i] + degree[i];
	for (ept i = 0; i < m / 2; i++)
	{
		ui a = edge_list[i].first, b = edge_list[i].second;
		edges[pstart[a]++] = b;
		edges[pstart[b]++] = a;
	}
	for (ui i = 0; i < n; i++)
		pstart[i] -= degree[i];
}

void Graph::extract_graph(ui n, ui m, ui *degree, ui *ids, ui &ids_n, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer)
{
	ids_n = 0;
	vp.clear();
	for (ui i = 0; i < n; ++i)
	{
		if (degree[i])
		{
			ids[ids_n] = i;
			rid[i] = ids_n++;
		}
	}
	for (ui i = 0; i < ids_n; i++)
	{
		ui u = ids[i];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]] && u < edges[j])
			{
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
	}
}

void Graph::extract_subgraph(ui u, ui *ids, ui &ids_n, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer)
{
	ids_n = 0;
	vp.clear();
	ids[ids_n++] = u;
	exists[u] = 1;
	rid[u] = 0;
	ui u_n = pstart[u];
	for (ept i = pstart[u]; i < pend[u]; i++)
		if (!deleted[edgelist_pointer[i]])
		{
			edges[u_n] = edges[i];
			edgelist_pointer[u_n++] = edgelist_pointer[i];
			ui v = edges[i];
			rid[v] = ids_n;
			ids[ids_n++] = v;
			exists[v] = 1;
		}
	pend[u] = u_n;
	ui old_size = ids_n;
	for (ui i = 1; i < old_size; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				ui v = edges[j];
				if (exists[v])
					continue;
				rid[v] = ids_n;
				ids[ids_n++] = v;
				exists[v] = 1;
			}
		pend[u] = u_n;
	}
	for (ui i = 0; i < old_size; i++)
	{
		u = ids[i];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (edges[j] > u)
			{
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
	}
	for (ui i = old_size; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (edges[j] > u && exists[edges[j]])
					vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
		pend[u] = u_n;
	}
	for (ui i = 0; i < ids_n; i++)
		exists[ids[i]] = 0;
}

void Graph::extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, ui *rid, vector<pair<ui, ui>> &vp, ui *Q, ui *degree, char *exists, ept *pend, char *deleted, ui *edgelist_pointer)
{
	vp.clear();
	ids_n = 0;
	ids[ids_n++] = u;
	exists[u] = 1;
	ui u_n = pstart[u];
	for (ept i = pstart[u]; i < pend[u]; i++)
		if (!deleted[edgelist_pointer[i]])
		{
			edges[u_n] = edges[i];
			edgelist_pointer[u_n++] = edgelist_pointer[i];
			ui v = edges[i];
			ids[ids_n++] = v;
			exists[v] = 2;
		}
	pend[u] = u_n;

	ui Q_n = 0;
	for (ui i = 1; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		degree[u] = 0;
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (exists[edges[j]] == 2)
					++degree[u];
			}
		pend[u] = u_n;
		if (degree[u] + 2 * K <= kplex.size())
			Q[Q_n++] = u;
	}
	for (ui i = 0; i < Q_n; i++)
	{
		u = Q[i];
		exists[u] = 10;
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (exists[edges[j]] == 2)
			{
				if ((degree[edges[j]]--) + 2 * K == kplex.size() + 1)
				{
					assert(Q_n < m / 2);
					Q[Q_n++] = edges[j];
				}
			}
	}
	assert(Q_n <= ids_n);
	if (ids_n - 1 - Q_n + K <= kplex.size())
	{
		for (ui i = 0; i < ids_n; i++)
			exists[ids[i]] = 0;
		ids_n = 0;
		return;
	}

	ui nr_size = ids_n;
	for (ui i = 1; i < nr_size; i++)
		if (exists[ids[i]] == 2)
		{
			u = ids[i];
			for (ept j = pstart[u]; j < pend[u]; j++)
			{
				if (!exists[edges[j]])
				{
					ids[ids_n++] = edges[j];
					exists[edges[j]] = 3;
					degree[edges[j]] = 1;
				}
				else if (exists[edges[j]] == 3)
					++degree[edges[j]];
			}
		}

#ifndef NDEBUG
	// printf("Entire list: ");
	// for(ui i = 0;i < nr_size;i ++) printf(" %u", ids[i]);
	// printf("\n");
#endif

	ui new_size = 1;
	for (ui i = 1; i < nr_size; i++)
	{
		if (exists[ids[i]] == 10)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
#ifndef NDEBUG
	if (new_size + Q_n != nr_size)
	{
		printf("new_size: %u, Q_n: %u, nr_size: %u\n", new_size, Q_n, nr_size);
		printf("New list: ");
		for (ui i = 0; i < new_size; i++)
			printf(" %u", ids[i]);
		printf("\n");
		printf("Pruned list: ");
		for (ui i = 0; i < Q_n; i++)
			printf(" %u", Q[i]);
		printf("\n");
	}
#endif
	assert(new_size + Q_n == nr_size);
	ui old_nr_size = nr_size;
	nr_size = new_size;
	for (ui i = old_nr_size; i < ids_n; i++)
	{
		if (degree[ids[i]] + 2 * K <= kplex.size() + 2)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	ids_n = new_size;
#ifndef NDEBUG
	assert(exists[ids[0]] == 1);
	for (ui i = 1; i < nr_size; i++)
		assert(exists[ids[i]] == 2);
	for (ui i = nr_size; i < ids_n; i++)
		assert(exists[ids[i]] == 3);
#endif

	// for(ui i = 0;i < ids_n;i ++) printf(" %u", ids[i]);
	// printf("\n");

	for (ui i = 0; i < ids_n; i++)
	{
		assert(exists[ids[i]]);
		rid[ids[i]] = i;
	}

	for (ui i = 0; i < nr_size; i++)
	{
		u = ids[i];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (exists[edges[j]] && edges[j] > u)
			{
				assert(!deleted[edgelist_pointer[j]]);
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
	}
	for (ui i = nr_size; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (edges[j] > u && exists[edges[j]])
					vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
		pend[u] = u_n;
	}
	for (ui i = 0; i < ids_n; i++)
		exists[ids[i]] = 0;
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(exists[i] == 0);
#endif
}

// max-degree-based heuristic k-plex computation
void Graph::heuristic_kplex_max_degree(ui processed_threshold)
{
	Timer t;
	ui *head = new ui[n];
	ui *next = new ui[n];
	ui *degree = new ui[n];

	ui *vis = new ui[n];
	memset(vis, 0, sizeof(ui) * n);

	ui max_degree = 0;
	for (ui i = 0; i < n; i++)
		head[i] = n;
	for (ui i = 0; i < n; i++)
	{
		degree[i] = pstart[i + 1] - pstart[i];
		if (degree[i] > max_degree)
			max_degree = degree[i];
		next[i] = head[degree[i]];
		head[degree[i]] = i;
	}

	for (ui processed_vertices = 0; max_degree + K >= kplex.size() && processed_vertices < processed_threshold; processed_vertices++)
	{
		ui u = n;
		while (max_degree >= 0 && max_degree + K >= kplex.size() && u == n)
		{
			for (ui v = head[max_degree]; v != n;)
			{
				ui tmp = next[v];
				if (degree[v] == max_degree)
				{
					u = v;
					head[max_degree] = tmp;
					break;
				}
				else if (degree[v] + K >= kplex.size())
				{
					next[v] = head[degree[v]];
					head[degree[v]] = v;
				}
				v = tmp;
			}
			if (u == n)
			{
				head[max_degree] = n;
				--max_degree;
			}
		}
		if (u == n)
			break;

		vis[u] = 1;
		for (ui k = pstart[u]; k < pstart[u + 1]; k++)
			if (!vis[edges[k]])
				--degree[edges[k]];

		vector<ui> vs;
		for (ui j = pstart[u]; j < pstart[u + 1]; j++)
			if (!vis[edges[j]])
				vs.pb(edges[j]);

		vector<ui> vs_deg(vs.size());
		for (ui j = 0; j < vs.size(); j++)
			vis[vs[j]] = 2;
		for (ui j = 0; j < vs.size(); j++)
		{
			ui v = vs[j], d = 0;
			for (ui k = pstart[v]; k < pstart[v + 1]; k++)
			{
				if (vis[edges[k]] == 2)
					++d;
			}
			vs_deg[j] = d;
		}
		for (ui j = 0; j < vs.size(); j++)
			vis[vs[j]] = 0;

		vector<ui> res;
		res.pb(u);
		ui vs_size = vs.size();
		while (vs_size > 0 && res.size() + vs_size + K - 1 > kplex.size())
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

			ui new_size = 0;
			for (ui j = pstart[u]; j < pstart[u + 1]; j++)
				if (!vis[edges[j]])
					vis[edges[j]] = 2;
			for (ui j = 0; j < vs_size; j++)
				if (vis[vs[j]])
				{
					if (j != new_size)
						swap(vs[new_size], vs[j]);
					vs_deg[new_size] = vs_deg[j];
					++new_size;
				}
			for (ui j = pstart[u]; j < pstart[u + 1]; j++)
				if (vis[edges[j]] == 2)
					vis[edges[j]] = 0;

			res.pb(u);
			for (ui k = 0; k < new_size; k++)
				vis[vs[k]] = k + 2;
			for (ui j = new_size; j < vs_size; j++)
			{
				ui v = vs[j];
				for (ui k = pstart[v]; k < pstart[v + 1]; k++)
				{
					if (vis[edges[k]] >= 2)
						--vs_deg[vis[edges[k]] - 2];
				}
			}
			for (ui k = 0; k < new_size; k++)
				vis[vs[k]] = 0;

			vs_size = new_size;
		}

		// TO DO: extend res to be a maximal k-plex

		if (res.size() > kplex.size())
			kplex = res;
	}

	delete[] vis;
	delete[] head;
	delete[] next;
	delete[] degree;

	printf("*** Heuristic kplex size: %lu, time: %s (microseconds)\n", kplex.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

// degeneracy-based k-plex
// return an upper bound of the maximum k-plex size
ui Graph::degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output)
{
	Timer t;

	ui threshold = (kplex.size() + 1 > K ? kplex.size() + 1 - K : 0);

	for (ui i = 0; i < n; i++)
		degree[i] = pstart[i + 1] - pstart[i];

	ui queue_n = 0, new_size = 0;
	for (ui i = 0; i < n; i++)
		if (degree[i] < threshold)
			peel_sequence[queue_n++] = i;
	for (ui i = 0; i < queue_n; i++)
	{
		ui u = peel_sequence[i];
		degree[u] = 0;
		for (ept j = pstart[u]; j < pstart[u + 1]; j++)
			if (degree[edges[j]] > 0)
			{
				if ((degree[edges[j]]--) == threshold)
					peel_sequence[queue_n++] = edges[j];
			}
	}
	ui UB = n;
	if (queue_n == n)
		UB = kplex.size();

	memset(vis, 0, sizeof(char) * n);
	for (ui i = 0; i < n; i++)
	{
		if (degree[i] >= threshold)
			peel_sequence[queue_n + (new_size++)] = i;
		else
		{
			vis[i] = 1;
			core[i] = 0;
		}
	}
	vector<ui> kcore_sizes;
	assert(queue_n + new_size == n);
	if (new_size != 0)
	{
		heap->init(new_size, new_size - 1, peel_sequence + queue_n, degree);
		ui max_core = 0;
		ui idx = n;
		UB = 0;
		for (ui i = 0; i < new_size; i++)
		{
			ui u, key;
			heap->pop_min(u, key);

			if (key > max_core)
			{
				max_core = key;
				kcore_sizes.push_back(heap->sz);
			}
			core[u] = max_core;
			peel_sequence[queue_n + i] = u;

			ui t_UB = core[u] + K;
			if (new_size - i < t_UB)
				t_UB = new_size - i;
			if (t_UB > UB)
				UB = t_UB;

			if (idx == n && key + K >= new_size - i)
				idx = i;
			vis[u] = 1;

			for (ept j = pstart[u]; j < pstart[u + 1]; j++)
				if (!vis[edges[j]])
					heap->decrement(edges[j], 1);
		}

		if (!dense_search && new_size - idx > kplex.size())
		{
			kplex.clear();
			for (ui i = idx; i < new_size; i++)
				kplex.pb(peel_sequence[queue_n + i]);
			if (output)
				printf("Find a k-plex of size: %u\n", new_size - idx);
		}
		printf("*** max_core: %u, UB: %u\n", max_core, UB);
		if (output)
		{
			printf("*** max_core_residual: %u\n", max_core);
			cout << "kcore_size: ";

			for (ui i : kcore_sizes)
				cout << i << " ";
			cout << endl;
		}
	}
	else
	{
		cout << "kcore_size: " << 0 << endl;
	}
	return UB;
}

// in_mapping and out_mapping can be the same array
// note that core is not maintained, and is assumed to not be used anymore
void Graph::shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *pstart, ui *edges, bool output)
{
	ui cnt = 0;
	for (ui i = 0; i < n; i++)
		if (core[i] + K > kplex.size())
		{
			rid[i] = cnt;
			if (in_mapping == NULL)
				out_mapping[cnt] = i;
			else
				out_mapping[cnt] = in_mapping[i];
			++cnt;
		}

	if (cnt != n)
	{
		cnt = 0;
		ept pos = 0;
		for (ui i = 0; i < n; i++)
			if (core[i] + K > kplex.size())
			{
				ept t_start = pstart[i];
				pstart[cnt] = pos;
				for (ept j = t_start; j < pstart[i + 1]; j++)
					if (core[edges[j]] + K > kplex.size())
					{
						edges[pos++] = rid[edges[j]];
					}
				++cnt;
			}
		pstart[cnt] = pos;

		// printf("%u %u %u %u\n", n, cnt, core[peel_sequence[n-cnt-1]], core[peel_sequence[n-cnt]]);
		assert(core[peel_sequence[n - cnt - 1]] == 0 || core[peel_sequence[n - cnt - 1]] + K <= kplex.size());
		assert(cnt == 0 || core[peel_sequence[n - cnt]] + K > kplex.size());
		for (ui i = 0; i < cnt; i++)
		{
			peel_sequence[i] = rid[peel_sequence[n - cnt + i]];
			// core[i] = core[out_mapping[i]];
		}

		n = cnt;
		m = pos;
	}

	if (output)
		printf("*** After core shrink: n = %s, m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m / 2).c_str());
}

// orient graph and triangle counting
void Graph::oriented_triangle_counting(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj)
{
	ui *rid = adj;
	for (ui i = 0; i < n; i++)
		rid[peel_sequence[i]] = i;
	for (ui i = 0; i < n; i++)
	{
		ept &end = pend[i] = pstart[i];
		for (ept j = pstart[i]; j < pstart[i + 1]; j++)
			if (rid[edges[j]] > rid[i])
				edges[end++] = edges[j];
	}

#ifndef NDEBUG
	long long sum = 0;
	for (ui i = 0; i < n; i++)
		sum += pend[i] - pstart[i];
	// printf("%lld %lld\n", sum, m);
	assert(sum * 2 == m);
#endif

	memset(adj, 0, sizeof(ui) * n);
	long long cnt = 0;
	memset(tri_cnt, 0, sizeof(ui) * m);
	for (ui u = 0; u < n; u++)
	{
		for (ept j = pstart[u]; j < pend[u]; j++)
			adj[edges[j]] = j + 1;

		for (ept j = pstart[u]; j < pend[u]; j++)
		{
			ui v = edges[j];
			for (ept k = pstart[v]; k < pend[v]; k++)
				if (adj[edges[k]])
				{
					++tri_cnt[j];
					++tri_cnt[k];
					++tri_cnt[adj[edges[k]] - 1];
					++cnt;
				}
		}

		for (ept j = pstart[u]; j < pend[u]; j++)
			adj[edges[j]] = 0;
	}
#ifndef NDEBUG
	// printf("*** Total number of triangles: %s\n", Utility::integer_to_string(cnt).c_str());
#endif
}

// reorganize the adjacency lists
// and sort each adjacency list to be in increasing order
void Graph::reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf)
{
	for (ui i = 0; i < n; i++)
		pend2[i] = pend[i];
	ept pos = 0;
	for (ui i = 0; i < n; i++)
	{
		for (ept j = pstart[i]; j < pend[i]; j++)
		{
			tri_cnt[pos >> 1] = edgelist_pointer[j];
			edge_list[pos++] = i;
			edge_list[pos++] = edges[j];

			ept &k = pend2[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j] = (pos >> 1) - 1;
			edges[k++] = i;
		}
	}

#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend2[i] == pstart[i + 1]);
#endif

	for (ui i = 0; i < n; i++)
	{
		pend2[i] = pend[i];
		pend[i] = pstart[i];
	}
	for (ui i = 0; i < n; i++)
	{
		for (ept j = pend2[i]; j < pstart[i + 1]; j++)
		{
			ept &k = pend[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j];
			edges[k++] = i;
		}
	}

	ept *ids = pend2;
	for (ui i = 0; i < n; i++)
	{
		if (pend[i] == pstart[i] || pend[i] == pstart[i + 1])
			continue;
		ept j = pstart[i], k = pend[i], pos = 0;
		while (j < pend[i] && k < pstart[i + 1])
		{
			if (edges[j] < edges[k])
			{
				ids[pos] = edges[j];
				buf[pos++] = edgelist_pointer[j++];
			}
			else
			{
				ids[pos] = edges[k];
				buf[pos++] = edgelist_pointer[k++];
			}
		}
		while (j < pend[i])
		{
			ids[pos] = edges[j];
			buf[pos++] = edgelist_pointer[j++];
		}
		while (k < pstart[i + 1])
		{
			ids[pos] = edges[k];
			buf[pos++] = edgelist_pointer[k++];
		}
		for (ept j = 0; j < pos; j++)
		{
			edges[pstart[i] + j] = ids[j];
			edgelist_pointer[pstart[i] + j] = buf[j];
		}
	}
}

char Graph::find(ui u, ui w, ept &b, ept e, char *deleted, ept &idx, ui *edgelist_pointer, ui *edges)
{
	if (b >= e)
		return 0;

	while (b + 1 < e)
	{
		idx = b + (e - b) / 2;
		if (edges[idx] > w)
			e = idx;
		else
			b = idx;
	}

	if (edges[b] == w)
	{
		idx = edgelist_pointer[b];
		if (!deleted[idx])
			return 1;
	}

	return 0;
}

// return the number of peeled edges
ept Graph::peeling(ui critical_vertex, ListLinearHeap *linear_heap, ui *Qv, ui &Qv_n, ui d_threshold, ui *Qe, bool initialize_Qe, ui t_threshold, ui *tri_cnt, ui *active_edgelist, ui &active_edgelist_n, ui *edge_list, ui *edgelist_pointer, char *deleted, ui *degree, ept *pstart, ept *pend, ui *edges, char *exists)
{
	ept Qe_n = 0;
#ifndef NO_TRUSS_PRUNE
	if (initialize_Qe)
	{
		ept active_edgelist_newn = 0;
		for (ept j = 0; j < active_edgelist_n; j++)
			if (!deleted[active_edgelist[j]])
			{
				if (tri_cnt[active_edgelist[j]] < t_threshold)
					Qe[Qe_n++] = active_edgelist[j];
				else
					active_edgelist[active_edgelist_newn++] = active_edgelist[j];
			}
		active_edgelist_n = active_edgelist_newn;
	}
#endif

	// printf("%lu\n", Qe_n);

	ept deleted_edges_n = 0;
	ui Qv_idx = 0;
	while (Qv_idx < Qv_n || Qe_n)
	{
		if (Qe_n == 0)
		{
			// printf("hit\n");
			ui u = Qv[Qv_idx++]; // delete u from the graph due to have a degree < d_threshold
			ept u_n = pstart[u];
			for (ept k = pstart[u]; k < pend[u]; k++)
				if (!deleted[edgelist_pointer[k]])
				{
					edges[u_n] = edges[k];
					edgelist_pointer[u_n++] = edgelist_pointer[k];
					exists[edges[k]] = 1;
				}
			pend[u] = u_n;

			for (ept k = pstart[u]; k < pend[u]; k++)
				deleted[edgelist_pointer[k]] = 1;
			deleted_edges_n += pend[u] - pstart[u];
			degree[u] = 0;
			if (linear_heap != NULL)
				linear_heap->del(u);
			// printf("Removed %u\n", u);

			for (ept k = pstart[u]; k < pend[u]; k++)
			{
				ui v = edges[k];
				ept v_n = pstart[v];
				for (ept x = pstart[v]; x < pend[v]; x++)
					if (!deleted[edgelist_pointer[x]])
					{
						edges[v_n] = edges[x];
						edgelist_pointer[v_n++] = edgelist_pointer[x];
						if (edges[x] > v && exists[edges[x]])
						{
							if ((tri_cnt[edgelist_pointer[x]]--) == t_threshold)
								Qe[Qe_n++] = edgelist_pointer[x];
						}
					}
				pend[v] = v_n;

				if ((degree[v]--) == d_threshold)
				{
					Qv[Qv_n++] = v;
					if (v == critical_vertex)
					{
						for (ept k = pstart[u]; k < pend[u]; k++)
							exists[edges[k]] = 0;
						return 0;
					}
				}
				if (linear_heap != NULL)
					linear_heap->decrement(v, 1);
			}

			for (ept k = pstart[u]; k < pend[u]; k++)
				exists[edges[k]] = 0;
		}
#ifdef NO_TRUSS_PRUNE
		Qe_n = 0;
#endif
		for (ept j = 0; j < Qe_n; j++)
		{
			ept idx = Qe[j];
			ui u = edge_list[idx << 1], v = edge_list[(idx << 1) + 1];
			ui tri_n = tri_cnt[idx];
			// printf("remove %u %u\n", u, v);
			deleted[idx] = 1;
			if ((degree[u]--) == d_threshold)
			{
				Qv[Qv_n++] = u;
				if (u == critical_vertex)
					return 0;
			}
			if ((degree[v]--) == d_threshold)
			{
				Qv[Qv_n++] = v;
				if (v == critical_vertex)
					return 0;
			}
			// printf("before\n");
			if (linear_heap != NULL)
			{
				linear_heap->decrement(u, 1);
				linear_heap->decrement(v, 1);
			}
			// printf("after\n");
			deleted_edges_n++;

			if (degree[u] < degree[v])
				swap(u, v);
			// printf("here\n");

			if (degree[u] > degree[v] * 2)
			{ // binary search
				// if(false) {
				ept v_n = pstart[v], start = pstart[u];
				for (ept k = pstart[v]; k < pend[v]; k++)
					if (!deleted[edgelist_pointer[k]])
					{
						edges[v_n] = edges[k];
						edgelist_pointer[v_n++] = edgelist_pointer[k];

						if (tri_n && find(u, edges[k], start, pend[u], deleted, idx, edgelist_pointer, edges))
						{
							--tri_n;
							if ((tri_cnt[idx]--) == t_threshold)
								Qe[Qe_n++] = idx;
							if ((tri_cnt[edgelist_pointer[k]]--) == t_threshold)
								Qe[Qe_n++] = edgelist_pointer[k];
						}
					}
				pend[v] = v_n;
				assert(tri_n == 0);
			}
			else
			{ // sorted_merge
				ept ii = pstart[u], jj = pstart[v];
				ept u_n = pstart[u], v_n = pstart[v];

				while (true)
				{
					while (ii < pend[u] && deleted[edgelist_pointer[ii]])
						++ii;
					while (jj < pend[v] && deleted[edgelist_pointer[jj]])
						++jj;
					if (ii >= pend[u] || jj >= pend[v])
						break;

					if (edges[ii] == edges[jj])
					{
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];

						if ((tri_cnt[edgelist_pointer[ii]]--) == t_threshold)
							Qe[Qe_n++] = edgelist_pointer[ii];
						if ((tri_cnt[edgelist_pointer[jj]]--) == t_threshold)
							Qe[Qe_n++] = edgelist_pointer[jj];

						++ii;
						++jj;
					}
					else if (edges[ii] < edges[jj])
					{
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
						++ii;
					}
					else
					{
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];
						++jj;
					}
				}
				while (ii < pend[u])
				{
					if (!deleted[edgelist_pointer[ii]])
					{
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
					}
					++ii;
				}
				while (jj < pend[v])
				{
					if (!deleted[edgelist_pointer[jj]])
					{
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];
					}
					++jj;
				}
				pend[u] = u_n;
				pend[v] = v_n;
			}
		}
		Qe_n = 0;
	}
	return deleted_edges_n;
}
void Graph::ego_degen(ui n, ui m, ui *peel_sequence, ept *pstart, ui *edges, ui *degree, ui *rid, char *vis, ListLinearHeap *heap, ui *edgelist_pointer, bool output)
{
	Timer t;
	if (pend == nullptr)
		pend = new ept[n + 1];
	for (ui i = 0; i < n; i++)
		rid[peel_sequence[i]] = i;
	for (ui i = 0; i < n; i++)
	{
		ept &end = pend[i] = pstart[i];
		for (ept j = pstart[i]; j < pstart[i + 1]; j++)
			if (rid[edges[j]] > rid[i])
				edges[end++] = edges[j];
	}

	if (pend_buf == nullptr)
		pend_buf = new ept[n + 1];
	if (edgelist_pointer == nullptr)
		edgelist_pointer = new ui[m];
	ept *pstart_s = pend_buf;
	ui *pend_s = rid;
	ui *edges_s = edgelist_pointer;

	vector<ui> Q;
	vector<ui> vs;
	memset(vis, 0, sizeof(char) * n);
	for (ui i = n; i > 0; i--)
	{
		ui u = peel_sequence[i - 1];
		if (pend[u] - pstart[u] < kplex.size())
			continue;

		vs.clear();
		for (ui j = pstart[u]; j < pend[u]; j++)
		{
			vs.push_back(edges[j]);
			vis[edges[j]] = 1;
			degree[edges[j]] = 0;
		}
		for (ui j = 0; j < vs.size(); j++)
			for (ui k = pstart[vs[j]]; k < pend[vs[j]]; k++)
				if (vis[edges[k]])
				{
					++degree[vs[j]];
					++degree[edges[k]];
				}
		pend_s[vs[0]] = pstart_s[vs[0]] = 0;
		for (ui j = 1; j < vs.size(); j++)
			pend_s[vs[j]] = pstart_s[vs[j]] = pstart_s[vs[j - 1]] + degree[vs[j - 1]];
		for (ui j = 0; j < vs.size(); j++)
			for (ui k = pstart[vs[j]]; k < pend[vs[j]]; k++)
				if (vis[edges[k]])
				{
					edges_s[pend_s[vs[j]]++] = edges[k];
					edges_s[pend_s[edges[k]]++] = vs[j];
				}

		ui threshold = (kplex.size() > K ? kplex.size() - K : 0); // all vertices with degree < threshold can be pruned
		Q.clear();
		for (ui j = 0; j < vs.size(); j++)
			if (degree[vs[j]] < threshold)
			{
				Q.push_back(vs[j]);
				vis[vs[j]] = 0;
			}
		for (ui j = 0; j < Q.size(); j++)
			for (ui k = pstart_s[Q[j]]; k < pend_s[Q[j]]; k++)
				if (vis[edges_s[k]])
				{
					if ((degree[edges_s[k]]--) == threshold)
					{
						Q.push_back(edges_s[k]);
						vis[edges_s[k]] = 0;
					}
				}
		ui cnt = 0;
		for (ui j = 0; j < vs.size(); j++)
			if (vis[vs[j]])
				vs[cnt++] = vs[j];
		assert(cnt + Q.size() == vs.size());
		vs.resize(cnt);
		if (cnt == 0)
			continue;

		heap->init(vs.size(), vs.size() - 1, vs.data(), degree);
		bool found = false;
		for (ui ii = 0; ii < vs.size(); ii++)
		{
			ui v, key;
			heap->pop_min(v, key);
			if (found)
			{
				kplex.push_back(v);
				continue;
			}

			if (vs.size() - ii + 1 <= kplex.size())
				break;

			if (key + K >= vs.size() - ii)
			{
				kplex.clear();
				kplex.push_back(u);
				kplex.push_back(v);
				found = true;
				continue;
			}

			vis[v] = 0;
			for (ept j = pstart_s[v]; j < pend_s[v]; j++)
				if (vis[edges_s[j]])
					heap->decrement(edges_s[j], 1);
		}
		for (ui j = 0; j < vs.size(); j++)
			vis[vs[j]] = 0;
	}
	for (ui i = 0; i < n; i++)
		pend_buf[i] = pend[i];
	for (ui i = 0; i < n; i++)
		for (ept j = pstart[i]; j < pend[i]; j++)
			edges[pend_buf[edges[j]]++] = i;
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend_buf[i] == pstart[i + 1]);
#endif

	if (output)
		printf("*** EGo-Degen kPlex size: %lu, Time: %s (microseconds)\n", kplex.size(), Utility::integer_to_string(t.elapsed()).c_str());
}
