#ifndef _KPLEX_BB_MATRIX_
#define _KPLEX_BB_MATRIX_

#include "Utility.h"
#include "Timer.h"
// #define _SECOND_ORDER_PRUNING_

// pruning switches
#define S2RULE
#define ALTRB
// SR_BRANCHING can take values S_branching, R_branching, SR_branching
#define SR_BRANCHING S_branching
// if PART_BRANCH is false, then pivot branch gets executed...
// #define PART_BRANCH (false)
#define PART_BRANCH (K <= 5 && sparse)

// Upper bounding switches...
// #define SEESAW
// #define COLORBOUND
// #define PART_BOUND

#define CSIZE (R_end - S_end)
class KPLEX_BB_MATRIX
{
private:
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
		ui *peel_sequence;

	ui K;
	ui *best_solution;
	ui best_solution_size;
	ui _UB_;

	ui *neighbors;
	ui *nonneighbors;
	ui *S2;

	ui *SR;		// union of S and R, where S is at the front
	ui *SR_rid; // reverse ID for SR
	std::queue<ui> Qv;
	ui *level_id;
	ui max_level;

	std::vector<std::pair<ui, ui>> vp;
	std::vector<ui> non_adj;

	// std::vector<std::pair<ui,ui> > vp2;

	bool sparse;
	vector<ui> B, PI, PIMax, ISc, peelOrder, psz;
	ui *LPI;
	MBitSet bmp;
	ui sz1h;
	bool found_larger = false;
	bool ctcp_enabled = false;
	bool dense_search, forward_sol = false;

public:
	ui best_n_edges;
	ui nmkp = 0, color, part;
	KPLEX_BB_MATRIX(bool _ds = false)
	{
		n = 0;
		matrix = nullptr;
		matrix_size = 0;

#ifdef _SECOND_ORDER_PRUNING_
		cn = nullptr;
		removed_edges_n = 0;
#endif

		degree = degree_in_S = nullptr;

		best_solution = nullptr;
		K = best_solution_size = _UB_ = 0;

		neighbors = nonneighbors = nullptr;
		S2 = nullptr;

		SR = SR_rid = nullptr;
		level_id = nullptr;
		max_level = 0;
		best_n_edges = 0;
		dense_search = _ds;
	}

	~KPLEX_BB_MATRIX()
	{
		if (matrix != NULL)
		{
			delete[] matrix;
			matrix = NULL;
		}
#ifdef _SECOND_ORDER_PRUNING_
		if (cn != NULL)
		{
			delete[] cn;
			cn = NULL;
		}
#endif
		if (degree != NULL)
		{
			delete[] degree;
			degree = NULL;
		}
		if (degree_in_S != NULL)
		{
			delete[] degree_in_S;
			degree_in_S = NULL;
		}

		if (best_solution != NULL)
		{
			delete[] best_solution;
			best_solution = NULL;
		}
		if (SR != NULL)
		{
			delete[] SR;
			SR = NULL;
		}
		if (SR_rid != NULL)
		{
			delete[] SR_rid;
			SR_rid = NULL;
		}
		if (neighbors != NULL)
		{
			delete[] neighbors;
			neighbors = NULL;
		}
		if (nonneighbors != NULL)
		{
			delete[] nonneighbors;
			nonneighbors = NULL;
		}
		if (level_id != NULL)
		{
			delete[] level_id;
			level_id = NULL;
		}
	}

	void allocateMemory(ui n)
	{
		if (n <= 0)
			return;

		matrix_size = 1;
		// printf("matrix size: %lldMB\n", (((long long)UB)*v_n*4)/1024/1024);
		matrix = new char[matrix_size];
#ifdef _SECOND_ORDER_PRUNING_
		cn = new ui[matrix_size];
#endif

		degree = new ui[n];
		degree_in_S = new ui[n];
		best_solution = new ui[n];
		SR = new ui[n];
		SR_rid = new ui[n];
		neighbors = new ui[n];
		nonneighbors = new ui[n];
		peel_sequence = new ui[n];
		S2 = new ui[n];
		level_id = new ui[n];
		B.reserve(n);
		LPI = new ui[matrix_size];

		peelOrder.resize(n);
		psz.resize(n);
		PI.reserve(n);
		PIMax.reserve(n);
		ISc.reserve(n);
		bmp.init(n);
	}

	void load_graph(ui _n, const std::vector<std::pair<ui, ui>> &vp)
	{
		n = _n;
		if (((long long)n) * n > matrix_size)
		{
			do
			{
				matrix_size *= 2;
			} while (((long long)n) * n > matrix_size);
			delete[] matrix;
			matrix = new char[matrix_size];
			delete[] LPI;
			LPI = new ui[matrix_size];
#ifdef _SECOND_ORDER_PRUNING_
			delete[] cn;
			cn = new ui[matrix_size];
#endif
		}

#ifdef _SECOND_ORDER_PRUNING_
		memset(cn, 0, sizeof(ui) * ((long long)n) * n);
#endif
		sparse = 2.0 * vp.size() / n / (n - 1) < 0.9;
		memset(matrix, 0, sizeof(char) * ((long long)n) * n);
		for (ui i = 0; i < n; i++)
			degree[i] = 0;
		for (ui i = 0; i < vp.size(); i++)
		{
			assert(vp[i].first >= 0 && vp[i].first < n && vp[i].second >= 0 && vp[i].second < n);
			ui a = vp[i].first, b = vp[i].second;
			degree[a]++;
			degree[b]++;
			if (matrix[a * n + b])
				printf("Duplicate edge in KPLEX_BB_matrix.load_graph()\n");
			matrix[a * n + b] = matrix[b * n + a] = 1;
		}

#ifndef NDEBUG
		// printf("load graph of size n=%lld, m=%lu\n", n, vp.size());
		// for(ui i = 0;i < vp.size();i ++) printf("%d %d\n", vp[i].first, vp[i].second);
#endif
	}
	ui create_initial_kplex(ui R_end, bool must_include_0){
		ui S_end=0;
		fill(degree_in_S, degree_in_S+n, 0);
		if(must_include_0)
			move_to_S(S_end, R_end, 0);
		for(ui i=n-1; i>0; i--){
			ui u=peel_sequence[i];
			if(SR_rid[u]>=n or u==0) continue;
			if(can_move_to_S(S_end, R_end, u))
			move_to_S(S_end, R_end, u);
		}
		cout<<"S_end: "<<S_end<<" ";
		if(S_end > best_solution_size) store_solution(S_end);
		return S_end;
	}
	void kPlex(ui K_, ui UB_, std::vector<ui> &kplex, bool must_include_0)
	{
		K = K_;
		_UB_ = UB_;
		if (K <= 1)
		{
			printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
			return;
		}
		best_solution_size = kplex.size();
		ui n_edges = best_n_edges;
		ui R_end;
		initialization(R_end, must_include_0);

		// ui S_end = create_initial_kplex(R_end, must_include_0);

		if (R_end && best_solution_size < _UB_)
			BB_search(0, R_end, true, true);

		if (dense_search && best_n_edges > n_edges)
		{
			kplex.clear();
			for (int i = 0; i < best_solution_size + 1; i++)
				kplex.push_back(best_solution[i]);
		}
		if (best_solution_size > kplex.size())
		{
			kplex.clear();
			for (int i = 0; i < best_solution_size; i++)
				kplex.push_back(best_solution[i]);
		}
	}

	int main(int argc, char *argv[])
	{
		if (argc < 3)
		{
			printf("Usage: [1]exe [2]dir [3]k\n");
			return 0;
		}
		readGraph_binary(argv[1]);
		printf("Finish reading graph\n");
		K = atoi(argv[2]);
		if (K == 1)
		{
			printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
			return 0;
		}
		best_solution_size = 1;
		_UB_ = n;
		ui R_end;
		Timer t;
		initialization(R_end, false);
		if (R_end)
			BB_search(0, R_end, 1, 0);
		printf("Maximum %u-plex size: %u, time excluding reading: %s (micro seconds)\n", K, best_solution_size, Utility::integer_to_string(t.elapsed()).c_str());
		return 0;
	}

private:
	void readGraph_binary(char *dir)
	{
		FILE *f = Utility::open_file((std::string(dir) + std::string("/b_degree.bin")).c_str(), "rb");

		int tt;
		fread(&tt, sizeof(int), 1, f);
		if (tt != sizeof(int))
		{
			printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
			return;
		}
		ui m;
		fread(&n, sizeof(int), 1, f);
		fread(&m, sizeof(int), 1, f);
		printf("n = %u, m = %u\n", n, m / 2);

		allocateMemory(n);

		if (((long long)n) * n > matrix_size)
		{
			do
			{
				matrix_size *= 2;
			} while (((long long)n) * n > matrix_size);
			delete[] matrix;
			matrix = new char[matrix_size];
#ifdef _SECOND_ORDER_PRUNING_
			delete[] cn;
			cn = new ui[matrix_size];
#endif
		}

		fread(degree, sizeof(ui), n, f);
		fclose(f);

		f = Utility::open_file((std::string(dir) + std::string("/b_adj.bin")).c_str(), "rb");

#ifdef _SECOND_ORDER_PRUNING_
		memset(cn, 0, sizeof(ui) * n * n);
#endif
		memset(matrix, 0, sizeof(char) * n * n);
		for (ui i = 0; i < n; i++)
			if (degree[i] > 0)
			{
				fread(neighbors, sizeof(ui), degree[i], f);
				for (ui j = 0; j < degree[i]; j++)
					matrix[i * n + neighbors[j]] = 1;
				ui d = 0;
				for (ui j = 0; j < n; j++)
					if (matrix[i * n + j])
						++d;
				if (d != degree[i])
				{
					printf("%u may have duplicate edges\n", i);
					degree[i] = d;
				}
			}
		fclose(f);
	}

	void initialization(ui &R_end, bool must_include_0)
	{
		// the following computes a degeneracy ordering and a heuristic solution
		ui *core = nonneighbors;
		ui *vis = SR;
		memset(vis, 0, sizeof(ui) * n);
		ui max_core = 0, idx = n;
		for (ui i = 0; i < n; i++)
		{
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

			if (idx == n && min_degree + K >= n - i)
				idx = i;

			for (ui j = 0; j < n; j++)
				if (!vis[j] && matrix[u * n + j])
					--degree[j];
		}
		// _UB_=max_core+K;
		if (!dense_search && (n - idx > best_solution_size))
		{
			best_solution_size = n - idx;
			for (ui i = idx; i < n; i++)
				best_solution[i - idx] = peel_sequence[i];
			printf("Degen find a solution of size %u\n", best_solution_size);
		}

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

		if ((must_include_0 && SR_rid[0] == n) || best_solution_size >= _UB_)
		{
			R_end = 0;
			return;
		}

		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i];
			degree[u] = degree_in_S[u] = 0;
			for (ui j = 0; j < R_end; j++)
				if (matrix[u * n + SR[j]])
					++degree[u];
		}

		memset(level_id, 0, sizeof(ui) * n);
		for (ui i = 0; i < R_end; i++)
			level_id[SR[i]] = n;

		if (!Qv.empty())
			printf("!!! Something wrong. Qv must be empty in initialization\n");

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

		if (!remove_vertices_and_edges_with_prune(0, R_end, 0))
			R_end = 0;
	}

	void store_solution(ui size)
	{
		// if(!is_kplex(size)) return;
		if (size <= best_solution_size)
		{
			printf("!!! the solution to store is no larger than the current best solution!");
			return;
		}
		ui n_edges = 0;
		for (ui i = 0; i < size; i++)
		{
			if (forward_sol)
				n_edges += degree[SR[i]];
			else
				n_edges += degree_in_S[SR[i]];
		}

		// verification of kplex
		// for (ui i = 0; i < best_solution_size; i++)
		// {
		// 	if (forward_sol)
		// 	{
		// 		if (degree[SR[i]] + K < best_solution_size)
		// 			std::cout << degree[SR[i]] + K << " E: ";
		// 	}
		// 	else
		// 	{
		// 		if (degree_in_S[SR[i]] + K < best_solution_size)
		// 			std::cout << degree_in_S[SR[i]] + K << " E: ";
		// 	}
		// }

		forward_sol = false;
		nmkp++;
		printf("!!! BB_Search found a kplex of size: %u, n_edges: %u \n", size, n_edges);

		if (dense_search)
		{
			if (n_edges > best_n_edges)
			{
				best_n_edges = n_edges;
				for (ui i = 0; i < size; i++)
					best_solution[i] = SR[i];
			}
		}
		else
		{
			best_solution_size = size;
			found_larger = true;
			for (ui i = 0; i < best_solution_size; i++)
				best_solution[i] = SR[i];
			best_n_edges = n_edges;
		}
	}

	bool is_kplex(ui R_end)
	{
		for (ui i = 0; i < R_end; i++)
			if (degree[SR[i]] + K < R_end)
				return false;
		forward_sol = true;
		return true;
	}

	void BB_search(ui S_end, ui R_end, ui level, bool choose_zero, bool root_level = true)
	{
		if (S_end > best_solution_size)
			store_solution(S_end);
		if (R_end > best_solution_size && is_kplex(R_end))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + K >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + K > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + K == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#ifdef _SECOND_ORDER_PRUNING_
		for (ui i = 0; i < R_end; i++)
			for (ui j = i + 1; j < R_end; j++)
			{
				ui v = SR[i], w = SR[j];
				ui common_neighbors = 0;
				for (ui k = S_end; k < R_end; k++)
					if (matrix[SR[k] * n + v] && matrix[SR[k] * n + w])
						++common_neighbors;
				assert(cn[v * n + w] == common_neighbors);
				assert(cn[w * n + v] == common_neighbors);
			}
#endif
#endif

		ui old_removed_edges_n = 0, old_S_end = S_end, old_R_end = R_end;
		assert(Qv.empty());
#ifdef _SECOND_ORDER_PRUNING_
		while (!Qe.empty())
			Qe.pop();
		old_removed_edges_n = removed_edges_n;
#endif

		// choosing 0 as branching vertex
		if (choose_zero && SR_rid[0] < R_end && !move_u_to_S_with_prune(0, S_end, R_end, level))
		{
			// printf("here1\n");
			restore_SR_and_edges(S_end, R_end, old_S_end, old_R_end, level, old_removed_edges_n);
			return;
		}
#ifdef S2RULE
		ui S2_n = 0;
		for (ui i = 0; i < S_end; i++)
			if (R_end - degree[SR[i]] > K)
				S2[S2_n++] = SR[i];

		if (S2_n >= 2)
		{
			collect_removable_vertices_based_on_total_edges(S2_n, S_end, R_end, level);
			if (!remove_vertices_and_edges_with_prune(S_end, R_end, level))
			{
				restore_SR_and_edges(S_end, R_end, old_S_end, old_R_end, level, old_removed_edges_n);
				return;
			}
		}
#endif

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + K >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + K > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + K == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif

		// greedily add vertices to S
		if (!greedily_add_vertices_to_S(S_end, R_end, level))
		{
			// printf("here2\n");
			restore_SR_and_edges(S_end, R_end, old_S_end, old_R_end, level, old_removed_edges_n);
			return;
		}

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + K >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + K > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + K == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif

		if (S_end > best_solution_size)
			store_solution(S_end);
		if (R_end > best_solution_size && is_kplex(R_end))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
		{
			// printf("here3\n");
			restore_SR_and_edges(S_end, R_end, old_S_end, old_R_end, level, old_removed_edges_n);
			return;
		}
		ui beta = best_solution_size - S_end;
#ifdef PART_BOUND
		if (S_bound(S_end, R_end, level))
		{
			restore_SR_and_edges(S_end, R_end, old_S_end, old_R_end, level, old_removed_edges_n);
			return;
		}
#endif

#ifdef SEESAW
		if (seesawUB(S_end, R_end) <= best_solution_size)
		{
			// if (seesawUB(S_end, R_end)<=best_solution_size) {
			restore_SR_and_edges(S_end, R_end, old_S_end, old_R_end, level, old_removed_edges_n);
			return;
		}

#endif

#ifdef COLORBOUND
		// if (CSIZE>3*beta && seesawUB(S_end, R_end)<=best_solution_size) {
		if (colorUB(S_end, R_end) <= best_solution_size)
		{
			restore_SR_and_edges(S_end, R_end, old_S_end, old_R_end, level, old_removed_edges_n);
			return;
		}
#endif

// Adding kpex stuff here...
#ifdef ALTRB
		bounding.tick();
		if (!alt_RB(S_end, R_end, level))
		{
			// printf("here2\n");
			restore_SR_and_edges(S_end, R_end, old_S_end, old_R_end, level, old_removed_edges_n);
			return;
		}
		bounding.tock();
#endif

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + K >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + K > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + K == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif
		if (PART_BRANCH)
		{

			// ******************* Adding our branching stuff here...
			ui t_R_end = R_end;

			R_end = SR_BRANCHING(S_end, R_end, level);
			while (R_end < t_R_end)
			{
				// branching vertices are now in R_end to t_R_end, and they are already sorted in peelOrder
				// move branching vertex back to C
				ui u = SR[R_end];
				assert(level_id[u] == level && SR_rid[u] == R_end);
				R_end++;
				level_id[u] = n;
				char *t_matrix = matrix + u * n;
				degree[u] = degree_in_S[u] = 0;
				for (ui i = 0; i < R_end; i++)
					if (t_matrix[SR[i]])
					{
						ui w = SR[i];
						++degree[w];
						++degree[u];
						if (i < S_end)
							++degree_in_S[u];
					}

				if (best_solution_size >= _UB_)
					return;
				if (root_level)
					found_larger = false;
				if (found_larger)
					continue;

				ui pre_best_solution_size = best_solution_size, t_old_S_end = S_end, t_old_R_end = R_end, t_old_removed_edges_n = 0;
#ifdef _SECOND_ORDER_PRUNING_
				if (ctcp_enabled)
				{
					while (!Qe.empty())
						Qe.pop();
					t_old_removed_edges_n = removed_edges_n;
				}
#endif
				if (move_u_to_S_with_prune(u, S_end, R_end, level))
					BB_search(S_end, R_end, level + 1, false, false);
				restore_SR_and_edges(S_end, R_end, t_old_S_end, t_old_R_end, level, t_old_removed_edges_n);
			}

			restore_SR_and_edges(S_end, R_end, old_S_end, old_R_end, level, old_removed_edges_n);
		}
		else
		{ // pivot based branching
			if (B.empty() || SR_rid[B.back()] >= R_end || SR_rid[B.back()] < S_end)
				branch(S_end, R_end);
			ui u = B.back();
			B.pop_back();

			// First branch moves u to S
			ui pre_best_solution_size = best_solution_size, t_old_S_end = S_end, t_old_R_end = R_end, t_old_removed_edges_n = 0;
#ifdef _SECOND_ORDER_PRUNING_
			if (ctcp_enabled)
			{
				while (!Qe.empty())
					Qe.pop();
				t_old_removed_edges_n = removed_edges_n;
			}
#endif
			if (move_u_to_S_with_prune(u, S_end, R_end, level))
				BB_search(S_end, R_end, level + 1, false);

			// This version of 2nd branch restores u through remove_u function
			{
				restore_SR_and_edges(S_end, R_end, S_end, t_old_R_end, level, t_old_removed_edges_n);
				while (!Qv.empty())
				{
					ui v = Qv.front();
					Qv.pop();
					level_id[v] = n;
				}
				B.clear();
#ifdef _SECOND_ORDER_PRUNING_
				if (ctcp_enabled)
				{
					while (!Qe.empty())
						Qe.pop();
					t_old_removed_edges_n = removed_edges_n;
				}
#endif
				bool succeed = remove_u_from_S_with_prune(S_end, R_end, level);
				if (succeed && best_solution_size > pre_best_solution_size)
					succeed = collect_removable_vertices_and_edges(S_end, R_end, level);
				if (remove_vertices_and_edges_with_prune(S_end, R_end, level))
					BB_search(S_end, R_end, level + 1, false);
			}
			restore_SR_and_edges(S_end, R_end, old_S_end, old_R_end, level, old_removed_edges_n);
		}
	}

	void collect_removable_vertices_based_on_total_edges(ui S2_n, ui S_end, ui R_end, ui level)
	{
		vp.resize(R_end - S_end);
		ui max_nn = 0;
		for (ui i = S_end; i < R_end; i++)
		{
			ui nn = 0;
			if (S2_n != S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = 0; j < S2_n; j++)
					if (!t_matrix[S2[j]])
						++nn;
			}
			else
				nn = S_end - degree_in_S[SR[i]];
			if (nn > max_nn)
				max_nn = nn;
			vp[i - S_end].first = SR[i];
			vp[i - S_end].second = nn;
		}
		ui *cnt = neighbors;
		for (ui i = 0; i <= max_nn; i++)
			cnt[i] = 0;
		for (ui i = 0; i < vp.size(); i++)
			++cnt[vp[i].second];
		for (ui i = 0; i < max_nn; i++)
			cnt[i + 1] += cnt[i];
		for (ui i = max_nn; i > 0; i--)
			cnt[i] = cnt[i - 1];
		cnt[0] = 0;
		ui *ids = nonneighbors;
		for (ui i = 0; i < vp.size(); i++)
			ids[cnt[vp[i].second]++] = i;

		ui total_support = 0;
		for (ui i = 0; i < S2_n; i++)
			total_support += K - S_end + degree_in_S[S2[i]];

		ui new_n = 0;
		while (!Qv.empty())
			Qv.pop();
		for (ui i = 0; i < vp.size(); i++)
		{
			ui idx = ids[i], v = vp[ids[i]].first;
			ui t_support = total_support - vp[idx].second;
			char *t_matrix = matrix + v * n;
			ui j = 0, v_support = K - 1 - S_end + degree_in_S[v], ub = S_end + 1;
			// ISc.clear();
			while (true)
			{
				if (j == new_n)
					j = i + 1;
				if (j >= vp.size() || ub > best_solution_size || ub + vp.size() - j <= best_solution_size)
					break;
				ui u = vp[ids[j]].first, nn = vp[ids[j]].second;
				if (t_support < nn)
					break;
				if (t_matrix[u])
				{
					t_support -= nn;
					++ub;
					// ISc.push_back(u);
				}
				else if (v_support > 0)
				{
					--v_support;
					t_support -= nn;
					++ub;
					// ISc.push_back(u);
				}
				++j;
			}
			// if(ub > best_solution_size) {
			// 	ub=bound(S_end, R_end, v, ISc);
			// 	if(flag and ub<=best_solution_size) cout<<ub<<" "<<best_solution_size<<endl;
			// }
			if (ub <= best_solution_size)
			{
				level_id[v] = level;
				Qv.push(v);
			}
			else
				ids[new_n++] = ids[i];
		}
	}
	ui support(ui S_end, ui u)
	{
		return K - (S_end - degree_in_S[u]);
	}
	// ui bound(ui S_end, ui R_end, ui u, std::vector<ui>& R) {
	// 	char *t_matrix=matrix+u*n;
	// 	vp2.clear();
	// 	vp2.reserve(S_end);
	// 	for(ui i = 0;i < S_end;i ++) vp2.push_back(std::make_pair(support(S_end, SR[i]), SR[i]));
	// 	for(ui i=0;i<S_end;i++)if(!t_matrix[SR[i]])vp2[i].first--;
	// 	// for(ui i = 0;i < S_end;i ++) vp.push_back(std::make_pair(-(degree_in_S[SR[i]]-neiInP[SR[i]]), SR[i]));
	// 	sort(vp2.begin(), vp2.end());
	// 	ui UB = S_end+1, cursor = 0;
	// 	for(ui i = 0;i < (ui)vp2.size(); i++) {
	// 		ui u = vp2[i].second;
	// 		if(vp2[i].first == 0) continue;// boundary vertex
	// 		ui count = 0;
	// 		char *t_matrix = matrix + u*n;
	// 		for(ui j = cursor;j < R.size();j ++) if(!t_matrix[R[j]]) {
	// 			if(j != cursor + count) std::swap(R[j], R[cursor+count]);
	// 			++ count;
	// 		}
	// 		UB += std::min(count, vp2[i].first);

	// 		if(UB <= best_solution_size)
	// 			cursor += count;
	// 		else
	// 			return UB;
	// 	}
	// 	UB+=(R.size()-cursor);
	// 	return UB;
	// }
	pair<ui, ui> partition_left_right(ui S_end, ui R_end, ui level)
	{
		ui S_l = 0, R_l = S_end;
		ui u_s, nnc_s;
		ui *nn_u_s = neighbors, *nn_u = nonneighbors;

		while (true)
		{
			double nnr_max = 1;
			for (ui i = S_l; i < S_end; i++)
			{
				ui u = SR[i], nnc = 0;
				char *t_matrix = matrix + u * n;
				if (support(S_end, u) == 0)
					continue;
				for (ui i = R_l; i < R_end; i++)
					if (!t_matrix[SR[i]])
						nn_u[nnc++] = SR[i];
				double nnr = (double)nnc / support(S_end, u);
				if (nnr > nnr_max)
				{
					nnr_max = nnr;
					u_s = u;
					std::swap(nn_u, nn_u_s);
					nnc_s = nnc;
				}
			}
			if (nnr_max == 1)
				break;
			for (ui i = 0; i < nnc_s; i++)
				swap_pos(SR_rid[nn_u_s[i]], R_l++);
			swap_pos(SR_rid[u_s], S_l++);
		}
		return {S_l, R_l};
	}
	enum Partition
	{
		left,
		right
	};
	bool alt_reduction_rules(ui S_end, ui &R_end, ui &R_l, ui LB_l, ui LB_r, Partition side, ui level)
	{
		ui Ra, Rb, R_end_temp = R_end;
		if (side == left)
			Ra = S_end, Rb = R_l;
		else
			Ra = R_l, Rb = R_end;

		while (!Qv.empty())
			Qv.pop();

		for (ui i = Ra; i < Rb; i++)
		{
			ui u = SR[i];
			char *t_matrix = matrix + u * n;
			ui neigh_l = degree_in_S[u], neigh_r = degree_in_S[u];
			for (ui j = S_end; j < R_end_temp; j++)
			{
				if (t_matrix[SR[j]])
					if (j < R_l)
						neigh_l++;
					else
						neigh_r++;
			}
			if (side == left)
			{
				if (neigh_l + K < LB_l + S_end or neigh_r + K < LB_r + S_end + 1)
				{
					swap_pos(i, --R_l);
					swap_pos(R_l, --R_end_temp);
					Rb--;
					level_id[u] = level;
					Qv.push(u);
				}
			}
			else
			{
				if (neigh_l + K < LB_l + S_end + 1 or neigh_r + K < LB_r + S_end)
				{
					swap_pos(i, --R_end_temp);
					Rb--;
					level_id[u] = level;
					Qv.push(u);
				}
			}
		}
		return remove_vertices_and_edges_with_prune_left(S_end, R_end, R_l, level);
	}
	bool can_move_to_S(ui S_end, ui R_end, ui u)
	{
		if (support(S_end, u) == 0)
			return false;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < S_end; i++)
		{
			if (!t_matrix[SR[i]] && support(S_end, SR[i]) == 0)
				return false; // SR[i] is non-neighbor of u in S and it's saturated
		}
		return true;
	}

	void move_to_S(ui &S_end, ui R_end, ui u)
	{
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
		{
			if (t_matrix[SR[i]])
				degree_in_S[SR[i]]++;
		}
		swap_pos(SR_rid[u], S_end++);
	}
	bool RR2(ui &S_end, ui R_l, ui &R_end, ui level, Partition side)
	{
		ui old_S_end = S_end;
		bool ret = true;
		ui Ra, Rb;
		if (side == left)
			Ra = S_end, Rb = R_l;
		else
			Ra = R_l, Rb = R_end;
		for (ui i = Ra; i < Rb; i++)
		{
			if (!can_move_to_S(S_end, R_end, SR[i]))
			{
				ret = false;
				break;
			}
			move_to_S(S_end, R_end, SR[i]);
		}

		if (ret)
			ret = prune_R(S_end, R_end, level);
		return ret;
	}

	bool prune_R(ui S_end, ui &R_end, ui level)
	{
		while (!Qv.empty())
			Qv.pop();
		// prune vertices in remaining R
		for (ui i = S_end; i < R_end; i++)
			if (!can_move_to_S(S_end, R_end, SR[i]))
				Qv.push(SR[i]);
		return remove_vertices_and_edges_with_prune(S_end, R_end, level);
	}

	bool alt_RB(ui &S_end, ui &R_end, ui level)
	{

		ui old_S_end = S_end;
		// Saving S, as it should be recovered with same location when returning from this function.
		for (ui i = 0; i < old_S_end; i++)
			S2[i] = SR[i];

		auto SR_left = partition_left_right(S_end, R_end);
		bool ret = true;
		ui S_l = SR_left.first, R_l = SR_left.second;
		ui UB_l = R_l - S_end, LB_l = 0;

		while (true)
		{
			ui UB = compute_UB(S_end, R_end, S_l, R_l, left);
			if (UB == UB_l)
				break;
			else
				UB_l = UB;
			ui LB_r = S_end + UB_l < best_solution_size + 1 ? best_solution_size + 1 - S_end - UB_l : 0;

			// RR1 on C_r
			if (!alt_reduction_rules(S_end, R_end, R_l, LB_l, LB_r, right, level))
			{
				ret = false;
				break;
			}

			ui UB_r = compute_UB(S_end, R_end, S_l, R_l, right);
			LB_l = (S_end + UB_r < best_solution_size + 1) ? best_solution_size + 1 - S_end - UB_r : 0;

			//  RR1 on C_l
			if (!alt_reduction_rules(S_end, R_end, R_l, LB_l, LB_r, left, level))
			{
				ret = false;
				break;
			}

			if (UB_r + UB_l + S_end == best_solution_size + 1)
			{
				// RR2 on C_l
				if (UB_l == R_l - S_end)
				{
					ret = RR2(S_end, R_l, R_end, level, left);
					break;
				}
				// RR2 on C_r
				else if (UB_r == R_end - R_l)
				{
					ret = RR2(S_end, R_l, R_end, level, right);
					break;
				}
			}
		}
		// restoring S
		for (ui i = 0; i < old_S_end; i++)
			SR[i] = S2[i], SR_rid[S2[i]] = i;
		return ret;
	}
	bool greedily_add_vertices_to_S(ui &S_end, ui &R_end, ui level)
	{
		while (true)
		{
			ui *candidates = S2;
			ui candidates_n = 0;
			ui skip;
			if (dense_search)
				skip = 1;
			else
				skip = 2;
			for (ui i = S_end; i < R_end; i++)
			{
				ui u = SR[i];
				if (R_end - degree[u] > K)
					continue;

				if (degree[u] >= R_end - skip)
				{
					candidates[candidates_n++] = u;
					continue;
				}

				char *t_matrix = matrix + u * n;
				bool OK = true;
				for (ui j = 0; j < R_end; j++)
					if (j != i && !t_matrix[SR[j]] && R_end - degree[SR[j]] > K)
					{
						OK = false;
						break;
					}
				if (OK)
					candidates[candidates_n++] = u;
			}

			if (!candidates_n)
				break;

			while (candidates_n)
			{
				ui u = candidates[--candidates_n];
				assert(SR_rid[u] >= S_end);
				if (SR_rid[u] >= R_end)
					return false;

				if (!move_u_to_S_with_prune(u, S_end, R_end, level))
					return false;
			}
		}
		return true;
	}

	void check_degrees(ui S_end, ui R_end)
	{
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
	}

	bool greedily_add_nonneighbors(ui *candidates, ui candidates_n, ui &S_end, ui &R_end, ui level)
	{
		while (candidates_n)
		{
			ui u = candidates[--candidates_n];
			assert(SR_rid[u] >= S_end);
			if (SR_rid[u] >= R_end || !move_u_to_S_with_prune(u, S_end, R_end, level))
				return false;
		}
		return true;
	}

	void get_neighbors_and_nonneighbors(ui u, ui R_end, ui &neighbors_n, ui &nonneighbors_n)
	{
		neighbors_n = 0;
		nonneighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
			if (SR[i] != u)
			{
				if (t_matrix[SR[i]])
					neighbors[neighbors_n++] = SR[i];
				else
					nonneighbors[nonneighbors_n++] = SR[i];
			}
	}

	bool move_u_to_S_with_prune(ui u, ui &S_end, ui &R_end, ui level)
	{
		assert(SR_rid[u] >= S_end && SR_rid[u] < R_end && SR[SR_rid[u]] == u);
		assert(degree_in_S[u] + K > S_end);
#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + K >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + K > S_end);
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif
		if (SR_rid[u] != S_end)
			swap_pos(S_end, SR_rid[u]);
		++S_end;

		ui neighbors_n = 0, nonneighbors_n = 0;
		get_neighbors_and_nonneighbors(u, R_end, neighbors_n, nonneighbors_n);
		assert(neighbors_n + nonneighbors_n == R_end - 1);
		for (ui i = 0; i < neighbors_n; i++)
			++degree_in_S[neighbors[i]];
#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + K >= S_end);
#endif

		while (!Qv.empty())
			Qv.pop();
		// reduction rules based on the fact that each vertex can have at most k-1 nonneighbors
		for (ui i = 0; i < nonneighbors_n; i++)
		{
			ui v = nonneighbors[i];
			if (SR_rid[v] >= S_end)
			{
				if (level_id[v] == level)
					continue;
				if (S_end - degree_in_S[v] >= K || S_end - degree_in_S[u] == K)
				{
					level_id[v] = level;
					Qv.push(v);
				}
			}
			else if (S_end - degree_in_S[v] == K)
			{
				char *tt_matrix = matrix + v * n;
				for (ui j = S_end; j < R_end; j++)
					if (level_id[SR[j]] > level && !tt_matrix[SR[j]])
					{
						level_id[SR[j]] = level;
						Qv.push(SR[j]);
					}
			}
		}

#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + K == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(level_id[SR[j]] == level || t_matrix[SR[j]]);
			}
#endif

#ifdef _SECOND_ORDER_PRUNING_
		for (ui i = 0; i < nonneighbors_n; i++)
		{
			int v = nonneighbors[i];
			if (SR_rid[v] < S_end || level_id[v] == level)
				continue;
			if (upper_bound_based_prune(S_end, u, v))
			{
				level_id[v] = level;
				Qv.push(v);
			}
		}

		// update cn(.,.) for pairs of u's neighbors
		for (ui i = 0; i < neighbors_n; i++)
			for (ui j = i + 1; j < neighbors_n; j++)
			{
				assert(cn[neighbors[i] * n + neighbors[j]]);
				--cn[neighbors[i] * n + neighbors[j]];
				--cn[neighbors[j] * n + neighbors[i]];
			}

		while (!Qe.empty())
			Qe.pop();
		ui new_n = 0;
		for (ui i = 0; i < nonneighbors_n; i++)
			if (level_id[nonneighbors[i]] > level)
				nonneighbors[new_n++] = nonneighbors[i];
		nonneighbors_n = new_n;
		for (ui i = 1; i < nonneighbors_n; i++)
		{ // process pairs of u's non-neighbors
			ui w = nonneighbors[i];
			for (ui j = 0; j < i; j++)
			{
				ui v = nonneighbors[j];
				if (!upper_bound_based_prune(S_end, v, w))
					continue;
				assert(SR_rid[w] > SR_rid[v]);
				if (SR_rid[w] < S_end)
					return false; // v, w \in S --- pruned
				else if (SR_rid[v] >= S_end)
				{ // v, w, \in R --- the edge (v,w) can be removed
					if (matrix[v * n + w])
						Qe.push(std::make_pair(v, w));
				}
				else
				{
					assert(level_id[w] > level);
					level_id[w] = level;
					Qv.push(w);
					break;
				}
			}
		}
#endif

		return remove_vertices_and_edges_with_prune(S_end, R_end, level);
	}

	bool remove_vertices_and_edges_with_prune(ui S_end, ui &R_end, ui level)
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
					return false;
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
					return false;
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

		return true;
	}

	bool remove_vertices_and_edges_with_prune_left(ui S_end, ui &R_end, ui &R_l, ui level)
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
				if (SR_rid[u] < R_l)
					swap_pos(SR_rid[u], --R_l);
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
					return false;
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
					return false;
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
		return true;
	}

	void restore_SR_and_edges(ui &S_end, ui &R_end, ui old_S_end, ui old_R_end, ui level, ui old_removed_edges_n)
	{
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
#endif
		while (!Qv.empty())
		{
			ui u = Qv.front();
			Qv.pop();
			assert(level_id[u] == level);
			assert(SR_rid[u] < R_end);
			level_id[u] = n;
		}

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			if (level_id[SR[i]] <= level)
				printf("level_id[%u] = %u, level = %u, n = %u\n", SR[i], level_id[SR[i]], level, n);
			assert(level_id[SR[i]] > level);
		}
#endif

		for (; R_end < old_R_end; R_end++)
		{ // insert u back into R
			ui u = SR[R_end];
			assert(level_id[u] == level && SR_rid[u] == R_end);
			level_id[u] = n;

			ui neighbors_n = 0;
			char *t_matrix = matrix + u * n;
			degree[u] = degree_in_S[u] = 0;

			for (ui i = 0; i < R_end; i++)
				if (t_matrix[SR[i]])
				{
					ui w = SR[i];
					neighbors[neighbors_n++] = w;
					++degree[w];
					++degree[u];
					if (i < S_end)
						++degree_in_S[u];
				}
#ifdef _SECOND_ORDER_PRUNING_
			for (ui i = 0; i < neighbors_n; i++)
			{
				ui v = neighbors[i];
				for (ui j = i + 1; j < neighbors_n; j++)
				{
					ui w = neighbors[j];
					++cn[v * n + w];
					++cn[w * n + v];
				}
			}
			ui *t_cn = cn + u * n;
			for (ui i = 0; i < R_end; i++)
				t_cn[SR[i]] = 0;
			for (ui i = 0; i < neighbors_n; i++)
				if (SR_rid[neighbors[i]] >= S_end)
				{
					ui v = neighbors[i];
					char *t_matrix = matrix + v * n;
					for (ui j = 0; j < R_end; j++)
						if (t_matrix[SR[j]])
							++t_cn[SR[j]];
				}
			for (ui i = 0; i < R_end; i++)
			{
				cn[SR[i] * n + u] = t_cn[SR[i]];
#ifndef NDEBUG
				ui common_neighbors = 0, v = SR[i], w = u;
				for (ui k = S_end; k < R_end; k++)
					if (matrix[SR[k] * n + v] && matrix[SR[k] * n + w])
						++common_neighbors;
				if (t_cn[SR[i]] != common_neighbors)
					printf("t_cn[SR[i]] = %u, comon_neighbors = %u\n", t_cn[SR[i]], common_neighbors);
				assert(t_cn[SR[i]] == common_neighbors);
#endif
			}
#endif
		}

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
#endif

		for (; S_end > old_S_end; S_end--)
		{ // move u from S to R
			ui u = SR[S_end - 1];
			assert(SR_rid[u] == S_end - 1);

			ui neighbors_n = 0;
			char *t_matrix = matrix + u * n;
			for (ui i = 0; i < R_end; i++)
				if (t_matrix[SR[i]])
				{
					ui w = SR[i];
					neighbors[neighbors_n++] = w;
					--degree_in_S[w];
				}
#ifdef _SECOND_ORDER_PRUNING_
			for (ui i = 0; i < neighbors_n; i++)
			{
				ui v = neighbors[i];
				for (ui j = i + 1; j < neighbors_n; j++)
				{
					ui w = neighbors[j];
					++cn[v * n + w];
					++cn[w * n + v];
				}
			}
			ui *t_cn = cn + u * n;
			for (ui i = 0; i < R_end; i++)
				t_cn[SR[i]] = 0;
			for (ui i = 0; i < neighbors_n; i++)
				if (SR_rid[neighbors[i]] >= S_end)
				{
					ui v = neighbors[i];
					char *t_matrix = matrix + v * n;
					for (ui j = 0; j < R_end; j++)
						if (t_matrix[SR[j]])
							++t_cn[SR[j]];
				}
			for (ui i = 0; i < R_end; i++)
			{
				cn[SR[i] * n + u] = t_cn[SR[i]];
#ifndef NDEBUG
				ui common_neighbors = 0, v = SR[i], w = u;
				for (ui k = S_end; k < R_end; k++)
					if (matrix[SR[k] * n + v] && matrix[SR[k] * n + w])
						++common_neighbors;
				if (t_cn[SR[i]] != common_neighbors)
					printf("t_cn[SR[i]] = %u, comon_neighbors = %u\n", t_cn[SR[i]], common_neighbors);
				assert(t_cn[SR[i]] == common_neighbors);
#endif
			}
#endif
		}

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
#endif

#ifdef _SECOND_ORDER_PRUNING_
		for (ui i = old_removed_edges_n; i < removed_edges_n; i++)
		{ // insert edge back into matrix
			ui v = removed_edges[i].first, w = removed_edges[i].second;
			assert(SR_rid[v] >= S_end && SR_rid[v] < R_end && SR_rid[w] >= S_end && SR_rid[w] < R_end);
			if (matrix[v * n + w])
				continue;

#ifndef NDEBUG
			// printf("restore edge between %u and %u\n", v, w);
#endif
			matrix[v * n + w] = matrix[w * n + v] = 1;
			++degree[v];
			++degree[w];

			char *t_matrix = matrix + v * n;
			for (ui i = 0; i < R_end; i++)
				if (t_matrix[SR[i]])
				{
					++cn[w * n + SR[i]];
					++cn[SR[i] * n + w];
				}
			t_matrix = matrix + w * n;
			for (ui i = 0; i < R_end; i++)
				if (t_matrix[SR[i]])
				{
					++cn[v * n + SR[i]];
					++cn[SR[i] * n + v];
				}
		}
		removed_edges_n = old_removed_edges_n;
#endif

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + K >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(old_S_end == S_end || degree_in_S[SR[i]] + K > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + K == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(old_S_end == S_end || t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif
	}

	bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level)
	{
		assert(S_end);
		ui u = SR[S_end - 1];
		--S_end;
		--R_end;
		swap_pos(S_end, R_end);
		level_id[u] = level;

		bool terminate = false;
		ui neighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
			if (t_matrix[SR[i]])
				neighbors[neighbors_n++] = SR[i];
		for (ui i = 0; i < neighbors_n; i++)
		{
			ui v = neighbors[i];
			--degree_in_S[v];
			--degree[v];
			if (degree[v] + K <= best_solution_size)
			{
				if (SR_rid[v] < S_end)
					terminate = true;
				else
				{
					assert(level_id[v] > level);
					level_id[v] = level;
					Qv.push(v);
				}
			}
		}
		if (terminate)
			return false;

#ifdef _SECOND_ORDER_PRUNING_
		for (ui i = 1; i < neighbors_n; i++)
			if (level_id[neighbors[i]] > level)
			{
				ui w = neighbors[i];
				for (ui j = 0; j < i; j++)
				{
					ui v = neighbors[j];
					if (!upper_bound_based_prune(S_end, v, w))
						continue;

					assert(SR_rid[v] < SR_rid[w]);
					if (SR_rid[w] < S_end)
						return false; // v, w \in S
					else if (SR_rid[v] >= S_end)
					{ // v, w, \in R
						if (matrix[v * n + w])
							Qe.push(std::make_pair(v, w));
					}
					else
					{
						assert(level_id[w] > level);
						level_id[w] = level;
						Qv.push(w);
						break;
					}
				}
			}
#endif
		return true;
	}

	bool collect_removable_vertices_and_edges(ui S_end, ui R_end, ui level)
	{
		for (ui i = 0; i < S_end; i++)
			if (degree[SR[i]] + K <= best_solution_size)
				return false;

#ifdef _SECOND_ORDER_PRUNING_
		for (ui i = 0; i < S_end; i++)
			for (ui j = i + 1; j < S_end; j++)
			{
				if (upper_bound_based_prune(S_end, SR[i], SR[j]))
					return false;
			}
#endif

		for (ui i = S_end; i < R_end; i++)
			if (level_id[SR[i]] > level)
			{
				ui v = SR[i];
				if (S_end - degree_in_S[v] >= K || degree[v] + K <= best_solution_size)
				{
					level_id[v] = level;
					Qv.push(v);
					continue;
				}
				char *t_matrix = matrix + v * n;
				for (ui j = 0; j < S_end; j++)
				{
#ifdef _SECOND_ORDER_PRUNING_
					if ((S_end - degree_in_S[SR[j]] == K && !t_matrix[SR[j]]) || upper_bound_based_prune(S_end, v, SR[j]))
#else
					if (S_end - degree_in_S[SR[j]] == K && !t_matrix[SR[j]])
#endif
					{
						level_id[v] = level;
						Qv.push(v);
						break;
					}
				}
			}

#ifdef _SECOND_ORDER_PRUNING_
		for (ui i = S_end; i < R_end; i++)
			if (level_id[SR[i]] > level)
			{
				for (ui j = i + 1; j < R_end; j++)
					if (level_id[SR[j]] > level && matrix[SR[i] * n + SR[j]])
					{
						if (upper_bound_based_prune(S_end, SR[i], SR[j]))
							Qe.push(std::make_pair(SR[i], SR[j]));
					}
			}
#endif

		return true;
	}

#ifdef _SECOND_ORDER_PRUNING_
	bool upper_bound_based_prune(ui S_end, ui u, ui v)
	{
		// ui ub = S_end + 2*K - (S_end - degree_in_S[u]) - (S_end - degree_in_S[v]) + cn[u*n + v];
		ui ub = 2 * K + degree_in_S[u] - S_end + degree_in_S[v] + cn[u * n + v];
		if (SR_rid[u] >= S_end)
		{
			--ub; // S_end ++
			if (matrix[u * n + v])
				++ub; // degree_in_S[v] ++
		}
		if (SR_rid[v] >= S_end)
		{
			--ub;
			if (matrix[v * n + u])
				++ub;
		}
		return ub <= best_solution_size;
	}
#endif

	void swap_pos(ui i, ui j)
	{
		std::swap(SR[i], SR[j]);
		SR_rid[SR[i]] = i;
		SR_rid[SR[j]] = j;
	}

	ui choose_branch_vertex(ui S_end, ui R_end)
	{
		ui *D = neighbors;
		ui D_n = 0;
		for (ui i = 0; i < R_end; i++)
			if (R_end - degree[SR[i]] > K)
				D[D_n++] = SR[i];
		assert(D_n != 0);

		ui min_degree_in_S = n;
		for (ui i = 0; i < D_n; i++)
			if (degree_in_S[D[i]] < min_degree_in_S)
				min_degree_in_S = degree_in_S[D[i]];
		ui u = n, min_degree = n;
		for (ui i = 0; i < D_n; i++)
			if (degree_in_S[D[i]] == min_degree_in_S && degree[D[i]] < min_degree)
			{
				min_degree = degree[D[i]];
				u = D[i];
			}
		assert(u != n);

		if (SR_rid[u] < S_end)
		{
			ui max_degree = 0, b = n;
			char *t_matrix = matrix + u * n;
			for (ui i = S_end; i < R_end; i++)
				if (!t_matrix[SR[i]] && degree[SR[i]] > max_degree)
				{
					max_degree = degree[SR[i]];
					b = SR[i];
				}
			return b;
		}
		else if (degree_in_S[u] < S_end || R_end - degree[u] > K + 1)
			return u;
		else
		{
			ui max_degree = 0, w = n;
			for (ui i = S_end; i < R_end; i++)
				if (degree[SR[i]] > max_degree)
				{
					max_degree = degree[SR[i]];
					w = SR[i];
				}
			if (degree[w] + 1 >= R_end)
			{
				printf("!!! WA degree[w]: %u, R_end: %u\n", degree[w], R_end);
			}
			assert(degree[w] + 1 < R_end);
			if (R_end - degree[w] == 2)
				return w;
			char *t_matrix = matrix + w * n;
			for (ui i = S_end; i < R_end; i++)
				if (!t_matrix[SR[i]] && R_end - degree[SR[i]] == K + 1)
					return SR[i];
		}

		printf("!!! WA in choose_branch_vertex\n");
		return n;
	}

	void branch(ui S_end, ui R_end)
	{
		B.clear();
		ui minnei = 0x3f3f3f3f;
		ui pivot; // should it be 0xffffffff?
		char *t_matrix = matrix + 0 * n;
		for (ui i = S_end; i < R_end; i++)
		{
			ui v = SR[i];
			if (!t_matrix[v] && // HOP2 first
				(degree[v] + K <= best_solution_size + 1 ||
				 support(S_end, v) == 1 ||
				 support(S_end, 0) == 1))
			{
				B.push_back(v);
				return;
			}
			if (degree[v] < minnei)
			{
				minnei = degree[v];
				pivot = v;
			}
		}

		t_matrix = matrix + pivot * n;
		for (ui i = S_end; i < R_end; i++)
			if (!t_matrix[SR[i]])
				B.push_back(SR[i]);

		auto comp = [&](int a, int b)
		{ return degree[a] > degree[b]; };
		std::sort(B.begin(), B.end(), comp);
	}

	ui SR_branching(ui S_end, ui R_end, ui level)
	{
		ui cend = R_end;
		ui beta = best_solution_size - S_end;
		ui ub = 0;
		bool flag = true;
		do
		{
			double ubp = 0;
			if (flag)
			{
				ubp = tryPartition(S_end, cend);
				if (ubp == 0)
					flag = false;
			}
			double ubc = tryColor(S_end, cend);

			if (ubp == 0 or
				(ISc.size() / ubc > PIMax.size() / ubp) or
				((ISc.size() / ubc == PIMax.size() / ubp) and (ISc.size() > PIMax.size())))

			{
				if (ubc <= beta)
					beta -= ubc;
				else
					break;
				for (ui v : ISc)
					swap_pos(v, --cend);
			}
			else
			{
				if (ubp <= beta)
					beta -= ubp;
				else
					break;
				for (ui v : PIMax)
					swap_pos(v, --cend);
			}
		} while (beta > 0 && cend > S_end);
		if (beta > 0)
			cend -= min(beta, cend - S_end);
		return move_candidates_to_end(S_end, cend, R_end, level);
	}
	ui R_branching(ui S_end, ui R_end, ui level)
	{

		ui cend = R_end;
		ui beta = best_solution_size - S_end;
		while (beta > 0 && cend > S_end)
		{
			ui ub = tryColor(S_end, cend);
			if (ub <= beta)
			{
				for (ui i : ISc)
				{
					swap_pos(i, --cend);
				}
				beta -= ub;
			}
			else
				break;
		}
		if (beta > 0)
			cend -= min(beta, cend - S_end);

		return move_candidates_to_end(S_end, cend, R_end, level);
	}

	ui S_branching(ui S_end, ui R_end, ui level)
	{
		for (ui i = 0; i < S_end; i++)
		{
			ui u = SR[i];
			psz[i] = 0;
			if (support(S_end, u) == 0)
				continue;
			// skipping it, because this is a boundary vertex, and it can't have any non-neighbor candidate
			ui *t_LPI = LPI + i * n;
			for (ui j = S_end; j < R_end; j++)
			{
				ui v = SR[j];
				if (!matrix[u * n + v])
					// PI[u].push_back(v);
					t_LPI[psz[i]++] = v;
			}
		}
		ui beta = best_solution_size - S_end;
		ui cend = R_end;

		while (S_end < cend)
		{
			ui maxpi = -1;
			double maxdise = 0;
			for (ui i = 0; i < S_end; i++)
			{
				ui u = SR[i];
				if (psz[i] == 0)
					continue;
				double cost = min(support(S_end, u), psz[i]);
				double dise = psz[i] / cost;
				if (cost <= beta and dise > maxdise)
					maxpi = i, maxdise = dise;
			}
			if (maxpi == -1)
				break;
			else
			{
				// remove pi* from C
				ui *t_LPI = LPI + maxpi * n;
				for (ui i = 0; i < psz[maxpi]; i++)
				{
					// removing from C
					ui v = t_LPI[i];
					swap_pos(SR_rid[v], --cend);
				}
				// beta-=cost(pi*)
				beta -= min(support(S_end, SR[maxpi]), psz[maxpi]);

				// remove maxpi from every pi
				psz[maxpi] = 0;
				for (ui i = 0; i < S_end; i++)
				{
					ui j = 0;
					ui *t_LPI = LPI + i * n;
					for (ui k = 0; k < psz[i]; k++)
						if (SR_rid[t_LPI[k]] < cend)
							t_LPI[j++] = t_LPI[k];
					psz[i] = j;
				}
			}
			if (beta == 0)
				break;
		}
		if (beta > 0)
			cend -= min(beta, cend - S_end);
		return move_candidates_to_end(S_end, cend, R_end, level);
	}
	ui compute_UB(ui S_end, ui R_end, ui S_l, ui R_l, Partition side)
	{
		// when computing for left side
		if(side==right) return color_ub_partial_C(R_l, R_end);
		ui Sa, Sb, Ra, Rb;
		if (side == left)
			Sa = 0, Sb = S_l, Ra = S_end, Rb = R_l;
		else
			Sa = S_l, Sb = S_end, Ra = R_l, Rb = R_end;

		for (ui i = Sa; i < Sb; i++)
		{
			ui u = SR[i];
			psz[i] = 0;
			if (support(S_end, u) == 0)
				continue;
			// skipping it, because this is a boundary vertex, and it can't have any non-neighbor candidate
			ui *t_LPI = LPI + i * n;
			for (ui j = Ra; j < Rb; j++)
			{
				ui v = SR[j];
				if (!matrix[u * n + v])
					t_LPI[psz[i]++] = v;
			}
		}

		ui cend = Rb;
		ui UB = 0;
		while (true)
		{
			ui maxpi = -1, mincost = 0;
			double maxdise = 0;
			for (ui i = Sa; i < Sb; i++)
			{
				ui u = SR[i];
				if (psz[i] == 0)
					continue;
				double cost = min(support(S_end, u), psz[i]);
				double dise = psz[i] / cost;
				if (dise > maxdise)
					maxpi = i, maxdise = dise, mincost = cost;
			}
			if (maxpi == -1)
				break;
			// remove pi* from C
			UB += mincost;
			ui *t_LPI = LPI + maxpi * n;
			for (ui i = 0; i < psz[maxpi]; i++)
			{
				// removing from C
				ui v = t_LPI[i];
				swap_pos(SR_rid[v], --cend);
			}

			// remove maxpi from every pi
			psz[maxpi] = 0;
			for (ui i = Sa; i < Sb; i++)
			{
				ui j = 0;
				ui *t_LPI = LPI + i * n;
				for (ui k = 0; k < psz[i]; k++)
					if (SR_rid[t_LPI[k]] < cend)
						t_LPI[j++] = t_LPI[k];
				psz[i] = j;
			}
		}
		assert(cend >= Ra);
		return UB + cend - Ra;
	}

	pair<ui, ui> partition_left_right(ui S_end, ui R_end)
	{
		for (ui i = 0; i < S_end; i++)
		{
			ui u = SR[i];
			psz[u] = 0;
			if (support(S_end, u) == 0)
				continue;
			// skipping it, because this is a boundary vertex, and it can't have any non-neighbor candidate
			ui *t_LPI = LPI + u * n;
			char *t_matrix = matrix + u * n;
			for (ui j = S_end; j < R_end; j++)
			{
				ui v = SR[j];
				if (!t_matrix[v])
					t_LPI[psz[u]++] = v;
			}
		}
		ui S_l = 0, R_l = S_end;
		ui UB = 0;

		while (true)
		{
			ui maxpi = -1;
			double maxdise = 1;
			for (ui i = S_l; i < S_end; i++)
			{
				ui u = SR[i];
				if (psz[u] == 0 or support(S_end, u) == 0)
					continue;
				// double cost = min(support(S_end, u), psz[u]);
				// double dise = psz[u] / min(support(S_end, u), psz[u]);
				double dise = (double)psz[u] / support(S_end, u);
				if (dise > maxdise)
					maxpi = u, maxdise = dise;
			}
			if (maxdise == 1)
				break;
			else
			{
				// remove pi* from C
				ui *t_LPI = LPI + maxpi * n;
				for (ui i = 0; i < psz[maxpi]; i++)
				{
					// removing from C
					ui v = t_LPI[i];
					swap_pos(SR_rid[v], R_l++);
				}
				// remove maxpi from every pi
				psz[maxpi] = 0;
				swap_pos(SR_rid[maxpi], S_l++);
				for (ui i = S_l; i < S_end; i++)
				{
					ui j = 0;
					ui u = SR[i];
					ui *t_LPI = LPI + u * n;
					for (ui k = 0; k < psz[u]; k++)
						if (SR_rid[t_LPI[k]] >= R_l)
							t_LPI[j++] = t_LPI[k];
					psz[u] = j;
				}
			}
		}
		return {S_l, R_l};
	}
	ui S_bound(ui S_end, ui R_end, ui level)
	{
		for (ui i = 0; i < S_end; i++)
		{
			ui u = SR[i];
			psz[i] = 0;
			if (support(S_end, u) == 0)
				continue;
			// skipping it, because this is a boundary vertex, and it can't have any non-neighbor candidate
			ui *t_LPI = LPI + i * n;
			for (ui j = S_end; j < R_end; j++)
			{
				ui v = SR[j];
				if (!matrix[u * n + v])
					// PI[u].push_back(v);
					t_LPI[psz[i]++] = v;
			}
		}
		ui beta = best_solution_size - S_end;
		ui cend = R_end;

		while (S_end < cend)
		{
			ui maxpi = -1;
			double maxdise = 0;
			for (ui i = 0; i < S_end; i++)
			{
				ui u = SR[i];
				if (psz[i] == 0)
					continue;
				double cost = min(support(S_end, u), psz[i]);
				double dise = psz[i] / cost;
				if (cost <= beta and dise > maxdise)
					maxpi = i, maxdise = dise;
			}
			if (maxpi == -1)
				break;
			else
			{
				// remove pi* from C
				ui *t_LPI = LPI + maxpi * n;
				for (ui i = 0; i < psz[maxpi]; i++)
				{
					// removing from C
					ui v = t_LPI[i];
					swap_pos(SR_rid[v], --cend);
				}
				// beta-=cost(pi*)
				beta -= min(support(S_end, SR[maxpi]), psz[maxpi]);

				// remove maxpi from every pi
				psz[maxpi] = 0;
				for (ui i = 0; i < S_end; i++)
				{
					ui j = 0;
					ui *t_LPI = LPI + i * n;
					for (ui k = 0; k < psz[i]; k++)
						if (SR_rid[t_LPI[k]] < cend)
							t_LPI[j++] = t_LPI[k];
					psz[i] = j;
				}
			}
			if (beta == 0)
				break;
		}
		if (beta > 0)
			cend -= min(beta, cend - S_end);
		return S_end == cend;
	}

	ui S_branching_backup(ui S_end, ui R_end, ui level)
	{
		ui cend = R_end;
		ui beta = best_solution_size - S_end;
		ui ub = 0;
		do
		{
			ub = tryPartition(S_end, cend);
			if (ub <= beta)
			{
				for (ui i : PIMax)
					swap_pos(i, --cend);
				beta -= ub;
			}
			else
				break;
		} while (beta > 0 && cend > S_end && ub != 0);

		if (beta > 0)
			cend -= min(beta, cend - S_end);

		for (ui i = S_end; i < cend; i++)
		{
			// get a vertex with lowest peelOrder
			ui u = SR[i], ind = i;
			for (ui j = i + 1; j < cend; j++)
			{
				ui v = SR[j];
				if (peelOrder[v] < peelOrder[u])
					ind = j, u = v;
			}
			swap_pos(i, ind);
			swap_pos(i, --R_end);
			level_id[u] = level;
			char *t_matrix = matrix + u * n;
			degree[u] = degree_in_S[u] = 0;
			for (ui i = 0; i < R_end; i++)
			{
				ui w = SR[i];
				// if(level_id[w]==level) continue;
				if (t_matrix[w])
					--degree[w];
			}
		}
		return R_end;
	}
	// ui SR_branching_backup(ui S_end, ui R_end, ui level)
	// {
	//     for (ui i = 0; i < S_end; i++)
	//     {
	//         ui u = SR[i];
	//         psz[i] = 0;
	//         if (support(S_end, u) == 0)
	//             continue;
	//         // skipping it, because this is a boundary vertex, and it can't have any non-neighbor candidate
	// 		ui* t_LPI = LPI+i*n;
	//         for (ui j = S_end; j < R_end; j++)
	//         {
	//             ui v = SR[j];
	//             if (!matrix[u * n + v])
	//                 // PI[u].push_back(v);
	//                 t_LPI[psz[i]++] = v;
	//         }
	//     }
	//     ui beta = best_solution_size - S_end;
	//     ui cend = R_end;

	//     while (S_end<cend)
	//     {

	//         ui maxpi = -1;
	//         double maxdise = 0;
	//         for (ui i = 0; i < S_end; i++)
	//         {
	//             ui u = SR[i];
	//             if (psz[i] == 0)
	//                 continue;
	//             double cost = min(support(S_end, u), psz[i]);
	//             double dise = psz[i] / cost;
	//             if (cost <= beta and dise > maxdise)
	//                 maxpi = i, maxdise = dise;
	//         }

	// 		double ubc = tryColor(S_end, cend);
	// 		double coldise=ISc.size()/ubc;
	// 		if((ubc<=beta)&&(maxpi==-1||coldise>maxdise||(coldise==maxdise&&ISc.size()>psz[maxpi]))){j
	// 			beta-=ubc;
	// 			for(ui i: ISc)
	// 				swap_pos(i, --cend);
	// 			for (ui i = 0; i < S_end; i++)
	//             {
	//                 ui j = 0;
	// 				ui* t_LPI = LPI+i*n;
	//                 for (ui k = 0; k < psz[i]; k++)
	//                     if (SR_rid[t_LPI[k]]<cend) t_LPI[j++] = t_LPI[k];
	//                 psz[i] = j;
	//             }
	// 		}

	//         else
	//         {
	//             // remove pi* from C
	// 			ui* t_LPI = LPI+maxpi*n;
	//             for (ui i = 0; i < psz[maxpi]; i++)
	//             {
	// 				// removing from C
	//                 ui v = t_LPI[i];
	//                 swap_pos(SR_rid[v], --cend);
	//             }
	//             // beta-=cost(pi*)
	//             beta-=min(support(S_end, SR[maxpi]), psz[maxpi]);

	//             // remove maxpi from every pi
	//             psz[maxpi] = 0;
	//             for (ui i = 0; i < S_end; i++)
	//             {
	//                 ui j = 0;
	// 				ui* t_LPI = LPI+i*n;
	//                 for (ui k = 0; k < psz[i]; k++)
	//                     if (SR_rid[t_LPI[k]]<cend) t_LPI[j++] = t_LPI[k];
	//                 psz[i] = j;
	//             }
	//         }

	//         if (beta == 0)
	//             break;
	//     }
	//     if (beta > 0)
	//         cend -= min(beta, cend-S_end);
	//     return move_candidates_to_end(S_end, cend, R_end, level);
	// }

	ui move_candidates_to_end(ui S_end, ui cend, ui R_end, ui level)
	{
		for (ui i = S_end; i < cend; i++)
		{
			// get a vertex with lowest peelOrder at location i
			ui u = SR[i], ind = i;
			for (ui j = i + 1; j < cend; j++)
			{
				ui v = SR[j];
				if (peelOrder[v] < peelOrder[u])
					ind = j, u = v;
			}

			swap_pos(i, ind);
			swap_pos(i, --R_end);

			level_id[u] = level;
			char *t_matrix = matrix + u * n;
			degree[u] = degree_in_S[u] = 0;
			for (ui i = 0; i < R_end; i++)
			{
				ui w = SR[i];
				// if(level_id[w]==level) continue;
				if (t_matrix[w])
					--degree[w];
			}
		}
		return R_end;
	}
	ui getBranchings2(ui S_end, ui R_end, ui level)
	{
		ui cend = bound(S_end, R_end);
		for (ui i = cend; i < R_end; i++)
		{
			// get a vertex with highest peelOrder at location i
			ui u = SR[i], ind = i;
			for (ui j = i + 1; j < R_end; j++)
			{
				ui v = SR[j];
				if (peelOrder[v] > peelOrder[u])
					ind = j, u = v;
			}
			if (i != ind)
				swap_pos(i, ind);
		}
		// remove vertex at i location
		// assert(level_id[u] == level&&SR_rid[u] == R_end);
		while (R_end > cend)
		{
			ui u = SR[--R_end];
			level_id[u] = level;
			char *t_matrix = matrix + u * n;
			degree[u] = degree_in_S[u] = 0;
			for (ui i = 0; i < R_end; i++)
			{
				ui w = SR[i];
				// if(level_id[w]==level) continue;
				if (t_matrix[w])
					--degree[w];
			}
		}

		return cend;
	}

	ui colorUB(ui S_end, ui R_end)
	{
		ui UB = S_end;
		while (R_end > S_end && UB <= best_solution_size)
		{
			double ubc = tryColor(S_end, R_end);
			for (ui v : ISc)
				swap_pos(v, --R_end);
			UB += ubc;
		}
		return UB;
	}


	ui color_ub_partial_C(ui S_end, ui R_end)
	{
		ui UB = 0;
		while (R_end > S_end)
		{
			double ubc = tryColor(S_end, R_end);
			for (ui v : ISc)
				swap_pos(v, --R_end);
			UB += ubc;
		}
		return UB;
	}


	ui seesawUB(ui S_end, ui R_end)
	{
		ui UB = S_end;
		while (R_end > S_end)
		{
			double ubp = tryPartition(S_end, R_end);
			double ubc = tryColor(S_end, R_end);
			if (ubp == 0 or
				(ISc.size() / ubc > PIMax.size() / ubp) or
				((ISc.size() / ubc == PIMax.size() / ubp) and (ISc.size() > PIMax.size())))

			{
				color++;
				for (ui v : ISc)
					swap_pos(v, --R_end);
				UB += ubc;
			}
			else
			{
				part++;
				for (ui v : PIMax)
					swap_pos(v, --R_end);
				UB += ubp;
			}
		}
		return UB;
	}



	void createIS(ui S_end, ui R_end)
	{
		ISc.clear();
		for (ui i = S_end; i < R_end; i++)
		{
			bool flag = true;
			for (ui j : ISc)
				if (is_neigh(i, j))
				{
					flag = false;
					break;
				}
			if (flag)
				ISc.push_back(i);
		}
	}

	ui TISUB(ui S_end)
	{
		ui maxsup = 0;
		for (ui i = 0; i < ISc.size(); i++)
		{
			for (ui j = i + 1; j < ISc.size(); j++)
			{
				if (support(S_end, SR[ISc[j]]) > support(S_end, SR[ISc[i]]))
					std::swap(ISc[i], ISc[j]);
			}
			// not using <= condition because i is starting from 0...
			if (support(S_end, SR[ISc[i]]) > i)
				maxsup++;
			else
				break;
		}
		return maxsup;
	}
	ui tryColor(ui S_end, ui R_end)
	{
		createIS(S_end, R_end);
		ui ub = TISUB(S_end);
		return ub;
		ui vlc = 0;
		// collect loose vertices i.e. v \in ISc | support(v) > ub
		for (ui i = 0; i < ISc.size(); i++)
		{
			if (support(S_end, SR[ISc[i]]) > ub)
			{
				std::swap(ISc[i], ISc[vlc]);
				vlc++;
			}
		}
		// ISc[0... vlc) we have loose vertices

		// Lookup inIS(&lookup, &ISc, true);
		bmp.setup(ISc, n);
		for (ui i = S_end; vlc < ub and i < R_end; i++)
		{
			if (bmp.test(i)) // this loop running for C\ISc
				continue;
			ui vc = 0;
			for (ui j = vlc; j < ISc.size(); j++) // this loop runs in ISc\LC
			{
				if (is_neigh(i, ISc[j]))
				{
					std::swap(ISc[vlc + vc], ISc[j]);
					vc++;
				}
			}
			if (vlc + vc + 1 <= ub)
			{
				vlc += vc;
				ISc.push_back(i);
				std::swap(ISc.back(), ISc[vlc++]);
				bmp.set(i);
			}
		}

		for (ui i = S_end; i < R_end; i++)
		{
			if (bmp.test(i) or support(S_end, SR[i]) >= ub) // this loop running for C\ISc
				continue;
			ui nv = 0;
			for (ui j : ISc)
			{
				if (is_neigh(i, j))
					nv++;
			}

			if (nv + support(S_end, SR[i]) <= ub)
			{
				ISc.push_back(i);
				bmp.set(i);
			}
		}
		return ub;
	}
	bool is_neigh(ui i, ui j)
	{
		return matrix[SR[i] * n + SR[j]];
	}

	ui tryPartition(ui S_end, ui R_end)
	{
		double maxdise = 0;
		ui ub = 0;
		PIMax.clear();
		for (ui i = 0; i < S_end; i++)
		{
			if (support(S_end, SR[i]) == 0)
				continue;
			PI.clear();
			for (ui j = S_end; j < R_end; j++)
			{
				if (!is_neigh(i, j))
					PI.push_back(j);
			}
			if (PI.empty())
				continue;
			ui cost = min(support(S_end, SR[i]), (ui)PI.size());
			double dise = (double)PI.size() / (double)cost;
			if (dise > maxdise or (dise == maxdise and PI.size() > PIMax.size()))
			{
				maxdise = dise, ub = cost;
				std::swap(PI, PIMax);
			}
		}
		return ub;
	}
	ui bound(ui S_end, ui R_end)
	{
		vp.clear();
		for (ui i = 0; i < S_end; i++)
			vp.push_back(std::make_pair(support(S_end, SR[i]), SR[i]));
		// for(ui i = 0;i < S_end;i ++) vp.push_back(std::make_pair(-(degree_in_S[SR[i]]-neiInP[SR[i]]), SR[i]));
		sort(vp.begin(), vp.end());
		ui UB = S_end, cursor = S_end;
		for (ui i = 0; i < (ui)vp.size(); i++)
		{
			ui u = vp[i].second;
			if (vp[i].first == 0)
				continue; // boundary vertex
			ui nn = 0;
			char *t_matrix = matrix + u * n;
			for (ui j = cursor; j < R_end; j++)
				if (!t_matrix[SR[j]])
				{
					if (j != cursor + nn)
						swap_pos(j, cursor + nn);
					++nn;
				}
			ui t_ub = min(nn, vp[i].first);
			// ui t_ub = nn;
			// if(vp[i].first < t_ub) t_ub = vp[i].first;
			if (UB + t_ub <= best_solution_size)
			{
				UB += t_ub;
				cursor += nn;
			}
			else
			{
				return cursor + (best_solution_size - UB);
			}
		}
		cursor += (best_solution_size - UB);
		if (cursor > R_end)
			cursor = R_end;
		return cursor;
	}
	void print_array(const char *str, const ui *array, ui idx_start, ui idx_end, ui l)
	{
		for (ui i = 0; i < l; i++)
			printf(" ");
		printf("%s:", str);
		for (ui i = idx_start; i < idx_end; i++)
			printf(" %u", array[i]);
		printf("\n");
	}
};
#endif
