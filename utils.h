#ifndef UTILS
#define UTILS
#include "graph.h"
enum RecLevel{FIRST, OTHER};

class ListLinearHeap
{
private:
	ui n;		// number vertices
	ui key_cap; // the maximum allowed key value

	ui min_key; // possible min key
	ui max_key; // possible max key

	ui *key_s; // key of vertices

	ui *head_s; // head of doubly-linked list for a specific weight

	ui *pre_s;	// pre for doubly-linked list
	ui *next_s; // next for doubly-linked list

public:
	ListLinearHeap(ui _n, ui _key_cap)
	{
		n = _n;
		key_cap = _key_cap;

		min_key = max_key = key_cap;

		head_s = key_s = pre_s = next_s = nullptr;
	}

	ListLinearHeap() : ListLinearHeap(0, 0)
	{
		head_s = key_s = pre_s = next_s = nullptr;
	}
	void init(ui _n)
	{
		n = _n;
		key_cap = n - 1;
		min_key = max_key = key_cap;
		if (key_s == nullptr)
			key_s = new ui[n];
		if (pre_s == nullptr)
			pre_s = new ui[n];
		if (next_s == nullptr)
			next_s = new ui[n];
		if (head_s == nullptr)
			head_s = new ui[key_cap + 1];
		for (ui i = 0; i <= key_cap; i++)
			head_s[i] = n;
	}
	void init(Graph &g)
	{
		n = g.V;
		assert(n);
		key_cap = n - 1;
		min_key = max_key = key_cap;
		if (key_s == nullptr)
			key_s = new ui[n];
		if (pre_s == nullptr)
			pre_s = new ui[n];
		if (next_s == nullptr)
			next_s = new ui[n];
		if (head_s == nullptr)
			head_s = new ui[key_cap + 1];
		for (ui i = 0; i <= key_cap; i++)
			head_s[i] = n;

		for (ui i = 0; i < n; i++)
		{
			ui id = i;
			ui key = g.adjList[i].size();
			assert(id < n);
			assert(key <= key_cap);

			key_s[id] = key;
			pre_s[id] = n;
			next_s[id] = head_s[key];
			if (head_s[key] != n)
				pre_s[head_s[key]] = id;
			head_s[key] = id;

			if (key < min_key)
				min_key = key;
		}
	}
	~ListLinearHeap()
	{
		if (head_s != nullptr)
		{
			delete[] head_s;
			head_s = nullptr;
		}
		if (pre_s != nullptr)
		{
			delete[] pre_s;
			pre_s = nullptr;
		}
		if (next_s != nullptr)
		{
			delete[] next_s;
			next_s = nullptr;
		}
		if (key_s != nullptr)
		{
			delete[] key_s;
			key_s = nullptr;
		}
	}

	void init(ui _n, ui _key_cap, ui *_id_s, ui *_key_s)
	{
		if (key_s == nullptr)
			key_s = new ui[n];
		if (pre_s == nullptr)
			pre_s = new ui[n];
		if (next_s == nullptr)
			next_s = new ui[n];
		if (head_s == nullptr)
			head_s = new ui[key_cap + 1];

		// assert(_key_cap <= key_cap);
		min_key = max_key = _key_cap;
		for (ui i = 0; i <= _key_cap; i++)
			head_s[i] = n;

		for (ui i = 0; i < _n; i++)
		{
			ui id = _id_s[i];
			ui key = _key_s[id];
			// assert(id < n); assert(key <= _key_cap);

			key_s[id] = key;
			pre_s[id] = n;
			next_s[id] = head_s[key];
			if (head_s[key] != n)
				pre_s[head_s[key]] = id;
			head_s[key] = id;

			if (key < min_key)
				min_key = key;
		}
	}

	ui get_key(ui id) { return key_s[id]; }

	void get_ids(ui *vs, ui &vs_size)
	{
		for (ui i = min_key; i <= max_key; i++)
		{
			for (ui id = head_s[i]; id != n; id = next_s[id])
			{
				vs[vs_size++] = id;
			}
		}
	}

	bool get_min(ui &id, ui &key)
	{ // return true if success, return false otherwise
		while (min_key <= max_key && head_s[min_key] == n)
			++min_key;
		if (min_key > max_key)
			return false;

		id = head_s[min_key];
		key = min_key;

		// assert(key_s[id] == key);

		return true;
	}

	bool pop_min(ui &id, ui &key)
	{ // return true if success, return false otherwise
		while (min_key <= max_key && head_s[min_key] == n)
			++min_key;
		if (min_key > max_key)
			return false;

		id = head_s[min_key];
		key = min_key;

		key_s[id] = key_cap + 1;
		// assert(key_s[id] == key);

		head_s[min_key] = next_s[id];
		if (head_s[min_key] != n)
			pre_s[head_s[min_key]] = n;
		return true;
	}

	ui decrement(ui id, ui dec)
	{

		assert(key_s[id] >= dec);
		if (key_s[id] > key_cap)
			return 0;

		if (pre_s[id] == n)
		{
			// assert(head_s[key_s[id]] == id);
			head_s[key_s[id]] = next_s[id];
			if (next_s[id] != n)
				pre_s[next_s[id]] = n;
		}
		else
		{
			ui pid = pre_s[id];
			next_s[pid] = next_s[id];
			if (next_s[id] != n)
				pre_s[next_s[id]] = pid;
		}

		ui &key = key_s[id];
		key -= dec;
		pre_s[id] = n;
		next_s[id] = head_s[key];
		if (head_s[key] != n)
			pre_s[head_s[key]] = id;
		head_s[key] = id;

		if (key < min_key)
			min_key = key;
		return key;
	}
};

class Timer
{
#define TIME_NOW chrono::steady_clock::now()

public:
	Timer() : m_start(TIME_NOW) {}
	void restart() { m_start = TIME_NOW; }
	long long elapsed()
	{
		return chrono::duration_cast<chrono::microseconds>(TIME_NOW - m_start).count();
	}
	void tick(){
		tic = TIME_NOW;
	}
	void tock(){
		toc+=chrono::duration_cast<chrono::nanoseconds>(TIME_NOW - tic).count();
	}
	double ticktock(){
		return toc/(1000'000'000.0);
	}
private:
	std::chrono::steady_clock::time_point m_start, tic;
	unsigned long long toc=0;
};

class CommandLine
{
public:
	int argc;
	char **argv;

	CommandLine(int _argc, char **_argv) : argc(_argc), argv(_argv) {}

	void BadArgument()
	{
		std::cout << "usage: " << argv[0] << " bad argument" << std::endl;
		abort();
	}

	char *GetOptionValue(const std::string &option)
	{
		for (int i = 1; i < argc; i++)
			if ((std::string)argv[i] == option)
				return argv[i + 1];
		return NULL;
	}

	std::string GetOptionValue(const std::string &option, std::string defaultValue)
	{
		for (int i = 1; i < argc; i++)
			if ((std::string)argv[i] == option)
				return (std::string)argv[i + 1];
		return defaultValue;
	}

	int GetOptionIntValue(const std::string &option, int defaultValue)
	{
		for (int i = 1; i < argc; i++)
			if ((std::string)argv[i] == option)
			{
				int r = atoi(argv[i + 1]);
				return r;
			}
		return defaultValue;
	}

	long GetOptionLongValue(const std::string &option, long defaultValue)
	{
		for (int i = 1; i < argc; i++)
			if ((std::string)argv[i] == option)
			{
				long r = atol(argv[i + 1]);
				return r;
			}
		return defaultValue;
	}

	double GetOptionDoubleValue(const std::string &option, double defaultValue)
	{
		for (int i = 1; i < argc; i++)
			if ((std::string)argv[i] == option)
			{
				double val;
				if (sscanf(argv[i + 1], "%lf", &val) == EOF)
				{
					BadArgument();
				}
				return val;
			}
		return defaultValue;
	}
};
Timer check2;
class MBitSet
{
private:
public:
	ui n;
	ui cap;
	ui *buf;

	MBitSet()
	{
		buf = nullptr;
		cap = n = 0;
	}
	MBitSet(ui _cap)
	{
		cap = _cap;
		n = (cap >> 5) + 1;
		buf = new ui[n];
		fill(buf, buf + n, 0);
		// for (ui i = 0; i < n; ++i)
		// 	buf[i] = 0;
	}
	void init(ui _cap)
	{
		cap = _cap;
		n = (cap >> 5) + 1;
		buf = new ui[n];
		fill(buf, buf + n, 0);
	}
	~MBitSet()
	{
		// todo doing double free, see what causing it and fix

		if (buf != nullptr)
			// delete[] buf;
			buf = nullptr;
	}
	void reset()
	{
		fill(buf, buf + n, 0);
	}

	void setAll()
	{
		fill(buf, buf + n, 0xffffffff);
	}
	// FLIP all the bits
	void flip()
	{
		for (ui i = 0; i < n; ++i)
			buf[i] = ~buf[i];
	}
	void set(ui x)
	{
		// assert(x < cap);
		buf[x >> 5] |= (ui)1 << (x & 31);
	}

	bool test(ui x)
	{
		// cout << x << " " << n << " " << cap << endl;
		return buf[x >> 5] >> (x & 31) & 1;
	}

	bool empty()
	{
		for (ui i = 0; i < n; ++i)
			if (buf[i])
				return false;
		return true;
	}
	void setup(vecui &adj, ui m)
	{
		m = (m >> 5) + 1;
		fill(buf, buf + m, 0);
		
		for (ui u : adj)
			set(u);
		
	}
};

class Lookup
{
	vector<ui> *lookup;
	vector<ui> *data;
	bool mode;

public:
	Lookup(vector<ui> *_lookup, vector<ui> *_data, bool _mode = false) : lookup(_lookup), data(_data)
	{
		mode = _mode;
		for (ui ind = 0; ind < data->size(); ind++)
		{
			ui u = data->at(ind);
			if (mode)
				lookup->at(u) = 1;
			else
				lookup->at(u) = ind + 1;
		}
	}

	~Lookup()
	{
		erase();
	}

	void erase()
	{
		for (const ui &u : *data)
		{
			lookup->at(u) = 0;
		}
	}

	ui &operator[](ui ind)
	{
		return lookup->at(ind);
	}
};

class RandList
{
private:
	ui *vlist;
	ui *vpos;
	ui vnum;
	ui cap;

public:
	RandList()
	{
		vlist = vpos = nullptr;
		vnum = cap = 0;
	};
	RandList(int _cap)
	{
		init(_cap);
		// cap = _cap;
		// vlist = new ui[cap];
		// vpos = new ui[cap];
		// vnum = 0;
		// for (ui i = 0; i < cap; i++)
		// {
		// 	vpos[i] = cap;
		// }
	}
	RandList(const RandList &rl)
	{
		cap = rl.cap;
		vlist = new ui[cap];
		vpos = new ui[cap];
		vnum = rl.vnum;
		memcpy(vlist, rl.vlist, sizeof(ui) * cap);
		memcpy(vpos, rl.vpos, sizeof(ui) * cap);

		// for (ui i = 0; i < cap; i++)
		// {
		// 	vpos[i] = rl.vpos[i];
		// 	vlist[i] = rl.vlist[i];
		// }
	}
	void init(int _cap)
	{
		cap = _cap;
		vlist = new ui[cap];
		vpos = new ui[cap];
		vnum = 0;
		for (ui i = 0; i < cap; i++)
		{
			vpos[i] = cap;
		}
	}
	void add(int vid)
	{
		assert(vpos[vid] == cap);
		vlist[vnum] = vid;
		vpos[vid] = vnum;
		vnum++;
	};
	void remove(int vid)
	{
		assert(vpos[vid] < vnum);
		ui last_id = vlist[vnum - 1];
		ui id_pos = vpos[vid];
		vlist[id_pos] = last_id;
		vpos[last_id] = id_pos;
		vnum--;
		vpos[vid] = cap; /*set as visited*/
	}

	void clear()
	{
		for (ui i = 0; i < vnum; i++)
			vpos[vlist[i]] = cap;
		vnum = 0;
	}
	ui get(ui i)
	{
		return vlist[i];
	}

	ui getIndex(ui vid)
	{
		assert(contains(vid));
		return vpos[vid];
	}

	void swapElements(ui i, ui j)
	{
		// swaps element at i index with j index
		ui u = vlist[i];
		ui v = vlist[j];
		swap(vlist[i], vlist[j]);
		swap(vpos[u], vpos[v]);
	}
	ui operator[](ui i)
	{
		assert(i < vnum);
		return vlist[i];
	}
	bool contains(int vid)
	{
		if (vid >= cap)
			return false;
		return vpos[vid] != cap;
	}

	bool empty() { return vnum == 0; }
	ui size() { return vnum; }
	ui getCap() { return cap; }
	vector<ui> getData()
	{
		vector<ui> data;
		data.insert(data.begin(), vlist, vlist + vnum);
		return data;
	}
	void copyDataTo(vector<ui> &data)
	{
		data.insert(data.begin(), vlist, vlist + vnum);
	}
	ui top()
	{
		assert(!empty());
		return vlist[vnum - 1];
	}
	ui pop()
	{
		assert(!empty());
		vnum--;
		ui v = vlist[vnum];
		vpos[v] = cap;
		return v;
	}
	ui fakePop()
	{
		assert(!empty());
		return vlist[--vnum];
	}
	ui fakeRecPop()
	{
		return vlist[vnum++];
	}
	void fakeRemove(ui v)
	{
		assert(!empty());
		vnum--;
		ui idx = vpos[v];
		ui u = vlist[vnum]; // last element
		swap(vpos[u], vpos[v]);
		swap(vlist[idx], vlist[vnum]);
	}
	void fakeRecover(ui sz)
	{
		vnum += sz;
	}
	void loadData(vector<ui> data)
	{
		clear();
		for (ui u : data)
		{
			add(u);
		}
	}
	void dispose()
	{
		if (vlist != nullptr)
		{
			delete[] vlist;
			vlist = nullptr;
		}
		if (vpos != nullptr)
		{
			delete[] vpos;
			vpos = nullptr;
		}
	}
	~RandList()
	{
		dispose();
	}
	void print()
	{
		if (vnum == 0)
			return;
		for (ui i = 0; i < vnum; i++)
			cout << vlist[i] << " ";
		cout << "; " << endl;
	}
#ifdef DBGMOD
	void printList(FILE *f = stdout)
	{
		fprintf(f, "Total %d: ", vnum);
		int *tmp_lst = new int[cap];
		memcpy(tmp_lst, vlist, vnum * sizeof(int));
		std::sort(tmp_lst, tmp_lst + vnum);
		// qsort(tmp_lst, vnum, sizeof(int), cmpfunc);
		for (ui i = 0; i < vnum; i++)
		{
			fprintf(f, "%d ", tmp_lst[i]);
		}
		fprintf(f, "\n");
	};
#else
	void printList(FILE *f = stdout){};
#endif
};

inline ui getLowerBound(auto &vec, ui x)
{
	auto it = lower_bound(vec.begin(), vec.end(), x);

	if (it == vec.end() || *it != x)
		return vec.size();
	return it - vec.begin();
}
inline ui isNeighbor(auto &vec, ui x)
{
	return binary_search(vec.begin(), vec.end(), x);
}
void printvec(string st, vecui vec)
{
	cout << st << " : ";
	for (ui v : vec)
		cout << v << " ";
	cout << endl;
}


#endif