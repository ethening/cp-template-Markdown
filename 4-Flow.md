# Flow

## MaxFlow

```cpp
constexpr int inf = 1E9;

template<class T>
struct MaxFlow {
	struct _Edge {
		int to;
		T cap;
		_Edge(int to, T cap) : to(to), cap(cap) {}
	};

	int n;
	std::vector<_Edge> e;
	std::vector<std::vector<int>> g;
	std::vector<int> cur, h;

	MaxFlow() {}
	MaxFlow(int n) {
		init(n);
	}

	void init(int n) {
		this->n = n;
		e.clear();
		g.assign(n, {});
		cur.resize(n);
		h.resize(n);
	}

	bool bfs(int s, int t) {
		h.assign(n, -1);
		std::queue<int> que;
		h[s] = 0;
		que.push(s);
		while (!que.empty()) {
			const int u = que.front();
			que.pop();
			for (int i : g[u]) {
				auto [v, c] = e[i];
				if (c > 0 && h[v] == -1) {
					h[v] = h[u] + 1;
					if (v == t) {
						return true;
					}
					que.push(v);
				}
			}
		}
		return false;
	}

	T dfs(int u, int t, T f) {
		if (u == t) {
			return f;
		}
		auto r = f;
		for (int &i = cur[u]; i < int(g[u].size()); ++i) {
			const int j = g[u][i];
			auto [v, c] = e[j];
			if (c > 0 && h[v] == h[u] + 1) {
				auto a = dfs(v, t, std::min(r, c));
				e[j].cap -= a;
				e[j ^ 1].cap += a;
				r -= a;
				if (r == 0) {
					return f;
				}
			}
		}
		return f - r;
	}
	void addEdge(int u, int v, T c) {
		g[u].push_back(e.size());
		e.emplace_back(v, c);
		g[v].push_back(e.size());
		e.emplace_back(u, 0);
	}
	T flow(int s, int t) {
		T ans = 0;
		while (bfs(s, t)) {
			cur.assign(n, 0);
			ans += dfs(s, t, std::numeric_limits<T>::max());
		}
		return ans;
	}

	std::vector<bool> minCut() {
		std::vector<bool> c(n);
		for (int i = 0; i < n; i++) {
			c[i] = (h[i] != -1);
		}
		return c;
	}

	struct Edge {
		int from;
		int to;
		T cap;
		T flow;
	};
	std::vector<Edge> edges() {
		std::vector<Edge> a;
		for (int i = 0; i < e.size(); i += 2) {
			Edge x;
			x.from = e[i + 1].to;
			x.to = e[i].to;
			x.cap = e[i].cap + e[i + 1].cap;
			x.flow = e[i + 1].cap;
			a.push_back(x);
		}
		return a;
	}
};
```

## Min Cost Flow

```c++
constexpr int inf = 1E9;

template<class T>
struct MinCostFlow {
	struct _Edge {
		int to;
		T cap;
		T cost;
		_Edge(int to_, T cap_, T cost_) : to(to_), cap(cap_), cost(cost_) {}
	};
	int n;
	std::vector<_Edge> e;
	std::vector<std::vector<int>> g;
	std::vector<T> h, dis;
	std::vector<int> pre;
	bool dijkstra(int s, int t) {
		dis.assign(n, std::numeric_limits<T>::max());
		pre.assign(n, -1);
		std::priority_queue<std::pair<T, int>, std::vector<std::pair<T, int>>, std::greater<std::pair<T, int>>> que;
		dis[s] = 0;
		que.emplace(0, s);
		while (!que.empty()) {
			T d = que.top().first;
			int u = que.top().second;
			que.pop();
			if (dis[u] != d) {
				continue;
			}
			for (int i : g[u]) {
				int v = e[i].to;
				T cap = e[i].cap;
				T cost = e[i].cost;
				if (cap > 0 && dis[v] > d + h[u] - h[v] + cost) {
					dis[v] = d + h[u] - h[v] + cost;
					pre[v] = i;
					que.emplace(dis[v], v);
				}
			}
		}
		return dis[t] != std::numeric_limits<T>::max();
	}
	MinCostFlow() {}
	MinCostFlow(int n_) {
		init(n_);
	}
	void init(int n_) {
		n = n_;
		e.clear();
		g.assign(n, {});
	}
	void addEdge(int u, int v, T cap, T cost) {
		g[u].push_back(e.size());
		e.emplace_back(v, cap, cost);
		g[v].push_back(e.size());
		e.emplace_back(u, 0, -cost);
	}
	std::pair<T, T> flow(int s, int t) {
		T flow = 0;
		T cost = 0;
		h.assign(n, 0);
		while (dijkstra(s, t)) {
			for (int i = 0; i < n; ++i) {
				h[i] += dis[i];
			}
			T aug = std::numeric_limits<int>::max();
			for (int i = t; i != s; i = e[pre[i] ^ 1].to) {
				aug = std::min(aug, e[pre[i]].cap);
			}
			for (int i = t; i != s; i = e[pre[i] ^ 1].to) {
				e[pre[i]].cap -= aug;
				e[pre[i] ^ 1].cap += aug;
			}
			flow += aug;
			cost += aug * h[t];
		}
		return std::make_pair(flow, cost);
	}
	struct Edge {
		int from;
		int to;
		T cap;
		T cost;
		T flow;
	};
	std::vector<Edge> edges() {
		std::vector<Edge> a;
		for (int i = 0; i < e.size(); i += 2) {
			Edge x;
			x.from = e[i + 1].to;
			x.to = e[i].to;
			x.cap = e[i].cap + e[i + 1].cap;
			x.cost = e[i].cost;
			x.flow = e[i + 1].cap;
			a.push_back(x);
		}
		return a;
	}
};
```
