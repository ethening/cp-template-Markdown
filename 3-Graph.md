# Graph

## SCC

```c++
struct SCC {
    int n;
    std::vector<std::vector<int>> adj;
    std::vector<int> stk;
    std::vector<int> dfn, low, bel;
    int cur, cnt;

    SCC() {}
    SCC(int n) {
        init(n);
    }

    void init(int n) {
        this->n = n;
        adj.assign(n, {});
        dfn.assign(n, -1);
        low.resize(n);
        bel.assign(n, -1);
        stk.clear();
        cur = cnt = 0;
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
    }

    void dfs(int x) {
        dfn[x] = low[x] = cur++;
        stk.push_back(x);

        for (auto y : adj[x]) {
            if (dfn[y] == -1) {
                dfs(y);
                low[x] = std::min(low[x], low[y]);
            } else if (bel[y] == -1) {
                low[x] = std::min(low[x], dfn[y]);
            }
        }

        if (dfn[x] == low[x]) {
            int y;
            do {
                y = stk.back();
                bel[y] = cnt;
                stk.pop_back();
            } while (y != x);
            cnt++;
        }
    }

    std::vector<int> work() {
        for (int i = 0; i < n; i++) {
            if (dfn[i] == -1) {
                dfs(i);
            }
        }
        return bel;
    }
};
```

## Centroid Decomposition

```c++
void solve() {
	int n, T;
	cin >> n >> T;

	vector<vector<array<int, 2>>> g(n);

	for (int i = 1; i < n; i++) {
		int u, v, w;
		cin >> u >> v >> w;
		u--, v--;
		g[u].push_back({v, w});
		g[v].push_back({u, w});
	}

	int q;
	cin >> q;

	vector<vector<array<int, 3>>> qry(n);
	for (int i = 0; i < q; i++) {
		int a, b;
		cin >> a >> b;
		a--, b--;
		qry[0].push_back({a, b, i});
	}
	vector<ll> ans(q, inf);

	vector<bool> vis(n);
	vector<int> siz(n), bel(n);

	auto dfs = [&](auto self, int x, int p) -> void {
		siz[x] = 1;
		for (auto [y, w] : g[x]) {
			if (y == p) {
				continue;
			}
			self(self, y, x);
			siz[x] += siz[y];
		}
	};
	dfs(dfs, 0, -1);

	vector<ll> dep(n);
	auto solve = [&](auto &&self, int r) -> void {
		auto Q = std::move(qry[r]);

		auto find = [&](auto self, int x, int p, int s) -> int {
			for (auto [y, _] : g[x]) {
				if (y == p || vis[y] || 2 * siz[y] <= s) {
					continue;
				}
				return self(self, y, x, s);
			}
			return x;
		};
		r = find(find, r, -1, siz[r]);
		vis[r] = true;

		auto dfs = [&](auto self, int x, int p) -> void {
			siz[x] = 1;
			for (auto [y, w] : g[x]) {
				if (y == p || vis[y]) {
					continue;
				}
				dep[y] = dep[x] + w;
				bel[y] = x == r ? y : bel[x];
				self(self, y, x);
				siz[x] += siz[y];
			}
		};
		dfs(dfs, r, -1);

		for (auto [a, b, i] : Q) {
			ans[i] = min(ans[i], dep[a] + dep[b]);
			if (bel[a] == bel[b]) {
				qry[bel[a]].push_back({a, b, i});
			}
		}
		for (auto [y, _] : g[r]) {
			if (!vis[y]) {
				self(self, y);
			}
		}
	};
	solve(solve, 0);

	for (int i = 0; i < q; i++) {
		cout << ans[i] << "\n";
	}
}

```

## DSU On Tree

```c++
struct FreqBuckets {
	vector<int> occ;
	vector<int> freq;
	FreqBuckets(int n, int maxC) : occ(maxC + 1, 0), freq(n + 1) { }

	void add(int x, int mul) {
		if (mul == +1) {
			occ[x]++;
			freq[occ[x]]++;
		}
		else if (mul == -1) {
			freq[occ[x]]--;
			occ[x]--;
		}
		else assert(false);
	}
};

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);

	int n,m; cin >> n >> m;

	vector<int> c(n);
	for (int i = 0; i < n; i++)
		cin >> c[i];

	vector<vector<int>> g(n);
	for (int i = 0; i + 1 < n; i++) {
		int u,v; cin >> u >> v; u--; v--;
		g[u].emplace_back(v);
		g[v].emplace_back(u);
	}

	vector<int> sz(n, 1);
	function<void(int, int)> dfs_hld = [&](int u, int p) {
		if (p != -1) {
			auto it = find(g[u].begin(), g[u].end(), p);
			assert(it != g[u].end());
			g[u].erase(it);
		}

		for (auto& v : g[u]) {
			dfs_hld(v, u);
			sz[u] += sz[v];
			if (sz[v] > sz[g[u][0]])
				swap(v, g[u][0]);
		}
	};
	dfs_hld(0, -1);

	vector<vector<pair<int, int>>> qry(n);
	for (int i = 0; i < m; i++) {
		int v,k; cin >> v >> k; v--;
		qry[v].emplace_back(k, i);
	}

	const int maxC = 100000;
	FreqBuckets cnt(n, maxC);
	vector<int> ans(m, -1);

	function<void(int, int)> dfs_addonly = [&](int u, int mul) {
		cnt.add(c[u], mul);
		for (auto& v : g[u])
			dfs_addonly(v, mul);
	};

	function<void(int)> dfs_solve = [&](int u) {
		for (auto& v : g[u]) {
			if (v == g[u][0]) continue;
			dfs_solve(v);
			dfs_addonly(v, -1);
		}

		if (!g[u].empty())
			dfs_solve(g[u][0]);
		cnt.add(c[u], +1);
		for (auto& v : g[u]) {
			if (v == g[u][0]) continue;
			dfs_addonly(v, +1);
		}

		for (auto& [k, i] : qry[u])
			ans[i] = (k <= n ? cnt.freq[k] : 0);
	};
	dfs_solve(0);

	for (int i = 0; i < m; i++)
		cout << ans[i] << "\n";

	return 0;
}
```
