# Geometry

## Planar graph to dual

```cpp
#include "bits/stdc++.h"
#include <algorithm>
#include <ostream>
#include <random>
using namespace std;

using ll = long long;
using LL = long long;
using pii = pair<int, int>;

using ull = unsigned long long;
using int128 = __int128;

ostream& operator<<(ostream& os, int128 p) {
	if (p == 0) return os << 0;
	string s;
	bool flag = 0;
	if (p < 0) { flag = 1; p *= -1; }
	while (p) { s += '0' + p % 10; p /= 10; }
	if (flag) s += '-';
	reverse(s.begin(), s.end());
	return os << s;
}

template<class T> int sgn(T x) { return (x > 0) - (x < 0); }
template<class T>
struct Point {
	typedef Point P;
	T x, y;
	Point(T x = 0, T y = 0) : x(x), y(y) {}
	bool operator<(P p) { return tie(x, y) < tie(p.x, p.y); }
	bool operator==(P p) { return tie(x, y) == tie(p.x, p.y); }

	P operator+(P p) { return P(x+p.x, y+p.y); }
	P operator-(P p) { return P(x-p.x, y-p.y); }
	P operator*(T d) { return P(x*d, y*d); }
	P operator/(T d) { return P(x/d, y/d); }

	T dot(P p) { return x*p.x + y*p.y; }
	T cross(P p) { return x*p.y - y*p.x; }
	T cross(P a, P b) { return (a-*this).cross(b-*this); }
	T dist2() { return x*x + y*y; }
	double dist() { return sqrt(double(dist2())); }
	double angle() { return atan2(ll(y), ll(x)); }

	P unit() { return *this / dist(); }
	P perp() { return P(-y, x); }
	P normal() { return perp().unit(); }

	int phase() {
		if (y != 0) return y > 0 ? 0 : 1;
		return x > 0 ? 0 : 1;
	}

	friend ostream& operator<<(ostream& os, P p) {
		return os << "(" << p.x << "," << p.y << ")";
	}
};

using P = Point<int128>;
double eps = 1e-9;

bool onSegment(P s, P e, P p) {
	return p.cross(s, e) == 0 && (s - p).dot(e - p) <= 0;
}

struct Edge {
	int id;
	int fr, to;

	bool operator==(const Edge &o) {
		return id == o.id;
	}
};

void solve(int TC) {
	int n, m, e;
	cin >> n >> m >> e;
	vector<P> baseP(n), sourceP(m);

	for (int i = 0; i < n; i++) {
		ll x, y; cin >> x >> y;
		x *= 2, y *= 2;
		baseP[i] = {x, y};
	}
	for (int i = 0; i < m; i++) {
		ll x, y; cin >> x >> y;
		x *= 2, y *= 2;
		sourceP[i] = {x, y};
	}

	vector<vector<Edge>> g(n);
	vector<Edge> edge(2 * e);
	for (int i = 0; i < e; i++) {
		int u, v; cin >> u >> v;
		--u, --v;
		edge[2 * i] = {2 * i, u, v};
		g[u].push_back(edge[2 * i]);
		edge[2 * i + 1] = {2 * i + 1, v, u};
		g[v].push_back(edge[2 * i + 1]);
	}

	vector<int> pos(2 * e);
	for (int i = 0; i < n; i++) {
		sort(g[i].begin(), g[i].end(), [&](Edge& a, Edge& b) {
			P u = baseP[a.to] - baseP[a.fr], v = baseP[b.to] - baseP[b.fr];
			if (u.phase() != v.phase()) return u.phase() < v.phase();
			return u.cross(v) > 0;
		});
		for (int j = 0; j < g[i].size(); j++) {
			pos[g[i][j].id] = j;
		}
	}

	int face = 0;
	vector<vector<int>> facePoly;
	vector<int128> faceArea;
	vector<int> belong(2 * e, -1);
	for (int i = 0; i < 2 * e; i++) {
		if (belong[i] != -1) continue;
		vector<int> poly;
		int j = i;
		int128 area = 0;
		while (belong[j] == -1) {
			poly.push_back(j);
			belong[j] = face;
			area += baseP[edge[j].fr].cross(baseP[edge[j].to]);

			int nxt = edge[j].to;
			j = g[nxt][(pos[j ^ 1] - 1 + g[nxt].size()) % g[nxt].size()].id;
		}
		facePoly.push_back(poly);
		faceArea.push_back(area);

		face++;
	}
	vector<P> intervalP(e);
	for (int i = 0; i < e; i++) {
		auto [id, u, v] = edge[2 * i];
		intervalP[i] = (baseP[u] + baseP[v]) / 2;
	}

	auto inPolygon = [&](vector<int> &p, P a, bool strict = true) -> bool {
		int cnt = 0, n = int(size(p));
		for (int i = 0; i < n; i++) {
			P fr = baseP[edge[p[i]].fr];
			P to = baseP[edge[p[i]].to];
			if (onSegment(fr, to, a)) return !strict;
			cnt ^= ((a.y < fr.y) - (a.y < to.y)) * a.cross(fr, to) > 0;
		}
		return cnt;
	};

	mt19937_64 rng(42);
	vector<ull> sourcePHash(m), intervalPHash(2 * e);
	for (int i = 0; i < face; i++) {
		if (faceArea[i] <= 0) continue;

		// cout << "face: " << i << endl;
		ull hash = rng();
		for (int j = 0; j < m; j++) {
			if (inPolygon(facePoly[i], sourceP[j])) {
				// cout << "source point: " << j << " poly: " << i << endl;
				sourcePHash[j] ^= hash;
			}
		}
		for (int j = 0; j < e; j++) {
			bool online = false;
			if (belong[2 * j] == i) {
				online = true;
				// cout << "interval mp: " << 2 * j << " poly: " << i << endl;
				intervalPHash[2 * j] ^= hash;
			}
			if (belong[2 * j + 1] == i) {
				online = true;
				// cout << "interval mp: " << 2 * j + 1 << " poly: " << i << endl;
				intervalPHash[2 * j + 1] ^= hash;
			}
			if (!online && inPolygon(facePoly[i], intervalP[j])) {
				// cout << "interval mp: " << 2 * j << " " << 2 * j + 1 << " poly: " << i << endl;
				intervalPHash[2 * j] ^= hash;
				intervalPHash[2 * j + 1] ^= hash;
			}
		}
	}

	// for (int i = 0; i < e; i++) {
	// 	cout << intervalPHash[2 * i] << " " << intervalPHash[2 * i + 1] << endl;
	// }

	unordered_set<ull> S;
	S.reserve(m * 2);
	for (int i = 0; i < m; i++) S.insert(sourcePHash[i]);

	for (int i = 0; i < e; i++) {
		if (S.count(intervalPHash[2 * i]) || S.count(intervalPHash[2 * i + 1])) {
			cout << "1";
		}
		else {
			cout << "0";
		}
	}
	cout << "\n";
}

int32_t main() {
	cin.tie(0)->sync_with_stdio(0);
	cout << fixed << setprecision(10);

	int t = 1;
	// cin >> t;

	for (int i = 1; i <= t; i++) {
		solve(i);
	}
}
```


