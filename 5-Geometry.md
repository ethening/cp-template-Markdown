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

## Macau F

```cpp
#include <bits/stdc++.h>
#pragma GCC optimize("unroll-loops, Ofast")
#define sz(a) (int) (a).size()
using namespace std;

typedef double LD;

long double eps = 1e-6;
struct P{
    P(LD x, LD y) : x(x), y(y){};
    P(int x, int y) : x(x), y(y){};
    bool operator== (P p) const { return fabs(x - p.x) < eps && fabs(y - p.y) < eps; }
    P operator- (P p) const { return P {x - p.x, y - p.y}; }
    LD cross(P p) const { return x * p.y - y * p.x; }
    LD cross(P a, P b) const { return (a - *this).cross(b - *this); }
    LD dot(P p) const { return x * p.x + y * p.y; }
    LD x, y;
};

int sgn(LD x){
    if(fabs(x) < eps) return 0;
    return x > 0 ? 1 : -1;
}

string substring(const string& s, int l, int r){
    return s.substr(l, r - l + 1);
}

struct HP{
    HP(int a, int b, int c) : a(a), b(b), c(c){};
    HP(string s) {
        s = substring(s, 1, sz(s) - 2);
        vector<int> sep;
        sep.push_back(-1);
        for(int i = 0; i < sz(s); i++) if(s[i] ==',') sep.push_back(i);
        sep.push_back(sz(s));
        vector<int> val;
        for(int i = 0; i < 3; i++){
            string t = substring(s, sep[i] + 1, sep[i + 1] - 1);
            val.push_back(stoi(t));
        }
        HP::a = val[0], HP::b = val[1], HP::c = val[2];
    }
    P at(LD t) const {
        if(sgn(a) && (!sgn(b) || fabs(a) < fabs(b)))
            return {-c / a + t * b, -t * a};
        else
            return {-t * b, -c / b + t * a};
    }
    LD a, b, c;
    bool in(const P& p){
        return a * p.x + b * p.y + c >= 0;
    }
};

typedef P V;

typedef vector<P> S;

const int LEAF = 0;
const int AND = 1;
const int OR = 2;
const int XOR = 3;
const int NOT = 4;
const int INF = 1 << 30;


string STR;

struct node{
    int type = -1;
    HP plane {INF, INF, INF};
    vector<node*> ch;

    bool in(const P& x){
        if(type == LEAF){
            return plane.in(x);

        }else if(type == AND){
            return ch[0]->in(x) && ch[1]->in(x);

        }else if(type == OR){
            return ch[0]->in(x) || ch[1]->in(x);

        }else if(type == XOR){
            return ch[0]->in(x) ^ ch[1]->in(x);

        }else if(type == NOT){
            return !ch[0]->in(x);

        }
    }
    node(int l, int r){
        // cout << l << ' ' << r << endl;
        if(STR[l] == '['){
            // atomic
            type = LEAF;
            plane = HP(substring(STR, l, r - 1));
            return;
        }else{
            int net = 0;
            // cout << "HELLO" << ' ' << substring(STR, l, r - 1) << endl;
            for(int i = l + 1; i < r; i++){
                if(STR[i] == '(') ++net;
                else if(STR[i] == ')') --net;

                else if(net == 0){
                    if (STR[i] == '&') {
                        type = AND;
                        ch = vector<node*> {
                            new node(l + 1, i),
                            new node(i + 1, r - 1)
                        };
                    } else if (STR[i] =='|') {
                        type = OR;
                        ch = vector<node*> {
                            new node(l + 1, i),
                            new node(i + 1, r - 1)
                        };
                    } else if (STR[i] =='^') {
                        type = XOR;
                        ch = vector<node*> {
                            new node(l + 1, i),
                            new node(i + 1, r - 1)
                        };
                    } else if (STR[i] =='!') {
                        type = NOT;
                        ch = vector<node*> {
                            new node(i + 1, r - 1)
                        };
                    }
                }
            }
        }
    }
};


P operator+ (const P& a, const P& b){
    return P(a.x + b.x, a.y + b.y);
}

P operator* (const P& a, const LD& c){
    return P(a.x * c, a.y * c);
}

P operator/ (const P& a, const LD& c){
    return P(a.x / c, a.y / c);
}

bool onSegment(P s, P e, P p) {
    return sgn(p.cross(s, e)) == 0 && sgn((s - p).dot(e - p)) <= 0;
}

LD cross(P a, P b){ return a.x * b.y - a.y * b.x; }
P center(const S& v) {
    P res(0, 0); double A = 0;
    for (int i = 0, j = sz(v) - 1; i < sz(v); j = i++) {
        res = res + (v[i] + v[j]) * cross(v[j], v[i]);
        A += cross(v[j], v[i]);
    }
    return res / A / 3;
}

vector<S> PS;

P lineInter(const P& s1, const P& e1, const P& s2, const P& e2) {
    auto d = (e1 - s1).cross(e2 - s2);
    if (fabs(d) < eps)
        return {INF, INF};
    auto p = s2.cross(e1, e2), q = s2.cross(e2, s1);
    return (s1 * p + e1 * q) / d;
}

P segInter(const P& a, const P& b, const P& c, const P& d) {
    auto oa = c.cross(d, a), ob = c.cross(d, b),
         oc = a.cross(b, c), od = a.cross(b, d);
    if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0)
        return (a * ob - b * oa) / (ob - oa);
    return {INF, INF};
}

void print(S v){
    for(auto u : v){
        printf("(%lf, %lf) ", u.x, u.y);
    }
    printf("\n");
    fflush(stdout);
}

vector<S> split(const S& v, const HP& p){
    // assert(v.size() >= 3);
    P s2 = p.at(-2e6), e2 = p.at(2e6);

    vector<pair<int, P>> inter;
    for(int i = 0; i < sz(v); i++){
        P s1 = v[i], e1 = v[(i + 1) % sz(v)];
        P a = segInter(s1, e1, s2, e2);
        if(fabs(a.x) <= 2e3 && fabs(a.y) <= 2e3){
            inter.push_back({i, a});
        }
        // cout << a.x << ' ' << a.y << endl;
    }

    // assert(inter.size() <= 2);

    if(inter.empty()){
        for(int i = 0; i < sz(v); i++){
            P s1 = v[i], e1 = v[(i + 1) % sz(v)];
            P a = lineInter(s1, e1, s2, e2);
            if(a == e1){
                inter.push_back({i, a});
            }
        }
        if(inter.size() == 2 && !(inter[0].first + 1 == inter[1].first || inter[0].first == 0 && inter[1].first == sz(v) - 1)){
            vector<S> res;
            {
                S s {inter[0].second};
                for(int id = inter[1].first; id != inter[0].first; id = (id + 1) % sz(v)){
                    int vid = (id + 1) % sz(v);
                    s.push_back(v[vid]);
                }
                res.push_back(s);
                // assert(s.size() >= 3);
            }
            {
                S s {inter[1].second};
                for(int id = inter[0].first; id != inter[1].first; id = (id + 1) % sz(v)){
                    int vid = (id + 1) % sz(v);
                    s.push_back(v[vid]);
                }
                res.push_back(s);
                // assert(s.size() >= 3);
            }
            return res;
        }else{
            return {v};
        }
    }else if(inter.size() == 1){
        for(int i = 0; i < sz(v); i++){
            P s1 = v[i], e1 = v[(i + 1) % sz(v)];
            P a = lineInter(s1, e1, s2, e2);
            if(a == e1){
                inter.push_back({i, a});
                break;
            }
        }
        // assert(inter.size() == 2);
        if(inter.size() == 2){
            vector<S> res;
            {
                S s {inter[0].second};
                for(int id = inter[1].first; id != inter[0].first; id = (id + 1) % sz(v)){
                    int vid = (id + 1) % sz(v);
                    s.push_back(v[vid]);
                }
                res.push_back(s);
                // assert(s.size() >= 3); OK
            }
            {
                S s {inter[0].second};
                for(int id = inter[0].first; id != (inter[1].first + 1) % sz(v); id = (id + 1) % sz(v)){
                    int vid = (id + 1) % sz(v);
                    s.push_back(v[vid]);
                }
                res.push_back(s);
            }
            return res;
        }else{
            return {v};
        }

    }else{
        vector<S> res;
        {
            S s {inter[0].second, inter[1].second};
            for(int id = inter[1].first; id != inter[0].first; id = (id + 1) % sz(v)){
                int vid = (id + 1) % sz(v);
                s.push_back(v[vid]);
            }
            res.push_back(s);
        }
        {
            S s {inter[1].second, inter[0].second};
            for(int id = inter[0].first; id != inter[1].first; id = (id + 1) % sz(v)){
                int vid = (id + 1) % sz(v);
                s.push_back(v[vid]);
            }
            res.push_back(s);
        }
        return res;
    }
}

void dfs0(node* v){
    for(auto u : v->ch){
        dfs0(u);
    }
    if(v->type == LEAF){
        auto& plane = v->plane;
        vector<S> NPS;
        for(S& v : PS){
            for(const S& u : split(v, plane)){
                // assert(u.size() >= 3);
                if(u.size() < 3) continue;
                NPS.push_back(u);
                // print(u);
            }
        }
        PS = NPS;
    }
}

LD polygonArea(const S& v){
    LD a = v.back().cross(v[0]);
    for(int i = 0; i < sz(v) - 1; i++)
        a += v[i].cross(v[i + 1]);
    return a / 2;
}

void solve(){
    int x0, x1, y0, y1; cin >> x0 >> x1 >> y0 >> y1;
    // cout << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1 << endl;
    {
        S init = S{P{x0, y0}, P{x1, y0}, P{x1, y1}, P{x0, y1}};
        // print(init);
        PS.push_back(init);
    }
    cin >> STR;

    node root(0, sz(STR));
    dfs0(&root);

    LD ans = 0;
    for(auto v : PS){
        LD area = polygonArea(v);
        // print(v); cout << area << endl;
        P g = center(v);
        bool black = root.in(g);
        // print(v); cout << area << endl;
        // cout << (black ? "BLACK" : "WHITE") << endl;
        if(black) ans += area;
    }
    cout << ans << '\n';
}

int32_t main(){
#ifndef LOCAL
    cin.tie(0)->sync_with_stdio(false);
#endif
    cout << fixed << setprecision(12);
    int t = 1; // cin >> t;
    while(t--){
        solve();
    }
}
```
