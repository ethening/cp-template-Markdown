# Math


## Matrix Inverse

Given matrix A, find inverse B. You are given matrix C which is matrix B edited at most 12 position.

jiangly's team solution

```cpp
#include <bits/stdc++.h>

using i64 = long long;

constexpr int P = 1000000007;

std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

int power(int a, int b) {
    int res = 1;
    for (; b; b /= 2, a = 1LL * a * a % P) {
        if (b % 2) {
            res = 1LL * res * a % P;
        }
    }
    return res;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n;
    std::cin >> n;

    std::vector A(n, std::vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> A[i][j];
            // A[i][j] = (i == j);
        }
    }
    std::vector C(n, std::vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> C[i][j];
            // C[i][j] = (i == j);
        }
    }
    // for (int i = 0; i < 12; i++) {
    //     C[rng() % n][rng() % n] = rng() % P;
    // }

    std::vector<int> idr(n), idc(n);
    for (int t = 0; t < 10; t++) {
        std::vector<int> v(n);
        for (int i = 0; i < n; i++) {
            v[i] = rng() % P;
        }
        std::vector<int> vA(n), vAC(n), Av(n), CAv(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                vA[j] = (vA[j] + 1LL * v[i] * A[i][j]) % P;
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                vAC[j] = (vAC[j] + 1LL * vA[i] * C[i][j]) % P;
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Av[i] = (Av[i] + 1LL * A[i][j] * v[j]) % P;
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                CAv[i] = (CAv[i] + 1LL * C[i][j] * Av[j]) % P;
            }
        }
        for (int i = 0; i < n; i++) {
            if (vAC[i] != v[i]) {
                idc[i] = 1;
            }
            if (CAv[i] != v[i]) {
                idr[i] = 1;
            }
        }
    }

    // for (int i = 0; i < n; i++) {
    //     std::cout << idr[i];
    // }
    // std::cout << "\n";
    // for (int i = 0; i < n; i++) {
    //     std::cout << idc[i];
    // }
    // std::cout << "\n";

    std::vector<int> row, col;
    int nr = 0, nc = 0;
    for (int i = 0; i < n; i++) {
        if (idr[i]) {
            idr[i] = nr++;
            row.push_back(i);
        } else {
            idr[i] = -1;
        }
        if (idc[i]) {
            idc[i] = nc++;
            col.push_back(i);
        } else {
            idc[i] = -1;
        }
    }

    int tot = nr * nc;
    int rank = 0;
    std::vector f(tot, std::vector<int>(tot + 1));
    while (tot > rank) {
        // std::cerr << "rank : " << rank << "\n";
        std::vector<int> v(n);
        for (int i = 0; i < n; i++) {
            v[i] = rng() % P;
        }
        std::vector<int> vA(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                vA[j] = (vA[j] + 1LL * v[i] * A[i][j]) % P;
            }
        }
        for (int j = 0; j < n; j++) {
            if (idc[j] == -1) {
                continue;
            }
            int res = 0;
            std::vector<int> g(tot + 1);
            for (int i = 0; i < n; i++) {
                if (idr[i] == -1) {
                    res = (res + 1LL * vA[i] * C[i][j]) % P;
                } else {
                    g[idr[i] * nc + idc[j]] = vA[i];
                }
            }
            g[tot] = (v[j] - res + P) % P;
            // for (int i = 0; i <= tot; i++) {
            //     std::cerr << g[i] << " \n"[i == tot];
            // }
            for (int i = 0; i < tot; i++) {
                if (g[i] != 0) {
                    if (f[i][i] == 0) {
                        int v = power(g[i], P - 2);
                        for (int j = 0; j <= tot; j++) {
                            g[j] = 1LL * g[j] * v % P;
                        }
                        f[i] = g;
                        for (int j = 0; j < i; j++) {
                            int x = f[j][i];
                            for (int k = j; k <= tot; k++) {
                                f[j][k] = (f[j][k] + 1LL * (P - x) * g[k]) % P;
                            }
                        }
                        rank++;
                        break;
                    }
                    int x = g[i];
                    for (int j = i; j <= tot; j++) {
                        g[j] = (g[j] + 1LL * (P - x) * f[i][j]) % P;
                    }
                }
            }
        }
    }

    std::vector<std::array<int, 3>> ans;
    for (int i = 0; i < tot; i++) {
        if (f[i][tot] != C[row[i / nc]][col[i % nc]]) {
            ans.push_back({row[i / nc] + 1, col[i % nc] + 1, f[i][tot]});
        }
    }
    std::cout << ans.size() << "\n";
    for (auto [x, y, z] : ans) {
        std::cout << x << " " << y << " " << z << "\n";
    }

    return 0;
}
```

## Energy Distribution

There are n planets in the galaxy. Some undirected tunnels connect planets. There exists at most one tunnel connecting each pair of planets. So these tunnels can be described as an n x n matrix W(n x n). Specifically, the tunnel connecting planet i and j has a width of w(i, j) (If there is no tunnel between planet i and j, then wi, j=0).

Now, you want to distribute exactly 1.0 unit of energy among the n planets. Suppose that planet i is distributed ei (a real number) unit of energy (sum of ei = 1), these planets will bring E magical value, where E = sum ei x ej x wi,j

Please distribute the energy and maximize the magical value.

```cpp
#include "bits/stdc++.h"
#define N 15
using namespace std;

using ll = long long;

#define vi vector<int>
#define sz(x) (int)size((x))
#define all(x) x.begin(), x.end()
typedef vector<double> vd;
const double eps=1e-12;
int solveLinear(vector<vd>& A,vd& b, vd& x){
	for(int i=1;i<=sz(b);i++) x.push_back(0);
	int n=sz(A),m=sz(x),rk=0,br,bc;
	if(n) assert(sz(A[0])==m);
	vi col(m);iota(all(col),0);
	for(int i=0;i<n;i++){
		double v,bv=0;
		for(int r=i;r<n;r++) for(int c=i;c<m;c++)
		if((v=fabs(A[r][c]))>bv) br=r,bc=c,bv=v;
		if(bv<=eps){
			for(int j=i;j<n;j++) if(fabs(b[j])>eps) return -1;
			break;
		}
		swap(A[i],A[br]);
		swap(b[i],b[br]);
		swap(col[i],col[bc]);
		for(int j=0;j<n;j++) swap(A[j][i],A[j][bc]);
		bv=1/A[i][i];for(int j=i+1;j<n;j++){
			double fac=A[j][i]*bv;
			b[j]-=fac*b[i];
			for(int k=i+1;k<m;k++) A[j][k]-=fac*A[i][k];
		}
		rk++;
	}
	x.assign(m,0);
	for(int i=rk;i--;){
		b[i]/=A[i][i];
		x[col[i]]=b[i];
		for(int j=0;j<i;j++) b[j]-=A[j][i]*b[i];
	}
	return rk;
}
int n,w[N][N];double e[N];bool ok[N];double Ans;
void solve(int TC) {
	cin>>n;for(int i=1;i<=n;i++) for(int j=1;j<=n;j++) cin>>w[i][j];
	for(int i=0;i<(1<<n);i++){
		vector<vd> A;vd b,ans;int tot=0;
		for(int j=1;j<=n;j++) if(!(i&(1<<(j-1)))) ok[j]=1,tot++;else ok[j]=0;
		if(tot<=1) continue;int nn=0;
		for(int j=n;j;j--) if(ok[j]){nn=j;break;}
		for(int j=1;j<nn;j++) if(ok[j]){
			b.clear();for(int k=1;k<nn;k++) if(ok[k]){
				if(j!=k) b.push_back(w[j][k]-w[j][nn]-w[k][nn]);
				else b.push_back(-2.0*w[j][nn]);
			}
			A.push_back(b);
		}
		b.clear();
		for(int j=1;j<nn;j++) if(ok[j]) b.push_back(-w[j][nn]);
		// cout<<b.size()<<' '<<A[0].size()<<'\n';
		solveLinear(A, b, ans);bool flg=0;double sum=0;
		// cout<<tot<<' '<<sz(b)<<' '<<sz(ans)<<'\n';
		for(int j=0;j<tot-1;j++) if(ans[j]<0){flg=1;break;}else sum+=ans[j];
		if(sum>1.0) flg=1;if(!flg){
			int tt=0;
			for(int j=1;j<nn;j++) if(ok[j]) e[j]=ans[tt],tt++;
			e[nn]=1.0-sum;for(int j=1;j<=n;j++) if(!ok[j]) e[j]=0;sum=0;
			for(int j=1;j<=n;j++) for(int k=j+1;k<=n;k++) sum+=e[j]*e[k]*w[j][k];
			Ans=max(Ans,sum);
		}
	}
	printf("%.9f\n",Ans);
}

int main() {
	cin.tie(0)->sync_with_stdio(0);
	cout << fixed << setprecision(9);
	int t = 1;
	// cin >> t;
	for (int i = 1; i <= t; i++) {
		solve(i);
	}
}
```

## Nash Equalibrium

```cpp
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <bitset>
#include <complex>
#include <deque>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

using Int = long long;

template <class T1, class T2> ostream &operator<<(ostream &os, const pair<T1, T2> &a) { return os << "(" << a.first << ", " << a.second << ")"; };
template <class T> ostream &operator<<(ostream &os, const vector<T> &as) { const int sz = as.size(); os << "["; for (int i = 0; i < sz; ++i) { if (i >= 256) { os << ", ..."; break; } if (i > 0) { os << ", "; } os << as[i]; } return os << "]"; }
template <class T> void pv(T a, T b) { for (T i = a; i != b; ++i) cerr << *i << " "; cerr << endl; }
template <class T> bool chmin(T &t, const T &f) { if (t > f) { t = f; return true; } return false; }
template <class T> bool chmax(T &t, const T &f) { if (t < f) { t = f; return true; } return false; }
#define COLOR(s) ("\x1b[" s "m")


using Double = double;
constexpr int ITER = 10000;

int N;
vector<int> A, B;
int S;

int dist[110][110];
int d[110][110];

int main() {
  for (; ~scanf("%d", &N); ) {
    A.resize(N - 1);
    B.resize(N - 1);
    for (int i = 0; i < N - 1; ++i) {
      scanf("%d%d", &A[i], &B[i]);
      --A[i];
      --B[i];
    }
    scanf("%d", &S);
    --S;

    for (int u = 0; u < N; ++u) for (int v = 0; v < N; ++v) {
      dist[u][v] = (u == v) ? 0 : N;
    }
    for (int i = 0; i < N - 1; ++i) {
      chmin(dist[A[i]][B[i]], 1);
      chmin(dist[B[i]][A[i]], 1);
    }
    for (int w = 0; w < N; ++w) for (int u = 0; u < N; ++u) for (int v = 0; v < N; ++v) {
      chmin(dist[u][v], dist[u][w] + dist[w][v]);
    }

    vector<int> deg(N, 0);
    for (int i = 0; i < N - 1; ++i) {
      ++deg[A[i]];
      ++deg[B[i]];
    }
    vector<int> ls;
    for (int u = 0; u < N; ++u) if (deg[u] == 1) {
      ls.push_back(u);
    }
    const int lsLen = ls.size();
    for (int i = 0; i < lsLen; ++i) for (int j = 0; j < lsLen; ++j) {
      d[i][j] = dist[ls[i]][ls[j]];
    }
    vector<Double> es(lsLen, 1.0), fs(lsLen);
    for (int iter = 1; iter <= ITER; ++iter) {
      Double sum = 0.0;
      for (int i = 0; i < lsLen; ++i) sum += es[i] = 1.0 / es[i];
      fill(fs.begin(), fs.end(), 0.0);
      for (int i = 0; i < lsLen; ++i) {
        Double t = 0.0;
        for (int j = 0; j < lsLen; ++j) if (i != j) {
          t += es[j] * d[i][j];
        }
        t += (lsLen - 2);
        fs[i] = t / (sum - es[i]);
      }
      es.swap(fs);
// if(!(iter&(iter-1)))cerr<<"iter = "<<iter<<", es = "<<es<<endl;
    }

    Double ans = 0.0;
    if (deg[S] == 1) {
      for (int i = 0; i < lsLen; ++i) if (S == ls[i]) {
        ans = es[i];
      }
    } else {
      Double sum = 0.0;
      for (int i = 0; i < lsLen; ++i) sum += es[i] = 1.0 / es[i];
      Double t = 0.0;
      for (int j = 0; j < lsLen; ++j) {
        t += es[j] * dist[S][ls[j]];
      }
      t += (lsLen - 1);
      ans = t / sum;
    }
    printf("%.12f\n", ans);
  }
  return 0;
}
```

## Markov Chain

```cpp
#include <cstdio>
#include <algorithm>

#define ll long long
#define db double
#define ull unsigned long long
#define uint unsigned int
#define FIO ""
#define dbug(...) fprintf(stderr, __VA_ARGS__)

template <typename Y> inline bool updmin(Y &a, Y b){if (a > b) {a = b; return 1;} return 0;}
template <typename Y> inline bool updmax(Y &a, Y b){if (a < b) {a = b; return 1;} return 0;}
template <typename Y> inline Y abs(Y a){if (a < 0) a = -a; return a;}
template <typename Y> inline Y sqr(Y a){return a * a;}

typedef std::pair<int, int> par;
#define fx first
#define fy second
#define mpar std::make_pair
#define pb push_back

int read() {
  int w = 1, q = 0, ch = ' ';
  for (; ch < '0' || ch > '9'; ch = getchar()) if (ch == '-') w = -1;
  for (; ch >= '0' && ch <= '9'; ch = getchar()) q = q * 10 + ch - 48;
  return q * w;
}

inline void FileIO(){freopen(FIO".in", "r", stdin); freopen(FIO".out", "w", stdout);}

const int N = 833;
double a[N][N], e[N][N], p[N], g[N], f[N][N];
int n, m, q, deg[N];
const double eps = 1e-9;
inline int sgn(const double &x) {
  if (x < -eps) return -1;
  if (x > eps) return 1;
  return 0;
}

void Gauss(int n, int m) {
  for (int i = 0, j; i < n; i++) {
    for (j = i; j < n && !sgn(a[j][i]); j++);
    if (j != i) {
      for (int k = 0; k < m; k++) {
        std::swap(a[i][k], a[j][k]);
      }
    }
    double d = a[i][i];
    if (!sgn(d)) break;
    for (int k = i; k < m; k++) {
      a[i][k] /= d;
    }
    for (j = 0; j < n; j++) if (j ^ i) {
      double mul = a[j][i];
      for (int k = i; k < m; k++) {
        a[j][k] -= mul * a[i][k];
      }
    }
  }
}

int main() {
#ifndef ONLINE_JUDGE
  freopen("101981F.in", "r", stdin);
#endif
  n = read();
  m = read();
  q = read();
  for (int i = 1; i <= m; i++) {
    int x = read(), y = read();
    e[x][y]++;
    deg[x]++;
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      e[i][j] /= deg[i];
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= n; j++) {
      a[i][j] = 0;
    }
  }
  for (int i = 0; i <= n; i++) {
    a[0][i] = 1;
  }
  for (int i = 0; i < n; i++) {
    a[i][i] = 1;
  }
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < n; j++) {
      a[i][j] -= e[j][i];
    }
  }
  Gauss(n, n + 1);
  for (int i = 0; i < n; i++) {
    p[i] = a[i][n];
    g[i] = 1 / p[i];
    //printf("g[%d] = %.3lf\n", i, g[i]);
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n << 1; j++) {
      a[i][j] = 0;
    }
  }
  for (int i = 0; i < n; i++) {
    a[i][i]++;
    for (int j = 0; j < n; j++) {
      a[i][j] -= e[i][j];
    }
  }
  for (int i = 0; i < n; i++) {
    a[i][i + n] -= g[i];
    for (int j = 0; j < n; j++) {
      a[i][j + n]++;
    }
  }
  Gauss(n, n << 1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      f[i][j] = a[i][j + n] - a[j][j + n];
    }
  }
  while (q--) {
    int k = read() - 1, x = read();
    if (!k) {
      puts("1");
    } else {
      double ans = 0;
      while (k--) {
        int y = read();
        ans += f[x][y];
        x = y;
      }
      printf("%.9lf\n", ans);
    }
  }
  return 0;
}
```
