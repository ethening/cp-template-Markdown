# Graph

## LGV Lemma

## 简介

Lindström–Gessel–Viennot lemma，即 LGV 引理，可以用来处理有向无环图上不相交路径计数等问题。

LGV 引理仅适用于 **有向无环图**。

题意：有一个 $n\times m$ 的格点棋盘，其中某些格子可走，某些格子不可走。有一只海龟从 $(x, y)$ 只能走到 $(x+1, y)$ 和 $(x, y+1)$ 的位置，求海龟从 $(1, 1)$ 到 $(n, m)$ 的不相交路径数对 $10^9+7$ 取模之后的结果。$2\le n,m\le3000$。

比较直接的 LGV 引理的应用。考虑所有合法路径，发现从 $(1,1)$ 出发一定要经过 $A=\{(1,2), (2,1)\}$，而到达终点一定要经过 $B=\{(n-1, m), (n, m-1)\}$，则 $A, B$ 可立即选定。应用 LGV 引理可得答案为：

$$
\begin{vmatrix}
f(a_1, b_1) & f(a_1, b_2) \\
f(a_2, b_1) & f(a_2, b_2)
\end{vmatrix} = f(a_1, b_1)\times f(a_2, b_2) - f(a_1, b_2)\times f(a_2, b_1)
$$

其中 $f(a, b)$ 为图上 $a\rightarrow b$ 的路径数，带有障碍格点的路径计数问题可以直接做一个 $O(nm)$ 的 dp，则 $f$ 易求。最终复杂度 $O(nm)$。

```cpp
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;

using ll = long long;
const int MOD = 1e9 + 7;
const int SIZE = 3010;

char board[SIZE][SIZE];
int dp[SIZE][SIZE];

int f(int x1, int y1, int x2, int y2) {
  memset(dp, 0, sizeof dp);

  dp[x1][y1] = board[x1][y1] == '.';
  for (int i = 1; i <= x2; i++) {
    for (int j = 1; j <= y2; j++) {
      if (board[i][j] == '#') {
        continue;
      }
      dp[i][j] = (dp[i][j] + dp[i - 1][j]) % MOD;
      dp[i][j] = (dp[i][j] + dp[i][j - 1]) % MOD;
    }
  }
  return dp[x2][y2] % MOD;
}

int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);
  cout.tie(nullptr);

  int n, m;
  cin >> n >> m;

  for (int i = 1; i <= n; i++) {
    cin >> (board[i] + 1);
  }

  ll f11 = f(1, 2, n - 1, m);
  ll f12 = f(1, 2, n, m - 1);
  ll f21 = f(2, 1, n - 1, m);
  ll f22 = f(2, 1, n, m - 1);

  ll ans = ((f11 * f22) % MOD - (f12 * f21) % MOD + MOD) % MOD;
  cout << ans << '\n';

  return 0;
}
```

题意：有一个 $n\times n$ 的棋盘，一个棋子从 $(x, y)$ 只能走到 $(x, y+1)$ 或 $(x + 1, y)$，有 $k$ 个棋子，一开始第 $i$ 个棋子放在 $(1, a_i)$，最终要到 $(n, b_i)$，路径要两两不相交，求方案数对 $10^9+7$ 取模。$1\le n\le 10^5$，$1\le k\le 100$，保证 $1\le a_1<a_2<\dots<a_n\le n$，$1\le b_1<b_2<\dots<b_n\le n$。

观察到如果路径不相交就一定是 $a_i$ 到 $b_i$，因此 LGV 引理中一定有 $\sigma(S)_i=i$，不需要考虑符号问题。边权设为 $1$，直接套用引理即可。

从 $(1, a_i)$ 到 $(n, b_j)$ 的路径条数相当于从 $n-1+b_j-a_i$ 步中选 $n-1$ 步向下走，所以 $e(A_i, B_j)=\binom{n-1+b_j-a_i}{n-1}$。

行列式可以使用高斯消元求。

复杂度为 $O(n+k(k^2 + \log p))$，其中 $\log p$ 是求逆元复杂度。

```cpp
#include <algorithm>
#include <cstdio>

typedef long long ll;

const int K = 105;
const int N = 100005;
const int mod = 1e9 + 7;

int T, n, k, a[K], b[K], fact[N << 1], m[K][K];

int qpow(int x, int y) {
  int out = 1;
  while (y) {
    if (y & 1) out = (ll)out * x % mod;
    x = (ll)x * x % mod;
    y >>= 1;
  }
  return out;
}

int c(int x, int y) {
  return (ll)fact[x] * qpow(fact[y], mod - 2) % mod *
         qpow(fact[x - y], mod - 2) % mod;
}

int main() {
  fact[0] = 1;
  for (int i = 1; i < N * 2; ++i) fact[i] = (ll)fact[i - 1] * i % mod;

  scanf("%d", &T);

  while (T--) {
    scanf("%d%d", &n, &k);

    for (int i = 1; i <= k; ++i) scanf("%d", a + i);
    for (int i = 1; i <= k; ++i) scanf("%d", b + i);

    for (int i = 1; i <= k; ++i) {
      for (int j = 1; j <= k; ++j) {
        if (a[i] <= b[j])
          m[i][j] = c(b[j] - a[i] + n - 1, n - 1);
        else
          m[i][j] = 0;
      }
    }

    for (int i = 1; i < k; ++i) {
      if (!m[i][i]) {
        for (int j = i + 1; j <= k; ++j) {
          if (m[j][i]) {
            std::swap(m[i], m[j]);
            break;
          }
        }
      }
      if (!m[i][i]) continue;
      int inv = qpow(m[i][i], mod - 2);
      for (int j = i + 1; j <= k; ++j) {
        if (!m[j][i]) continue;
        int mul = (ll)m[j][i] * inv % mod;
        for (int p = i; p <= k; ++p) {
          m[j][p] = (m[j][p] - (ll)m[i][p] * mul % mod + mod) % mod;
        }
      }
    }

    int ans = 1;

    for (int i = 1; i <= k; ++i) ans = (ll)ans * m[i][i] % mod;

    printf("%d\n", ans);
  }

  return 0;
}
```


## Flow with bound

在阅读这篇文章之前请先阅读 [最大流](./max-flow.md) 并确保自己熟练掌握最大流算法。

### 概述

上下界网络流本质是给流量网络的每一条边设置了流量上界 $c(u,v)$ 和流量下界 $b(u,v)$。也就是说，一种可行的流必须满足 $b(u,v) \leq f(u,v) \leq c(u,v)$。同时必须满足除了源点和汇点之外的其余点流量平衡。

根据题目要求，我们可以使用上下界网络流解决不同问题。

### 无源汇上下界可行流

给定无源汇流量网络 $G$。询问是否存在一种标定每条边流量的方式，使得每条边流量满足上下界同时每一个点流量平衡。

不妨假设每条边已经流了 $b(u,v)$ 的流量，设其为初始流。同时我们在新图中加入 $u$ 连向 $v$ 的流量为 $c(u,v) - b(u,v)$ 的边。考虑在新图上进行调整。

由于最大流需要满足初始流量平衡条件（最大流可以看成是下界为 $0$ 的上下界最大流），但是构造出来的初始流很有可能不满足初始流量平衡。假设一个点初始流入流量减初始流出流量为 $M$。

若 $M=0$，此时流量平衡，不需要附加边。

若 $M>0$，此时入流量过大，需要新建附加源点 $S'$，$S'$ 向其连流量为 $M$ 的附加边。

若 $M<0$，此时出流量过大，需要新建附加汇点 $T'$，其向 $T'$ 连流量为 $-M$ 的附加边。

如果附加边满流，说明这一个点的流量平衡条件可以满足，否则这个点的流量平衡条件不满足。（因为原图加上附加流之后才会满足原图中的流量平衡。）

在建图完毕之后跑 $S'$ 到 $T'$ 的最大流，若 $S'$ 连出去的边全部满流，则存在可行流，否则不存在。

### 有源汇上下界可行流

给定有源汇流量网络 $G$。询问是否存在一种标定每条边流量的方式，使得每条边流量满足上下界同时除了源点和汇点每一个点流量平衡。

假设源点为 $S$，汇点为 $T$。

则我们可以加入一条 $T$ 到 $S$ 的上界为 $\infty$，下界为 $0$ 的边转化为无源汇上下界可行流问题。

若有解，则 $S$ 到 $T$ 的可行流流量等于 $T$ 到 $S$ 的附加边的流量。

### 有源汇上下界最大流

给定有源汇流量网络 $G$。询问是否存在一种标定每条边流量的方式，使得每条边流量满足上下界同时除了源点和汇点每一个点流量平衡。如果存在，询问满足标定的最大流量。

我们找到网络上的任意一个可行流。如果找不到解就可以直接结束。

否则我们考虑删去所有附加边之后的残量网络并且在网络上进行调整。

我们在残量网络上再跑一次 $S$ 到 $T$ 的最大流，将可行流流量和最大流流量相加即为答案。

"一个非常易错的问题"
$S$ 到 $T$ 的最大流直接在跑完有源汇上下界可行的残量网络上跑。
千万不可以在原来的流量网络上跑。

### 有源汇上下界最小流

给定有源汇流量网络 $G$。询问是否存在一种标定每条边流量的方式，使得每条边流量满足上下界同时除了源点和汇点每一个点流量平衡。如果存在，询问满足标定的最小流量。

类似的，我们考虑将残量网络中不需要的流退掉。

我们找到网络上的任意一个可行流。如果找不到解就可以直接结束。

否则我们考虑删去所有附加边之后的残量网络。

我们在残量网络上再跑一次 $T$ 到 $S$ 的最大流，将可行流流量减去最大流流量即为答案。

[AHOI 2014 支线剧情](https://loj.ac/problem/2226)
    对于每条 $x$ 到 $y$ 花费 $v$ 的剧情边设上界为 $\infty$, 下界为 $1$。
    对于每个点，向 $T$ 连边权 $c$, 上界 $\infty$, 下界为 $1$。
    $S$ 点为 $1$ 号节点。
    跑一次 上下界带源汇最小费用可行流 即可。
    因为最小费用可行流解法与最小可行流类似，这里不再展开。
