# String

## GSAM

### 约定

字符串个数为 $k$ 个，即 $S_1, S_2, S_3 \dots S_k$

约定字典树和广义后缀自动机的根节点为 $0$ 号节点

### 所有字符中不同子串个数

可以根据后缀自动机的性质得到，以点 $i$ 为结束节点的子串个数等于 $len[i] - len[link[i]]$

所以可以遍历所有的节点求和得到

例题：[【模板】广义后缀自动机（广义 SAM）](https://www.luogu.com.cn/problem/P6139)

```cpp
#include <bits/stdc++.h>
using namespace std;
const int MAXN = 2000000;  // 双倍字符串长度
const int CHAR_NUM = 30;   // 字符集个数，注意修改下方的 (-'a')

struct exSAM {
  int len[MAXN];             // 节点长度
  int link[MAXN];            // 后缀链接，link
  int next[MAXN][CHAR_NUM];  // 转移
  int tot;                   // 节点总数：[0, tot)

  void init() {  // 初始化函数
    tot = 1;
    link[0] = -1;
  }

  int insertSAM(int last, int c) {  // last 为父 c 为子
    int cur = next[last][c];
    if (len[cur]) return cur;
    len[cur] = len[last] + 1;
    int p = link[last];
    while (p != -1) {
      if (!next[p][c])
        next[p][c] = cur;
      else
        break;
      p = link[p];
    }
    if (p == -1) {
      link[cur] = 0;
      return cur;
    }
    int q = next[p][c];
    if (len[p] + 1 == len[q]) {
      link[cur] = q;
      return cur;
    }
    int clone = tot++;
    for (int i = 0; i < CHAR_NUM; ++i)
      next[clone][i] = len[next[q][i]] != 0 ? next[q][i] : 0;
    len[clone] = len[p] + 1;
    while (p != -1 && next[p][c] == q) {
      next[p][c] = clone;
      p = link[p];
    }
    link[clone] = link[q];
    link[cur] = clone;
    link[q] = clone;
    return cur;
  }

  int insertTrie(int cur, int c) {
    if (next[cur][c]) return next[cur][c];  // 已有该节点 直接返回
    return next[cur][c] = tot++;            // 无该节点 建立节点
  }

  void insert(const string &s) {
    int root = 0;
    for (auto ch : s) root = insertTrie(root, ch - 'a');
  }

  void insert(const char *s, int n) {
    int root = 0;
    for (int i = 0; i < n; ++i)
      root =
          insertTrie(root, s[i] - 'a');  // 一边插入一边更改所插入新节点的父节点
  }

  void build() {
    queue<pair<int, int>> q;
    for (int i = 0; i < 26; ++i)
      if (next[0][i]) q.push({i, 0});
    while (!q.empty()) {  // 广搜遍历
      auto item = q.front();
      q.pop();
      auto last = insertSAM(item.second, item.first);
      for (int i = 0; i < 26; ++i)
        if (next[last][i]) q.push({i, last});
    }
  }
} exSam;

char s[1000100];

int main() {
  int n;
  cin >> n;
  exSam.init();
  for (int i = 0; i < n; ++i) {
    cin >> s;
    int len = strlen(s);
    exSam.insert(s, len);
  }
  exSam.build();
  long long ans = 0;
  for (int i = 1; i < exSam.tot; ++i) {
    ans += exSam.len[i] - exSam.len[exSam.link[i]];
  }
  cout << ans << endl;
}
```

### 多个字符串间的最长公共子串

我们需要对每个节点建立一个长度为 $k$ 的数组 `flag`（对于本题而言，可以仅为标记数组，若需要求出此子串的个数，则需要改成计数数组）

在字典树插入字符串时，对所有节点进行计数，保存在当前字符串所在的数组

然后按照 `len` 递减的顺序遍历，通过后缀链接将当前节点的 `flag` 与其他节点的合并

遍历所有的节点，找到一个 `len` 最大且满足对于所有的 `k`，其 `flag` 的值均为非 $0$ 的节点，此节点的 $len$ 即为解

例题：[SPOJ Longest Common Substring II](https://www.spoj.com/problems/LCS2/)

```c++
#include <bits/stdc++.h>
using namespace std;

const int MAXN = 2000000;  // 双倍字符串长度
const int CHAR_NUM = 30;   // 字符集个数，注意修改下方的 (-'a')
const int NUM = 15;        // 字符串个数

struct exSAM {
  int len[MAXN];             // 节点长度
  int link[MAXN];            // 后缀链接，link
  int next[MAXN][CHAR_NUM];  // 转移
  int tot;                   // 节点总数：[0, tot)
  int lenSorted[MAXN];   // 按照 len 排序后的数组，仅排序 [1, tot)
                         // 部分，最终下标范围 [0, tot - 1)
  int sizeC[MAXN][NUM];  // 表示某个字符串的子串个数
  int curString;         // 字符串实际个数
  /**
   * 计数排序使用的辅助空间数组
   */
  int lc[MAXN];  // 统计个数

  void init() {
    tot = 1;
    link[0] = -1;
  }

  int insertSAM(int last, int c) {
    int cur = next[last][c];
    len[cur] = len[last] + 1;
    int p = link[last];
    while (p != -1) {
      if (!next[p][c])
        next[p][c] = cur;
      else
        break;
      p = link[p];
    }
    if (p == -1) {
      link[cur] = 0;
      return cur;
    }
    int q = next[p][c];
    if (len[p] + 1 == len[q]) {
      link[cur] = q;
      return cur;
    }
    int clone = tot++;
    for (int i = 0; i < CHAR_NUM; ++i)
      next[clone][i] = len[next[q][i]] != 0 ? next[q][i] : 0;
    len[clone] = len[p] + 1;
    while (p != -1 && next[p][c] == q) {
      next[p][c] = clone;
      p = link[p];
    }
    link[clone] = link[q];
    link[cur] = clone;
    link[q] = clone;
    return cur;
  }

  int insertTrie(int cur, int c) {
    if (!next[cur][c]) next[cur][c] = tot++;
    sizeC[next[cur][c]][curString]++;
    return next[cur][c];
  }

  void insert(const string &s) {
    int root = 0;
    for (auto ch : s) root = insertTrie(root, ch - 'a');
    curString++;
  }

  void insert(const char *s, int n) {
    int root = 0;
    for (int i = 0; i < n; ++i) root = insertTrie(root, s[i] - 'a');
    curString++;
  }

  void build() {
    queue<pair<int, int>> q;
    for (int i = 0; i < 26; ++i)
      if (next[0][i]) q.push({i, 0});
    while (!q.empty()) {  // 广搜遍历
      auto item = q.front();
      q.pop();
      auto last = insertSAM(item.second, item.first);
      for (int i = 0; i < 26; ++i)
        if (next[last][i]) q.push({i, last});
    }
  }

  void sortLen() {
    for (int i = 1; i < tot; ++i) lc[i] = 0;
    for (int i = 1; i < tot; ++i) lc[len[i]]++;
    for (int i = 2; i < tot; ++i) lc[i] += lc[i - 1];
    for (int i = 1; i < tot; ++i) lenSorted[--lc[len[i]]] = i;
  }

  void getSizeLen() {
    for (int i = tot - 2; i >= 0; --i)
      for (int j = 0; j < curString; ++j)
        sizeC[link[lenSorted[i]]][j] += sizeC[lenSorted[i]][j];
  }
} exSam;

int main() {
  exSam.init();  // 初始化
  string s;
  while (cin >> s) exSam.insert(s);
  exSam.build();
  exSam.sortLen();
  exSam.getSizeLen();
  int ans = 0;
  for (int i = 0; i < exSam.tot; ++i) {
    bool flag = true;
    for (int j = 0; j < exSam.curString; ++j) {
      if (!exSam.sizeC[i][j]) {
        flag = false;
        break;
      }
    }
    if (flag) ans = max(ans, exSam.len[i]);
  }
  cout << ans << endl;
}
```

## Lyndon

首先我们介绍 Lyndon 分解的概念。

Lyndon 串：对于字符串 $s$，如果 $s$ 的字典序严格小于 $s$ 的所有后缀的字典序，我们称 $s$ 是简单串，或者 **Lyndon 串**。举一些例子，`a`,`b`,`ab`,`aab`,`abb`,`ababb`,`abcd` 都是 Lyndon 串。当且仅当 $s$ 的字典序严格小于它的所有非平凡的（非平凡：非空且不同于自身）循环同构串时，$s$ 才是 Lyndon 串。

Lyndon 分解：串 $s$ 的 Lyndon 分解记为 $s=w_1w_2\cdots w_k$，其中所有 $w_i$ 为简单串，并且他们的字典序按照非严格单减排序，即 $w_1\ge w_2\ge\cdots\ge w_k$。可以发现，这样的分解存在且唯一。

## Duval 算法

### 解释

Duval 可以在 $O(n)$ 的时间内求出一个串的 Lyndon 分解。

首先我们介绍另外一个概念：如果一个字符串 $t$ 能够分解为 $t=ww\cdots\overline{w}$ 的形式，其中 $w$ 是一个 Lyndon 串，而 $\overline{w}$ 是 $w$ 的前缀（$\overline{w}$ 可能是空串），那么称 $t$ 是近似简单串（pre-simple），或者近似 Lyndon 串。一个 Lyndon 串也是近似 Lyndon 串。

Duval 算法运用了贪心的思想。算法过程中我们把串 $s$ 分成三个部分 $s=s_1s_2s_3$，其中 $s_1$ 是一个 Lyndon 串，它的 Lyndon 分解已经记录；$s_2$ 是一个近似 Lyndon 串；$s_3$ 是未处理的部分。

### 过程

整体描述一下，该算法每一次尝试将 $s_3$ 的首字符添加到 $s_2$ 的末尾。如果 $s_2$ 不再是近似 Lyndon 串，那么我们就可以将 $s_2$ 截出一部分前缀（即 Lyndon 分解）接在 $s_1$ 末尾。

我们来更详细地解释一下算法的过程。定义一个指针 $i$ 指向 $s_2$ 的首字符，则 $i$ 从 $1$ 遍历到 $n$（字符串长度）。在循环的过程中我们定义另一个指针 $j$ 指向 $s_3$ 的首字符，指针 $k$ 指向 $s_2$ 中我们当前考虑的字符（意义是 $j$ 在 $s_2$ 的上一个循环节中对应的字符）。我们的目标是将 $s[j]$ 添加到 $s_2$ 的末尾，这就需要将 $s[j]$ 与 $s[k]$ 做比较：

1. 如果 $s[j]=s[k]$，则将 $s[j]$ 添加到 $s_2$ 末尾不会影响它的近似简单性。于是我们只需要让指针 $j,k$ 自增（移向下一位）即可。
2. 如果 $s[j]>s[k]$，那么 $s_2s[j]$ 就变成了一个 Lyndon 串，于是我们将指针 $j$ 自增，而让 $k$ 指向 $s_2$ 的首字符，这样 $s_2$ 就变成了一个循环次数为 1 的新 Lyndon 串了。
3. 如果 $s[j]<s[k]$，则 $s_2s[j]$ 就不是一个近似简单串了，那么我们就要把 $s_2$ 分解出它的一个 Lyndon 子串，这个 Lyndon 子串的长度将是 $j-k$，即它的一个循环节。然后把 $s_2$ 变成分解完以后剩下的部分，继续循环下去（注意，这个情况下我们没有改变指针 $j,k$），直到循环节被截完。对于剩余部分，我们只需要将进度「回退」到剩余部分的开头即可。

### 实现

下面的代码返回串 $s$ 的 Lyndon 分解方案。

```cpp
    // duval_algorithm
    vector<string> duval(string const& s) {
      int n = s.size(), i = 0;
      vector<string> factorization;
      while (i < n) {
        int j = i + 1, k = i;
        while (j < n && s[k] <= s[j]) {
          if (s[k] < s[j])
            k = i;
          else
            k++;
          j++;
        }
        while (i <= k) {
          factorization.push_back(s.substr(i, j - k));
          i += j - k;
        }
      }
      return factorization;
    }
```

### 复杂度分析

接下来我们证明一下这个算法的复杂度。

外层的循环次数不超过 $n$，因为每一次 $i$ 都会增加。第二个内层循环也是 $O(n)$ 的，因为它只记录 Lyndon 分解的方案。接下来我们分析一下内层循环。很容易发现，每一次在外层循环中找到的 Lyndon 串是比我们所比较过的剩余的串要长的，因此剩余的串的长度和要小于 $n$，于是我们最多在内层循环 $O(n)$ 次。事实上循环的总次数不超过 $4n-3$，时间复杂度为 $O(n)$。

### 最小表示法（Finding the smallest cyclic shift）

对于长度为 $n$ 的串 $s$，我们可以通过上述算法寻找该串的最小表示法。

我们构建串 $ss$ 的 Lyndon 分解，然后寻找这个分解中的一个 Lyndon 串 $t$，使得它的起点小于 $n$ 且终点大于等于 $n$。可以很容易地使用 Lyndon 分解的性质证明，子串 $t$ 的首字符就是 $s$ 的最小表示法的首字符，即我们沿着 $t$ 的开头往后 $n$ 个字符组成的串就是 $s$ 的最小表示法。

于是我们在分解的过程中记录每一次的近似 Lyndon 串的开头即可。

```cpp
// smallest_cyclic_string
string min_cyclic_string(string s) {
    s += s;
    int n = s.size();
    int i = 0, ans = 0;
    while (i < n / 2) {
        ans = i;
        int j = i + 1, k = i;
        while (j < n && s[k] <= s[j]) {
            if (s[k] < s[j]) k = i;
            else k++;
            j++;
        }
        while (i <= k) i += j - k;
    }
    return s.substr(ans, n / 2);
}
```

## Qinhuangdao C Palindrome

Erase l to r of string S costing r - l + 1, min cost to make it palindrome, and numbers of ways
```cpp
#include "bits/stdc++.h"
#include <iostream>
using namespace std;


using ll = long long;
using LL = long long;

using pii = pair<int, int>;
using pil = pair<int, ll>;
using pli = pair<ll, int>;
using pll = pair<ll, ll>;

inline int lg(int __n) {
	return sizeof(int) * __CHAR_BIT__	- 1 - __builtin_clz(__n);
}

template<size_t CHARSETSZ = 26>
struct PalindromesAutomaton {
	/**
	 * sz: 节点数
	 * root0: 偶数回文串的根
	 * root1：奇数根
	 * link: 最大回文后缀指针。也叫后缀指针或者失配指针。link指针反向构成一棵树
	 * len: 每个节点对应回文串的长度
	 * sum: 每个节点对应回文串包含的回文后缀数目
	 * ch: 儿子指针
	 */
	int sz, root0, root1;
	vector<int> link, len, rep;
	vector<array<int, CHARSETSZ>> ch;

	vector<vector<int>> sp;

	PalindromesAutomaton(int n) :sz(0), ch(n), link(n), len(n), rep(n), sp(20) {
		root0 = sz++;
		root1 = sz++;
		len[root0] = 0;
		// 奇数根的长度初始化为-1
		len[root1] = -1;
		// 两个根的link都指向奇数根
		link[root0] = root1;
		link[root1] = root1;
	}

	int getId(char c) {
		return c - 'a';
	}

	inline int getLink(const string&s, int i, int u){
		while (i - len[u] - 1 < 0 || s[i] != s[i - len[u] - 1]) u = link[u];
		return u;
	}

	void build(string& s) {
		int ans = 0, last = root1;
		for (int i = 0; i < s.length(); i++) {
			int c = getId(s[i]);
			int u = getLink(s, i, last);

			if (!ch[u][c]) {
				int cur = sz++;
				int v = getLink(s, i, link[u]);
				len[cur] = len[u] + 2;
				link[cur] = ch[v][c];
				ch[u][c] = cur;
			}

			last = ch[u][c];
			rep[i] = last;
			// cout << "i: " << i << " rep: " << rep[i] << " " << " len: " << len[rep[i]] << " " << link[rep[i]] << " " << len[link[rep[i]]] << endl;
		}
	}

	void buildsp() {
		for (int i = 0; i < 20; i++) sp[i].resize(sz);
		for (int i = 0; i < sz; i++) {
			sp[0][i] = link[i];
		}
		for (int i = 1; i < 20; i++) {
			for (int j = 0; j < sz; j++) {
				sp[i][j] = sp[i - 1][sp[i - 1][j]];
			}
		}
	}

	int find(int endpos, int mxlen) {
		int x = rep[endpos];
		// cout << "*" << endpos << " " << x << " " << len[x] << endl;
		for (int i = 19; i >= 0; i--) {
			if (len[sp[i][x]] >= mxlen) {
				x = sp[i][x];
				// cout << "#" << x << " " << len[x] << endl;
			}
		}
		if (len[x] > mxlen) {
			x = sp[0][x];
		}
		// cout << "**" << len[x] << endl;
		return len[x];
	}
};

typedef uint64_t ull;
struct H {
	ull x; H(ull x=0) : x(x) {}
	H operator+(H o) { return x + o.x + (x + o.x < x); }
	H operator-(H o) { return *this + ~o.x; }
	H operator*(H o) { auto m = (__uint128_t)x * o.x;
		return H((ull)m) + (ull)(m >> 64); }
	ull get() const { return x + !~x; }
	bool operator==(H o) const { return get() == o.get(); }
	bool operator<(H o) const { return get() < o.get(); }
};
static const H C = (ll)1e11+3; // (order ~ 3e9; random also ok)

struct HashInterval {
	vector<H> ha, pw;
	HashInterval(string& str) : ha(size(str)+1), pw(ha) {
		pw[0] = 1;
		for (int i = 0; i < size(str); i++) {
			ha[i+1] = ha[i] * C + str[i],
			pw[i+1] = pw[i] * C;
		}
	}
	H hashInterval(int a, int b) { // hash [a, b)
		return ha[b] - ha[a] * pw[b - a];
	}
};

void solve(int TC) {
	int n;
	cin >> n;
	string s;
	cin >> s;
	string r = s;
	reverse(r.begin(), r.end());

	HashInterval S(s), R(r);

	PalindromesAutomaton SP(n + 5), RP(n + 5);
	SP.build(s);
	SP.buildsp();
	RP.build(r);
	RP.buildsp();

	int q;
	cin >> q;
	for (int i = 0; i < q; i++) {
		int a, b;
		cin >> a >> b;
		--a, --b;

		auto check = [&](int len) -> bool {
			if (S.hashInterval(a, a + len) == R.hashInterval(n - 1 - b, n - 1 - b + len)) return true;
			else return false;
		};

		int len{};
		{
			int l = 1;
			int r = (b - a + 1);
			while (l <= r) {
				int mid = (l + r) / 2;
				if (check(mid)) {
					l = mid + 1;
				}
				else {
					r = mid - 1;
				}
			}

			if (r == (b - a + 1)) {
				cout << 0 << " " << 0 << "\n";
				continue;
			}
			len = r;
		}

		// [a + len .. b - len]

		int rsufplen = SP.find(b - len, (b - len) - (a + len) + 1);
		// [a + len .. b - len - rsufplen]

		int lpreplen = RP.find((n - 1) - (a + len), (b - len) - (a + len) + 1);
		// [a + len + lpreplen .. b - len]

		// cout << a + len << " - " << b - len - rsufplen << endl;
		// cout << a + len + lpreplen << " - " << b - len << endl;

		int len1 = (b - len - rsufplen) - (a + len) + 1;
		int len2 = (b - len) - (a + len + lpreplen) + 1;

		int alen = min(len1, len2);

		int ways = 0;

		if (len1 == alen) {
			// find lcp(a + len - 1 .. a + len - k, b - len - rsufplen .. b - len - rsufplen - k + 1)

			{
				auto check = [&](int l) -> bool {
					if (R.hashInterval(n - 1 - (a + len - 1), n - 1 - (a + len - 1) + l) == R.hashInterval(n - 1 - (b - len - rsufplen), n - 1 - (b - len - rsufplen) + l)) return true;
					else return false;
				};
				int l = 1;
				int r = len;
				while (l <= r) {
					int mid = (l + r) / 2;
					if (check(mid)) {
						l = mid + 1;
					}
					else {
						r = mid - 1;
					}
				}

				ways += r + 1;
			}

		}

		if (len2 == alen) {

			// find lcp(b - len + 1 .. b - len + k, a + len + lpreplen .. a + len + lpreplen + k - 1)
			{
				auto check = [&](int l) -> bool {
					if (S.hashInterval(b - len + 1, b - len + 1 + l) == S.hashInterval(a + len + lpreplen, a + len + lpreplen + l)) return true;
					else return false;
				};
				int l = 1;
				int r = len;
				while (l <= r) {
					int mid = (l + r) / 2;
					if (check(mid)) {
						l = mid + 1;
					}
					else {
						r = mid - 1;
					}
				}

				ways += r + 1;
			}
		}

		// cout << "ANS: ";
		cout << alen << " " << ways << "\n";
	}
}


int main() {
	cout << fixed << setprecision(9);
	cin.tie(0)->sync_with_stdio(0);
	int t = 1;
	// cin >> t;
	for (int i = 1; i <= t; i++) {
		solve(i);
	}
}
```
