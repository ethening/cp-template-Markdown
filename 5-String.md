# String

## KMP

```c++
int main(){
	string s;
	cin >> s;
	int n = (int)s.size();
	s = " " + s;

	vector<int> nxt(n + 1);
	for (int i = 2, p = 0; i <= n; i++) {
		while (p && s[p + 1] != s[i]) p = nxt[p];
		p = nxt[i] = p + (s[p + 1] == s[i]);
	}


	string t;
	cin >> t;
	int m = (int)t.size();
	t = " " + t;
	for (int i = 1, p = 0; i <= m; i++) {
		while (p && s[p + 1] != t[i]) p = nxt[p];
		p += (s[p + 1] == t[i]);
		if (p == n) {
			cout << i - n + 1 << "\n";
			p = nxt[p];
		}
	}

	for (int i = 1; i <= n; i++) {
		cout << nxt[i] << "\n";
	}

	return 0;
}
```

## Trie

```c++
struct Trie {
	const int S = 26;
	int N;

	int cnt = 0; // root = 0

	vector<vector<int>> tr;
	vector<int> pos;

	Trie() : N(0) {}
	Trie(int n) : N(2 * n + 5), tr(N, vector<int>(S)), pos(N) {}

	void insert(string s, int id) {
		int p = 0;
		for (char it: s) {
			if (!tr[p][it - 'a']) tr[p][it - 'a'] = ++cnt;
			p = tr[p][it - 'a'];
		}
	}
};
```

## GSAM

```c++
struct GSAM {
	const int S = 26;
	int N;

	vector<vector<int>> son;
    vector<int> len, fa;

	int cnt = 1; // root = 1

    GSAM() : N(0) {}
    GSAM(int n) : N(2 * n + 5), son(N, vector<int>(S)), len(N), fa(N) {}

	int insert(int p, int it) {
		int cur = ++cnt;
		len[cur] = len[p] + 1;
		while (!son[p][it]) son[p][it] = cur, p = fa[p];
		if (!p) return fa[cur] = 1, cur;
		int q = son[p][it];
		if (len[p] + 1 == len[q]) return fa[cur] = q, cur;
		int cl = ++cnt;
		son[cl] = son[q];
		len[cl] = len[p] + 1, fa[cl] = fa[q], fa[q] = fa[cur] = cl;
		while (son[p][it] == q) son[p][it] = cl, p = fa[p];
		return cur;
	}

	void build(Trie &T) {
		queue<int> Q;
		Q.push(0), T.pos[0] = 1;
		while (!Q.empty()) {
			int cur = Q.front();
			Q.pop();
			for (int i = 0; i < S; i++) {
				int p = T.tr[cur][i];
				if (!p) continue;
				T.pos[p] = insert(T.pos[cur], i), Q.push(p);
			}
		}
	}
};
```

## ACA

```c++
#include <bits/stdc++.h>

using i64 = long long;

struct AhoCorasick {
    static constexpr int ALPHABET = 26;
    struct Node {
        int len;
        int link;
        std::array<int, ALPHABET> next;
        Node() : link{}, next{} {}
    };

    std::vector<Node> t;

    AhoCorasick() {
        init();
    }

    void init() {
        t.assign(2, Node());
        t[0].next.fill(1);
        t[0].len = -1;
    }

    int newNode() {
        t.emplace_back();
        return t.size() - 1;
    }

    int add(const std::vector<int> &a) {
        int p = 1;
        for (auto x : a) {
            if (t[p].next[x] == 0) {
                t[p].next[x] = newNode();
                t[t[p].next[x]].len = t[p].len + 1;
            }
            p = t[p].next[x];
        }
        return p;
    }

    int add(const std::string &a, char offset = 'a') {
        std::vector<int> b(a.size());
        for (int i = 0; i < a.size(); i++) {
            b[i] = a[i] - offset;
        }
        return add(b);
    }

    void work() {
        std::queue<int> q;
        q.push(1);

        while (!q.empty()) {
            int x = q.front();
            q.pop();

            for (int i = 0; i < ALPHABET; i++) {
                if (t[x].next[i] == 0) {
                    t[x].next[i] = t[t[x].link].next[i];
                } else {
                    t[t[x].next[i]].link = t[t[x].link].next[i];
                    q.push(t[x].next[i]);
                }
            }
        }
    }

    int next(int p, int x) {
        return t[p].next[x];
    }

    int next(int p, char c, char offset = 'a') {
        return next(p, c - 'a');
    }

    int link(int p) {
        return t[p].link;
    }

    int len(int p) {
        return t[p].len;
    }

    int size() {
        return t.size();
    }
};
```

## Manacher

```c++
using i64 = long long;
std::vector<int> manacher(std::string s) {
    std::string t = "#";
    for (auto c : s) {
        t += c;
        t += '#';
    }
    int n = t.size();
    std::vector<int> r(n);
    for (int i = 0, j = 0; i < n; i++) {
        if (2 * j - i >= 0 && j + r[j] > i) {
            r[i] = std::min(r[2 * j - i], j + r[j] - i);
        }
        while (i - r[i] >= 0 && i + r[i] < n && t[i - r[i]] == t[i + r[i]]) {
            r[i] += 1;
        }
        if (i + r[i] > j + r[j]) {
            j = i;
        }
    }
    return r;
}
```


## SAM

```cpp
struct SuffixAutomaton {
    static constexpr int ALPHABET_SIZE = 26, N = 1e6;
    struct Node {
        int len;
        int link;
        int next[ALPHABET_SIZE];
        Node() : len(0), link(0), next{} {}
    } t[2 * N];
    int cntNodes;
    SuffixAutomaton() {
        cntNodes = 1;
        std::fill(t[0].next, t[0].next + ALPHABET_SIZE, 1);
        t[0].len = -1;
    }
    int extend(int p, int c) {
        if (t[p].next[c]) {
            int q = t[p].next[c];
            if (t[q].len == t[p].len + 1)
                return q;
            int r = ++cntNodes;
            t[r].len = t[p].len + 1;
            t[r].link = t[q].link;
            std::copy(t[q].next, t[q].next + ALPHABET_SIZE, t[r].next);
            t[q].link = r;
            while (t[p].next[c] == q) {
                t[p].next[c] = r;
                p = t[p].link;
            }
            return r;
        }
        int cur = ++cntNodes;
        t[cur].len = t[p].len + 1;
        while (!t[p].next[c]) {
            t[p].next[c] = cur;
            p = t[p].link;
        }
        t[cur].link = extend(p, c);
        return cur;
    }
};
```

## Suffix Array

```cpp

#include <bits/stdc++.h>

using i64 = long long;
struct SuffixArray {
    int n;
    std::vector<int> sa, rk, lc;
    SuffixArray(const std::string &s) {
        n = s.length();
        sa.resize(n);
        lc.resize(n - 1);
        rk.resize(n);
        std::iota(sa.begin(), sa.end(), 0);
        std::sort(sa.begin(), sa.end(), [&](int a, int b) {return s[a] < s[b];});
        rk[sa[0]] = 0;
        for (int i = 1; i < n; ++i)
            rk[sa[i]] = rk[sa[i - 1]] + (s[sa[i]] != s[sa[i - 1]]);
        int k = 1;
        std::vector<int> tmp, cnt(n);
        tmp.reserve(n);
        while (rk[sa[n - 1]] < n - 1) {
            tmp.clear();
            for (int i = 0; i < k; ++i)
                tmp.push_back(n - k + i);
            for (auto i : sa)
                if (i >= k)
                    tmp.push_back(i - k);
            std::fill(cnt.begin(), cnt.end(), 0);
            for (int i = 0; i < n; ++i)
                ++cnt[rk[i]];
            for (int i = 1; i < n; ++i)
                cnt[i] += cnt[i - 1];
            for (int i = n - 1; i >= 0; --i)
                sa[--cnt[rk[tmp[i]]]] = tmp[i];
            std::swap(rk, tmp);
            rk[sa[0]] = 0;
            for (int i = 1; i < n; ++i)
                rk[sa[i]] = rk[sa[i - 1]] + (tmp[sa[i - 1]] < tmp[sa[i]] || sa[i - 1] + k == n || tmp[sa[i - 1] + k] < tmp[sa[i] + k]);
            k *= 2;
        }
        for (int i = 0, j = 0; i < n; ++i) {
            if (rk[i] == 0) {
                j = 0;
            } else {
                for (j -= j > 0; i + j < n && sa[rk[i] - 1] + j < n && s[i + j] == s[sa[rk[i] - 1] + j]; )
                    ++j;
                lc[rk[i] - 1] = j;
            }
        }
    }
};
```
