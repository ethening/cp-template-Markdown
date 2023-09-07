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
