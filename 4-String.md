# String

## GSAM

### 约定

字符串个数为 $k$ 个，即 $S_1, S_2, S_3 \dots S_k$

约定字典树和广义后缀自动机的根节点为 $0$ 号节点

## 应用

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
