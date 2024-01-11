# Misc

## Eval expression

```cpp
bool delim(char c) { return c == ' '; }

bool is_op(char c) { return c == '+' || c == '-' || c == '*' || c == '/'; }

bool is_unary(char c) { return c == '+' || c == '-'; }

int priority(char op) {
  if (op < 0)  // unary operator
    return 3;
  if (op == '+' || op == '-') return 1;
  if (op == '*' || op == '/') return 2;
  return -1;
}

void process_op(stack<int>& st, char op) {
  if (op < 0) {
    int l = st.top();
    st.pop();
    switch (-op) {
      case '+':
        st.push(l);
        break;
      case '-':
        st.push(-l);
        break;
    }
  } else {  // 取出栈顶元素，注意顺序
    int r = st.top();
    st.pop();
    int l = st.top();
    st.pop();
    switch (op) {
      case '+':
        st.push(l + r);
        break;
      case '-':
        st.push(l - r);
        break;
      case '*':
        st.push(l * r);
        break;
      case '/':
        st.push(l / r);
        break;
    }
  }
}

int evaluate(string& s) {
  stack<int> st;
  stack<char> op;
  bool may_be_unary = true;
  for (int i = 0; i < (int)s.size(); i++) {
    if (delim(s[i])) continue;

    if (s[i] == '(') {
      op.push('(');  // 2. 如果遇到左括号，那么将其放在运算符栈上
      may_be_unary = true;
    } else if (s[i] == ')') {  // 3. 如果遇到右括号，执行一对括号内的所有运算符
      while (op.top() != '(') {
        process_op(st, op.top());
        op.pop();  // 不断输出栈顶元素，直至遇到左括号
      }
      op.pop();  // 左括号出栈
      may_be_unary = false;
    } else if (is_op(s[i])) {  // 4. 如果遇到其他运算符
      char cur_op = s[i];
      if (may_be_unary && is_unary(cur_op)) cur_op = -cur_op;
      while (!op.empty() &&
             ((cur_op >= 0 && priority(op.top()) >= priority(cur_op)) ||
              (cur_op < 0 && priority(op.top()) > priority(cur_op)))) {
        process_op(st, op.top());
        op.pop();  // 不断输出所有运算优先级大于等于当前运算符的运算符
      }
      op.push(cur_op);  // 新的运算符入运算符栈
      may_be_unary = true;
    } else {  // 1. 如果遇到数字，直接输出该数字
      int number = 0;
      while (i < (int)s.size() && isalnum(s[i]))
        number = number * 10 + s[i++] - '0';
      --i;
      st.push(number);
      may_be_unary = false;
    }
  }

  while (!op.empty()) {
    process_op(st, op.top());
    op.pop();
  }
  return st.top();
}
```

## Digit DP

题目大意：给定一个区间 [l,r]，求其中满足条件 不含前导 0 且相邻两个数字相差至少为 2 的数字个数。

```cpp
int dfs(int x, int st, int op)  // op=1 =;op=0 <
{
  if (!x) return 1;
  if (!op && ~f[x][st]) return f[x][st];
  int maxx = op ? dim[x] : 9, ret = 0;
  for (int i = 0; i <= maxx; i++) {
    if (abs(st - i) < 2) continue;
    if (st == 11 && i == 0)
      ret += dfs(x - 1, 11, op & (i == maxx));
    else
      ret += dfs(x - 1, i, op & (i == maxx));
  }
  if (!op) f[x][st] = ret;
  return ret;
}

int solve(int x) {
  memset(f, -1, sizeof f);
  dim.clear();
  dim.push_back(-1);
  int t = x;
  while (x) {
    dim.push_back(x % 10);
    x /= 10;
  }
  return dfs(dim.size() - 1, 11, 1);
}
```

## Nanjing 2018 B Alien Trick

```cpp
#include "bits/stdc++.h"
#include <iomanip>
#define int long long
using namespace std;

using ll = long long;
using LL = long long;

const int N = 3e5 + 11;
const int INF = 1LL << 60;

int a[N], dp[N], cnt[N];

// 1-based


int ps1[N];
// ps2[N];

int s0(int l, int r){
	return r - l + 1;
}
int s1(int l, int r){
	return ps1[r] - ps1[l - 1];
}
// int s2(int l, int r){
// 	return ps2[r] - ps2[l - 1];
// }

int cost_L(int l, int m){
	return a[m] * s0(l, m - 1) - s1(l, m - 1);
}

int cost_R(int m, int r){
	return s1(m + 1, r) - a[m] * s0(m + 1, r);
}

int Cost(int l, int r){
	int m = (l + r) / 2;
	// cout << cost_L(l, m) << ' ' << cost_R(m, r) << '\n';
	return cost_L(l, m) + cost_R(m, r);
}

int n, k;

struct OptDS {
	struct Node {
		int l, r;
		int id;
	};

	vector<Node> v;

	int getOpt(int x) {
		int n = v.size();
		int L = 0;
		int R = n - 1;
		while (L <= R) {
			int mid = (L + R) / 2;
			if (v[mid].l <= x && x <= v[mid].r) {
				return v[mid].id;
			}
			else if (v[mid].l > x) {
				R = mid - 1;
			}
			else {
				L = mid + 1;
			}
		}
		return -1;
	}

	bool better(int id1, int id2, int x) {
		ll c1 = dp[id1] + Cost(id1 + 1, x);
		ll c2 = dp[id2] + Cost(id2 + 1, x);
		if (c1 != c2) return c1 <= c2;
		return cnt[id1] < cnt[id2];
	}

	void insert(int x) {
		while (!v.empty() && v.back().l > x && better(x, v.back().id, v.back().l)) v.pop_back();

		if (v.empty()) {
			v.push_back({x + 1, n, x});
			return;
		}

		if (better(v.back().id, x, n)) return;

		int L = max(x + 1, v.back().l);
		int R = n;

		while (L <= R) {
			int mid = (L + R) / 2;
			if (better(x, v.back().id, mid)) {
				R = mid - 1;
			}
			else {
				L = mid + 1;
			}
		}
		if (L <= n) {
			v.back().r = L - 1;
			v.push_back({L, n, x});
		}
	}

	void print() {
		for (auto [l, r, id]: v) {
			cout << "l: " << l << ", r: " << r << ", id: " << id << endl;
		}
	}
};

pair<ll, ll> calc(ll lambda){
	fill(dp + 1, dp + n + 1, INF);
	fill(cnt, cnt + n + 1, 0);

	OptDS DS;
	DS.insert(0);
	for (int i = 1; i <= n; i++) {
		int opt = DS.getOpt(i);
		dp[i] = dp[opt] + Cost(opt + 1, i) + lambda;
		cnt[i] = cnt[opt] + 1;

		DS.insert(i);

		// cout << "i: " << i << endl;
		// DS.print();
	}

	// for (int i = 1; i <= n; i++) {
	// 	// find opt
	// 	int opt = 0, cur = INF;
	// 	for (int j = 0; j < i; j++) {
	// 		int cost = dp[j] + Cost(j + 1, i) + lambda;
	// 		if (cost < cur) {
	// 			opt = j;
	// 			cur = cost;
	// 		}
	// 	}
	// 	dp[i] = dp[opt] + Cost(opt + 1, i) + lambda;
	// 	cnt[i] = cnt[opt] + 1;
	// }
	return {dp[n], cnt[n]};
}

void solve(int TC) {
	cin >> n >> k;
	for(int i = 1; i <= n; i++){
		cin >> a[i];
		ps1[i] = ps1[i - 1] + a[i];
		// ps2[i] = ps2[i - 1] + a[i] * i;
	}

	// printf("(%lld %lld) => %lld\n", 1, 3, Cost(1, 3));
	// printf("(%lld %lld) => %lld\n", 1, 4, Cost(1, 4));
	ll L = 0, R = 3 * 1E15;
	while (L <= R) {
		int mid = (L + R) / 2;
		auto [c, photo] = calc(mid);
		if (photo > k) {
			L = mid + 1;
		}
		else {
			R = mid - 1;
		}
	}
	auto [c, ph] = calc(L);
	// cout << c << ' ' << ph << endl;
	ll ans = c - k * L;
	cout << ans << '\n';
	return;
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
