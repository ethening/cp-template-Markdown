# Start

## Makefile

```Makefile
# ICPC
% : %.cpp
	g++ -Wall -Wfatal-errors -Wshadow -g -std=c++2a $< -o $@ -O2
```

## Code Runner

```sh
g++ -g -std=c++20 -fsanitize=address,undefined -fno-sanitize-recover -fstack-protector-all -O2 $fileName -o $fileNameWithoutExt
```

## Template

```cpp
#include "bits/stdc++.h"
using namespace std;
using ll = long long;
using pii = pair<int, int>;
using pil = pair<int, ll>;
using pli = pair<ll, int>;
using pll = pair<ll, ll>;

inline void yes() { cout  << "YES" << "\n"; return; }
inline void no() { cout << "NO" << "\n"; return; }

template<typename T>
bool chmin(T &x, T val) { if (val < x) { x = val; return true; } return false; }
template<typename T>
bool chmax(T &x, T val) { if (x < val) { x = val; return true; } return false; }

#define DEBUG 1
#define MULTI_TEST 1

void solve(int TC) {
	int n;
	cin >> n;
	vector<int> a(n);
	for (int i = 0; i < n; i++) {
		cin >> a[i];
	}
}

int main() {
	cin.tie(0)->sync_with_stdio(0);
	cout << fixed << setprecision(9);
	int t = 1;
	if (MULTI_TEST) cin >> t;
	for (int _ = 1; _ <= t; _++) {
		solve(_);
	}
}
```

## Pragma

```c++
#pragma GCC optimize("Ofast")
#pragma GCC target("avx2")
```
