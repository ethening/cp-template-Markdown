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
