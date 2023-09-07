# Data Structure

## Segment Tree

```cpp
template<class Info>
struct SegmentTree {
	int n;
	std::vector<Info> info;

	SegmentTree() : n(0) {}
	SegmentTree(int n_, Info v_ = Info()) {
		init(n_, v_);
	}
	template<class T>
	SegmentTree(std::vector<T> init_) {
		init(init_);
	}

	void init(int n_, Info v_ = Info()) {
		init(std::vector(n_, v_));
	}
	template<class T>
	void init(std::vector<T> init_) {
		n = init_.size();
		info.assign(4 << std::__lg(n), Info());
		std::function<void(int, int, int)> build = [&](int p, int l, int r) {
			if (r - l == 1) {
				info[p] = init_[l];
				return;
			}
			int m = (l + r) / 2;
			build(2 * p, l, m);
			build(2 * p + 1, m, r);
			pull(p);
		};
		build(1, 0, n);
	}
	void pull(int p) {
		info[p] = info[2 * p] + info[2 * p + 1];
	}
	void modify(int p, int l, int r, int x, const Info &v) {
		if (r - l == 1) {
			info[p] = v;
			return;
		}
		int m = (l + r) / 2;
		if (x < m) {
			modify(2 * p, l, m, x, v);
		} else {
			modify(2 * p + 1, m, r, x, v);
		}
		pull(p);
	}
	void modify(int p, const Info &v) {
		modify(1, 0, n, p, v);
	}
	Info rangeQuery(int p, int l, int r, int x, int y) {
		if (l >= y || r <= x) {
			return Info();
		}
		if (l >= x && r <= y) {
			return info[p];
		}
		int m = (l + r) / 2;
		return rangeQuery(2 * p, l, m, x, y) + rangeQuery(2 * p + 1, m, r, x, y);
	}
	Info rangeQuery(int l, int r) {
		return rangeQuery(1, 0, n, l, r);
	}
	template<class F>
	int findFirst(int p, int l, int r, int x, int y, F pred) {
		if (l >= y || r <= x || !pred(info[p])) {
			return -1;
		}
		if (r - l == 1) {
			return l;
		}
		int m = (l + r) / 2;
		int res = findFirst(2 * p, l, m, x, y, pred);
		if (res == -1) {
			res = findFirst(2 * p + 1, m, r, x, y, pred);
		}
		return res;
	}
	template<class F>
	int findFirst(int l, int r, F pred) {
		return findFirst(1, 0, n, l, r, pred);
	}
	template<class F>
	int findLast(int p, int l, int r, int x, int y, F pred) {
		if (l >= y || r <= x || !pred(info[p])) {
			return -1;
		}
		if (r - l == 1) {
			return l;
		}
		int m = (l + r) / 2;
		int res = findLast(2 * p + 1, m, r, x, y, pred);
		if (res == -1) {
			res = findLast(2 * p, l, m, x, y, pred);
		}
		return res;
	}
	template<class F>
	int findLast(int l, int r, F pred) {
		return findLast(1, 0, n, l, r, pred);
	}
};

constexpr ll inf = 1E18;

struct Info {
	ll cnt = 0;
	ll sum = 0;
	ll min = inf;
};

Info operator+(Info a, Info b) {
	Info c;
	c.cnt = a.cnt + b.cnt;
	c.sum = a.sum + b.sum;
	c.min = std::min(a.min, b.min);
	return c;
}
```

## Lazy Segment Tree

```cpp
template<class Info, class Tag>
struct LazySegmentTree {
	int n;
	std::vector<Info> info;
	std::vector<Tag> tag;
	LazySegmentTree() : n(0) {}
	LazySegmentTree(int n_, Info v_ = Info()) {
		init(n_, v_);
	}
	template<class T>
	LazySegmentTree(std::vector<T> init_) {
		init(init_);
	}
	void init(int n_, Info v_ = Info()) {
		init(std::vector(n_, v_));
	}
	template<class T>
	void init(std::vector<T> init_) {
		n = init_.size();
		info.assign(4 << std::__lg(n), Info());
		tag.assign(4 << std::__lg(n), Tag());
		std::function<void(int, int, int)> build = [&](int p, int l, int r) {
			if (r - l == 1) {
				info[p] = init_[l];
				return;
			}
			int m = (l + r) / 2;
			build(2 * p, l, m);
			build(2 * p + 1, m, r);
			pull(p);
		};
		build(1, 0, n);
	}
	void pull(int p) {
		info[p] = info[2 * p] + info[2 * p + 1];
	}
	void apply(int p, const Tag &v) {
		info[p].apply(v);
		tag[p].apply(v);
	}
	void push(int p) {
		apply(2 * p, tag[p]);
		apply(2 * p + 1, tag[p]);
		tag[p] = Tag();
	}
	void modify(int p, int l, int r, int x, const Info &v) {
		if (r - l == 1) {
			info[p] = v;
			return;
		}
		int m = (l + r) / 2;
		push(p);
		if (x < m) {
			modify(2 * p, l, m, x, v);
		} else {
			modify(2 * p + 1, m, r, x, v);
		}
		pull(p);
	}
	void modify(int p, const Info &v) {
		modify(1, 0, n, p, v);
	}
	Info rangeQuery(int p, int l, int r, int x, int y) {
		if (l >= y || r <= x) {
			return Info();
		}
		if (l >= x && r <= y) {
			return info[p];
		}
		int m = (l + r) / 2;
		push(p);
		return rangeQuery(2 * p, l, m, x, y) + rangeQuery(2 * p + 1, m, r, x, y);
	}
	Info rangeQuery(int l, int r) {
		return rangeQuery(1, 0, n, l, r);
	}
	void rangeApply(int p, int l, int r, int x, int y, const Tag &v) {
		if (l >= y || r <= x) {
			return;
		}
		if (l >= x && r <= y) {
			apply(p, v);
			return;
		}
		int m = (l + r) / 2;
		push(p);
		rangeApply(2 * p, l, m, x, y, v);
		rangeApply(2 * p + 1, m, r, x, y, v);
		pull(p);
	}
	void rangeApply(int l, int r, const Tag &v) {
		return rangeApply(1, 0, n, l, r, v);
	}
	template<class F>
	int findFirst(int p, int l, int r, int x, int y, F pred) {
		if (l >= y || r <= x || !pred(info[p])) {
			return -1;
		}
		if (r - l == 1) {
			return l;
		}
		int m = (l + r) / 2;
		push(p);
		int res = findFirst(2 * p, l, m, x, y, pred);
		if (res == -1) {
			res = findFirst(2 * p + 1, m, r, x, y, pred);
		}
		return res;
	}
	template<class F>
	int findFirst(int l, int r, F pred) {
		return findFirst(1, 0, n, l, r, pred);
	}
	template<class F>
	int findLast(int p, int l, int r, int x, int y, F pred) {
		if (l >= y || r <= x || !pred(info[p])) {
			return -1;
		}
		if (r - l == 1) {
			return l;
		}
		int m = (l + r) / 2;
		push(p);
		int res = findLast(2 * p + 1, m, r, x, y, pred);
		if (res == -1) {
			res = findLast(2 * p, l, m, x, y, pred);
		}
		return res;
	}
	template<class F>
	int findLast(int l, int r, F pred) {
		return findLast(1, 0, n, l, r, pred);
	}
};

struct Tag {
	ll add = 0;

	void apply(Tag t) {
		add += t.add;
	}
};

struct Info {
	ll mn = 1E18;

	void apply(Tag t) {
		mn += t.add;
	}
};

Info operator+(Info a, Info b) {
	return {std::min(a.mn, b.mn)};
}
```

## BIT

```c++
template <typename T>
struct BIT {
	int n;
	vector<T> a;

	BIT(int n = 0) {
		init(n);
	}

	void init(int n) {
		this->n = n;
		a.assign(n, T());
	}

	void add(int x, int v) {
		while (x < n) {
			a[x] += v;
			x += x & -x;
		}
	}

	T sum(int x) {
		auto ret = T();
		while (x > 0) {
			ret += a[x];
			x -= x & -x;
		}
		return ret;
	}

	T rangeSum(int l, int r) {
		if (l == 0) return sum(r);
		else return sum(r) - sum(l - 1);
	}

};
```

## DSU

```c++
struct DSU {
	vector<int> f, siz, rnk;
	int cc;

	DSU() {}
	DSU(int n) {
		init(n);
	}

	void init(int n) {
		f.resize(n);
		iota(f.begin(), f.end(), 0);
		siz.assign(n, 1);
		rnk.assign(n, 0);
		cc = n;
	}

	int find(int x) {
		while (x != f[x]) {
			x = f[x] = f[f[x]];
		}
		return x;
	}

	bool same(int x, int y) {
		return find(x) == find(y);
	}

	int merge(int x, int y) {
		x = find(x);
		y = find(y);
		if (x == y) return -1;
		--cc;
		if (rnk[x] > rnk[y]) swap(x, y);
		siz[y] += siz[x];
		f[x] = y;
		if (rnk[x] == rnk[y]) rnk[y]++;
		return y;
	}

	int size(int x) {
		return siz[find(x)];
	}
};
```

## KD Tree

```c++
template <class T> int sgn(T x) { return (x > 0) - (x < 0); }
template<class T>
struct Point {
	typedef Point P;
	T x, y;
	explicit Point(T x=0, T y=0) : x(x), y(y) {}
	bool operator<(P p) const { return tie(x,y) < tie(p.x,p.y); }
	bool operator==(P p) const { return tie(x,y)==tie(p.x,p.y); }
	P operator+(P p) const { return P(x+p.x, y+p.y); }
	P operator-(P p) const { return P(x-p.x, y-p.y); }
	P operator*(T d) const { return P(x*d, y*d); }
	P operator/(T d) const { return P(x/d, y/d); }
	T dot(P p) const { return x*p.x + y*p.y; }
	T cross(P p) const { return x*p.y - y*p.x; }
	T cross(P a, P b) const { return (a-*this).cross(b-*this); }
	T dist2() const { return x*x + y*y; }
	double dist() const { return sqrt((double)dist2()); }
	// angle to x-axis in interval [-pi, pi]
	double angle() const { return atan2(y, x); }
	P unit() const { return *this/dist(); } // makes dist()=1
	P perp() const { return P(-y, x); } // rotates +90 degrees
	P normal() const { return perp().unit(); }
	// returns point rotated 'a' radians ccw around the origin
	P rotate(double a) const {
		return P(x*cos(a)-y*sin(a),x*sin(a)+y*cos(a)); }
	friend ostream& operator<<(ostream& os, P p) {
		return os << "(" << p.x << "," << p.y << ")"; }
};

using T = ll;
using P = Point<T>;
const T INF = numeric_limits<T>::max();

bool on_x(const P& a, const P& b) { return a.x < b.x; }
bool on_y(const P& a, const P& b) { return a.y < b.y; }

struct Node {
	P pt; // if this is a leaf, the single point in it
	T x0 = INF, x1 = -INF, y0 = INF, y1 = -INF; // bounds
	Node *first = 0, *second = 0;

	T distance(const P& p) { // min squared distance to a point
		T x = (p.x < x0 ? x0 : p.x > x1 ? x1 : p.x);
		T y = (p.y < y0 ? y0 : p.y > y1 ? y1 : p.y);
		return (P(x,y) - p).dist2();
	}

	Node(vector<P>&& vp) : pt(vp[0]) {
		for (P p : vp) {
			x0 = min(x0, p.x); x1 = max(x1, p.x);
			y0 = min(y0, p.y); y1 = max(y1, p.y);
		}
		if (vp.size() > 1) {
			// split on x if width >= height (not ideal...)
			sort(begin(vp), end(vp), x1 - x0 >= y1 - y0 ? on_x : on_y);
			// divide by taking half the array for each child (not
			// best performance with many duplicates in the middle)
			int half = (int)size(vp)/2;
			first = new Node({vp.begin(), vp.begin() + half});
			second = new Node({vp.begin() + half, vp.end()});
		}
	}
};

struct KDTree {
	Node* root;
	KDTree(const vector<P>& vp) : root(new Node({begin(vp), end(vp)})) {}

	pair<T, P> search(Node *node, const P& p) {
		if (!node->first) {
			// uncomment if we should not find the point itself:
			// if (p == node->pt) return {INF, P()};
			return make_pair((p - node->pt).dist2(), node->pt);
		}

		Node *f = node->first, *s = node->second;
		T bfirst = f->distance(p), bsec = s->distance(p);
		if (bfirst > bsec) swap(bsec, bfirst), swap(f, s);

		// search closest side first, other side if needed
		auto best = search(f, p);
		if (bsec < best.first)
			best = min(best, search(s, p));
		return best;
	}

	// find nearest point to a point, and its squared distance
	// (requires an arbitrary operator< for Point)
	pair<T, P> nearest(const P& p) {
		return search(root, p);
	}
};
```

## Segment Tree Beats

```c++
#define REP(i, n) for (int i = 0; (i) < (int)(n); ++ (i))
#define REP3(i, m, n) for (int i = (m); (i) < (int)(n); ++ (i))
#define REP_R(i, n) for (int i = (int)(n) - 1; (i) >= 0; -- (i))
#define REP3R(i, m, n) for (int i = (int)(n) - 1; (i) >= (int)(m); -- (i))
#define ALL(x) std::begin(x), std::end(x)

/**
 * @brief a segment tree beats
 */
class segment_tree_beats {
    // MEMO: values for queries (max, min, lazy_add, and lazy_update) already apply to the current node; but not for children
    typedef struct {
        int64_t max;
        int64_t max_second;
        int max_count;
        int64_t min;
        int64_t min_second;
        int min_count;
        int64_t lazy_add;
        int64_t lazy_update;
        int64_t sum;
    } value_type;

    int n;
    std::vector<value_type> a;

public:
    segment_tree_beats() = default;
    segment_tree_beats(int n_) {
        n = 1; while (n < n_) n *= 2;
        a.resize(2 * n - 1);
        tag<UPDATE>(0, 0);
    }
    template <class InputIterator>
    segment_tree_beats(InputIterator first, InputIterator last) {
        int n_ = std::distance(first, last);
        n = 1; while (n < n_) n *= 2;
        a.resize(2 * n - 1);
        REP (i, n_) {
            tag<UPDATE>(n - 1 + i, *(first + i));
        }
        REP3 (i, n_, n) {
            tag<UPDATE>(n - 1 + i, 0);
        }
        REP_R (i, n - 1) {
            update(i);
        }
    }

    void range_chmin(int l, int r, int64_t value) {  // 0-based, [l, r)
        assert (0 <= l and l <= r and r <= n);
        range_apply<CHMIN>(0, 0, n, l, r, value);
    }
    void range_chmax(int l, int r, int64_t value) {  // 0-based, [l, r)
        assert (0 <= l and l <= r and r <= n);
        range_apply<CHMAX>(0, 0, n, l, r, value);
    }
    void range_add(int l, int r, int64_t value) {  // 0-based, [l, r)
        assert (0 <= l and l <= r and r <= n);
        range_apply<ADD>(0, 0, n, l, r, value);
    }
    void range_update(int l, int r, int64_t value) {  // 0-based, [l, r)
        assert (0 <= l and l <= r and r <= n);
        range_apply<UPDATE>(0, 0, n, l, r, value);
    }

    int64_t range_min(int l, int r) {  // 0-based, [l, r)
        assert (0 <= l and l <= r and r <= n);
        return range_get<MIN>(0, 0, n, l, r);
    }
    int64_t range_max(int l, int r) {  // 0-based, [l, r)
        assert (0 <= l and l <= r and r <= n);
        return range_get<MAX>(0, 0, n, l, r);
    }
    int64_t range_sum(int l, int r) {  // 0-based, [l, r)
        assert (0 <= l and l <= r and r <= n);
        return range_get<SUM>(0, 0, n, l, r);
    }

private:
    static constexpr char CHMIN = 0;
    static constexpr char CHMAX = 1;
    static constexpr char ADD = 2;
    static constexpr char UPDATE = 3;
    static constexpr char MIN = 10;
    static constexpr char MAX = 11;
    static constexpr char SUM = 12;

    template <char TYPE>
    void range_apply(int i, int il, int ir, int l, int r, int64_t g) {
        if (ir <= l or r <= il or break_condition<TYPE>(i, g)) {
            // break
        } else if (l <= il and ir <= r and tag_condition<TYPE>(i, g)) {
            tag<TYPE>(i, g);
        } else {
            pushdown(i);
            range_apply<TYPE>(2 * i + 1, il, (il + ir) / 2, l, r, g);
            range_apply<TYPE>(2 * i + 2, (il + ir) / 2, ir, l, r, g);
            update(i);
        }
    }
    template <char TYPE>
    inline bool break_condition(int i, int64_t g) {
        switch (TYPE) {
            case CHMIN: return a[i].max <= g;
            case CHMAX: return g <= a[i].min;
            case ADD: return false;
            case UPDATE: return false;
            default: assert (false);
        }
    }
    template <char TYPE>
    inline bool tag_condition(int i, int64_t g) {
        switch (TYPE) {
            case CHMIN: return a[i].max_second < g and g < a[i].max;
            case CHMAX: return a[i].min < g and g < a[i].min_second;
            case ADD: return true;
            case UPDATE: return true;
            default: assert (false);
        }
    }
    template <char TYPE>
    inline void tag(int i, int64_t g) {
        int length = n >> (32 - __builtin_clz(i + 1) - 1);
        if (TYPE == CHMIN) {
            if (a[i].max == a[i].min or g <= a[i].min) {
                tag<UPDATE>(i, g);
                return;
            }
            if (a[i].max != INT64_MIN) {
                a[i].sum -= a[i].max * a[i].max_count;
            }
            a[i].max = g;
            a[i].min_second = std::min(a[i].min_second, g);
            if (a[i].lazy_update != INT64_MAX) {
                a[i].lazy_update = std::min(a[i].lazy_update, g);
            }
            a[i].sum += g * a[i].max_count;
        } else if (TYPE == CHMAX) {
            if (a[i].max == a[i].min or a[i].max <= g) {
                tag<UPDATE>(i, g);
                return;
            }
            if (a[i].min != INT64_MAX) {
                a[i].sum -= a[i].min * a[i].min_count;
            }
            a[i].min = g;
            a[i].max_second = std::max(a[i].max_second, g);
            if (a[i].lazy_update != INT64_MAX) {
                a[i].lazy_update = std::max(a[i].lazy_update, g);
            }
            a[i].sum += g * a[i].min_count;
        } else if (TYPE == ADD) {
            if (a[i].max != INT64_MAX) {
                a[i].max += g;
            }
            if (a[i].max_second != INT64_MIN) {
                a[i].max_second += g;
            }
            if (a[i].min != INT64_MIN) {
                a[i].min += g;
            }
            if (a[i].min_second != INT64_MAX) {
                a[i].min_second += g;
            }
            a[i].lazy_add += g;
            if (a[i].lazy_update != INT64_MAX) {
                a[i].lazy_update += g;
            }
            a[i].sum += g * length;
        } else if (TYPE == UPDATE) {
            a[i].max = g;
            a[i].max_second = INT64_MIN;
            a[i].max_count = length;
            a[i].min = g;
            a[i].min_second = INT64_MAX;
            a[i].min_count = length;
            a[i].lazy_add = 0;
            a[i].lazy_update = INT64_MAX;
            a[i].sum = g * length;
        } else {
            assert (false);
        }
    }
    void pushdown(int i) {
        int l = 2 * i + 1;
        int r = 2 * i + 2;
        // update
        if (a[i].lazy_update != INT64_MAX) {
            tag<UPDATE>(l, a[i].lazy_update);
            tag<UPDATE>(r, a[i].lazy_update);
            a[i].lazy_update = INT64_MAX;
            return;
        }
        // add
        if (a[i].lazy_add != 0) {
            tag<ADD>(l, a[i].lazy_add);
            tag<ADD>(r, a[i].lazy_add);
            a[i].lazy_add = 0;
        }
        // chmin
        if (a[i].max < a[l].max) {
            tag<CHMIN>(l, a[i].max);
        }
        if (a[i].max < a[r].max) {
            tag<CHMIN>(r, a[i].max);
        }
        // chmax
        if (a[l].min < a[i].min) {
            tag<CHMAX>(l, a[i].min);
        }
        if (a[r].min < a[i].min) {
            tag<CHMAX>(r, a[i].min);
        }
    }
    void update(int i) {
        int l = 2 * i + 1;
        int r = 2 * i + 2;
        // chmin
        std::vector<int64_t> b { a[l].max, a[l].max_second, a[r].max, a[r].max_second };
        std::sort(b.rbegin(), b.rend());
        b.erase(std::unique(ALL(b)), b.end());
        a[i].max = b[0];
        a[i].max_second = b[1];
        a[i].max_count = (b[0] == a[l].max ? a[l].max_count : 0) + (b[0] == a[r].max ? a[r].max_count : 0);
        // chmax
        std::vector<int64_t> c { a[l].min, a[l].min_second, a[r].min, a[r].min_second };
        std::sort(ALL(c));
        c.erase(std::unique(ALL(c)), c.end());
        a[i].min = c[0];
        a[i].min_second = c[1];
        a[i].min_count = (c[0] == a[l].min ? a[l].min_count : 0) + (c[0] == a[r].min ? a[r].min_count : 0);
        // add
        a[i].lazy_add = 0;
        // update
        a[i].lazy_update = INT64_MAX;
        // sum
        a[i].sum = a[l].sum + a[r].sum;
    }

    template <char TYPE>
    int64_t range_get(int i, int il, int ir, int l, int r) {
        if (ir <= l or r <= il) {
            return 0;
        } else if (l <= il and ir <= r) {
            // base
            switch (TYPE) {
                case MIN: return a[i].min;
                case MAX: return a[i].max;
                case SUM: return a[i].sum;
                default: assert (false);
            }
        } else {
            pushdown(i);
            int64_t value_l = range_get<TYPE>(2 * i + 1, il, (il + ir) / 2, l, r);
            int64_t value_r = range_get<TYPE>(2 * i + 2, (il + ir) / 2, ir, l, r);
            // mult
            switch (TYPE) {
                case MIN: return std::min(value_l, value_r);
                case MAX: return std::max(value_l, value_r);
                case SUM: return value_l + value_r;
                default: assert (false);
            }
        }
    }
};

int main() {
    int n, q; scanf("%d%d", &n, &q);

    std::vector<long long> a(n);
    for (int i = 0; i < n; i++) {
        scanf("%lld", &a[i]);
    }
    segment_tree_beats beats(ALL(a));

    for (int ph = 0; ph < q; ph++) {
        int ty, l, r; scanf("%d%d%d", &ty, &l, &r);
        if (ty == 0) {
            long long b; scanf("%lld", &b);
            beats.range_chmin(l, r, b);
        } else if (ty == 1) {
            long long b; scanf("%lld", &b);
            beats.range_chmax(l, r, b);
        } else if (ty == 2) {
            long long b; scanf("%lld", &b);
            beats.range_add(l, r, b);
        } else {
            long long sum = beats.range_sum(l, r);
            printf("%lld\n", sum);
        }
    }
    return 0;
}
```
