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
