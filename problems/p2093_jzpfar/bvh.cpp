#include <bits/stdc++.h>

struct Point {
    int x, y, id;
};

class bvhTree {
public:
    struct Candidate {
        long long d;
        int id;
    };

    struct CandidateCmp {
        bool operator () (const Candidate& c1, const Candidate& c2) {
            return c1.d == c2.d ? c1.id < c2.id : c1.d > c2.d;
        }
    };

    struct Aabb {
        int min_x = 1, min_y = 1;
        int max_x = -1, max_y = -1;
    };

    struct Node {
        int l, r;
        int pl, pr;
        int min_id;
        Aabb box;
    };
    
    using Heap = std::priority_queue<Candidate, std::vector<Candidate>, CandidateCmp>;

    bvhTree(const std::vector<Point>& p, const int leaf_size): 
        m_points(p), m_leaf_size(leaf_size) {
        int n = p.size();
        m_tree.reserve(n << 1);
        m_root = build (0, n);
    }

    int query(int x, int y, int k) {
        Heap best_k;
        query (m_root, x, y, k, best_k);
        return best_k.top().id;
    }

private:
    struct PointCmp {
        int x_axis;
        PointCmp(int _axis): x_axis(_axis) {}
        bool operator () (const Point& p1, const Point& p2) {
            const auto p1_v = x_axis ? p1.x : p1.y;
            const auto p2_v = x_axis ? p2.x : p2.y;
            return p1_v == p2_v ? p1.id < p2.id : p1_v < p2_v;
        }
    };

    int build (int l, int r) {
        int mid = (l + r) >> 1;
        const auto aabb = mergeAabbs(l, r);
        const auto min_id_iter = std::min_element(m_points.begin() + l, m_points.begin() + r, [](const Point& p1, const Point& p2) {
            return p1.id < p2.id;
        });
        const int id = m_tree.size();
        m_tree.emplace_back(Node {-1, -1, l, r, min_id_iter->id, aabb});
        if (r - l <= m_leaf_size) {
            return id;
        }
        const int x_axis = aabb.max_x - aabb.min_x >= aabb.max_y - aabb.min_y;
        std::nth_element(m_points.begin() + l, m_points.begin() + mid, m_points.begin() + r, PointCmp(x_axis));
        m_tree[id].l = build (l, mid + 1);
        m_tree[id].r = build (mid + 1, r);
        return id;
    }

    void query (int cur, int x, int y, int k, Heap& best_k) {
        if (cur == -1) return ;
        const Candidate c0{point2AabbDistance(x, y, m_tree[cur].box), m_tree[cur].min_id};
        if (best_k.size() == k && CandidateCmp()(best_k.top(), c0)) return ;
        if (isLeaf (cur)) {
            for (int i = m_tree[cur].pl; i < m_tree[cur].pr; ++ i) {
                const Candidate c{point2PointDistance(x, y, m_points[i]), m_points[i].id};
                if (best_k.size() < k) best_k.push(c);
                else if (CandidateCmp()(c, best_k.top())) best_k.pop(), best_k.push(c);
            }
            return ;
        }
        const int ls = m_tree[cur].l;
        const int rs = m_tree[cur].r;
        const Candidate cl = ls == -1 ? Candidate{-1, 0} : Candidate{point2AabbDistance(x, y, m_tree[ls].box), m_tree[ls].min_id};
        const Candidate cr = rs == -1 ? Candidate{-1, 0} : Candidate{point2AabbDistance(x, y, m_tree[rs].box), m_tree[rs].min_id};

        std::vector<int> children = {rs, ls};
        if (CandidateCmp()(cl, cr)) std::swap(children[0], children[1]);
        for (int child : children) query (child, x, y, k, best_k);
    }

    Aabb mergeAabbs(int l, int r) {
        assert (l < r);
        const auto& p = m_points[l];
        Aabb aabb{p.x, p.y, p.x, p.y};
        for (int i = l + 1; i < r; i ++) {
            aabb.min_x = std::min (aabb.min_x, m_points[i].x);
            aabb.max_x = std::max (aabb.max_x, m_points[i].x);
            aabb.min_y = std::min (aabb.min_y, m_points[i].y);
            aabb.max_y = std::max (aabb.max_y, m_points[i].y);
        }
        return aabb;
    }

    long long point2PointDistance(int x, int y, const Point& p) {
        int dx = x - p.x;
        int dy = y - p.y;
        return (long long)dx * dx + (long long)dy * dy;
    }

    long long point2AabbDistance(int x, int y, const Aabb& aabb) {
        int dx = std::max (std::abs (x - aabb.min_x), std::abs (x - aabb.max_x));
        int dy = std::max (std::abs (y - aabb.min_y), std::abs (y - aabb.max_y));
        return (long long)dx * dx + (long long)dy * dy;
    }

    bool isLeaf (int x) {
        return m_tree[x].l == -1 && m_tree[x].r == -1;
    }

    int m_root;
    int m_leaf_size;
    std::vector<Point> m_points;
    std::vector<Node> m_tree;
};

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n; std::cin >> n;
    std::vector<Point> p(n);
    for (int i = 0; i < n; i ++) {
        std::cin >> p[i].x >> p[i].y;
        p[i].id = i + 1;
    }

    bvhTree bvh_tree(p, 20);

    int m; std::cin >> m;
    while (m --) {
        int x, y, k; std::cin >> x >> y >> k;
        std::cout << bvh_tree.query (x, y, k) << "\n";
    }
    return 0;
}