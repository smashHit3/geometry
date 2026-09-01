#include <bits/stdc++.h>

struct Point {
    int x, y, id;
};

struct Aabb {
    int min_x, min_y;
    int max_x, max_y;
};

struct Node {
    int l, r;
    int min_id;
    Point p;
    Aabb box;
};

struct Candidate {
    long long d;
    int id;
};

struct PointCmp{
    int axis;

    PointCmp(int _axis) : axis(_axis) {}
    bool operator()(const Point& p1, const Point& p2) const{
        int p1_v = axis ? p1.x : p1.y;
        int p2_v = axis ? p2.x : p2.y;
        return p1_v == p2_v ? p1.id < p2.id : p1_v < p2_v;
    }
};

struct CandidateCmp {
    bool operator()(const Candidate& c1, const Candidate& c2) const{
        return c1.d == c2.d ? c1.id < c2.id : c1.d > c2.d;
    }
};

using Heap = std::priority_queue<Candidate, std::vector<Candidate>, CandidateCmp>;

Aabb merge (const Aabb& b1, const Aabb& b2) {
    return Aabb{std::min (b1.min_x, b2.min_x), std::min (b1.min_y, b2.min_y), 
                std::max (b1.max_x, b2.max_x), std::max (b1.max_y, b2.max_y)};
}

long long point2AabbDistance(const Point& p, const Aabb& box) {
    int dx = std::max (std::abs (p.x - box.min_x), std::abs (p.x - box.max_x));
    int dy = std::max (std::abs (p.y - box.min_y), std::abs (p.y - box.max_y));
    return (long long)dx * dx + (long long)dy * dy;
}

long long point2PointDistance(const Point& p1, const Point& p2) {
    int dx = p1.x - p2.x;
    int dy = p1.y - p2.y;
    return (long long)dx * dx + (long long)dy * dy;
}

int build (
    // inputs
    std::vector<Point>& p, int l, int r, int axis,
    // output
    std::vector<Node>& tree) {
    if (l == r) {
        return -1; // empty set
    }

    int mid = (l + r) >> 1; // find the middle number
    std::nth_element(p.begin() + l, p.begin() + mid, p.begin() + r, PointCmp(axis));

    const int id = tree.size();
    tree.emplace_back(Node{-1, -1, p[mid].id, p[mid], Aabb{p[mid].x, p[mid].y, p[mid].x, p[mid].y}});
    const int ls = build (p, l, mid, axis ^ 1, tree);
    const int rs = build (p, mid + 1, r, axis ^ 1, tree);

    tree[id].l = ls, tree[id].r = rs;
    if (ls != -1) {
        tree[id].box = merge (tree[id].box, tree[ls].box);
        tree[id].min_id = std::min (tree[id].min_id, tree[ls].min_id);
    }
    if (rs != -1) {
        tree[id].box = merge (tree[id].box, tree[rs].box);
        tree[id].min_id = std::min (tree[id].min_id, tree[rs].min_id);
    }
    return id;
}

void query (const std::vector<Node>& tree, int cur, const Point& p, int k, Heap& best_k) {
    if (cur == -1) return ;

    Candidate c{point2PointDistance(p, tree[cur].p), tree[cur].p.id};
    
    if (best_k.size() < k) best_k.push(c);
    else if (CandidateCmp()(c, best_k.top())) best_k.pop(), best_k.push(c);

    auto query_son = [&tree, &p, k, &best_k](const Candidate& c, int son) {
        if (son == -1) return ;
        if (best_k.size() < k || CandidateCmp()(c, best_k.top())) 
            query (tree, son, p, k, best_k);
    };

    const int ls = tree[cur].l;
    const int rs = tree[cur].r;
    Candidate cl = ls == -1 ? Candidate{-1, 0} : Candidate{point2AabbDistance(p, tree[ls].box), tree[ls].min_id};
    Candidate cr = rs == -1 ? Candidate{-1, 0} : Candidate{point2AabbDistance(p, tree[rs].box), tree[rs].min_id};

    if (CandidateCmp()(cl, cr)) {
        query_son (cl, ls);
        query_son (cr, rs);
    }else {
        query_son (cr, rs);
        query_son (cl, ls);
    }
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n;
    if (!(std::cin >> n)) {
        return 0;
    }
    std::vector<Point> points;
    points.resize(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> points[i].x >> points[i].y;
        points[i].id = i + 1;
    }
    std::vector<Node> tree;
    tree.reserve(n + 1);
    const int root = build(points, 0, n, 0, tree);

    int queries;
    std::cin >> queries;
    while (queries-- > 0) {
        int x, y, k;
        std::cin >> x >> y >> k;
        Heap best;
        query(tree, root, Point{x, y, 0}, k, best);
        std::cout << best.top().id << '\n';
    }
    return 0;
}