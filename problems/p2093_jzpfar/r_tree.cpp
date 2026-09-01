#include <bits/stdc++.h>

struct Point {
    int x, y, id;
};

class rTree {
public:
    struct Aabb {
        int min_x, min_y;
        int max_x, max_y;
    };

    struct Node {
        Aabb box;
        std::vector<int> entries;
        int min_id;
        bool is_leaf;
    };

    struct Candidate {
        long long d;
        int id;
    };

    struct CandidateCmp {
        bool operator()(const Candidate& l, const Candidate& r) {
            return l.d == r.d ? l.id < r.id : l.d > r.d;
        }
    };

    using Heap = std::priority_queue<Candidate, std::vector<Candidate>, CandidateCmp>;

    rTree(const std::vector<Point>& p, const int k_fanout)
        : m_points(p), m_k_fanout(k_fanout) {
        int level_size = static_cast<int>(m_points.size());
        int total_nodes = 0;
        while (level_size > 1) {
            level_size = (level_size + m_k_fanout - 1) / m_k_fanout;
            total_nodes += level_size;
        }
        m_tree.reserve(total_nodes);
        m_root = build();
    }

    int query (int x, int y, int k) {
        Heap best;
        query (m_root, x, y, k, best);
        return best.top().id;
    }

private:
    int build () {
        std::vector<int> level (m_points.size());
        std::iota (level.begin(), level.end(), 0);
        bool is_leaf = true;
        while (level.size() != 1) {
            level = packLevel(std::move(level), is_leaf);
            is_leaf = false;
        }
        if (m_tree.empty()) {
            level = packLevel(std::move(level), true);
        }
        return level.front();
    }

    void query (int cur, int x, int y, int k, Heap& best) {
        if (m_tree[cur].is_leaf) {
            const Node& node = m_tree[cur];
            for (const int entry : node.entries) {
                const Candidate c{point2PointDistance(x, y, m_points[entry]), m_points[entry].id};
                if (best.size() < k) best.push (c);
                else if (CandidateCmp()(c, best.top())) best.pop(), best.push (c);
            }
            return ;
        }
        std::vector<std::pair<Candidate, int>> children;
        const Node& node = m_tree[cur];
        children.reserve(node.entries.size());
        for (const int entry : node.entries) {
            const Candidate c{point2AabbDistance(x, y, m_tree[entry].box), m_tree[entry].min_id};
            children.emplace_back(c, entry);
        }
        std::sort (children.begin(), children.end(), 
            [](const std::pair<Candidate, int>& l, const std::pair<Candidate, int>& r){
            return CandidateCmp()(l.first, r.first);
        });
        for (const auto child : children) {
            if (best.size() < k || CandidateCmp()(child.first, best.top())) 
                query (child.second, x, y, k, best);
        }
    }

    long long point2PointDistance(const int x, const int y, const Point& p2) {
        int dx = x - p2.x, dy = y - p2.y;
        return (long long)dx * dx + (long long)dy * dy;
    }

    long long point2AabbDistance(const int x, const int y, const Aabb& b2) {
        int dx = std::max (std::abs(x - b2.min_x), std::abs(x - b2.max_x));
        int dy = std::max (std::abs(y - b2.min_y), std::abs(y - b2.max_y));
        return (long long)dx * dx + (long long)dy * dy;
    }

    Aabb convertPointToAabb(const int i) {
        const Point& p = m_points[i];
        return Aabb{p.x, p.y, p.x, p.y};
    }

    Node genNode(const std::vector<int>::iterator l, const std::vector<int>::iterator r, const bool is_leaf) {
        Aabb aabb = is_leaf ? convertPointToAabb(*l) : m_tree[*l].box;
        int min_id = is_leaf ? m_points[*l].id : m_tree[*l].min_id;
        std::vector<int> entries;
        entries.reserve(m_k_fanout);
        entries.emplace_back(*l);
        auto i = l;
        for (++ i; i != r; i ++) {
            aabb = mergeAabb (aabb, is_leaf ? convertPointToAabb(*i) : m_tree[*i].box);
            min_id = std::min (min_id, is_leaf ? m_points[*i].id : m_tree[*i].min_id);
            entries.emplace_back(*i);
        }
        return Node{aabb, std::move(entries), min_id, is_leaf};
    }

    Aabb mergeAabb(const Aabb& b1, const Aabb& b2) {
        return Aabb{std::min (b1.min_x, b2.min_x), std::min (b1.min_y, b2.min_y), 
                    std::max (b1.max_x, b2.max_x), std::max (b1.max_y, b2.max_y)};
    }

    std::vector<int> packLevel(std::vector<int> entries, bool is_leaf) {
        std::sort (entries.begin(), entries.end(), [&](const int lhs, const int rhs){
            const auto lhs_point = is_leaf ? m_points[lhs] : Point{m_tree[lhs].box.min_x, m_tree[lhs].box.min_y, m_tree[lhs].min_id}; 
            const auto rhs_point = is_leaf ? m_points[rhs] : Point{m_tree[rhs].box.min_x, m_tree[rhs].box.min_y, m_tree[rhs].min_id};
            if (lhs_point.x != rhs_point.x) return lhs_point.x < rhs_point.x;
            if (lhs_point.y != rhs_point.y) return lhs_point.y < rhs_point.y;
            return lhs_point.id < rhs_point.id;
        });
        
        const int group_count = (entries.size() + m_k_fanout - 1) / m_k_fanout;
        int slice_count = 1;
        for (; (long long)slice_count * slice_count < group_count;) slice_count ++;
        const int slice_size = (group_count + slice_count - 1) / slice_count * m_k_fanout;

        for (int l = 0; l < entries.size(); l += slice_size) {
            int r = std::min (l + slice_size, (int)entries.size());
            std::sort (entries.begin() + l, entries.begin() + r, [&](const int lhs, const int rhs){
                const auto lhs_point = is_leaf ? m_points[lhs] : Point{m_tree[lhs].box.min_x, m_tree[lhs].box.min_y, m_tree[lhs].min_id}; 
                const auto rhs_point = is_leaf ? m_points[rhs] : Point{m_tree[rhs].box.min_x, m_tree[rhs].box.min_y, m_tree[rhs].min_id};
                if (lhs_point.y != rhs_point.y) return lhs_point.y < rhs_point.y;
                if (lhs_point.x != rhs_point.x) return lhs_point.x < rhs_point.x;
                return lhs_point.id < rhs_point.id; 
            });
        }

        std::vector<int> parents;
        parents.reserve(group_count);
        for (int l = 0; l < entries.size(); l += m_k_fanout) {
            const int r = std::min (l + m_k_fanout, (int)entries.size());
            parents.emplace_back(m_tree.size());
            m_tree.emplace_back(genNode(entries.begin() + l, entries.begin() + r, is_leaf));
        }

        return parents;
    }

    int m_root;
    int m_k_fanout;
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

    rTree r_tree(p, 32);

    int m; std::cin >> m;
    while (m --) {
        int x, y, k; std::cin >> x >> y >> k;
        std::cout << r_tree.query (x, y, k) << "\n";
    }
    return 0;
}