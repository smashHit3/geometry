#include <algorithm>
#include <array>
#include <iostream>
#include <queue>
#include <vector>

using Int128 = __int128_t;

struct Point {
    long long x;
    long long y;
    int id;
};

struct Node {
    Point point;
    std::array<long long, 2> minimum;
    std::array<long long, 2> maximum;
    int minimumId;
    int left = 0;
    int right = 0;
};

struct Candidate {
    Int128 distance;
    int id;
};

// A larger distance is farther; at equal distance, a smaller original id is farther.
struct Farther {
    bool operator()(const Candidate& a, const Candidate& b) const {
        return a.distance != b.distance ? a.distance > b.distance : a.id < b.id;
    }
};

std::vector<Point> points;
std::vector<Node> tree(1);

void pull(int index) {
    Node& node = tree[index];
    node.minimum = {node.point.x, node.point.y};
    node.maximum = {node.point.x, node.point.y};
    node.minimumId = node.point.id;
    for (const int child : {node.left, node.right}) {
        if (child == 0) {
            continue;
        }
        for (int axis = 0; axis < 2; ++axis) {
            node.minimum[axis] = std::min(node.minimum[axis], tree[child].minimum[axis]);
            node.maximum[axis] = std::max(node.maximum[axis], tree[child].maximum[axis]);
        }
        node.minimumId = std::min(node.minimumId, tree[child].minimumId);
    }
}

int build(int left, int right, int depth) {
    if (left >= right) {
        return 0;
    }
    const int axis = depth % 2;
    const int middle = left + (right - left) / 2;
    std::nth_element(points.begin() + left, points.begin() + middle, points.begin() + right,
                     [axis](const Point& a, const Point& b) {
                         const long long av = axis == 0 ? a.x : a.y;
                         const long long bv = axis == 0 ? b.x : b.y;
                         if (av != bv) {
                             return av < bv;
                         }
                         return a.id < b.id;
                     });

    const int index = static_cast<int>(tree.size());
    tree.push_back({points[middle], {}, {}, points[middle].id});
    tree[index].left = build(left, middle, depth + 1);
    tree[index].right = build(middle + 1, right, depth + 1);
    pull(index);
    return index;
}

Int128 squaredDistance(long long x1, long long y1, long long x2, long long y2) {
    const Int128 dx = static_cast<Int128>(x1) - x2;
    const Int128 dy = static_cast<Int128>(y1) - y2;
    return dx * dx + dy * dy;
}

Int128 maximumDistance(const Node& node, long long x, long long y) {
    const auto absolute = [](Int128 value) { return value < 0 ? -value : value; };
    const long long farX =
        absolute(static_cast<Int128>(x) - node.minimum[0]) >
                absolute(static_cast<Int128>(x) - node.maximum[0])
            ? node.minimum[0]
            : node.maximum[0];
    const long long farY =
        absolute(static_cast<Int128>(y) - node.minimum[1]) >
                absolute(static_cast<Int128>(y) - node.maximum[1])
            ? node.minimum[1]
            : node.maximum[1];
    return squaredDistance(x, y, farX, farY);
}

using Heap = std::priority_queue<Candidate, std::vector<Candidate>, Farther>;

void query(int index, long long x, long long y, int k, Heap& best) {
    if (index == 0) {
        return;
    }

    const Node& node = tree[index];
    const Candidate current{squaredDistance(x, y, node.point.x, node.point.y), node.point.id};
    if (static_cast<int>(best.size()) < k) {
        best.push(current);
    } else if (Farther{}(current, best.top())) {
        best.pop();
        best.push(current);
    }

    const int first = node.left;
    const int second = node.right;
    const Candidate firstBound =
        first == 0 ? Candidate{-1, 0} : Candidate{maximumDistance(tree[first], x, y), tree[first].minimumId};
    const Candidate secondBound =
        second == 0 ? Candidate{-1, 0}
                    : Candidate{maximumDistance(tree[second], x, y), tree[second].minimumId};

    auto visit = [&](int child, const Candidate& bound) {
        if (child != 0 &&
            (static_cast<int>(best.size()) < k || Farther{}(bound, best.top()))) {
            query(child, x, y, k, best);
        }
    };

    if (Farther{}(firstBound, secondBound)) {
        visit(first, firstBound);
        visit(second, secondBound);
    } else {
        visit(second, secondBound);
        visit(first, firstBound);
    }
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n;
    if (!(std::cin >> n)) {
        return 0;
    }
    points.resize(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> points[i].x >> points[i].y;
        points[i].id = i + 1;
    }
    tree.reserve(n + 1);
    const int root = build(0, n, 0);

    int queries;
    std::cin >> queries;
    while (queries-- > 0) {
        long long x;
        long long y;
        int k;
        std::cin >> x >> y >> k;
        Heap best;
        query(root, x, y, k, best);
        std::cout << best.top().id << '\n';
    }
    return 0;
}
