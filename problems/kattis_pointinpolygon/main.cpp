#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <unordered_set>
#include <utility>
#include <vector>

using Int128 = __int128_t;

struct Point {
    long long x;
    long long y;
};

struct Segment {
    Point left;
    Point right;
};

class TrapezoidalMap {
public:
    explicit TrapezoidalMap(const std::vector<Point>& polygon) {
        const int n = static_cast<int>(polygon.size());
        edges_.reserve(n);
        xCoordinates_.reserve(n);

        for (int i = 0; i < n; ++i) {
            const Point a = polygon[i];
            const Point b = polygon[(i + 1) % n];
            vertices_.insert(pointKey(a));
            xCoordinates_.push_back(a.x);

            if (a.x == b.x) {
                verticalEdges_[a.x].push_back({std::min(a.y, b.y), std::max(a.y, b.y)});
            } else if (a.x < b.x) {
                edges_.push_back({a, b});
            } else {
                edges_.push_back({b, a});
            }
        }

        std::sort(xCoordinates_.begin(), xCoordinates_.end());
        xCoordinates_.erase(
            std::unique(xCoordinates_.begin(), xCoordinates_.end()), xCoordinates_.end());

        for (auto& [x, intervals] : verticalEdges_) {
            (void)x;
            std::sort(intervals.begin(), intervals.end());
        }

        std::vector<std::vector<int>> starting(xCoordinates_.size());
        std::vector<std::vector<int>> ending(xCoordinates_.size());
        for (int edge = 0; edge < static_cast<int>(edges_.size()); ++edge) {
            const int start = xIndex(edges_[edge].left.x);
            const int end = xIndex(edges_[edge].right.x);
            starting[start].push_back(edge);
            ending[end].push_back(edge);
        }

        nodes_.reserve(edges_.size() * 32);
        slabRoots_.resize(xCoordinates_.size() - 1, -1);
        int root = -1;
        for (int i = 0; i < static_cast<int>(slabRoots_.size()); ++i) {
            const long long x = xCoordinates_[i];
            for (const int edge : ending[i]) {
                root = erase(root, edge, x, -1);
            }
            for (const int edge : starting[i]) {
                root = insert(root, edge, x, 1);
            }
            slabRoots_[i] = root;
        }
    }

    const char* locate(const Point query) const {
        if (vertices_.count(pointKey(query)) != 0 || onVerticalEdge(query)) {
            return "on";
        }
        if (query.x < xCoordinates_.front() || query.x > xCoordinates_.back()) {
            return "out";
        }

        int slab = static_cast<int>(
                       std::upper_bound(xCoordinates_.begin(), xCoordinates_.end(), query.x) -
                       xCoordinates_.begin()) -
                   1;
        if (slab == static_cast<int>(slabRoots_.size())) {
            --slab;
        }

        int crossingsBelow = 0;
        int node = slabRoots_[slab];
        while (node != -1) {
            const Node& current = nodes_[node];
            const int comparison = compareEdgeWithY(edges_[current.edge], query.x, query.y);
            if (comparison == 0) {
                return "on";
            }
            if (comparison < 0) {
                crossingsBelow += nodeSize(current.left) + 1;
                node = current.right;
            } else {
                node = current.left;
            }
        }
        return crossingsBelow % 2 == 1 ? "in" : "out";
    }

private:
    struct Node {
        int edge;
        int left;
        int right;
        int height;
        int size;
    };

    static std::uint64_t pointKey(const Point point) {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(point.x)) << 32) |
               static_cast<std::uint32_t>(point.y);
    }

    static Int128 yNumerator(const Segment& edge, long long x) {
        const long long dx = edge.right.x - edge.left.x;
        const long long dy = edge.right.y - edge.left.y;
        return static_cast<Int128>(edge.left.y) * dx +
               static_cast<Int128>(dy) * (x - edge.left.x);
    }

    static int compareEdgeWithY(const Segment& edge, long long x, long long y) {
        const Int128 edgeY = yNumerator(edge, x);
        const Int128 pointY = static_cast<Int128>(y) * (edge.right.x - edge.left.x);
        return (edgeY > pointY) - (edgeY < pointY);
    }

    int compareEdges(int first, int second, long long x, int side) const {
        if (first == second) {
            return 0;
        }

        const Segment& a = edges_[first];
        const Segment& b = edges_[second];
        const long long aDx = a.right.x - a.left.x;
        const long long bDx = b.right.x - b.left.x;
        const Int128 lhs = yNumerator(a, x) * bDx;
        const Int128 rhs = yNumerator(b, x) * aDx;
        if (lhs != rhs) {
            return lhs < rhs ? -1 : 1;
        }

        const Int128 aSlope = static_cast<Int128>(a.right.y - a.left.y) * bDx;
        const Int128 bSlope = static_cast<Int128>(b.right.y - b.left.y) * aDx;
        if (aSlope != bSlope) {
            const int slopeOrder = aSlope < bSlope ? -1 : 1;
            return side > 0 ? slopeOrder : -slopeOrder;
        }
        return first < second ? -1 : 1;
    }

    int xIndex(long long x) const {
        return static_cast<int>(
            std::lower_bound(xCoordinates_.begin(), xCoordinates_.end(), x) -
            xCoordinates_.begin());
    }

    int nodeHeight(int node) const {
        return node == -1 ? 0 : nodes_[node].height;
    }

    int nodeSize(int node) const {
        return node == -1 ? 0 : nodes_[node].size;
    }

    int makeNode(int edge, int left, int right) {
        nodes_.push_back(
            {edge,
             left,
             right,
             1 + std::max(nodeHeight(left), nodeHeight(right)),
             1 + nodeSize(left) + nodeSize(right)});
        return static_cast<int>(nodes_.size()) - 1;
    }

    int balance(int edge, int left, int right) {
        const int balanceFactor = nodeHeight(left) - nodeHeight(right);
        if (balanceFactor > 1) {
            const Node leftNode = nodes_[left];
            if (nodeHeight(leftNode.left) >= nodeHeight(leftNode.right)) {
                const int newRight = makeNode(edge, leftNode.right, right);
                return makeNode(leftNode.edge, leftNode.left, newRight);
            }
            const Node pivot = nodes_[leftNode.right];
            const int newLeft = makeNode(leftNode.edge, leftNode.left, pivot.left);
            const int newRight = makeNode(edge, pivot.right, right);
            return makeNode(pivot.edge, newLeft, newRight);
        }
        if (balanceFactor < -1) {
            const Node rightNode = nodes_[right];
            if (nodeHeight(rightNode.right) >= nodeHeight(rightNode.left)) {
                const int newLeft = makeNode(edge, left, rightNode.left);
                return makeNode(rightNode.edge, newLeft, rightNode.right);
            }
            const Node pivot = nodes_[rightNode.left];
            const int newLeft = makeNode(edge, left, pivot.left);
            const int newRight = makeNode(rightNode.edge, pivot.right, rightNode.right);
            return makeNode(pivot.edge, newLeft, newRight);
        }
        return makeNode(edge, left, right);
    }

    int insert(int node, int edge, long long x, int side) {
        if (node == -1) {
            return makeNode(edge, -1, -1);
        }
        const Node current = nodes_[node];
        if (compareEdges(edge, current.edge, x, side) < 0) {
            return balance(current.edge, insert(current.left, edge, x, side), current.right);
        }
        return balance(current.edge, current.left, insert(current.right, edge, x, side));
    }

    int eraseMinimum(int node) {
        const Node current = nodes_[node];
        if (current.left == -1) {
            return current.right;
        }
        return balance(current.edge, eraseMinimum(current.left), current.right);
    }

    int erase(int node, int edge, long long x, int side) {
        const Node current = nodes_[node];
        const int comparison = compareEdges(edge, current.edge, x, side);
        if (comparison < 0) {
            return balance(current.edge, erase(current.left, edge, x, side), current.right);
        }
        if (comparison > 0) {
            return balance(current.edge, current.left, erase(current.right, edge, x, side));
        }
        if (current.left == -1) {
            return current.right;
        }
        if (current.right == -1) {
            return current.left;
        }
        int successor = current.right;
        while (nodes_[successor].left != -1) {
            successor = nodes_[successor].left;
        }
        return balance(
            nodes_[successor].edge, current.left, eraseMinimum(current.right));
    }

    bool onVerticalEdge(const Point query) const {
        const auto found = verticalEdges_.find(query.x);
        if (found == verticalEdges_.end()) {
            return false;
        }

        const std::vector<std::pair<long long, long long>>& intervals = found->second;
        const auto next = std::upper_bound(
            intervals.begin(), intervals.end(), std::pair<long long, long long>{query.y, query.y});
        if (next == intervals.begin()) {
            return false;
        }
        const auto& [low, high] = *std::prev(next);
        return low <= query.y && query.y <= high;
    }

    std::vector<Segment> edges_;
    std::vector<long long> xCoordinates_;
    std::vector<int> slabRoots_;
    std::vector<Node> nodes_;
    std::map<long long, std::vector<std::pair<long long, long long>>> verticalEdges_;
    std::unordered_set<std::uint64_t> vertices_;
};

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int vertexCount;
    while (std::cin >> vertexCount && vertexCount != 0) {
        std::vector<Point> polygon(vertexCount);
        for (Point& point : polygon) {
            std::cin >> point.x >> point.y;
        }

        const TrapezoidalMap map(polygon);

        int queryCount;
        std::cin >> queryCount;
        while (queryCount-- > 0) {
            Point query{};
            std::cin >> query.x >> query.y;
            std::cout << map.locate(query) << '\n';
        }
    }
    return 0;
}
