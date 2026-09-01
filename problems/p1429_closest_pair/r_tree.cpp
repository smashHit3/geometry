#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <vector>

using Distance = unsigned __int128;

// Static STR bulk-loaded R-tree with fixed fanout.
struct Point {
    long long x;
    long long y;
};

struct Box {
    long long minX;
    long long minY;
    long long maxX;
    long long maxY;
};

constexpr int kFanout = 32;

struct Node {
    Box box;
    std::array<int, kFanout> entries{};
    int count = 0;
    bool leaf = false;
};

const Distance kInfinity = ~Distance{0};

std::vector<Point> points;
std::vector<Node> nodes;

Distance squareDifference(long long a, long long b) {
    const __int128 difference = static_cast<__int128>(a) - static_cast<__int128>(b);
    const Distance magnitude =
        difference < 0 ? static_cast<Distance>(-difference) : static_cast<Distance>(difference);
    return magnitude * magnitude;
}

Distance squaredDistance(const Point& a, const Point& b) {
    return squareDifference(a.x, b.x) + squareDifference(a.y, b.y);
}

Distance squaredDistanceToBox(const Point& point, const Box& box) {
    Distance result = 0;
    if (point.x < box.minX) {
        result += squareDifference(point.x, box.minX);
    } else if (point.x > box.maxX) {
        result += squareDifference(point.x, box.maxX);
    }
    if (point.y < box.minY) {
        result += squareDifference(point.y, box.minY);
    } else if (point.y > box.maxY) {
        result += squareDifference(point.y, box.maxY);
    }
    return result;
}

Box pointBox(int pointIndex) {
    const Point& point = points[pointIndex];
    return Box{point.x, point.y, point.x, point.y};
}

Box unite(const Box& a, const Box& b) {
    return Box{std::min(a.minX, b.minX), std::min(a.minY, b.minY),
               std::max(a.maxX, b.maxX), std::max(a.maxY, b.maxY)};
}

Box entryBox(int entry, bool leaf) {
    return leaf ? pointBox(entry) : nodes[entry].box;
}

long long entryMinX(int entry, bool leaf) {
    return entryBox(entry, leaf).minX;
}

long long entryMinY(int entry, bool leaf) {
    return entryBox(entry, leaf).minY;
}

std::vector<int> packLevel(std::vector<int> entries, bool leaf) {
    std::sort(entries.begin(), entries.end(), [leaf](int lhs, int rhs) {
        const long long lhsX = entryMinX(lhs, leaf);
        const long long rhsX = entryMinX(rhs, leaf);
        if (lhsX != rhsX) {
            return lhsX < rhsX;
        }
        const long long lhsY = entryMinY(lhs, leaf);
        const long long rhsY = entryMinY(rhs, leaf);
        return lhsY != rhsY ? lhsY < rhsY : lhs < rhs;
    });

    const int groupCount = static_cast<int>((entries.size() + kFanout - 1) / kFanout);
    int sliceCount = 1;
    while (sliceCount * sliceCount < groupCount) {
        ++sliceCount;
    }
    const int nodesPerSlice = (groupCount + sliceCount - 1) / sliceCount;
    const int sliceSize = nodesPerSlice * kFanout;

    std::vector<int> ordered;
    ordered.reserve(entries.size());
    for (int begin = 0; begin < static_cast<int>(entries.size()); begin += sliceSize) {
        const int end = std::min(begin + sliceSize, static_cast<int>(entries.size()));
        std::sort(entries.begin() + begin, entries.begin() + end, [leaf](int lhs, int rhs) {
            const long long lhsY = entryMinY(lhs, leaf);
            const long long rhsY = entryMinY(rhs, leaf);
            if (lhsY != rhsY) {
                return lhsY < rhsY;
            }
            const long long lhsX = entryMinX(lhs, leaf);
            const long long rhsX = entryMinX(rhs, leaf);
            return lhsX != rhsX ? lhsX < rhsX : lhs < rhs;
        });
        ordered.insert(ordered.end(), entries.begin() + begin, entries.begin() + end);
    }

    std::vector<int> parents;
    for (int begin = 0; begin < static_cast<int>(ordered.size()); begin += kFanout) {
        const int end = std::min(begin + kFanout, static_cast<int>(ordered.size()));
        Node node;
        node.leaf = leaf;
        node.count = end - begin;
        node.box = entryBox(ordered[begin], leaf);
        for (int i = begin; i < end; ++i) {
            node.entries[i - begin] = ordered[i];
            node.box = unite(node.box, entryBox(ordered[i], leaf));
        }
        parents.push_back(static_cast<int>(nodes.size()));
        nodes.push_back(node);
    }
    return parents;
}

int build() {
    std::vector<int> level(points.size());
    std::iota(level.begin(), level.end(), 0);
    bool leaf = true;
    while (level.size() != 1) {
        level = packLevel(std::move(level), leaf);
        leaf = false;
    }
    if (nodes.empty()) {
        level = packLevel(std::move(level), true);
    }
    return level.front();
}

Distance nearestNeighbor(int root, int queryIndex, Distance best) {
    using QueueEntry = std::pair<Distance, int>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> pending;
    pending.push({0, root});

    while (!pending.empty()) {
        const auto [lowerBound, nodeIndex] = pending.top();
        pending.pop();
        if (lowerBound >= best) {
            break;
        }

        const Node& node = nodes[nodeIndex];
        if (node.leaf) {
            for (int i = 0; i < node.count; ++i) {
                const int candidate = node.entries[i];
                if (candidate != queryIndex) {
                    best = std::min(best, squaredDistance(points[queryIndex], points[candidate]));
                }
            }
            continue;
        }
        for (int i = 0; i < node.count; ++i) {
            const int child = node.entries[i];
            const Distance childBound = squaredDistanceToBox(points[queryIndex], nodes[child].box);
            if (childBound < best) {
                pending.push({childBound, child});
            }
        }
    }
    return best;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n;
    if (!(std::cin >> n)) {
        return 0;
    }
    points.resize(n);
    for (Point& point : points) {
        std::cin >> point.x >> point.y;
    }
    if (n < 2) {
        std::cout << "0.0000\n";
        return 0;
    }

    nodes.reserve((n + kFanout - 1) / kFanout * 2);
    const int root = build();
    Distance answer = kInfinity;
    for (int i = 0; i < n && answer != 0; ++i) {
        answer = nearestNeighbor(root, i, answer);
    }
    std::cout << std::fixed << std::setprecision(4)
              << std::sqrt(static_cast<long double>(answer)) << '\n';
    return 0;
}
