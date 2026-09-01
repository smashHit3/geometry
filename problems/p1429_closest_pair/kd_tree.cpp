#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

using Distance = unsigned __int128;

// Median-split K-D tree; each node stores one point.
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

struct Node {
    Box box;
    int pointIndex = -1;
    int left = -1;
    int right = -1;
};

const Distance kInfinity = ~Distance{0};

std::vector<Point> points;
std::vector<int> pointIndices;
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

Box bounds(int begin, int end) {
    Box result{points[pointIndices[begin]].x, points[pointIndices[begin]].y,
               points[pointIndices[begin]].x, points[pointIndices[begin]].y};
    for (int i = begin + 1; i < end; ++i) {
        const Point& point = points[pointIndices[i]];
        result.minX = std::min(result.minX, point.x);
        result.minY = std::min(result.minY, point.y);
        result.maxX = std::max(result.maxX, point.x);
        result.maxY = std::max(result.maxY, point.y);
    }
    return result;
}

int build(int begin, int end) {
    const Box box = bounds(begin, end);
    const Distance width = squareDifference(box.maxX, box.minX);
    const Distance height = squareDifference(box.maxY, box.minY);
    const bool splitX = width >= height;
    const int middle = begin + (end - begin) / 2;
    std::nth_element(pointIndices.begin() + begin, pointIndices.begin() + middle,
                     pointIndices.begin() + end, [splitX](int lhs, int rhs) {
                         const Point& a = points[lhs];
                         const Point& b = points[rhs];
                         if (splitX) {
                             return a.x != b.x ? a.x < b.x
                                               : (a.y != b.y ? a.y < b.y : lhs < rhs);
                         }
                         return a.y != b.y ? a.y < b.y
                                           : (a.x != b.x ? a.x < b.x : lhs < rhs);
                     });

    const int nodeIndex = static_cast<int>(nodes.size());
    nodes.push_back(Node{box, pointIndices[middle]});
    const int left = begin < middle ? build(begin, middle) : -1;
    const int right = middle + 1 < end ? build(middle + 1, end) : -1;
    nodes[nodeIndex].left = left;
    nodes[nodeIndex].right = right;
    return nodeIndex;
}

Distance nearestNeighbor(int queryIndex, Distance best) {
    std::vector<int> pending;
    pending.push_back(0);

    while (!pending.empty()) {
        const int nodeIndex = pending.back();
        pending.pop_back();
        const Node& node = nodes[nodeIndex];
        if (squaredDistanceToBox(points[queryIndex], node.box) >= best) {
            continue;
        }
        if (node.pointIndex != queryIndex) {
            best = std::min(best, squaredDistance(points[queryIndex], points[node.pointIndex]));
        }

        if (node.left == -1 && node.right == -1) {
            continue;
        }
        const Distance leftDistance =
            node.left == -1 ? kInfinity : squaredDistanceToBox(points[queryIndex], nodes[node.left].box);
        const Distance rightDistance =
            node.right == -1 ? kInfinity : squaredDistanceToBox(points[queryIndex], nodes[node.right].box);
        const int nearChild = leftDistance <= rightDistance ? node.left : node.right;
        const int farChild = leftDistance <= rightDistance ? node.right : node.left;
        const Distance nearDistance = std::min(leftDistance, rightDistance);
        const Distance farDistance = std::max(leftDistance, rightDistance);
        if (farChild != -1 && farDistance < best) {
            pending.push_back(farChild);
        }
        if (nearChild != -1 && nearDistance < best) {
            pending.push_back(nearChild);
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

    pointIndices.resize(n);
    std::iota(pointIndices.begin(), pointIndices.end(), 0);
    nodes.reserve(n);
    build(0, n);

    Distance answer = kInfinity;
    for (int i = 0; i < n && answer != 0; ++i) {
        answer = nearestNeighbor(i, answer);
    }
    std::cout << std::fixed << std::setprecision(4)
              << std::sqrt(static_cast<long double>(answer)) << '\n';
    return 0;
}
