#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

struct Point {
    long double x;
    long double y;
};

std::vector<Point> points;
std::vector<Point> buffer;
std::vector<Point> strip;

long double squaredDistance(const Point& a, const Point& b) {
    const long double dx = a.x - b.x;
    const long double dy = a.y - b.y;
    return dx * dx + dy * dy;
}

long double closestPair(int left, int right) {
    if (right - left <= 3) {
        long double best = std::numeric_limits<long double>::infinity();
        for (int i = left; i < right; ++i) {
            for (int j = i + 1; j < right; ++j) {
                best = std::min(best, squaredDistance(points[i], points[j]));
            }
        }
        std::sort(points.begin() + left, points.begin() + right,
                  [](const Point& a, const Point& b) { return a.y < b.y; });
        return best;
    }

    const int middle = left + (right - left) / 2;
    const long double middleX = points[middle].x;
    long double best = std::min(closestPair(left, middle), closestPair(middle, right));

    std::merge(points.begin() + left, points.begin() + middle,
               points.begin() + middle, points.begin() + right,
               buffer.begin() + left,
               [](const Point& a, const Point& b) { return a.y < b.y; });
    std::copy(buffer.begin() + left, buffer.begin() + right, points.begin() + left);

    strip.clear();
    for (int i = left; i < right; ++i) {
        const long double dx = points[i].x - middleX;
        if (dx * dx < best) {
            for (int j = static_cast<int>(strip.size()) - 1; j >= 0; --j) {
                const long double dy = points[i].y - strip[j].y;
                if (dy * dy >= best) {
                    break;
                }
                best = std::min(best, squaredDistance(points[i], strip[j]));
            }
            strip.push_back(points[i]);
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
    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a.x != b.x ? a.x < b.x : a.y < b.y;
    });
    buffer.resize(n);
    strip.reserve(n);

    const long double answer = closestPair(0, n);
    std::cout << std::fixed << std::setprecision(4) << static_cast<double>(std::sqrt(answer))
              << '\n';
    return 0;
}
