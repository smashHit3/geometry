#include <algorithm>
#include <iostream>
#include <vector>

using Int128 = __int128_t;

struct Point {
    long long x;
    long long y;
};

struct Event {
    Int128 x;
    Int128 y;
    int query;
    int coefficient;
};

class FenwickTree {
public:
    explicit FenwickTree(int size) : values_(size + 1) {}

    void add(int index, int value) {
        for (; index < static_cast<int>(values_.size()); index += index & -index) {
            values_[index] += value;
        }
    }

    int sum(int index) const {
        int result = 0;
        for (; index > 0; index -= index & -index) {
            result += values_[index];
        }
        return result;
    }

private:
    std::vector<int> values_;
};

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n;
    int m;
    if (!(std::cin >> n >> m)) {
        return 0;
    }

    std::vector<Point> points(n);
    std::vector<long long> ys;
    ys.reserve(n);
    for (Point& point : points) {
        std::cin >> point.x >> point.y;
        ys.push_back(point.y);
    }
    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a.x != b.x ? a.x < b.x : a.y < b.y;
    });
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());

    std::vector<Event> events;
    events.reserve(4 * m);
    for (int i = 0; i < m; ++i) {
        long long x1;
        long long y1;
        long long x2;
        long long y2;
        std::cin >> x1 >> y1 >> x2 >> y2;
        // [x1, x2] x [y1, y2] = F(x2,y2)-F(x1-1,y2)-F(x2,y1-1)+F(x1-1,y1-1).
        events.push_back({x2, y2, i, 1});
        events.push_back({static_cast<Int128>(x1) - 1, y2, i, -1});
        events.push_back({x2, static_cast<Int128>(y1) - 1, i, -1});
        events.push_back(
            {static_cast<Int128>(x1) - 1, static_cast<Int128>(y1) - 1, i, 1});
    }
    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        return a.x < b.x;
    });

    FenwickTree bit(static_cast<int>(ys.size()));
    std::vector<long long> answers(m);
    int pointIndex = 0;
    for (const Event& event : events) {
        while (pointIndex < n && static_cast<Int128>(points[pointIndex].x) <= event.x) {
            const int yIndex = static_cast<int>(
                                   std::lower_bound(ys.begin(), ys.end(), points[pointIndex].y) -
                                   ys.begin()) +
                               1;
            bit.add(yIndex, 1);
            ++pointIndex;
        }
        const int yCount = static_cast<int>(
            std::upper_bound(ys.begin(), ys.end(), event.y,
                             [](const Int128 value, long long coordinate) {
                                 return value < static_cast<Int128>(coordinate);
                             }) -
            ys.begin());
        answers[event.query] += static_cast<long long>(event.coefficient) * bit.sum(yCount);
    }

    for (const long long answer : answers) {
        std::cout << answer << '\n';
    }
    return 0;
}
