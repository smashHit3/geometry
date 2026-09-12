#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace {

using Integer = long long;
using WideInteger = __int128_t;
using Real = long double;

constexpr Real EPSILON = 1e-12L;

struct IntegerPoint {
    Integer x;
    Integer y;
};

struct Point {
    Real x;
    Real y;
};

Point operator+(const Point& left, const Point& right)
{
    return {left.x + right.x, left.y + right.y};
}

Point operator-(const Point& left, const Point& right)
{
    return {left.x - right.x, left.y - right.y};
}

Point operator*(const Point& point, Real scalar)
{
    return {point.x * scalar, point.y * scalar};
}

Real cross(const Point& left, const Point& right)
{
    return left.x * right.y - left.y * right.x;
}

Real cross(const Point& origin, const Point& first, const Point& second)
{
    return cross(first - origin, second - origin);
}

Real dot(const Point& left, const Point& right)
{
    return left.x * right.x + left.y * right.y;
}

Real squaredDistance(const Point& first, const Point& second)
{
    const Point difference = first - second;
    return dot(difference, difference);
}

WideInteger cross(const IntegerPoint& origin, const IntegerPoint& first,
                  const IntegerPoint& second)
{
    return static_cast<WideInteger>(first.x - origin.x) *
               (second.y - origin.y) -
           static_cast<WideInteger>(first.y - origin.y) *
               (second.x - origin.x);
}

WideInteger incircleDeterminant(const IntegerPoint& first,
                                const IntegerPoint& second,
                                const IntegerPoint& third,
                                const IntegerPoint& point)
{
    const WideInteger ax = first.x - point.x;
    const WideInteger ay = first.y - point.y;
    const WideInteger bx = second.x - point.x;
    const WideInteger by = second.y - point.y;
    const WideInteger cx = third.x - point.x;
    const WideInteger cy = third.y - point.y;

    return
        (ax * ax + ay * ay) * (bx * cy - by * cx) +
        (bx * bx + by * by) * (cx * ay - cy * ax) +
        (cx * cx + cy * cy) * (ax * by - ay * bx);
}

class DelaunayTriangulation {
public:
    explicit DelaunayTriangulation(
        const std::vector<IntegerPoint>& points)
        : points_(points)
    {
    }

    std::vector<std::vector<int>> build()
    {
        std::vector<int> order(points_.size());
        for (int index = 0; index < static_cast<int>(points_.size());
             ++index) {
            order[index] = index;
        }
        std::sort(order.begin(), order.end(),
                  [&](int left, int right) {
                      if (points_[left].x != points_[right].x) {
                          return points_[left].x < points_[right].x;
                      }
                      return points_[left].y < points_[right].y;
                  });

        buildRange(order, 0, static_cast<int>(order.size()) - 1);

        std::vector<std::set<int>> neighbor_sets(points_.size());
        for (QuadEdge* edge : primal_edges_) {
            if (edge->removed) {
                continue;
            }
            const int first = edge->origin;
            const int second = edge->destination();
            neighbor_sets[first].insert(second);
            neighbor_sets[second].insert(first);
        }

        std::vector<std::vector<int>> neighbors(points_.size());
        for (std::size_t index = 0; index < points_.size(); ++index) {
            neighbors[index].assign(neighbor_sets[index].begin(),
                                    neighbor_sets[index].end());
        }
        return neighbors;
    }

private:
    struct QuadEdge {
        int origin = -1;
        QuadEdge* rotated = nullptr;
        QuadEdge* next = nullptr;
        bool removed = false;

        QuadEdge* reversed() const { return rotated->rotated; }
        int destination() const { return reversed()->origin; }
        QuadEdge* leftNext() const
        {
            return rotated->reversed()->next->rotated;
        }
        QuadEdge* originPrevious() const
        {
            return rotated->next->rotated;
        }
    };

    const std::vector<IntegerPoint>& points_;
    std::vector<std::unique_ptr<QuadEdge>> storage_;
    std::vector<QuadEdge*> primal_edges_;

    QuadEdge* allocateEdge()
    {
        storage_.push_back(std::make_unique<QuadEdge>());
        return storage_.back().get();
    }

    QuadEdge* makeEdge(int from, int to)
    {
        QuadEdge* forward = allocateEdge();
        QuadEdge* backward = allocateEdge();
        QuadEdge* dual_forward = allocateEdge();
        QuadEdge* dual_backward = allocateEdge();

        forward->origin = from;
        backward->origin = to;
        forward->rotated = dual_forward;
        backward->rotated = dual_backward;
        dual_forward->rotated = backward;
        dual_backward->rotated = forward;
        forward->next = forward;
        backward->next = backward;
        dual_forward->next = dual_backward;
        dual_backward->next = dual_forward;
        primal_edges_.push_back(forward);
        return forward;
    }

    static void splice(QuadEdge* first, QuadEdge* second)
    {
        std::swap(first->next->rotated->next,
                  second->next->rotated->next);
        std::swap(first->next, second->next);
    }

    void removeEdge(QuadEdge* edge)
    {
        splice(edge, edge->originPrevious());
        splice(edge->reversed(), edge->reversed()->originPrevious());
        edge->removed = true;
        edge->reversed()->removed = true;
    }

    QuadEdge* connect(QuadEdge* first, QuadEdge* second)
    {
        QuadEdge* edge =
            makeEdge(first->destination(), second->origin);
        splice(edge, first->leftNext());
        splice(edge->reversed(), second);
        return edge;
    }

    bool leftOf(int point, QuadEdge* edge) const
    {
        return cross(points_[point], points_[edge->origin],
                     points_[edge->destination()]) > 0;
    }

    bool rightOf(int point, QuadEdge* edge) const
    {
        return cross(points_[point], points_[edge->origin],
                     points_[edge->destination()]) < 0;
    }

    bool inCircle(int first, int second, int third, int fourth) const
    {
        return incircleDeterminant(
                   points_[first], points_[second],
                   points_[third], points_[fourth]) > 0;
    }

    std::pair<QuadEdge*, QuadEdge*> buildRange(
        const std::vector<int>& order, int left, int right)
    {
        if (right - left == 1) {
            QuadEdge* edge = makeEdge(order[left], order[right]);
            return {edge, edge->reversed()};
        }

        if (right - left == 2) {
            QuadEdge* first =
                makeEdge(order[left], order[left + 1]);
            QuadEdge* second =
                makeEdge(order[left + 1], order[right]);
            splice(first->reversed(), second);

            const WideInteger orientation =
                cross(points_[order[left]], points_[order[left + 1]],
                      points_[order[right]]);
            if (orientation == 0) {
                return {first, second->reversed()};
            }

            QuadEdge* closing = connect(second, first);
            if (orientation > 0) {
                return {first, second->reversed()};
            }
            return {closing->reversed(), closing};
        }

        const int middle = (left + right) / 2;
        auto [left_outer, left_inner] =
            buildRange(order, left, middle);
        auto [right_inner, right_outer] =
            buildRange(order, middle + 1, right);

        while (true) {
            if (leftOf(right_inner->origin, left_inner)) {
                left_inner = left_inner->leftNext();
            } else if (rightOf(left_inner->origin, right_inner)) {
                right_inner = right_inner->reversed()->next;
            } else {
                break;
            }
        }

        QuadEdge* base = connect(right_inner->reversed(), left_inner);
        if (left_inner->origin == left_outer->origin) {
            left_outer = base->reversed();
        }
        if (right_inner->origin == right_outer->origin) {
            right_outer = base;
        }

        const auto valid = [&](QuadEdge* edge) {
            return rightOf(edge->destination(), base);
        };

        while (true) {
            QuadEdge* left_candidate = base->reversed()->next;
            if (valid(left_candidate)) {
                while (inCircle(
                    base->destination(), base->origin,
                    left_candidate->destination(),
                    left_candidate->next->destination())) {
                    QuadEdge* next = left_candidate->next;
                    removeEdge(left_candidate);
                    left_candidate = next;
                }
            }

            QuadEdge* right_candidate = base->originPrevious();
            if (valid(right_candidate)) {
                while (inCircle(
                    base->destination(), base->origin,
                    right_candidate->destination(),
                    right_candidate->originPrevious()->destination())) {
                    QuadEdge* previous =
                        right_candidate->originPrevious();
                    removeEdge(right_candidate);
                    right_candidate = previous;
                }
            }

            const bool left_valid = valid(left_candidate);
            const bool right_valid = valid(right_candidate);
            if (!left_valid && !right_valid) {
                break;
            }
            if (!left_valid ||
                (right_valid &&
                 inCircle(left_candidate->destination(),
                          left_candidate->origin,
                          right_candidate->origin,
                          right_candidate->destination()))) {
                base = connect(right_candidate, base->reversed());
            } else {
                base =
                    connect(base->reversed(),
                            left_candidate->reversed());
            }
        }
        return {left_outer, right_outer};
    }
};

std::vector<std::vector<int>> delaunayNeighbors(
    const std::vector<IntegerPoint>& points)
{
    return DelaunayTriangulation(points).build();
}

Real halfPlaneValue(const Point& point, const Point& normal, Real offset)
{
    return dot(point, normal) - offset;
}

std::vector<Point> clipPolygon(const std::vector<Point>& polygon,
                               const Point& normal, Real offset)
{
    std::vector<Point> clipped;
    if (polygon.empty()) {
        return clipped;
    }

    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const Point& current = polygon[index];
        const Point& next = polygon[(index + 1) % polygon.size()];
        const Real current_value =
            halfPlaneValue(current, normal, offset);
        const Real next_value = halfPlaneValue(next, normal, offset);
        const bool current_inside = current_value <= EPSILON;
        const bool next_inside = next_value <= EPSILON;

        if (current_inside) {
            clipped.push_back(current);
        }
        if (current_inside == next_inside) {
            continue;
        }

        const Real ratio = current_value / (current_value - next_value);
        clipped.push_back(current + (next - current) * ratio);
    }
    return clipped;
}

std::vector<Point> buildVoronoiCell(
    int site_index, const std::vector<Point>& sites,
    const std::vector<int>& neighbors, Real minimum_x, Real minimum_y,
    Real maximum_x, Real maximum_y)
{
    std::vector<Point> cell{
        {minimum_x, minimum_y},
        {maximum_x, minimum_y},
        {maximum_x, maximum_y},
        {minimum_x, maximum_y},
    };

    const Point& site = sites[site_index];
    for (const int neighbor_index : neighbors) {
        const Point& neighbor = sites[neighbor_index];
        const Point normal = neighbor - site;
        const Real offset =
            (dot(neighbor, neighbor) - dot(site, site)) / 2;
        cell = clipPolygon(cell, normal, offset);
    }
    return cell;
}

bool pointOnSegment(const Point& point, const Point& first,
                    const Point& second)
{
    if (std::abs(cross(first, second, point)) > EPSILON) {
        return false;
    }
    return dot(point - first, point - second) <= EPSILON;
}

bool pointInPolygon(const Point& point, const std::vector<Point>& polygon)
{
    bool inside = false;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const Point& first = polygon[index];
        const Point& second = polygon[(index + 1) % polygon.size()];
        if (pointOnSegment(point, first, second)) {
            return true;
        }

        const bool crosses =
            (first.y > point.y) != (second.y > point.y);
        if (crosses) {
            const Real intersection_x =
                first.x + (second.x - first.x) *
                              (point.y - first.y) /
                              (second.y - first.y);
            if (intersection_x > point.x) {
                inside = !inside;
            }
        }
    }
    return inside;
}

bool segmentIntersection(const Point& first_start, const Point& first_end,
                         const Point& second_start, const Point& second_end,
                         Point& intersection)
{
    const Point first_direction = first_end - first_start;
    const Point second_direction = second_end - second_start;
    const Real denominator = cross(first_direction, second_direction);
    if (std::abs(denominator) <= EPSILON) {
        return false;
    }

    const Point difference = second_start - first_start;
    const Real first_parameter =
        cross(difference, second_direction) / denominator;
    const Real second_parameter =
        cross(difference, first_direction) / denominator;
    if (first_parameter < -EPSILON ||
        first_parameter > 1 + EPSILON ||
        second_parameter < -EPSILON ||
        second_parameter > 1 + EPSILON) {
        return false;
    }

    intersection = first_start + first_direction * first_parameter;
    return true;
}

} // namespace

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int vertex_count;
    std::cin >> vertex_count;

    std::vector<IntegerPoint> integer_vertices(vertex_count);
    std::vector<Point> polygon(vertex_count);
    Real minimum_x = 0;
    Real minimum_y = 0;
    Real maximum_x = 0;
    Real maximum_y = 0;
    for (int index = 0; index < vertex_count; ++index) {
        std::cin >> integer_vertices[index].x >> integer_vertices[index].y;
        polygon[index] = {
            static_cast<Real>(integer_vertices[index].x),
            static_cast<Real>(integer_vertices[index].y)};
        if (index == 0) {
            minimum_x = maximum_x = polygon[index].x;
            minimum_y = maximum_y = polygon[index].y;
        } else {
            minimum_x = std::min(minimum_x, polygon[index].x);
            minimum_y = std::min(minimum_y, polygon[index].y);
            maximum_x = std::max(maximum_x, polygon[index].x);
            maximum_y = std::max(maximum_y, polygon[index].y);
        }
    }

    const std::vector<std::vector<int>> neighbors =
        delaunayNeighbors(integer_vertices);
    const Real padding =
        std::max(maximum_x - minimum_x, maximum_y - minimum_y) + 1;

    Real answer_squared = 0;
    for (int site_index = 0; site_index < vertex_count; ++site_index) {
        const std::vector<Point> cell = buildVoronoiCell(
            site_index, polygon, neighbors[site_index],
            minimum_x - padding, minimum_y - padding,
            maximum_x + padding, maximum_y + padding);

        for (const Point& vertex : cell) {
            if (pointInPolygon(vertex, polygon)) {
                answer_squared =
                    std::max(answer_squared,
                             squaredDistance(vertex, polygon[site_index]));
            }
        }

        for (std::size_t cell_index = 0; cell_index < cell.size();
             ++cell_index) {
            const Point& cell_start = cell[cell_index];
            const Point& cell_end =
                cell[(cell_index + 1) % cell.size()];
            for (int polygon_index = 0; polygon_index < vertex_count;
                 ++polygon_index) {
                Point intersection;
                if (segmentIntersection(
                        cell_start, cell_end, polygon[polygon_index],
                        polygon[(polygon_index + 1) % vertex_count],
                        intersection)) {
                    answer_squared = std::max(
                        answer_squared,
                        squaredDistance(intersection,
                                        polygon[site_index]));
                }
            }
        }
    }

    std::cout << std::fixed << std::setprecision(10)
              << std::sqrt(answer_squared) << '\n';
    return 0;
}