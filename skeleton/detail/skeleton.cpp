#include "skeleton/skeleton.h"

#include "common/BasePointUtil.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

namespace skeleton::detail {
namespace {

constexpr double kEpsilon = 1e-10;

double cross(const Point& origin, const Point& first, const Point& second) {
    return common::cross(first - origin, second - origin);
}

bool onSegment(const Point& point, const Point& first, const Point& second) {
    return std::abs(cross(first, second, point)) <= kEpsilon &&
           point.x() >= std::min(first.x(), second.x()) - kEpsilon &&
           point.x() <= std::max(first.x(), second.x()) + kEpsilon &&
           point.y() >= std::min(first.y(), second.y()) - kEpsilon &&
           point.y() <= std::max(first.y(), second.y()) + kEpsilon;
}

bool intersects(const Point& a, const Point& b, const Point& c, const Point& d) {
    const double ab_c = cross(a, b, c);
    const double ab_d = cross(a, b, d);
    const double cd_a = cross(c, d, a);
    const double cd_b = cross(c, d, b);
    return (std::abs(ab_c) <= kEpsilon && onSegment(c, a, b)) ||
           (std::abs(ab_d) <= kEpsilon && onSegment(d, a, b)) ||
           (std::abs(cd_a) <= kEpsilon && onSegment(a, c, d)) ||
           (std::abs(cd_b) <= kEpsilon && onSegment(b, c, d)) ||
           ((ab_c > 0.0) != (ab_d > 0.0) && (cd_a > 0.0) != (cd_b > 0.0));
}

bool contains(const std::vector<Point>& polygon, const Point& point) {
    bool inside = false;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const Point& a = polygon[i];
        const Point& b = polygon[(i + 1) % polygon.size()];
        if (onSegment(point, a, b)) return true;
        if ((a.y() > point.y()) != (b.y() > point.y()) &&
            a.x() + (b.x() - a.x()) * (point.y() - a.y()) / (b.y() - a.y()) > point.x()) {
            inside = !inside;
        }
    }
    return inside;
}

void validate(const std::vector<Point>& polygon) {
    if (polygon.size() < 3) throw std::invalid_argument("A polygon needs at least three vertices");
    double area = 0.0;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const Point& a = polygon[i];
        const Point& b = polygon[(i + 1) % polygon.size()];
        if (common::distance(a, b) <= kEpsilon)
            throw std::invalid_argument("A polygon cannot have zero-length edges");
        area += a.x() * b.y() - a.y() * b.x();
    }
    if (std::abs(area) <= kEpsilon) throw std::invalid_argument("A polygon must have non-zero area");
    for (std::size_t i = 0; i < polygon.size(); ++i)
        for (std::size_t j = i + 1; j < polygon.size(); ++j)
            if ((i + 1) % polygon.size() != j && (j + 1) % polygon.size() != i &&
                intersects(polygon[i], polygon[(i + 1) % polygon.size()],
                           polygon[j], polygon[(j + 1) % polygon.size()]))
                throw std::invalid_argument("The polygon must be simple");
}

struct Triangle { std::array<std::size_t, 3> vertices; Point center; };

Point circumcenter(const Point& a, const Point& b, const Point& c) {
    const double divisor = 2.0 * (a.x() * (b.y() - c.y()) + b.x() * (c.y() - a.y()) +
                                  c.x() * (a.y() - b.y()));
    if (std::abs(divisor) <= kEpsilon) throw std::domain_error("Collinear Delaunay points");
    const double aa = common::squaredLength(a), bb = common::squaredLength(b),
                 cc = common::squaredLength(c);
    return {(aa * (b.y() - c.y()) + bb * (c.y() - a.y()) + cc * (a.y() - b.y())) / divisor,
            (aa * (c.x() - b.x()) + bb * (a.x() - c.x()) + cc * (b.x() - a.x())) / divisor};
}

bool inCircle(const Triangle& triangle, const std::vector<Point>& points, const Point& p) {
    const Point& a = points[triangle.vertices[0]];
    const Point& b = points[triangle.vertices[1]];
    const Point& c = points[triangle.vertices[2]];
    const double ax = a.x() - p.x(), ay = a.y() - p.y();
    const double bx = b.x() - p.x(), by = b.y() - p.y();
    const double cx = c.x() - p.x(), cy = c.y() - p.y();
    return (ax * ax + ay * ay) * (bx * cy - by * cx) -
               (bx * bx + by * by) * (ax * cy - ay * cx) +
               (cx * cx + cy * cy) * (ax * by - ay * bx) > kEpsilon;
}

std::vector<Triangle> triangulate(std::vector<Point> points) {
    const std::size_t count = points.size();
    double min_x = points[0].x(), max_x = min_x, min_y = points[0].y(), max_y = min_y;
    for (const auto& p : points) { min_x = std::min(min_x, p.x()); max_x = std::max(max_x, p.x()); min_y = std::min(min_y, p.y()); max_y = std::max(max_y, p.y()); }
    const double scale = std::max({max_x - min_x, max_y - min_y, 1.0}) * 32.0;
    const Point center{(min_x + max_x) / 2.0, (min_y + max_y) / 2.0};
    points.push_back({center.x() - 2 * scale, center.y() - scale});
    points.push_back({center.x(), center.y() + 2 * scale});
    points.push_back({center.x() + 2 * scale, center.y() - scale});
    std::vector<Triangle> triangles{{{count, count + 2, count + 1},
                                      circumcenter(points[count], points[count + 2], points[count + 1])}};
    for (std::size_t p = 0; p < count; ++p) {
        std::map<std::pair<std::size_t, std::size_t>, int> edges;
        std::vector<Triangle> kept;
        for (const auto& t : triangles) {
            if (!inCircle(t, points, points[p])) { kept.push_back(t); continue; }
            for (std::size_t e = 0; e < 3; ++e) ++edges[std::minmax(t.vertices[e], t.vertices[(e + 1) % 3])];
        }
        triangles = std::move(kept);
        for (const auto& [edge, uses] : edges) if (uses == 1) {
            std::array<std::size_t, 3> v{edge.first, edge.second, p};
            if (cross(points[v[0]], points[v[1]], points[v[2]]) < 0) std::swap(v[0], v[1]);
            triangles.push_back({v, circumcenter(points[v[0]], points[v[1]], points[v[2]])});
        }
    }
    triangles.erase(std::remove_if(triangles.begin(), triangles.end(), [count](const Triangle& t) {
        return std::any_of(t.vertices.begin(), t.vertices.end(), [count](std::size_t v) { return v >= count; });
    }), triangles.end());
    return triangles;
}

} // namespace

std::vector<MedialAxisEdge> medialAxis(const std::vector<Point>& polygon,
                                       const MedialAxisOptions& options) {
    if (options.max_boundary_segment_length < 0.0) throw std::invalid_argument("The boundary sampling length cannot be negative");
    detail::validate(polygon);
    double min_x = polygon[0].x(), max_x = min_x, min_y = polygon[0].y(), max_y = min_y;
    for (const auto& p : polygon) { min_x = std::min(min_x, p.x()); max_x = std::max(max_x, p.x()); min_y = std::min(min_y, p.y()); max_y = std::max(max_y, p.y()); }
    const double spacing = options.max_boundary_segment_length > 0 ? options.max_boundary_segment_length : std::hypot(max_x - min_x, max_y - min_y) / 64.0;
    std::vector<Point> samples; std::vector<std::size_t> sources;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const Point& a = polygon[i]; const Point& b = polygon[(i + 1) % polygon.size()];
        const auto n = static_cast<std::size_t>(std::ceil(common::distance(a, b) / spacing));
        for (std::size_t j = 0; j < n; ++j) { samples.push_back(a + (b - a) * (static_cast<double>(j) / n)); sources.push_back(i); }
    }
    const auto triangles = detail::triangulate(samples);
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> adjacent;
    std::vector<MedialAxisEdge> result;
    for (std::size_t i = 0; i < triangles.size(); ++i) for (std::size_t e = 0; e < 3; ++e) {
        const auto edge = std::minmax(triangles[i].vertices[e], triangles[i].vertices[(e + 1) % 3]);
        const auto [it, new_edge] = adjacent.emplace(edge, i);
        if (!new_edge && sources[edge.first] != sources[edge.second]) {
            const Point middle = (triangles[i].center + triangles[it->second].center) / 2.0;
            if (detail::contains(polygon, middle))
                result.push_back({triangles[it->second].center, triangles[i].center});
        }
    }
    return result;
}

} // namespace skeleton::detail
