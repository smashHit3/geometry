#ifndef SKELETON_POLYGON_FORTUNE_VORONOI_H
#define SKELETON_POLYGON_FORTUNE_VORONOI_H

#include "common/geometry/Aabb.h"
#include "skeleton/SegmentFortuneVoronoi.h"
#include "skeleton/Types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace skeleton {

struct PolygonFeatureRef {
    enum class Kind { Vertex, Edge };

    Kind kind = Kind::Vertex;
    std::size_t inputIndex = 0;
};

struct PolygonVoronoiEdge {
    PolygonFeatureRef firstFeature;
    PolygonFeatureRef secondFeature;
    bool curved = false;
    std::vector<Point> vertices;
};

struct PolygonVoronoiDiagram {
    std::vector<PolygonVoronoiEdge> edges;
};

namespace detail {

inline bool polygonFeatureBefore(const PolygonFeatureRef& left,
                                 const PolygonFeatureRef& right) {
    if (left.kind != right.kind) {
        return left.kind == PolygonFeatureRef::Kind::Vertex;
    }
    return left.inputIndex < right.inputIndex;
}

class PolygonBoundaryIndex {
public:
    using Real = long double;

    struct Vector {
        Real x = 0.0L;
        Real y = 0.0L;
    };

    struct Bounds {
        Real minX = 0.0L;
        Real minY = 0.0L;
        Real maxX = 0.0L;
        Real maxY = 0.0L;
    };

    struct NearestFeature {
        PolygonFeatureRef feature;
        Real distance = std::numeric_limits<Real>::infinity();
    };

    static constexpr std::size_t invalidIndex() {
        return std::numeric_limits<std::size_t>::max();
    }

    PolygonBoundaryIndex(const std::vector<Point>& vertices,
                         common::geometry::Aabb<double> bounds,
                         double maximumError)
        : m_publicBounds(bounds) {
        Real minimumX = bounds.minX();
        Real minimumY = bounds.minY();
        Real maximumX = bounds.maxX();
        Real maximumY = bounds.maxY();
        for (const Point& vertex : vertices) {
            minimumX = std::min(minimumX, static_cast<Real>(vertex.x()));
            minimumY = std::min(minimumY, static_cast<Real>(vertex.y()));
            maximumX = std::max(maximumX, static_cast<Real>(vertex.x()));
            maximumY = std::max(maximumY, static_cast<Real>(vertex.y()));
        }
        m_origin = {minimumX * 0.5L + maximumX * 0.5L,
                    minimumY * 0.5L + maximumY * 0.5L};
        m_scale = std::max(
            {std::abs(static_cast<Real>(bounds.minX()) - m_origin.x),
             std::abs(static_cast<Real>(bounds.maxX()) - m_origin.x),
             std::abs(static_cast<Real>(bounds.minY()) - m_origin.y),
             std::abs(static_cast<Real>(bounds.maxY()) - m_origin.y)});
        for (const Point& vertex : vertices) {
            m_scale = std::max(
                {m_scale,
                 std::abs(static_cast<Real>(vertex.x()) - m_origin.x),
                 std::abs(static_cast<Real>(vertex.y()) - m_origin.y)});
        }
        if (!(m_scale > 0.0L) || !std::isfinite(m_scale)) {
            throw std::invalid_argument(
                "Polygon Voronoi normalization radius must be finite and positive");
        }

        m_bounds = {normalizeX(bounds.minX()), normalizeY(bounds.minY()),
                    normalizeX(bounds.maxX()), normalizeY(bounds.maxY())};
        m_vertices.reserve(vertices.size());
        for (const Point& vertex : vertices) {
            const Vector local{
                normalizeX(vertex.x()), normalizeY(vertex.y())};
            if (!finite(local)) {
                throw std::invalid_argument(
                    "Polygon Voronoi normalized coordinates must be finite");
            }
            m_vertices.push_back(local);
        }
        m_error = maximumError == 0.0
                      ? 1e-4L
                      : static_cast<Real>(maximumError) / m_scale;
        if (!std::isfinite(m_error)) {
            m_error = std::numeric_limits<Real>::max();
        }
        m_distanceTolerance =
            std::max(1e-10L, std::min(1e-3L, m_error * 8.0L));

        m_order.resize(m_vertices.size());
        std::iota(m_order.begin(), m_order.end(), 0U);
        m_nodes.reserve(m_vertices.size() * 2U);
        if (!m_order.empty()) buildNode(0U, m_order.size());
    }

    std::size_t size() const { return m_vertices.size(); }
    const Vector& vertex(std::size_t index) const {
        return m_vertices[index];
    }
    const Bounds& bounds() const { return m_bounds; }
    Real error() const { return m_error; }

    Vector localPoint(const Point& point) const {
        const Vector result{
            normalizeX(point.x()), normalizeY(point.y())};
        if (!finite(result)) {
            throw std::invalid_argument(
                "Polygon Voronoi output contains a non-finite coordinate");
        }
        return result;
    }

    Point publicPoint(Vector point) const {
        point.x = std::clamp(point.x, m_bounds.minX, m_bounds.maxX);
        point.y = std::clamp(point.y, m_bounds.minY, m_bounds.maxY);
        const Real x = std::clamp(
            point.x * m_scale + m_origin.x,
            static_cast<Real>(m_publicBounds.minX()),
            static_cast<Real>(m_publicBounds.maxX()));
        const Real y = std::clamp(
            point.y * m_scale + m_origin.y,
            static_cast<Real>(m_publicBounds.minY()),
            static_cast<Real>(m_publicBounds.maxY()));
        return {static_cast<double>(x), static_cast<double>(y)};
    }

    std::optional<Real> featureDistance(
        const PolygonFeatureRef& feature, const Vector& point) const {
        if (feature.inputIndex >= m_vertices.size()) return std::nullopt;
        if (feature.kind == PolygonFeatureRef::Kind::Vertex) {
            const std::size_t previous =
                (feature.inputIndex + m_vertices.size() - 1U) %
                m_vertices.size();
            const Vector& pointValue = m_vertices[feature.inputIndex];
            const Vector& previousValue = m_vertices[previous];
            const Vector& nextValue =
                m_vertices[(feature.inputIndex + 1U) % m_vertices.size()];
            const Real previousParameter =
                dot(subtract(point, previousValue),
                    subtract(pointValue, previousValue)) /
                squaredDistance(pointValue, previousValue);
            const Real nextParameter =
                dot(subtract(point, pointValue),
                    subtract(nextValue, pointValue)) /
                squaredDistance(nextValue, pointValue);
            const Real domainTolerance =
                64.0L * std::numeric_limits<Real>::epsilon();
            if (previousParameter < 1.0L - domainTolerance ||
                nextParameter > domainTolerance) {
                return std::nullopt;
            }
            return std::sqrt(
                squaredDistance(point, m_vertices[feature.inputIndex]));
        }

        const Vector& first = m_vertices[feature.inputIndex];
        const Vector& second =
            m_vertices[(feature.inputIndex + 1U) % m_vertices.size()];
        const Vector direction = subtract(second, first);
        const Real lengthSquared = squaredLength(direction);
        const Real parameter =
            dot(subtract(point, first), direction) / lengthSquared;
        const Real domainTolerance =
            64.0L * std::numeric_limits<Real>::epsilon();
        if (parameter < -domainTolerance ||
            parameter > 1.0L + domainTolerance) {
            return std::nullopt;
        }
        return std::abs(cross(direction, subtract(point, first))) /
               std::sqrt(lengthSquared);
    }

    Real nearestDistance(const Vector& point) const {
        Real bestSquared = std::numeric_limits<Real>::infinity();
        std::size_t bestEdge = 0;
        Real bestParameter = 0.0L;
        nearestNode(0U, point, invalidIndex(),
                    invalidIndex(), bestSquared,
                    bestEdge, bestParameter);
        return std::sqrt(bestSquared);
    }

    NearestFeature nearestFeature(
        const Vector& point,
        std::size_t excludedFirst = invalidIndex(),
        std::size_t excludedSecond = invalidIndex()) const {
        Real bestSquared = std::numeric_limits<Real>::infinity();
        std::size_t bestEdge = invalidIndex();
        Real bestParameter = 0.0L;
        nearestNode(0U, point, excludedFirst, excludedSecond,
                    bestSquared, bestEdge, bestParameter);
        NearestFeature result;
        result.distance = std::sqrt(bestSquared);
        if (bestParameter <= 0.0L) {
            result.feature = {
                PolygonFeatureRef::Kind::Vertex, bestEdge};
        } else if (bestParameter >= 1.0L) {
            result.feature = {
                PolygonFeatureRef::Kind::Vertex,
                (bestEdge + 1U) % m_vertices.size()};
        } else {
            result.feature = {
                PolygonFeatureRef::Kind::Edge, bestEdge};
        }
        return result;
    }

    bool validVoronoiPoint(const PolygonFeatureRef& first,
                           const PolygonFeatureRef& second,
                           const Vector& point) const {
        const auto firstDistance = featureDistance(first, point);
        const auto secondDistance = featureDistance(second, point);
        if (!firstDistance.has_value() ||
            !secondDistance.has_value()) {
            return false;
        }
        const Real scale =
            std::max({1.0L, *firstDistance, *secondDistance});
        if (std::abs(*firstDistance - *secondDistance) >
            m_distanceTolerance * scale) {
            return false;
        }
        const Real nearest = nearestDistance(point);
        return std::max(*firstDistance, *secondDistance) <=
               nearest + m_distanceTolerance * scale;
    }

private:
    struct Node {
        Bounds bounds;
        std::size_t begin = 0;
        std::size_t end = 0;
        std::size_t lower = PolygonBoundaryIndex::invalidIndex();
        std::size_t upper = PolygonBoundaryIndex::invalidIndex();

        bool leaf() const {
            return lower == PolygonBoundaryIndex::invalidIndex();
        }
    };

    std::vector<Vector> m_vertices;
    std::vector<std::size_t> m_order;
    std::vector<Node> m_nodes;
    common::geometry::Aabb<double> m_publicBounds;
    Vector m_origin;
    Bounds m_bounds;
    Real m_scale = 1.0L;
    Real m_error = 1e-4L;
    Real m_distanceTolerance = 1e-9L;

    Real normalizeX(double value) const {
        return (static_cast<Real>(value) - m_origin.x) / m_scale;
    }
    Real normalizeY(double value) const {
        return (static_cast<Real>(value) - m_origin.y) / m_scale;
    }

    static bool finite(const Vector& point) {
        return std::isfinite(point.x) && std::isfinite(point.y);
    }
    static Vector add(const Vector& left, const Vector& right) {
        return {left.x + right.x, left.y + right.y};
    }
    static Vector subtract(const Vector& left, const Vector& right) {
        return {left.x - right.x, left.y - right.y};
    }
    static Vector multiply(const Vector& vector, Real scalar) {
        return {vector.x * scalar, vector.y * scalar};
    }
    static Real dot(const Vector& left, const Vector& right) {
        return left.x * right.x + left.y * right.y;
    }
    static Real cross(const Vector& left, const Vector& right) {
        return left.x * right.y - left.y * right.x;
    }
    static Real squaredLength(const Vector& vector) {
        return dot(vector, vector);
    }
    static Real squaredDistance(const Vector& left, const Vector& right) {
        return squaredLength(subtract(left, right));
    }

    Bounds edgeBounds(std::size_t edge) const {
        const Vector& first = m_vertices[edge];
        const Vector& second =
            m_vertices[(edge + 1U) % m_vertices.size()];
        return {std::min(first.x, second.x),
                std::min(first.y, second.y),
                std::max(first.x, second.x),
                std::max(first.y, second.y)};
    }

    static void mergeBounds(Bounds& destination,
                            const Bounds& source) {
        destination.minX = std::min(destination.minX, source.minX);
        destination.minY = std::min(destination.minY, source.minY);
        destination.maxX = std::max(destination.maxX, source.maxX);
        destination.maxY = std::max(destination.maxY, source.maxY);
    }

    std::size_t buildNode(std::size_t begin, std::size_t end) {
        const std::size_t nodeIndex = m_nodes.size();
        m_nodes.push_back({});
        Bounds bounds = edgeBounds(m_order[begin]);
        for (std::size_t index = begin + 1U; index < end; ++index) {
            mergeBounds(bounds, edgeBounds(m_order[index]));
        }
        m_nodes[nodeIndex].bounds = bounds;
        m_nodes[nodeIndex].begin = begin;
        m_nodes[nodeIndex].end = end;
        if (end - begin <= 4U) return nodeIndex;

        const bool splitX =
            bounds.maxX - bounds.minX >= bounds.maxY - bounds.minY;
        const std::size_t middle = begin + (end - begin) / 2U;
        std::nth_element(
            m_order.begin() + static_cast<std::ptrdiff_t>(begin),
            m_order.begin() + static_cast<std::ptrdiff_t>(middle),
            m_order.begin() + static_cast<std::ptrdiff_t>(end),
            [&](std::size_t left, std::size_t right) {
                const Bounds leftBounds = edgeBounds(left);
                const Bounds rightBounds = edgeBounds(right);
                const Real leftCenter =
                    splitX ? leftBounds.minX * 0.5L +
                                 leftBounds.maxX * 0.5L
                           : leftBounds.minY * 0.5L +
                                 leftBounds.maxY * 0.5L;
                const Real rightCenter =
                    splitX ? rightBounds.minX * 0.5L +
                                 rightBounds.maxX * 0.5L
                           : rightBounds.minY * 0.5L +
                                 rightBounds.maxY * 0.5L;
                if (leftCenter != rightCenter) {
                    return leftCenter < rightCenter;
                }
                return left < right;
            });
        const std::size_t lower = buildNode(begin, middle);
        const std::size_t upper = buildNode(middle, end);
        m_nodes[nodeIndex].lower = lower;
        m_nodes[nodeIndex].upper = upper;
        return nodeIndex;
    }

    static Real distanceSquaredToBounds(const Bounds& bounds,
                                        const Vector& point) {
        const Real x =
            point.x < bounds.minX
                ? bounds.minX - point.x
                : (point.x > bounds.maxX ? point.x - bounds.maxX : 0.0L);
        const Real y =
            point.y < bounds.minY
                ? bounds.minY - point.y
                : (point.y > bounds.maxY ? point.y - bounds.maxY : 0.0L);
        return x * x + y * y;
    }

    Real distanceSquaredToEdge(std::size_t edge,
                               const Vector& point,
                               Real& parameter) const {
        const Vector& first = m_vertices[edge];
        const Vector& second =
            m_vertices[(edge + 1U) % m_vertices.size()];
        const Vector direction = subtract(second, first);
        parameter = std::clamp(
            dot(subtract(point, first), direction) /
                squaredLength(direction),
            0.0L, 1.0L);
        return squaredDistance(
            point, add(first, multiply(direction, parameter)));
    }

    void nearestNode(std::size_t nodeIndex, const Vector& point,
                     std::size_t excludedFirst,
                     std::size_t excludedSecond,
                     Real& bestSquared, std::size_t& bestEdge,
                     Real& bestParameter) const {
        const Node& node = m_nodes[nodeIndex];
        if (distanceSquaredToBounds(node.bounds, point) >= bestSquared) {
            return;
        }
        if (node.leaf()) {
            for (std::size_t index = node.begin;
                 index < node.end; ++index) {
                const std::size_t edge = m_order[index];
                if (edge == excludedFirst || edge == excludedSecond) {
                    continue;
                }
                Real parameter = 0.0L;
                const Real distance =
                    distanceSquaredToEdge(edge, point, parameter);
                if (distance < bestSquared ||
                    (distance == bestSquared && edge < bestEdge)) {
                    bestSquared = distance;
                    bestEdge = edge;
                    bestParameter = parameter;
                }
            }
            return;
        }
        const Real lowerDistance =
            distanceSquaredToBounds(m_nodes[node.lower].bounds, point);
        const Real upperDistance =
            distanceSquaredToBounds(m_nodes[node.upper].bounds, point);
        if (lowerDistance <= upperDistance) {
            nearestNode(node.lower, point, excludedFirst, excludedSecond,
                        bestSquared, bestEdge, bestParameter);
            nearestNode(node.upper, point, excludedFirst, excludedSecond,
                        bestSquared, bestEdge, bestParameter);
        } else {
            nearestNode(node.upper, point, excludedFirst, excludedSecond,
                        bestSquared, bestEdge, bestParameter);
            nearestNode(node.lower, point, excludedFirst, excludedSecond,
                        bestSquared, bestEdge, bestParameter);
        }
    }
};

class PolygonInteriorLocator {
public:
    using Vector = PolygonBoundaryIndex::Vector;

    struct Query {
        Vector point;
        std::size_t index = 0;
    };

    explicit PolygonInteriorLocator(
        const PolygonBoundaryIndex& polygon)
        : m_polygon(polygon), m_nodes(polygon.size()) {
        for (std::size_t edge = 0; edge < polygon.size(); ++edge) {
            m_nodes[edge].edge = edge;
            m_nodes[edge].priority = priority(edge);
        }
    }

    std::vector<bool> classify(const std::vector<Query>& queries) {
        enum class EventKind { Remove, Insert, Query };
        struct Event {
            long double y = 0.0L;
            EventKind kind = EventKind::Query;
            std::size_t value = 0;
            long double x = 0.0L;
        };

        std::vector<Event> events;
        events.reserve(m_polygon.size() * 2U + queries.size());
        for (std::size_t edge = 0; edge < m_polygon.size(); ++edge) {
            const Vector& first = m_polygon.vertex(edge);
            const Vector& second =
                m_polygon.vertex((edge + 1U) % m_polygon.size());
            if (first.y == second.y) continue;
            const Vector& lower = first.y < second.y ? first : second;
            const Vector& upper = first.y < second.y ? second : first;
            events.push_back(
                {upper.y, EventKind::Remove, edge, upper.x});
            events.push_back(
                {lower.y, EventKind::Insert, edge, lower.x});
        }
        for (const Query& query : queries) {
            events.push_back(
                {query.point.y, EventKind::Query,
                 query.index, query.point.x});
        }
        std::sort(events.begin(), events.end(),
                  [](const Event& left, const Event& right) {
                      if (left.y != right.y) return left.y < right.y;
                      if (left.kind != right.kind) {
                          return static_cast<int>(left.kind) <
                                 static_cast<int>(right.kind);
                      }
                      if (left.x != right.x) return left.x < right.x;
                      return left.value < right.value;
                  });

        std::vector<bool> result(queries.size(), false);
        for (const Event& event : events) {
            m_y = event.y;
            if (event.kind == EventKind::Remove) {
                erase(event.value);
            } else if (event.kind == EventKind::Insert) {
                insert(event.value);
            } else {
                result[event.value] =
                    (rankBefore(event.x) & 1U) != 0U;
            }
        }
        return result;
    }

private:
    struct Node {
        std::size_t edge = 0;
        std::size_t lower = PolygonBoundaryIndex::invalidIndex();
        std::size_t upper = PolygonBoundaryIndex::invalidIndex();
        std::size_t parent = PolygonBoundaryIndex::invalidIndex();
        std::size_t subtreeSize = 1;
        std::uint64_t priority = 0;
        bool active = false;
    };

    const PolygonBoundaryIndex& m_polygon;
    std::vector<Node> m_nodes;
    std::size_t m_root = PolygonBoundaryIndex::invalidIndex();
    long double m_y = 0.0L;

    static std::uint64_t priority(std::size_t value) {
        std::uint64_t result =
            static_cast<std::uint64_t>(value) +
            0x9e3779b97f4a7c15ULL;
        result = (result ^ (result >> 30U)) *
                 0xbf58476d1ce4e5b9ULL;
        result = (result ^ (result >> 27U)) *
                 0x94d049bb133111ebULL;
        return result ^ (result >> 31U);
    }

    std::size_t size(std::size_t node) const {
        return node == PolygonBoundaryIndex::invalidIndex()
                   ? 0U
                   : m_nodes[node].subtreeSize;
    }

    void update(std::size_t node) {
        if (node == PolygonBoundaryIndex::invalidIndex()) return;
        m_nodes[node].subtreeSize =
            1U + size(m_nodes[node].lower) +
            size(m_nodes[node].upper);
    }

    void updateAncestors(std::size_t node) {
        while (node != PolygonBoundaryIndex::invalidIndex()) {
            update(node);
            node = m_nodes[node].parent;
        }
    }

    long double edgeX(std::size_t edge) const {
        const Vector& first = m_polygon.vertex(edge);
        const Vector& second =
            m_polygon.vertex((edge + 1U) % m_polygon.size());
        const long double parameter =
            (m_y - first.y) / (second.y - first.y);
        return first.x * (1.0L - parameter) +
               second.x * parameter;
    }

    long double edgeSlope(std::size_t edge) const {
        const Vector& first = m_polygon.vertex(edge);
        const Vector& second =
            m_polygon.vertex((edge + 1U) % m_polygon.size());
        return (second.x - first.x) / (second.y - first.y);
    }

    bool edgeBefore(std::size_t left, std::size_t right) const {
        const long double leftX = edgeX(left);
        const long double rightX = edgeX(right);
        if (leftX != rightX) return leftX < rightX;
        const long double leftSlope = edgeSlope(left);
        const long double rightSlope = edgeSlope(right);
        if (leftSlope != rightSlope) return leftSlope < rightSlope;
        return left < right;
    }

    void replaceParent(std::size_t node, std::size_t replacement) {
        const std::size_t parent = m_nodes[node].parent;
        if (parent == PolygonBoundaryIndex::invalidIndex()) {
            m_root = replacement;
        } else if (m_nodes[parent].lower == node) {
            m_nodes[parent].lower = replacement;
        } else {
            m_nodes[parent].upper = replacement;
        }
        if (replacement != PolygonBoundaryIndex::invalidIndex()) {
            m_nodes[replacement].parent = parent;
        }
    }

    void rotateLeft(std::size_t node) {
        const std::size_t pivot = m_nodes[node].upper;
        const std::size_t parent = m_nodes[node].parent;
        m_nodes[node].upper = m_nodes[pivot].lower;
        if (m_nodes[node].upper !=
            PolygonBoundaryIndex::invalidIndex()) {
            m_nodes[m_nodes[node].upper].parent = node;
        }
        m_nodes[pivot].lower = node;
        m_nodes[node].parent = pivot;
        m_nodes[pivot].parent = parent;
        if (parent == PolygonBoundaryIndex::invalidIndex()) {
            m_root = pivot;
        } else if (m_nodes[parent].lower == node) {
            m_nodes[parent].lower = pivot;
        } else {
            m_nodes[parent].upper = pivot;
        }
        update(node);
        update(pivot);
    }

    void rotateRight(std::size_t node) {
        const std::size_t pivot = m_nodes[node].lower;
        const std::size_t parent = m_nodes[node].parent;
        m_nodes[node].lower = m_nodes[pivot].upper;
        if (m_nodes[node].lower !=
            PolygonBoundaryIndex::invalidIndex()) {
            m_nodes[m_nodes[node].lower].parent = node;
        }
        m_nodes[pivot].upper = node;
        m_nodes[node].parent = pivot;
        m_nodes[pivot].parent = parent;
        if (parent == PolygonBoundaryIndex::invalidIndex()) {
            m_root = pivot;
        } else if (m_nodes[parent].lower == node) {
            m_nodes[parent].lower = pivot;
        } else {
            m_nodes[parent].upper = pivot;
        }
        update(node);
        update(pivot);
    }

    void insert(std::size_t edge) {
        Node& inserted = m_nodes[edge];
        inserted.lower = PolygonBoundaryIndex::invalidIndex();
        inserted.upper = PolygonBoundaryIndex::invalidIndex();
        inserted.parent = PolygonBoundaryIndex::invalidIndex();
        inserted.subtreeSize = 1U;
        inserted.active = true;
        if (m_root == PolygonBoundaryIndex::invalidIndex()) {
            m_root = edge;
            return;
        }

        std::size_t current = m_root;
        while (true) {
            if (edgeBefore(edge, current)) {
                if (m_nodes[current].lower ==
                    PolygonBoundaryIndex::invalidIndex()) {
                    m_nodes[current].lower = edge;
                    inserted.parent = current;
                    break;
                }
                current = m_nodes[current].lower;
            } else {
                if (m_nodes[current].upper ==
                    PolygonBoundaryIndex::invalidIndex()) {
                    m_nodes[current].upper = edge;
                    inserted.parent = current;
                    break;
                }
                current = m_nodes[current].upper;
            }
        }
        updateAncestors(inserted.parent);
        while (inserted.parent !=
                   PolygonBoundaryIndex::invalidIndex() &&
               inserted.priority <
                   m_nodes[inserted.parent].priority) {
            const std::size_t parent = inserted.parent;
            if (m_nodes[parent].lower == edge) {
                rotateRight(parent);
            } else {
                rotateLeft(parent);
            }
        }
        updateAncestors(inserted.parent);
    }

    void erase(std::size_t edge) {
        Node& removed = m_nodes[edge];
        if (!removed.active) return;
        while (removed.lower !=
                   PolygonBoundaryIndex::invalidIndex() ||
               removed.upper !=
                   PolygonBoundaryIndex::invalidIndex()) {
            if (removed.upper ==
                    PolygonBoundaryIndex::invalidIndex() ||
                (removed.lower !=
                     PolygonBoundaryIndex::invalidIndex() &&
                 m_nodes[removed.lower].priority <
                     m_nodes[removed.upper].priority)) {
                rotateRight(edge);
            } else {
                rotateLeft(edge);
            }
        }
        const std::size_t parent = removed.parent;
        replaceParent(edge, PolygonBoundaryIndex::invalidIndex());
        removed.parent = PolygonBoundaryIndex::invalidIndex();
        removed.active = false;
        updateAncestors(parent);
    }

    std::size_t rankBefore(long double x) const {
        std::size_t result = 0;
        std::size_t node = m_root;
        while (node != PolygonBoundaryIndex::invalidIndex()) {
            if (edgeX(node) < x) {
                result += size(m_nodes[node].lower) + 1U;
                node = m_nodes[node].upper;
            } else {
                node = m_nodes[node].lower;
            }
        }
        return result;
    }
};

inline bool adjacentPolygonEdges(const PolygonFeatureRef& first,
                                 const PolygonFeatureRef& second,
                                 std::size_t edgeCount) {
    if (first.kind != PolygonFeatureRef::Kind::Edge ||
        second.kind != PolygonFeatureRef::Kind::Edge) {
        return false;
    }
    return (first.inputIndex + 1U) % edgeCount == second.inputIndex ||
           (second.inputIndex + 1U) % edgeCount == first.inputIndex;
}

inline PolygonVoronoiDiagram publicPolygonDiagram(
    PolygonSweepDiagram sweepDiagram,
    const std::vector<Point>& vertices,
    common::geometry::Aabb<double> bounds, double maximumError) {
    PolygonBoundaryIndex index(vertices, bounds, maximumError);
    PolygonVoronoiDiagram result;
    result.edges.reserve(sweepDiagram.edges.size());
    for (PolygonSweepEdge& edge : sweepDiagram.edges) {
        PolygonFeatureRef first{
            edge.firstFeature.vertex ? PolygonFeatureRef::Kind::Vertex
                                     : PolygonFeatureRef::Kind::Edge,
            edge.firstFeature.inputIndex};
        PolygonFeatureRef second{
            edge.secondFeature.vertex ? PolygonFeatureRef::Kind::Vertex
                                      : PolygonFeatureRef::Kind::Edge,
            edge.secondFeature.inputIndex};
        if (polygonFeatureBefore(second, first)) std::swap(first, second);
        if (adjacentPolygonEdges(first, second, vertices.size())) continue;

        std::vector<Point> current;
        auto flush = [&]() {
            if (current.size() >= 2U) {
                result.edges.push_back(
                    {first, second, edge.curved, std::move(current)});
                current.clear();
            }
        };
        for (std::size_t vertex = 1U;
             vertex < edge.vertices.size(); ++vertex) {
            const auto firstPoint =
                index.localPoint(edge.vertices[vertex - 1U]);
            const auto secondPoint =
                index.localPoint(edge.vertices[vertex]);
            const PolygonBoundaryIndex::Vector middle{
                firstPoint.x * 0.5L + secondPoint.x * 0.5L,
                firstPoint.y * 0.5L + secondPoint.y * 0.5L};
            if (index.validVoronoiPoint(first, second, middle)) {
                if (current.empty()) {
                    current.push_back(edge.vertices[vertex - 1U]);
                }
                current.push_back(edge.vertices[vertex]);
            } else {
                flush();
            }
        }
        flush();
    }

    const auto subtract = [](const PolygonBoundaryIndex::Vector& left,
                             const PolygonBoundaryIndex::Vector& right) {
        return PolygonBoundaryIndex::Vector{
            left.x - right.x, left.y - right.y};
    };
    const auto add = [](const PolygonBoundaryIndex::Vector& left,
                        const PolygonBoundaryIndex::Vector& right) {
        return PolygonBoundaryIndex::Vector{
            left.x + right.x, left.y + right.y};
    };
    const auto multiply = [](const PolygonBoundaryIndex::Vector& value,
                             long double scalar) {
        return PolygonBoundaryIndex::Vector{
            value.x * scalar, value.y * scalar};
    };
    const auto length = [](const PolygonBoundaryIndex::Vector& value) {
        return std::sqrt(value.x * value.x + value.y * value.y);
    };
    const auto dot = [](const PolygonBoundaryIndex::Vector& left,
                        const PolygonBoundaryIndex::Vector& right) {
        return left.x * right.x + left.y * right.y;
    };
    const auto cross = [](const PolygonBoundaryIndex::Vector& left,
                          const PolygonBoundaryIndex::Vector& right) {
        return left.x * right.y - left.y * right.x;
    };
    const auto rayBoundsInterval =
        [&](const PolygonBoundaryIndex::Vector& origin,
            const PolygonBoundaryIndex::Vector& direction)
        -> std::optional<std::pair<long double, long double>> {
        long double minimum = 0.0L;
        long double maximum =
            std::numeric_limits<long double>::infinity();
        const auto update =
            [&](long double coordinate, long double component,
                long double lower, long double upper) {
                if (std::abs(component) <= 1e-18L) {
                    return coordinate >= lower && coordinate <= upper;
                }
                long double first = (lower - coordinate) / component;
                long double second = (upper - coordinate) / component;
                if (first > second) std::swap(first, second);
                minimum = std::max(minimum, first);
                maximum = std::min(maximum, second);
                return minimum <= maximum;
            };
        const auto localBounds = index.bounds();
        if (!update(origin.x, direction.x,
                    localBounds.minX, localBounds.maxX) ||
            !update(origin.y, direction.y,
                    localBounds.minY, localBounds.maxY)) {
            return std::nullopt;
        }
        return std::pair<long double, long double>{minimum, maximum};
    };

    for (std::size_t vertex = 0; vertex < index.size(); ++vertex) {
        const std::size_t previous =
            (vertex + index.size() - 1U) % index.size();
        const std::size_t next = (vertex + 1U) % index.size();
        const auto center = index.vertex(vertex);
        const auto firstVector =
            subtract(index.vertex(previous), center);
        const auto secondVector =
            subtract(index.vertex(next), center);
        const long double firstLength = length(firstVector);
        const long double secondLength = length(secondVector);
        const auto firstUnit =
            multiply(firstVector, 1.0L / firstLength);
        const auto secondUnit =
            multiply(secondVector, 1.0L / secondLength);
        const std::array<PolygonBoundaryIndex::Vector, 4> candidates{
            add(firstUnit, secondUnit),
            subtract(firstUnit, secondUnit),
            multiply(add(firstUnit, secondUnit), -1.0L),
            multiply(subtract(firstUnit, secondUnit), -1.0L)};

        PolygonFeatureRef firstFeature{
            PolygonFeatureRef::Kind::Edge, previous};
        PolygonFeatureRef secondFeature{
            PolygonFeatureRef::Kind::Edge, vertex};
        if (polygonFeatureBefore(secondFeature, firstFeature)) {
            std::swap(firstFeature, secondFeature);
        }

        std::vector<PolygonBoundaryIndex::Vector> usedDirections;
        for (auto direction : candidates) {
            const long double directionLength = length(direction);
            if (directionLength <= 1e-15L) continue;
            direction = multiply(direction, 1.0L / directionLength);
            bool duplicate = false;
            for (const auto& used : usedDirections) {
                const long double differenceX = direction.x - used.x;
                const long double differenceY = direction.y - used.y;
                if (differenceX * differenceX +
                        differenceY * differenceY <=
                    1e-24L) {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate) continue;
            usedDirections.push_back(direction);

            const long double probe =
                std::max(1e-12L,
                         std::min(firstLength, secondLength) * 1e-7L);
            if (!index.validVoronoiPoint(
                    firstFeature, secondFeature,
                    add(center, multiply(direction, probe)))) {
                continue;
            }

            const auto interval = rayBoundsInterval(center, direction);
            if (!interval.has_value() ||
                interval->second - interval->first <= 1e-13L) {
                continue;
            }
            const long double minimum = interval->first;
            const long double maximum = interval->second;

            const long double intervalProbe =
                std::min(maximum,
                         minimum + std::max(
                                       1e-12L,
                                       (maximum - minimum) * 1e-9L));
            if (!index.validVoronoiPoint(
                    firstFeature, secondFeature,
                    add(center, multiply(direction, intervalProbe)))) {
                continue;
            }

            long double end = maximum;
            if (!index.validVoronoiPoint(
                    firstFeature, secondFeature,
                    add(center, multiply(direction, maximum)))) {
                long double valid = intervalProbe;
                long double invalid = maximum;
                const long double target =
                    std::max(1e-12L, index.error() * 0.25L);
                for (int iteration = 0;
                     iteration < 96 && invalid - valid > target;
                     ++iteration) {
                    const long double middle =
                        valid + (invalid - valid) * 0.5L;
                    if (index.validVoronoiPoint(
                            firstFeature, secondFeature,
                            add(center,
                                multiply(direction, middle)))) {
                        valid = middle;
                    } else {
                        invalid = middle;
                    }
                }
                end = valid;
            }
            if (end - minimum <= 1e-12L) continue;
            result.edges.push_back(
                {firstFeature, secondFeature, false,
                 {index.publicPoint(
                      add(center, multiply(direction, minimum))),
                  index.publicPoint(
                      add(center, multiply(direction, end)))}});
        }
    }

    long double twiceArea = 0.0L;
    for (std::size_t vertex = 0; vertex < index.size(); ++vertex) {
        twiceArea += cross(
            index.vertex(vertex),
            index.vertex((vertex + 1U) % index.size()));
    }
    const int winding = twiceArea > 0.0L ? 1 : -1;

    const auto ownersOverlap =
        [&](const PolygonFeatureRef& first,
            const PolygonFeatureRef& second) {
            const auto owners = [&](const PolygonFeatureRef& feature) {
                if (feature.kind == PolygonFeatureRef::Kind::Edge) {
                     return std::array<std::size_t, 2>{
                         feature.inputIndex,
                         PolygonBoundaryIndex::invalidIndex()};
                }
                return std::array<std::size_t, 2>{
                     (feature.inputIndex + index.size() - 1U) %
                         index.size(),
                     feature.inputIndex};
            };
            const auto firstOwners = owners(first);
            const auto secondOwners = owners(second);
            for (std::size_t left : firstOwners) {
                if (left == PolygonBoundaryIndex::invalidIndex()) continue;
                for (std::size_t right : secondOwners) {
                     if (left == right) return true;
                }
            }
            return false;
        };

    for (std::size_t vertex = 0; vertex < index.size(); ++vertex) {
        const std::size_t previous =
            (vertex + index.size() - 1U) % index.size();
        const std::size_t next = (vertex + 1U) % index.size();
        const auto center = index.vertex(vertex);
        const auto incoming =
            subtract(center, index.vertex(previous));
        const auto outgoing =
            subtract(index.vertex(next), center);
        if (cross(incoming, outgoing) * winding >= -1e-13L) {
            continue;
        }

        const auto previousRay =
            multiply(incoming, -1.0L / length(incoming));
        const auto nextRay =
            multiply(outgoing, 1.0L / length(outgoing));
        auto interiorDirection = add(previousRay, nextRay);
        const long double bisectorLength = length(interiorDirection);
        if (bisectorLength <= 1e-15L) continue;
        interiorDirection =
            multiply(interiorDirection, -1.0L / bisectorLength);
        const auto rayInterval =
            rayBoundsInterval(center, interiorDirection);
        if (!rayInterval.has_value() ||
            rayInterval->second <= 1e-12L) {
            continue;
        }

        std::vector<PolygonFeatureRef> candidates;
        for (int sample = 1; sample <= 32; ++sample) {
            const long double fraction =
                static_cast<long double>(sample) / 32.0L;
            const long double parameter =
                rayInterval->second *
                fraction * fraction;
            const auto nearest = index.nearestFeature(
                add(center,
                     multiply(interiorDirection, parameter)),
                previous, vertex);
            PolygonFeatureRef feature = nearest.feature;
            const PolygonFeatureRef reflexFeature{
                PolygonFeatureRef::Kind::Vertex, vertex};
            if (feature.kind != PolygonFeatureRef::Kind::Edge ||
                ownersOverlap(reflexFeature, feature)) {
                continue;
            }
            if (std::find_if(
                     candidates.begin(), candidates.end(),
                     [&](const PolygonFeatureRef& existing) {
                         return existing.kind == feature.kind &&
                                existing.inputIndex == feature.inputIndex;
                     }) == candidates.end()) {
                candidates.push_back(feature);
            }
        }

        for (PolygonFeatureRef edgeFeature : candidates) {
            PolygonFeatureRef vertexFeature{
                PolygonFeatureRef::Kind::Vertex, vertex};
            const auto lineStart =
                index.vertex(edgeFeature.inputIndex);
            const auto lineEnd = index.vertex(
                (edgeFeature.inputIndex + 1U) % index.size());
            const auto lineVector = subtract(lineEnd, lineStart);
            const long double lineLength = length(lineVector);
            const auto tangent =
                multiply(lineVector, 1.0L / lineLength);
            auto normal = PolygonBoundaryIndex::Vector{
                -tangent.y, tangent.x};
            auto relative = subtract(center, lineStart);
            long double focusNormal = dot(relative, normal);
            if (focusNormal < 0.0L) {
                normal = multiply(normal, -1.0L);
                focusNormal = -focusNormal;
            }
            if (focusNormal <= 1e-13L) continue;
            const long double focusTangent = dot(relative, tangent);
            const auto evaluate =
                [&](long double parameter) {
                     const long double normalCoordinate =
                         ((parameter - focusTangent) *
                              (parameter - focusTangent) +
                          focusNormal * focusNormal) /
                         (2.0L * focusNormal);
                     return add(
                         add(lineStart,
                             multiply(tangent, parameter)),
                         multiply(normal, normalCoordinate));
                };
            const auto validParameter =
                [&](long double parameter) {
                     const auto point = evaluate(parameter);
                     const auto localBounds = index.bounds();
                     if (point.x < localBounds.minX ||
                         point.x > localBounds.maxX ||
                         point.y < localBounds.minY ||
                         point.y > localBounds.maxY) {
                         return false;
                     }
                     return index.validVoronoiPoint(
                         vertexFeature, edgeFeature, point);
                };

            constexpr int sampleCount = 192;
            std::array<bool, sampleCount> validCells{};
            for (int sample = 0; sample < sampleCount; ++sample) {
                const long double middle =
                     lineLength *
                     (static_cast<long double>(sample) + 0.5L) /
                     static_cast<long double>(sampleCount);
                validCells[static_cast<std::size_t>(sample)] =
                     validParameter(middle);
            }
            for (int begin = 0; begin < sampleCount;) {
                while (begin < sampleCount &&
                        !validCells[static_cast<std::size_t>(begin)]) {
                     ++begin;
                }
                if (begin == sampleCount) break;
                int end = begin + 1;
                while (end < sampleCount &&
                        validCells[static_cast<std::size_t>(end)]) {
                     ++end;
                }

                long double firstParameter =
                     lineLength * static_cast<long double>(begin) /
                     static_cast<long double>(sampleCount);
                long double secondParameter =
                     lineLength * static_cast<long double>(end) /
                     static_cast<long double>(sampleCount);
                const long double target =
                     std::max(1e-12L, index.error() * 0.25L);
                if (begin > 0) {
                     long double invalid =
                         lineLength *
                         (static_cast<long double>(begin) - 0.5L) /
                         static_cast<long double>(sampleCount);
                     long double valid =
                         lineLength *
                         (static_cast<long double>(begin) + 0.5L) /
                         static_cast<long double>(sampleCount);
                     for (int iteration = 0;
                          iteration < 80 && valid - invalid > target;
                          ++iteration) {
                         const long double middle =
                             invalid + (valid - invalid) * 0.5L;
                         if (validParameter(middle)) {
                             valid = middle;
                         } else {
                             invalid = middle;
                         }
                     }
                     firstParameter = valid;
                }
                if (end < sampleCount) {
                     long double valid =
                         lineLength *
                         (static_cast<long double>(end) - 0.5L) /
                         static_cast<long double>(sampleCount);
                     long double invalid =
                         lineLength *
                         (static_cast<long double>(end) + 0.5L) /
                         static_cast<long double>(sampleCount);
                     for (int iteration = 0;
                          iteration < 80 && invalid - valid > target;
                          ++iteration) {
                         const long double middle =
                             valid + (invalid - valid) * 0.5L;
                         if (validParameter(middle)) {
                             valid = middle;
                         } else {
                             invalid = middle;
                         }
                     }
                     secondParameter = valid;
                }
                if (secondParameter - firstParameter <= 1e-12L) {
                     begin = end;
                     continue;
                }

                std::vector<PolygonBoundaryIndex::Vector> localVertices{
                     evaluate(firstParameter)};
                std::function<void(long double,
                                    const PolygonBoundaryIndex::Vector&,
                                    long double,
                                    const PolygonBoundaryIndex::Vector&,
                                    int)>
                     appendAdaptive;
                appendAdaptive =
                     [&](long double first,
                         const PolygonBoundaryIndex::Vector& firstPoint,
                         long double second,
                         const PolygonBoundaryIndex::Vector& secondPoint,
                         int depth) {
                         const long double middleParameter =
                             first + (second - first) * 0.5L;
                         const auto middlePoint =
                             evaluate(middleParameter);
                         const auto chord =
                             subtract(secondPoint, firstPoint);
                         const long double chordSquared =
                             dot(chord, chord);
                         long double chordDistance = 0.0L;
                         if (chordSquared > 0.0L) {
                             const long double projection =
                                 std::clamp(
                                     dot(subtract(middlePoint, firstPoint),
                                         chord) /
                                         chordSquared,
                                     0.0L, 1.0L);
                             chordDistance = length(
                                 subtract(
                                     middlePoint,
                                     add(firstPoint,
                                         multiply(chord,
                                                  projection))));
                         }
                         if (depth < 32 &&
                             chordDistance > index.error()) {
                             appendAdaptive(
                                 first, firstPoint, middleParameter,
                                 middlePoint, depth + 1);
                             appendAdaptive(
                                 middleParameter, middlePoint, second,
                                 secondPoint, depth + 1);
                         } else {
                             localVertices.push_back(secondPoint);
                         }
                     };
                appendAdaptive(
                     firstParameter, localVertices.front(),
                     secondParameter, evaluate(secondParameter), 0);

                PolygonFeatureRef outputFirst = vertexFeature;
                PolygonFeatureRef outputSecond = edgeFeature;
                if (polygonFeatureBefore(outputSecond, outputFirst)) {
                     std::swap(outputFirst, outputSecond);
                }
                PolygonVoronoiEdge output{
                     outputFirst, outputSecond, true, {}};
                output.vertices.reserve(localVertices.size());
                for (const auto& point : localVertices) {
                     output.vertices.push_back(index.publicPoint(point));
                }
                result.edges.push_back(std::move(output));
                begin = end;
            }
        }
    }
    return result;
}

inline PolygonVoronoiDiagram retainPolygonInterior(
    PolygonVoronoiDiagram diagram,
    const std::vector<Point>& vertices,
    common::geometry::Aabb<double> bounds) {
    PolygonBoundaryIndex index(vertices, bounds, 0.0);
    std::vector<PolygonInteriorLocator::Query> queries;
    std::vector<std::vector<std::size_t>> queryIndexes(
        diagram.edges.size());
    for (std::size_t edge = 0; edge < diagram.edges.size(); ++edge) {
        const auto& polyline = diagram.edges[edge].vertices;
        queryIndexes[edge].reserve(
            polyline.empty() ? 0U : polyline.size() - 1U);
        for (std::size_t vertex = 1U;
             vertex < polyline.size(); ++vertex) {
            const auto first = index.localPoint(polyline[vertex - 1U]);
            const auto second = index.localPoint(polyline[vertex]);
            const std::size_t query = queries.size();
            queries.push_back(
                {{first.x * 0.5L + second.x * 0.5L,
                  first.y * 0.5L + second.y * 0.5L},
                 query});
            queryIndexes[edge].push_back(query);
        }
    }

    PolygonInteriorLocator locator(index);
    const std::vector<bool> inside = locator.classify(queries);
    PolygonVoronoiDiagram result;
    for (std::size_t edge = 0; edge < diagram.edges.size(); ++edge) {
        PolygonVoronoiEdge& source = diagram.edges[edge];
        std::vector<Point> current;
        auto flush = [&]() {
            if (current.size() >= 2U) {
                result.edges.push_back(
                    {source.firstFeature, source.secondFeature,
                     source.curved, std::move(current)});
                current.clear();
            }
        };
        for (std::size_t segment = 0;
             segment < queryIndexes[edge].size(); ++segment) {
            if (inside[queryIndexes[edge][segment]]) {
                if (current.empty()) {
                    current.push_back(source.vertices[segment]);
                }
                current.push_back(source.vertices[segment + 1U]);
            } else {
                flush();
            }
        }
        flush();
    }
    return result;
}

} // namespace detail

// Builds the nearest-boundary-feature Voronoi diagram of a simple polygon.
// Vertices are canonical point features shared by their two incident edges;
// edge interiors are separate features. Decomposition edges whose feature
// owners overlap are omitted, while bisectors between adjacent edge interiors
// remain. The result is clipped to bounds and curved branches are tessellated
// to maximumError. Generic sweep construction is O(n log n); validation and
// nearest-boundary filtering use balanced sweep/BVH structures, and emitted
// polyline work is output-sensitive.
template <typename coordinate_type, typename bounds_type>
PolygonVoronoiDiagram polygonBoundaryVoronoiDiagram(
    const std::vector<common::geometry::BasePoint<coordinate_type>>& vertices,
    const common::geometry::Aabb<bounds_type>& bounds, double maximumError = 0.0) {
    std::vector<Point> doubleVertices;
    doubleVertices.reserve(vertices.size());
    for (const common::geometry::BasePoint<coordinate_type>& vertex : vertices) {
        doubleVertices.emplace_back(
            detail::checkedSegmentVoronoiCoordinate(vertex.x()),
            detail::checkedSegmentVoronoiCoordinate(vertex.y()));
    }
    const common::geometry::Aabb<double> doubleBounds{
        detail::checkedSegmentVoronoiCoordinate(bounds.minX()),
        detail::checkedSegmentVoronoiCoordinate(bounds.minY()),
        detail::checkedSegmentVoronoiCoordinate(bounds.maxX()),
        detail::checkedSegmentVoronoiCoordinate(bounds.maxY())};
    const auto sweepDiagram =
        detail::SegmentVoronoiBuilder(
            doubleVertices, doubleBounds, maximumError,
            detail::PolygonSweepTag{})
            .buildPolygonDiagram();
    return detail::publicPolygonDiagram(
        std::move(sweepDiagram), doubleVertices, doubleBounds,
        maximumError);
}

// Returns the portions of the polygon boundary-feature diagram inside the
// polygon. Point location is performed for all tessellated intervals in one
// offline order-statistic sweep, O((n + k) log n) for k output intervals,
// rather than ray-casting against every polygon edge for every sample.
template <typename coordinate_type>
PolygonVoronoiDiagram polygonMedialAxis(
    const std::vector<common::geometry::BasePoint<coordinate_type>>& vertices,
    double maximumError = 0.0) {
    if (vertices.size() < 3U) {
        throw std::invalid_argument(
            "Polygon medial axis requires at least three vertices");
    }
    std::vector<Point> doubleVertices;
    doubleVertices.reserve(vertices.size());
    for (const common::geometry::BasePoint<coordinate_type>& vertex : vertices) {
        doubleVertices.emplace_back(
            detail::checkedSegmentVoronoiCoordinate(vertex.x()),
            detail::checkedSegmentVoronoiCoordinate(vertex.y()));
    }

    double minimumX = doubleVertices.front().x();
    double minimumY = doubleVertices.front().y();
    double maximumX = minimumX;
    double maximumY = minimumY;
    for (const Point& vertex : doubleVertices) {
        minimumX = std::min(minimumX, vertex.x());
        minimumY = std::min(minimumY, vertex.y());
        maximumX = std::max(maximumX, vertex.x());
        maximumY = std::max(maximumY, vertex.y());
    }
    const common::geometry::Aabb<double> bounds{
        minimumX, minimumY, maximumX, maximumY};
    const auto sweepDiagram =
        detail::SegmentVoronoiBuilder(
            doubleVertices, bounds, maximumError,
            detail::PolygonSweepTag{})
            .buildPolygonDiagram();
    auto boundaryDiagram = detail::publicPolygonDiagram(
        std::move(sweepDiagram), doubleVertices, bounds,
        maximumError);
    return detail::retainPolygonInterior(
        std::move(boundaryDiagram), doubleVertices, bounds);
}

} // namespace skeleton

#endif // SKELETON_POLYGON_FORTUNE_VORONOI_H
