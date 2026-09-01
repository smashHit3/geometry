#ifndef SKELETON_SEGMENT_FORTUNE_VORONOI_H
#define SKELETON_SEGMENT_FORTUNE_VORONOI_H

#include "common/geometry/Aabb.h"
#include "common/geometryutil/BasePointUtil.h"
#include "skeleton/Types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <list>
#include <memory>
#include <optional>
#include <queue>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace skeleton {

struct SegmentVoronoiEdge {
    std::size_t firstSite = 0;
    std::size_t secondSite = 0;
    bool curved = false;
    std::vector<Point> vertices;
};

struct SegmentVoronoiCell {
    Segment site;
    std::vector<std::size_t> edges;
};

struct SegmentVoronoiDiagram {
    std::vector<SegmentVoronoiEdge> edges;
    std::vector<SegmentVoronoiCell> cells;
};

namespace detail {

template <typename coordinate_type>
double checkedSegmentVoronoiCoordinate(coordinate_type value) {
    const double converted = static_cast<double>(value);
    if (!std::isfinite(converted)) {
        throw std::invalid_argument(
            "Segment Voronoi coordinates must be exactly representable as double");
    }

    using Value = std::remove_cv_t<coordinate_type>;
    if constexpr (std::is_integral_v<Value>) {
        // A floating-to-integral conversion is only defined in the target
        // range.  Check that range in double first, then perform the
        // type-preserving round trip.  This does not rely on long double
        // having more precision than double.
        const double limit =
            std::ldexp(1.0, std::numeric_limits<Value>::digits);
        if constexpr (std::is_signed_v<Value>) {
            if (converted < -limit || converted >= limit) {
                throw std::invalid_argument(
                    "Segment Voronoi coordinates must be exactly representable as double");
            }
        } else if (converted < 0.0 || converted >= limit) {
            throw std::invalid_argument(
                "Segment Voronoi coordinates must be exactly representable as double");
        }
        if (static_cast<Value>(converted) != value) {
            throw std::invalid_argument(
                "Segment Voronoi coordinates must be exactly representable as double");
        }
    } else {
        if (!std::isfinite(value) ||
            static_cast<Value>(converted) != value) {
            throw std::invalid_argument(
                "Segment Voronoi coordinates must be exactly representable as double");
        }
    }
    return converted;
}

struct PolygonSweepTag {};

struct PolygonSweepFeatureRef {
    bool vertex = false;
    std::size_t inputIndex = 0;
};

struct PolygonSweepEdge {
    PolygonSweepFeatureRef firstFeature;
    PolygonSweepFeatureRef secondFeature;
    bool curved = false;
    std::vector<Point> vertices;
};

struct PolygonSweepDiagram {
    std::vector<PolygonSweepEdge> edges;
};

class SegmentVoronoiBuilder {
public:
    SegmentVoronoiBuilder(std::vector<Segment> sites,
                          common::geometry::Aabb<double> bounds,
                          double maximumError);
    SegmentVoronoiBuilder(std::vector<Point> vertices,
                          common::geometry::Aabb<double> bounds,
                          double maximumError, PolygonSweepTag);

    SegmentVoronoiDiagram build();
    PolygonSweepDiagram buildPolygonDiagram();

private:
    using Real = long double;

    struct Vector {
        Real x = 0.0L;
        Real y = 0.0L;
    };

    struct LocalSegment {
        Vector start;
        Vector end;
    };

    struct SweepSegmentEndpoints {
        Vector first;
        Vector second;
    };

    enum class FeatureKind { Point, Segment };

    struct FeatureDomain {
        std::size_t segment = 0;
        int endpoint = -1;
    };

    struct Feature {
        FeatureKind kind = FeatureKind::Point;
        std::size_t inputIndex = 0;
        std::array<std::size_t, 2> owners{};
        std::size_t ownerCount = 0;
        std::array<FeatureDomain, 2> domains{};
        std::size_t domainCount = 0;
        Vector point;
        Vector sweepStart;
        Vector sweepEnd;
        Real activationX = 0.0L;
    };

    struct BeachSite {
        std::size_t feature = 0;
        bool inverse = false;
    };

    // The two slots are the two half-breakpoints born for a sweep edge. A
    // completed slot is a DCEL vertex; an incomplete slot is an open ray.
    struct EdgeEndpoint {
        Vector point;
        BeachSite opposite;
        bool completed = false;
    };

    struct EdgeRecord;

    struct HalfEdge {
        EdgeRecord* edge = nullptr;
        unsigned int endpoint = 0;
    };

    struct EdgeRecord {
        BeachSite first;
        BeachSite second;
        std::array<EdgeEndpoint, 2> endpoints;
        // A directed half-edge belongs to first (slot 0) or second (slot 1)
        // and leaves its stored vertex in this curve-parameter direction.
        std::array<int, 2> continuationDirections{};
        std::array<HalfEdge, 2> previous;
        std::array<HalfEdge, 2> next;
        bool removed = false;
    };

    struct SegmentState {
        std::size_t activationPoint = 0;
        std::size_t retirementPoint = 0;
        std::size_t interior = 0;
    };

    struct SiteEvent {
        BeachSite site;
        Vector point;
    };

    struct CircleCandidate {
        Vector center;
        Real radius = 0.0L;
        Real eventX = 0.0L;
    };

    struct CircleEvent;

    struct BeachNode {
        BeachSite left;
        BeachSite right;
        HalfEdge edge;
        CircleEvent* circle = nullptr;
        BeachNode* previous = nullptr;
        BeachNode* next = nullptr;
        BeachNode* lower = nullptr;
        BeachNode* upper = nullptr;
        BeachNode* parent = nullptr;
        int height = 1;
        int rootBranch = 0;
        Real lastTransverse = 0.0L;
        bool hasPosition = false;
        bool attached = false;
    };

    struct CircleEvent {
        CircleCandidate circle;
        BeachNode* node = nullptr;
        bool valid = true;
        std::uint64_t sequence = 0;
    };

    struct CircleAfter {
        bool operator()(const CircleEvent* left,
                        const CircleEvent* right) const {
            if (left->circle.eventX != right->circle.eventX) {
                return left->circle.eventX > right->circle.eventX;
            }
            if (left->circle.center.y != right->circle.center.y) {
                return left->circle.center.y > right->circle.center.y;
            }
            return left->sequence > right->sequence;
        }
    };

    // Beach breakpoints have an order determined by their neighboring arcs,
    // not by a geometric comparator whose result changes as the directrix
    // advances.  This intrusive AVL tree keeps that order structurally while
    // retaining predecessor/successor links for circle-event updates.
    class BeachTree {
    public:
        BeachNode* first() const { return m_first; }
        BeachNode* last() const { return m_last; }
        BeachNode* root() const { return m_root; }

        BeachNode* insertBefore(BeachNode* next, BeachSite left,
                                BeachSite right, HalfEdge edge) {
            auto storage = std::make_unique<BeachNode>();
            BeachNode* node = storage.get();
            node->left = left;
            node->right = right;
            node->edge = edge;
            node->attached = true;
            m_nodes.push_back(std::move(storage));

            node->next = next;
            node->previous = next == nullptr ? m_last : next->previous;
            if (node->previous != nullptr) {
                node->previous->next = node;
            } else {
                m_first = node;
            }
            if (next != nullptr) {
                next->previous = node;
            } else {
                m_last = node;
            }

            insertNode(node);
            return node;
        }

        void erase(BeachNode* node) {
            if (node == nullptr || !node->attached) return;

            if (node->previous != nullptr) {
                node->previous->next = node->next;
            } else {
                m_first = node->next;
            }
            if (node->next != nullptr) {
                node->next->previous = node->previous;
            } else {
                m_last = node->previous;
            }
            node->previous = nullptr;
            node->next = nullptr;
            node->attached = false;

            BeachNode* rebalanceStart = nullptr;
            if (node->lower == nullptr || node->upper == nullptr) {
                BeachNode* replacement =
                    node->lower != nullptr ? node->lower : node->upper;
                rebalanceStart =
                    node->parent != nullptr ? node->parent : replacement;
                replaceParent(node, replacement);
            } else {
                BeachNode* successor = node->upper;
                while (successor->lower != nullptr) {
                    successor = successor->lower;
                }
                if (successor->parent != node) {
                    BeachNode* successorParent = successor->parent;
                    replaceParent(successor, successor->upper);
                    successor->upper = node->upper;
                    successor->upper->parent = successor;
                    rebalanceStart = successorParent;
                } else {
                    rebalanceStart = successor;
                }
                replaceParent(node, successor);
                successor->lower = node->lower;
                successor->lower->parent = successor;
                updateHeight(successor);
            }

            node->lower = nullptr;
            node->upper = nullptr;
            node->parent = nullptr;
            node->height = 1;
            rebalanceFrom(rebalanceStart);
        }

        void clear() {
            m_nodes.clear();
            m_root = nullptr;
            m_first = nullptr;
            m_last = nullptr;
        }

    private:
        std::vector<std::unique_ptr<BeachNode>> m_nodes;
        BeachNode* m_root = nullptr;
        BeachNode* m_first = nullptr;
        BeachNode* m_last = nullptr;

        static int height(const BeachNode* node) {
            return node == nullptr ? 0 : node->height;
        }

        static int balanceFactor(const BeachNode* node) {
            return height(node->lower) - height(node->upper);
        }

        static void updateHeight(BeachNode* node) {
            if (node != nullptr) {
                node->height =
                    1 + std::max(height(node->lower), height(node->upper));
            }
        }

        void replaceParent(BeachNode* node, BeachNode* replacement) {
            BeachNode* parent = node->parent;
            if (parent == nullptr) {
                m_root = replacement;
            } else if (parent->lower == node) {
                parent->lower = replacement;
            } else {
                parent->upper = replacement;
            }
            if (replacement != nullptr) replacement->parent = parent;
        }

        void rotateLeft(BeachNode* node) {
            BeachNode* pivot = node->upper;
            BeachNode* parent = node->parent;
            node->upper = pivot->lower;
            if (node->upper != nullptr) node->upper->parent = node;
            pivot->lower = node;
            node->parent = pivot;
            pivot->parent = parent;
            if (parent == nullptr) {
                m_root = pivot;
            } else if (parent->lower == node) {
                parent->lower = pivot;
            } else {
                parent->upper = pivot;
            }
            updateHeight(node);
            updateHeight(pivot);
        }

        void rotateRight(BeachNode* node) {
            BeachNode* pivot = node->lower;
            BeachNode* parent = node->parent;
            node->lower = pivot->upper;
            if (node->lower != nullptr) node->lower->parent = node;
            pivot->upper = node;
            node->parent = pivot;
            pivot->parent = parent;
            if (parent == nullptr) {
                m_root = pivot;
            } else if (parent->lower == node) {
                parent->lower = pivot;
            } else {
                parent->upper = pivot;
            }
            updateHeight(node);
            updateHeight(pivot);
        }

        void rebalanceFrom(BeachNode* node) {
            while (node != nullptr) {
                updateHeight(node);
                if (balanceFactor(node) > 1) {
                    if (balanceFactor(node->lower) < 0) {
                        rotateLeft(node->lower);
                    }
                    rotateRight(node);
                } else if (balanceFactor(node) < -1) {
                    if (balanceFactor(node->upper) > 0) {
                        rotateRight(node->upper);
                    }
                    rotateLeft(node);
                }
                node = node->parent;
            }
        }

        void insertNode(BeachNode* node) {
            if (m_root == nullptr) {
                m_root = node;
                return;
            }

            BeachNode* parent = nullptr;
            if (node->previous != nullptr) {
                parent = node->previous;
                if (parent->upper != nullptr) {
                    parent = parent->upper;
                    while (parent->lower != nullptr) {
                        parent = parent->lower;
                    }
                    parent->lower = node;
                } else {
                    parent->upper = node;
                }
            } else {
                parent = node->next;
                if (parent->lower != nullptr) {
                    parent = parent->lower;
                    while (parent->upper != nullptr) {
                        parent = parent->upper;
                    }
                    parent->upper = node;
                } else {
                    parent->lower = node;
                }
            }
            node->parent = parent;
            rebalanceFrom(parent);
        }
    };

    struct Polynomial {
        Real quadratic = 0.0L;
        Real linear = 0.0L;
        Real constant = 0.0L;
    };

    struct Curve {
        bool curved = false;
        Vector quadratic;
        Vector linear;
        Vector constant;
        Vector lineOrigin;
        Vector lineDirection;
        Real parameterSign = 1.0L;
        Real segmentLength = 0.0L;
    };

    struct IntersectionEvent {
        Real x = 0.0L;
        Real y = 0.0L;
        std::size_t segment = 0;
        bool start = false;
    };

    struct SweepOutputEdge {
        std::size_t firstFeature = 0;
        std::size_t secondFeature = 0;
        bool curved = false;
        std::vector<Point> vertices;
    };

    struct ActiveSegmentLess {
        const SegmentVoronoiBuilder* owner = nullptr;
        const Real* sweepX = nullptr;

        bool operator()(std::size_t left, std::size_t right) const {
            if (left == right) return false;
            const Real leftY = owner->segmentYAt(left, *sweepX);
            const Real rightY = owner->segmentYAt(right, *sweepX);
            if (leftY != rightY) return leftY < rightY;
            const Real leftSlope = owner->segmentSweepSlope(left);
            const Real rightSlope = owner->segmentSweepSlope(right);
            if (leftSlope != rightSlope) return leftSlope < rightSlope;
            return left < right;
        }
    };

    std::vector<Segment> m_originalSites;
    std::vector<Point> m_polygonVertices;
    common::geometry::Aabb<double> m_inputBounds;
    std::vector<LocalSegment> m_localSites;
    std::vector<Vector> m_eventPoints;
    std::vector<Feature> m_features;
    std::vector<SegmentState> m_segmentStates;
    Vector m_sweepDirection{0.0L, -1.0L};
    Vector m_sweepTransverse{1.0L, 0.0L};
    Real m_originX = 0.0L;
    Real m_originY = 0.0L;
    Real m_minX = 0.0L;
    Real m_minY = 0.0L;
    Real m_maxX = 0.0L;
    Real m_maxY = 0.0L;
    Real m_scale = 1.0L;
    Real m_lengthTolerance = 1e-12L;
    Real m_areaTolerance = 1e-24L;
    Real m_eventTolerance = 1e-13L;
    Real m_distanceTolerance = 1e-20L;
    Real m_maximumError = 1e-4L;
    Real m_currentX = -std::numeric_limits<Real>::infinity();
    std::list<EdgeRecord> m_edges;
    std::vector<SiteEvent> m_siteEvents;
    std::size_t m_nextSiteEvent = 0;
    BeachTree m_beach;
    std::vector<BeachNode*> m_temporaryNodes;
    std::vector<std::unique_ptr<CircleEvent>> m_circleEvents;
    std::priority_queue<CircleEvent*, std::vector<CircleEvent*>, CircleAfter>
        m_circleQueue;
    std::uint64_t m_nextSequence = 0;
    bool m_polygonMode = false;

    static Vector add(const Vector& left, const Vector& right) {
        return {left.x + right.x, left.y + right.y};
    }

    static Vector subtract(const Vector& left, const Vector& right) {
        return {left.x - right.x, left.y - right.y};
    }

    static Vector multiply(const Vector& point, Real scalar) {
        return {point.x * scalar, point.y * scalar};
    }

    static Vector divide(const Vector& point, Real scalar) {
        return {point.x / scalar, point.y / scalar};
    }

    static Real dot(const Vector& left, const Vector& right) {
        return left.x * right.x + left.y * right.y;
    }

    static Real cross(const Vector& left, const Vector& right) {
        return left.x * right.y - left.y * right.x;
    }

    static Real squaredLength(const Vector& point) {
        return dot(point, point);
    }

    static Real length(const Vector& point) {
        return std::sqrt(squaredLength(point));
    }

    static Real squaredDistance(const Vector& left, const Vector& right) {
        return squaredLength(subtract(left, right));
    }

    static Vector normalized(const Vector& point) {
        return divide(point, length(point));
    }

    static Vector leftNormal(const Vector& vector) {
        return {-vector.y, vector.x};
    }

    static Vector rightNormal(const Vector& vector) {
        return {vector.y, -vector.x};
    }

    static bool finite(const Vector& point) {
        return std::isfinite(point.x) && std::isfinite(point.y);
    }

    static bool sameBeachSite(const BeachSite& left,
                              const BeachSite& right) {
        return left.feature == right.feature && left.inverse == right.inverse;
    }

    void validateScalars(double maximumError) const {
        if (!std::isfinite(maximumError) || maximumError < 0.0) {
            throw std::invalid_argument(
                "Segment Voronoi maximum error must be finite and nonnegative");
        }
        if (m_inputBounds.minX() >= m_inputBounds.maxX() ||
            m_inputBounds.minY() >= m_inputBounds.maxY()) {
            throw std::invalid_argument(
                "Segment Voronoi bounds must have positive area");
        }
        if (!std::isfinite(m_inputBounds.minX()) ||
            !std::isfinite(m_inputBounds.minY()) ||
            !std::isfinite(m_inputBounds.maxX()) ||
            !std::isfinite(m_inputBounds.maxY())) {
            throw std::invalid_argument("Segment Voronoi bounds must be finite");
        }
        for (const Segment& site : m_originalSites) {
            if (!std::isfinite(site.start.x()) ||
                !std::isfinite(site.start.y()) ||
                !std::isfinite(site.end.x()) ||
                !std::isfinite(site.end.y())) {
                throw std::invalid_argument(
                    "Segment Voronoi sites must be finite");
            }
            if (site.start.x() == site.end.x() &&
                site.start.y() == site.end.y()) {
                throw std::invalid_argument(
                    "Segment Voronoi sites must have positive length");
            }
        }
    }

    void initializeLocalCoordinates() {
        Real minimumX = m_inputBounds.minX();
        Real minimumY = m_inputBounds.minY();
        Real maximumX = m_inputBounds.maxX();
        Real maximumY = m_inputBounds.maxY();
        for (const Segment& site : m_originalSites) {
            minimumX = std::min(
                {minimumX, static_cast<Real>(site.start.x()),
                 static_cast<Real>(site.end.x())});
            minimumY = std::min(
                {minimumY, static_cast<Real>(site.start.y()),
                 static_cast<Real>(site.end.y())});
            maximumX = std::max(
                {maximumX, static_cast<Real>(site.start.x()),
                 static_cast<Real>(site.end.x())});
            maximumY = std::max(
                {maximumY, static_cast<Real>(site.start.y()),
                 static_cast<Real>(site.end.y())});
        }
        m_originX = minimumX * 0.5L + maximumX * 0.5L;
        m_originY = minimumY * 0.5L + maximumY * 0.5L;
        if (!std::isfinite(m_originX) || !std::isfinite(m_originY)) {
            throw std::invalid_argument(
                "Segment Voronoi normalization origin is not finite");
        }

        const auto radiusFromOrigin = [](Real value, Real origin) {
            return std::abs(value - origin);
        };
        m_scale = std::max(
            {radiusFromOrigin(static_cast<Real>(m_inputBounds.minX()),
                              m_originX),
             radiusFromOrigin(static_cast<Real>(m_inputBounds.maxX()),
                              m_originX),
             radiusFromOrigin(static_cast<Real>(m_inputBounds.minY()),
                              m_originY),
             radiusFromOrigin(static_cast<Real>(m_inputBounds.maxY()),
                              m_originY)});
        for (const Segment& site : m_originalSites) {
            m_scale = std::max(
                {m_scale,
                 radiusFromOrigin(static_cast<Real>(site.start.x()),
                                  m_originX),
                 radiusFromOrigin(static_cast<Real>(site.end.x()),
                                  m_originX),
                 radiusFromOrigin(static_cast<Real>(site.start.y()),
                                  m_originY),
                 radiusFromOrigin(static_cast<Real>(site.end.y()),
                                  m_originY)});
        }
        if (!(m_scale > 0.0L) || !std::isfinite(m_scale)) {
            throw std::invalid_argument(
                "Segment Voronoi normalization radius must be finite and positive");
        }

        const auto normalizeCoordinate = [&](Real value, Real origin) {
            const Real normalizedValue = (value - origin) / m_scale;
            if (!std::isfinite(normalizedValue)) {
                throw std::invalid_argument(
                    "Segment Voronoi normalized coordinates must be finite");
            }
            return normalizedValue;
        };
        m_minX = normalizeCoordinate(
            static_cast<Real>(m_inputBounds.minX()), m_originX);
        m_minY = normalizeCoordinate(
            static_cast<Real>(m_inputBounds.minY()), m_originY);
        m_maxX = normalizeCoordinate(
            static_cast<Real>(m_inputBounds.maxX()), m_originX);
        m_maxY = normalizeCoordinate(
            static_cast<Real>(m_inputBounds.maxY()), m_originY);

        m_localSites.reserve(m_originalSites.size());
        m_eventPoints.reserve(m_originalSites.size() * 2U);
        for (const Segment& site : m_originalSites) {
            const LocalSegment local{
                {normalizeCoordinate(static_cast<Real>(site.start.x()),
                                     m_originX),
                 normalizeCoordinate(static_cast<Real>(site.start.y()),
                                     m_originY)},
                {normalizeCoordinate(static_cast<Real>(site.end.x()),
                                     m_originX),
                 normalizeCoordinate(static_cast<Real>(site.end.y()),
                                     m_originY)}};
            if (!finite(local.start) || !finite(local.end)) {
                throw std::invalid_argument(
                    "Segment Voronoi normalized coordinates must be finite");
            }
            m_localSites.push_back(local);
            m_eventPoints.push_back(local.start);
            m_eventPoints.push_back(local.end);
        }
        std::sort(m_eventPoints.begin(), m_eventPoints.end(),
                  [](const Vector& left, const Vector& right) {
                      if (left.x != right.x) return left.x < right.x;
                      return left.y < right.y;
                  });
        m_eventPoints.erase(
            std::unique(m_eventPoints.begin(), m_eventPoints.end(),
                        [](const Vector& left, const Vector& right) {
                            return left.x == right.x && left.y == right.y;
                        }),
            m_eventPoints.end());
        m_lengthTolerance =
            256.0L * std::numeric_limits<Real>::epsilon();
        m_lengthTolerance = std::max(m_lengthTolerance, 1e-13L);
        m_areaTolerance = m_lengthTolerance;
        m_eventTolerance = 2e-13L;
        m_distanceTolerance = 2e-11L;
    }

    Vector toSweep(const Vector& point) const {
        return {dot(point, m_sweepDirection),
                dot(point, m_sweepTransverse)};
    }

    Vector fromSweep(const Vector& point) const {
        return add(multiply(m_sweepDirection, point.x),
                   multiply(m_sweepTransverse, point.y));
    }

    static Real normalizedAngle(Real angle) {
        const Real pi = std::acos(-1.0L);
        while (angle < 0.0L) angle += pi;
        while (angle >= pi) angle -= pi;
        return angle;
    }

    bool usableSweepDirection(const Vector& direction) const {
        std::vector<Real> eventCoordinates;
        eventCoordinates.reserve(m_eventPoints.size());
        for (const LocalSegment& segment : m_localSites) {
            const Vector delta = subtract(segment.end, segment.start);
            const Real segmentLength = length(delta);
            const Real projection =
                std::abs(dot(delta, direction)) / segmentLength;
            if (projection <= 1e-8L) return false;
            const Vector transverse{-direction.y, direction.x};
            const Real tangentY = dot(delta, transverse) / segmentLength;
            const Real forwardDenominator = 1.0L + tangentY;
            const Real reverseDenominator = 1.0L - tangentY;
            if (forwardDenominator <= 1e-10L ||
                reverseDenominator <= 1e-10L) {
                return false;
            }
        }
        for (const Vector& point : m_eventPoints) {
            eventCoordinates.push_back(dot(point, direction));
        }
        std::sort(eventCoordinates.begin(), eventCoordinates.end());
        for (std::size_t index = 1; index < eventCoordinates.size(); ++index) {
            if (eventCoordinates[index] - eventCoordinates[index - 1] <=
                m_eventTolerance) {
                return false;
            }
        }
        return true;
    }

    bool trySweepAngle(Real angle) {
        angle = normalizedAngle(angle);
        Vector direction{std::cos(angle), std::sin(angle)};
        if (direction.y > 0.0L ||
            (direction.y == 0.0L && direction.x < 0.0L)) {
            direction = multiply(direction, -1.0L);
        }
        if (!usableSweepDirection(direction)) return false;
        m_sweepDirection = direction;
        m_sweepTransverse = {-direction.y, direction.x};
        return true;
    }

    struct AngularInterval {
        Real first = 0.0L;
        Real second = 0.0L;
    };

    void addForbiddenInterval(std::vector<AngularInterval>& intervals,
                              Real center, Real halfWidth) const {
        const Real pi = std::acos(-1.0L);
        center = normalizedAngle(center);
        if (halfWidth >= pi * 0.5L) {
            intervals.push_back({0.0L, pi});
            return;
        }
        const Real first = center - halfWidth;
        const Real second = center + halfWidth;
        if (first < 0.0L) {
            intervals.push_back({0.0L, second});
            intervals.push_back({first + pi, pi});
        } else if (second > pi) {
            intervals.push_back({first, pi});
            intervals.push_back({0.0L, second - pi});
        } else {
            intervals.push_back({first, second});
        }
    }

    bool chooseSweepBasisFromForbiddenIntervals() {
        const Real pi = std::acos(-1.0L);
        constexpr Real denominatorTolerance = 1e-10L;
        const Real segmentProjectionTolerance = std::max(
            1e-8L,
            std::sqrt(2.0L * denominatorTolerance -
                      denominatorTolerance * denominatorTolerance));
        std::vector<AngularInterval> intervals;
        const std::size_t pointCount = m_eventPoints.size();
        if (pointCount > 1U &&
            pointCount >
                std::numeric_limits<std::size_t>::max() /
                    (pointCount - 1U)) {
            throw std::invalid_argument(
                "Segment Voronoi sweep direction search is too large");
        }
        intervals.reserve(m_localSites.size() * 2U +
                          pointCount * (pointCount - 1U));

        for (const LocalSegment& segment : m_localSites) {
            const Vector delta = subtract(segment.end, segment.start);
            const Real angle =
                std::atan2(delta.y, delta.x) + pi * 0.5L;
            addForbiddenInterval(
                intervals, angle,
                std::asin(std::min(1.0L, segmentProjectionTolerance)));
        }
        for (std::size_t first = 0; first < pointCount; ++first) {
            for (std::size_t second = first + 1U;
                 second < pointCount; ++second) {
                const Vector delta =
                    subtract(m_eventPoints[second], m_eventPoints[first]);
                const Real separation = length(delta);
                if (!(separation > m_eventTolerance)) return false;
                const Real ratio =
                    std::min(1.0L, m_eventTolerance / separation);
                addForbiddenInterval(
                    intervals,
                    std::atan2(delta.y, delta.x) + pi * 0.5L,
                    std::asin(ratio));
            }
        }

        std::sort(intervals.begin(), intervals.end(),
                  [](const AngularInterval& left,
                     const AngularInterval& right) {
                      if (left.first != right.first) {
                          return left.first < right.first;
                      }
                      return left.second < right.second;
                  });
        std::vector<AngularInterval> merged;
        for (const AngularInterval& interval : intervals) {
            if (merged.empty() ||
                interval.first > merged.back().second) {
                merged.push_back(interval);
            } else {
                merged.back().second =
                    std::max(merged.back().second, interval.second);
            }
        }

        std::vector<AngularInterval> gaps;
        Real previous = 0.0L;
        for (const AngularInterval& interval : merged) {
            if (interval.first > previous) {
                gaps.push_back({previous, interval.first});
            }
            previous = std::max(previous, interval.second);
        }
        if (previous < pi) gaps.push_back({previous, pi});
        std::sort(gaps.begin(), gaps.end(),
                  [](const AngularInterval& left,
                     const AngularInterval& right) {
                      return left.second - left.first >
                             right.second - right.first;
                  });
        for (const AngularInterval& gap : gaps) {
            const Real width = gap.second - gap.first;
            for (Real fraction : {0.5L, 0.25L, 0.75L}) {
                if (trySweepAngle(gap.first + width * fraction)) {
                    return true;
                }
            }
        }
        return false;
    }

    void chooseSweepBasis() {
        const Real pi = std::acos(-1.0L);
        if (trySweepAngle(-pi * 0.5L)) return;

        const Real inverseGoldenRatio =
            0.61803398874989484820458683436563811772L;
        constexpr std::size_t genericAttempts = 32U;
        Real fraction = 0.5L;
        for (std::size_t attempt = 0;
             attempt < genericAttempts; ++attempt) {
            fraction += inverseGoldenRatio;
            fraction -= std::floor(fraction);
            if (trySweepAngle(fraction * pi)) return;
        }
        if (chooseSweepBasisFromForbiddenIntervals()) return;
        throw std::invalid_argument(
            "Segment Voronoi input has no numerically stable sweep direction");
    }

    static Real orientation(const Vector& first, const Vector& second,
                            const Vector& third) {
        return cross(subtract(second, first), subtract(third, first));
    }

    bool pointOnSegment(const Vector& point,
                        const LocalSegment& segment) const {
        if (point.x < std::min(segment.start.x, segment.end.x) -
                          m_lengthTolerance ||
            point.x > std::max(segment.start.x, segment.end.x) +
                          m_lengthTolerance ||
            point.y < std::min(segment.start.y, segment.end.y) -
                          m_lengthTolerance ||
            point.y > std::max(segment.start.y, segment.end.y) +
                          m_lengthTolerance) {
            return false;
        }
        return std::abs(orientation(segment.start, segment.end, point)) <=
               m_areaTolerance;
    }

    bool segmentsIntersect(const LocalSegment& first,
                           const LocalSegment& second) const {
        const Real o1 =
            orientation(first.start, first.end, second.start);
        const Real o2 = orientation(first.start, first.end, second.end);
        const Real o3 =
            orientation(second.start, second.end, first.start);
        const Real o4 = orientation(second.start, second.end, first.end);
        const auto sign = [&](Real value) {
            if (value > m_areaTolerance) return 1;
            if (value < -m_areaTolerance) return -1;
            return 0;
        };
        const int s1 = sign(o1);
        const int s2 = sign(o2);
        const int s3 = sign(o3);
        const int s4 = sign(o4);
        if (s1 * s2 < 0 && s3 * s4 < 0) return true;
        return (s1 == 0 && pointOnSegment(second.start, first)) ||
               (s2 == 0 && pointOnSegment(second.end, first)) ||
               (s3 == 0 && pointOnSegment(first.start, second)) ||
               (s4 == 0 && pointOnSegment(first.end, second));
    }

    SweepSegmentEndpoints orderedSweepEndpoints(
        std::size_t segment) const {
        Vector first = toSweep(m_localSites[segment].start);
        Vector second = toSweep(m_localSites[segment].end);
        if (first.x > second.x ||
            (first.x == second.x && first.y > second.y)) {
            std::swap(first, second);
        }
        return {first, second};
    }

    Real segmentYAt(std::size_t segment, Real sweepX) const {
        const auto [first, second] = orderedSweepEndpoints(segment);
        const Real parameter =
            (sweepX - first.x) / (second.x - first.x);
        return first.y + (second.y - first.y) * parameter;
    }

    Real segmentSweepSlope(std::size_t segment) const {
        const auto [first, second] = orderedSweepEndpoints(segment);
        return (second.y - first.y) / (second.x - first.x);
    }

    void validateDisjointSites() const {
        std::vector<IntersectionEvent> events;
        events.reserve(m_localSites.size() * 2);
        for (std::size_t index = 0; index < m_localSites.size(); ++index) {
            const auto [first, second] = orderedSweepEndpoints(index);
            events.push_back({first.x, first.y, index, true});
            events.push_back({second.x, second.y, index, false});
        }
        std::sort(events.begin(), events.end(),
                  [](const IntersectionEvent& left,
                     const IntersectionEvent& right) {
                      if (left.x != right.x) return left.x < right.x;
                      if (left.start != right.start) return left.start;
                      if (left.y != right.y) return left.y < right.y;
                      return left.segment < right.segment;
                  });

        Real sweepX = events.empty() ? 0.0L : events.front().x;
        using ActiveSet = std::set<std::size_t, ActiveSegmentLess>;
        ActiveSet active(ActiveSegmentLess{this, &sweepX});
        std::vector<typename ActiveSet::iterator> positions(
            m_localSites.size(), active.end());
        const auto rejectIfIntersecting = [&](std::size_t first,
                                              std::size_t second) {
            if (segmentsIntersect(m_localSites[first],
                                  m_localSites[second])) {
                throw std::invalid_argument(
                    "Segment Voronoi sites must be pairwise disjoint");
            }
        };

        for (const IntersectionEvent& event : events) {
            sweepX = event.x;
            if (event.start) {
                auto inserted = active.insert(event.segment).first;
                positions[event.segment] = inserted;
                if (inserted != active.begin()) {
                    rejectIfIntersecting(*std::prev(inserted), *inserted);
                }
                const auto next = std::next(inserted);
                if (next != active.end()) {
                    rejectIfIntersecting(*inserted, *next);
                }
            } else {
                const auto position = positions[event.segment];
                if (position == active.end()) {
                    throw std::invalid_argument(
                        "Segment Voronoi intersection validation failed");
                }
                const auto next = std::next(position);
                if (position != active.begin() && next != active.end()) {
                    rejectIfIntersecting(*std::prev(position), *next);
                }
                active.erase(position);
                positions[event.segment] = active.end();
            }
        }
    }

    bool polygonEdgesAdjacent(std::size_t first,
                              std::size_t second) const {
        const std::size_t count = m_localSites.size();
        return (first + 1U) % count == second ||
               (second + 1U) % count == first;
    }

    bool validAdjacentPolygonContact(std::size_t first,
                                     std::size_t second) const {
        if (!polygonEdgesAdjacent(first, second)) return false;
        const std::size_t count = m_localSites.size();
        if ((second + 1U) % count == first) std::swap(first, second);
        const Vector shared = m_localSites[first].end;
        if (squaredDistance(shared, m_localSites[second].start) >
            m_lengthTolerance * m_lengthTolerance) {
            return false;
        }
        const Vector firstDirection =
            subtract(m_localSites[first].start, shared);
        const Vector secondDirection =
            subtract(m_localSites[second].end, shared);
        if (std::abs(cross(firstDirection, secondDirection)) >
            m_areaTolerance) {
            return true;
        }
        return dot(firstDirection, secondDirection) <= 0.0L;
    }

    void validatePolygonVertices() const {
        const std::size_t count = m_polygonVertices.size();
        if (count < 3U) {
            throw std::invalid_argument(
                "Polygon Voronoi input must have at least three vertices");
        }

        std::vector<Point> ordered = m_polygonVertices;
        std::sort(ordered.begin(), ordered.end(),
                  [](const Point& left, const Point& right) {
                      if (left.x() != right.x()) {
                          return left.x() < right.x();
                      }
                      return left.y() < right.y();
                  });
        for (std::size_t index = 1; index < ordered.size(); ++index) {
            if (ordered[index - 1U].x() == ordered[index].x() &&
                ordered[index - 1U].y() == ordered[index].y()) {
                throw std::invalid_argument(
                    "Polygon Voronoi vertices must not repeat");
            }
        }

        Real twiceArea = 0.0L;
        for (std::size_t vertex = 0; vertex < count; ++vertex) {
            const LocalSegment& edge = m_localSites[vertex];
            twiceArea += cross(edge.start, edge.end);
            const std::size_t previous = (vertex + count - 1U) % count;
            const Vector incoming =
                subtract(m_localSites[previous].start,
                         m_localSites[previous].end);
            const Vector outgoing =
                subtract(m_localSites[vertex].end,
                         m_localSites[vertex].start);
            if (std::abs(cross(incoming, outgoing)) <= m_areaTolerance &&
                dot(incoming, outgoing) > 0.0L) {
                throw std::invalid_argument(
                    "Polygon Voronoi input contains collinear backtracking");
            }
        }
        if (std::abs(twiceArea) <= m_areaTolerance) {
            throw std::invalid_argument(
                "Polygon Voronoi input must enclose nonzero area");
        }
    }

    void validateSimplePolygon() const {
        validatePolygonVertices();

        std::vector<IntersectionEvent> events;
        events.reserve(m_localSites.size() * 2U);
        for (std::size_t index = 0; index < m_localSites.size(); ++index) {
            const auto [first, second] = orderedSweepEndpoints(index);
            events.push_back({first.x, first.y, index, true});
            events.push_back({second.x, second.y, index, false});
        }
        std::sort(events.begin(), events.end(),
                  [](const IntersectionEvent& left,
                     const IntersectionEvent& right) {
                      if (left.x != right.x) return left.x < right.x;
                      if (left.start != right.start) return left.start;
                      if (left.y != right.y) return left.y < right.y;
                      return left.segment < right.segment;
                  });

        Real sweepX = events.front().x;
        using ActiveSet = std::set<std::size_t, ActiveSegmentLess>;
        ActiveSet active(ActiveSegmentLess{this, &sweepX});
        std::vector<typename ActiveSet::iterator> positions(
            m_localSites.size(), active.end());
        const auto rejectIfInvalidIntersection =
            [&](std::size_t first, std::size_t second) {
                if (!segmentsIntersect(m_localSites[first],
                                       m_localSites[second])) {
                    return;
                }
                if (validAdjacentPolygonContact(first, second)) return;
                throw std::invalid_argument(
                    "Polygon Voronoi input must be a simple polygon");
            };

        for (const IntersectionEvent& event : events) {
            sweepX = event.x;
            if (event.start) {
                auto inserted = active.insert(event.segment).first;
                positions[event.segment] = inserted;
                if (inserted != active.begin()) {
                    rejectIfInvalidIntersection(
                        *std::prev(inserted), *inserted);
                }
                const auto next = std::next(inserted);
                if (next != active.end()) {
                    rejectIfInvalidIntersection(*inserted, *next);
                }
            } else {
                const auto position = positions[event.segment];
                if (position == active.end()) {
                    throw std::invalid_argument(
                        "Polygon Voronoi intersection validation failed");
                }
                const auto next = std::next(position);
                if (position != active.begin() && next != active.end()) {
                    rejectIfInvalidIntersection(
                        *std::prev(position), *next);
                }
                active.erase(position);
                positions[event.segment] = active.end();
            }
        }
    }

    std::size_t addPointFeature(
        std::size_t inputIndex, const Vector& point,
        std::array<std::size_t, 2> owners, std::size_t ownerCount,
        std::array<FeatureDomain, 2> domains, std::size_t domainCount) {
        const std::size_t index = m_features.size();
        Feature value;
        value.kind = FeatureKind::Point;
        value.inputIndex = inputIndex;
        value.owners = owners;
        value.ownerCount = ownerCount;
        value.domains = domains;
        value.domainCount = domainCount;
        value.point = point;
        value.activationX = point.x;
        m_features.push_back(value);
        return index;
    }

    std::size_t addSegmentFeature(std::size_t inputIndex,
                                  std::size_t segment,
                                  const Vector& sweepStart,
                                  const Vector& sweepEnd) {
        const std::size_t index = m_features.size();
        Feature value;
        value.kind = FeatureKind::Segment;
        value.inputIndex = inputIndex;
        value.owners = {segment, 0U};
        value.ownerCount = 1U;
        value.domains = {{{segment, -1}, {}}};
        value.domainCount = 1U;
        value.sweepStart = sweepStart;
        value.sweepEnd = sweepEnd;
        value.activationX = sweepStart.x;
        m_features.push_back(value);
        return index;
    }

    void initializeSegmentFeaturesAndEvents() {
        m_segmentStates.resize(m_localSites.size());
        for (std::size_t index = 0; index < m_localSites.size(); ++index) {
            const Vector publicStart = toSweep(m_localSites[index].start);
            const Vector publicEnd = toSweep(m_localSites[index].end);
            Vector activation = publicStart;
            Vector retirement = publicEnd;
            int activationEndpoint = 0;
            int retirementEndpoint = 1;
            if (activation.x > retirement.x ||
                (activation.x == retirement.x &&
                 activation.y > retirement.y)) {
                std::swap(activation, retirement);
                std::swap(activationEndpoint, retirementEndpoint);
            }

            SegmentState& state = m_segmentStates[index];
            state.activationPoint = addPointFeature(
                index, activation, {index, 0U}, 1U,
                {{{index, activationEndpoint}, {}}}, 1U);
            state.retirementPoint = addPointFeature(
                index, retirement, {index, 0U}, 1U,
                {{{index, retirementEndpoint}, {}}}, 1U);
            state.interior =
                addSegmentFeature(index, index, activation, retirement);

            m_siteEvents.push_back(
                {{state.activationPoint, false}, activation});
            m_siteEvents.push_back(
                {{state.retirementPoint, false}, retirement});
            m_siteEvents.push_back(
                {{state.interior, false}, activation});
        }
    }

    void initializePolygonFeaturesAndEvents() {
        const std::size_t count = m_localSites.size();
        m_segmentStates.resize(count);
        std::vector<std::size_t> vertexFeatures(count);
        for (std::size_t vertex = 0; vertex < count; ++vertex) {
            const std::size_t previous = (vertex + count - 1U) % count;
            const Vector point = toSweep(m_localSites[vertex].start);
            vertexFeatures[vertex] = addPointFeature(
                vertex, point, {previous, vertex}, 2U,
                {{{previous, 1}, {vertex, 0}}}, 2U);
            m_siteEvents.push_back(
                {{vertexFeatures[vertex], false}, point});
        }

        for (std::size_t edge = 0; edge < count; ++edge) {
            const std::size_t nextVertex = (edge + 1U) % count;
            Vector activation = toSweep(m_localSites[edge].start);
            Vector retirement = toSweep(m_localSites[edge].end);
            std::size_t activationVertex = edge;
            std::size_t retirementVertex = nextVertex;
            if (activation.x > retirement.x ||
                (activation.x == retirement.x &&
                 activation.y > retirement.y)) {
                std::swap(activation, retirement);
                std::swap(activationVertex, retirementVertex);
            }

            SegmentState& state = m_segmentStates[edge];
            state.activationPoint = vertexFeatures[activationVertex];
            state.retirementPoint = vertexFeatures[retirementVertex];
            state.interior =
                addSegmentFeature(edge, edge, activation, retirement);
            m_siteEvents.push_back(
                {{state.interior, false}, activation});
        }
    }

    const Feature& feature(const BeachSite& site) const {
        return m_features[site.feature];
    }

    SweepSegmentEndpoints directedSegment(
        const BeachSite& site) const {
        const Feature& value = feature(site);
        if (!site.inverse) return {value.sweepStart, value.sweepEnd};
        return {value.sweepEnd, value.sweepStart};
    }

    Vector sitePoint(const BeachSite& site) const {
        return feature(site).point;
    }

    Vector segmentNormal(const BeachSite& site) const {
        const auto [first, second] = directedSegment(site);
        return rightNormal(normalized(subtract(second, first)));
    }

    Real signedLineDistance(const BeachSite& site,
                            const Vector& point) const {
        const auto [first, second] = directedSegment(site);
        const Vector normal =
            rightNormal(normalized(subtract(second, first)));
        return dot(normal, subtract(point, first));
    }

    std::pair<HalfEdge, HalfEdge> createEdge(BeachSite first,
                                             BeachSite second) {
        m_edges.emplace_back();
        EdgeRecord* edge = &m_edges.back();
        edge->first = first;
        edge->second = second;
        edge->continuationDirections = {1, -1};
        return {{edge, 0U}, {edge, 1U}};
    }

    void completeBreakpoint(const HalfEdge& halfEdge, const Vector& point,
                            BeachSite opposite) {
        if (halfEdge.edge == nullptr) {
            return;
        }
        EdgeEndpoint& endpoint =
            halfEdge.edge->endpoints[halfEdge.endpoint];
        if (endpoint.completed) {
            if (squaredDistance(endpoint.point, point) <=
                m_lengthTolerance * m_lengthTolerance) {
                return;
            }
            throw std::invalid_argument(
                "Segment Voronoi breakpoint completed more than once");
        }
        endpoint.point = point;
        endpoint.opposite = opposite;
        endpoint.completed = true;
    }

    static HalfEdge twin(const HalfEdge& edge) {
        return {edge.edge, edge.endpoint ^ 1U};
    }

    static bool sameHalfEdge(const HalfEdge& left, const HalfEdge& right) {
        return left.edge == right.edge && left.endpoint == right.endpoint;
    }

    void linkNext(const HalfEdge& first, const HalfEdge& second) {
        if (first.edge == nullptr || second.edge == nullptr) return;
        first.edge->next[first.endpoint] = second;
        second.edge->previous[second.endpoint] = first;
    }

    void assignVertex(const HalfEdge& edge, const Vector& point) {
        if (edge.edge == nullptr) return;
        EdgeEndpoint& endpoint = edge.edge->endpoints[edge.endpoint];
        endpoint.point = point;
        endpoint.completed = true;
    }

    std::optional<HalfEdge> rotationNext(const HalfEdge& edge) const {
        if (edge.edge == nullptr) return std::nullopt;
        const HalfEdge previous = edge.edge->previous[edge.endpoint];
        if (previous.edge == nullptr) return std::nullopt;
        return twin(previous);
    }

    std::optional<HalfEdge> rotationPrevious(const HalfEdge& edge) const {
        if (edge.edge == nullptr) return std::nullopt;
        const HalfEdge next = edge.edge->next[edge.endpoint ^ 1U];
        if (next.edge == nullptr) return std::nullopt;
        return next;
    }

    void collapseDegenerateEdge(EdgeRecord& edge) {
        if (edge.removed || !edge.endpoints[0].completed ||
            !edge.endpoints[1].completed ||
            squaredDistance(edge.endpoints[0].point,
                            edge.endpoints[1].point) >
                m_lengthTolerance * m_lengthTolerance) {
            return;
        }
        const HalfEdge first{&edge, 0U};
        const HalfEdge second{&edge, 1U};
        const auto firstPrevious = rotationPrevious(first);
        const auto firstNext = rotationNext(first);
        const auto secondPrevious = rotationPrevious(second);
        const auto secondNext = rotationNext(second);
        if (firstPrevious.has_value() && firstNext.has_value() &&
            secondPrevious.has_value() && secondNext.has_value()) {
            HalfEdge current = *secondNext;
            const HalfEdge stop = twin(first);
            const std::size_t maximumSteps = m_edges.size() * 2U + 1U;
            for (std::size_t step = 0;
                 step < maximumSteps && !sameHalfEdge(current, stop);
                 ++step) {
                assignVertex(current, edge.endpoints[0].point);
                const auto next = rotationNext(current);
                if (!next.has_value()) break;
                current = *next;
            }

            linkNext(twin(*firstNext), *secondPrevious);
            linkNext(twin(*secondNext), *firstPrevious);
        }
        edge.removed = true;
    }

    void collapseDegenerateEdges() {
        for (EdgeRecord& edge : m_edges) collapseDegenerateEdge(edge);
    }

    Real arcX(const BeachSite& site, Real directrix,
              Real transverse) const {
        const Feature& value = feature(site);
        if (value.kind == FeatureKind::Point) {
            const Real denominator =
                2.0L * (directrix - value.point.x);
            if (denominator <= m_eventTolerance) {
                throw std::invalid_argument(
                    "Unsupported coincident segment Voronoi site events");
            }
            const Real delta = transverse - value.point.y;
            return (directrix + value.point.x) * 0.5L -
                   delta * delta / denominator;
        }

        const auto [first, second] = directedSegment(site);
        const Vector tangent = normalized(subtract(second, first));
        const Vector normal = rightNormal(tangent);
        const Real denominator = 1.0L + normal.x;
        if (denominator <= 1e-10L) {
            throw std::invalid_argument(
                "Segment Voronoi sweep direction is numerically singular");
        }
        return (directrix - normal.y * transverse +
                dot(normal, first)) /
               denominator;
    }

    bool isPointSite(const BeachSite& site) const {
        return feature(site).kind == FeatureKind::Point;
    }

    bool isSegmentSite(const BeachSite& site) const {
        return feature(site).kind == FeatureKind::Segment;
    }

    Vector siteEventPoint0(const BeachSite& site) const {
        if (isPointSite(site)) return sitePoint(site);
        return directedSegment(site).first;
    }

    Vector siteEventPoint1(const BeachSite& site) const {
        if (isPointSite(site)) return sitePoint(site);
        return directedSegment(site).second;
    }

    BeachSite inverseSite(BeachSite site) const {
        if (isSegmentSite(site)) site.inverse = !site.inverse;
        return site;
    }

    int orientationSign(const Vector& first, const Vector& second,
                        const Vector& third) const {
        const Real value = orientation(first, second, third);
        if (value < 0.0L) return -1;
        if (value > 0.0L) return 1;
        return 0;
    }

    bool siteEventBefore(const SiteEvent& left,
                         const SiteEvent& right) const {
        if (left.point.x != right.point.x) {
            return left.point.x < right.point.x;
        }
        if (isPointSite(left.site)) {
            if (isPointSite(right.site)) {
                return left.point.y < right.point.y;
            }
            return true;
        }
        if (isPointSite(right.site)) return false;
        if (left.point.y != right.point.y) {
            return left.point.y < right.point.y;
        }
        const int angularOrder =
            orientationSign(siteEventPoint1(left.site),
                            siteEventPoint0(left.site),
                            siteEventPoint1(right.site));
        if (angularOrder != 0) {
            return angularOrder > 0;
        }
        return left.site.feature < right.site.feature;
    }

    static std::vector<Real> solveQuadratic(Real quadratic, Real linear,
                                            Real constant) {
        constexpr Real inputEpsilon =
            static_cast<Real>(std::numeric_limits<double>::epsilon());
        const Real scale =
            std::max({1.0L, std::abs(quadratic), std::abs(linear),
                      std::abs(constant)});
        const Real epsilon =
            1024.0L * inputEpsilon * scale;
        if (std::abs(quadratic) <= epsilon) {
            if (std::abs(linear) <= epsilon) return {};
            return {-constant / linear};
        }
        Real discriminant =
            linear * linear - 4.0L * quadratic * constant;
        const Real discriminantTolerance =
            4096.0L * inputEpsilon *
            std::max({1.0L, linear * linear,
                      std::abs(4.0L * quadratic * constant)});
        if (discriminant < -discriminantTolerance) return {};
        if (discriminant < 0.0L) discriminant = 0.0L;
        const Real root = std::sqrt(discriminant);
        const Real stable =
            -0.5L * (linear + std::copysign(root, linear));
        if (stable == 0.0L) {
            return {-linear / (2.0L * quadratic)};
        }
        const Real first = stable / quadratic;
        const Real second = constant / stable;
        if (std::abs(first - second) <=
            1024.0L * inputEpsilon *
                std::max({1.0L, std::abs(first), std::abs(second)})) {
            return {(first + second) * 0.5L};
        }
        return first < second ? std::vector<Real>{first, second}
                              : std::vector<Real>{second, first};
    }

    Polynomial beachFront(const BeachSite& site, Real directrix) const {
        const Feature& value = feature(site);
        if (value.kind == FeatureKind::Point) {
            const Real delta = directrix - value.point.x;
            if (delta <= 0.0L) {
                throw std::invalid_argument(
                    "Segment Voronoi beach query reached an unprocessed point");
            }
            return {-1.0L / (2.0L * delta), value.point.y / delta,
                    (directrix + value.point.x) * 0.5L -
                        value.point.y * value.point.y / (2.0L * delta)};
        }

        const auto [first, second] = directedSegment(site);
        const Vector normal =
            rightNormal(normalized(subtract(second, first)));
        const Real denominator = 1.0L + normal.x;
        if (denominator <= 1e-10L) {
            throw std::invalid_argument(
                "Segment Voronoi sweep direction is numerically singular");
        }
        return {0.0L, -normal.y / denominator,
                (directrix + dot(normal, first)) / denominator};
    }

    std::vector<Real> breakpointRoots(const BeachNode& node,
                                      Real directrix) const {
        if (node.left.feature == node.right.feature) {
            const Feature& value = feature(node.left);
            if (value.kind != FeatureKind::Segment) return {};
            const Vector direction =
                subtract(value.sweepEnd, value.sweepStart);
            if (std::abs(direction.x) <= m_eventTolerance) return {};
            return {value.sweepStart.y +
                    (directrix - value.sweepStart.x) * direction.y /
                        direction.x};
        }

        const Polynomial left = beachFront(node.left, directrix);
        const Polynomial right = beachFront(node.right, directrix);
        return solveQuadratic(left.quadratic - right.quadratic,
                              left.linear - right.linear,
                              left.constant - right.constant);
    }

    Real breakpointTransverse(BeachNode& node, Real directrix) {
        const std::vector<Real> roots = breakpointRoots(node, directrix);
        if (roots.empty()) {
            if (node.hasPosition) return node.lastTransverse;
            throw std::invalid_argument(
                "Segment Voronoi beach breakpoint has no real branch");
        }

        Real selected = roots.front();
        if (roots.size() > 1U) {
            if (node.rootBranch < 0) {
                selected = roots.front();
            } else if (node.rootBranch > 0) {
                selected = roots.back();
            } else if (node.hasPosition) {
                if (std::abs(roots.front() - node.lastTransverse) <=
                    std::abs(roots.back() - node.lastTransverse)) {
                    selected = roots.front();
                    node.rootBranch = -1;
                } else {
                    selected = roots.back();
                    node.rootBranch = 1;
                }
            }
        }
        node.lastTransverse = selected;
        node.hasPosition = true;
        return selected;
    }

    void seedBreakpoint(BeachNode* node, int rootBranch, Real transverse) {
        node->rootBranch = rootBranch;
        node->lastTransverse = transverse;
        node->hasPosition = true;
    }

    void orientCircleBreakpoint(BeachNode* node, BeachNode* before,
                                BeachNode* after, Real eventDirectrix,
                                Real eventTransverse) {
        const Real probe = std::max(
            1e-17L,
            128.0L * std::numeric_limits<Real>::epsilon() *
                std::max(1.0L, std::abs(eventDirectrix)));
        const Real directrix = eventDirectrix + probe;
        const std::vector<Real> roots =
            breakpointRoots(*node, directrix);
        if (roots.size() < 2U) {
            seedBreakpoint(node, 0,
                           roots.empty() ? eventTransverse : roots.front());
            return;
        }

        Real lower = -std::numeric_limits<Real>::infinity();
        Real upper = std::numeric_limits<Real>::infinity();
        if (before != nullptr && before->attached) {
            lower = breakpointTransverse(*before, directrix);
        }
        if (after != nullptr && after->attached) {
            upper = breakpointTransverse(*after, directrix);
        }
        const Real tolerance = m_lengthTolerance * 64.0L;
        std::size_t selected =
            std::abs(roots.front() - eventTransverse) <=
                    std::abs(roots.back() - eventTransverse)
                ? 0U
                : 1U;
        bool foundInNeighbors = false;
        for (std::size_t index = 0; index < roots.size(); ++index) {
            if (roots[index] >= lower - tolerance &&
                roots[index] <= upper + tolerance) {
                if (!foundInNeighbors ||
                    std::abs(roots[index] - eventTransverse) <
                        std::abs(roots[selected] - eventTransverse)) {
                    selected = index;
                    foundInNeighbors = true;
                }
            }
        }
        seedBreakpoint(node, selected == 0U ? -1 : 1,
                       roots[selected]);
    }

    BeachNode* firstBreakpointAtOrAfter(const SiteEvent& event) {
        const Real probe = std::max(
            1e-17L,
            128.0L * std::numeric_limits<Real>::epsilon() *
                std::max(1.0L, std::abs(event.point.x)));
        const Real directrix = event.point.x + probe;
        BeachNode* current = m_beach.root();
        BeachNode* result = nullptr;
        while (current != nullptr) {
            const Real transverse =
                breakpointTransverse(*current, directrix);
            if (transverse < event.point.y) {
                current = current->upper;
            } else {
                result = current;
                current = current->lower;
            }
        }
        return result;
    }

    void appendPointLineCandidates(const Vector& point,
                                   const BeachSite& segment,
                                   const Vector& equationNormal,
                                   Real equationOffset,
                                   int selectedRoot,
                                   std::vector<CircleCandidate>& output) const {
        const Real normalLengthSquared = squaredLength(equationNormal);
        if (normalLengthSquared <= m_areaTolerance) return;
        const Vector linePoint =
            multiply(equationNormal,
                     equationOffset / normalLengthSquared);
        const Vector lineDirection =
            normalized(leftNormal(equationNormal));
        const Vector segmentNormalValue = segmentNormal(segment);
        const auto directed = directedSegment(segment);
        const Real segmentOffset =
            dot(segmentNormalValue, directed.first);
        const Real radiusBase =
            dot(segmentNormalValue, linePoint) - segmentOffset;
        const Real radiusSlope =
            dot(segmentNormalValue, lineDirection);
        const Vector pointOffset = subtract(linePoint, point);
        const Real quadratic =
            squaredLength(lineDirection) - radiusSlope * radiusSlope;
        const Real linear =
            2.0L * (dot(pointOffset, lineDirection) -
                    radiusBase * radiusSlope);
        const Real constant =
            squaredLength(pointOffset) - radiusBase * radiusBase;
        const std::vector<Real> roots =
            solveQuadratic(quadratic, linear, constant);
        std::size_t begin = 0;
        std::size_t end = roots.size();
        if (roots.size() > 1 && selectedRoot != 0) {
            begin = selectedRoot < 0 ? 0 : roots.size() - 1;
            end = begin + 1;
        }
        for (std::size_t index = begin; index < end; ++index) {
            const Real parameter = roots[index];
            const Vector center =
                add(linePoint, multiply(lineDirection, parameter));
            const Real radius = radiusBase + radiusSlope * parameter;
            output.push_back({center, radius, center.x + radius});
        }
    }

    void formPPP(const BeachSite& first, const BeachSite& second,
                 const BeachSite& third,
                 std::vector<CircleCandidate>& output) const {
        const Vector p1 = sitePoint(first);
        const Vector p2 = sitePoint(second);
        const Vector p3 = sitePoint(third);
        const Vector d12 = subtract(p2, p1);
        const Vector d13 = subtract(p3, p1);
        const Real determinant = 2.0L * cross(d12, d13);
        if (std::abs(determinant) <= m_areaTolerance) return;
        const Real q12 = squaredLength(d12);
        const Real q13 = squaredLength(d13);
        const Vector center{
            p1.x + (q12 * d13.y - q13 * d12.y) / determinant,
            p1.y + (q13 * d12.x - q12 * d13.x) / determinant};
        const Real radius = std::sqrt(squaredDistance(center, p1));
        output.push_back({center, radius, center.x + radius});
    }

    void formPPS(const std::array<BeachSite, 3>& sites,
                 std::vector<CircleCandidate>& output) const {
        std::array<BeachSite, 2> points{};
        BeachSite segment;
        std::size_t pointCount = 0;
        int segmentIndex = 0;
        for (std::size_t index = 0; index < sites.size(); ++index) {
            const BeachSite site = sites[index];
            if (feature(site).kind == FeatureKind::Point) {
                points[pointCount++] = site;
            } else {
                segment = site;
                segmentIndex = static_cast<int>(index) + 1;
            }
        }
        const Vector first = sitePoint(points[0]);
        const Vector second = sitePoint(points[1]);
        const Vector normal = multiply(subtract(second, first), 2.0L);
        const Real offset =
            squaredLength(second) - squaredLength(first);
        appendPointLineCandidates(
            first, segment, normal, offset,
            segmentIndex == 2 ? 1 : -1, output);
    }

    void formPSS(const std::array<BeachSite, 3>& sites,
                 std::vector<CircleCandidate>& output) const {
        BeachSite point;
        std::array<BeachSite, 2> segments{};
        std::size_t segmentCount = 0;
        int pointIndex = 0;
        for (std::size_t index = 0; index < sites.size(); ++index) {
            const BeachSite site = sites[index];
            if (feature(site).kind == FeatureKind::Point) {
                point = site;
                pointIndex = static_cast<int>(index) + 1;
            } else {
                segments[segmentCount++] = site;
            }
        }
        const Vector firstNormal = segmentNormal(segments[0]);
        const Vector secondNormal = segmentNormal(segments[1]);
        const auto firstDirected = directedSegment(segments[0]);
        const auto secondDirected = directedSegment(segments[1]);
        const Real firstOffset =
            dot(firstNormal, firstDirected.first);
        const Real secondOffset =
            dot(secondNormal, secondDirected.first);
        appendPointLineCandidates(
            sitePoint(point), segments[0],
            subtract(firstNormal, secondNormal),
            firstOffset - secondOffset,
            pointIndex == 2 ? -1 : 1, output);
    }

    static bool solveThreeByThree(
        std::array<std::array<Real, 4>, 3> matrix,
        std::array<Real, 3>& solution) {
        for (std::size_t column = 0; column < 3; ++column) {
            std::size_t pivot = column;
            for (std::size_t row = column + 1; row < 3; ++row) {
                if (std::abs(matrix[row][column]) >
                    std::abs(matrix[pivot][column])) {
                    pivot = row;
                }
            }
            const Real pivotScale =
                std::max({1.0L, std::abs(matrix[pivot][0]),
                          std::abs(matrix[pivot][1]),
                          std::abs(matrix[pivot][2])});
            if (std::abs(matrix[pivot][column]) <=
                2048.0L * std::numeric_limits<Real>::epsilon() *
                    pivotScale) {
                return false;
            }
            if (pivot != column) std::swap(matrix[pivot], matrix[column]);
            for (std::size_t row = column + 1; row < 3; ++row) {
                const Real factor =
                    matrix[row][column] / matrix[column][column];
                for (std::size_t entry = column; entry < 4; ++entry) {
                    matrix[row][entry] -=
                        factor * matrix[column][entry];
                }
            }
        }
        for (int row = 2; row >= 0; --row) {
            Real value = matrix[static_cast<std::size_t>(row)][3];
            for (std::size_t column =
                     static_cast<std::size_t>(row) + 1;
                 column < 3; ++column) {
                value -= matrix[static_cast<std::size_t>(row)][column] *
                         solution[column];
            }
            solution[static_cast<std::size_t>(row)] =
                value /
                matrix[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(row)];
        }
        return true;
    }

    void formSSS(const std::array<BeachSite, 3>& sites,
                 std::vector<CircleCandidate>& output) const {
        std::array<std::array<Real, 4>, 3> matrix{};
        for (std::size_t index = 0; index < 3; ++index) {
            const Vector normal = segmentNormal(sites[index]);
            const auto directed = directedSegment(sites[index]);
            matrix[index] = {normal.x, normal.y, -1.0L,
                             dot(normal, directed.first)};
        }
        std::array<Real, 3> solution{};
        if (!solveThreeByThree(matrix, solution)) return;
        const Vector center{solution[0], solution[1]};
        const Real radius = solution[2];
        output.push_back({center, radius, center.x + radius});
    }

    Real projectionParameter(std::size_t segmentIndex,
                             const Vector& sweepPoint) const {
        const LocalSegment& segment = m_localSites[segmentIndex];
        const Vector point = fromSweep(sweepPoint);
        const Vector direction =
            subtract(segment.end, segment.start);
        return dot(subtract(point, segment.start), direction) /
               squaredLength(direction);
    }

    bool featureDomainContains(const BeachSite& site,
                               const Vector& sweepPoint,
                               Real toleranceMultiplier = 1.0L) const {
        const Feature& value = feature(site);
        if (value.kind == FeatureKind::Point) {
            for (std::size_t index = 0; index < value.domainCount; ++index) {
                const FeatureDomain& domain = value.domains[index];
                const LocalSegment& segment =
                    m_localSites[domain.segment];
                const Real segmentLength =
                    length(subtract(segment.end, segment.start));
                const Real tolerance =
                    toleranceMultiplier * m_lengthTolerance /
                    std::max(segmentLength, m_lengthTolerance);
                const Real parameter =
                    projectionParameter(domain.segment, sweepPoint);
                if ((domain.endpoint == 0 && parameter > tolerance) ||
                    (domain.endpoint == 1 &&
                     parameter < 1.0L - tolerance)) {
                    return false;
                }
            }
            return true;
        }

        const std::size_t segmentIndex = value.domains[0].segment;
        const LocalSegment& segment = m_localSites[segmentIndex];
        const Real segmentLength =
            length(subtract(segment.end, segment.start));
        const Real tolerance =
            toleranceMultiplier * m_lengthTolerance /
            std::max(segmentLength, m_lengthTolerance);
        const Real parameter =
            projectionParameter(segmentIndex, sweepPoint);
        if (parameter < -tolerance || parameter > 1.0L + tolerance) {
            return false;
        }
        return signedLineDistance(site, sweepPoint) >=
               -toleranceMultiplier * m_lengthTolerance;
    }

    Real primitiveDistanceSquared(const BeachSite& site,
                                  const Vector& sweepPoint) const {
        const Feature& value = feature(site);
        if (value.kind == FeatureKind::Point) {
            return squaredDistance(sweepPoint, value.point);
        }
        const Real distance = signedLineDistance(site, sweepPoint);
        return distance * distance;
    }

    bool validCircleCandidate(const CircleCandidate& candidate,
                              const std::array<BeachSite, 3>& sites,
                              Real currentX) const {
        if (!finite(candidate.center) ||
            !std::isfinite(candidate.radius) ||
            !std::isfinite(candidate.eventX) ||
            candidate.radius <= m_lengthTolerance ||
            candidate.eventX + m_eventTolerance < currentX) {
            return false;
        }

        Real maximumActivation =
            -std::numeric_limits<Real>::infinity();
        Real referenceDistance = -1.0L;
        for (const BeachSite& site : sites) {
            maximumActivation =
                std::max(maximumActivation, feature(site).activationX);
            const Real distance =
                primitiveDistanceSquared(site, candidate.center);
            if (referenceDistance < 0.0L) {
                referenceDistance = distance;
            } else if (std::abs(distance - referenceDistance) >
                       std::max(m_distanceTolerance,
                                referenceDistance * 2e-10L)) {
                return false;
            }
        }
        if (candidate.eventX <=
            maximumActivation + m_eventTolerance) {
            return false;
        }

        const Real eventFirstX =
            arcX(sites[0], candidate.eventX, candidate.center.y);
        const Real eventMiddleX =
            arcX(sites[1], candidate.eventX, candidate.center.y);
        const Real eventLastX =
            arcX(sites[2], candidate.eventX, candidate.center.y);
        const Real eventTolerance = m_lengthTolerance * 512.0L;
        const bool result =
            std::abs(eventFirstX - candidate.center.x) <= eventTolerance &&
            std::abs(eventMiddleX - candidate.center.x) <= eventTolerance &&
            std::abs(eventLastX - candidate.center.x) <= eventTolerance;
        return result;
    }

    bool ppsCircleExists(const BeachSite& firstPoint,
                         const BeachSite& secondPoint,
                         const BeachSite& segment,
                         int segmentIndex) const {
        if (segmentIndex == 2) {
            return !(siteEventPoint0(segment).x ==
                         sitePoint(firstPoint).x &&
                     siteEventPoint0(segment).y ==
                         sitePoint(firstPoint).y &&
                     siteEventPoint1(segment).x ==
                         sitePoint(secondPoint).x &&
                     siteEventPoint1(segment).y ==
                         sitePoint(secondPoint).y);
        }
        const int firstOrientation =
            orientationSign(sitePoint(firstPoint),
                            sitePoint(secondPoint),
                            siteEventPoint0(segment));
        const int secondOrientation =
            orientationSign(sitePoint(firstPoint),
                            sitePoint(secondPoint),
                            siteEventPoint1(segment));
        if (segmentIndex == 1 &&
            sitePoint(firstPoint).x >= sitePoint(secondPoint).x) {
            return firstOrientation == -1;
        }
        if (segmentIndex == 3 &&
            sitePoint(secondPoint).x >= sitePoint(firstPoint).x) {
            return secondOrientation == -1;
        }
        return firstOrientation == -1 || secondOrientation == -1;
    }

    bool pssCircleExists(const BeachSite& point,
                         const BeachSite& firstSegment,
                         const BeachSite& secondSegment,
                         int pointIndex) const {
        if (firstSegment.feature == secondSegment.feature) return false;
        if (pointIndex == 2) {
            if (!firstSegment.inverse && secondSegment.inverse) {
                return false;
            }
            if (firstSegment.inverse == secondSegment.inverse &&
                orientationSign(siteEventPoint0(firstSegment),
                                sitePoint(point),
                                siteEventPoint1(secondSegment)) != -1) {
                return false;
            }
        }
        return true;
    }

    bool circleExists(BeachSite first, BeachSite middle,
                      BeachSite last) const {
        const bool firstPoint = isPointSite(first);
        const bool middlePoint = isPointSite(middle);
        const bool lastPoint = isPointSite(last);
        if (firstPoint && middlePoint && lastPoint) {
            return orientationSign(sitePoint(first), sitePoint(middle),
                                   sitePoint(last)) == -1;
        }
        if (firstPoint && middlePoint) {
            return ppsCircleExists(first, middle, last, 3);
        }
        if (firstPoint && lastPoint) {
            return ppsCircleExists(first, last, middle, 2);
        }
        if (middlePoint && lastPoint) {
            return ppsCircleExists(middle, last, first, 1);
        }
        if (firstPoint) {
            return pssCircleExists(first, middle, last, 1);
        }
        if (middlePoint) {
            return pssCircleExists(middle, first, last, 2);
        }
        if (lastPoint) {
            return pssCircleExists(last, first, middle, 3);
        }
        return first.feature != middle.feature &&
               middle.feature != last.feature;
    }

    std::optional<CircleCandidate> formCircle(
        BeachSite first, BeachSite middle, BeachSite last,
        Real currentX) const {
        if (sameBeachSite(first, middle) ||
            sameBeachSite(middle, last) ||
            sameBeachSite(first, last) ||
            !circleExists(first, middle, last)) {
            return std::nullopt;
        }
        const std::array<BeachSite, 3> sites{first, middle, last};
        std::size_t pointCount = 0;
        for (const BeachSite& site : sites) {
            if (feature(site).kind == FeatureKind::Point) ++pointCount;
        }

        std::vector<CircleCandidate> candidates;
        if (pointCount == 3) {
            formPPP(first, middle, last, candidates);
        } else if (pointCount == 2) {
            formPPS(sites, candidates);
        } else if (pointCount == 1) {
            formPSS(sites, candidates);
        } else {
            formSSS(sites, candidates);
        }

        std::optional<CircleCandidate> selected;
        for (const CircleCandidate& candidate : candidates) {
            const bool valid =
                validCircleCandidate(candidate, sites, currentX);
            if (!valid) continue;
            if (!selected.has_value() ||
                candidate.eventX < selected->eventX) {
                selected = candidate;
            } else if (std::abs(candidate.eventX - selected->eventX) <=
                           m_eventTolerance &&
                       squaredDistance(candidate.center, selected->center) >
                           m_lengthTolerance * m_lengthTolerance * 64.0L) {
                throw std::invalid_argument(
                    "Unsupported ambiguous segment Voronoi tangent event");
            }
        }
        return selected;
    }

    void deactivateCircle(BeachNode* node) {
        if (node != nullptr && node->circle != nullptr) {
            node->circle->valid = false;
            node->circle = nullptr;
        }
    }

    void activateCircle(BeachSite first, BeachSite middle, BeachSite last,
                        BeachNode* bisectorNode) {
        if (bisectorNode == nullptr ||
            !bisectorNode->attached ||
            bisectorNode->edge.edge == nullptr) {
            return;
        }
        const auto candidate =
            formCircle(first, middle, last, m_currentX);
        if (!candidate.has_value()) return;
        m_circleEvents.push_back(std::make_unique<CircleEvent>());
        CircleEvent* event = m_circleEvents.back().get();
        event->circle = *candidate;
        event->node = bisectorNode;
        event->sequence = m_nextSequence++;
        bisectorNode->circle = event;
        m_circleQueue.push(event);
    }

    BeachNode* insertBeachArc(BeachSite leftArc, BeachSite rightArc,
                              const SiteEvent& siteEvent,
                              BeachNode* position) {
        BeachSite leftNodeLeft = leftArc;
        BeachSite leftNodeRight = siteEvent.site;
        BeachSite rightNodeLeft = siteEvent.site;
        BeachSite rightNodeRight = rightArc;
        if (isSegmentSite(siteEvent.site)) {
            rightNodeLeft = inverseSite(rightNodeLeft);
        }

        const auto edges = createEdge(rightArc, siteEvent.site);
        BeachNode* rightNode =
            m_beach.insertBefore(position, rightNodeLeft, rightNodeRight,
                                 edges.second);
        seedBreakpoint(rightNode, 1, siteEvent.point.y);
        BeachNode* beforeRight = rightNode;
        if (isSegmentSite(siteEvent.site)) {
            BeachNode* temporary = m_beach.insertBefore(
                rightNode, siteEvent.site, inverseSite(siteEvent.site),
                {});
            seedBreakpoint(temporary, 0, siteEvent.point.y);
            const Feature& segmentFeature = feature(siteEvent.site);
            m_temporaryNodes[segmentFeature.domains[0].segment] = temporary;
            beforeRight = temporary;
        }
        BeachNode* leftNode = m_beach.insertBefore(
            beforeRight, leftNodeLeft, leftNodeRight, edges.first);
        seedBreakpoint(leftNode, -1, siteEvent.point.y);
        return leftNode;
    }

    void initializeBeach() {
        std::sort(m_siteEvents.begin(), m_siteEvents.end(),
                  [&](const SiteEvent& left, const SiteEvent& right) {
                      return siteEventBefore(left, right);
                  });
        if (m_siteEvents.size() < 2) {
            m_nextSiteEvent = m_siteEvents.size();
            return;
        }
        m_temporaryNodes.resize(m_localSites.size());

        const SiteEvent& first = m_siteEvents[0];
        const SiteEvent& second = m_siteEvents[1];
        if (!isPointSite(first.site)) {
            throw std::invalid_argument(
                "Segment Voronoi sweep must begin with an endpoint");
        }
        m_currentX = second.point.x;
        insertBeachArc(first.site, first.site, second, nullptr);
        m_nextSiteEvent = 2;
        if (m_polygonMode) {
            while (m_nextSiteEvent < m_siteEvents.size() &&
                   m_siteEvents[m_nextSiteEvent].point.x == first.point.x &&
                   m_siteEvents[m_nextSiteEvent].point.y == first.point.y) {
                processSiteEvent();
            }
        }
    }

    void removeRetiredTemporaryNode(const SiteEvent& siteEvent) {
        const Feature& pointFeature = feature(siteEvent.site);
        for (std::size_t index = 0;
             index < pointFeature.domainCount; ++index) {
            const std::size_t segment =
                pointFeature.domains[index].segment;
            SegmentState& state = m_segmentStates[segment];
            if (siteEvent.site.feature != state.retirementPoint) continue;
            BeachNode*& stored = m_temporaryNodes[segment];
            if (stored == nullptr || !stored->attached) {
                stored = nullptr;
                continue;
            }
            if (stored->edge.edge == nullptr) m_beach.erase(stored);
            stored = nullptr;
        }
    }

    void processSiteEvent() {
        const SiteEvent siteEvent = m_siteEvents[m_nextSiteEvent];
        std::size_t last = m_nextSiteEvent + 1;
        if (!isPointSite(siteEvent.site)) {
            while (last < m_siteEvents.size() &&
                   isSegmentSite(m_siteEvents[last].site) &&
                   m_siteEvents[last].point.x == siteEvent.point.x &&
                   m_siteEvents[last].point.y == siteEvent.point.y) {
                ++last;
            }
        }
        if (isPointSite(siteEvent.site)) {
            removeRetiredTemporaryNode(siteEvent);
        }

        BeachNode* right = firstBreakpointAtOrAfter(siteEvent);
        for (; m_nextSiteEvent < last; ++m_nextSiteEvent) {
            const SiteEvent current = m_siteEvents[m_nextSiteEvent];
            BeachNode* left = right;
            if (right == nullptr) {
                left = m_beach.last();
                const BeachSite arc = left->right;
                right = insertBeachArc(arc, arc, current, nullptr);
                activateCircle(left->left, left->right, current.site, right);
            } else if (right == m_beach.first()) {
                const BeachSite arc = right->left;
                left = insertBeachArc(arc, arc, current, right);
                BeachSite circleSite = current.site;
                if (isSegmentSite(circleSite)) {
                    circleSite = inverseSite(circleSite);
                }
                activateCircle(circleSite, right->left, right->right, right);
                right = left;
            } else {
                const BeachSite secondArc = right->left;
                const BeachSite third = right->right;
                deactivateCircle(right);
                left = left->previous;
                const BeachSite firstArc = left->right;
                const BeachSite first = left->left;
                BeachNode* inserted =
                    insertBeachArc(firstArc, secondArc, current, right);
                activateCircle(first, firstArc, current.site, inserted);
                BeachSite circleSite = current.site;
                if (isSegmentSite(circleSite)) {
                    circleSite = inverseSite(circleSite);
                }
                activateCircle(circleSite, secondArc, third, right);
                right = inserted;
            }
        }
    }

    void processCircleEvent() {
        CircleEvent* event = m_circleQueue.top();
        m_circleQueue.pop();
        if (!event->valid) return;
        BeachNode* last = event->node;
        if (last == nullptr || !last->attached || last->circle != event ||
            last == m_beach.first()) {
            return;
        }
        BeachNode* first = last->previous;
        const BeachSite site1 = first->left;
        const BeachSite site2 = first->right;
        BeachSite site3 = last->right;
        if (isPointSite(site1) && isSegmentSite(site3) &&
            siteEventPoint1(site3).x == sitePoint(site1).x &&
            siteEventPoint1(site3).y == sitePoint(site1).y) {
            site3 = inverseSite(site3);
        }

        const HalfEdge firstBisector = first->edge;
        const HalfEdge secondBisector = last->edge;
        deactivateCircle(first);
        last->circle = nullptr;
        BeachNode* after = last->next;
        if (after != nullptr) deactivateCircle(after);
        BeachNode* before = first->previous;

        const Vector vertex = fromSweep(event->circle.center);
        completeBreakpoint(firstBisector, vertex, site3);
        completeBreakpoint(secondBisector, vertex, site1);
        const auto newEdge = createEdge(site1, site3);
        completeBreakpoint(newEdge.second, vertex, site2);
        linkNext(newEdge.first, firstBisector);
        linkNext(twin(firstBisector), secondBisector);
        linkNext(twin(secondBisector), newEdge.second);

        // Reinsert the replacement between the unchanged neighbors.  The
        // AVL stores positional order, so no geometric key is ever mutated
        // while it is resident in the tree.
        m_beach.erase(last);
        m_beach.erase(first);
        BeachNode* inserted =
            m_beach.insertBefore(after, site1, site3, newEdge.first);
        orientCircleBreakpoint(inserted, before, after,
                               event->circle.eventX,
                               event->circle.center.y);

        if (before != nullptr) {
            activateCircle(before->left, site1, site3, inserted);
        }
        if (after != nullptr) {
            activateCircle(site1, site3, after->right, after);
        }
    }

    void discardInactiveCircles() {
        while (!m_circleQueue.empty() && !m_circleQueue.top()->valid) {
            m_circleQueue.pop();
        }
    }

    void processPolygonVertexEvent() {
        const Vector vertex = m_siteEvents[m_nextSiteEvent].point;
        processSiteEvent();
        while (m_nextSiteEvent < m_siteEvents.size() &&
               m_siteEvents[m_nextSiteEvent].point.x == vertex.x &&
               m_siteEvents[m_nextSiteEvent].point.y == vertex.y) {
            processSiteEvent();
        }
    }

    void processSweep() {
        initializeBeach();
        discardInactiveCircles();
        while (m_nextSiteEvent < m_siteEvents.size() ||
               !m_circleQueue.empty()) {
            const bool haveSite =
                m_nextSiteEvent < m_siteEvents.size();
            const bool processSite =
                haveSite &&
                (m_circleQueue.empty() ||
                 m_siteEvents[m_nextSiteEvent].point.x <=
                     m_circleQueue.top()->circle.eventX);
            if (processSite) {
                m_currentX =
                    m_siteEvents[m_nextSiteEvent].point.x;
                if (m_polygonMode) {
                    processPolygonVertexEvent();
                } else {
                    processSiteEvent();
                }
            } else {
                m_currentX =
                    m_circleQueue.top()->circle.eventX;
                processCircleEvent();
            }
            discardInactiveCircles();
        }
        collapseDegenerateEdges();
        m_beach.clear();
        std::fill(m_temporaryNodes.begin(), m_temporaryNodes.end(), nullptr);
    }

    Vector primitiveGradient(const BeachSite& site,
                             const Vector& sweepPoint) const {
        const Feature& value = feature(site);
        if (value.kind == FeatureKind::Point) {
            return multiply(subtract(sweepPoint, value.point), 2.0L);
        }
        const Vector normal = segmentNormal(site);
        return multiply(normal,
                        2.0L * signedLineDistance(site, sweepPoint));
    }

    std::optional<Curve> curveFor(const EdgeRecord& edge) const {
        const Feature& first = feature(edge.first);
        const Feature& second = feature(edge.second);
        Curve curve;

        if (first.kind == FeatureKind::Point &&
            second.kind == FeatureKind::Point) {
            const Vector firstPoint = first.point;
            const Vector secondPoint = second.point;
            const Vector origin =
                multiply(add(firstPoint, secondPoint), 0.5L);
            Vector direction =
                normalized(leftNormal(subtract(secondPoint, firstPoint)));
            curve.lineOrigin = fromSweep(origin);
            curve.lineDirection = fromSweep(direction);
            curve.linear = curve.lineDirection;
            curve.constant = curve.lineOrigin;
            return curve;
        }

        if (first.kind == FeatureKind::Segment &&
            second.kind == FeatureKind::Segment) {
            const Vector firstNormal = segmentNormal(edge.first);
            const Vector secondNormal = segmentNormal(edge.second);
            const auto firstDirected = directedSegment(edge.first);
            const auto secondDirected = directedSegment(edge.second);
            const Vector normal =
                subtract(firstNormal, secondNormal);
            const Real normalSquared = squaredLength(normal);
            if (normalSquared <= m_areaTolerance) return std::nullopt;
            const Real offset =
                dot(firstNormal, firstDirected.first) -
                dot(secondNormal, secondDirected.first);
            const Vector origin =
                multiply(normal, offset / normalSquared);
            const Vector direction = normalized(leftNormal(normal));
            curve.lineOrigin = fromSweep(origin);
            curve.lineDirection = fromSweep(direction);
            curve.linear = curve.lineDirection;
            curve.constant = curve.lineOrigin;
            return curve;
        }

        const BeachSite pointSite =
            first.kind == FeatureKind::Point ? edge.first : edge.second;
        const BeachSite segmentSite =
            first.kind == FeatureKind::Segment ? edge.first : edge.second;
        const Vector focus = sitePoint(pointSite);
        const auto [lineStart, lineEnd] = directedSegment(segmentSite);
        const Vector tangent =
            normalized(subtract(lineEnd, lineStart));
        const Vector normal = rightNormal(tangent);
        const Vector relative = subtract(focus, lineStart);
        const Real focusTangent = dot(relative, tangent);
        const Real focusNormal = dot(relative, normal);
        if (focusNormal <= m_lengthTolerance) {
            return std::nullopt;
        }

        Vector quadratic =
            divide(normal, 2.0L * focusNormal);
        Vector linear =
            subtract(tangent,
                     multiply(normal, focusTangent / focusNormal));
        Vector constant = add(
            lineStart,
            multiply(normal,
                     (focusTangent * focusTangent +
                      focusNormal * focusNormal) /
                         (2.0L * focusNormal)));

        const Real sampleQ =
            length(subtract(lineEnd, lineStart)) * 0.5L;
        const Vector sampleSweep =
            add(add(multiply(quadratic, sampleQ * sampleQ),
                    multiply(linear, sampleQ)),
                constant);
        const Vector derivativeSweep =
            add(multiply(quadratic, 2.0L * sampleQ), linear);
        const Vector gradient = subtract(
            primitiveGradient(edge.first, sampleSweep),
            primitiveGradient(edge.second, sampleSweep));
        Real parameterSign = 1.0L;
        if (dot(gradient, leftNormal(derivativeSweep)) >= 0.0L) {
            parameterSign = -1.0L;
            linear = multiply(linear, -1.0L);
        }

        curve.curved = true;
        curve.quadratic = fromSweep(quadratic);
        curve.linear = fromSweep(linear);
        curve.constant = fromSweep(constant);
        curve.parameterSign = parameterSign;
        curve.segmentLength = length(subtract(lineEnd, lineStart));
        return curve;
    }

    Vector evaluate(const Curve& curve, Real parameter) const {
        return add(add(multiply(curve.quadratic, parameter * parameter),
                       multiply(curve.linear, parameter)),
                   curve.constant);
    }

    bool clipLineToBounds(const Curve& curve, Real& minimum,
                          Real& maximum) const {
        minimum = -std::numeric_limits<Real>::infinity();
        maximum = std::numeric_limits<Real>::infinity();
        const auto update = [&](Real origin, Real direction, Real lower,
                                Real upper) {
            if (std::abs(direction) <= m_lengthTolerance) {
                return origin >= lower - m_lengthTolerance &&
                       origin <= upper + m_lengthTolerance;
            }
            Real first = (lower - origin) / direction;
            Real second = (upper - origin) / direction;
            if (first > second) std::swap(first, second);
            minimum = std::max(minimum, first);
            maximum = std::min(maximum, second);
            return minimum <= maximum + m_lengthTolerance;
        };
        return update(curve.lineOrigin.x, curve.lineDirection.x,
                      m_minX, m_maxX) &&
               update(curve.lineOrigin.y, curve.lineDirection.y,
                      m_minY, m_maxY);
    }

    bool featureDomainContainsLocal(const BeachSite& site,
                                    const Vector& localPoint,
                                    Real multiplier = 1.0L) const {
        return featureDomainContains(site, toSweep(localPoint), multiplier);
    }

    void addPolynomialRoots(const Polynomial& polynomial, Real minimum,
                            Real maximum,
                            std::vector<Real>& breakpoints) const {
        for (Real root :
             solveQuadratic(polynomial.quadratic, polynomial.linear,
                            polynomial.constant)) {
            if (root > minimum + m_lengthTolerance &&
                root < maximum - m_lengthTolerance) {
                breakpoints.push_back(root);
            }
        }
    }

    Polynomial coordinatePolynomial(const Curve& curve, bool xCoordinate,
                                    Real offset = 0.0L) const {
        return {xCoordinate ? curve.quadratic.x : curve.quadratic.y,
                xCoordinate ? curve.linear.x : curve.linear.y,
                (xCoordinate ? curve.constant.x : curve.constant.y) -
                    offset};
    }

    Polynomial projectionPolynomial(const Curve& curve,
                                    const LocalSegment& segment,
                                    Real offset) const {
        const Vector direction =
            subtract(segment.end, segment.start);
        return {dot(curve.quadratic, direction),
                dot(curve.linear, direction),
                dot(subtract(curve.constant, segment.start), direction) -
                    offset};
    }

    Polynomial signedDistancePolynomial(const Curve& curve,
                                        const BeachSite& site) const {
        const auto [firstSweep, secondSweep] = directedSegment(site);
        const Vector normal =
            rightNormal(normalized(subtract(secondSweep, firstSweep)));
        const Vector firstLocal = fromSweep(firstSweep);
        const Vector normalLocal = fromSweep(normal);
        return {dot(curve.quadratic, normalLocal),
                dot(curve.linear, normalLocal),
                dot(subtract(curve.constant, firstLocal), normalLocal)};
    }

    bool insideBounds(const Vector& point) const {
        return point.x >= m_minX - m_lengthTolerance &&
               point.x <= m_maxX + m_lengthTolerance &&
               point.y >= m_minY - m_lengthTolerance &&
               point.y <= m_maxY + m_lengthTolerance;
    }

    bool validOutputPoint(const EdgeRecord& edge,
                          const Vector& point) const {
        return insideBounds(point) &&
               featureDomainContainsLocal(edge.first, point, 32.0L) &&
               featureDomainContainsLocal(edge.second, point, 32.0L);
    }

    std::vector<std::pair<Real, Real>> visibleIntervals(
        const EdgeRecord& edge, const Curve& curve,
        Real minimum, Real maximum) const {
        std::vector<Real> breakpoints{minimum, maximum};
        if (curve.curved) {
            addPolynomialRoots(
                coordinatePolynomial(curve, true, m_minX),
                minimum, maximum, breakpoints);
            addPolynomialRoots(
                coordinatePolynomial(curve, true, m_maxX),
                minimum, maximum, breakpoints);
            addPolynomialRoots(
                coordinatePolynomial(curve, false, m_minY),
                minimum, maximum, breakpoints);
            addPolynomialRoots(
                coordinatePolynomial(curve, false, m_maxY),
                minimum, maximum, breakpoints);
        }

        for (const BeachSite& site : {edge.first, edge.second}) {
            const Feature& value = feature(site);
            for (std::size_t index = 0;
                 index < value.domainCount; ++index) {
                const LocalSegment& segment =
                    m_localSites[value.domains[index].segment];
                const Real segmentLengthSquared =
                    squaredLength(subtract(segment.end, segment.start));
                addPolynomialRoots(
                    projectionPolynomial(curve, segment, 0.0L),
                    minimum, maximum, breakpoints);
                addPolynomialRoots(
                    projectionPolynomial(curve, segment,
                                         segmentLengthSquared),
                    minimum, maximum, breakpoints);
            }
            if (value.kind == FeatureKind::Segment) {
                addPolynomialRoots(
                    signedDistancePolynomial(curve, site),
                    minimum, maximum, breakpoints);
            }
        }

        std::sort(breakpoints.begin(), breakpoints.end());
        breakpoints.erase(
            std::unique(breakpoints.begin(), breakpoints.end(),
                        [&](Real left, Real right) {
                            return std::abs(left - right) <=
                                   m_lengthTolerance;
                        }),
            breakpoints.end());

        std::vector<std::pair<Real, Real>> intervals;
        for (std::size_t index = 1; index < breakpoints.size(); ++index) {
            const Real first = breakpoints[index - 1];
            const Real second = breakpoints[index];
            if (second - first <= m_lengthTolerance) continue;
            const Real middle = first + (second - first) * 0.5L;
            if (!validOutputPoint(edge, evaluate(curve, middle))) continue;
            if (!intervals.empty() &&
                std::abs(intervals.back().second - first) <=
                    m_lengthTolerance * 8.0L) {
                intervals.back().second = second;
            } else {
                intervals.push_back({first, second});
            }
        }
        return intervals;
    }

    Vector clampToBounds(Vector point) const {
        point.x = std::clamp(point.x, m_minX, m_maxX);
        point.y = std::clamp(point.y, m_minY, m_maxY);
        return point;
    }

    static Real distanceToChord(const Vector& point,
                                const Vector& first,
                                const Vector& second) {
        const Vector direction = subtract(second, first);
        const Real lengthSquared = squaredLength(direction);
        if (lengthSquared == 0.0L) {
            return std::sqrt(squaredDistance(point, first));
        }
        const Real parameter =
            std::clamp(dot(subtract(point, first), direction) /
                           lengthSquared,
                       0.0L, 1.0L);
        const Vector closest =
            add(first, multiply(direction, parameter));
        return std::sqrt(squaredDistance(point, closest));
    }

    void appendAdaptive(const Curve& curve, Real firstParameter,
                        const Vector& first, Real secondParameter,
                        const Vector& second, int depth,
                        std::vector<Vector>& output) const {
        const Real middleParameter =
            firstParameter + (secondParameter - firstParameter) * 0.5L;
        const Vector middle = evaluate(curve, middleParameter);
        if (depth < 32 &&
            distanceToChord(middle, first, second) >
                m_maximumError) {
            appendAdaptive(curve, firstParameter, first,
                           middleParameter, middle, depth + 1, output);
            appendAdaptive(curve, middleParameter, middle,
                           secondParameter, second, depth + 1, output);
        } else {
            output.push_back(second);
        }
    }

    Point publicPoint(Vector point) const {
        point = clampToBounds(point);
        const Real x = std::clamp(point.x * m_scale + m_originX,
                                  static_cast<Real>(m_inputBounds.minX()),
                                  static_cast<Real>(m_inputBounds.maxX()));
        const Real y = std::clamp(point.y * m_scale + m_originY,
                                  static_cast<Real>(m_inputBounds.minY()),
                                  static_cast<Real>(m_inputBounds.maxY()));
        return {static_cast<double>(x), static_cast<double>(y)};
    }

    bool featureOwnersOverlap(const Feature& first,
                              const Feature& second) const {
        for (std::size_t left = 0; left < first.ownerCount; ++left) {
            for (std::size_t right = 0; right < second.ownerCount; ++right) {
                if (first.owners[left] == second.owners[right]) return true;
            }
        }
        return false;
    }

    void appendResultEdge(const EdgeRecord& source, const Curve& curve,
                          Real firstParameter, Real secondParameter,
                          std::vector<SweepOutputEdge>& result) const {
        Vector first = clampToBounds(evaluate(curve, firstParameter));
        Vector second = clampToBounds(evaluate(curve, secondParameter));
        if (squaredDistance(first, second) <=
            m_lengthTolerance * m_lengthTolerance) {
            return;
        }

        std::vector<Vector> localVertices{first};
        if (curve.curved) {
            appendAdaptive(curve, firstParameter, first,
                           secondParameter, second, 0, localVertices);
        } else {
            localVertices.push_back(second);
        }

        std::vector<Point> vertices;
        vertices.reserve(localVertices.size());
        for (const Vector& vertex : localVertices) {
            const Point converted = publicPoint(vertex);
            const Real deltaX =
                vertices.empty()
                    ? 0.0L
                    : (static_cast<Real>(vertices.back().x()) -
                       static_cast<Real>(converted.x())) /
                          m_scale;
            const Real deltaY =
                vertices.empty()
                    ? 0.0L
                    : (static_cast<Real>(vertices.back().y()) -
                       static_cast<Real>(converted.y())) /
                          m_scale;
            if (vertices.empty() ||
                deltaX * deltaX + deltaY * deltaY >
                    m_lengthTolerance * m_lengthTolerance) {
                vertices.push_back(converted);
            }
        }
        if (vertices.size() < 2) return;

        result.push_back(
            {source.first.feature, source.second.feature,
             curve.curved, std::move(vertices)});
    }

    std::optional<Real> parameterAt(const Curve& curve,
                                    const Vector& point) const {
        if (!curve.curved) {
            const Real directionLengthSquared =
                squaredLength(curve.lineDirection);
            if (directionLengthSquared <=
                m_lengthTolerance * m_lengthTolerance) {
                return std::nullopt;
            }
            return dot(subtract(point, curve.lineOrigin),
                       curve.lineDirection) /
                   directionLengthSquared;
        }

        const bool useX =
            std::abs(curve.quadratic.x) + std::abs(curve.linear.x) >=
            std::abs(curve.quadratic.y) + std::abs(curve.linear.y);
        const Polynomial primary =
            coordinatePolynomial(curve, useX,
                                 useX ? point.x : point.y);
        const Polynomial secondary =
            coordinatePolynomial(curve, !useX,
                                 useX ? point.y : point.x);
        std::vector<Real> candidates =
            solveQuadratic(primary.quadratic, primary.linear,
                           primary.constant);
        for (Real candidate :
             solveQuadratic(secondary.quadratic, secondary.linear,
                            secondary.constant)) {
            candidates.push_back(candidate);
        }

        std::optional<Real> result;
        Real bestError = std::numeric_limits<Real>::infinity();
        for (Real candidate : candidates) {
            const Real error =
                squaredDistance(evaluate(curve, candidate), point);
            if (error < bestError) {
                bestError = error;
                result = candidate;
            }
        }
        const Real tolerance = m_lengthTolerance * 4096.0L;
        if (!result.has_value() || bestError > tolerance * tolerance) {
            return std::nullopt;
        }
        return result;
    }

    std::optional<std::pair<Real, Real>> topologyRange(
        const EdgeRecord& edge, const Curve& curve) const {
        Real minimum = -std::numeric_limits<Real>::infinity();
        Real maximum = std::numeric_limits<Real>::infinity();
        if (curve.curved) {
            minimum = std::min(0.0L,
                               curve.parameterSign * curve.segmentLength);
            maximum = std::max(0.0L,
                               curve.parameterSign * curve.segmentLength);
        }

        std::array<Real, 2> parameters{};
        std::size_t completed = 0;
        for (std::size_t index = 0; index < edge.endpoints.size(); ++index) {
            const EdgeEndpoint& endpoint = edge.endpoints[index];
            if (!endpoint.completed) continue;
            const auto parameter = parameterAt(curve, endpoint.point);
            if (!parameter.has_value()) return std::nullopt;
            parameters[index] = *parameter;
            ++completed;
        }

        if (completed == 2U) {
            minimum = std::max(minimum,
                               std::min(parameters[0], parameters[1]));
            maximum = std::min(maximum,
                               std::max(parameters[0], parameters[1]));
        } else if (completed == 1U) {
            const std::size_t index =
                edge.endpoints[0].completed ? 0U : 1U;
            const int direction = edge.continuationDirections[index];
            if (direction < 0) {
                maximum = std::min(maximum, parameters[index]);
            } else if (direction > 0) {
                minimum = std::max(minimum, parameters[index]);
            } else {
                return std::nullopt;
            }
        }

        if (maximum - minimum <= m_lengthTolerance) return std::nullopt;
        return std::pair<Real, Real>{minimum, maximum};
    }

    void emitSweepEdge(const EdgeRecord& edge,
                       std::vector<SweepOutputEdge>& result) const {
        if (edge.removed) return;
        if (featureOwnersOverlap(feature(edge.first),
                                 feature(edge.second))) {
            return;
        }

        // A point on another segment's supporting line can generate a
        // zero-width primitive pair during the sweep.  It has no parabolic
        // branch for a disjoint finite segment, so omit only that secondary
        // pair instead of rejecting the complete input.
        const auto curveValue = curveFor(edge);
        if (!curveValue.has_value()) return;
        const Curve& curve = *curveValue;
        const auto topology = topologyRange(edge, curve);
        if (!topology.has_value()) return;

        Real minimum = topology->first;
        Real maximum = topology->second;
        if (!curve.curved) {
            Real boundsMinimum;
            Real boundsMaximum;
            if (!clipLineToBounds(curve, boundsMinimum, boundsMaximum)) {
                return;
            }
            minimum = std::max(minimum, boundsMinimum);
            maximum = std::min(maximum, boundsMaximum);
        }
        if (!std::isfinite(minimum) || !std::isfinite(maximum) ||
            maximum - minimum <= m_lengthTolerance) {
            return;
        }

        for (const auto& interval :
             visibleIntervals(edge, curve, minimum, maximum)) {
            appendResultEdge(edge, curve, interval.first, interval.second,
                             result);
        }
    }

    std::vector<SweepOutputEdge> buildSweepEdges() {
        std::vector<SweepOutputEdge> result;
        processSweep();
        for (const EdgeRecord& edge : m_edges) {
            emitSweepEdge(edge, result);
        }
        return result;
    }

};

inline SegmentVoronoiDiagram SegmentVoronoiBuilder::build() {
    SegmentVoronoiDiagram result;
    result.cells.reserve(m_originalSites.size());
    for (const Segment& site : m_originalSites) {
        result.cells.push_back({site, {}});
    }
    if (m_originalSites.size() < 2U) return result;

    for (SweepOutputEdge& edge : buildSweepEdges()) {
        const std::size_t firstSite =
            m_features[edge.firstFeature].owners[0];
        const std::size_t secondSite =
            m_features[edge.secondFeature].owners[0];
        const std::size_t lower = std::min(firstSite, secondSite);
        const std::size_t upper = std::max(firstSite, secondSite);
        const std::size_t edgeIndex = result.edges.size();
        result.edges.push_back(
            {lower, upper, edge.curved, std::move(edge.vertices)});
        result.cells[lower].edges.push_back(edgeIndex);
        result.cells[upper].edges.push_back(edgeIndex);
    }
    return result;
}

inline PolygonSweepDiagram SegmentVoronoiBuilder::buildPolygonDiagram() {
    if (!m_polygonMode) {
        throw std::logic_error(
            "Polygon sweep output requested from a segment builder");
    }
    PolygonSweepDiagram result;
    for (SweepOutputEdge& edge : buildSweepEdges()) {
        const Feature& first = m_features[edge.firstFeature];
        const Feature& second = m_features[edge.secondFeature];
        result.edges.push_back(
            {{first.kind == FeatureKind::Point, first.inputIndex},
             {second.kind == FeatureKind::Point, second.inputIndex},
             edge.curved, std::move(edge.vertices)});
    }
    return result;
}

inline SegmentVoronoiBuilder::SegmentVoronoiBuilder(
    std::vector<Segment> sites, common::geometry::Aabb<double> bounds,
    double maximumError)
    : m_originalSites(std::move(sites)), m_inputBounds(bounds) {
    validateScalars(maximumError);
    if (m_originalSites.size() < 2U) return;
    initializeLocalCoordinates();
    chooseSweepBasis();
    validateDisjointSites();
    initializeSegmentFeaturesAndEvents();
    m_maximumError =
        maximumError == 0.0
            ? 1e-4L
            : static_cast<Real>(maximumError) / m_scale;
}

inline SegmentVoronoiBuilder::SegmentVoronoiBuilder(
    std::vector<Point> vertices, common::geometry::Aabb<double> bounds,
    double maximumError, PolygonSweepTag)
    : m_polygonVertices(std::move(vertices)),
      m_inputBounds(bounds), m_polygonMode(true) {
    if (m_polygonVertices.size() < 3U) {
        throw std::invalid_argument(
            "Polygon Voronoi input must have at least three vertices");
    }
    m_originalSites.reserve(m_polygonVertices.size());
    for (std::size_t index = 0; index < m_polygonVertices.size(); ++index) {
        m_originalSites.push_back(
            {m_polygonVertices[index],
             m_polygonVertices[(index + 1U) % m_polygonVertices.size()]});
    }

    validateScalars(maximumError);
    initializeLocalCoordinates();
    chooseSweepBasis();
    validateSimplePolygon();
    initializePolygonFeaturesAndEvents();
    m_maximumError =
        maximumError == 0.0
            ? 1e-4L
            : static_cast<Real>(maximumError) / m_scale;
}

} // namespace detail

template <typename coordinate_type, typename bounds_type>
SegmentVoronoiDiagram segmentVoronoiDiagram(
    const std::vector<BaseSegment<coordinate_type>>& sites,
    const common::geometry::Aabb<bounds_type>& bounds, double maximumError = 0.0) {
    std::vector<Segment> doubleSites;
    doubleSites.reserve(sites.size());
    for (const BaseSegment<coordinate_type>& site : sites) {
        doubleSites.push_back(
            {{detail::checkedSegmentVoronoiCoordinate(site.start.x()),
              detail::checkedSegmentVoronoiCoordinate(site.start.y())},
             {detail::checkedSegmentVoronoiCoordinate(site.end.x()),
              detail::checkedSegmentVoronoiCoordinate(site.end.y())}});
    }
    return detail::SegmentVoronoiBuilder(
               std::move(doubleSites),
               {detail::checkedSegmentVoronoiCoordinate(bounds.minX()),
                detail::checkedSegmentVoronoiCoordinate(bounds.minY()),
                detail::checkedSegmentVoronoiCoordinate(bounds.maxX()),
                detail::checkedSegmentVoronoiCoordinate(bounds.maxY())},
               maximumError)
        .build();
}

// Builds the bounded nearest-segment diagram with a generalized Fortune sweep.
// Each finite segment contributes two endpoint point features and two directed
// views of its interior feature. The interior is activated at the first
// endpoint in a deterministically rotated top-to-bottom sweep and its temporary
// beach-line separator is retired at the other endpoint. PPP, PPS, PSS, and SSS
// tangent events update an explicit AVL beach tree through a priority queue.
// The tree stores breakpoint order structurally and tracks each branch of a
// point/line wavefront intersection; no mutable geometric comparator is used.
// Circle events complete DCEL edge endpoints directly.  Finished EdgeRecords
// are then clipped to the requested AABB, and only then are parabolic edges
// tessellated. Only boundaries between different input segments are returned.
// The sweep and validation are O(n log n) for generic disjoint input; clipping
// and parabolic discretization are output-sensitive.
//
// Sites must be finite, nonzero, pairwise disjoint (shared endpoints are not
// accepted), and exactly representable as doubles. The floating-point
// implementation rejects inputs for which no stable sweep direction exists,
// and ambiguous simultaneous tangent events. A point/interior feature pair
// that lies on a different segment's supporting line has no contributing
// parabolic branch and is skipped without rejecting the valid input.
template <typename coordinate_type, typename bounds_type>
SegmentVoronoiDiagram voronoiDiagram(
    const std::vector<BaseSegment<coordinate_type>>& sites,
    const common::geometry::Aabb<bounds_type>& bounds, double maximumError = 0.0) {
    return segmentVoronoiDiagram(sites, bounds, maximumError);
}

} // namespace skeleton

#endif // SKELETON_SEGMENT_FORTUNE_VORONOI_H
