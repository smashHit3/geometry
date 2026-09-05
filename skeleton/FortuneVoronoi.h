#ifndef SKELETON_FORTUNE_VORONOI_H
#define SKELETON_FORTUNE_VORONOI_H

#include "common/geometry/Aabb.h"
#include "common/geometryutil/BasePointUtil.h"
#include "skeleton/PolygonFortuneVoronoi.h"
#include "skeleton/SegmentFortuneVoronoi.h"
#include "skeleton/Types.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace skeleton {

struct VoronoiCell {
    Point site;
    std::vector<Point> vertices;
};

struct VoronoiEdge {
    Point firstSite;
    Point secondSite;
    Point start;
    Point end;
};

struct VoronoiDiagram {
    std::vector<VoronoiEdge> edges;
    std::vector<VoronoiCell> cells;
};

namespace detail {

class FortuneVoronoiBuilder {
public:
    FortuneVoronoiBuilder(std::vector<Point> sites, common::geometry::Aabb<double> bounds)
        : m_sites(std::move(sites)), m_bounds(bounds) {
        validateInput();
        m_scale = coordinateScale();
        m_lengthEpsilon = m_scale * 1e-11;
        m_areaEpsilon = m_scale * m_scale * 1e-12;
        initializeSweepSites();
    }

    VoronoiDiagram build() {
        if (m_sites.empty()) return {};
        if (allSitesCollinear()) {
            buildCollinearEdges();
        } else if (Point center; allSitesCocircular(center)) {
            buildCocircularEdges(center);
        } else {
            initializeSiteEvents();
            processEvents();
        }

        VoronoiDiagram diagram;
        std::vector<ClippedEdge> clipped_edges;
        clipped_edges.reserve(m_edges.size());
        for (const EdgeRecord& edge : m_edges) {
            const auto clipped = clipEdge(edge);
            if (!clipped.has_value()) continue;

            const ClippedEdge& clipped_edge = *clipped;
            diagram.edges.push_back(
                {m_sites[edge.first], m_sites[edge.second],
                 clipped_edge.start, clipped_edge.end});
            clipped_edges.push_back(clipped_edge);
        }
        diagram.cells = buildCells(clipped_edges);
        return diagram;
    }

private:
    struct Event;

    struct EdgeVertex {
        Point point;
        std::size_t opposite;
    };

    struct EdgeRecord {
        std::size_t first;
        std::size_t second;
        std::vector<EdgeVertex> vertices;
        unsigned int completedBreakpoints = 0;
    };

    struct Arc {
        std::size_t site;
        Arc* previous = nullptr;
        Arc* next = nullptr;
        Arc* left = nullptr;
        Arc* right = nullptr;
        Arc* parent = nullptr;
        int height = 1;
        EdgeRecord* edgeToNext = nullptr;

        Event* circle = nullptr;
    };

    enum class EventKind { Site, Circle };

    struct Event {
        EventKind kind;
        double y;
        double x;
        std::size_t site = 0;
        Arc* arc = nullptr;
        Point center;
        bool valid = true;
        std::uint64_t sequence = 0;
    };

    struct EventBefore {
        bool operator()(const Event* left, const Event* right) const {
            if (left->y != right->y) return left->y < right->y;
            if (left->kind != right->kind) {
                return left->kind == EventKind::Circle;
            }
            if (left->x != right->x) return left->x > right->x;
            return left->sequence > right->sequence;
        }
    };

    struct ClippedEdge {
        std::size_t first;
        std::size_t second;
        Point start;
        Point end;
    };

    std::vector<Point> m_sites;
    std::vector<Point> m_sweepSites;
    common::geometry::Aabb<double> m_bounds;
    double m_scale = 1.0;
    double m_lengthEpsilon = 1e-11;
    double m_areaEpsilon = 1e-12;
    std::list<EdgeRecord> m_edges;
    std::map<std::uint64_t, EdgeRecord*> m_edgeIndexes;
    std::vector<std::unique_ptr<Arc>> m_arcs;
    std::vector<std::unique_ptr<Event>> m_events;
    std::priority_queue<Event*, std::vector<Event*>, EventBefore> m_queue;
    Arc* m_root = nullptr;
    Arc* m_first = nullptr;
    Arc* m_last = nullptr;
    std::uint64_t m_nextEventSequence = 0;

    void validateInput() const {
        if (m_bounds.minX() >= m_bounds.maxX() ||
            m_bounds.minY() >= m_bounds.maxY()) {
            throw std::invalid_argument("Voronoi bounds must have positive area");
        }
        if (!std::isfinite(m_bounds.minX()) || !std::isfinite(m_bounds.minY()) ||
            !std::isfinite(m_bounds.maxX()) || !std::isfinite(m_bounds.maxY())) {
            throw std::invalid_argument("Voronoi bounds must be finite");
        }

        std::vector<Point> sorted_sites = m_sites;
        for (const Point& site : sorted_sites) {
            if (!std::isfinite(site.x()) || !std::isfinite(site.y())) {
                throw std::invalid_argument("Voronoi sites must be finite");
            }
        }
        std::sort(sorted_sites.begin(), sorted_sites.end(),
                  [](const Point& left, const Point& right) {
                      if (left.x() != right.x()) return left.x() < right.x();
                      return left.y() < right.y();
                  });
        for (std::size_t i = 1; i < sorted_sites.size(); ++i) {
            if (sorted_sites[i - 1].x() == sorted_sites[i].x() &&
                sorted_sites[i - 1].y() == sorted_sites[i].y()) {
                throw std::invalid_argument("Voronoi sites must be distinct");
            }
        }
    }

    double coordinateScale() const {
        double scale = std::max(
            {1.0, std::abs(m_bounds.minX()), std::abs(m_bounds.minY()),
             std::abs(m_bounds.maxX()), std::abs(m_bounds.maxY())});
        for (const Point& site : m_sites) {
            scale = std::max(scale, std::abs(site.x()));
            scale = std::max(scale, std::abs(site.y()));
        }
        return scale;
    }

    void initializeSweepSites() {
        m_sweepSites = m_sites;
        std::vector<std::size_t> ordered(m_sites.size());
        for (std::size_t index = 0; index < ordered.size(); ++index) {
            ordered[index] = index;
        }
        std::sort(ordered.begin(), ordered.end(), [&](std::size_t left,
                                                       std::size_t right) {
            if (m_sites[left].y() != m_sites[right].y()) {
                return m_sites[left].y() > m_sites[right].y();
            }
            return m_sites[left].x() < m_sites[right].x();
        });

        for (std::size_t begin = 0; begin < ordered.size();) {
            std::size_t end = begin + 1;
            while (end < ordered.size() &&
                   m_sites[ordered[begin]].y() == m_sites[ordered[end]].y()) {
                ++end;
            }
            if (end - begin > 1) {
                double nearest_level_distance = std::numeric_limits<double>::infinity();
                if (begin > 0) {
                    nearest_level_distance =
                        m_sites[ordered[begin - 1]].y() -
                        m_sites[ordered[begin]].y();
                }
                if (end < ordered.size()) {
                    nearest_level_distance =
                        std::min(nearest_level_distance,
                                 m_sites[ordered[begin]].y() -
                                     m_sites[ordered[end]].y());
                }
                const double minimum_nudge =
                    32.0 * std::numeric_limits<double>::epsilon() *
                    std::max(1.0, std::abs(m_sites[ordered[begin]].y()));
                const double nudge =
                    std::max(nearest_level_distance * 1e-8, minimum_nudge);
                if (std::isfinite(nearest_level_distance) &&
                    nudge * static_cast<double>(end - begin) <
                        nearest_level_distance / 4.0) {
                    for (std::size_t index = begin; index < end; ++index) {
                        m_sweepSites[ordered[index]].setY(
                            m_sites[ordered[index]].y() +
                            nudge * static_cast<double>(end - index - 1));
                    }
                }
            }
            begin = end;
        }
    }

    static std::uint64_t edgeKey(std::size_t first, std::size_t second) {
        const auto lower = static_cast<std::uint64_t>(std::min(first, second));
        const auto upper = static_cast<std::uint64_t>(std::max(first, second));
        return (upper << 32U) | lower;
    }

    EdgeRecord* edgeFor(std::size_t first, std::size_t second) {
        const std::uint64_t key = edgeKey(first, second);
        const auto existing = m_edgeIndexes.find(key);
        if (existing != m_edgeIndexes.end()) return existing->second;

        m_edges.push_back({first, second, {}, 0});
        EdgeRecord* edge = &m_edges.back();
        m_edgeIndexes.emplace(key, edge);
        return edge;
    }

    Arc* makeArc(std::size_t site) {
        m_arcs.push_back(std::make_unique<Arc>());
        Arc* arc = m_arcs.back().get();
        arc->site = site;
        return arc;
    }

    Event* makeEvent(Event event) {
        event.sequence = m_nextEventSequence++;
        m_events.push_back(std::make_unique<Event>(std::move(event)));
        return m_events.back().get();
    }

    void initializeSiteEvents() {
        for (std::size_t index = 0; index < m_sites.size(); ++index) {
            m_queue.push(makeEvent(
                {EventKind::Site, m_sweepSites[index].y(), m_sweepSites[index].x(),
                 index}));
        }
    }

    void replaceParent(Arc* old_child, Arc* new_child) {
        Arc* parent = old_child->parent;
        if (parent == nullptr) {
            m_root = new_child;
        } else if (parent->left == old_child) {
            parent->left = new_child;
        } else {
            parent->right = new_child;
        }
        if (new_child != nullptr) new_child->parent = parent;
    }

    void rotateLeft(Arc* arc) {
        Arc* pivot = arc->right;
        arc->right = pivot->left;
        if (arc->right != nullptr) arc->right->parent = arc;
        replaceParent(arc, pivot);
        pivot->left = arc;
        arc->parent = pivot;
        updateHeight(arc);
        updateHeight(pivot);
    }

    void rotateRight(Arc* arc) {
        Arc* pivot = arc->left;
        arc->left = pivot->right;
        if (arc->left != nullptr) arc->left->parent = arc;
        replaceParent(arc, pivot);
        pivot->right = arc;
        arc->parent = pivot;
        updateHeight(arc);
        updateHeight(pivot);
    }

    static int height(const Arc* arc) {
        return arc == nullptr ? 0 : arc->height;
    }

    static int balanceFactor(const Arc* arc) {
        return height(arc->left) - height(arc->right);
    }

    static void updateHeight(Arc* arc) {
        if (arc != nullptr) {
            arc->height = 1 + std::max(height(arc->left), height(arc->right));
        }
    }

    void rebalanceFrom(Arc* arc) {
        while (arc != nullptr) {
            updateHeight(arc);
            if (balanceFactor(arc) > 1) {
                if (balanceFactor(arc->left) < 0) rotateLeft(arc->left);
                rotateRight(arc);
            } else if (balanceFactor(arc) < -1) {
                if (balanceFactor(arc->right) > 0) rotateRight(arc->right);
                rotateLeft(arc);
            }
            arc = arc->parent;
        }
    }

    void insertBetween(Arc* previous, Arc* next, Arc* arc) {
        if (m_root == nullptr) {
            m_root = arc;
            return;
        }

        Arc* parent = nullptr;
        if (previous != nullptr) {
            parent = previous;
            if (parent->right != nullptr) {
                parent = parent->right;
                while (parent->left != nullptr) parent = parent->left;
                parent->left = arc;
            } else {
                parent->right = arc;
            }
            arc->parent = parent;
        } else {
            parent = next;
            if (parent->left != nullptr) {
                parent = parent->left;
                while (parent->right != nullptr) parent = parent->right;
                parent->right = arc;
            } else {
                parent->left = arc;
            }
            arc->parent = parent;
        }
        rebalanceFrom(parent);
    }

    void removeFromTree(Arc* arc) {
        Arc* rebalance_start = nullptr;
        if (arc->left == nullptr) {
            rebalance_start = arc->parent;
            replaceParent(arc, arc->right);
            if (rebalance_start == nullptr) rebalance_start = arc->right;
        } else if (arc->right == nullptr) {
            rebalance_start = arc->parent;
            replaceParent(arc, arc->left);
            if (rebalance_start == nullptr) rebalance_start = arc->left;
        } else {
            Arc* successor = arc->right;
            while (successor->left != nullptr) successor = successor->left;
            if (successor->parent != arc) {
                rebalance_start = successor->parent;
                replaceParent(successor, successor->right);
                successor->right = arc->right;
                successor->right->parent = successor;
            } else {
                rebalance_start = successor;
            }
            replaceParent(arc, successor);
            successor->left = arc->left;
            successor->left->parent = successor;
            updateHeight(successor);
        }

        arc->left = nullptr;
        arc->right = nullptr;
        arc->parent = nullptr;
        arc->height = 1;
        rebalanceFrom(rebalance_start);
    }

    double breakpointX(std::size_t left_index, std::size_t right_index,
                       double directrix) const {
        const Point& left = m_sweepSites[left_index];
        const Point& right = m_sweepSites[right_index];
        if (left.y() == right.y()) return (left.x() + right.x()) / 2.0;
        if (left.y() == directrix) return left.x();
        if (right.y() == directrix) return right.x();

        const double left_denominator = 2.0 * (left.y() - directrix);
        const double right_denominator = 2.0 * (right.y() - directrix);
        const double a = 1.0 / left_denominator - 1.0 / right_denominator;
        const double b = -2.0 * (left.x() / left_denominator -
                                 right.x() / right_denominator);
        const double c =
            (left.x() * left.x() + left.y() * left.y() -
             directrix * directrix) /
                left_denominator -
            (right.x() * right.x() + right.y() * right.y() -
             directrix * directrix) /
                right_denominator;
        if (a == 0.0) return -c / b;

        double discriminant = b * b - 4.0 * a * c;
        if (discriminant < 0.0 && discriminant > -m_areaEpsilon) {
            discriminant = 0.0;
        }
        if (discriminant < 0.0) return (left.x() + right.x()) / 2.0;

        const double root = std::sqrt(discriminant);
        const double stable_term = -0.5 * (b + std::copysign(root, b));
        const double first_root = stable_term / a;
        const double second_root =
            stable_term == 0.0 ? -b / (2.0 * a) : c / stable_term;
        const double lower = std::min(first_root, second_root);
        const double upper = std::max(first_root, second_root);
        return left.y() > right.y() ? lower : upper;
    }

    Arc* findArcAbove(double x, double directrix) const {
        Arc* current = m_root;
        Arc* candidate = nullptr;
        while (current != nullptr) {
            if (current->next == nullptr ||
                x <= breakpointX(current->site, current->next->site, directrix)) {
                candidate = current;
                current = current->left;
            } else {
                current = current->right;
            }
        }
        return candidate == nullptr ? m_last : candidate;
    }

    void invalidateCircle(Arc* arc) {
        if (arc->circle != nullptr) {
            arc->circle->valid = false;
            arc->circle = nullptr;
        }
    }

    static double orientation(const Point& first, const Point& second,
                              const Point& third) {
        return common::geometryutil::cross(second - first, third - second);
    }

    bool circumcenter(const Point& first, const Point& second, const Point& third,
                      Point& center) const {
        const Point second_offset = second - first;
        const Point third_offset = third - first;
        const double denominator = 2.0 * common::geometryutil::cross(second_offset, third_offset);
        if (std::abs(denominator) <= m_areaEpsilon) return false;
        center = {
            first.x() +
                (common::geometryutil::squaredLength(second_offset) * third_offset.y() -
                 common::geometryutil::squaredLength(third_offset) * second_offset.y()) /
                    denominator,
            first.y() +
                (common::geometryutil::squaredLength(third_offset) * second_offset.x() -
                 common::geometryutil::squaredLength(second_offset) * third_offset.x()) /
                    denominator};
        return std::isfinite(center.x()) && std::isfinite(center.y());
    }

    void scheduleCircle(Arc* middle, double directrix) {
        if (middle == nullptr || middle->previous == nullptr ||
            middle->next == nullptr) {
            return;
        }
        const Point& first = m_sweepSites[middle->previous->site];
        const Point& second = m_sweepSites[middle->site];
        const Point& third = m_sweepSites[middle->next->site];
        if (orientation(first, second, third) >= -m_areaEpsilon) return; // why
        if (std::abs(orientation(m_sites[middle->previous->site],
                                 m_sites[middle->site],
                                 m_sites[middle->next->site])) <= m_areaEpsilon) {
            return;
        }

        Point center;
        if (!circumcenter(first, second, third, center)) return;
        const double event_y =
            center.y() - std::hypot(center.x() - first.x(), center.y() - first.y());
        if (!std::isfinite(event_y) || event_y >= directrix - m_lengthEpsilon) {
            return;
        }

        Event* event = makeEvent(
            {EventKind::Circle, event_y, center.x(), 0, middle, center});
        middle->circle = event;
        m_queue.push(event);
    }

    void addVertex(EdgeRecord* edge, const Point& point, std::size_t opposite) {
        ++edge->completedBreakpoints;
        for (const EdgeVertex& existing : edge->vertices) {
            if (common::geometryutil::squaredDistance(existing.point, point) <=
                m_lengthEpsilon * m_lengthEpsilon) {
                return;
            }
        }
        edge->vertices.push_back({point, opposite});
    }

    void processSite(Event* event) {
        if (m_root == nullptr) {
            Arc* arc = makeArc(event->site);
            m_root = m_first = m_last = arc;
            return;
        }

        Arc* split = findArcAbove(m_sweepSites[event->site].x(), event->y);
        invalidateCircle(split);

        Arc* previous = split->previous;
        Arc* next = split->next;
        EdgeRecord* split_edge = edgeFor(split->site, event->site);
        Arc* left = makeArc(split->site);
        Arc* middle = makeArc(event->site);
        Arc* right = makeArc(split->site);
        left->previous = previous;
        left->next = middle;
        left->edgeToNext = split_edge;
        middle->previous = left;
        middle->next = right;
        middle->edgeToNext = split_edge;
        right->previous = middle;
        right->next = next;
        right->edgeToNext = split->edgeToNext;
        if (previous != nullptr) {
            previous->next = left;
        } else {
            m_first = left;
        }
        if (next != nullptr) {
            next->previous = right;
        } else {
            m_last = right;
        }

        removeFromTree(split);
        insertBetween(previous, next, left);
        insertBetween(left, next, middle);
        insertBetween(middle, next, right);

        scheduleCircle(left, event->y);
        scheduleCircle(right, event->y);
    }

    void processCircle(Event* event) {
        Arc* middle = event->arc;
        if (!event->valid || middle == nullptr || middle->circle != event ||
            middle->previous == nullptr || middle->next == nullptr) {
            return;
        }

        Arc* first = middle->previous;
        Arc* third = middle->next;
        middle->circle = nullptr;
        invalidateCircle(first);
        invalidateCircle(third);
        Point center;
        if (!circumcenter(m_sites[first->site], m_sites[middle->site],
                          m_sites[third->site], center)) {
            return;
        }
        addVertex(first->edgeToNext, center, third->site);
        addVertex(middle->edgeToNext, center, first->site);
        EdgeRecord* new_edge = edgeFor(first->site, third->site);
        addVertex(new_edge, center, middle->site);

        first->next = third;
        first->edgeToNext = new_edge;
        third->previous = first;
        removeFromTree(middle);

        scheduleCircle(first, event->y);
        scheduleCircle(third, event->y);
    }

    void processEvents() {
        while (!m_queue.empty()) {
            Event* event = m_queue.top();
            m_queue.pop();
            if (event->kind == EventKind::Site) {
                processSite(event);
            } else {
                processCircle(event);
            }
        }
    }

    bool allSitesCollinear() const {
        if (m_sites.size() < 3) return true;
        const Point& first = m_sites[0];
        std::size_t second_index = 1;
        while (second_index < m_sites.size() &&
               common::geometryutil::squaredDistance(first, m_sites[second_index]) <=
                   m_lengthEpsilon * m_lengthEpsilon) {
            ++second_index;
        }
        if (second_index == m_sites.size()) return true;
        const Point direction = m_sites[second_index] - first;
        for (std::size_t index = second_index + 1; index < m_sites.size(); ++index) {
            if (std::abs(common::geometryutil::cross(direction, m_sites[index] - first)) >
                m_areaEpsilon) {
                return false;
            }
        }
        return true;
    }

    bool allSitesCocircular(Point& center) const {
        std::size_t second_index = 1;
        std::size_t third_index = 2;
        while (third_index < m_sites.size() &&
               std::abs(orientation(m_sites[0], m_sites[second_index],
                                    m_sites[third_index])) <= m_areaEpsilon) {
            ++third_index;
        }
        if (third_index == m_sites.size()) return false;

        if (!circumcenter(m_sites[0], m_sites[second_index],
                          m_sites[third_index], center)) {
            return false;
        }
        const double radius_squared = common::geometryutil::squaredDistance(m_sites[0], center);
        for (const Point& site : m_sites) {
            if (std::abs(common::geometryutil::squaredDistance(site, center) - radius_squared) >
                m_areaEpsilon) {
                return false;
            }
        }
        return true;
    }

    void buildCollinearEdges() {
        if (m_sites.size() < 2) return;
        const Point direction = m_sites[1] - m_sites[0];
        std::vector<std::size_t> ordered(m_sites.size());
        for (std::size_t index = 0; index < ordered.size(); ++index) {
            ordered[index] = index;
        }
        std::sort(ordered.begin(), ordered.end(), [&](std::size_t left,
                                                       std::size_t right) {
            return common::geometryutil::dot(m_sites[left], direction) <
                   common::geometryutil::dot(m_sites[right], direction);
        });
        for (std::size_t index = 1; index < ordered.size(); ++index) {
            edgeFor(ordered[index - 1], ordered[index]);
        }
    }

    void buildCocircularEdges(const Point& center) {
        std::vector<std::size_t> ordered(m_sites.size());
        for (std::size_t index = 0; index < ordered.size(); ++index) {
            ordered[index] = index;
        }
        std::sort(ordered.begin(), ordered.end(), [&](std::size_t left,
                                                       std::size_t right) {
            return std::atan2(m_sites[left].y() - center.y(),
                              m_sites[left].x() - center.x()) <
                   std::atan2(m_sites[right].y() - center.y(),
                              m_sites[right].x() - center.x());
        });
        for (std::size_t index = 0; index < ordered.size(); ++index) {
            EdgeRecord* edge =
                edgeFor(ordered[index], ordered[(index + 1) % ordered.size()]);
            addVertex(edge, center, ordered[(index + 2) % ordered.size()]);
        }
    }

    Point clampToBounds(Point point) const {
        point.setX(std::clamp(point.x(), m_bounds.minX(), m_bounds.maxX()));
        point.setY(std::clamp(point.y(), m_bounds.minY(), m_bounds.maxY()));
        return point;
    }

    bool clipParameterRange(const Point& origin, const Point& direction,
                            double minimum_parameter, double maximum_parameter,
                            Point& start, Point& end) const {
        const auto update_range = [&](double origin_coordinate,
                                      double direction_coordinate,
                                      double lower, double upper,
                                      double& minimum, double& maximum) {
            if (direction_coordinate == 0.0) {
                return origin_coordinate >= lower - m_lengthEpsilon &&
                       origin_coordinate <= upper + m_lengthEpsilon;
            }
            double first = (lower - origin_coordinate) / direction_coordinate;
            double second = (upper - origin_coordinate) / direction_coordinate;
            if (first > second) std::swap(first, second);
            minimum = std::max(minimum, first);
            maximum = std::min(maximum, second);
            return minimum <= maximum + m_lengthEpsilon;
        };

        if (!update_range(origin.x(), direction.x(), m_bounds.minX(), m_bounds.maxX(),
                          minimum_parameter, maximum_parameter) ||
            !update_range(origin.y(), direction.y(), m_bounds.minY(), m_bounds.maxY(),
                          minimum_parameter, maximum_parameter) ||
            minimum_parameter > maximum_parameter + m_lengthEpsilon) {
            return false;
        }
        start = clampToBounds(origin + direction * minimum_parameter);
        end = clampToBounds(origin + direction * maximum_parameter);
        return common::geometryutil::squaredDistance(start, end) >
               m_lengthEpsilon * m_lengthEpsilon;
    }

    std::optional<ClippedEdge> clipEdge(const EdgeRecord& edge) const {
        const Point& first_site = m_sites[edge.first];
        const Point& second_site = m_sites[edge.second];
        const Point direction{
            first_site.y() - second_site.y(), second_site.x() - first_site.x()};
        Point start;
        Point end;
        if (edge.completedBreakpoints >= 2U) {
            if (edge.vertices.size() < 2U) return std::nullopt;
            const Point segment_direction =
                edge.vertices[1].point - edge.vertices[0].point;
            if (!clipParameterRange(edge.vertices[0].point, segment_direction, 0.0,
                                    1.0, start, end)) {
                return std::nullopt;
            }
        } else if (edge.completedBreakpoints == 1U) {
            const EdgeVertex& vertex = edge.vertices.front();
            Point ray_direction = direction;
            const Point to_opposite = m_sites[vertex.opposite] - first_site;
            if (common::geometryutil::dot(to_opposite, ray_direction) > 0.0) {
                ray_direction = ray_direction * -1.0;
            }
            if (!clipParameterRange(vertex.point, ray_direction, 0.0,
                                    std::numeric_limits<double>::infinity(), start,
                                    end)) {
                return std::nullopt;
            }
        } else if (!clipParameterRange(
                       (first_site + second_site) / 2.0, direction,
                       -std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity(), start, end)) {
            return std::nullopt;
        }
        return ClippedEdge{edge.first, edge.second, start, end};
    }

    std::vector<VoronoiCell> buildCells(
        const std::vector<ClippedEdge>& edges) const {
        std::vector<std::vector<Point>> vertices(m_sites.size());
        for (const ClippedEdge& edge : edges) {
            vertices[edge.first].push_back(edge.start);
            vertices[edge.first].push_back(edge.end);
            vertices[edge.second].push_back(edge.start);
            vertices[edge.second].push_back(edge.end);
        }

        const std::vector<Point> corners{
            {m_bounds.minX(), m_bounds.minY()},
            {m_bounds.maxX(), m_bounds.minY()},
            {m_bounds.maxX(), m_bounds.maxY()},
            {m_bounds.minX(), m_bounds.maxY()},
        };
        for (const Point& corner : corners) {
            double closest_distance = std::numeric_limits<double>::infinity();
            for (const Point& site : m_sites) {
                closest_distance =
                    std::min(closest_distance, common::geometryutil::squaredDistance(site, corner));
            }
            for (std::size_t index = 0; index < m_sites.size(); ++index) {
                if (std::abs(common::geometryutil::squaredDistance(m_sites[index], corner) -
                             closest_distance) <= m_areaEpsilon) {
                    vertices[index].push_back(corner);
                }
            }
        }

        std::vector<VoronoiCell> cells;
        cells.reserve(m_sites.size());
        for (std::size_t index = 0; index < m_sites.size(); ++index) {
            std::vector<Point>& cell_vertices = vertices[index];
            std::sort(cell_vertices.begin(), cell_vertices.end(),
                      [&](const Point& left, const Point& right) {
                          const double left_angle = std::atan2(
                              left.y() - m_sites[index].y(),
                              left.x() - m_sites[index].x());
                          const double right_angle = std::atan2(
                              right.y() - m_sites[index].y(),
                              right.x() - m_sites[index].x());
                          if (left_angle != right_angle) {
                              return left_angle < right_angle;
                          }
                          return common::geometryutil::squaredDistance(left, m_sites[index]) <
                                 common::geometryutil::squaredDistance(right, m_sites[index]);
                      });
            cell_vertices.erase(
                std::unique(cell_vertices.begin(), cell_vertices.end(),
                            [&](const Point& left, const Point& right) {
                                return common::geometryutil::squaredDistance(left, right) <=
                                       m_lengthEpsilon * m_lengthEpsilon;
                            }),
                cell_vertices.end());
            cells.push_back({m_sites[index], std::move(cell_vertices)});
        }
        return cells;
    }
};

inline VoronoiDiagram fortuneVoronoiDiagram(const std::vector<Point>& sites,
                                            const common::geometry::Aabb<double>& bounds) {
    return FortuneVoronoiBuilder(sites, bounds).build();
}

} // namespace detail

// Builds the point-site Voronoi diagram using a Fortune sweep line. For finite,
// nondegenerate input, event processing and clipped-output construction are
// O(n log n). Returned edges and cell boundaries use doubles and are clipped to
// bounds. Collinear sites are handled separately; this uses scale-relative
// floating-point predicates, not exact-arithmetic degeneracy handling.
template <typename site_type, typename bounds_type>
VoronoiDiagram voronoiDiagram(
    const std::vector<common::geometry::BasePoint<site_type>>& sites,
    const common::geometry::Aabb<bounds_type>& bounds) {
    std::vector<Point> double_sites;
    double_sites.reserve(sites.size());
    for (const auto& site : sites) {
        double_sites.emplace_back(static_cast<double>(site.x()),
                                  static_cast<double>(site.y()));
    }
    return detail::fortuneVoronoiDiagram(
        double_sites,
        {static_cast<double>(bounds.minX()), static_cast<double>(bounds.minY()),
         static_cast<double>(bounds.maxX()), static_cast<double>(bounds.maxY())});
}

template <typename site_type, typename bounds_type>
std::vector<VoronoiCell> halfPlaneVoronoiCells(
    const std::vector<common::geometry::BasePoint<site_type>>& sites,
    const common::geometry::Aabb<bounds_type>& bounds) {
    return voronoiDiagram(sites, bounds).cells;
}

} // namespace skeleton

#endif // SKELETON_FORTUNE_VORONOI_H
