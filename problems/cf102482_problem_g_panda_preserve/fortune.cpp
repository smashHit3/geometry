#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using Real = long double;

constexpr Real EPSILON = 1e-9L;

class Point {
public:
    constexpr Point() = default;
    constexpr Point(Real x, Real y) : x_(x), y_(y) {}

    constexpr Real& x() { return x_; }
    constexpr const Real& x() const { return x_; }
    constexpr Real& y() { return y_; }
    constexpr const Real& y() const { return y_; }
    constexpr void setX(Real x) { x_ = x; }
    constexpr void setY(Real y) { y_ = y; }

private:
    Real x_ = 0.0;
    Real y_ = 0.0;
};

constexpr Point operator+(const Point& left, const Point& right)
{
    return {left.x() + right.x(), left.y() + right.y()};
}

constexpr Point operator-(const Point& left, const Point& right)
{
    return {left.x() - right.x(), left.y() - right.y()};
}

constexpr Point operator*(const Point& point, Real scalar)
{
    return {point.x() * scalar, point.y() * scalar};
}

constexpr Point operator/(const Point& point, Real scalar)
{
    return {point.x() / scalar, point.y() / scalar};
}

constexpr Real dot(const Point& left, const Point& right)
{
    return left.x() * right.x() + left.y() * right.y();
}

constexpr Real vectorCross(const Point& left, const Point& right)
{
    return left.x() * right.y() - left.y() * right.x();
}

constexpr Real squaredLength(const Point& point)
{
    return dot(point, point);
}

constexpr Real squaredDistance(const Point& first, const Point& second)
{
    return squaredLength(first - second);
}

class Aabb {
public:
    constexpr Aabb(Real min_x, Real min_y, Real max_x, Real max_y)
        : min_x_(min_x), min_y_(min_y), max_x_(max_x), max_y_(max_y)
    {
    }

    constexpr Real minX() const { return min_x_; }
    constexpr Real minY() const { return min_y_; }
    constexpr Real maxX() const { return max_x_; }
    constexpr Real maxY() const { return max_y_; }

private:
    Real min_x_;
    Real min_y_;
    Real max_x_;
    Real max_y_;
};

struct Segment {
    Point start;
    Point end;
};

struct VoronoiEdge {
    std::size_t firstSite;
    std::size_t secondSite;
    Segment segment;
    std::vector<Point> vertices;
};

class FortuneVoronoiBuilder {
public:
    FortuneVoronoiBuilder(std::vector<Point> sites, Aabb bounds)
        : sites_(std::move(sites)), bounds_(bounds)
    {
        validateInput();
        const Real scale = coordinateScale();
        length_epsilon_ = scale * 1e-11;
        area_epsilon_ = scale * scale * 1e-12;
    }

    std::vector<VoronoiEdge> build()
    {
        if (sites_.empty()) {
            return {};
        }
        initializeSiteEvents();
        processEvents();

        std::vector<VoronoiEdge> result;
        result.reserve(edges_.size());
        for (const EdgeRecord& edge : edges_) {
            const auto clipped = clipEdge(edge);
            if (clipped.has_value()) {
                std::vector<Point> vertices;
                vertices.reserve(edge.vertices.size());
                for (const EdgeVertex& vertex : edge.vertices) {
                    vertices.push_back(vertex.point);
                }
                result.push_back({
                    edge.first, edge.second, *clipped,
                    std::move(vertices)});
            }
        }
        return result;
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
        Real y;
        Real x;
        std::size_t site = 0;
        Arc* arc = nullptr;
        bool valid = true;
        std::uint64_t sequence = 0;
    };

    struct EventBefore {
        bool operator()(const Event* left, const Event* right) const
        {
            if (left->y != right->y) {
                return left->y < right->y;
            }
            if (left->kind != right->kind) {
                return left->kind == EventKind::Circle;
            }
            if (left->x != right->x) {
                return left->x > right->x;
            }
            return left->sequence > right->sequence;
        }
    };

    std::vector<Point> sites_;
    Aabb bounds_;
    Real length_epsilon_ = 1e-11;
    Real area_epsilon_ = 1e-12;
    std::list<EdgeRecord> edges_;
    std::map<std::uint64_t, EdgeRecord*> edge_indexes_;
    std::vector<std::unique_ptr<Arc>> arcs_;
    std::vector<std::unique_ptr<Event>> events_;
    std::priority_queue<Event*, std::vector<Event*>, EventBefore> queue_;
    Arc* root_ = nullptr;
    Arc* last_ = nullptr;
    std::uint64_t next_event_sequence_ = 0;

    void validateInput() const
    {
        if (bounds_.minX() >= bounds_.maxX() ||
            bounds_.minY() >= bounds_.maxY()) {
            throw std::invalid_argument("Voronoi bounds must have positive area");
        }
        std::vector<Point> sorted_sites = sites_;
        std::sort(sorted_sites.begin(), sorted_sites.end(),
                  [](const Point& left, const Point& right) {
                      if (left.x() != right.x()) {
                          return left.x() < right.x();
                      }
                      return left.y() < right.y();
                  });
        for (std::size_t index = 1; index < sorted_sites.size(); ++index) {
            if (sorted_sites[index - 1].x() == sorted_sites[index].x() &&
                sorted_sites[index - 1].y() == sorted_sites[index].y()) {
                throw std::invalid_argument("Voronoi sites must be distinct");
            }
        }
    }

    Real coordinateScale() const
    {
        Real scale =
            std::max({1.0L, std::abs(bounds_.minX()), std::abs(bounds_.minY()),
                      std::abs(bounds_.maxX()), std::abs(bounds_.maxY())});
        for (const Point& site : sites_) {
            scale = std::max(scale, std::abs(site.x()));
            scale = std::max(scale, std::abs(site.y()));
        }
        return scale;
    }

    static std::uint64_t edgeKey(std::size_t first, std::size_t second)
    {
        const auto lower =
            static_cast<std::uint64_t>(std::min(first, second));
        const auto upper =
            static_cast<std::uint64_t>(std::max(first, second));
        return (upper << 32U) | lower;
    }

    EdgeRecord* edgeFor(std::size_t first, std::size_t second)
    {
        const std::uint64_t key = edgeKey(first, second);
        const auto existing = edge_indexes_.find(key);
        if (existing != edge_indexes_.end()) {
            return existing->second;
        }
        edges_.push_back({first, second, {}, 0});
        EdgeRecord* edge = &edges_.back();
        edge_indexes_.emplace(key, edge);
        return edge;
    }

    Arc* makeArc(std::size_t site)
    {
        arcs_.push_back(std::make_unique<Arc>());
        Arc* arc = arcs_.back().get();
        arc->site = site;
        return arc;
    }

    Event* makeEvent(Event event)
    {
        event.sequence = next_event_sequence_++;
        events_.push_back(std::make_unique<Event>(std::move(event)));
        return events_.back().get();
    }

    void initializeSiteEvents()
    {
        for (std::size_t index = 0; index < sites_.size(); ++index) {
            queue_.push(makeEvent(
                {EventKind::Site, sites_[index].y(),
                 sites_[index].x(), index, nullptr, true, 0}));
        }
    }

    void replaceParent(Arc* old_child, Arc* new_child)
    {
        Arc* parent = old_child->parent;
        if (parent == nullptr) {
            root_ = new_child;
        } else if (parent->left == old_child) {
            parent->left = new_child;
        } else {
            parent->right = new_child;
        }
        if (new_child != nullptr) {
            new_child->parent = parent;
        }
    }

    static int height(const Arc* arc)
    {
        return arc == nullptr ? 0 : arc->height;
    }

    static int balanceFactor(const Arc* arc)
    {
        return height(arc->left) - height(arc->right);
    }

    static void updateHeight(Arc* arc)
    {
        if (arc != nullptr) {
            arc->height = 1 + std::max(height(arc->left), height(arc->right));
        }
    }

    void rotateLeft(Arc* arc)
    {
        Arc* pivot = arc->right;
        arc->right = pivot->left;
        if (arc->right != nullptr) {
            arc->right->parent = arc;
        }
        replaceParent(arc, pivot);
        pivot->left = arc;
        arc->parent = pivot;
        updateHeight(arc);
        updateHeight(pivot);
    }

    void rotateRight(Arc* arc)
    {
        Arc* pivot = arc->left;
        arc->left = pivot->right;
        if (arc->left != nullptr) {
            arc->left->parent = arc;
        }
        replaceParent(arc, pivot);
        pivot->right = arc;
        arc->parent = pivot;
        updateHeight(arc);
        updateHeight(pivot);
    }

    void rebalanceFrom(Arc* arc)
    {
        while (arc != nullptr) {
            updateHeight(arc);
            if (balanceFactor(arc) > 1) {
                if (balanceFactor(arc->left) < 0) {
                    rotateLeft(arc->left);
                }
                rotateRight(arc);
            } else if (balanceFactor(arc) < -1) {
                if (balanceFactor(arc->right) > 0) {
                    rotateRight(arc->right);
                }
                rotateLeft(arc);
            }
            arc = arc->parent;
        }
    }

    void insertBetween(Arc* previous, Arc* next, Arc* arc)
    {
        if (root_ == nullptr) {
            root_ = arc;
            return;
        }

        Arc* parent = nullptr;
        if (previous != nullptr) {
            parent = previous;
            if (parent->right != nullptr) {
                parent = parent->right;
                while (parent->left != nullptr) {
                    parent = parent->left;
                }
                parent->left = arc;
            } else {
                parent->right = arc;
            }
        } else {
            parent = next;
            if (parent->left != nullptr) {
                parent = parent->left;
                while (parent->right != nullptr) {
                    parent = parent->right;
                }
                parent->right = arc;
            } else {
                parent->left = arc;
            }
        }
        arc->parent = parent;
        rebalanceFrom(parent);
    }

    void removeFromTree(Arc* arc)
    {
        Arc* rebalance_start = nullptr;
        if (arc->left == nullptr) {
            rebalance_start = arc->parent;
            replaceParent(arc, arc->right);
            if (rebalance_start == nullptr) {
                rebalance_start = arc->right;
            }
        } else if (arc->right == nullptr) {
            rebalance_start = arc->parent;
            replaceParent(arc, arc->left);
            if (rebalance_start == nullptr) {
                rebalance_start = arc->left;
            }
        } else {
            Arc* successor = arc->right;
            while (successor->left != nullptr) {
                successor = successor->left;
            }
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

    Real breakpointX(std::size_t left_index, std::size_t right_index,
                       Real directrix) const
    {
        const Point& left = sites_[left_index];
        const Point& right = sites_[right_index];
        if (left.y() == right.y()) {
            return (left.x() + right.x()) / 2.0;
        }
        if (left.y() == directrix) {
            return left.x();
        }
        if (right.y() == directrix) {
            return right.x();
        }

        const Real left_denominator = 2.0 * (left.y() - directrix);
        const Real right_denominator = 2.0 * (right.y() - directrix);
        const Real a =
            1.0 / left_denominator - 1.0 / right_denominator;
        const Real b =
            -2.0 * (left.x() / left_denominator -
                    right.x() / right_denominator);
        const Real c =
            (left.x() * left.x() + left.y() * left.y() -
             directrix * directrix) /
                left_denominator -
            (right.x() * right.x() + right.y() * right.y() -
             directrix * directrix) /
                right_denominator;
        if (a == 0.0) {
            return -c / b;
        }

        Real discriminant = b * b - 4.0 * a * c;
        if (discriminant < 0.0 && discriminant > -area_epsilon_) {
            discriminant = 0.0;
        }
        if (discriminant < 0.0) {
            return (left.x() + right.x()) / 2.0;
        }
        const Real root = std::sqrt(discriminant);
        const Real stable_term = -0.5 * (b + std::copysign(root, b));
        const Real first_root = stable_term / a;
        const Real second_root =
            stable_term == 0.0 ? -b / (2.0 * a) : c / stable_term;
        const Real lower = std::min(first_root, second_root);
        const Real upper = std::max(first_root, second_root);
        return left.y() > right.y() ? lower : upper;
    }

    Arc* findArcAbove(Real x, Real directrix) const
    {
        Arc* current = root_;
        Arc* candidate = nullptr;
        while (current != nullptr) {
            if (current->next == nullptr ||
                x <= breakpointX(current->site, current->next->site,
                                 directrix)) {
                candidate = current;
                current = current->left;
            } else {
                current = current->right;
            }
        }
        return candidate == nullptr ? last_ : candidate;
    }

    void invalidateCircle(Arc* arc)
    {
        if (arc->circle != nullptr) {
            arc->circle->valid = false;
            arc->circle = nullptr;
        }
    }

    static Real orientation(const Point& first, const Point& second,
                              const Point& third)
    {
        return vectorCross(second - first, third - second);
    }

    bool circumcenter(const Point& first, const Point& second,
                      const Point& third, Point& center) const
    {
        const Point second_offset = second - first;
        const Point third_offset = third - first;
        const Real denominator =
            2.0 * vectorCross(second_offset, third_offset);
        if (std::abs(denominator) <= area_epsilon_) {
            return false;
        }
        center = {
            first.x() +
                (squaredLength(second_offset) * third_offset.y() -
                 squaredLength(third_offset) * second_offset.y()) /
                    denominator,
            first.y() +
                (squaredLength(third_offset) * second_offset.x() -
                 squaredLength(second_offset) * third_offset.x()) /
                    denominator};
        return std::isfinite(center.x()) && std::isfinite(center.y());
    }

    void scheduleCircle(Arc* middle, Real directrix)
    {
        if (middle == nullptr || middle->previous == nullptr ||
            middle->next == nullptr) {
            return;
        }
        const Point& first = sites_[middle->previous->site];
        const Point& second = sites_[middle->site];
        const Point& third = sites_[middle->next->site];
        if (orientation(first, second, third) >= -area_epsilon_) {
            return;
        }
        if (std::abs(orientation(sites_[middle->previous->site],
                                 sites_[middle->site],
                                 sites_[middle->next->site])) <=
            area_epsilon_) {
            return;
        }

        Point center;
        if (!circumcenter(first, second, third, center)) {
            return;
        }
        Real event_y =
            center.y() -
            std::hypot(center.x() - first.x(), center.y() - first.y());
        if (!std::isfinite(event_y) ||
            event_y > directrix + length_epsilon_) {
            return;
        }

        Event* event = makeEvent(
            {EventKind::Circle, event_y, center.x(), 0, middle, true, 0});
        middle->circle = event;
        queue_.push(event);
    }

    void addVertex(EdgeRecord* edge, const Point& point,
                   std::size_t opposite)
    {
        ++edge->completedBreakpoints;
        for (const EdgeVertex& existing : edge->vertices) {
            if (squaredDistance(existing.point, point) <=
                length_epsilon_ * length_epsilon_) {
                return;
            }
        }
        edge->vertices.push_back({point, opposite});
    }

    void processSite(Event* event)
    {
        if (root_ == nullptr) {
            Arc* arc = makeArc(event->site);
            root_ = last_ = arc;
            return;
        }

        Arc* split = findArcAbove(sites_[event->site].x(), event->y);
        invalidateCircle(split);

        if (sites_[split->site].y() == event->y &&
            sites_[split->site].x() < sites_[event->site].x()) {
            Arc* next = split->next;
            if (next != nullptr) {
                invalidateCircle(next);
            }

            Arc* inserted = makeArc(event->site);
            EdgeRecord* edge = edgeFor(split->site, event->site);
            inserted->previous = split;
            inserted->next = next;
            inserted->edgeToNext =
                next == nullptr
                    ? nullptr
                    : edgeFor(inserted->site, next->site);
            split->next = inserted;
            split->edgeToNext = edge;
            if (next != nullptr) {
                next->previous = inserted;
            } else {
                last_ = inserted;
            }
            insertBetween(split, next, inserted);

            scheduleCircle(split, event->y);
            scheduleCircle(inserted, event->y);
            scheduleCircle(next, event->y);
            return;
        }

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
        }
        if (next != nullptr) {
            next->previous = right;
        } else {
            last_ = right;
        }

        removeFromTree(split);
        insertBetween(previous, next, left);
        insertBetween(left, next, middle);
        insertBetween(middle, next, right);

        scheduleCircle(left, event->y);
        scheduleCircle(right, event->y);
    }

    void processCircle(Event* event)
    {
        Arc* middle = event->arc;
        if (!event->valid || middle == nullptr ||
            middle->circle != event || middle->previous == nullptr ||
            middle->next == nullptr) {
            return;
        }

        Arc* first = middle->previous;
        Arc* third = middle->next;
        middle->circle = nullptr;
        invalidateCircle(first);
        invalidateCircle(third);
        Point center;
        if (!circumcenter(sites_[first->site], sites_[middle->site],
                          sites_[third->site], center)) {
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

    void processEvents()
    {
        while (!queue_.empty()) {
            Event* event = queue_.top();
            queue_.pop();
            if (event->kind == EventKind::Site) {
                processSite(event);
            } else {
                processCircle(event);
            }
        }
    }

    Point clampToBounds(Point point) const
    {
        point.setX(
            std::clamp(point.x(), bounds_.minX(), bounds_.maxX()));
        point.setY(
            std::clamp(point.y(), bounds_.minY(), bounds_.maxY()));
        return point;
    }

    bool clipParameterRange(const Point& origin, const Point& direction,
                            Real minimum_parameter,
                            Real maximum_parameter, Point& start,
                            Point& end) const
    {
        const auto update_range =
            [&](Real origin_coordinate, Real direction_coordinate,
                Real lower, Real upper, Real& minimum,
                Real& maximum) {
                if (direction_coordinate == 0.0) {
                    return origin_coordinate >= lower - length_epsilon_ &&
                           origin_coordinate <= upper + length_epsilon_;
                }
                Real first =
                    (lower - origin_coordinate) / direction_coordinate;
                Real second =
                    (upper - origin_coordinate) / direction_coordinate;
                if (first > second) {
                    std::swap(first, second);
                }
                minimum = std::max(minimum, first);
                maximum = std::min(maximum, second);
                return minimum <= maximum + length_epsilon_;
            };

        if (!update_range(origin.x(), direction.x(), bounds_.minX(),
                          bounds_.maxX(), minimum_parameter,
                          maximum_parameter) ||
            !update_range(origin.y(), direction.y(), bounds_.minY(),
                          bounds_.maxY(), minimum_parameter,
                          maximum_parameter) ||
            minimum_parameter > maximum_parameter + length_epsilon_) {
            return false;
        }
        start = clampToBounds(origin + direction * minimum_parameter);
        end = clampToBounds(origin + direction * maximum_parameter);
        return squaredDistance(start, end) >
               length_epsilon_ * length_epsilon_;
    }

    std::optional<Segment> clipEdge(const EdgeRecord& edge) const
    {
        const Point& first_site = sites_[edge.first];
        const Point& second_site = sites_[edge.second];
        const Point direction{first_site.y() - second_site.y(),
                              second_site.x() - first_site.x()};
        Point start;
        Point end;
        if (edge.completedBreakpoints >= 2U) {
            if (edge.vertices.size() < 2U) {
                return std::nullopt;
            }
            const Point segment_direction =
                edge.vertices[1].point - edge.vertices[0].point;
            if (!clipParameterRange(edge.vertices[0].point,
                                    segment_direction, 0.0, 1.0,
                                    start, end)) {
                return std::nullopt;
            }
        } else if (edge.completedBreakpoints == 1U) {
            const EdgeVertex& vertex = edge.vertices.front();
            Point ray_direction = direction;
            const Point to_opposite =
                sites_[vertex.opposite] - first_site;
            if (dot(to_opposite, ray_direction) > 0.0) {
                ray_direction = ray_direction * -1.0;
            }
            if (!clipParameterRange(
                    vertex.point, ray_direction, 0.0,
                    std::numeric_limits<Real>::infinity(), start, end)) {
                return std::nullopt;
            }
        } else {
            return std::nullopt;
        }
        return Segment{start, end};
    }
};

std::optional<Point> segmentIntersection(const Segment& first,
                                         const Segment& second)
{
    const Point first_direction = first.end - first.start;
    const Point second_direction = second.end - second.start;
    const Real denominator =
        vectorCross(first_direction, second_direction);
    if (std::abs(denominator) <= EPSILON) {
        return std::nullopt;
    }

    const Point offset = second.start - first.start;
    const Real first_parameter =
        vectorCross(offset, second_direction) / denominator;
    const Real second_parameter =
        vectorCross(offset, first_direction) / denominator;
    if (first_parameter < -EPSILON ||
        first_parameter > 1.0L + EPSILON ||
        second_parameter < -EPSILON ||
        second_parameter > 1.0L + EPSILON) {
        return std::nullopt;
    }
    return first.start + first_direction * first_parameter;
}

bool pointOnSegment(const Point& point, const Segment& segment)
{
    const Point direction = segment.end - segment.start;
    if (std::abs(vectorCross(direction, point - segment.start)) >
        EPSILON) {
        return false;
    }
    return dot(point - segment.start, point - segment.end) <= EPSILON;
}

bool pointInPolygon(const Point& point,
                    const std::vector<Point>& polygon)
{
    bool inside = false;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const Segment edge{
            polygon[index],
            polygon[(index + 1) % polygon.size()]};
        if (pointOnSegment(point, edge)) {
            return true;
        }
        const Point& first = edge.start;
        const Point& second = edge.end;
        if ((first.y() > point.y()) != (second.y() > point.y())) {
            const Real crossing_x =
                first.x() +
                (second.x() - first.x()) *
                    (point.y() - first.y()) /
                    (second.y() - first.y());
            if (crossing_x > point.x()) {
                inside = !inside;
            }
        }
    }
    return inside;
}

std::vector<VoronoiEdge> buildFortuneEdges(
    const std::vector<Point>& sites)
{
    Real minimum_x = std::numeric_limits<Real>::infinity();
    Real minimum_y = std::numeric_limits<Real>::infinity();
    Real maximum_x = -std::numeric_limits<Real>::infinity();
    Real maximum_y = -std::numeric_limits<Real>::infinity();
    for (const Point& site : sites) {
        minimum_x = std::min(minimum_x, site.x());
        minimum_y = std::min(minimum_y, site.y());
        maximum_x = std::max(maximum_x, site.x());
        maximum_y = std::max(maximum_y, site.y());
    }
    return FortuneVoronoiBuilder(
               sites, Aabb{minimum_x, minimum_y, maximum_x, maximum_y})
        .build();
}

Real clipPolygon(const std::vector<VoronoiEdge>& voronoi_edges,
                 const std::vector<Point>& polygon)
{
    Real answer_squared = 0.0L;
    for (const VoronoiEdge& edge : voronoi_edges) {
        for (const Point& vertex : edge.vertices) {
            if (pointInPolygon(vertex, polygon)) {
                answer_squared = std::max(
                    answer_squared,
                    squaredDistance(vertex, polygon[edge.firstSite]));
            }
        }

        for (std::size_t index = 0; index < polygon.size(); ++index) {
            const Segment boundary{
                polygon[index],
                polygon[(index + 1) % polygon.size()]};
            const auto intersection =
                segmentIntersection(edge.segment, boundary);
            if (intersection.has_value()) {
                answer_squared = std::max(
                    answer_squared,
                    squaredDistance(*intersection,
                                    polygon[edge.firstSite]));
            }
        }
    }
    return answer_squared;
}

} // namespace

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int vertex_count;
    std::cin >> vertex_count;

    std::vector<Point> polygon(vertex_count);
    for (Point& point : polygon) {
        std::cin >> point.x() >> point.y();
    }

    const std::vector<VoronoiEdge> voronoi_edges =
        buildFortuneEdges(polygon);

    const Real maximum_squared_radius =
        clipPolygon(voronoi_edges, polygon);

    std::cout << std::setprecision(12)
              << std::sqrt(maximum_squared_radius) << '\n';
    return 0;
}
