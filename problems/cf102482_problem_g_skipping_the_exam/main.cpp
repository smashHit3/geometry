#include <bits/stdc++.h>

constexpr double EPSILON = 1e-9;

class VoronoiEngine {
public:
    struct Point {
        double x, y;
        Point operator-(const Point& other) const 
        {
            return {x - other.x, y - other.y};
        }

        Point operator*(double scalar) const 
        {
            return {x * scalar, y * scalar};
        }

        Point operator+(const Point& other) const 
        {
            return {x + other.x, y + other.y};
        }
    };

    struct Aabb {
        Point min, max;
    };

    struct VoronoiCell {
        Point site; // The site point that generates this cell
        std::vector<Point> vertices; // The vertices of the Voronoi cell in counter-clockwise order
    };

    struct VoronoiEdge {
        Point first_site, second_site; // The two sites that generate this edge
        Point start, end; // The start and end points of the edge (may be at infinity)
    };

    struct VoronoiDiagram {
        std::vector<VoronoiCell> cells;
        std::vector<VoronoiEdge> edges;
    };

    VoronoiEngine() = default;
    VoronoiEngine(std::vector<Point> sites, Aabb bounding_box) :
        m_sites(std::move(sites)), m_bounding_box(std::move(bounding_box))
    {
        initializeSweepSites();
    }

    VoronoiDiagram computeVoronoiDiagram();

private:
    struct Event;

    struct EdgeVertex {
        Point point; // The point of the edge vertex
        std::size_t opposite; // The index of the opposite vertex in the Voronoi diagram
    };

    struct EdgeRecord {
        std::size_t first_site_index; // The index of the first site in the Voronoi diagram
        std::size_t second_site_index; // The index of the second site in the Voronoi diagram
        std::vector<EdgeVertex> vertices; // The vertices of the edge in counter-clockwise order
        unsigned int completed_break_points = 0; // The number of completed break points for this edge
    };

    struct Arc {
        std::size_t site_index; // The index of the site that generates this arc
        Arc* prev = nullptr; // Pointer to the previous arc in the beach line
        Arc* next = nullptr; // Pointer to the next arc in the beach line
        Arc* left_child = nullptr; // Pointer to the left child in the binary search tree
        Arc* right_child = nullptr; // Pointer to the right child in the binary search tree
        Arc* parent = nullptr; // Pointer to the parent in the binary search tree
        int height = 1; // The height of the arc in the binary search tree
        EdgeRecord* edge_to_next = nullptr; // Pointer to the edge record that connects this arc to the next arc

        Event* circle_event = nullptr; // Pointer to the circle event associated with this arc (if any)
    };
    
    enum class EventType {
        SITE_EVENT,
        CIRCLE_EVENT
    };

    struct Event {
        EventType type; // The type of the event (site or circle)
        double y; // The y-coordinate of the event
        double x; // The x-coordinate of the event
        std::size_t site_index; // The index of the site associated with this event (for site events)
        Arc* arc; // Pointer to the arc associated with this event (for circle events)
        Point circle_center; // The center of the circle associated with this event (for circle events)
        bool valid; // Whether the event is still valid (for circle events)
        std::uint64_t id; // A unique identifier for the event (for circle events)
    };

    struct EventComparator {
        bool operator()(const Event* lhs, const Event* rhs) const 
        {
            if (lhs->y != rhs->y) {
                return lhs->y < rhs->y; // Sort by y-coordinate
            }
            if (lhs->type != rhs->type) {
                return lhs->type == EventType::CIRCLE_EVENT; // Circle events have higher priority than site events
            }
            if (lhs->x != rhs->x) {
                return lhs->x > rhs->x; // Process equal-height sites left-to-right
            }
            return lhs->id > rhs->id; // Process otherwise equal events in insertion order
        }
    };

    struct ClippedEdge {
        std::size_t first_site_index; // The index of the first site in the Voronoi diagram
        std::size_t second_site_index; // The index of the second site in the Voronoi diagram
        Point start; // The start point of the clipped edge
        Point end; // The end point of the clipped edge
    };

    double dot(const Point& a, const Point& b) const 
    {
        return a.x * b.x + a.y * b.y;
    }

    double cross(const Point& a, const Point& b) const 
    {
        return a.x * b.y - a.y * b.x;
    }

    bool allSitesCollinear() const;
    void computeCollinearEdges();
    bool allSitesCoCircular(Point& center) const;
    void computeCoCircularEdges(const Point& center);
    bool circumCircle(const Point& a, const Point& b, const Point& c, Point& center) const;
    void initializeSweepSites();
    void initializeEventQueue();
    void processEventQueue();
    void processSiteEvent(Event* event);
    void processCircleEvent(Event* event);
    void scheduleCircleEvent(Arc* arc, double sweep_line_y);
    void invalidateCircleEvent(Arc* arc);

    std::optional<ClippedEdge> clipEdge(const EdgeRecord& edge) const;
    bool clipParameterRange(const Point& origin, const Point& direction, double min_param, double max_param, Point& start, Point& end) const;
    Point clampToBounds(const Point& point) const;
    std::vector<VoronoiCell> buildCells(const std::vector<ClippedEdge>& clipped_edges);

    std::pair<std::size_t, std::size_t> makeOrderedPair(std::size_t first, std::size_t second) const;
    Event* makeEvent(EventType type, const Point& point, std::size_t site_index);
    EdgeRecord* findEdgeRecord(std::size_t first_site_index, std::size_t second_site_index);
    void addEdgeVertex(EdgeRecord* edge_record, const Point& vertex_point, std::size_t opposite_vertex_index);

    Arc* makeArc(std::size_t site_index);
    void rotateLeft(Arc* node);
    void rotateRight(Arc* node);
    void replaceArc(Arc* old_arc, Arc* new_arc);
    void rebalanceFrom(Arc* node);
    void insertBetween(Arc* left_arc, Arc* right_arc, Arc* new_arc);
    void removeArc(Arc* arc);
    double breakPointX(std::size_t left_site_index, std::size_t right_site_index, double sweep_line_y) const;
    Arc* findArcAbove(double x, double sweep_line_y) const;

    static int getHeight(Arc* node) 
    {
        return node ? node->height : 0;
    }

    static int getBalanceFactor(Arc* node) 
    {
        return node ? getHeight(node->left_child) - getHeight(node->right_child) : 0;
    }

    static void updateHeight(Arc* node) 
    {
        if (node) {
            node->height = 1 + std::max(getHeight(node->left_child), getHeight(node->right_child));
        }
    }
    
    std::vector<Point> m_sites; // The input sites for the Voronoi diagram
    std::vector<Point> m_sweep_sites; // Slightly perturbed sites used only by the sweep
    Aabb m_bounding_box; // The bounding box for the Voronoi diagram
    std::list<EdgeRecord> m_edge_records; // The list of edge records for the Voronoi diagram
    std::map<std::pair<std::size_t, std::size_t>, EdgeRecord*> m_edge_map; // A map from site index pairs to edge records
    std::vector<std::unique_ptr<Arc>> m_arcs; // The arcs in the beach line
    std::vector<std::unique_ptr<Event>> m_events; // The events in the event queue
    std::priority_queue<Event*, std::vector<Event*>, EventComparator> m_event_queue; // The event queue for Fortune's algorithm
    Arc* m_beach_line_root = nullptr; // The root of the beach line binary search tree
    Arc* m_last_arc = nullptr; // The last arc in the beach line
    Arc* m_first_arc = nullptr; // The first arc in the beach line
    std::uint64_t m_event_id_counter = 0; // A counter for generating unique event identifiers
};

VoronoiEngine::VoronoiDiagram VoronoiEngine::computeVoronoiDiagram() 
{
    // Implementation of Fortune's algorithm to compute the Voronoi diagram
    // This is a placeholder for the actual implementation
    VoronoiDiagram diagram;
    if (m_sites.empty()) {
        return diagram; // Return an empty diagram if there are no sites
    }
    if (allSitesCollinear()) {
        // Handle the case where all sites are collinear
        computeCollinearEdges();
    }else if (Point center; allSitesCoCircular(center)) {
        // Handle the case where all sites are co-circular
        computeCoCircularEdges(center);
    } else {
        // Handle the general case using Fortune's algorithm
        initializeEventQueue();
        processEventQueue();
    }

    std::vector<ClippedEdge> clipped_edges;
    clipped_edges.reserve(m_edge_records.size());
    for (const EdgeRecord& edge : m_edge_records) {
        const std::optional<ClippedEdge> clipped_edge = clipEdge(edge);
        if (!clipped_edge.has_value()) {
            continue;
        }
        clipped_edges.push_back(*clipped_edge);
        diagram.edges.push_back({
            m_sites[clipped_edge->first_site_index],
            m_sites[clipped_edge->second_site_index],
            clipped_edge->start,
            clipped_edge->end,
        });
    }
    diagram.cells = buildCells(clipped_edges);

    return diagram;
}

bool VoronoiEngine::allSitesCollinear() const 
{
    if (m_sites.size() < 3) {
        return true; // Less than 3 points are always collinear
    }
    const auto& p0 = m_sites[0];
    std::size_t second_index = 1;
    // Find the second point that is not the same as the first point
    while (second_index < m_sites.size() && std::abs(m_sites[second_index].x - p0.x) < EPSILON && std::abs(m_sites[second_index].y - p0.y) < EPSILON) {
        ++second_index;
    }
    if (second_index == m_sites.size()) {
        return true; // All points are the same
    }
    const auto& p1 = m_sites[second_index];
    for (std::size_t i = second_index + 1; i < m_sites.size(); ++i) {
        const auto& p2 = m_sites[i];
        // Check if the area of the triangle formed by p0, p1, and p2 is zero
        if (std::abs((p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y)) > EPSILON) {
            return false; // Found a point that is not collinear
        }
    }
    return true; // All points are collinear
}

void VoronoiEngine::computeCollinearEdges()
{
    if (m_sites.size() < 2) {
        return; // Not enough sites to form a Voronoi diagram
    }
    const auto& p0 = m_sites[0];
    std::size_t second_index = 1;
    // Find the second point that is not the same as the first point
    while (second_index < m_sites.size() && std::abs(m_sites[second_index].x - p0.x) < EPSILON && std::abs(m_sites[second_index].y - p0.y) < EPSILON) {
        ++second_index;
    }
    const auto& direction = m_sites[second_index] - p0;
    std::vector<std::size_t> sorted_indices(m_sites.size());
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    std::sort(sorted_indices.begin(), sorted_indices.end(), [&](std::size_t i, std::size_t j) {
        return dot(m_sites[i] - p0, direction) < dot(m_sites[j] - p0, direction);
    });
    for (std::size_t i = 0; i + 1 < sorted_indices.size(); ++i) {
        findEdgeRecord(sorted_indices[i], sorted_indices[i + 1]);
    }
}

bool VoronoiEngine::allSitesCoCircular(Point& center) const 
{
    if (m_sites.size() < 3) {
        return false;
    }

    std::size_t second_index = 1;
    // Find the second and third points that are not the same as the first point
    while (second_index < m_sites.size() && std::abs(m_sites[second_index].x - m_sites[0].x) < EPSILON && std::abs(m_sites[second_index].y - m_sites[0].y) < EPSILON) {
        ++second_index;
    }
    if (second_index == m_sites.size()) {
        return false;
    }

    std::size_t third_index = second_index + 1;
    while (third_index < m_sites.size() && std::abs(cross(m_sites[second_index] - m_sites[0], m_sites[third_index] - m_sites[0])) < EPSILON) {
        ++third_index;
    }
    if (third_index == m_sites.size()) {
        return false; // All points are collinear, cannot form a circumcircle
    }

    if (!circumCircle(m_sites[0], m_sites[second_index], m_sites[third_index], center)) {
        return false;
    }

    const double radius_sq = dot(center - m_sites[0], center - m_sites[0]);
    for (const Point& site : m_sites) {
        const Point offset = site - center;
        const double distance_sq = dot(offset, offset);
        const double tolerance = EPSILON * std::max({1.0, distance_sq});
        if (std::abs(distance_sq - radius_sq) > tolerance) {
            return false;
        }
    }
    return true;
}

void VoronoiEngine::computeCoCircularEdges(const Point& center) 
{
    std::vector<std::pair<double, std::size_t>> angles;
    for (std::size_t i = 0; i < m_sites.size(); ++i) {
        const Point offset = m_sites[i] - center;
        const double angle = std::atan2(offset.y, offset.x);
        angles.emplace_back(angle, i);
    }
    std::sort(angles.begin(), angles.end());
    for (std::size_t i = 0; i < angles.size(); ++i) {
        const std::size_t first_index = angles[i].second;
        const std::size_t second_index = angles[(i + 1) % angles.size()].second;
        EdgeRecord* edge = findEdgeRecord(first_index, second_index);
        addEdgeVertex(edge, center, angles[(i + 2) % angles.size()].second);
    }
}

bool VoronoiEngine::circumCircle(const Point& a, const Point& b, const Point& c, Point& center) const 
{
    // For center = a + u, equal distances to a, b, and c give:
    // 2 * dot(u, b - a) = |b - a|^2 and
    // 2 * dot(u, c - a) = |c - a|^2.
    // Solving this system yields the offsets below. Coincident or nearly
    // collinear points do not define a stable circumcircle.
    const Point ab_dir = b - a;
    const Point ac_dir = c - a;
    const double ab_len_sq = dot(ab_dir, ab_dir);
    const double ac_len_sq = dot(ac_dir, ac_dir);
    const double length_product = std::sqrt(ab_len_sq * ac_len_sq);
    const double cross_product = cross(ab_dir, ac_dir);
    if (length_product == 0.0 ||
        std::abs(cross_product) <= EPSILON * length_product) {
        return false;
    }

    const double d = 2.0 * cross_product;
    center.x = a.x + (ac_dir.y * ab_len_sq - ab_dir.y * ac_len_sq) / d;
    center.y = a.y + (ab_dir.x * ac_len_sq - ac_dir.x * ab_len_sq) / d;
    return std::isfinite(center.x) && std::isfinite(center.y);
}

void VoronoiEngine::initializeSweepSites()
{
    m_sweep_sites = m_sites;
    std::vector<std::size_t> ordered(m_sites.size());
    std::iota(ordered.begin(), ordered.end(), 0);
    std::sort(ordered.begin(), ordered.end(),
              [&](std::size_t first, std::size_t second) {
                  if (m_sites[first].y != m_sites[second].y) {
                      return m_sites[first].y > m_sites[second].y;
                  }
                  return m_sites[first].x < m_sites[second].x;
              });

    for (std::size_t begin = 0; begin < ordered.size();) {
        std::size_t end = begin + 1;
        while (end < ordered.size() &&
               m_sites[ordered[begin]].y == m_sites[ordered[end]].y) {
            ++end;
        }
        if (end - begin > 1) {
            double nearest_level_distance =
                std::numeric_limits<double>::infinity();
            if (begin > 0) {
                nearest_level_distance =
                    m_sites[ordered[begin - 1]].y -
                    m_sites[ordered[begin]].y;
            }
            if (end < ordered.size()) {
                nearest_level_distance = std::min(
                    nearest_level_distance,
                    m_sites[ordered[begin]].y - m_sites[ordered[end]].y);
            }
            const double minimum_nudge =
                32.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, std::abs(m_sites[ordered[begin]].y));
            const double nudge =
                std::max(nearest_level_distance * 1e-8, minimum_nudge);
            if (std::isfinite(nearest_level_distance) &&
                nudge * static_cast<double>(end - begin) <
                    nearest_level_distance / 4.0) {
                for (std::size_t i = begin; i < end; ++i) {
                    m_sweep_sites[ordered[i]].y +=
                        nudge * static_cast<double>(i - begin);
                }
            }
        }
        begin = end;
    }
}

void VoronoiEngine::initializeEventQueue() 
{
    const double highest_y = std::max_element(
        m_sites.begin(), m_sites.end(),
        [](const Point& first, const Point& second) {
            return first.y < second.y;
        })->y;
    std::vector<std::size_t> highest_sites;
    for (std::size_t i = 0; i < m_sites.size(); ++i) {
        if (m_sites[i].y == highest_y) {
            highest_sites.push_back(i);
        } else {
            m_event_queue.push(
                makeEvent(EventType::SITE_EVENT, m_sweep_sites[i], i));
        }
    }
    std::sort(highest_sites.begin(), highest_sites.end(),
              [&](std::size_t first, std::size_t second) {
                  return m_sites[first].x < m_sites[second].x;
              });

    Arc* previous = nullptr;
    for (std::size_t site_index : highest_sites) {
        Arc* arc = makeArc(site_index);
        arc->prev = previous;
        if (previous != nullptr) {
            previous->next = arc;
            previous->edge_to_next =
                findEdgeRecord(previous->site_index, site_index);
        } else {
            m_first_arc = arc;
        }
        insertBetween(previous, nullptr, arc);
        previous = arc;
    }
    m_last_arc = previous;
}

void VoronoiEngine::processEventQueue() 
{
    while (!m_event_queue.empty()) {
        Event* event = m_event_queue.top();
        m_event_queue.pop();
        if (!event->valid) {
            continue; // Skip invalidated events
        }
        // Process the event here
        if (event->type == EventType::SITE_EVENT) {
            processSiteEvent(event);
        } else if (event->type == EventType::CIRCLE_EVENT) {
            processCircleEvent(event);
        }
    }
}

void VoronoiEngine::processSiteEvent(Event* event) 
{
    if (m_beach_line_root == nullptr) {
        m_beach_line_root = makeArc(event->site_index);
        m_first_arc = m_beach_line_root;
        m_last_arc = m_beach_line_root;
        return;
    }
    Arc* spilt =
        findArcAbove(m_sweep_sites[event->site_index].x, event->y);
    invalidateCircleEvent(spilt);
    // Split the arc and create new arcs for the beach line
    Arc* left_arc = makeArc(spilt->site_index);
    Arc* center_arc = makeArc(event->site_index);
    Arc* right_arc = makeArc(spilt->site_index);

    Arc* prev_arc = spilt->prev;
    Arc* next_arc = spilt->next;
    EdgeRecord* edge_record = findEdgeRecord(spilt->site_index, event->site_index);

    left_arc->prev = prev_arc;
    if (prev_arc) {
        prev_arc->next = left_arc;
    } else {
        m_first_arc = left_arc;
    }

    left_arc->next = center_arc;
    center_arc->prev = left_arc;
    center_arc->next = right_arc;
    right_arc->prev = center_arc;

    right_arc->next = next_arc;
    if (next_arc) {
        next_arc->prev = right_arc;
    } else {
        m_last_arc = right_arc;
    }

    left_arc->edge_to_next = edge_record;
    center_arc->edge_to_next = edge_record;
    right_arc->edge_to_next = spilt->edge_to_next;

    removeArc(spilt);
    insertBetween(prev_arc, next_arc, left_arc);
    insertBetween(left_arc, next_arc, center_arc);
    insertBetween(center_arc, next_arc, right_arc);

    scheduleCircleEvent(left_arc, event->y);
    scheduleCircleEvent(right_arc, event->y);
}

void VoronoiEngine::processCircleEvent(Event* event) 
{
    Arc* middle = event->arc;
    if (middle == nullptr || middle->prev == nullptr || 
        middle->next == nullptr || middle->circle_event != event) {
        return;
    }

    Arc* first = middle->prev;
    Arc* third = middle->next;
    invalidateCircleEvent(first);
    invalidateCircleEvent(third);

    Point center;
    if (!circumCircle(m_sites[first->site_index], m_sites[middle->site_index],
                      m_sites[third->site_index], center)) {
        return;
    }
    addEdgeVertex(first->edge_to_next, center, third->site_index);
    addEdgeVertex(middle->edge_to_next, center, first->site_index);
    EdgeRecord* new_edge = findEdgeRecord(first->site_index, third->site_index);
    addEdgeVertex(new_edge, center, middle->site_index);

    first->next = third;
    first->edge_to_next = new_edge;
    third->prev = first;
    removeArc(middle);

    scheduleCircleEvent(first, event->y);
    scheduleCircleEvent(third, event->y);
}

void VoronoiEngine::scheduleCircleEvent(Arc* arc, double sweep_line_y) 
{
    if (arc == nullptr || arc->prev == nullptr || arc->next == nullptr) {
        return;
    }
    // Implementation for scheduling a circle event for the given arc
    // This typically involves computing the circumcircle of the arc and its neighbors
    // and adding a new event to the event queue if a valid circle is found.
    const Point& a = m_sweep_sites[arc->prev->site_index];
    const Point& b = m_sweep_sites[arc->site_index];
    const Point& c = m_sweep_sites[arc->next->site_index];

    if (cross(b - a, c - a) >= -EPSILON) {
        return; // Points are collinear or clockwise, no valid circle event
    }
    if (std::abs(cross(
            m_sites[arc->site_index] - m_sites[arc->prev->site_index],
            m_sites[arc->next->site_index] -
                m_sites[arc->prev->site_index])) <= EPSILON) {
        return;
    }

    // Compute the circumcircle of points a, b, c
    Point circumcenter;
    if (!circumCircle(a, b, c, circumcenter)) {
        return; // Failed to compute circumcenter
    }
    const double event_y = circumcenter.y - std::sqrt((b.x - circumcenter.x) * (b.x - circumcenter.x) + 
                                                      (b.y - circumcenter.y) * (b.y - circumcenter.y));
    if (!std::isfinite(event_y) || event_y >= sweep_line_y - EPSILON) {
        return; // Invalid event_y, do not schedule the circle event
    }
    Event* event = makeEvent(
        EventType::CIRCLE_EVENT,
        Point{circumcenter.x, event_y},
        arc->site_index);
    event->arc = arc;
    event->circle_center = circumcenter;
    arc->circle_event = event;
    m_event_queue.push(event);
}

void VoronoiEngine::invalidateCircleEvent(Arc* arc) 
{
    if (arc && arc->circle_event) {
        arc->circle_event->valid = false;
        arc->circle_event = nullptr;
    }
}

std::optional<VoronoiEngine::ClippedEdge> VoronoiEngine::clipEdge(const VoronoiEngine::EdgeRecord& edge) const
{
    const auto& first_site = m_sites[edge.first_site_index];
    const auto& second_site = m_sites[edge.second_site_index];
    const Point direction{first_site.y - second_site.y, second_site.x - first_site.x};
    
    Point start, end;
    if (edge.completed_break_points >= 2) {
        if (edge.vertices.size() < 2) {
            return std::nullopt; // Not enough vertices to form a clipped edge
        }
        const Point segment_direction = edge.vertices[1].point - edge.vertices[0].point;
        if (!clipParameterRange(edge.vertices[0].point, segment_direction, 0.0, 1.0, start, end)) {
            return std::nullopt; // Segment does not intersect the bounding box
        }
    }else if (edge.completed_break_points == 1) {
        if (edge.vertices.empty()) {
            return std::nullopt; // Not enough vertices to form a clipped edge
        }
        const auto& vertex = edge.vertices[0];
        Point ray_direction = direction;
        const Point to_opposite = m_sites[vertex.opposite] - first_site;
        if (dot(to_opposite, ray_direction) > EPSILON) {
            ray_direction = ray_direction * (-1.0); // Reverse the ray direction to point towards the opposite site
        }
        if (!clipParameterRange(vertex.point, ray_direction, 0.0, std::numeric_limits<double>::infinity(), start, end)) {
            return std::nullopt; // Ray does not intersect the bounding box
        }
    }else if (!clipParameterRange((first_site + second_site) * 0.5, direction, 
              -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), start, end)) {
        return std::nullopt; // Segment does not intersect the bounding box
    }

    const Point clipped_direction = end - start;
    double minimum_parameter = 0.0;
    double maximum_parameter = 1.0;
    for (const Point& site : m_sites) {
        const Point normal = site - first_site;
        const double coefficient =
            2.0 * dot(clipped_direction, normal);
        const double limit =
            dot(site, site) - dot(first_site, first_site) -
            2.0 * dot(start, normal);
        if (std::abs(coefficient) <= EPSILON) {
            if (limit < -EPSILON) {
                return std::nullopt;
            }
            continue;
        }

        const double parameter = limit / coefficient;
        if (coefficient > 0.0) {
            maximum_parameter = std::min(maximum_parameter, parameter);
        } else {
            minimum_parameter = std::max(minimum_parameter, parameter);
        }
        if (minimum_parameter > maximum_parameter + EPSILON) {
            return std::nullopt;
        }
    }
    start = clampToBounds(
        start + clipped_direction * std::max(0.0, minimum_parameter));
    end = clampToBounds(
        start + clipped_direction *
                    (std::min(1.0, maximum_parameter) -
                     std::max(0.0, minimum_parameter)));
    if (dot(end - start, end - start) <= EPSILON) {
        return std::nullopt;
    }
    return ClippedEdge{edge.first_site_index, edge.second_site_index, start, end};
}

bool VoronoiEngine::clipParameterRange(const Point& origin, const Point& direction, 
                                       double min_param, double max_param, Point& start, Point& end) const 
{
    const auto update_range = [&](double origin_coordinate, double direction_coordinate, 
                                  double lower, double upper, double& min, double& max) {
        if (std::abs(direction_coordinate) < EPSILON) {
            return origin_coordinate >= lower - EPSILON && origin_coordinate <= upper + EPSILON;
        }
        auto first = (lower - origin_coordinate) / direction_coordinate;
        auto second = (upper - origin_coordinate) / direction_coordinate;
        if (first > second) std::swap(first, second);
        min = std::max(min, first);
        max = std::min(max, second);
        return min <= max + EPSILON;
    };

    if (!update_range(origin.x, direction.x, m_bounding_box.min.x, m_bounding_box.max.x, min_param, max_param)) {
        return false;
    }
    if (!update_range(origin.y, direction.y, m_bounding_box.min.y, m_bounding_box.max.y, min_param, max_param)) {
        return false;
    }
    start = clampToBounds(origin + direction * min_param);
    end = clampToBounds(origin + direction * max_param);
    return dot(end - start, end - start) > EPSILON;
}

VoronoiEngine::Point VoronoiEngine::clampToBounds(const Point& point) const 
{
    Point clamped_point = point;
    clamped_point.x = std::clamp(point.x, m_bounding_box.min.x, m_bounding_box.max.x);
    clamped_point.y = std::clamp(point.y, m_bounding_box.min.y, m_bounding_box.max.y);
    return clamped_point;
}

std::vector<VoronoiEngine::VoronoiCell> VoronoiEngine::buildCells(const std::vector<ClippedEdge>& clipped_edges)
{
    std::vector<std::vector<Point>> vertices(m_sites.size());
    for (const auto& clipped_edge : clipped_edges) {
        vertices[clipped_edge.first_site_index].push_back(clipped_edge.start);
        vertices[clipped_edge.first_site_index].push_back(clipped_edge.end);
        vertices[clipped_edge.second_site_index].push_back(clipped_edge.start);
        vertices[clipped_edge.second_site_index].push_back(clipped_edge.end);
    }
    const std::vector<Point> corner_points = {
        m_bounding_box.min,
        {m_bounding_box.min.x, m_bounding_box.max.y},
        m_bounding_box.max,
        {m_bounding_box.max.x, m_bounding_box.min.y}
    };
    for (const auto& corner_point : corner_points) {
        double closest_distance = std::numeric_limits<double>::infinity();
        for (const auto& site : m_sites) {
            closest_distance = std::min(closest_distance, dot(site - corner_point, site - corner_point));
        }
        for (std::size_t site_index = 0; site_index < m_sites.size(); ++site_index) {
            if (std::abs(dot(m_sites[site_index] - corner_point, m_sites[site_index] - corner_point) - closest_distance) <= EPSILON) {
                vertices[site_index].push_back(corner_point);
            }
        }
    }

    std::vector<VoronoiCell> cells(m_sites.size());
    for (std::size_t site_index = 0; site_index < m_sites.size(); ++site_index) {
        std::vector<Point> unique_vertices = vertices[site_index];
        std::sort(unique_vertices.begin(), unique_vertices.end(), [&](const Point& a, const Point& b) {
            const auto angle_a = atan2(a.y - m_sites[site_index].y, a.x - m_sites[site_index].x);
            const auto angle_b = atan2(b.y - m_sites[site_index].y, b.x - m_sites[site_index].x);
            return angle_a != angle_b ? angle_a < angle_b : 
                        dot(a - m_sites[site_index], a - m_sites[site_index]) < 
                        dot(b - m_sites[site_index], b - m_sites[site_index]);
        });
        unique_vertices.erase(std::unique(unique_vertices.begin(), unique_vertices.end(), [](const Point& a, const Point& b) {
            return std::abs(a.x - b.x) < EPSILON && std::abs(a.y - b.y) < EPSILON;
        }), unique_vertices.end());
        cells[site_index].site = m_sites[site_index];
        cells[site_index].vertices = std::move(unique_vertices);
    }

    return cells;
}

std::pair<std::size_t, std::size_t> VoronoiEngine::makeOrderedPair(std::size_t first, std::size_t second) const 
{
    return (first < second) ? std::make_pair(first, second) : std::make_pair(second, first);
}

VoronoiEngine::Event* VoronoiEngine::makeEvent(EventType type, const Point& point, std::size_t site_index) {
    m_events.emplace_back(std::make_unique<Event>(Event{type, point.y, point.x, site_index, nullptr, {}, true, m_event_id_counter++}));
    return m_events.back().get();
}

VoronoiEngine::EdgeRecord* VoronoiEngine::findEdgeRecord(std::size_t first_site_index, std::size_t second_site_index) 
{
    auto ordered_pair = makeOrderedPair(first_site_index, second_site_index);
    auto it = m_edge_map.find(ordered_pair);
    if (it != m_edge_map.end()) {
        return it->second;
    }
    m_edge_records.emplace_back(
        EdgeRecord{first_site_index, second_site_index, {}, 0});
    EdgeRecord* new_edge_record = &m_edge_records.back();
    m_edge_map[ordered_pair] = new_edge_record;
    return new_edge_record;
}

void VoronoiEngine::addEdgeVertex(EdgeRecord* edge_record, const Point& vertex_point, std::size_t opposite_vertex_index) 
{
    assert (edge_record != nullptr);
    ++ edge_record->completed_break_points;
    for (const auto& vertex : edge_record->vertices) {
        if (std::abs(vertex.point.x - vertex_point.x) < EPSILON && std::abs(vertex.point.y - vertex_point.y) < EPSILON) {
            return; // Vertex already exists
        }
    }
    edge_record->vertices.emplace_back(EdgeVertex{vertex_point, opposite_vertex_index});
}

VoronoiEngine::Arc* VoronoiEngine::makeArc(std::size_t site_index) 
{
    m_arcs.emplace_back(std::make_unique<Arc>());
    m_arcs.back()->site_index = site_index;
    return m_arcs.back().get();
}

void VoronoiEngine::replaceArc(Arc* old_arc, Arc* new_arc) 
{
    // Replace old_arc with new_arc in the beach line binary search tree
    if (old_arc->parent) {
        if (old_arc->parent->left_child == old_arc) {
            old_arc->parent->left_child = new_arc;
        } else {
            old_arc->parent->right_child = new_arc;
        }
    } else {
        m_beach_line_root = new_arc;
    }
    if (new_arc) {
        new_arc->parent = old_arc->parent;
    }
}

void VoronoiEngine::rotateLeft(Arc* arc) 
{
    Arc* right_child = arc->right_child;
    replaceArc(arc, right_child);
    arc->right_child = right_child->left_child;
    if (right_child->left_child) {
        right_child->left_child->parent = arc;
    }
    right_child->left_child = arc;
    arc->parent = right_child;
    updateHeight(arc);
    updateHeight(right_child);
}

void VoronoiEngine::rotateRight(Arc* arc) 
{
    Arc* left_child = arc->left_child;
    replaceArc(arc, left_child);
    arc->left_child = left_child->right_child;
    if (left_child->right_child) {
        left_child->right_child->parent = arc;
    }
    left_child->right_child = arc;
    arc->parent = left_child;
    updateHeight(arc);
    updateHeight(left_child);
}

void VoronoiEngine::rebalanceFrom(Arc* node) 
{
    while (node) {
        updateHeight(node);
        int balance = getBalanceFactor(node);
        if (balance > 1) {
            if (getBalanceFactor(node->left_child) < 0) {
                rotateLeft(node->left_child);
            }
            rotateRight(node);
        } else if (balance < -1) {
            if (getBalanceFactor(node->right_child) > 0) {
                rotateRight(node->right_child);
            }
            rotateLeft(node);
        }
        node = node->parent;
    }
}

void VoronoiEngine::insertBetween(Arc* left_arc, Arc* right_arc, Arc* new_arc) 
{
    if (m_beach_line_root == nullptr) {
        m_beach_line_root = new_arc;
        new_arc->parent = nullptr;
        return;
    }
    assert(left_arc || right_arc); // Ensure at least one of the arcs is not null
    Arc* current = nullptr;
    if (left_arc){
        current = left_arc;
        if (current->right_child != nullptr) {
            current = current->right_child;
            while (current->left_child != nullptr) {
                current = current->left_child;
            }
            current->left_child = new_arc;
        }else {
            current->right_child = new_arc;
        }
    }else {
        current = right_arc;
        if (current->left_child != nullptr) {
            current = current->left_child;
            while (current->right_child != nullptr) {
                current = current->right_child;
            }
            current->right_child = new_arc;
        }else {
            current->left_child = new_arc;
        }
    }
    new_arc->parent = current;
    rebalanceFrom(current);
}

void VoronoiEngine::removeArc(Arc* arc) 
{
    Arc* current = nullptr;
    if (arc->left_child == nullptr) {
        current = arc->parent;
        replaceArc(arc, arc->right_child);
        if (current == nullptr) {
            current = arc->right_child;
        }
    }else if (arc->right_child == nullptr) {
        current = arc->parent;
        replaceArc(arc, arc->left_child);
        if (current == nullptr) {
            current = arc->left_child;
        }
    }else {
        Arc* successor = arc->right_child;
        while (successor->left_child != nullptr) {
            successor = successor->left_child;
        }
        if (successor->parent != arc) {
            current = successor->parent;
            replaceArc(successor, successor->right_child);
            successor->right_child = arc->right_child;
            successor->right_child->parent = successor;
        }else {
            current = successor;
        }
        replaceArc(arc, successor);
        successor->left_child = arc->left_child;
        successor->left_child->parent = successor;
        updateHeight(successor);
        rebalanceFrom(current);
    }

    arc->parent = nullptr;
    arc->left_child = nullptr;
    arc->right_child = nullptr;
    arc->height = 1;
    rebalanceFrom(current);
}

double VoronoiEngine::breakPointX(std::size_t left_site_index, std::size_t right_site_index, double sweep_line_y) const 
{
    const Point& left_site = m_sweep_sites[left_site_index];
    const Point& right_site = m_sweep_sites[right_site_index];

    if (std::abs(left_site.y - right_site.y) < EPSILON) {
        return (left_site.x + right_site.x) / 2.0;
    }
    if (std::abs(left_site.y - sweep_line_y) < EPSILON) {
        return left_site.x;
    }
    if (std::abs(right_site.y - sweep_line_y) < EPSILON) {
        return right_site.x;
    }
    const double d_left = 2.0 * (left_site.y - sweep_line_y);
    const double d_right = 2.0 * (right_site.y - sweep_line_y);
    const double a = 1.0 / d_left - 1.0 / d_right;
    const double b = -2.0 * (left_site.x / d_left - right_site.x / d_right);
    const double c = (left_site.x * left_site.x + left_site.y * left_site.y - sweep_line_y * sweep_line_y) / d_left
                   - (right_site.x * right_site.x + right_site.y * right_site.y - sweep_line_y * sweep_line_y) / d_right;
    if (a == 0.0) return -c / b;
    double discriminant = b * b - 4 * a * c;
    if (discriminant < -EPSILON) {
        return (left_site.x + right_site.x) / 2.0;
    }
    if (discriminant < EPSILON) {
        return (-b) / (2 * a);
    }
    const double sqrt_discriminant = std::sqrt(discriminant);
    const double stable_term = -0.5 * (b + std::copysign(sqrt_discriminant, b));
    const double root1 = stable_term / a;
    const double root2 = stable_term == 0.0 ? -b / (2 * a) : c / stable_term;
    return (left_site.y > right_site.y) ? std::min(root1, root2) : std::max(root1, root2);
}

VoronoiEngine::Arc* VoronoiEngine::findArcAbove(double x, double sweep_line_y) const 
{
    Arc* current = m_beach_line_root;
    Arc* candidate = nullptr;
    while (current) {
        if (current->next == nullptr || x < breakPointX(current->site_index, current->next->site_index, sweep_line_y)) {
            candidate = current;
            current = current->left_child;
        } else {
            current = current->right_child;
        }
    }
    return candidate == nullptr ? m_last_arc : candidate;
}

double squaredDistance(const VoronoiEngine::Point& first,
                       const VoronoiEngine::Point& second)
{
    const double dx = first.x - second.x;
    const double dy = first.y - second.y;
    return dx * dx + dy * dy;
}

double cross(const VoronoiEngine::Point& first,
             const VoronoiEngine::Point& second,
             const VoronoiEngine::Point& third)
{
    return (second.x - first.x) * (third.y - first.y) -
           (second.y - first.y) * (third.x - first.x);
}

bool pointOnSegment(const VoronoiEngine::Point& point,
                    const VoronoiEngine::Point& start,
                    const VoronoiEngine::Point& end)
{
    return std::abs(cross(start, end, point)) <= EPSILON &&
           point.x >= std::min(start.x, end.x) - EPSILON &&
           point.x <= std::max(start.x, end.x) + EPSILON &&
           point.y >= std::min(start.y, end.y) - EPSILON &&
           point.y <= std::max(start.y, end.y) + EPSILON;
}

bool pointInPolygon(const VoronoiEngine::Point& point,
                    const std::vector<VoronoiEngine::Point>& polygon)
{
    bool inside = false;
    for (std::size_t i = 0, previous = polygon.size() - 1;
         i < polygon.size(); previous = i++) {
        const auto& first = polygon[previous];
        const auto& second = polygon[i];
        if (pointOnSegment(point, first, second)) {
            return true;
        }
        if ((first.y > point.y) != (second.y > point.y)) {
            const double intersection_x =
                first.x + (second.x - first.x) *
                              (point.y - first.y) /
                              (second.y - first.y);
            if (point.x < intersection_x) {
                inside = !inside;
            }
        }
    }
    return inside;
}

std::vector<VoronoiEngine::Point> segmentIntersections(
    const VoronoiEngine::Point& first_start,
    const VoronoiEngine::Point& first_end,
    const VoronoiEngine::Point& second_start,
    const VoronoiEngine::Point& second_end)
{
    const auto first_direction = first_end - first_start;
    const auto second_direction = second_end - second_start;
    const double denominator =
        first_direction.x * second_direction.y -
        first_direction.y * second_direction.x;
    if (std::abs(denominator) <= EPSILON) {
        std::vector<VoronoiEngine::Point> intersections;
        if (pointOnSegment(first_start, second_start, second_end)) {
            intersections.push_back(first_start);
        }
        if (pointOnSegment(first_end, second_start, second_end)) {
            intersections.push_back(first_end);
        }
        if (pointOnSegment(second_start, first_start, first_end)) {
            intersections.push_back(second_start);
        }
        if (pointOnSegment(second_end, first_start, first_end)) {
            intersections.push_back(second_end);
        }
        return intersections;
    }

    const double offset_x = second_start.x - first_start.x;
    const double offset_y = second_start.y - first_start.y;
    const double first_parameter =
        (offset_x * second_direction.y -
         offset_y * second_direction.x) /
        denominator;
    const double second_parameter =
        (offset_x * first_direction.y -
         offset_y * first_direction.x) /
        denominator;
    if (first_parameter < -EPSILON || first_parameter > 1.0 + EPSILON ||
        second_parameter < -EPSILON || second_parameter > 1.0 + EPSILON) {
        return {};
    }
    return {{
        first_start.x + first_direction.x * first_parameter,
        first_start.y + first_direction.y * first_parameter,
    }};
}

double boundaryMaximumSquared(
    const std::vector<VoronoiEngine::Point>& polygon)
{
    struct Line {
        double slope;
        double intercept;
        double start;
    };

    double result = 0.0;
    for (std::size_t edge_index = 0;
         edge_index < polygon.size(); ++edge_index) {
        const auto& start = polygon[edge_index];
        const auto& end = polygon[(edge_index + 1) % polygon.size()];
        const auto direction = end - start;
        const double quadratic = squaredDistance(start, end);

        std::vector<Line> lines;
        lines.reserve(polygon.size());
        for (const auto& site : polygon) {
            const auto offset = start - site;
            lines.push_back({
                2.0 * (direction.x * offset.x +
                       direction.y * offset.y),
                offset.x * offset.x + offset.y * offset.y,
                0.0,
            });
        }
        std::sort(lines.begin(), lines.end(),
                  [](const Line& first, const Line& second) {
                      if (first.slope != second.slope) {
                          return first.slope > second.slope;
                      }
                      return first.intercept < second.intercept;
                  });

        std::vector<Line> envelope;
        envelope.reserve(lines.size());
        for (Line line : lines) {
            if (!envelope.empty() &&
                line.slope == envelope.back().slope) {
                continue;
            }
            line.start = -std::numeric_limits<double>::infinity();
            while (!envelope.empty()) {
                const Line& previous = envelope.back();
                const double intersection =
                    (line.intercept - previous.intercept) /
                    (previous.slope - line.slope);
                if (intersection > previous.start) {
                    line.start = intersection;
                    break;
                }
                envelope.pop_back();
            }
            envelope.push_back(line);
        }

        for (const Line& line : envelope) {
            if (line.start <= 0.0 || line.start >= 1.0) {
                continue;
            }
            const double parameter = line.start;
            result = std::max(
                result,
                quadratic * parameter * parameter +
                    line.slope * parameter + line.intercept);
        }
    }
    return result;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);
    int n; std::cin >> n;
    std::vector<VoronoiEngine::Point> sites(n);
    VoronoiEngine::Aabb bounding_box{
        VoronoiEngine::Point{
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()},
        VoronoiEngine::Point{
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()}};
    for (int i = 0; i < n; ++i) {
        std::cin >> sites[i].x >> sites[i].y;
        bounding_box.min.x = std::min(bounding_box.min.x, sites[i].x);
        bounding_box.min.y = std::min(bounding_box.min.y, sites[i].y);
        bounding_box.max.x = std::max(bounding_box.max.x, sites[i].x);
        bounding_box.max.y = std::max(bounding_box.max.y, sites[i].y);
    }
    VoronoiEngine engine(sites, bounding_box);
    const auto result = engine.computeVoronoiDiagram();
    double max_squared_radius = boundaryMaximumSquared(sites);
    for (const auto& edge : result.edges) {
        if (pointInPolygon(edge.start, sites)) {
            max_squared_radius = std::max(
                max_squared_radius,
                squaredDistance(edge.start, edge.first_site));
        }
        if (pointInPolygon(edge.end, sites)) {
            max_squared_radius = std::max(
                max_squared_radius,
                squaredDistance(edge.end, edge.first_site));
        }
        for (std::size_t i = 0; i < sites.size(); ++i) {
            const auto& fence_start = sites[i];
            const auto& fence_end = sites[(i + 1) % sites.size()];
            for (const auto& intersection :
                 segmentIntersections(edge.start, edge.end,
                                      fence_start, fence_end)) {
                max_squared_radius = std::max(
                    max_squared_radius,
                    squaredDistance(intersection, edge.first_site));
            }
        }
    }
    std::cout << std::setprecision(12)
              << std::sqrt(max_squared_radius) << "\n";
    return 0;
}