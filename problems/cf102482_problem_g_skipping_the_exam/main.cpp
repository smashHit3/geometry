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
        m_sites(std::move(sites)), m_bounding_box(std::move(bounding_box)) {}

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
            if (std::abs(lhs->y - rhs->y) > EPSILON) {
                return lhs->y < rhs->y; // Sort by y-coordinate
            }
            if (lhs->type != rhs->type) {
                return lhs->type == EventType::CIRCLE_EVENT; // Circle events have higher priority than site events
            }
            if (std::abs(lhs->x - rhs->x) > EPSILON) {
                return lhs->x < rhs->x; // Sort by x-coordinate
            }
            return lhs->id < rhs->id; // Sort by unique identifier
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
    void initializeEventQueue();
    void processEventQueue();
    void processSiteEvent(Event* event);
    void processCircleEvent(Event* event);

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
    for (std::size_t i = 1; i < sorted_indices.size() - 1; ++i) {
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
        findEdgeRecord(first_index, second_index);
        addEdgeVertex(findEdgeRecord(first_index, second_index), center, (i + 2) % angles.size());
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

void VoronoiEngine::initializeEventQueue() {
    for (std::size_t i = 0; i < m_sites.size(); ++i) {
        m_event_queue.push(makeEvent(EventType::SITE_EVENT, m_sites[i], i));
    }
}

void VoronoiEngine::processEventQueue() {
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

void VoronoiEngine::processSiteEvent(Event* event) {
    // Implementation for processing a site event
}

void VoronoiEngine::processCircleEvent(Event* event) {
    // Implementation for processing a circle event
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
    m_edge_records.emplace_back(EdgeRecord{first_site_index, second_site_index});
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
        if (current) {
            current = arc->right_child;
        }
    }else if (arc->right_child == nullptr) {
        current = arc->parent;
        replaceArc(arc, arc->left_child);
        if (current) {
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
        rebalanceFrom(current);
    }

    arc->parent = nullptr;
    arc->left_child = nullptr;
    arc->right_child = nullptr;
    arc->height = 0;
    rebalanceFrom(current);
}

double VoronoiEngine::breakPointX(std::size_t left_site_index, std::size_t right_site_index, double sweep_line_y) const 
{
    const Point& left_site = m_sites[left_site_index];
    const Point& right_site = m_sites[right_site_index];

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
    return (left_site.y < right_site.y) ? std::min(root1, root2) : std::max(root1, root2);
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

int main() {

    return 0;
}