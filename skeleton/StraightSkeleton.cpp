#include "skeleton/StraightSkeleton.h"

#include "common/BasePointUtil.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <queue>
#include <stdexcept>
#include <vector>

namespace skeleton::detail {

namespace {

inline Point unit(const Point& p) {
    const double len = std::hypot(static_cast<double>(p.x()), static_cast<double>(p.y()));
    if (len == 0.0) return {0.0, 0.0};
    return {p.x() / len, p.y() / len};
}

} // namespace

struct StraightSkeletonBuilder::Node {
    Point pos;
    double time = 0.0;
    Point bisector{0.0, 0.0};
    bool reflex = false;
    int lav = 0;
    std::size_t origin = 0; // index into StraightSkeleton::nodes
    bool alive = true;
    Node* prev = nullptr;
    Node* next = nullptr;
    Event* edge_event = nullptr; // forward edge (this, next)
    Event* split_event = nullptr;
};

struct StraightSkeletonBuilder::Event {
    enum class Kind { Edge, Split } kind;
    double time = 0.0;
    bool valid = true;
    Point point;
    Node* v = nullptr;
    Node* w = nullptr; // edge: v.next; split: unused
    Node* a = nullptr; // split: edge start
    Node* b = nullptr; // split: edge end
};

StraightSkeletonBuilder::StraightSkeletonBuilder(std::vector<Point> polygon)
    : m_polygon(std::move(polygon)) {
    double scale = 1.0;
    for (const Point& v : m_polygon) {
        scale = std::max(scale, std::abs(static_cast<double>(v.x())));
        scale = std::max(scale, std::abs(static_cast<double>(v.y())));
    }
    m_scale = std::max(scale, 1.0);
    m_lengthEpsilon = m_scale * 1e-9;
    m_areaEpsilon = m_scale * m_scale * 1e-18;
}

StraightSkeleton StraightSkeletonBuilder::build() {
    StraightSkeleton result;
    if (m_polygon.size() < 3U) return result;

    // Orient CCW.
    double signed_area = 0.0;
    for (std::size_t i = 0; i < m_polygon.size(); ++i) {
        const std::size_t j = (i + 1U) % m_polygon.size();
        signed_area += static_cast<double>(m_polygon[i].x()) *
                           static_cast<double>(m_polygon[j].y()) -
                       static_cast<double>(m_polygon[i].y()) *
                           static_cast<double>(m_polygon[j].x());
    }
    if (signed_area < 0.0) std::reverse(m_polygon.begin(), m_polygon.end());

    result.nodes = m_polygon;

    std::vector<std::unique_ptr<Node>> nodes;
    std::vector<std::unique_ptr<Event>> events;
    struct EventLater {
        bool operator()(const Event* left, const Event* right) const {
            return left->time > right->time;
        }
    };
    std::priority_queue<Event*, std::vector<Event*>, EventLater> queue;

    const auto make_node = [&](std::size_t origin) {
        nodes.push_back(std::make_unique<Node>());
        Node* n = nodes.back().get();
        n->origin = origin;
        return n;
    };
    const auto make_event = [&] {
        events.push_back(std::make_unique<Event>());
        return events.back().get();
    };

    const auto position_at = [](const Node* n, double t) {
        return n->pos + (t - n->time) * n->bisector;
    };

    // Reuse an existing skeleton node when an event lands on an already-known
    // point (common for symmetric inputs where several events collapse to the
    // same location). Keeps the skeleton a clean planar graph.
    const auto add_node = [&](const Point& point) -> std::size_t {
        for (std::size_t i = 0; i < result.nodes.size(); ++i) {
            if (common::squaredDistance(result.nodes[i], point) <=
                m_lengthEpsilon * m_lengthEpsilon) {
                return i;
            }
        }
        result.nodes.push_back(point);
        return result.nodes.size() - 1U;
    };
    const auto push_arc = [&](std::size_t from, std::size_t to) {
        if (from == to) return;
        if (common::squaredDistance(result.nodes[from], result.nodes[to]) <=
            m_lengthEpsilon * m_lengthEpsilon) {
            return;
        }
        result.arcs.push_back({from, to});
    };

    const auto compute_bisector = [&](Node* n) {
        // Edge directions are constant, so evaluate every incident vertex at
        // this node's creation time to keep the directions consistent.
        const Point prev_at = position_at(n->prev, n->time);
        const Point next_at = position_at(n->next, n->time);
        const Point prev_dir = unit(n->pos - prev_at);
        const Point next_dir = unit(next_at - n->pos);
        // Inward (left) normals of the two incident wavefront edges.
        const Point n_prev{-prev_dir.y(), prev_dir.x()};
        const Point n_next{-next_dir.y(), next_dir.x()};
        n->reflex = common::cross(prev_dir, next_dir) < -m_areaEpsilon;
        // Vertex velocity v with v.n_prev = v.n_next = 1: each wavefront edge
        // advances at unit speed along its inward normal, so the vertex (the
        // intersection of the two moving edges) satisfies both. The unique
        // solution is v = (n_prev + n_next) / (1 + n_prev.n_next). Using the
        // true velocity (not the unit direction) is what makes adjacent edges
        // collapse consistently: v - w is along the shared edge, so the
        // perpendicular residual in meeting_time vanishes to floating-point
        // noise. The direction (n_prev + n_next) has non-negative dot with both
        // inward normals, so it points inward for convex and reflex vertices
        // alike; no sign flip is needed for reflex vertices.
        const double denom = 1.0 + common::dot(n_prev, n_next);
        if (denom < 1e-9) {
            // Degenerate (nearly straight) vertex: translate with the incident
            // edge at unit speed along its inward normal.
            n->bisector = n_prev;
        } else {
            const Point sum = n_prev + n_next;
            n->bisector = Point{sum.x() / denom, sum.y() / denom};
        }
    };

    const auto invalidate_node = [&](Node* n) {
        if (n->edge_event != nullptr) {
            n->edge_event->valid = false;
            n->edge_event = nullptr;
        }
        if (n->split_event != nullptr) {
            n->split_event->valid = false;
            n->split_event = nullptr;
        }
    };

    // Time at which two adjacent moving vertices meet (the forward edge
    // collapses). Returns -1 if they do not converge in the future.
    const auto meeting_time = [&](const Node* v, const Node* w, Point& meet) -> double {
        const Point base_v = v->pos - v->time * v->bisector;
        const Point base_w = w->pos - w->time * w->bisector;
        const Point c = base_v - base_w;
        const Point d = w->bisector - v->bisector;
        const double dx = static_cast<double>(d.x());
        const double dy = static_cast<double>(d.y());
        const double cx = static_cast<double>(c.x());
        const double cy = static_cast<double>(c.y());
        const bool use_x = std::fabs(dx) >= std::fabs(dy);
        const double denom = use_x ? dx : dy;
        if (std::fabs(denom) < m_areaEpsilon) return -1.0; // parallel bisectors
        const double comp = use_x ? cx : cy;
        const double t = comp / denom;
        if (!std::isfinite(t)) return -1.0;
        // Consistency in the other axis: for two adjacent wavefront vertices
        // the shared edge translates rigidly, so C and D are both along the
        // edge and the perpendicular residual is just floating-point noise.
        const double consistency = m_scale * 1e-6;
        const double other = use_x ? (cy - t * dy) : (cx - t * dx);
        if (std::fabs(other) > consistency) return -1.0;
        if (t < std::max(v->time, w->time) - m_lengthEpsilon) return -1.0;
        meet = v->pos + (t - v->time) * v->bisector;
        return t;
    };

    const auto schedule_edge_event = [&](Node* n) {
        if (!n->alive || n->next == nullptr || !n->next->alive) return;
        if (n->next == n) return; // single-node LAV
        Point meet;
        const double t = meeting_time(n, n->next, meet);
        if (t < 0.0) return;
        Event* e = make_event();
        e->kind = Event::Kind::Edge;
        e->time = t;
        e->point = meet;
        e->v = n;
        e->w = n->next;
        n->edge_event = e;
        queue.push(e);
    };

    // Earliest time a reflex vertex pierces an opposite edge of the same LAV.
    const auto schedule_split_event = [&](Node* r) {
        if (!r->alive || !r->reflex) return;
        const double t_birth = std::max({r->time, r->prev->time, r->next->time});
        double best_t = std::numeric_limits<double>::infinity();
        Node* best_a = nullptr;
        Node* best_b = nullptr;
        Point best_point{0.0, 0.0};
        for (const auto& candidate : nodes) {
            Node* a = candidate.get();
            if (!a->alive) continue;
            Node* b = a->next;
            if (b == nullptr || !b->alive) continue;
            if (a == r || b == r) continue;            // reflex's own edges
            if (a->lav != r->lav) continue;            // cross-loop => merge event
            const double t_dir = std::max({t_birth, a->time, b->time});
            const Point a_at = a->pos + (t_dir - a->time) * a->bisector;
            const Point b_at = b->pos + (t_dir - b->time) * b->bisector;
            const Point dir = unit(b_at - a_at);
            const Point r0 = r->pos - r->time * r->bisector;
            const Point a0 = a->pos - a->time * a->bisector;
            const Point c = r0 - a0;
            const Point v = r->bisector - a->bisector;
            const double cross_dv = common::cross(dir, v);
            if (std::fabs(cross_dv) < m_areaEpsilon) continue;
            const double t = -common::cross(dir, c) / cross_dv;
            if (!std::isfinite(t) || t <= t_birth + m_lengthEpsilon) continue;
            const Point p = r->pos + (t - r->time) * r->bisector;
            const Point a_t = a->pos + (t - a->time) * a->bisector;
            const Point b_t = b->pos + (t - b->time) * b->bisector;
            const double span = common::dot(b_t - a_t, dir);
            if (span <= m_lengthEpsilon) continue; // edge already collapsed
            const double s = common::dot(p - a_t, dir) / span;
            if (s < -m_lengthEpsilon || s > 1.0 + m_lengthEpsilon) continue;
            if (t < best_t) {
                best_t = t;
                best_a = a;
                best_b = b;
                best_point = p;
            }
        }
        if (best_a == nullptr) return;
        Event* e = make_event();
        e->kind = Event::Kind::Split;
        e->time = best_t;
        e->point = best_point;
        e->v = r;
        e->a = best_a;
        e->b = best_b;
        r->split_event = e;
        queue.push(e);
    };

    const auto recompute_all_splits = [&] {
        for (const auto& n : nodes) {
            if (!n->alive || !n->reflex) continue;
            if (n->split_event != nullptr) {
                n->split_event->valid = false;
                n->split_event = nullptr;
            }
            schedule_split_event(n.get());
        }
    };

    int next_lav = 0;
    std::size_t total_events = 0;
    const std::size_t event_cap = 32U * m_polygon.size() + 64U;

    const auto process_edge_event = [&](Event* e) {
        Node* v = e->v;
        Node* w = e->w;
        if (!v->alive || !w->alive || v->next != w) return;
        Node* p = v->prev;
        Node* q = w->next;

        const std::size_t s = add_node(e->point);
        push_arc(v->origin, s);
        push_arc(w->origin, s);

        Node* m = make_node(s);
        m->pos = e->point;
        m->time = e->time;
        m->lav = v->lav;
        m->prev = p;
        m->next = q;
        p->next = m;
        q->prev = m;
        v->alive = false;
        w->alive = false;

        invalidate_node(v);
        invalidate_node(w);
        invalidate_node(p);
        invalidate_node(q);
        compute_bisector(m);
        schedule_edge_event(p);
        schedule_edge_event(m);
        schedule_edge_event(q);
        recompute_all_splits();
    };

    const auto process_split_event = [&](Event* e) {
        Node* v = e->v;
        Node* a = e->a;
        Node* b = e->b;
        if (!v->alive || !a->alive || !b->alive || a->next != b || v->reflex == false) {
            return;
        }
        Node* p = v->prev;
        Node* q = v->next;

        const std::size_t s = add_node(e->point);
        push_arc(v->origin, s);

        Node* m1 = make_node(s);
        Node* m2 = make_node(s);
        m1->pos = e->point;
        m2->pos = e->point;
        m1->time = e->time;
        m2->time = e->time;
        m1->lav = v->lav;
        m2->lav = ++next_lav;

        m1->prev = a;
        m1->next = q;
        a->next = m1;
        q->prev = m1;
        m2->prev = p;
        m2->next = b;
        p->next = m2;
        b->prev = m2;
        v->alive = false;

        invalidate_node(v);
        invalidate_node(a);
        invalidate_node(b);
        invalidate_node(p);
        invalidate_node(q);
        compute_bisector(m1);
        compute_bisector(m2);
        // Relabel m2's sub-loop (m2 -> b -> ... -> p -> m2) to the new LAV id.
        // Before the split the whole loop shared v's LAV, so the new id is the
        // only thing that distinguishes the two resulting wavefronts.
        Node* cursor = m2;
        do {
            cursor->lav = m2->lav;
            cursor = cursor->next;
        } while (cursor != nullptr && cursor != m2);
        schedule_edge_event(a);
        schedule_edge_event(m1);
        schedule_edge_event(q);
        schedule_edge_event(p);
        schedule_edge_event(m2);
        schedule_edge_event(b);
        recompute_all_splits();
    };

    // Initialize the LAV and seed events.
    const std::size_t n = m_polygon.size();
    std::vector<Node*> ring(n);
    for (std::size_t i = 0; i < n; ++i) {
        ring[i] = make_node(i);
        ring[i]->pos = m_polygon[i];
        ring[i]->time = 0.0;
    }
    for (std::size_t i = 0; i < n; ++i) {
        ring[i]->prev = ring[(i + n - 1U) % n];
        ring[i]->next = ring[(i + 1U) % n];
        ring[i]->lav = 0;
    }
    for (Node* node : ring) compute_bisector(node);
    for (Node* node : ring) schedule_edge_event(node);
    for (Node* node : ring) if (node->reflex) schedule_split_event(node);

    while (!queue.empty() && total_events < event_cap) {
        Event* e = queue.top();
        queue.pop();
        if (!e->valid) continue;
        ++total_events;
        if (e->kind == Event::Kind::Edge) {
            process_edge_event(e);
        } else {
            process_split_event(e);
        }
    }

    return result;
}

} // namespace skeleton::detail
