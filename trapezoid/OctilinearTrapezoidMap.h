#ifndef TRAPEZOID_OCTILINEARTRAPEZOIDMAP_H
#define TRAPEZOID_OCTILINEARTRAPEZOIDMAP_H

#include "common/geometry/BasePoint.h"
#include "trapezoid/detail/DisjointSet.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace trapezoid {

// Trapezoidal map of an octilinear (VLSI-style) subdivision: every input edge
// is horizontal, vertical, or at +/- 45 degrees.  Built by a deterministic
// vertical sweep line in O(n log n).
//
// Vertical rays are shot from every endpoint; between two consecutive event
// x-coordinates the active non-vertical edges form a fixed stack of curves
// (edges never cross), so each band is a trapezoid whose top/bottom are
// input edges (horizontal or +/- 45 degree) and whose left/right sides are
// vertical walls.
//
// With integral coordinates every ordering predicate is an exact integer
// comparison: active curves evaluate to y = slope*x + intercept with
// slope in {-1, 0, 1}, so no floating point is used anywhere.
//
// Point location uses a uniform spatial hash (expected O(1) for uniformly
// distributed layouts); the decomposition itself does not depend on the
// distribution.  Face/band identity is carried by edge ids.
template <typename coordinate_type>
class OctilinearTrapezoidMap {
    static_assert(std::is_integral_v<coordinate_type>,
                  "OctilinearTrapezoidMap requires integral coordinates");

public:
    using Point = common::geometry::BasePoint<coordinate_type>;
    static constexpr std::size_t INVALID =
        std::numeric_limits<std::size_t>::max();

    struct Segment {
        Point start;
        Point end;
    };

    struct Cell {
        // Vertical sides.
        coordinate_type xmin = 0;
        coordinate_type xmax = 0;
        // Corners: lower/upper edge y at the left and right side.
        coordinate_type yLowerLeft = 0;
        coordinate_type yUpperLeft = 0;
        coordinate_type yLowerRight = 0;
        coordinate_type yUpperRight = 0;
        // Defining edges: curve index, kBottomSentinel (-1) or
        // kTopSentinel (-2).
        long long lowerEdge = 0;
        long long upperEdge = 0;
        std::size_t id = 0;
        std::size_t face = 0;
    };

    struct Vertex {
        Point point;
        std::size_t incidentHalfEdge = INVALID;
    };

    struct HalfEdge {
        std::size_t origin = INVALID;
        std::size_t twin = INVALID;
        std::size_t next = INVALID;
        std::size_t previous = INVALID;
        std::size_t face = INVALID;
        bool isInputEdge = false;
        bool isBoundingEdge = false;
    };

    struct Face {
        // INVALID for the exterior face.
        std::size_t cell = INVALID;
        // Original-subdivision face label shared by connected cells.
        std::size_t subdivisionFace = INVALID;
        std::size_t incidentHalfEdge = INVALID;
    };

    // margin < 0 selects an automatic bounding box margin.
    explicit OctilinearTrapezoidMap(
        const std::vector<Segment>& segments, long double margin = -1.0L) {
        if (segments.empty()) {
            throw std::invalid_argument("Segment set must not be empty");
        }
        classifySegments(segments);
        validateNoOverlap();
        computeBounds(segments, margin);
        buildEvents();
        sweep();
        labelFaces();
        buildDcel();
        buildQueryGrid();
    }

    const std::vector<Cell>& cells() const {
        return m_cells;
    }

    std::size_t faceCount() const {
        return m_face_count;
    }

    const std::vector<Vertex>& vertices() const {
        return m_vertices;
    }

    const std::vector<HalfEdge>& halfEdges() const {
        return m_half_edges;
    }

    const std::vector<Face>& faces() const {
        return m_faces;
    }

    std::size_t exteriorFace() const {
        return m_cells.size();
    }

    struct Bounds {
        coordinate_type xmin = 0;
        coordinate_type ymin = 0;
        coordinate_type xmax = 0;
        coordinate_type ymax = 0;
    };

    Bounds bounds() const {
        return {m_bounds_xmin, m_bounds_ymin, m_bounds_xmax, m_bounds_ymax};
    }

    // Returns the id of the cell containing (x, y), or -1 when the point is
    // outside the bounding box.  Boundary points go to one adjacent cell.
    std::int64_t locate(coordinate_type x, coordinate_type y) const {
        if (x < m_bounds_xmin || x > m_bounds_xmax ||
            y < m_bounds_ymin || y > m_bounds_ymax) {
            return -1;
        }
        const std::size_t ix = clampedIndex(x, m_bounds_xmin,
                                            m_cell_size_x, m_grid_nx);
        const std::size_t iy = clampedIndex(y, m_bounds_ymin,
                                            m_cell_size_y, m_grid_ny);
        for (const std::size_t id : m_grid[iy * m_grid_nx + ix]) {
            const Cell& cell = m_cells[id];
            if (x < cell.xmin || x > cell.xmax) {
                continue;
            }
            const coordinate_type y_lower = edgeY(cell.lowerEdge, x);
            const coordinate_type y_upper = edgeY(cell.upperEdge, x);
            if (y >= y_lower && y <= y_upper) {
                return static_cast<std::int64_t>(id);
            }
        }
        return -1; // unreachable for a point inside the bounding box
    }

    std::int64_t locateFace(coordinate_type x, coordinate_type y) const {
        const std::int64_t id = locate(x, y);
        return id < 0 ? -1
                      : static_cast<std::int64_t>(m_cells[id].face);
    }

private:
    using value_type = coordinate_type;

    static constexpr long long kBottom = -1;
    static constexpr long long kTop = -2;
    static constexpr long long kTransient = -3;

    struct Curve {
        value_type xmin = 0;
        value_type xmax = 0;
        long long slope = 0;       // -1, 0 or +1
        value_type intercept = 0; // y = slope * x + intercept
    };

    struct Vertical {
        value_type x = 0;
        value_type ymin = 0;
        value_type ymax = 0;
    };

    struct Event {
        std::vector<long long> starts;
        std::vector<long long> ends;
        std::vector<Vertical> verticals;
        std::vector<value_type> rayYs;
    };

    struct BandRef {
        long long lower = 0;
        long long upper = 0;
        value_type ylo = 0;
        value_type yhi = 0;
    };

    struct OpenCell {
        value_type xLeft = 0;
        std::size_t id = 0;
    };

    struct EventCell {
        std::size_t id = 0;
        value_type ylo = 0;
        value_type yhi = 0;
    };

    // Ordered set of active curves.  Ordering is evaluated at the current
    // event x; `side` selects the symbolic x - epsilon (before the event) or
    // x + epsilon (after), which shifts a slope-s curve by -s / +s.  This is
    // the standard endpoint-contact tie-break: the ordering of any pair of
    // curves both present on the same side never changes without a crossing,
    // and crossings are rejected.
    struct ActiveState {
        const std::vector<Curve>* curves = nullptr;
        value_type boundsYMin = 0;
        value_type boundsYMax = 0;
        value_type x = 0;
        int side = 1;          // -1 = before, +1 = after
        bool transientActive = false;
        value_type transientY = 0;
    };

    struct ActiveComparator {
        ActiveState* state = nullptr;

        // Order at x + side*epsilon:
        //   1. exact y at the event (differences are integers >= 1, so an
        //      infinitesimal perturbation can never reorder unequal ys);
        //   2. at a true endpoint contact (equal ys), side*slope decides
        //      which curve rises above the other toward the queried side;
        //   3. id as a final deterministic tie-break.
        bool operator()(long long a, long long b) const {
            const value_type ya = exactY(a);
            const value_type yb = exactY(b);
            if (ya != yb) {
                return ya < yb;
            }
            const int order_a = secondaryKey(a);
            const int order_b = secondaryKey(b);
            if (order_a != order_b) {
                return order_a < order_b;
            }
            return a < b;
        }

        int secondaryKey(long long id) const {
            if (id < 0) {
                return 0; // sentinels and the transient query point
            }
            return state->side * (*state->curves)[id].slope;
        }

        // True when curve id is strictly below y on the queried side, i.e.
        // its perturbed position is below the point.
        bool belowPoint(long long id, value_type y) const {
            const value_type exact = exactY(id);
            if (exact != y) {
                return exact < y;
            }
            return secondaryKey(id) < 0;
        }

        value_type exactY(long long id) const {
            if (id == kBottom) {
                return state->boundsYMin;
            }
            if (id == kTop) {
                return state->boundsYMax;
            }
            if (id == kTransient) {
                return state->transientY;
            }
            const Curve& curve = (*state->curves)[id];
            return static_cast<value_type>(
                static_cast<long long>(curve.slope) * state->x +
                curve.intercept);
        }
    };

    class ActiveSet {
    public:
        explicit ActiveSet(ActiveState* state)
            : m_state(state), m_set(ActiveComparator{state}) {}

        void clear() {
            m_set.clear();
            m_set.insert(kBottom);
            m_set.insert(kTop);
        }

        void setEvent(value_type x, int side) {
            m_state->x = x;
            m_state->side = side;
        }

        void insertCurve(long long id) {
            m_set.insert(id);
        }

        void eraseCurve(long long id) {
            m_set.erase(id);
        }

        // Curves immediately below/above the horizontal line y.  Either may
        // be absent when y is outside the active range.
        struct Neighbors {
            bool hasBelow = false;
            bool hasAbove = false;
            long long below = kBottom;
            long long above = kTop;
        };

        Neighbors neighbors(value_type y) const {
            m_state->transientY = y;
            m_state->transientActive = true;
            const auto inserted = m_set.insert(kTransient);
            auto it = inserted.first;
            Neighbors result;
            if (it != m_set.begin()) {
                result.hasBelow = true;
                result.below = *std::prev(it);
            }
            auto next_it = std::next(it);
            if (next_it != m_set.end()) {
                result.hasAbove = true;
                result.above = *next_it;
            }
            m_set.erase(it);
            m_state->transientActive = false;
            return result;
        }

        // Predecessor/successor of an active curve in the current ordering.
        bool predecessor(long long id, long long& result) const {
            auto it = m_set.find(id);
            if (it == m_set.end() || it == m_set.begin()) {
                return false;
            }
            result = *std::prev(it);
            return true;
        }

        bool successor(long long id, long long& result) const {
            auto it = m_set.find(id);
            if (it == m_set.end()) {
                return false;
            }
            auto next_it = std::next(it);
            if (next_it == m_set.end()) {
                return false;
            }
            result = *next_it;
            return true;
        }

        using Iterator = typename std::set<long long,
                                          ActiveComparator>::const_iterator;

        // First active curve strictly above y (symbolic side ordering).
        Iterator firstAbove(value_type y) const {
            const Neighbors near = neighbors(y);
            if (!near.hasAbove) {
                return m_set.end();
            }
            return m_set.find(near.above);
        }

        Iterator end() const {
            return m_set.end();
        }

        Iterator next(Iterator it) const {
            return std::next(it);
        }

        value_type exact(long long id) const {
            return ActiveComparator{m_state}.exactY(id);
        }

        bool belowPoint(long long id, value_type y) const {
            return ActiveComparator{m_state}.belowPoint(id, y);
        }

    private:
        ActiveState* m_state;
        mutable std::set<long long, ActiveComparator> m_set;
    };

    value_type edgeY(long long edge, value_type x) const {
        if (edge == kBottom) {
            return m_bounds_ymin;
        }
        if (edge == kTop) {
            return m_bounds_ymax;
        }
        const Curve& curve = m_curves[edge];
        return static_cast<value_type>(
            static_cast<long long>(curve.slope) * x + curve.intercept);
    }

    void classifySegments(const std::vector<Segment>& segments) {
        m_curves.reserve(segments.size());
        m_verticals.reserve(segments.size());
        for (const Segment& segment : segments) {
            const long long x1 = segment.start.x();
            const long long y1 = segment.start.y();
            const long long x2 = segment.end.x();
            const long long y2 = segment.end.y();
            if (x1 == x2 && y1 == y2) {
                throw std::invalid_argument("Zero-length segments are rejected");
            }
            if (x1 == x2) {
                m_verticals.push_back({static_cast<value_type>(x1),
                                       static_cast<value_type>(std::min(y1, y2)),
                                       static_cast<value_type>(std::max(y1, y2))});
            } else {
                const long long dx = x2 - x1;
                const long long dy = y2 - y1;
                long long slope = 0;
                if (dy == 0) {
                    slope = 0;
                } else if (dx == dy || dx == -dy) {
                    slope = dy / dx; // +1 or -1
                } else {
                    throw std::invalid_argument(
                        "Only horizontal, vertical and +/-45 degree segments "
                        "are accepted");
                }
                const long long left_x = std::min(x1, x2);
                const long long right_x = std::max(x1, x2);
                const long long left_y = (x1 < x2) ? y1 : y2;
                m_curves.push_back({
                    static_cast<value_type>(left_x),
                    static_cast<value_type>(right_x),
                    slope,
                    static_cast<value_type>(left_y - slope * left_x)});
            }
        }
        if (m_curves.empty() && m_verticals.empty()) {
            throw std::invalid_argument("Segment set must not be empty");
        }
    }

    void validateNoOverlap() const {
        std::vector<Curve> curves = m_curves;
        std::sort(curves.begin(), curves.end(), [](const Curve& a,
                                                  const Curve& b) {
            if (a.slope != b.slope) {
                return a.slope < b.slope;
            }
            if (a.intercept != b.intercept) {
                return a.intercept < b.intercept;
            }
            return a.xmin < b.xmin;
        });
        for (std::size_t i = 1; i < curves.size(); ++i) {
            const Curve& previous = curves[i - 1];
            const Curve& current = curves[i];
            if (previous.slope == current.slope &&
                previous.intercept == current.intercept &&
                current.xmin < previous.xmax) {
                throw std::invalid_argument(
                    "Overlapping collinear segments");
            }
        }
        std::vector<Vertical> verticals = m_verticals;
        std::sort(verticals.begin(), verticals.end(),
                  [](const Vertical& a, const Vertical& b) {
                      if (a.x != b.x) {
                          return a.x < b.x;
                      }
                      return a.ymin < b.ymin;
                  });
        for (std::size_t i = 1; i < verticals.size(); ++i) {
            if (verticals[i - 1].x == verticals[i].x &&
                verticals[i].ymin < verticals[i - 1].ymax) {
                throw std::invalid_argument(
                    "Overlapping collinear vertical segments");
            }
        }
    }

    void computeBounds(const std::vector<Segment>& segments,
                       long double configured_margin) {
        value_type min_x = segments.front().start.x();
        value_type max_x = min_x;
        value_type min_y = segments.front().start.y();
        value_type max_y = min_y;
        for (const Segment& segment : segments) {
            for (const Point* p : {&segment.start, &segment.end}) {
                min_x = std::min(min_x, p->x());
                max_x = std::max(max_x, p->x());
                min_y = std::min(min_y, p->y());
                max_y = std::max(max_y, p->y());
            }
        }
        long double span = std::max<long double>(
            static_cast<long double>(max_x) - min_x,
            static_cast<long double>(max_y) - min_y);
        if (span == 0.0L) {
            span = 1.0L;
        }
        long double delta = configured_margin;
        if (delta <= 0.0L) {
            delta = std::max<long double>(1.0L, span * 0.02L);
        }
        m_bounds_xmin = static_cast<value_type>(
            static_cast<long double>(min_x) - delta);
        m_bounds_xmax = static_cast<value_type>(
            static_cast<long double>(max_x) + delta);
        m_bounds_ymin = static_cast<value_type>(
            static_cast<long double>(min_y) - delta);
        m_bounds_ymax = static_cast<value_type>(
            static_cast<long double>(max_y) + delta);
    }

    void buildEvents() {
        for (long long id = 0;
             id < static_cast<long long>(m_curves.size()); ++id) {
            const Curve& curve = m_curves[id];
            m_events[curve.xmin].starts.push_back(id);
            m_events[curve.xmin].rayYs.push_back(
                edgeY(id, curve.xmin));
            m_events[curve.xmax].ends.push_back(id);
            m_events[curve.xmax].rayYs.push_back(
                edgeY(id, curve.xmax));
        }
        for (const Vertical& vertical : m_verticals) {
            Event& event = m_events[vertical.x];
            event.verticals.push_back(vertical);
            event.rayYs.push_back(vertical.ymin);
            event.rayYs.push_back(vertical.ymax);
        }
        for (auto& entry : m_events) {
            Event& event = entry.second;
            std::sort(event.rayYs.begin(), event.rayYs.end());
            event.rayYs.erase(
                std::unique(event.rayYs.begin(), event.rayYs.end()),
                event.rayYs.end());
            std::sort(event.verticals.begin(), event.verticals.end(),
                      [](const Vertical& a, const Vertical& b) {
                          return a.ymin < b.ymin;
                      });
        }
    }

    // Proper intersection of two curves at an x strictly greater than
    // `event_x` and strictly inside both x-ranges.  Exact rational test.
    bool crossAfter(long long a, long long b, value_type event_x) const {
        const Curve& ca = m_curves[a];
        const Curve& cb = m_curves[b];
        if (ca.slope == cb.slope) {
            return false; // parallel; overlap rejected separately
        }
        const __int128 denominator =
            static_cast<__int128>(ca.slope) - cb.slope;
        const __int128 numerator =
            static_cast<__int128>(cb.intercept) - ca.intercept;
        const __int128 x_event = event_x;
        const __int128 limit = std::min(
            static_cast<__int128>(ca.xmax),
            static_cast<__int128>(cb.xmax));
        // x_intersect = numerator / denominator, denominator is +-1 or +-2.
        if (denominator > 0) {
            return numerator > denominator * x_event &&
                   numerator < denominator * limit;
        }
        return numerator < denominator * x_event &&
               numerator > denominator * limit;
    }

    std::vector<BandRef> collectTouchedBands(const ActiveSet& active,
                                             const Event& event) const {
        std::set<std::pair<long long, long long>> seen;
        std::vector<BandRef> result;

        auto add_band = [&](long long lower, long long upper) {
            if (lower == upper) {
                return;
            }
            if (seen.insert({lower, upper}).second) {
                result.push_back({lower, upper, active.exact(lower),
                                  active.exact(upper)});
            }
        };

        for (const value_type ray_y : event.rayYs) {
            const typename ActiveSet::Neighbors near =
                active.neighbors(ray_y);

            // Collect the whole fan of curves passing exactly through the
            // ray point (several edges can meet at one vertex).  They form a
            // tie group around the transient position in both directions.
            std::vector<long long> group;
            long long cursor = near.hasBelow ? near.below : kBottom;
            while (cursor >= 0 && active.exact(cursor) == ray_y) {
                group.push_back(cursor);
                if (!active.predecessor(cursor, cursor)) {
                    break;
                }
            }
            std::reverse(group.begin(), group.end());
            cursor = near.hasAbove ? near.above : kTop;
            while (cursor >= 0 && active.exact(cursor) == ray_y) {
                group.push_back(cursor);
                if (!active.successor(cursor, cursor)) {
                    break;
                }
            }

            if (group.empty()) {
                // Endpoint in the interior of a band: the ray crosses it.
                if (near.hasBelow && near.hasAbove) {
                    add_band(near.below, near.above);
                }
                continue;
            }

            // Bands between the fan and whatever lies outside it.
            long long outside = kBottom;
            if (active.predecessor(group.front(), outside)) {
                add_band(outside, group.front());
            }
            for (std::size_t k = 1; k < group.size(); ++k) {
                add_band(group[k - 1], group[k]);
            }
            if (active.successor(group.back(), outside)) {
                add_band(group.back(), outside);
            }
        }

        for (const Vertical& vertical : event.verticals) {
            const typename ActiveSet::Neighbors around =
                active.neighbors(vertical.ymin);
            if (!around.hasAbove) {
                continue;
            }
            // Band straddling the lower endpoint of the wall.
            if (around.hasBelow &&
                active.exact(around.below) < vertical.ymax) {
                add_band(around.below, around.above);
            }
            auto it = active.firstAbove(vertical.ymin);
            while (it != active.end() &&
                   active.belowPoint(*it, vertical.ymax)) {
                auto next_it = active.next(it);
                if (next_it != active.end()) {
                    add_band(*it, *next_it);
                }
                it = next_it;
            }
        }

        std::sort(result.begin(), result.end(),
                  [](const BandRef& a, const BandRef& b) {
                      if (a.ylo != b.ylo) {
                          return a.ylo < b.ylo;
                      }
                      return a.yhi < b.yhi;
                  });
        return result;
    }

    void sweep() {
        ActiveState state;
        state.curves = &m_curves;
        state.boundsYMin = m_bounds_ymin;
        state.boundsYMax = m_bounds_ymax;
        ActiveSet active(&state);
        active.clear();

        std::map<std::pair<long long, long long>, OpenCell> open;
        openCell(open, {kBottom, kTop}, m_bounds_xmin);

        for (const auto& entry : m_events) {
            const value_type x = entry.first;
            const Event& event = entry.second;

            // ---- left side: close touched bands ----
            active.setEvent(x, -1);
            const std::vector<BandRef> left_bands =
                collectTouchedBands(active, event);
            std::vector<EventCell> closed;
            closed.reserve(left_bands.size());
            for (const BandRef& band : left_bands) {
                auto open_it = open.find({band.lower, band.upper});
                if (open_it == open.end()) {
                    continue;
                }
                Cell& cell = m_cells[open_it->second.id];
                cell.xmax = x;
                cell.yLowerRight = edgeY(band.lower, x);
                cell.yUpperRight = edgeY(band.upper, x);
                closed.push_back(
                    {open_it->second.id, band.ylo, band.yhi});
                open.erase(open_it);
            }

            // ---- advance the active set ----
            // Erase first, then inspect the new adjacency: the removed curve
            // must not appear as its own neighbor under the perturbation.
            for (const long long id : event.ends) {
                const Curve& curve = m_curves[id];
                const value_type end_y = edgeY(id, curve.xmax);
                active.eraseCurve(id);
                const typename ActiveSet::Neighbors near_before =
                    active.neighbors(end_y);
                if (near_before.hasBelow && near_before.hasAbove &&
                    near_before.below >= 0 && near_before.above >= 0 &&
                    crossAfter(near_before.below, near_before.above, x)) {
                    throw std::invalid_argument(
                        "Properly crossing non-vertical segments");
                }
            }

            active.setEvent(x, +1);

            for (const long long id : event.starts) {
                // Inspect future adjacency before the curve joins the set.
                const Curve& curve = m_curves[id];
                const value_type start_y = edgeY(id, curve.xmin);
                const typename ActiveSet::Neighbors near_after =
                    active.neighbors(start_y);
                if (near_after.hasBelow && near_after.below >= 0 &&
                    crossAfter(near_after.below, id, x)) {
                    throw std::invalid_argument(
                        "Properly crossing non-vertical segments");
                }
                if (near_after.hasAbove && near_after.above >= 0 &&
                    crossAfter(id, near_after.above, x)) {
                    throw std::invalid_argument(
                        "Properly crossing non-vertical segments");
                }
                active.insertCurve(id);
            }

            // Vertical walls vs active curves: exact y strictly inside.
            for (const Vertical& vertical : event.verticals) {
                auto it = active.firstAbove(vertical.ymin);
                while (it != active.end() &&
                       active.belowPoint(*it, vertical.ymax)) {
                    const long long id = *it;
                    const value_type y_here = active.exact(id);
                    if (y_here > vertical.ymin && y_here < vertical.ymax) {
                        throw std::invalid_argument(
                            "Properly crossing vertical and non-vertical "
                            "segments");
                    }
                    it = active.next(it);
                }
            }

            // ---- right side: open touched bands ----
            const std::vector<BandRef> right_bands =
                collectTouchedBands(active, event);
            std::vector<EventCell> opened;
            opened.reserve(right_bands.size());
            for (const BandRef& band : right_bands) {
                const std::pair<long long, long long> key = {band.lower,
                                                             band.upper};
                if (open.find(key) != open.end()) {
                    continue;
                }
                const std::size_t id = openCell(open, key, x);
                opened.push_back({id, band.ylo, band.yhi});
            }

            // Artificial rays never separate faces; only input vertical
            // segments do.
            connectAcrossEvent(closed, opened, event);
        }

        for (auto& entry : open) {
            Cell& cell = m_cells[entry.second.id];
            cell.xmax = m_bounds_xmax;
            cell.yLowerRight = edgeY(cell.lowerEdge, m_bounds_xmax);
            cell.yUpperRight = edgeY(cell.upperEdge, m_bounds_xmax);
        }
    }

    std::size_t openCell(
        std::map<std::pair<long long, long long>, OpenCell>& open,
        std::pair<long long, long long> key, value_type x_left) {
        const std::size_t id = m_cells.size();
        Cell cell;
        cell.xmin = x_left;
        cell.xmax = x_left;
        cell.lowerEdge = key.first;
        cell.upperEdge = key.second;
        cell.yLowerLeft = edgeY(key.first, x_left);
        cell.yUpperLeft = edgeY(key.second, x_left);
        cell.yLowerRight = cell.yLowerLeft;
        cell.yUpperRight = cell.yUpperLeft;
        cell.id = id;
        m_cells.push_back(cell);
        m_dsu.add();
        open.emplace(std::move(key), OpenCell{x_left, id});
        return id;
    }

    void connectAcrossEvent(const std::vector<EventCell>& closed,
                            const std::vector<EventCell>& opened,
                            const Event& event) {
        std::size_t i = 0;
        std::size_t j = 0;
        while (i < closed.size() && j < opened.size()) {
            const EventCell& left = closed[i];
            const EventCell& right = opened[j];
            const value_type lo = std::max(left.ylo, right.ylo);
            const value_type hi = std::min(left.yhi, right.yhi);
            if (hi > lo) {
                // Sum the blocked portions of [lo, hi) exactly.
                __int128 blocked = 0;
                for (const Vertical& vertical : event.verticals) {
                    const value_type blo = std::max(lo, vertical.ymin);
                    const value_type bhi = std::min(hi, vertical.ymax);
                    if (bhi > blo) {
                        blocked += static_cast<__int128>(bhi) - blo;
                    }
                }
                if (blocked < static_cast<__int128>(hi) - lo) {
                    m_dsu.unite(left.id, right.id);
                }
            }
            if (left.yhi < right.yhi) {
                ++i;
            } else if (right.yhi < left.yhi) {
                ++j;
            } else {
                ++i;
                ++j;
            }
        }
    }

    void labelFaces() {
        std::map<std::size_t, std::size_t> face_by_root;
        m_face_count = 0;
        for (Cell& cell : m_cells) {
            const std::size_t root = m_dsu.find(cell.id);
            auto found = face_by_root.find(root);
            if (found == face_by_root.end()) {
                found = face_by_root.emplace(root, m_face_count++).first;
            }
            cell.face = found->second;
        }
    }

    struct DcelLine {
        int kind = 0;
        value_type offset = 0;

        bool operator<(const DcelLine& other) const {
            return kind < other.kind ||
                   (kind == other.kind && offset < other.offset);
        }
    };

    struct DcelRawEdge {
        std::pair<value_type, value_type> start;
        std::pair<value_type, value_type> end;
        std::size_t face = INVALID;
        bool isInputEdge = false;
        bool isBoundingEdge = false;
    };

    static DcelLine dcelLine(
        const std::pair<value_type, value_type>& start,
        const std::pair<value_type, value_type>& end) {
        if (start.first == end.first) {
            return {0, start.first};
        }
        if (start.second == end.second) {
            return {1, start.second};
        }
        const value_type slope =
            (end.second > start.second) == (end.first > start.first) ? 1 : -1;
        return {slope > 0 ? 2 : 3,
                static_cast<value_type>(start.second - slope * start.first)};
    }

    static value_type dcelParameter(
        const DcelLine& line,
        const std::pair<value_type, value_type>& point) {
        return line.kind == 0 ? point.second : point.first;
    }

    static std::pair<value_type, value_type> dcelPoint(
        const DcelLine& line, value_type parameter) {
        if (line.kind == 0) {
            return {line.offset, parameter};
        }
        if (line.kind == 1) {
            return {parameter, line.offset};
        }
        const value_type slope = line.kind == 2 ? 1 : -1;
        return {parameter,
                static_cast<value_type>(slope * parameter + line.offset)};
    }

    bool isInputVertical(value_type x, value_type ymin,
                         value_type ymax) const {
        for (const Vertical& vertical : m_verticals) {
            if (vertical.x == x && vertical.ymin <= ymin &&
                vertical.ymax >= ymax) {
                return true;
            }
        }
        return false;
    }

    void buildDcel() {
        using CoordinatePair = std::pair<value_type, value_type>;
        std::vector<DcelRawEdge> raw_edges;
        std::map<DcelLine, std::set<value_type>> points_by_line;

        auto add_raw = [&](CoordinatePair start, CoordinatePair end,
                           std::size_t face, bool input, bool bounding) {
            if (start == end) {
                return;
            }
            const DcelLine line = dcelLine(start, end);
            points_by_line[line].insert(dcelParameter(line, start));
            points_by_line[line].insert(dcelParameter(line, end));
            raw_edges.push_back(
                {start, end, face, input, bounding});
        };

        for (const Cell& cell : m_cells) {
            const CoordinatePair lower_left{cell.xmin, cell.yLowerLeft};
            const CoordinatePair lower_right{cell.xmax, cell.yLowerRight};
            const CoordinatePair upper_right{cell.xmax, cell.yUpperRight};
            const CoordinatePair upper_left{cell.xmin, cell.yUpperLeft};
            add_raw(lower_left, lower_right, cell.id, cell.lowerEdge >= 0,
                    cell.lowerEdge == kBottom);
            add_raw(lower_right, upper_right, cell.id, false,
                    cell.xmax == m_bounds_xmax);
            add_raw(upper_right, upper_left, cell.id, cell.upperEdge >= 0,
                    cell.upperEdge == kTop);
            add_raw(upper_left, lower_left, cell.id, false,
                    cell.xmin == m_bounds_xmin);
        }

        m_faces.reserve(m_cells.size() + 1);
        for (const Cell& cell : m_cells) {
            m_faces.push_back({cell.id, cell.face, INVALID});
        }
        m_faces.push_back({INVALID, INVALID, INVALID});

        std::map<CoordinatePair, std::size_t> vertex_by_point;
        auto vertex_id = [&](const CoordinatePair& point) {
            const auto found = vertex_by_point.find(point);
            if (found != vertex_by_point.end()) {
                return found->second;
            }
            const std::size_t id = m_vertices.size();
            vertex_by_point.emplace(point, id);
            m_vertices.push_back({Point(point.first, point.second), INVALID});
            return id;
        };

        std::map<std::pair<std::size_t, std::size_t>, std::size_t>
            directed_edge;
        std::vector<std::vector<std::size_t>> face_edges(m_cells.size());
        for (const DcelRawEdge& raw : raw_edges) {
            const DcelLine line = dcelLine(raw.start, raw.end);
            const value_type first = dcelParameter(line, raw.start);
            const value_type last = dcelParameter(line, raw.end);
            std::vector<value_type> parameters;
            const value_type low = std::min(first, last);
            const value_type high = std::max(first, last);
            const auto& line_points = points_by_line[line];
            for (auto it = line_points.lower_bound(low);
                 it != line_points.end() && *it <= high; ++it) {
                parameters.push_back(*it);
            }
            if (first > last) {
                std::reverse(parameters.begin(), parameters.end());
            }

            for (std::size_t i = 1; i < parameters.size(); ++i) {
                const CoordinatePair start = dcelPoint(line, parameters[i - 1]);
                const CoordinatePair end = dcelPoint(line, parameters[i]);
                const std::size_t origin = vertex_id(start);
                const std::size_t destination = vertex_id(end);
                const bool vertical = start.first == end.first;
                const bool input =
                    raw.isInputEdge ||
                    (vertical &&
                     isInputVertical(start.first,
                                     std::min(start.second, end.second),
                                     std::max(start.second, end.second)));
                const std::size_t edge = m_half_edges.size();
                m_half_edges.push_back(
                    {origin, INVALID, INVALID, INVALID, raw.face, input,
                     raw.isBoundingEdge});
                if (!directed_edge.emplace(
                        std::make_pair(origin, destination), edge).second) {
                    throw std::logic_error(
                        "Invalid trapezoidal subdivision topology");
                }
                const auto twin = directed_edge.find(
                    std::make_pair(destination, origin));
                if (twin != directed_edge.end()) {
                    m_half_edges[edge].twin = twin->second;
                    m_half_edges[twin->second].twin = edge;
                    m_half_edges[edge].isInputEdge |=
                        m_half_edges[twin->second].isInputEdge;
                    m_half_edges[twin->second].isInputEdge =
                        m_half_edges[edge].isInputEdge;
                    m_half_edges[edge].isBoundingEdge |=
                        m_half_edges[twin->second].isBoundingEdge;
                    m_half_edges[twin->second].isBoundingEdge =
                        m_half_edges[edge].isBoundingEdge;
                }
                if (m_vertices[origin].incidentHalfEdge == INVALID) {
                    m_vertices[origin].incidentHalfEdge = edge;
                }
                if (m_faces[raw.face].incidentHalfEdge == INVALID) {
                    m_faces[raw.face].incidentHalfEdge = edge;
                }
                face_edges[raw.face].push_back(edge);
            }
        }

        for (const auto& edges : face_edges) {
            for (std::size_t i = 0; i < edges.size(); ++i) {
                const std::size_t edge = edges[i];
                const std::size_t next = edges[(i + 1) % edges.size()];
                m_half_edges[edge].next = next;
                m_half_edges[next].previous = edge;
            }
        }

        const std::size_t exterior = exteriorFace();
        const std::size_t interior_edge_count = m_half_edges.size();
        std::vector<std::size_t> exterior_edges;
        for (std::size_t edge = 0; edge < interior_edge_count; ++edge) {
            if (m_half_edges[edge].twin != INVALID) {
                continue;
            }
            const std::size_t destination =
                m_half_edges[m_half_edges[edge].next].origin;
            const std::size_t twin = m_half_edges.size();
            m_half_edges.push_back(
                {destination, edge, INVALID, INVALID, exterior, false, true});
            m_half_edges[edge].twin = twin;
            exterior_edges.push_back(twin);
            if (m_vertices[destination].incidentHalfEdge == INVALID) {
                m_vertices[destination].incidentHalfEdge = twin;
            }
        }

        std::map<std::size_t, std::size_t> exterior_from_origin;
        for (const std::size_t edge : exterior_edges) {
            if (!exterior_from_origin
                     .emplace(m_half_edges[edge].origin, edge)
                     .second) {
                throw std::logic_error(
                    "Invalid exterior boundary topology");
            }
        }
        for (const std::size_t edge : exterior_edges) {
            const std::size_t destination =
                m_half_edges[m_half_edges[edge].twin].origin;
            const auto next = exterior_from_origin.find(destination);
            if (next == exterior_from_origin.end()) {
                throw std::logic_error(
                    "Open exterior boundary in trapezoidal subdivision");
            }
            m_half_edges[edge].next = next->second;
            m_half_edges[next->second].previous = edge;
        }
        if (!exterior_edges.empty()) {
            m_faces[exterior].incidentHalfEdge = exterior_edges.front();
        }
    }

    void buildQueryGrid() {
        const long double span_x =
            static_cast<long double>(m_bounds_xmax) - m_bounds_xmin;
        const long double span_y =
            static_cast<long double>(m_bounds_ymax) - m_bounds_ymin;
        const long double target_cells =
            std::max<long double>(4.0L, 4.0L * m_cells.size());
        m_cell_size_x = static_cast<value_type>(
            std::sqrt(span_x * span_y / target_cells));
        if (static_cast<long double>(m_cell_size_x) <= 0.0L) {
            m_cell_size_x = static_cast<value_type>(1);
        }
        m_cell_size_y = m_cell_size_x;
        m_grid_nx = std::max<std::size_t>(
            1, static_cast<std::size_t>(
                   std::ceil(span_x / static_cast<long double>(
                                         m_cell_size_x))));
        m_grid_ny = std::max<std::size_t>(
            1, static_cast<std::size_t>(
                   std::ceil(span_y / static_cast<long double>(
                                         m_cell_size_y))));
        m_grid.assign(m_grid_nx * m_grid_ny, {});
        for (const Cell& cell : m_cells) {
            const value_type cell_ymin =
                std::min({cell.yLowerLeft, cell.yUpperLeft,
                          cell.yLowerRight, cell.yUpperRight});
            const value_type cell_ymax =
                std::max({cell.yLowerLeft, cell.yUpperLeft,
                          cell.yLowerRight, cell.yUpperRight});
            const std::size_t ix0 = clampedIndex(
                cell.xmin, m_bounds_xmin, m_cell_size_x, m_grid_nx);
            const std::size_t ix1 = clampedIndex(
                cell.xmax, m_bounds_xmin, m_cell_size_x, m_grid_nx);
            const std::size_t iy0 = clampedIndex(
                cell_ymin, m_bounds_ymin, m_cell_size_y, m_grid_ny);
            const std::size_t iy1 = clampedIndex(
                cell_ymax, m_bounds_ymin, m_cell_size_y, m_grid_ny);
            for (std::size_t iy = iy0; iy <= iy1; ++iy) {
                for (std::size_t ix = ix0; ix <= ix1; ++ix) {
                    m_grid[iy * m_grid_nx + ix].push_back(cell.id);
                }
            }
        }
    }

    static std::size_t clampedIndex(value_type value, value_type origin,
                                    value_type cell_size, std::size_t count) {
        long double index_f =
            (static_cast<long double>(value) - origin) /
            static_cast<long double>(cell_size);
        if (index_f < 0.0L) {
            index_f = 0.0L;
        }
        std::size_t index = static_cast<std::size_t>(std::floor(index_f));
        return std::min(index, count - 1);
    }

    std::vector<Curve> m_curves;
    std::vector<Vertical> m_verticals;
    std::map<value_type, Event> m_events;

    value_type m_bounds_xmin = 0;
    value_type m_bounds_xmax = 0;
    value_type m_bounds_ymin = 0;
    value_type m_bounds_ymax = 0;

    std::vector<Cell> m_cells;
    detail::DisjointSet m_dsu;
    std::size_t m_face_count = 0;
    std::vector<Vertex> m_vertices;
    std::vector<HalfEdge> m_half_edges;
    std::vector<Face> m_faces;

    value_type m_cell_size_x = 1;
    value_type m_cell_size_y = 1;
    std::size_t m_grid_nx = 1;
    std::size_t m_grid_ny = 1;
    std::vector<std::vector<std::size_t>> m_grid;
};

} // namespace trapezoid

#endif // TRAPEZOID_OCTILINEARTRAPEZOIDMAP_H
