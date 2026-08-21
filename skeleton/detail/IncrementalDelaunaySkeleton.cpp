#include "skeleton/detail/IncrementalDelaunaySkeleton.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace skeleton::detail {

namespace {

// 12-bit Hilbert curve key for a point on a 4096x4096 grid. Gives a
// cache- and locality-friendly insertion order that keeps the walking search
// near the previous triangle (expected O(log n) per insertion).
std::uint32_t hilbertKey(std::uint32_t x, std::uint32_t y) {
    constexpr std::uint32_t n = 4096; // grid side, must be a power of two
    std::uint32_t d = 0;
    for (std::uint32_t s = n / 2; s > 0; s /= 2) {
        const std::uint32_t rx = (x & s) ? 1U : 0U;
        const std::uint32_t ry = (y & s) ? 1U : 0U;
        d += s * s * ((3U * rx) ^ ry);
        if (ry == 0U) {
            if (rx == 1U) {
                x = (n - 1U) - x;
                y = (n - 1U) - y;
            }
            std::swap(x, y);
        }
    }
    return d;
}

} // namespace

IncrementalDelaunayBuilder::IncrementalDelaunayBuilder(std::vector<Point> sites)
    : m_sites(std::move(sites)) {
    double scale = 1.0;
    for (const Point& site : m_sites) {
        scale = std::max(scale, std::abs(static_cast<double>(site.x())));
        scale = std::max(scale, std::abs(static_cast<double>(site.y())));
    }
    m_scale = std::max(scale, 1.0);
    m_areaEpsilon = m_scale * m_scale * 1e-12;
}

double IncrementalDelaunayBuilder::orientation(std::size_t a, std::size_t b,
                                               std::size_t c) const {
    return common::cross(m_sites[b] - m_sites[a], m_sites[c] - m_sites[a]);
}

bool IncrementalDelaunayBuilder::inCircumcircle(
    std::size_t a, std::size_t b, std::size_t c, std::size_t p) const {
    const Point da = m_sites[a] - m_sites[p];
    const Point db = m_sites[b] - m_sites[p];
    const Point dc = m_sites[c] - m_sites[p];
    const double adx = static_cast<double>(da.x());
    const double ady = static_cast<double>(da.y());
    const double bdx = static_cast<double>(db.x());
    const double bdy = static_cast<double>(db.y());
    const double cdx = static_cast<double>(dc.x());
    const double cdy = static_cast<double>(dc.y());
    const double a2 = adx * adx + ady * ady;
    const double b2 = bdx * bdx + bdy * bdy;
    const double c2 = cdx * cdx + cdy * cdy;
    const double determinant =
        a2 * (bdx * cdy - bdy * cdx) -
        b2 * (adx * cdy - ady * cdx) +
        c2 * (adx * bdy - ady * bdx);
    // For a CCW triangle the determinant is positive iff p is inside the
    // circumcircle; for a CW triangle the sign flips.
    if (orientation(a, b, c) < -m_areaEpsilon) return determinant < -m_areaEpsilon;
    return determinant > m_areaEpsilon;
}

std::size_t IncrementalDelaunayBuilder::addTriangle(std::size_t a, std::size_t b,
                                                    std::size_t c) {
    m_triangles.push_back({{a, b, c}, {kNone, kNone, kNone}, true});
    return m_triangles.size() - 1U;
}

std::size_t IncrementalDelaunayBuilder::locate(std::size_t site) const {
    if (m_triangles.empty()) return kNone;

    // Visibility walking from the hint: move across any edge that the site
    // lies outside of until the site is inside the current triangle.
    auto walk_from = [&](std::size_t start) -> std::size_t {
        std::size_t current = start;
        std::size_t guard = 0;
        const std::size_t max_steps = m_triangles.size() * 2U + 8U;
        while (guard++ < max_steps) {
            const Triangle& triangle = m_triangles[current];
            bool moved = false;
            for (int edge = 0; edge < 3; ++edge) {
                const std::size_t v0 = triangle.vertices[(edge + 1U) % 3U];
                const std::size_t v1 = triangle.vertices[(edge + 2U) % 3U];
                if (triangle.neighbors[edge] == kNone) continue;
                const double cross = common::cross(m_sites[v1] - m_sites[v0],
                                                  m_sites[site] - m_sites[v0]);
                if (cross < -m_areaEpsilon) {
                    current = triangle.neighbors[edge];
                    moved = true;
                    break;
                }
            }
            if (!moved) return current;
        }
        return current;
    };

    std::size_t found = walk_from(m_hint);
    // Fallback: linear scan if walking failed to converge numerically.
    const Triangle& found_triangle = m_triangles[found];
    bool inside = true;
    for (int edge = 0; edge < 3; ++edge) {
        const std::size_t v0 = found_triangle.vertices[(edge + 1U) % 3U];
        const std::size_t v1 = found_triangle.vertices[(edge + 2U) % 3U];
        if (common::cross(m_sites[v1] - m_sites[v0],
                          m_sites[site] - m_sites[v0]) < -m_areaEpsilon) {
            inside = false;
            break;
        }
    }
    if (inside) return found;

    for (std::size_t t = 0; t < m_triangles.size(); ++t) {
        if (!m_triangles[t].alive) continue;
        const Triangle& triangle = m_triangles[t];
        bool ok = true;
        for (int edge = 0; edge < 3; ++edge) {
            const std::size_t v0 = triangle.vertices[(edge + 1U) % 3U];
            const std::size_t v1 = triangle.vertices[(edge + 2U) % 3U];
            if (common::cross(m_sites[v1] - m_sites[v0],
                              m_sites[site] - m_sites[v0]) < -m_areaEpsilon) {
                ok = false;
                break;
            }
        }
        if (ok) return t;
    }
    return found;
}

void IncrementalDelaunayBuilder::insertSite(std::size_t site) {
    const std::size_t start = locate(site);
    if (start == kNone) return;

    // Timestamp caches sized to the current triangle count let each insertion
    // run in time proportional to the cavity size rather than O(triangles).
    const std::size_t count = m_triangles.size();
    if (m_visit_stamp.size() < count) {
        m_visit_stamp.resize(count * 2U + 8U, 0);
        m_bad_stamp.resize(count * 2U + 8U, 0);
    }
    const std::size_t visit_gen = ++m_visit_counter;

    // Cavity: triangles whose circumcircle contains the site. The cavity of a
    // Delaunay triangulation around an inserted point is connected and
    // star-shaped, so its boundary is a single simple loop.
    std::vector<std::size_t> bad_triangles;
    std::vector<std::size_t> stack{start};
    m_visit_stamp[start] = visit_gen;
    while (!stack.empty()) {
        const std::size_t t = stack.back();
        stack.pop_back();
        const Triangle& triangle = m_triangles[t];
        if (!inCircumcircle(triangle.vertices[0], triangle.vertices[1],
                            triangle.vertices[2], site)) {
            continue;
        }
        bad_triangles.push_back(t);
        for (int edge = 0; edge < 3; ++edge) {
            const std::size_t neighbor = triangle.neighbors[edge];
            if (neighbor != kNone && m_visit_stamp[neighbor] != visit_gen) {
                m_visit_stamp[neighbor] = visit_gen;
                stack.push_back(neighbor);
            }
        }
    }

    if (bad_triangles.empty()) return;

    // Cavity membership is needed separately from the BFS visited set: a
    // neighbor can be visited but not in the cavity (its circumcircle does not
    // contain the site), and such a neighbor is exactly a boundary neighbor.
    const std::size_t bad_gen = ++m_bad_counter;
    for (const std::size_t t : bad_triangles) m_bad_stamp[t] = bad_gen;
    const auto is_bad = [&](std::size_t t) {
        return t != kNone && m_bad_stamp[t] == bad_gen;
    };

    // Collect boundary edges (edge of a bad triangle whose neighbor is alive
    // or absent), preserving each bad triangle's CCW orientation.
    struct Boundary {
        std::size_t a;
        std::size_t b;
        std::size_t neighbor;
    };
    std::vector<Boundary> boundary;
    for (const std::size_t t : bad_triangles) {
        Triangle& triangle = m_triangles[t];
        triangle.alive = false;
        for (int edge = 0; edge < 3; ++edge) {
            const std::size_t a = triangle.vertices[(edge + 1U) % 3U];
            const std::size_t b = triangle.vertices[(edge + 2U) % 3U];
            const std::size_t neighbor = triangle.neighbors[edge];
            if (!is_bad(neighbor)) {
                boundary.push_back({a, b, neighbor});
            }
        }
    }

    // Null out stale pointers from living neighbors to the dead cavity.
    for (const Boundary& edge : boundary) {
        if (edge.neighbor == kNone) continue;
        Triangle& neighbor_triangle = m_triangles[edge.neighbor];
        for (int edge_index = 0; edge_index < 3; ++edge_index) {
            const std::size_t neighbor_ref = neighbor_triangle.neighbors[edge_index];
            if (neighbor_ref != kNone && !m_triangles[neighbor_ref].alive) {
                neighbor_triangle.neighbors[edge_index] = kNone;
            }
        }
    }

    // Stitch the boundary into a single CCW loop around the site. Because each
    // boundary edge retains its bad triangle's CCW orientation, the loop is
    // already CCW around the cavity, so consecutive edges satisfy b == next.a.
    std::unordered_map<std::size_t, std::size_t> next_by_start;
    next_by_start.reserve(boundary.size() * 2);
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        next_by_start.emplace(boundary[i].a, i);
    }
    std::vector<Boundary> loop;
    loop.reserve(boundary.size());
    std::vector<char> visited(boundary.size(), 0);
    loop.push_back(boundary[0]);
    visited[0] = 1;
    while (loop.size() < boundary.size()) {
        auto it = next_by_start.find(loop.back().b);
        if (it == next_by_start.end() || visited[it->second]) break;
        visited[it->second] = 1;
        loop.push_back(boundary[it->second]);
    }

    // Verify CCW orientation around the site; reverse if necessary so that
    // every new triangle (a, b, site) is CCW.
    double signed_area = 0.0;
    for (const Boundary& edge : loop) {
        signed_area += static_cast<double>(m_sites[edge.a].x()) *
                           static_cast<double>(m_sites[edge.b].y()) -
                       static_cast<double>(m_sites[edge.a].y()) *
                           static_cast<double>(m_sites[edge.b].x());
    }
    if (signed_area < 0.0) {
        std::reverse(loop.begin(), loop.end());
        // Reverse each edge too so a->b stays consistent with the new order.
        for (Boundary& edge : loop) std::swap(edge.a, edge.b);
    }

    // Build the new triangle fan around the site.
    std::size_t first_new = kNone;
    std::size_t previous_new = kNone;
    for (const Boundary& edge : loop) {
        const std::size_t new_index = addTriangle(edge.a, edge.b, site);
        if (first_new == kNone) first_new = new_index;

        // neighbors[2] is across edge (a, b) -> the external neighbor.
        m_triangles[new_index].neighbors[2] = edge.neighbor;
        if (edge.neighbor != kNone) {
            Triangle& external = m_triangles[edge.neighbor];
            for (int edge_index = 0; edge_index < 3; ++edge_index) {
                const std::size_t na = external.vertices[(edge_index + 1U) % 3U];
                const std::size_t nb = external.vertices[(edge_index + 2U) % 3U];
                if ((na == edge.a && nb == edge.b) ||
                    (na == edge.b && nb == edge.a)) {
                    external.neighbors[edge_index] = new_index;
                    break;
                }
            }
        }

        // Shared edge with the previous new triangle is (b_prev, site) ==
        // (a_curr, site); it is opposite vertex a_prev (neighbors[0]) in the
        // previous triangle and opposite vertex b_curr (neighbors[1]) here.
        if (previous_new != kNone) {
            m_triangles[previous_new].neighbors[0] = new_index;
            m_triangles[new_index].neighbors[1] = previous_new;
        }
        previous_new = new_index;
    }
    if (first_new != kNone && previous_new != kNone && previous_new != first_new) {
        m_triangles[first_new].neighbors[1] = previous_new;
        m_triangles[previous_new].neighbors[0] = first_new;
    }
    m_hint = (first_new != kNone) ? first_new : m_hint;
}

DelaunayTriangulation IncrementalDelaunayBuilder::build() {
    DelaunayTriangulation result;
    result.sites = m_sites;
    if (m_sites.size() < 3U) return result;

    const std::size_t base = m_sites.size();
    const double span = m_scale * 1e6;
    const Point super_a{-span, -span};
    const Point super_b{span, -span};
    const Point super_c{0.0, span * 3.0};
    m_sites.push_back(super_a);
    m_sites.push_back(super_b);
    m_sites.push_back(super_c);

    addTriangle(base, base + 1U, base + 2U);
    m_hint = 0U;

    // Hilbert insertion order for locality-friendly walking.
    std::vector<std::size_t> order(base);
    for (std::size_t i = 0; i < base; ++i) order[i] = i;
    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < base; ++i) {
        min_x = std::min(min_x, static_cast<double>(m_sites[i].x()));
        min_y = std::min(min_y, static_cast<double>(m_sites[i].y()));
        max_x = std::max(max_x, static_cast<double>(m_sites[i].x()));
        max_y = std::max(max_y, static_cast<double>(m_sites[i].y()));
    }
    const double width = std::max(max_x - min_x, 1.0);
    const double height = std::max(max_y - min_y, 1.0);
    std::sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
        const auto key = [&](std::size_t i) {
            const double nx = (static_cast<double>(m_sites[i].x()) - min_x) / width;
            const double ny = (static_cast<double>(m_sites[i].y()) - min_y) / height;
            return hilbertKey(
                static_cast<std::uint32_t>(std::clamp(nx, 0.0, 1.0) * 4095.0),
                static_cast<std::uint32_t>(std::clamp(ny, 0.0, 1.0) * 4095.0));
        };
        return key(left) < key(right);
    });

    for (const std::size_t site : order) {
        insertSite(site);
    }

    result.triangles.clear();
    result.triangles.reserve(m_triangles.size());
    for (const Triangle& triangle : m_triangles) {
        if (!triangle.alive) continue;
        if (triangle.vertices[0] >= base || triangle.vertices[1] >= base ||
            triangle.vertices[2] >= base) {
            continue;
        }
        result.triangles.push_back({triangle.vertices[0], triangle.vertices[1],
                                   triangle.vertices[2]});
    }
    result.sites.resize(base);
    return result;
}

} // namespace skeleton::detail
