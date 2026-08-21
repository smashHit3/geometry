#ifndef SKELETON_DETAIL_INCREMENTAL_DELAUNAY_SKELETON_H
#define SKELETON_DETAIL_INCREMENTAL_DELAUNAY_SKELETON_H

#include "common/Aabb.h"
#include "common/BasePoint.h"
#include "common/BasePointUtil.h"
#include "skeleton/Types.h"

#include <cstddef>
#include <vector>

namespace skeleton {

struct DelaunayTriangle {
    std::size_t first;
    std::size_t second;
    std::size_t third;
};

struct DelaunayTriangulation {
    std::vector<Point> sites;
    std::vector<DelaunayTriangle> triangles;
};

} // namespace skeleton

namespace skeleton::detail {

// Builds a Delaunay triangulation of a planar point set via incremental
// Bowyer-Watson insertion. Insertion order is driven by a Hilbert curve to
// keep the walking locality search O(log n) expected per point; cavity
// retriangulation rebuilds only the triangles whose circumcircles contain the
// new site. The result is the Delaunay triangulation, whose dual approximates
// the point-site skeleton. Uses scale-relative floating-point predicates, not
// exact arithmetic.
class IncrementalDelaunayBuilder {
public:
    explicit IncrementalDelaunayBuilder(std::vector<Point> sites);

    DelaunayTriangulation build();

private:
    struct Triangle {
        std::size_t vertices[3];
        std::size_t neighbors[3]; // neighbor across the edge opposite vertices[i]
        bool alive = true;
    };

    std::vector<Point> m_sites;
    std::vector<Triangle> m_triangles;
    std::vector<std::size_t> m_visit_stamp; // timestamp cache for cavity BFS
    std::vector<std::size_t> m_bad_stamp;   // timestamp cache for cavity membership
    std::size_t m_visit_counter = 0;
    std::size_t m_bad_counter = 0;
    double m_scale = 1.0;
    double m_areaEpsilon = 1e-18;
    std::size_t m_hint = 0;

    static constexpr std::size_t kNone = static_cast<std::size_t>(-1);

    double orientation(std::size_t a, std::size_t b, std::size_t c) const;
    bool inCircumcircle(std::size_t a, std::size_t b, std::size_t c,
                        std::size_t p) const;
    std::size_t locate(std::size_t site) const;
    std::size_t addTriangle(std::size_t a, std::size_t b, std::size_t c);
    void insertSite(std::size_t site);
};

inline DelaunayTriangulation incrementalDelaunayTriangulation(
    const std::vector<Point>& sites) {
    return IncrementalDelaunayBuilder(sites).build();
}

} // namespace skeleton::detail

template <typename site_type>
skeleton::DelaunayTriangulation delaunayTriangulation(
    const std::vector<common::BasePoint<site_type>>& sites) {
    std::vector<skeleton::Point> double_sites;
    double_sites.reserve(sites.size());
    for (const auto& site : sites) {
        double_sites.emplace_back(static_cast<double>(site.x()),
                                  static_cast<double>(site.y()));
    }
    return skeleton::detail::incrementalDelaunayTriangulation(double_sites);
}

#endif // SKELETON_DETAIL_INCREMENTAL_DELAUNAY_SKELETON_H
