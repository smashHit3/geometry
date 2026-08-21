#ifndef SKELETON_STRAIGHT_SKELETON_H
#define SKELETON_STRAIGHT_SKELETON_H

#include "common/BasePoint.h"
#include "skeleton/Types.h"

#include <cstddef>
#include <vector>

namespace skeleton {

struct SkeletonArc {
    std::size_t from;
    std::size_t to;
};

struct StraightSkeleton {
    std::vector<Point> nodes;
    std::vector<SkeletonArc> arcs;
};

} // namespace skeleton

namespace skeleton::detail {

// Builds the straight skeleton of a simple polygon via wavefront propagation.
// A priority queue drives edge events (adjacent wavefront edges collapse) and
// split events (a reflex vertex's bisector pierces an opposite edge) in
// increasing time order. Each LAV is a circular doubly-linked list of moving
// vertices; events are invalidated lazily and recomputed locally, with split
// candidates re-scanned for affected reflex vertices. Uses scale-relative
// floating-point predicates, not exact arithmetic. Polygon vertices must be
// ordered counter-clockwise; clockwise input is reversed automatically.
class StraightSkeletonBuilder {
public:
    explicit StraightSkeletonBuilder(std::vector<Point> polygon);

    StraightSkeleton build();

private:
    struct Node;
    struct Event;

    std::vector<Point> m_polygon;
    double m_scale = 1.0;
    double m_lengthEpsilon = 1e-11;
    double m_areaEpsilon = 1e-18;
};

inline StraightSkeleton straightSkeleton(const std::vector<Point>& polygon) {
    return StraightSkeletonBuilder(polygon).build();
}

} // namespace skeleton::detail

template <typename vertex_type>
skeleton::StraightSkeleton straightSkeleton(
    const std::vector<common::BasePoint<vertex_type>>& polygon) {
    std::vector<skeleton::Point> double_polygon;
    double_polygon.reserve(polygon.size());
    for (const auto& vertex : polygon) {
        double_polygon.emplace_back(static_cast<double>(vertex.x()),
                                    static_cast<double>(vertex.y()));
    }
    return skeleton::detail::straightSkeleton(double_polygon);
}

#endif // SKELETON_STRAIGHT_SKELETON_H
