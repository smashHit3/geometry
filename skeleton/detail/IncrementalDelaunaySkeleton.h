#ifndef SKELETON_INCREMENTAL_DELAUNAY_SKELETON_H
#define SKELETON_INCREMENTAL_DELAUNAY_SKELETON_H

#include "skeleton/Types.h"

#include <vector>

namespace skeleton {

struct IncrementalDelaunayEdge {
    Point start;
    Point end;
};

struct IncrementalDelaunayOptions {
    // The boundary is sampled no farther apart than this value. A value of zero
    // selects one sixty-fourth of the polygon bounding-box diagonal.
    double max_boundary_segment_length = 0.0;
};

namespace detail {

std::vector<IncrementalDelaunayEdge> medialAxis(
    const std::vector<Point>& polygon,
    const IncrementalDelaunayOptions& options);

} // namespace detail

// Computes an approximation of the medial axis of a simple polygon. The
// boundary is sampled and a Delaunay/Voronoi dual is constructed; decreasing
// max_boundary_segment_length improves the approximation. Input coordinates
// may use any arithmetic BasePoint type; results use double precision.
template <typename coordinate_type>
std::vector<IncrementalDelaunayEdge> incrementalDelaunayMedialAxis(
    const std::vector<common::BasePoint<coordinate_type>>& polygon,
    const IncrementalDelaunayOptions& options = {}) {
    std::vector<Point> double_polygon;
    double_polygon.reserve(polygon.size());
    for (const auto& point : polygon) {
        double_polygon.emplace_back(
            static_cast<double>(point.x()),
            static_cast<double>(point.y()));
    }
    return detail::medialAxis(double_polygon, options);
}

} // namespace skeleton

#endif // SKELETON_INCREMENTAL_DELAUNAY_SKELETON_H
