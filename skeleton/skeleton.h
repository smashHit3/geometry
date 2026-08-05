#ifndef SKELETON_SKELETON_H
#define SKELETON_SKELETON_H

#include "common/BasePoint.h"

#include <vector>

namespace skeleton {

using Point = common::BasePoint<double>;

struct MedialAxisEdge {
    Point start;
    Point end;
};

struct MedialAxisOptions {
    // The boundary is sampled no farther apart than this value. A value of zero
    // selects one sixty-fourth of the polygon bounding-box diagonal.
    double max_boundary_segment_length = 0.0;
};

// Computes an approximation of the medial axis of a simple polygon. The
// boundary is sampled and a Delaunay/Voronoi dual is constructed; decreasing
// max_boundary_segment_length improves the approximation.
std::vector<MedialAxisEdge> medialAxis(
    const std::vector<Point>& polygon,
    const MedialAxisOptions& options = {});

} // namespace skeleton

#endif // SKELETON_SKELETON_H
