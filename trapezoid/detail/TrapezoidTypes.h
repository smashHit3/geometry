#ifndef TRAPEZOID_DETAIL_TRAPEZOID_TYPES_H
#define TRAPEZOID_DETAIL_TRAPEZOID_TYPES_H

#include "common/BasePoint.h"
#include "common/BasePointUtil.h"
#include "common/CommonUtil.h"

#include <cstddef>

namespace trapezoid::detail {

// A line segment with left endpoint having the smaller x-coordinate.
// If both endpoints share the same x, the one with smaller y is "left".
template <typename T>
struct Segment {
    common::BasePoint<T> left;
    common::BasePoint<T> right;
    std::size_t id = 0;

    // Linear-interpolate the segment's y value at the given x.
    T evalY(T x) const {
        const auto dx = right.x() - left.x();
        if (common::almostEqual(dx, T{0})) {
            return left.y();
        }
        const auto t = (x - left.x()) / dx;
        return left.y() + t * (right.y() - left.y());
    }
};

template <typename T>
struct DagNode;

// A trapezoid in the trapezoidal map.  Its left and right boundaries are
// vertical lines through *leftp and *rightp respectively.  Its top and
// bottom boundaries are the segments *top and *bottom (or the bounding-box
// edges, represented by sentinel segments).
template <typename T>
struct Trapezoid {
    Segment<T>* top = nullptr;
    Segment<T>* bottom = nullptr;
    common::BasePoint<T>* leftp = nullptr;
    common::BasePoint<T>* rightp = nullptr;

    Trapezoid<T>* upper_left = nullptr;
    Trapezoid<T>* lower_left = nullptr;
    Trapezoid<T>* upper_right = nullptr;
    Trapezoid<T>* lower_right = nullptr;

    DagNode<T>* dag_node = nullptr;
    bool removed = false;
};

enum class DagNodeType { X, Y, LEAF };

// A node in the search DAG.
//   X-node  — stores an endpoint *point; query compares q.x with point->x.
//   Y-node  — stores a segment *segment; query compares q against the segment.
//   LEAF    — stores a trapezoid *trapezoid; query terminates here.
template <typename T>
struct DagNode {
    DagNodeType type = DagNodeType::LEAF;
    common::BasePoint<T>* point = nullptr;
    Segment<T>* segment = nullptr;
    Trapezoid<T>* trapezoid = nullptr;
    DagNode<T>* left = nullptr;
    DagNode<T>* right = nullptr;
};

} // namespace trapezoid::detail

#endif // TRAPEZOID_DETAIL_TRAPEZOID_TYPES_H
