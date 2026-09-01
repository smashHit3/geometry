#ifndef SKELETON_TYPES_H
#define SKELETON_TYPES_H

#include "common/geometry/BasePoint.h"

#include <type_traits>

namespace skeleton {

using Point = common::geometry::BasePoint<double>;

template <typename coordinate_type,
          std::enable_if_t<std::is_arithmetic_v<coordinate_type>, int> = 0>
struct BaseSegment {
    common::geometry::BasePoint<coordinate_type> start;
    common::geometry::BasePoint<coordinate_type> end;
};

using Segment = BaseSegment<double>;

} // namespace skeleton

#endif // SKELETON_TYPES_H
