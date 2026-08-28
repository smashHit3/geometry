#ifndef SKELETON_TYPES_H
#define SKELETON_TYPES_H

#include "common/BasePoint.h"

#include <type_traits>

namespace skeleton {

using Point = common::BasePoint<double>;

template <typename coordinate_type,
          std::enable_if_t<std::is_arithmetic_v<coordinate_type>, int> = 0>
struct BaseSegment {
    common::BasePoint<coordinate_type> start;
    common::BasePoint<coordinate_type> end;
};

using Segment = BaseSegment<double>;

} // namespace skeleton

#endif // SKELETON_TYPES_H
