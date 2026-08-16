#ifndef BOOLEAN_DETAIL_RECT_CLIPPER_DETAIL_H
#define BOOLEAN_DETAIL_RECT_CLIPPER_DETAIL_H

#include "boolean/RectClipper.h"

#include <cstddef>
#include <vector>

namespace boolean::detail {

// A circular output list permits boundary overlaps to split and rejoin paths
// without copying every already-emitted vertex.
struct OutputVertex {
    IntegerPoint point;
    std::size_t owner = 0;
    OutputVertex* previous = nullptr;
    OutputVertex* next = nullptr;
};

using OutputPath = std::vector<OutputVertex*>;

} // namespace boolean::detail

#endif // BOOLEAN_DETAIL_RECT_CLIPPER_DETAIL_H
