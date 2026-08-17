#ifndef BOOLEAN_DETAIL_RECT_CLIPPER_H
#define BOOLEAN_DETAIL_RECT_CLIPPER_H

#include "common/Aabb.h"
#include "common/BasePoint.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace boolean {

using IntegerPoint = common::BasePoint<std::int64_t>;
using Path = std::vector<IntegerPoint>;
using Paths = std::vector<Path>;

// Clips closed polygon paths against an axis-aligned rectangle.
//
// A single input path may yield multiple output paths, so callers should use
// Paths for both input and output. Paths are not implicitly closed; the last
// vertex is connected to the first vertex during clipping.
class RectClipper {
public:
    using Rectangle = common::Aabb<std::int64_t>;

    explicit RectClipper(Rectangle rectangle)
        : m_rectangle(rectangle) {}

    const Rectangle& rectangle() const {
        return m_rectangle;
    }

    Paths clip(const Paths& subjects) const;

private:
    enum class Location {
        left,
        top,
        right,
        bottom,
        inside,
    };

    struct Intersection {
        IntegerPoint point;
        Location boundary;
    };

    Location locate(const IntegerPoint& point) const;
    bool intersectsRectangle(const Path& subject) const;
    std::vector<Intersection> findIntersections(
        const IntegerPoint& from, const IntegerPoint& to) const;

    Rectangle m_rectangle;
};

} // namespace boolean

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

#endif // BOOLEAN_DETAIL_RECT_CLIPPER_H
