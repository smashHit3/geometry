#ifndef COMMON_AABB_H
#define COMMON_AABB_H

#include <algorithm>
#include <stdexcept>
#include <type_traits>

namespace common {

// Axis-aligned bounding box represented by its minimum and maximum coordinates.
template <typename coordinate_type,
          std::enable_if_t<std::is_arithmetic_v<coordinate_type>, int> = 0>
class Aabb {
public:
    constexpr Aabb(coordinate_type min_x, coordinate_type min_y,
                   coordinate_type max_x, coordinate_type max_y)
        : m_min_x(min_x), m_min_y(min_y), m_max_x(max_x), m_max_y(max_y) {
        if (min_x > max_x || min_y > max_y) {
            throw std::invalid_argument("AABB minimum must not exceed maximum");
        }
    }

    constexpr coordinate_type minX() const { return m_min_x; }
    constexpr coordinate_type minY() const { return m_min_y; }
    constexpr coordinate_type maxX() const { return m_max_x; }
    constexpr coordinate_type maxY() const { return m_max_y; }

    constexpr bool intersects(const Aabb& other) const {
        return m_min_x <= other.m_max_x && m_max_x >= other.m_min_x &&
               m_min_y <= other.m_max_y && m_max_y >= other.m_min_y;
    }

    constexpr bool contains(const Aabb& other) const {
        return m_min_x <= other.m_min_x && m_min_y <= other.m_min_y &&
               m_max_x >= other.m_max_x && m_max_y >= other.m_max_y;
    }

    constexpr Aabb merged(const Aabb& other) const {
        return {std::min(m_min_x, other.m_min_x), std::min(m_min_y, other.m_min_y),
                std::max(m_max_x, other.m_max_x), std::max(m_max_y, other.m_max_y)};
    }

    constexpr auto area() const
        -> std::common_type_t<coordinate_type, double> {
        using result_type = std::common_type_t<coordinate_type, double>;
        return (static_cast<result_type>(m_max_x) - m_min_x) *
               (static_cast<result_type>(m_max_y) - m_min_y);
    }

private:
    coordinate_type m_min_x;
    coordinate_type m_min_y;
    coordinate_type m_max_x;
    coordinate_type m_max_y;
};

} // namespace common

#endif // COMMON_AABB_H
