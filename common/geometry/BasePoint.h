#ifndef COMMON_GEOMETRY_BASEPOINT_H
#define COMMON_GEOMETRY_BASEPOINT_H

#include <type_traits>

namespace common::geometry {

template <typename point_type,
          std::enable_if_t<std::is_arithmetic_v<point_type>, int> = 0>
class BasePoint {
public:
    constexpr BasePoint() = default;
    constexpr BasePoint(point_type x, point_type y)
        : m_x(x), m_y(y) {}
    ~BasePoint() = default;

    constexpr point_type& x() {
        return m_x;
    }

    constexpr const point_type& x() const {
        return m_x;
    }

    constexpr point_type& y() {
        return m_y;
    }

    constexpr const point_type& y() const {
        return m_y;
    }

    constexpr void setX(point_type x) {
        m_x = x;
    }

    constexpr void setY(point_type y) {
        m_y = y;
    }

private:
    point_type m_x = 0;
    point_type m_y = 0;
};

} // namespace common::geometry

#endif // COMMON_GEOMETRY_BASEPOINT_H
