#ifndef COMMON_BASE_POINT_UTIL_H
#define COMMON_BASE_POINT_UTIL_H

#include "common/BasePoint.h"

#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace common {

template <typename... value_types>
using calculation_type = std::conditional_t<
    (std::is_integral_v<value_types> && ...), long double,
    std::common_type_t<value_types...>>;

template <typename left_type, typename right_type>
constexpr auto add(const BasePoint<left_type>& left,
                   const BasePoint<right_type>& right)
    -> BasePoint<calculation_type<left_type, right_type>> {
    using result_type = calculation_type<left_type, right_type>;
    return {static_cast<result_type>(left.x()) + right.x(),
            static_cast<result_type>(left.y()) + right.y()};
}

template <typename left_type, typename right_type>
constexpr auto subtract(const BasePoint<left_type>& left,
                        const BasePoint<right_type>& right)
    -> BasePoint<calculation_type<left_type, right_type>> {
    using result_type = calculation_type<left_type, right_type>;
    return {static_cast<result_type>(left.x()) - right.x(),
            static_cast<result_type>(left.y()) - right.y()};
}

template <typename point_type, typename scalar_type,
          std::enable_if_t<std::is_arithmetic_v<scalar_type>, int> = 0>
constexpr auto multiply(const BasePoint<point_type>& point, scalar_type scalar)
    -> BasePoint<calculation_type<point_type, scalar_type>> {
    using result_type = calculation_type<point_type, scalar_type>;
    return {static_cast<result_type>(point.x()) * scalar,
            static_cast<result_type>(point.y()) * scalar};
}

template <typename point_type, typename scalar_type,
          std::enable_if_t<std::is_arithmetic_v<scalar_type>, int> = 0>
constexpr auto divide(const BasePoint<point_type>& point, scalar_type scalar)
    -> BasePoint<std::common_type_t<point_type, scalar_type, double>> {
    if (scalar == 0) {
        throw std::invalid_argument("Cannot divide a point by zero");
    }

    using result_type = std::common_type_t<point_type, scalar_type, double>;
    return {static_cast<result_type>(point.x()) / scalar,
            static_cast<result_type>(point.y()) / scalar};
}

template <typename left_type, typename right_type>
constexpr auto dot(const BasePoint<left_type>& left,
                   const BasePoint<right_type>& right)
    -> calculation_type<left_type, right_type> {
    using result_type = calculation_type<left_type, right_type>;
    return static_cast<result_type>(left.x()) * right.x() +
           static_cast<result_type>(left.y()) * right.y();
}

template <typename left_type, typename right_type>
constexpr auto cross(const BasePoint<left_type>& left,
                     const BasePoint<right_type>& right)
    -> calculation_type<left_type, right_type> {
    using result_type = calculation_type<left_type, right_type>;
    return static_cast<result_type>(left.x()) * right.y() -
           static_cast<result_type>(left.y()) * right.x();
}

template <typename point_type>
constexpr auto squaredLength(const BasePoint<point_type>& point)
    -> calculation_type<point_type, point_type> {
    return dot(point, point);
}

template <typename point_type>
auto length(const BasePoint<point_type>& point)
    -> decltype(std::sqrt(squaredLength(point))) {
    return std::sqrt(squaredLength(point));
}

template <typename left_type, typename right_type>
constexpr auto squaredDistance(const BasePoint<left_type>& left,
                              const BasePoint<right_type>& right)
    -> calculation_type<left_type, right_type> {
    return squaredLength(subtract(left, right));
}

template <typename left_type, typename right_type>
auto distance(const BasePoint<left_type>& left,
              const BasePoint<right_type>& right)
    -> decltype(std::sqrt(squaredDistance(left, right))) {
    return std::sqrt(squaredDistance(left, right));
}

template <typename point_type>
auto normalize(const BasePoint<point_type>& point)
    -> BasePoint<decltype(point.x() / length(point))> {
    const auto point_length = length(point);
    if (point_length == 0) {
        throw std::domain_error("Cannot normalize a zero-length point");
    }

    using result_type = decltype(point.x() / point_length);
    return {static_cast<result_type>(point.x()) / point_length,
            static_cast<result_type>(point.y()) / point_length};
}

template <typename point_type, typename angle_type,
          std::enable_if_t<std::is_arithmetic_v<angle_type>, int> = 0>
auto rotate(const BasePoint<point_type>& point, angle_type radians)
    -> BasePoint<std::common_type_t<point_type, angle_type, double>> {
    using result_type = std::common_type_t<point_type, angle_type, double>;
    const auto cosine = std::cos(static_cast<result_type>(radians));
    const auto sine = std::sin(static_cast<result_type>(radians));

    return {point.x() * cosine - point.y() * sine,
            point.x() * sine + point.y() * cosine};
}

} // namespace common

#endif // COMMON_BASE_POINT_UTIL_H
