#ifndef COMMON_COMMON_UTIL_H
#define COMMON_COMMON_UTIL_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace common {

inline constexpr double defaultTolerance = 1e-9;

template <typename left_type, typename right_type,
          typename value_type = std::common_type_t<left_type, right_type, double>,
          std::enable_if_t<std::is_arithmetic_v<left_type> &&
                               std::is_arithmetic_v<right_type>,
                           int> = 0>
constexpr bool almostEqual(
    left_type left, right_type right,
    value_type tolerance = static_cast<value_type>(defaultTolerance)) {
    if (tolerance < 0) {
        throw std::invalid_argument("Tolerance cannot be negative");
    }
    if (left == right) {
        return true;
    }

    const auto difference = std::fabs(static_cast<value_type>(left) - right);
    const auto scale = std::max(
        value_type{1}, std::max(std::fabs(static_cast<value_type>(left)),
                                std::fabs(static_cast<value_type>(right))));
    return difference <= tolerance * scale;
}

template <typename left_type, typename right_type,
          typename value_type = std::common_type_t<left_type, right_type, double>,
          std::enable_if_t<std::is_arithmetic_v<left_type> &&
                               std::is_arithmetic_v<right_type>,
                           int> = 0>
constexpr bool almostLess(
    left_type left, right_type right,
    value_type tolerance = static_cast<value_type>(defaultTolerance)) {
    return !almostEqual(left, right, tolerance) &&
           static_cast<value_type>(left) < static_cast<value_type>(right);
}

template <typename left_type, typename right_type,
          typename value_type = std::common_type_t<left_type, right_type, double>,
          std::enable_if_t<std::is_arithmetic_v<left_type> &&
                               std::is_arithmetic_v<right_type>,
                           int> = 0>
constexpr bool almostLessOrEqual(
    left_type left, right_type right,
    value_type tolerance = static_cast<value_type>(defaultTolerance)) {
    return almostEqual(left, right, tolerance) ||
           static_cast<value_type>(left) < static_cast<value_type>(right);
}

template <typename left_type, typename right_type,
          typename value_type = std::common_type_t<left_type, right_type, double>,
          std::enable_if_t<std::is_arithmetic_v<left_type> &&
                               std::is_arithmetic_v<right_type>,
                           int> = 0>
constexpr bool almostGreater(
    left_type left, right_type right,
    value_type tolerance = static_cast<value_type>(defaultTolerance)) {
    return !almostEqual(left, right, tolerance) &&
           static_cast<value_type>(left) > static_cast<value_type>(right);
}

template <typename left_type, typename right_type,
          typename value_type = std::common_type_t<left_type, right_type, double>,
          std::enable_if_t<std::is_arithmetic_v<left_type> &&
                               std::is_arithmetic_v<right_type>,
                           int> = 0>
constexpr bool almostGreaterOrEqual(
    left_type left, right_type right,
    value_type tolerance = static_cast<value_type>(defaultTolerance)) {
    return almostEqual(left, right, tolerance) ||
           static_cast<value_type>(left) > static_cast<value_type>(right);
}

template <typename value_type,
          std::enable_if_t<std::is_arithmetic_v<value_type>, int> = 0>
constexpr bool isNearZero(
    value_type value,
    std::common_type_t<value_type, double> tolerance =
        static_cast<std::common_type_t<value_type, double>>(
            defaultTolerance)) {
    return almostEqual(value, 0, tolerance);
}

template <typename value_type, typename step_type,
          std::enable_if_t<std::is_arithmetic_v<value_type> &&
                               std::is_arithmetic_v<step_type>,
                           int> = 0>
auto snap(value_type value, step_type step)
    -> std::common_type_t<value_type, step_type, double> {
    if (step <= 0) {
        throw std::invalid_argument("Snap step must be positive");
    }

    using result_type = std::common_type_t<value_type, step_type, double>;
    return std::round(static_cast<result_type>(value) / step) * step;
}

template <typename value_type,
          std::enable_if_t<std::is_arithmetic_v<value_type>, int> = 0>
auto roundToPrecision(value_type value, int decimal_places)
    -> std::common_type_t<value_type, double> {
    using result_type = std::common_type_t<value_type, double>;
    const auto scale = std::pow(result_type{10}, decimal_places);
    if (!std::isfinite(scale) || scale == result_type{0}) {
        throw std::out_of_range("Decimal precision is out of range");
    }

    return std::round(static_cast<result_type>(value) * scale) / scale;
}

} // namespace common

#endif // COMMON_COMMON_UTIL_H
