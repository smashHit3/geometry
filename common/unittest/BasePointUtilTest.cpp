#include "common/BasePointUtil.h"

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace common::unittest {

TEST(BasePointUtilTest, AddsAndSubtractsPoints) {
    const BasePoint<int> left{3, 4};
    const BasePoint<double> right{1.5, 2.0};

    const auto sum = add(left, right);
    const auto difference = subtract(left, right);

    static_assert(
        std::is_same_v<std::remove_cv_t<decltype(sum)>, BasePoint<double>>);
    EXPECT_DOUBLE_EQ(sum.x(), 4.5);
    EXPECT_DOUBLE_EQ(sum.y(), 6.0);
    EXPECT_DOUBLE_EQ(difference.x(), 1.5);
    EXPECT_DOUBLE_EQ(difference.y(), 2.0);
}

TEST(BasePointUtilTest, MultipliesAndDividesByScalars) {
    const BasePoint<int> point{3, -6};

    const auto product = multiply(point, 2);
    const auto quotient = divide(point, 2);

    EXPECT_EQ(product.x(), 6);
    EXPECT_EQ(product.y(), -12);
    EXPECT_DOUBLE_EQ(quotient.x(), 1.5);
    EXPECT_DOUBLE_EQ(quotient.y(), -3.0);
    EXPECT_THROW(divide(point, 0), std::invalid_argument);
}

TEST(BasePointUtilTest, DotReturnsScalarProjection) {
    const BasePoint<int> left{3, 4};
    const BasePoint<int> right{2, -1};

    EXPECT_EQ(dot(left, right), 2);
}

TEST(BasePointUtilTest, CrossReturnsSignedParallelogramArea) {
    const BasePoint<int> horizontal{3, 0};
    const BasePoint<int> vertical{0, 4};

    EXPECT_EQ(cross(horizontal, vertical), 12);
    EXPECT_EQ(cross(vertical, horizontal), -12);
}

TEST(BasePointUtilTest, SupportsMixedCoordinateTypes) {
    const BasePoint<int> integer_point{2, 3};
    const BasePoint<double> decimal_point{1.5, 2.0};

    EXPECT_DOUBLE_EQ(dot(integer_point, decimal_point), 9.0);
    EXPECT_DOUBLE_EQ(cross(integer_point, decimal_point), -0.5);
}

TEST(BasePointUtilTest, CalculatesLengthsAndDistances) {
    const BasePoint<int> point{3, 4};
    const BasePoint<int> other{6, 8};

    EXPECT_EQ(squared_length(point), 25);
    EXPECT_DOUBLE_EQ(length(point), 5.0);
    EXPECT_EQ(squared_distance(point, other), 25);
    EXPECT_DOUBLE_EQ(distance(point, other), 5.0);
}

TEST(BasePointUtilTest, NormalizesNonZeroPoints) {
    const auto normalized = normalize(BasePoint<int>{3, 4});

    static_assert(std::is_same_v<std::remove_cv_t<decltype(normalized)>,
                                 BasePoint<double>>);
    EXPECT_DOUBLE_EQ(normalized.x(), 0.6);
    EXPECT_DOUBLE_EQ(normalized.y(), 0.8);
    EXPECT_THROW(normalize(BasePoint<int>{}), std::domain_error);
}

TEST(BasePointUtilTest, RotatesPointsCounterclockwise) {
    const auto rotated = rotate(BasePoint<int>{1, 0}, std::acos(-1.0) / 2.0);

    EXPECT_NEAR(rotated.x(), 0.0, 1e-12);
    EXPECT_NEAR(rotated.y(), 1.0, 1e-12);
}

} // namespace common::unittest
