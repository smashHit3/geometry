#include "common/BasePoint.h"

#include <gtest/gtest.h>

#include <string>
#include <type_traits>

namespace {

template <typename type, typename = void>
struct is_point_coordinate_type : std::false_type {};

template <typename type>
struct is_point_coordinate_type<
    type, std::void_t<decltype(sizeof(common::BasePoint<type>))>>
    : std::true_type {};

static_assert(is_point_coordinate_type<int>::value);
static_assert(is_point_coordinate_type<double>::value);
static_assert(!is_point_coordinate_type<std::string>::value);

} // namespace

namespace common::unittest {

TEST(BasePointTest, DefaultConstructionInitializesCoordinatesToZero) {
    common::BasePoint<int> origin;
    EXPECT_EQ(origin.x(), 0);
    EXPECT_EQ(origin.y(), 0);
}

TEST(BasePointTest, ValueConstructionSetsCoordinates) {
    common::BasePoint<int> point{1, 2};
    EXPECT_EQ(point.x(), 1);
    EXPECT_EQ(point.y(), 2);
}

TEST(BasePointTest, SupportsFloatingPointCoordinates) {
    common::BasePoint<double> point{1.5, 2.25};
    EXPECT_DOUBLE_EQ(point.x(), 1.5);
    EXPECT_DOUBLE_EQ(point.y(), 2.25);

    point.setX(3.75);
    point.setY(4.5);
    EXPECT_DOUBLE_EQ(point.x(), 3.75);
    EXPECT_DOUBLE_EQ(point.y(), 4.5);
}

TEST(BasePointTest, RejectsNonNumericCoordinateTypes) {
    EXPECT_FALSE((is_point_coordinate_type<std::string>::value));
}

TEST(BasePointTest, CoordinatesCanBeMutated) {
    common::BasePoint<int> point;
    point.setX(3);
    point.y() = 4;
    EXPECT_EQ(point.x(), 3);
    EXPECT_EQ(point.y(), 4);
}

TEST(BasePointTest, ConstAccessorsReadCoordinates) {
    common::BasePoint<int> point{3, 4};
    const common::BasePoint<int>& constant_point = point;
    EXPECT_EQ(constant_point.x(), 3);
    EXPECT_EQ(constant_point.y(), 4);
}

} // namespace common::unittest