#include "common/geometryutil/CommonUtil.h"

#include <gtest/gtest.h>

#include <stdexcept>

namespace common::unittest {

using common::geometryutil::almostEqual;
using common::geometryutil::almostLess;
using common::geometryutil::almostLessOrEqual;
using common::geometryutil::almostGreater;
using common::geometryutil::almostGreaterOrEqual;
using common::geometryutil::isNearZero;
using common::geometryutil::snap;
using common::geometryutil::roundToPrecision;

TEST(CommonUtilTest, ComparesFloatingPointValuesWithRelativeTolerance) {
    EXPECT_TRUE(almostEqual(1.0, 1.0 + 1e-10));
    EXPECT_TRUE(almostEqual(1'000'000.0, 1'000'000.05, 1e-7));
    EXPECT_FALSE(almostEqual(1.0, 1.01));
    EXPECT_THROW(almostEqual(1.0, 1.0, -0.1), std::invalid_argument);
}

TEST(CommonUtilTest, ComparesValuesWithToleranceAwareOrdering) {
    constexpr auto tolerance = 1e-6;

    EXPECT_TRUE(almostLess(1.0, 2.0, tolerance));
    EXPECT_FALSE(almostLess(1.0, 1.0 + 1e-7, tolerance));
    EXPECT_TRUE(almostLessOrEqual(1.0, 1.0 + 1e-7, tolerance));
    EXPECT_FALSE(almostLessOrEqual(2.0, 1.0, tolerance));

    EXPECT_TRUE(almostGreater(2.0, 1.0, tolerance));
    EXPECT_FALSE(almostGreater(1.0 + 1e-7, 1.0, tolerance));
    EXPECT_TRUE(almostGreaterOrEqual(1.0 + 1e-7, 1.0, tolerance));
    EXPECT_FALSE(almostGreaterOrEqual(1.0, 2.0, tolerance));

    EXPECT_THROW(almostLess(1.0, 2.0, -tolerance), std::invalid_argument);
    EXPECT_THROW(almostLessOrEqual(1.0, 2.0, -tolerance),
                 std::invalid_argument);
    EXPECT_THROW(almostGreater(2.0, 1.0, -tolerance), std::invalid_argument);
    EXPECT_THROW(almostGreaterOrEqual(2.0, 1.0, -tolerance),
                 std::invalid_argument);
}

TEST(CommonUtilTest, DetectsValuesNearZero) {
    EXPECT_TRUE(isNearZero(1e-10));
    EXPECT_FALSE(isNearZero(1e-4));
    EXPECT_TRUE(isNearZero(0));
}

TEST(CommonUtilTest, SnapsValuesToNearestStep) {
    EXPECT_DOUBLE_EQ(snap(1.24, 0.5), 1.0);
    EXPECT_DOUBLE_EQ(snap(1.26, 0.5), 1.5);
    EXPECT_DOUBLE_EQ(snap(-1.26, 0.5), -1.5);
    EXPECT_THROW(snap(1.0, 0.0), std::invalid_argument);
}

TEST(CommonUtilTest, RoundsToRequestedDecimalPrecision) {
    EXPECT_DOUBLE_EQ(roundToPrecision(1.23456, 3), 1.235);
    EXPECT_DOUBLE_EQ(roundToPrecision(123.456, -1), 120.0);
    EXPECT_THROW(roundToPrecision(1.0, -1'000), std::out_of_range);
}

} // namespace common::unittest
