#include "common/CommonUtil.h"

#include <gtest/gtest.h>

#include <stdexcept>

namespace common::unittest {

TEST(CommonUtilTest, ComparesFloatingPointValuesWithRelativeTolerance) {
    EXPECT_TRUE(almost_equal(1.0, 1.0 + 1e-10));
    EXPECT_TRUE(almost_equal(1'000'000.0, 1'000'000.05, 1e-7));
    EXPECT_FALSE(almost_equal(1.0, 1.01));
    EXPECT_THROW(almost_equal(1.0, 1.0, -0.1), std::invalid_argument);
}

TEST(CommonUtilTest, DetectsValuesNearZero) {
    EXPECT_TRUE(is_near_zero(1e-10));
    EXPECT_FALSE(is_near_zero(1e-4));
    EXPECT_TRUE(is_near_zero(0));
}

TEST(CommonUtilTest, SnapsValuesToNearestStep) {
    EXPECT_DOUBLE_EQ(snap(1.24, 0.5), 1.0);
    EXPECT_DOUBLE_EQ(snap(1.26, 0.5), 1.5);
    EXPECT_DOUBLE_EQ(snap(-1.26, 0.5), -1.5);
    EXPECT_THROW(snap(1.0, 0.0), std::invalid_argument);
}

TEST(CommonUtilTest, RoundsToRequestedDecimalPrecision) {
    EXPECT_DOUBLE_EQ(round_to_precision(1.23456, 3), 1.235);
    EXPECT_DOUBLE_EQ(round_to_precision(123.456, -1), 120.0);
}

} // namespace common::unittest
