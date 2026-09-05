#include "common/geometry/Aabb.h"

#include <gtest/gtest.h>

#include <stdexcept>

namespace common::geometry::unittest {

TEST(AabbTest, TestsIntersectionContainmentAndMerge) {
    const Aabb<double> bounds{1.0, 2.0, 4.0, 6.0};
    const Aabb<double> overlapping{4.0, 5.0, 8.0, 9.0};
    const Aabb<double> disjoint{5.0, 7.0, 8.0, 9.0};

    EXPECT_TRUE(bounds.intersects(overlapping));
    EXPECT_FALSE(bounds.intersects(disjoint));
    EXPECT_TRUE(bounds.contains(Aabb<double>{2.0, 3.0, 3.0, 5.0}));
    EXPECT_EQ(bounds.merged(overlapping).maxX(), 8.0);
    EXPECT_DOUBLE_EQ(bounds.area(), 12.0);
    EXPECT_THROW((Aabb<double>{2.0, 1.0, 1.0, 3.0}),
                 std::invalid_argument);
}

} // namespace common::geometry::unittest
