#include "boolean/detail/RectClipper.h"

#include <gtest/gtest.h>

namespace boolean::unittest {

TEST(RectClipperTest, RetainsTheConfiguredRectangle) {
    const RectClipper::Rectangle rectangle{-10, -5, 20, 15};
    const RectClipper clipper{rectangle};

    EXPECT_EQ(clipper.rectangle().minX(), -10);
    EXPECT_EQ(clipper.rectangle().minY(), -5);
    EXPECT_EQ(clipper.rectangle().maxX(), 20);
    EXPECT_EQ(clipper.rectangle().maxY(), 15);
}

} // namespace boolean::unittest
