#include "skeleton/skeleton.h"

#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

namespace skeleton::unittest {

TEST(MedialAxisTest, ProducesInteriorEdgesForASquare) {
    const std::vector<common::BasePoint<double>> square{
        {0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0}};

    const auto edges = medialAxis(square, {.max_boundary_segment_length = 0.25});

    ASSERT_FALSE(edges.empty());
    for (const auto& edge : edges) {
        EXPECT_GE(edge.start.x(), 0.0);
        EXPECT_LE(edge.start.x(), 4.0);
        EXPECT_GE(edge.start.y(), 0.0);
        EXPECT_LE(edge.start.y(), 4.0);
        EXPECT_GE(edge.end.x(), 0.0);
        EXPECT_LE(edge.end.x(), 4.0);
        EXPECT_GE(edge.end.y(), 0.0);
        EXPECT_LE(edge.end.y(), 4.0);
    }
}

TEST(MedialAxisTest, SupportsConcavePolygons) {
    const std::vector<Point> polygon{
        {0.0, 0.0}, {5.0, 0.0}, {5.0, 2.0}, {2.0, 2.0}, {2.0, 5.0}, {0.0, 5.0}};

    const auto edges = medialAxis(polygon, {.max_boundary_segment_length = 0.25});

    EXPECT_FALSE(edges.empty());
}

TEST(MedialAxisTest, RejectsNonSimplePolygons) {
    const std::vector<common::BasePoint<double>> bow_tie{
        {0.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}, {2.0, 0.0}};

    EXPECT_THROW(medialAxis(bow_tie), std::invalid_argument);
}

} // namespace skeleton::unittest
