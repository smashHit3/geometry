#include "skeleton/detail/IncrementalDelaunaySkeleton.h"

#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

namespace skeleton::unittest {

TEST(IncrementalDelaunaySkeletonTest, ProducesInteriorEdgesForASquare) {
    const std::vector<common::BasePoint<int>> square{
        {0, 0}, {4, 0}, {4, 4}, {0, 4}};

    const auto edges = incrementalDelaunayMedialAxis(
        square, {.max_boundary_segment_length = 0.25});

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

TEST(IncrementalDelaunaySkeletonTest, SupportsConcavePolygons) {
    const std::vector<common::BasePoint<float>> polygon{
        {0.0F, 0.0F}, {5.0F, 0.0F}, {5.0F, 2.0F},
        {2.0F, 2.0F}, {2.0F, 5.0F}, {0.0F, 5.0F}};

    const auto edges = incrementalDelaunayMedialAxis(
        polygon, {.max_boundary_segment_length = 0.25});

    EXPECT_FALSE(edges.empty());
}

TEST(IncrementalDelaunaySkeletonTest, RejectsNonSimplePolygons) {
    const std::vector<common::BasePoint<double>> bow_tie{
        {0.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}, {2.0, 0.0}};

    EXPECT_THROW(incrementalDelaunayMedialAxis(bow_tie), std::invalid_argument);
}

} // namespace skeleton::unittest
