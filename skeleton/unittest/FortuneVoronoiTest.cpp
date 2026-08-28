#include "skeleton/FortuneVoronoi.h"

#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

namespace skeleton::unittest {

TEST(FortuneVoronoiTest, SplitsBoundsBetweenTwoSites) {
    const std::vector<common::BasePoint<int>> sites{{2, 5}, {8, 5}};

    const auto diagram = voronoiDiagram(sites, common::Aabb<int>{0, 0, 10, 10});
    const auto& cells = diagram.cells;

    ASSERT_EQ(cells.size(), 2U);
    ASSERT_EQ(diagram.edges.size(), 1U);
    EXPECT_EQ(halfPlaneVoronoiCells(sites, common::Aabb<int>{0, 0, 10, 10}).size(),
              cells.size());
    EXPECT_NEAR(diagram.edges.front().start.x(), 5.0, 1e-10);
    EXPECT_NEAR(diagram.edges.front().end.x(), 5.0, 1e-10);
    for (const auto& vertex : cells[0].vertices) EXPECT_LE(vertex.x(), 5.0);
    for (const auto& vertex : cells[1].vertices) EXPECT_GE(vertex.x(), 5.0);
}

TEST(FortuneVoronoiTest, ProducesBoundedCellsAndEdgesForTriangle) {
    const std::vector<common::BasePoint<float>> sites{
        {2.0F, 2.0F}, {8.0F, 2.0F}, {5.0F, 8.0F}};
    const common::Aabb<int> bounds{0, 0, 10, 10};

    const auto diagram = voronoiDiagram(sites, bounds);

    ASSERT_EQ(diagram.edges.size(), 3U);
    ASSERT_EQ(diagram.cells.size(), sites.size());
    for (const auto& edge : diagram.edges) {
        EXPECT_GE(edge.start.x(), 0.0);
        EXPECT_LE(edge.start.x(), 10.0);
        EXPECT_GE(edge.start.y(), 0.0);
        EXPECT_LE(edge.start.y(), 10.0);
        EXPECT_GE(edge.end.x(), 0.0);
        EXPECT_LE(edge.end.x(), 10.0);
        EXPECT_GE(edge.end.y(), 0.0);
        EXPECT_LE(edge.end.y(), 10.0);
    }
    for (const auto& cell : diagram.cells) {
        ASSERT_GE(cell.vertices.size(), 3U);
        for (const auto& vertex : cell.vertices) {
            const double cell_distance = common::squaredDistance(cell.site, vertex);
            for (const auto& site : sites) {
                EXPECT_LE(cell_distance,
                          common::squaredDistance(
                              Point{static_cast<double>(site.x()),
                                    static_cast<double>(site.y())},
                              vertex) +
                              1e-8);
            }
        }
    }
}

TEST(FortuneVoronoiTest, HandlesCollinearSitesWithoutUnboundedOutput) {
    const std::vector<common::BasePoint<int>> sites{{2, 5}, {6, 5}, {10, 5}};

    const auto diagram = voronoiDiagram(sites, common::Aabb<int>{0, 0, 12, 10});

    ASSERT_EQ(diagram.edges.size(), 2U);
    ASSERT_EQ(diagram.cells.size(), 3U);
    for (const auto& cell : diagram.cells) {
        EXPECT_EQ(cell.vertices.size(), 4U);
    }
}

TEST(FortuneVoronoiTest, HandlesCocircularSitesWithoutZeroLengthDiagonals) {
    const std::vector<common::BasePoint<int>> sites{
        {2, 2}, {8, 2}, {8, 8}, {2, 8}};

    const auto diagram = voronoiDiagram(sites, common::Aabb<int>{0, 0, 10, 10});

    ASSERT_EQ(diagram.edges.size(), 4U);
    ASSERT_EQ(diagram.cells.size(), 4U);
    for (const auto& cell : diagram.cells) {
        EXPECT_EQ(cell.vertices.size(), 4U);
    }
}

TEST(FortuneVoronoiTest, ProcessesSitesSharingSweepLevels) {
    const std::vector<common::BasePoint<double>> sites{
        {1.0, 9.0}, {5.0, 9.0}, {9.0, 9.0},
        {2.0, 2.0}, {8.0, 2.0}, {5.0, 5.0}};

    const auto diagram =
        voronoiDiagram(sites, common::Aabb<double>{0.0, 0.0, 10.0, 10.0});

    ASSERT_EQ(diagram.cells.size(), sites.size());
    EXPECT_FALSE(diagram.edges.empty());
    for (const auto& edge : diagram.edges) {
        for (const Point& endpoint : {edge.start, edge.end}) {
            const double edge_distance =
                common::squaredDistance(edge.firstSite, endpoint);
            EXPECT_NEAR(edge_distance,
                        common::squaredDistance(edge.secondSite, endpoint), 1e-8);
            for (const auto& site : sites) {
                EXPECT_LE(edge_distance,
                          common::squaredDistance(Point{site.x(), site.y()}, endpoint) +
                              1e-8);
            }
        }
    }
}

TEST(FortuneVoronoiTest, BuildsManyNondegenerateSites) {
    std::vector<common::BasePoint<double>> sites;
    sites.reserve(128);
    for (int index = 0; index < 128; ++index) {
        sites.emplace_back((index * 47) % 997 + 0.001 * index,
                           (index * index * 73 + index * 19) % 991 +
                               0.003 * index);
    }

    const auto diagram =
        voronoiDiagram(sites, common::Aabb<double>{0.0, 0.0, 1000.0, 1000.0});

    ASSERT_EQ(diagram.cells.size(), sites.size());
    EXPECT_FALSE(diagram.edges.empty());
    for (const auto& cell : diagram.cells) {
        EXPECT_GE(cell.vertices.size(), 3U);
        for (const auto& vertex : cell.vertices) {
            const double cell_distance = common::squaredDistance(cell.site, vertex);
            for (const auto& site : sites) {
                EXPECT_LE(cell_distance,
                          common::squaredDistance(
                              Point{site.x(), site.y()}, vertex) +
                              1e-4);
            }
        }
    }
    for (const auto& edge : diagram.edges) {
        EXPECT_NEAR(common::squaredDistance(edge.firstSite, edge.start),
                    common::squaredDistance(edge.secondSite, edge.start), 1e-5);
        EXPECT_NEAR(common::squaredDistance(edge.firstSite, edge.end),
                    common::squaredDistance(edge.secondSite, edge.end), 1e-5);
        for (const auto& site : sites) {
            const Point double_site{site.x(), site.y()};
            EXPECT_LE(common::squaredDistance(edge.firstSite, edge.start),
                      common::squaredDistance(double_site, edge.start) + 1e-4);
            EXPECT_LE(common::squaredDistance(edge.firstSite, edge.end),
                      common::squaredDistance(double_site, edge.end) + 1e-4);
        }
    }
}

TEST(FortuneVoronoiTest, RejectsDuplicateSites) {
    const std::vector<common::BasePoint<double>> sites{{1.0, 1.0}, {1.0, 1.0}};

    EXPECT_THROW(voronoiDiagram(sites, common::Aabb<double>{0.0, 0.0, 2.0, 2.0}),
                 std::invalid_argument);
}

} // namespace skeleton::unittest
