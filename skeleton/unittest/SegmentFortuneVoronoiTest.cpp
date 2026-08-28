#include "skeleton/SegmentFortuneVoronoi.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace skeleton::unittest {
namespace {

double squaredDistanceToSegment(const Point& point, const Segment& segment) {
    const Point direction = segment.end - segment.start;
    const double lengthSquared = common::squaredLength(direction);
    const double parameter =
        std::clamp(common::dot(point - segment.start, direction) / lengthSquared,
                   0.0, 1.0);
    return common::squaredDistance(
        point, segment.start + direction * parameter);
}

} // namespace

TEST(SegmentFortuneVoronoiTest, SplitsBoundsBetweenParallelSegments) {
    const std::vector<BaseSegment<int>> sites{
        {{2, 2}, {8, 2}},
        {{2, 8}, {8, 8}},
    };

    const auto diagram = voronoiDiagram(
        sites, common::Aabb<int>{0, 0, 10, 10});

    ASSERT_EQ(diagram.cells.size(), sites.size());
    ASSERT_FALSE(diagram.edges.empty());
    double minimumX = 10.0;
    double maximumX = 0.0;
    for (const SegmentVoronoiEdge& edge : diagram.edges) {
        EXPECT_FALSE(edge.curved);
        for (const Point& vertex : edge.vertices) {
            EXPECT_NEAR(vertex.y(), 5.0, 1e-6);
            minimumX = std::min(minimumX, vertex.x());
            maximumX = std::max(maximumX, vertex.x());
        }
    }
    EXPECT_NEAR(minimumX, 0.0, 1e-6);
    EXPECT_NEAR(maximumX, 10.0, 1e-6);
}

TEST(SegmentFortuneVoronoiTest, DiscretizesParabolicBisectors) {
    const std::vector<BaseSegment<double>> sites{
        {{1.0, 2.0}, {9.0, 2.0}},
        {{4.0, 8.0}, {6.0, 8.0}},
    };

    const auto diagram = segmentVoronoiDiagram(
        sites, common::Aabb<double>{0.0, 0.0, 10.0, 10.0}, 1e-4);

    ASSERT_EQ(diagram.cells.size(), sites.size());
    ASSERT_FALSE(diagram.edges.empty());
    EXPECT_TRUE(std::any_of(
        diagram.edges.begin(), diagram.edges.end(),
        [](const SegmentVoronoiEdge& edge) {
            return edge.curved && edge.vertices.size() > 2;
        }));

    for (const SegmentVoronoiEdge& edge : diagram.edges) {
        for (const Point& vertex : edge.vertices) {
            EXPECT_GE(vertex.x(), -1e-8);
            EXPECT_LE(vertex.x(), 10.0 + 1e-8);
            EXPECT_GE(vertex.y(), -1e-8);
            EXPECT_LE(vertex.y(), 10.0 + 1e-8);
            EXPECT_NEAR(squaredDistanceToSegment(vertex,
                                                diagram.cells[edge.firstSite].site),
                        squaredDistanceToSegment(
                            vertex, diagram.cells[edge.secondSite].site),
                        2e-3);
        }
    }
}

TEST(SegmentFortuneVoronoiTest, RecordsIncidentEdgesForEachCell) {
    const std::vector<BaseSegment<int>> sites{
        {{1, 1}, {4, 1}},
        {{6, 3}, {9, 3}},
        {{2, 8}, {8, 8}},
    };

    const auto diagram =
        segmentVoronoiDiagram(sites, common::Aabb<int>{0, 0, 10, 10});

    ASSERT_EQ(diagram.cells.size(), sites.size());
    ASSERT_FALSE(diagram.edges.empty());
    for (std::size_t index = 0; index < diagram.edges.size(); ++index) {
        const SegmentVoronoiEdge& edge = diagram.edges[index];
        EXPECT_NE(std::find(diagram.cells[edge.firstSite].edges.begin(),
                            diagram.cells[edge.firstSite].edges.end(), index),
                  diagram.cells[edge.firstSite].edges.end());
        EXPECT_NE(std::find(diagram.cells[edge.secondSite].edges.begin(),
                            diagram.cells[edge.secondSite].edges.end(), index),
                  diagram.cells[edge.secondSite].edges.end());
        for (const Point& vertex : edge.vertices) {
            const double edgeDistance = squaredDistanceToSegment(
                vertex, diagram.cells[edge.firstSite].site);
            EXPECT_NEAR(edgeDistance,
                        squaredDistanceToSegment(
                            vertex, diagram.cells[edge.secondSite].site),
                        2e-3);
            for (const SegmentVoronoiCell& cell : diagram.cells) {
                EXPECT_LE(edgeDistance,
                          squaredDistanceToSegment(vertex, cell.site) + 2e-3);
            }
        }
    }
}

TEST(SegmentFortuneVoronoiTest, RejectsDegenerateAndIntersectingSites) {
    EXPECT_THROW(
        segmentVoronoiDiagram(
            std::vector<BaseSegment<int>>{{{1, 1}, {1, 1}}},
            common::Aabb<int>{0, 0, 10, 10}),
        std::invalid_argument);

    EXPECT_THROW(
        segmentVoronoiDiagram(
            std::vector<BaseSegment<int>>{
                {{1, 1}, {9, 9}},
                {{1, 9}, {9, 1}},
            },
            common::Aabb<int>{0, 0, 10, 10}),
        std::invalid_argument);

    EXPECT_THROW(
        segmentVoronoiDiagram(
            std::vector<BaseSegment<int>>{
                {{1, 5}, {8, 5}},
                {{4, 5}, {8, 5}},
            },
            common::Aabb<int>{0, 0, 10, 10}),
        std::invalid_argument);

    EXPECT_THROW(
        segmentVoronoiDiagram(
            std::vector<BaseSegment<int>>{{{1, 1}, {2, 2}}},
            common::Aabb<int>{0, 0, 10, 10}, -1.0),
        std::invalid_argument);

}

TEST(SegmentFortuneVoronoiTest,
     HandlesEndpointOnAnotherSupportingLine) {
    const std::vector<BaseSegment<int>> sites{
        {{0, 0}, {2, 0}},
        {{3, 0}, {4, 1}},
    };

    const auto diagram =
        segmentVoronoiDiagram(sites, common::Aabb<int>{-1, -2, 6, 3});

    ASSERT_EQ(diagram.cells.size(), sites.size());
    ASSERT_FALSE(diagram.edges.empty());
    for (const SegmentVoronoiEdge& edge : diagram.edges) {
        for (const Point& vertex : edge.vertices) {
            const double edgeDistance = squaredDistanceToSegment(
                vertex, diagram.cells[edge.firstSite].site);
            EXPECT_NEAR(edgeDistance,
                        squaredDistanceToSegment(
                            vertex, diagram.cells[edge.secondSite].site),
                        2e-3);
        }
    }
}

TEST(SegmentFortuneVoronoiTest, HandlesEmptyAndSingleSiteInputs) {
    const auto empty = segmentVoronoiDiagram(
        std::vector<BaseSegment<int>>{}, common::Aabb<int>{0, 0, 10, 10});
    EXPECT_TRUE(empty.edges.empty());
    EXPECT_TRUE(empty.cells.empty());

    const std::vector<BaseSegment<int>> sites{{{1, 2}, {8, 7}}};
    const auto single =
        segmentVoronoiDiagram(sites, common::Aabb<int>{0, 0, 10, 10});
    EXPECT_TRUE(single.edges.empty());
    ASSERT_EQ(single.cells.size(), 1U);
    EXPECT_TRUE(single.cells.front().edges.empty());
}

TEST(SegmentFortuneVoronoiTest, RetainsNarrowVisibleBisectors) {
    const std::vector<BaseSegment<int>> sites{
        {{22, 90}, {55, 64}},
        {{51, 7}, {71, 86}},
        {{31, 0}, {7, 22}},
    };

    const auto diagram =
        segmentVoronoiDiagram(sites, common::Aabb<int>{0, 0, 100, 100});

    EXPECT_TRUE(std::any_of(
        diagram.edges.begin(), diagram.edges.end(),
        [](const SegmentVoronoiEdge& edge) {
            return edge.firstSite == 0 && edge.secondSite == 2;
        }));
}

TEST(SegmentFortuneVoronoiTest,
     MaintainsBeachOrderForInterleavedArcAndCircleEvents) {
    const std::vector<BaseSegment<double>> sites{
        {{1.57, 1.42}, {1.2044814371685406, 3.1416260279821584}},
        {{6.97, 0.22}, {8.2725551207105905, -1.0511216139735087}},
        {{11.44, 1.05}, {11.197165397784417, 0.73214570009704127}},
        {{0.27, 6.22}, {-0.8336507104173807, 5.1316732524672517}},
        {{6.86, 6.16}, {8.0294776402453074, 6.8141575108230965}},
        {{10.77, 5.2}, {8.781396000563575, 4.409396348714819}},
        {{1.47, 11.23}, {2.539169145425544, 12.577322284559273}},
        {{5.48, 10.08}, {6.988885162159968, 11.251053187269259}},
        {{10.81, 10.03}, {9.797260447654669, 11.163471922508615}},
    };

    const auto diagram = segmentVoronoiDiagram(
        sites, common::Aabb<double>{-4.0, -4.0, 17.0, 17.0}, 0.05);

    ASSERT_EQ(diagram.cells.size(), sites.size());
    ASSERT_FALSE(diagram.edges.empty());
    for (const SegmentVoronoiEdge& edge : diagram.edges) {
        ASSERT_LT(edge.firstSite, diagram.cells.size());
        ASSERT_LT(edge.secondSite, diagram.cells.size());
        for (const Point& vertex : edge.vertices) {
            const double edgeDistance = squaredDistanceToSegment(
                vertex, diagram.cells[edge.firstSite].site);
            EXPECT_NEAR(edgeDistance,
                        squaredDistanceToSegment(
                            vertex, diagram.cells[edge.secondSite].site),
                        2e-3);
            for (const SegmentVoronoiCell& cell : diagram.cells) {
                EXPECT_LE(edgeDistance,
                          squaredDistanceToSegment(vertex, cell.site) + 2e-3);
            }
        }
    }
}

TEST(SegmentFortuneVoronoiTest, HandlesOverlappingSiteAabbs) {
    const std::vector<BaseSegment<double>> sites{
        {{0.0, 0.0}, {10.0, 10.0}},
        {{0.0, 8.0}, {2.0, 10.0}},
        {{8.0, 0.0}, {10.0, 2.0}},
    };
    const common::Aabb<double> bounds{-2.0, -2.0, 12.0, 12.0};

    const auto diagram = segmentVoronoiDiagram(sites, bounds, 1e-4);

    ASSERT_EQ(diagram.cells.size(), sites.size());
    ASSERT_FALSE(diagram.edges.empty());
    for (const SegmentVoronoiEdge& edge : diagram.edges) {
        for (const Point& vertex : edge.vertices) {
            EXPECT_GE(vertex.x(), bounds.minX() - 1e-8);
            EXPECT_LE(vertex.x(), bounds.maxX() + 1e-8);
            EXPECT_GE(vertex.y(), bounds.minY() - 1e-8);
            EXPECT_LE(vertex.y(), bounds.maxY() + 1e-8);
            const double edgeDistance = squaredDistanceToSegment(
                vertex, diagram.cells[edge.firstSite].site);
            EXPECT_NEAR(edgeDistance,
                        squaredDistanceToSegment(
                            vertex, diagram.cells[edge.secondSite].site),
                        2e-3);
            for (const SegmentVoronoiCell& cell : diagram.cells) {
                EXPECT_LE(edgeDistance,
                          squaredDistanceToSegment(vertex, cell.site) + 2e-3);
            }
        }
    }
}

TEST(SegmentFortuneVoronoiTest, RetainsNearRepeatedPssCircleEvent) {
    const std::vector<BaseSegment<double>> sites{
        {{77.9652384404999, 73.5787644057206},
         {81.6543654453846, 68.9567735320739}},
        {{18.3989254667555, 30.5887526371546},
         {16.9279740315596, 33.3813585886033}},
        {{61.0893604840448, 6.9245679040964},
         {56.9530695747596, 4.7842058334866}},
    };

    const auto diagram = segmentVoronoiDiagram(
        sites, common::Aabb<double>{0.0, 0.0, 100.0, 100.0}, 1e-4);

    ASSERT_EQ(diagram.cells.size(), sites.size());
    for (const SegmentVoronoiEdge& edge : diagram.edges) {
        for (const Point& vertex : edge.vertices) {
            const double edgeDistance = squaredDistanceToSegment(
                vertex, diagram.cells[edge.firstSite].site);
            EXPECT_NEAR(edgeDistance,
                        squaredDistanceToSegment(
                            vertex, diagram.cells[edge.secondSite].site),
                        2e-3);
            for (const SegmentVoronoiCell& cell : diagram.cells) {
                EXPECT_LE(edgeDistance,
                          squaredDistanceToSegment(vertex, cell.site) + 2e-3);
            }
        }
    }
}

TEST(SegmentFortuneVoronoiTest, HandlesLargeCoordinateScale) {
    constexpr double scale = 1e7;
    const std::vector<BaseSegment<double>> sites{
        {{2.0 * scale, 2.0 * scale}, {8.0 * scale, 2.0 * scale}},
        {{2.0 * scale, 8.0 * scale}, {8.0 * scale, 8.0 * scale}},
    };

    const auto diagram = segmentVoronoiDiagram(
        sites, common::Aabb<double>{0.0, 0.0, 10.0 * scale, 10.0 * scale},
        100.0);

    ASSERT_FALSE(diagram.edges.empty());
    double minimumX = 10.0 * scale;
    double maximumX = 0.0;
    for (const SegmentVoronoiEdge& edge : diagram.edges) {
        for (const Point& vertex : edge.vertices) {
            EXPECT_NEAR(vertex.y(), 5.0 * scale, 1e-3);
            minimumX = std::min(minimumX, vertex.x());
            maximumX = std::max(maximumX, vertex.x());
        }
    }
    EXPECT_NEAR(minimumX, 0.0, 1e-3);
    EXPECT_NEAR(maximumX, 10.0 * scale, 1e-3);
}

TEST(SegmentFortuneVoronoiTest, NormalizesTinyDiagramsByTheirActualRadius) {
    constexpr double scale = 1e-200;
    const std::vector<BaseSegment<double>> sites{
        {{2.0 * scale, 2.0 * scale}, {8.0 * scale, 2.0 * scale}},
        {{2.0 * scale, 8.0 * scale}, {8.0 * scale, 8.0 * scale}},
    };

    const auto diagram = segmentVoronoiDiagram(
        sites,
        common::Aabb<double>{0.0, 0.0, 10.0 * scale, 10.0 * scale},
        1e-205);

    ASSERT_FALSE(diagram.edges.empty());
    double minimumX = 10.0 * scale;
    double maximumX = 0.0;
    for (const SegmentVoronoiEdge& edge : diagram.edges) {
        for (const Point& vertex : edge.vertices) {
            EXPECT_NEAR(vertex.y(), 5.0 * scale, 1e-210);
            minimumX = std::min(minimumX, vertex.x());
            maximumX = std::max(maximumX, vertex.x());
        }
    }
    EXPECT_NEAR(minimumX, 0.0, 1e-210);
    EXPECT_NEAR(maximumX, 10.0 * scale, 1e-210);
}

TEST(SegmentFortuneVoronoiTest, NormalizesFiniteDoubleExtremesSafely) {
    const double extent = std::numeric_limits<double>::max();
    const double half = extent * 0.5;
    const double quarter = extent * 0.25;
    const std::vector<BaseSegment<double>> sites{
        {{-half, -quarter}, {half, -quarter}},
        {{-half, quarter}, {half, quarter}},
    };

    const auto diagram = segmentVoronoiDiagram(
        sites, common::Aabb<double>{-extent, -extent, extent, extent});

    ASSERT_FALSE(diagram.edges.empty());
    for (const SegmentVoronoiEdge& edge : diagram.edges) {
        for (const Point& vertex : edge.vertices) {
            EXPECT_TRUE(std::isfinite(vertex.x()));
            EXPECT_TRUE(std::isfinite(vertex.y()));
            EXPECT_GE(vertex.x(), -extent);
            EXPECT_LE(vertex.x(), extent);
            EXPECT_GE(vertex.y(), -extent);
            EXPECT_LE(vertex.y(), extent);
        }
    }
}

TEST(SegmentFortuneVoronoiTest,
     SearchesBeyondCollidingCandidateSweepAngles) {
    const long double pi = std::acos(-1.0L);
    const std::array<long double, 9> collidingAngles{
        -pi * 0.5L,
        -pi * 0.5L + 0.173205080756887729L,
        -pi * 0.5L - 0.223606797749978969L,
        -pi * 0.5L + 0.414213562373095049L,
        -pi * 0.5L - 0.618033988749894848L,
        0.271828182845904523L,
        -0.732050807568877293L,
        1.141592653589793238L,
        -1.287901101718757732L,
    };
    std::vector<BaseSegment<double>> sites;
    for (std::size_t index = 0; index < collidingAngles.size(); ++index) {
        const long double angle = collidingAngles[index];
        const long double transverseX = -std::sin(angle);
        const long double transverseY = std::cos(angle);
        const double centerX =
            20.0 * static_cast<double>(index % 3U);
        const double centerY =
            20.0 * static_cast<double>(index / 3U);
        for (long double side : {-1.0L, 1.0L}) {
            const double startX = centerX +
                static_cast<double>(side * transverseX * 0.5L);
            const double startY = centerY +
                static_cast<double>(side * transverseY * 0.5L);
            const long double segmentAngle =
                angle + 0.37L + side * 0.03L;
            sites.push_back(
                {{startX, startY},
                 {startX + static_cast<double>(
                               std::cos(segmentAngle) * 0.2L),
                  startY + static_cast<double>(
                               std::sin(segmentAngle) * 0.2L)}});
        }
    }

    const auto diagram = segmentVoronoiDiagram(
        sites, common::Aabb<double>{-5.0, -5.0, 45.0, 65.0}, 0.05);

    EXPECT_EQ(diagram.cells.size(), sites.size());
    EXPECT_FALSE(diagram.edges.empty());
}

TEST(SegmentFortuneVoronoiTest, IsInvariantUnderLargeTranslation) {
    const std::vector<BaseSegment<double>> sites{
        {{6.0, 8.0}, {2.0, 3.0}},
        {{4.0, 7.0}, {3.0, 10.0}},
    };
    constexpr double offset = 1e15;
    const std::vector<BaseSegment<double>> translated{
        {{offset + 6.0, offset + 8.0}, {offset + 2.0, offset + 3.0}},
        {{offset + 4.0, offset + 7.0}, {offset + 3.0, offset + 10.0}},
    };

    const auto original =
        segmentVoronoiDiagram(sites, common::Aabb<double>{0.0, 0.0, 10.0, 10.0});
    const auto shifted = segmentVoronoiDiagram(
        translated,
        common::Aabb<double>{offset, offset, offset + 10.0, offset + 10.0});

    EXPECT_EQ(shifted.edges.size(), original.edges.size());
    for (const SegmentVoronoiEdge& edge : shifted.edges) {
        for (const Point& vertex : edge.vertices) {
            EXPECT_GE(vertex.x(), offset);
            EXPECT_LE(vertex.x(), offset + 10.0);
            EXPECT_GE(vertex.y(), offset);
            EXPECT_LE(vertex.y(), offset + 10.0);
        }
    }
}

TEST(SegmentFortuneVoronoiTest, RejectsSharedEndpoints) {
    EXPECT_THROW(
        segmentVoronoiDiagram(
            std::vector<BaseSegment<int>>{
                {{0, 0}, {10, 0}},
                {{0, 0}, {0, 10}},
            },
            common::Aabb<int>{-5, -5, 15, 15}),
        std::invalid_argument);
}

TEST(SegmentFortuneVoronoiTest, RejectsCoordinatesThatLosePrecision) {
    constexpr std::int64_t exactLimit = 9007199254740992LL;
    EXPECT_THROW(
        segmentVoronoiDiagram(
            std::vector<BaseSegment<std::int64_t>>{
                {{exactLimit, 0}, {exactLimit + 1, 1}},
            },
            common::Aabb<std::int64_t>{0, 0, exactLimit + 2, 10}),
        std::invalid_argument);
}

TEST(SegmentFortuneVoronoiTest, BuildsManySitesWithLinearOutputComplexity) {
    constexpr int siteCount = 320;
    std::vector<BaseSegment<double>> sites;
    sites.reserve(siteCount);
    for (int index = 0; index < siteCount; ++index) {
        const double x = static_cast<double>((index * 37) % 101) * 0.03;
        const double y = static_cast<double>(index) * 3.0;
        sites.push_back({{x, y}, {x + 1.0, y}});
    }

    const auto diagram = segmentVoronoiDiagram(
        sites, common::Aabb<double>{-10.0, -2.0, 20.0, 962.0}, 0.1);

    ASSERT_EQ(diagram.cells.size(), sites.size());
    EXPECT_GT(diagram.edges.size(), sites.size());
    EXPECT_LT(diagram.edges.size(), sites.size() * 12);
    for (const SegmentVoronoiCell& cell : diagram.cells) {
        EXPECT_FALSE(cell.edges.empty());
    }
}

} // namespace skeleton::unittest
