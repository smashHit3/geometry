#include "skeleton/PolygonFortuneVoronoi.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace skeleton::unittest {
namespace {

double featureDistance(const Point& point,
                       const std::vector<Point>& polygon,
                       const PolygonFeatureRef& feature) {
    if (feature.kind == PolygonFeatureRef::Kind::Vertex) {
        return common::geometryutil::distance(point, polygon[feature.inputIndex]);
    }
    const Point& first = polygon[feature.inputIndex];
    const Point& second =
        polygon[(feature.inputIndex + 1U) % polygon.size()];
    const Point direction = second - first;
    const double parameter = std::clamp(
        common::geometryutil::dot(point - first, direction) /
            common::geometryutil::squaredLength(direction),
        0.0, 1.0);
    return common::geometryutil::distance(
        point, first + direction * parameter);
}

bool insideOrOnPolygon(const Point& point,
                       const std::vector<Point>& polygon) {
    bool inside = false;
    for (std::size_t first = 0, second = polygon.size() - 1U;
         first < polygon.size(); second = first++) {
        const Point edge = polygon[first] - polygon[second];
        const Point relative = point - polygon[second];
        const double cross = common::geometryutil::cross(edge, relative);
        const double projection = common::geometryutil::dot(relative, edge);
        if (std::abs(cross) <= 1e-7 &&
            projection >= -1e-7 &&
            projection <= common::geometryutil::squaredLength(edge) + 1e-7) {
            return true;
        }
        const bool crosses =
            (polygon[first].y() > point.y()) !=
            (polygon[second].y() > point.y());
        if (crosses) {
            const double intersectionX =
                polygon[second].x() +
                (point.y() - polygon[second].y()) *
                    (polygon[first].x() - polygon[second].x()) /
                    (polygon[first].y() - polygon[second].y());
            if (point.x() < intersectionX) inside = !inside;
        }
    }
    return inside;
}

void expectValidFeatureEdges(const PolygonVoronoiDiagram& diagram,
                             const std::vector<Point>& polygon,
                             double tolerance) {
    for (const PolygonVoronoiEdge& edge : diagram.edges) {
        EXPECT_LT(edge.firstFeature.inputIndex, polygon.size());
        EXPECT_LT(edge.secondFeature.inputIndex, polygon.size());
        ASSERT_GE(edge.vertices.size(), 2U);
        for (const Point& vertex : edge.vertices) {
            EXPECT_NEAR(featureDistance(vertex, polygon, edge.firstFeature),
                        featureDistance(vertex, polygon, edge.secondFeature),
                        tolerance);
        }
    }
}

} // namespace

TEST(PolygonFortuneVoronoiTest, BuildsTriangleBoundaryFeatureDiagram) {
    const std::vector<Point> polygon{{0.0, 0.0}, {10.0, 0.0}, {5.0, 8.0}};
    const auto diagram = polygonBoundaryVoronoiDiagram(
        polygon, common::geometry::Aabb<double>{-2.0, -2.0, 12.0, 10.0}, 1e-4);

    ASSERT_GE(diagram.edges.size(), 3U);
    expectValidFeatureEdges(diagram, polygon, 2e-2);
    for (std::size_t edge = 0; edge < polygon.size(); ++edge) {
        const std::size_t previous =
            (edge + polygon.size() - 1U) % polygon.size();
        EXPECT_TRUE(std::any_of(
            diagram.edges.begin(), diagram.edges.end(),
            [&](const PolygonVoronoiEdge& value) {
                if (value.firstFeature.kind !=
                        PolygonFeatureRef::Kind::Edge ||
                    value.secondFeature.kind !=
                        PolygonFeatureRef::Kind::Edge) {
                    return false;
                }
                return (value.firstFeature.inputIndex == previous &&
                        value.secondFeature.inputIndex == edge) ||
                       (value.firstFeature.inputIndex == edge &&
                        value.secondFeature.inputIndex == previous);
            }));
    }
}

TEST(PolygonFortuneVoronoiTest, BuildsSquareMedialAxis) {
    const std::vector<Point> polygon{
        {0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}};
    const auto axis = polygonMedialAxis(polygon, 1e-4);

    ASSERT_EQ(axis.edges.size(), 4U);
    expectValidFeatureEdges(axis, polygon, 2e-2);
    for (const PolygonVoronoiEdge& edge : axis.edges) {
        EXPECT_EQ(edge.firstFeature.kind, PolygonFeatureRef::Kind::Edge);
        EXPECT_EQ(edge.secondFeature.kind, PolygonFeatureRef::Kind::Edge);
        for (const Point& vertex : edge.vertices) {
            EXPECT_TRUE(insideOrOnPolygon(vertex, polygon));
        }
        const Point& end = edge.vertices.back();
        EXPECT_NEAR(end.x(), 5.0, 1e-2);
        EXPECT_NEAR(end.y(), 5.0, 1e-2);
    }
}

TEST(PolygonFortuneVoronoiTest, ClipsBoundaryDiagramToSuppliedBounds) {
    const std::vector<Point> polygon{
        {0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}};
    const common::geometry::Aabb<double> bounds{2.0, 2.0, 8.0, 8.0};

    const auto diagram =
        polygonBoundaryVoronoiDiagram(polygon, bounds, 1e-4);

    ASSERT_EQ(diagram.edges.size(), 4U);
    for (const PolygonVoronoiEdge& edge : diagram.edges) {
        for (const Point& vertex : edge.vertices) {
            EXPECT_GE(vertex.x(), bounds.minX() - 1e-8);
            EXPECT_LE(vertex.x(), bounds.maxX() + 1e-8);
            EXPECT_GE(vertex.y(), bounds.minY() - 1e-8);
            EXPECT_LE(vertex.y(), bounds.maxY() + 1e-8);
        }
    }
}

TEST(PolygonFortuneVoronoiTest, RetainsReflexVertexParabola) {
    const std::vector<Point> polygon{
        {0.0, 0.0}, {8.0, 0.0}, {8.0, 8.0},
        {4.0, 3.0}, {0.0, 8.0}};
    const auto axis = polygonMedialAxis(polygon, 1e-4);

    expectValidFeatureEdges(axis, polygon, 2e-2);
    const auto reflex = std::find_if(
        axis.edges.begin(), axis.edges.end(),
        [](const PolygonVoronoiEdge& edge) {
            const auto isReflexVertex =
                [](const PolygonFeatureRef& feature) {
                    return feature.kind ==
                               PolygonFeatureRef::Kind::Vertex &&
                           feature.inputIndex == 3U;
                };
            return edge.curved &&
                   (isReflexVertex(edge.firstFeature) ||
                    isReflexVertex(edge.secondFeature));
        });
    ASSERT_NE(reflex, axis.edges.end());
    EXPECT_GT(reflex->vertices.size(), 2U);
    for (const Point& vertex : reflex->vertices) {
        EXPECT_TRUE(insideOrOnPolygon(vertex, polygon));
    }
}

TEST(PolygonFortuneVoronoiTest, IsIndependentOfWinding) {
    const std::vector<Point> counterClockwise{
        {0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}};
    std::vector<Point> clockwise = counterClockwise;
    std::reverse(clockwise.begin(), clockwise.end());

    const auto first = polygonMedialAxis(counterClockwise, 1e-4);
    const auto second = polygonMedialAxis(clockwise, 1e-4);

    ASSERT_EQ(first.edges.size(), 4U);
    ASSERT_EQ(second.edges.size(), first.edges.size());
    for (const PolygonVoronoiEdge& edge : second.edges) {
        for (const Point& vertex : edge.vertices) {
            EXPECT_TRUE(insideOrOnPolygon(vertex, clockwise));
        }
    }
}

TEST(PolygonFortuneVoronoiTest, RejectsNonSimplePolygons) {
    const common::geometry::Aabb<double> bounds{-1.0, -1.0, 6.0, 6.0};
    EXPECT_THROW(
        polygonBoundaryVoronoiDiagram(
            std::vector<Point>{{0.0, 0.0}, {4.0, 4.0},
                               {0.0, 4.0}, {4.0, 0.0}},
            bounds),
        std::invalid_argument);
    EXPECT_THROW(
        polygonBoundaryVoronoiDiagram(
            std::vector<Point>{{0.0, 0.0}, {4.0, 0.0},
                               {4.0, 4.0}, {0.0, 4.0},
                               {4.0, 0.0}},
            bounds),
        std::invalid_argument);
    EXPECT_THROW(
        polygonBoundaryVoronoiDiagram(
            std::vector<Point>{{0.0, 0.0}, {4.0, 0.0},
                               {4.0, 4.0}, {2.0, 0.0},
                               {0.0, 4.0}},
            bounds),
        std::invalid_argument);
    EXPECT_THROW(
        polygonBoundaryVoronoiDiagram(
            std::vector<Point>{{0.0, 0.0}, {4.0, 0.0},
                               {2.0, 0.0}, {3.0, 0.0},
                               {0.0, 4.0}},
            bounds),
        std::invalid_argument);
}

} // namespace skeleton::unittest
