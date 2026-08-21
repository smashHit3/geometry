#include "skeleton/StraightSkeleton.h"

#include "common/BasePointUtil.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

namespace skeleton::unittest {

namespace {

bool insidePolygon(const std::vector<Point>& polygon, const Point& p) {
    // Ray casting for a CCW simple polygon.
    bool result = false;
    const std::size_t n = polygon.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const Point a = polygon[i];
        const Point b = polygon[j];
        const bool intersects = (a.y() > p.y()) != (b.y() > p.y());
        if (!intersects) continue;
        const double x_intersect =
            a.x() + (p.y() - a.y()) / (b.y() - a.y()) * (b.x() - a.x());
        if (p.x() < x_intersect) result = !result;
    }
    return result;
}

double polygonArea(const std::vector<Point>& polygon) {
    double area = 0.0;
    const std::size_t n = polygon.size();
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        area += static_cast<double>(polygon[i].x()) *
                    static_cast<double>(polygon[j].y()) -
                static_cast<double>(polygon[i].y()) *
                    static_cast<double>(polygon[j].x());
    }
    return area * 0.5;
}

// The straight skeleton of a simple polygon is a tree whose leaves are the
// polygon vertices, so its number of arcs equals its number of nodes minus one.
// Union-find is used to verify the underlying graph is connected and acyclic.
struct DisjointSets {
    explicit DisjointSets(std::size_t n) : parent(n), rank(n, 0) {
        for (std::size_t i = 0; i < n; ++i) parent[i] = i;
    }
    std::size_t find(std::size_t i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    }
    bool unite(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) ++rank[a];
        return true;
    }
    std::vector<std::size_t> parent;
    std::vector<std::size_t> rank;
};

void verifySkeletonInvariants(const StraightSkeleton& skeleton,
                              const std::vector<Point>& polygon) {
    SCOPED_TRACE("Verifying skeleton invariants");
    ASSERT_FALSE(skeleton.nodes.empty());
    ASSERT_FALSE(skeleton.arcs.empty());

    // All arc endpoints must be valid node indices.
    for (const auto& arc : skeleton.arcs) {
        ASSERT_LT(arc.from, skeleton.nodes.size());
        ASSERT_LT(arc.to, skeleton.nodes.size());
        ASSERT_NE(arc.from, arc.to);
        const double dist =
            common::squaredDistance(skeleton.nodes[arc.from], skeleton.nodes[arc.to]);
        EXPECT_GT(dist, 1e-10);
    }

    // The skeleton is a connected forest (tree) on its nodes.
    DisjointSets ds(skeleton.nodes.size());
    std::size_t components = skeleton.nodes.size();
    for (const auto& arc : skeleton.arcs) {
        if (ds.unite(arc.from, arc.to)) {
            --components;
        }
    }
    EXPECT_EQ(components, 1U);
    EXPECT_EQ(skeleton.arcs.size(), skeleton.nodes.size() - 1U);

    // Every original polygon vertex must appear as a skeleton node.
    for (const Point& vertex : polygon) {
        bool found = false;
        for (const Point& node : skeleton.nodes) {
            if (common::squaredDistance(node, vertex) <= 1e-10) {
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found);
    }

    // Interior skeleton nodes (non-boundary) must lie inside the polygon.
    for (const Point& node : skeleton.nodes) {
        bool is_boundary = false;
        for (const Point& vertex : polygon) {
            if (common::squaredDistance(node, vertex) <= 1e-10) {
                is_boundary = true;
                break;
            }
        }
        if (!is_boundary) {
            EXPECT_TRUE(insidePolygon(polygon, node));
        }
    }
}

} // namespace

TEST(StraightSkeletonTest, DegenerateInputProducesEmptySkeleton) {
    std::vector<Point> polygon{{0.0, 0.0}, {5.0, 5.0}};
    const auto skeleton = straightSkeleton(polygon);
    EXPECT_TRUE(skeleton.nodes.empty());
    EXPECT_TRUE(skeleton.arcs.empty());
}

TEST(StraightSkeletonTest, TriangleSkeletonMeetsAtIncenter) {
    const std::vector<Point> polygon{
        {0.0, 0.0}, {10.0, 0.0}, {0.0, 10.0}};
    const auto skeleton = straightSkeleton(polygon);
    verifySkeletonInvariants(skeleton, polygon);

    // The straight skeleton of a triangle collapses at the incenter.
    // For a right triangle with legs 10 the inradius r = (a+b-c)/2 with
    // c=10*sqrt(2), giving r = 10/(2+sqrt(2)) ~ 2.9289.
    const double inradius = 10.0 / (2.0 + std::sqrt(2.0));
    const Point incenter{inradius, inradius};
    bool found = false;
    for (const Point& node : skeleton.nodes) {
        if (common::squaredDistance(node, incenter) <= 1e-6) {
            found = true;
            break;
        }
    }
    EXPECT_TRUE(found);
}

TEST(StraightSkeletonTest, SquareSkeletonMeetsAtCenter) {
    const std::vector<Point> polygon{
        {0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}};
    const auto skeleton = straightSkeleton(polygon);
    verifySkeletonInvariants(skeleton, polygon);

    const Point center{5.0, 5.0};
    bool found = false;
    for (const Point& node : skeleton.nodes) {
        if (common::squaredDistance(node, center) <= 1e-6) {
            found = true;
            break;
        }
    }
    EXPECT_TRUE(found);
}

TEST(StraightSkeletonTest, ConvexPentagonFormsTree) {
    const std::vector<Point> polygon{
        {0.0, 0.0}, {10.0, 0.0}, {12.0, 6.0}, {5.0, 10.0}, {0.0, 6.0}};
    ASSERT_GT(polygonArea(polygon), 0.0); // ensure CCW
    const auto skeleton = straightSkeleton(polygon);
    verifySkeletonInvariants(skeleton, polygon);
    // A generic convex polygon with n vertices has n-2 interior nodes.
    EXPECT_EQ(skeleton.nodes.size(), polygon.size() + polygon.size() - 2U);
}

TEST(StraightSkeletonTest, ReflexVertexTriggersSplitEvent) {
    // A hexagon with a single inward dent whose reflex vertex (7,5) has a
    // nearly-vertical (non-45-degree) bisector. The bisector pierces the bottom
    // wavefront edge at t~2.31, a time distinct from every adjacent edge
    // collapse, so the split event is genuinely non-degenerate.
    const std::vector<Point> polygon{
        {0.0, 0.0}, {12.0, 0.0}, {12.0, 8.0}, {7.0, 5.0},
        {2.0, 8.0}, {0.0, 8.0}};
    ASSERT_GT(polygonArea(polygon), 0.0);
    const auto skeleton = straightSkeleton(polygon);
    verifySkeletonInvariants(skeleton, polygon);

    // A generic (non-degenerate) straight skeleton is a tree whose n boundary
    // vertices are leaves and whose interior nodes all have degree 3, giving
    // exactly 2n-2 nodes total. The reflex vertex's split point is one of
    // those degree-3 interior nodes, so reaching this count with a valid tree
    // confirms the split event was handled correctly.
    EXPECT_EQ(skeleton.nodes.size(), 2 * polygon.size() - 2U);
}

TEST(StraightSkeletonTest, ClockwiseInputIsReoriented) {
    const std::vector<Point> cw_polygon{
        {0.0, 0.0}, {0.0, 10.0}, {10.0, 10.0}, {10.0, 0.0}};
    ASSERT_LT(polygonArea(cw_polygon), 0.0);
    const auto skeleton = straightSkeleton(cw_polygon);
    verifySkeletonInvariants(skeleton, cw_polygon);

    // The result must match the CCW square skeleton.
    const std::vector<Point> ccw_polygon{
        {0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}};
    const auto reference = straightSkeleton(ccw_polygon);
    EXPECT_EQ(skeleton.nodes.size(), reference.nodes.size());
    EXPECT_EQ(skeleton.arcs.size(), reference.arcs.size());
}

TEST(StraightSkeletonTest, LargerConcavePolygonStaysConsistent) {
    // A heptagon with two inward dents whose reflex vertices (14,6) and (5,5)
    // both have nearly-vertical bisectors piercing the bottom wavefront edge at
    // distinct times (t~2.19 and t~2.67) and distinct points, keeping both
    // split events genuinely non-degenerate.
    const std::vector<Point> polygon{
        {0.0, 0.0}, {18.0, 0.0}, {18.0, 9.0}, {14.0, 6.0},
        {10.0, 9.0}, {5.0, 5.0}, {0.0, 9.0}};
    ASSERT_GT(polygonArea(polygon), 0.0);
    const auto skeleton = straightSkeleton(polygon);
    verifySkeletonInvariants(skeleton, polygon);
}

} // namespace skeleton::unittest
