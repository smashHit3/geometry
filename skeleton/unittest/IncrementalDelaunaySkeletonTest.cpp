#include "skeleton/detail/IncrementalDelaunaySkeleton.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_set>
#include <vector>

namespace skeleton::unittest {

using skeleton::detail::incrementalDelaunayTriangulation;

namespace {

// Robust in-circumcircle predicate mirroring the builder's logic. Returns true
// when site p lies strictly inside the circumcircle of (a, b, c). Triangle
// orientation is handled by flipping the sign for clockwise inputs.
bool inCircumcircle(const Point& a, const Point& b, const Point& c,
                    const Point& p) {
    const Point da = a - p;
    const Point db = b - p;
    const Point dc = c - p;
    const double adx = da.x();
    const double ady = da.y();
    const double bdx = db.x();
    const double bdy = db.y();
    const double cdx = dc.x();
    const double cdy = dc.y();
    const double a2 = adx * adx + ady * ady;
    const double b2 = bdx * bdx + bdy * bdy;
    const double c2 = cdx * cdx + cdy * cdy;
    const double det = a2 * (bdx * cdy - bdy * cdx) -
                       b2 * (adx * cdy - ady * cdx) +
                       c2 * (adx * bdy - ady * bdx);
    const double orient = common::cross(b - a, c - a);
    const double eps = 1e-9;
    if (orient < -eps) return det < -eps;
    return det > eps;
}

bool isDelaunay(const DelaunayTriangulation& triangulation) {
    const auto& sites = triangulation.sites;
    for (const auto& triangle : triangulation.triangles) {
        const Point a = sites[triangle.first];
        const Point b = sites[triangle.second];
        const Point c = sites[triangle.third];
        if (std::fabs(common::cross(b - a, c - a)) <= 1e-9) continue;
        for (std::size_t i = 0; i < sites.size(); ++i) {
            if (i == triangle.first || i == triangle.second || i == triangle.third) {
                continue;
            }
            if (inCircumcircle(a, b, c, sites[i])) return false;
        }
    }
    return true;
}

std::size_t countTriangles(const DelaunayTriangulation& t) {
    return t.triangles.size();
}

} // namespace

TEST(IncrementalDelaunaySkeletonTest, EmptyInputProducesNoTriangles) {
    const std::vector<Point> sites;
    const auto triangulation = incrementalDelaunayTriangulation(sites);
    EXPECT_EQ(countTriangles(triangulation), 0U);
    EXPECT_TRUE(triangulation.sites.empty());
}

TEST(IncrementalDelaunaySkeletonTest, FewerThanThreeSitesProducesNoTriangles) {
    std::vector<Point> sites{{0.0, 0.0}, {5.0, 5.0}};
    const auto triangulation = incrementalDelaunayTriangulation(sites);
    EXPECT_EQ(countTriangles(triangulation), 0U);
    EXPECT_EQ(triangulation.sites.size(), 2U);
}

TEST(IncrementalDelaunaySkeletonTest, TriangleProducesSingleFace) {
    const std::vector<Point> sites{{0.0, 0.0}, {10.0, 0.0}, {0.0, 10.0}};
    const auto triangulation = incrementalDelaunayTriangulation(sites);
    ASSERT_EQ(countTriangles(triangulation), 1U);
    EXPECT_TRUE(isDelaunay(triangulation));
}

TEST(IncrementalDelaunaySkeletonTest, SquareProducesTwoTriangles) {
    const std::vector<Point> sites{
        {0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}};
    const auto triangulation = incrementalDelaunayTriangulation(sites);
    ASSERT_EQ(countTriangles(triangulation), 2U);
    EXPECT_TRUE(isDelaunay(triangulation));
}

TEST(IncrementalDelaunaySkeletonTest, CollinearSitesProduceNoTriangles) {
    const std::vector<Point> sites{{0.0, 5.0}, {5.0, 5.0}, {10.0, 5.0}};
    const auto triangulation = incrementalDelaunayTriangulation(sites);
    EXPECT_EQ(countTriangles(triangulation), 0U);
}

TEST(IncrementalDelaunaySkeletonTest, CocircularSitesStillDelaunay) {
    const std::vector<Point> sites{
        {2.0, 2.0}, {8.0, 2.0}, {8.0, 8.0}, {2.0, 8.0}};
    const auto triangulation = incrementalDelaunayTriangulation(sites);
    ASSERT_EQ(countTriangles(triangulation), 2U);
    EXPECT_TRUE(isDelaunay(triangulation));
}

TEST(IncrementalDelaunaySkeletonTest, RandomSitesSatisfyDelaunayProperty) {
    std::vector<Point> sites;
    sites.reserve(200);
    std::uint32_t state = 1234567U;
    auto lcg = [&]() {
        state = state * 1103515245U + 12345U;
        return (state >> 8) & 0x1FFFFU; // 0..131071
    };
    for (int i = 0; i < 200; ++i) {
        sites.emplace_back(static_cast<double>(lcg()) * 0.01,
                           static_cast<double>(lcg()) * 0.01);
    }

    const auto triangulation = incrementalDelaunayTriangulation(sites);
    EXPECT_EQ(triangulation.sites.size(), sites.size());
    EXPECT_FALSE(triangulation.triangles.empty());
    EXPECT_TRUE(isDelaunay(triangulation));
}

TEST(IncrementalDelaunaySkeletonTest, NoDuplicateOrDegenerateTriangles) {
    const std::vector<Point> sites{
        {0.0, 0.0}, {4.0, 0.0}, {4.0, 3.0}, {1.0, 1.0},
        {0.0, 3.0}, {2.0, 2.0}, {3.0, 1.5}};
    const auto triangulation = incrementalDelaunayTriangulation(sites);

    std::unordered_set<std::uint64_t> seen;
    for (const auto& triangle : triangulation.triangles) {
        std::array<std::size_t, 3> idx{triangle.first, triangle.second,
                                       triangle.third};
        std::sort(idx.begin(), idx.end());
        ASSERT_TRUE(idx[0] < idx[1] && idx[1] < idx[2]);
        const std::uint64_t key = (idx[0] << 42) | (idx[1] << 21) | idx[2];
        ASSERT_TRUE(seen.insert(key).second);
        const Point a = triangulation.sites[idx[0]];
        const Point b = triangulation.sites[idx[1]];
        const Point c = triangulation.sites[idx[2]];
        EXPECT_GT(std::fabs(common::cross(b - a, c - a)), 1e-9);
    }
}

TEST(IncrementalDelaunaySkeletonTest, SatisfiesEulerRelationForFullSet) {
    // For a triangulation of n points with h on the convex hull:
    //   triangles = 2n - 2 - h,  edges = 3n - 3 - h.
    std::vector<Point> sites{{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0},
                              {0.0, 10.0}, {3.0, 3.0}, {7.0, 3.0},
                              {5.0, 7.0}};
    const auto triangulation = incrementalDelaunayTriangulation(sites);
    const std::size_t n = sites.size();

    std::unordered_set<std::uint64_t> edges;
    for (const auto& triangle : triangulation.triangles) {
        const std::array<std::size_t, 3> v{triangle.first, triangle.second,
                                           triangle.third};
        for (int e = 0; e < 3; ++e) {
            std::size_t a = v[e];
            std::size_t b = v[(e + 1) % 3];
            if (a > b) std::swap(a, b);
            edges.insert((static_cast<std::uint64_t>(a) << 32) | b);
        }
    }
    // Hull edges are a subset of triangulation edges; count hull vertices via
    // edges that appear in exactly one triangle (boundary edges).
    std::unordered_map<std::uint64_t, int> edge_count;
    for (const auto& triangle : triangulation.triangles) {
        const std::array<std::size_t, 3> v{triangle.first, triangle.second,
                                            triangle.third};
        for (int e = 0; e < 3; ++e) {
            std::size_t a = v[e];
            std::size_t b = v[(e + 1) % 3];
            if (a > b) std::swap(a, b);
            edge_count[(static_cast<std::uint64_t>(a) << 32) | b]++;
        }
    }
    std::size_t boundary_edges = 0;
    for (const auto& pair : edge_count) {
        if (pair.second == 1) ++boundary_edges;
    }
    const std::size_t h = boundary_edges;
    EXPECT_EQ(triangulation.triangles.size(), 2 * n - 2 - h);
    EXPECT_EQ(edges.size(), 3 * n - 3 - h);
}

} // namespace skeleton::unittest
