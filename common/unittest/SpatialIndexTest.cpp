#include "common/Aabb.h"
#include "common/Bvh.h"
#include "common/RTree.h"
#include "common/UniformGrid.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace common::unittest {
namespace {

using Bounds = Aabb<double>;

std::vector<int> sorted(std::vector<int> ids) {
    std::sort(ids.begin(), ids.end());
    return ids;
}

} // namespace

TEST(AabbTest, TestsIntersectionContainmentAndMerge) {
    const Bounds bounds{1.0, 2.0, 4.0, 6.0};
    const Bounds overlapping{4.0, 5.0, 8.0, 9.0};
    const Bounds disjoint{5.0, 7.0, 8.0, 9.0};

    EXPECT_TRUE(bounds.intersects(overlapping));
    EXPECT_FALSE(bounds.intersects(disjoint));
    EXPECT_TRUE(bounds.contains(Bounds{2.0, 3.0, 3.0, 5.0}));
    EXPECT_EQ(bounds.merged(overlapping).maxX(), 8.0);
    EXPECT_DOUBLE_EQ(bounds.area(), 12.0);
    EXPECT_THROW((Bounds{2.0, 1.0, 1.0, 3.0}), std::invalid_argument);
}

TEST(UniformGridTest, QueriesAndUpdatesEntriesAcrossCells) {
    UniformGrid<int, double> grid{2.0};
    grid.insert(1, {0.0, 0.0, 1.0, 1.0});
    grid.insert(2, {1.5, 0.0, 3.0, 1.0});
    grid.insert(3, {8.0, 8.0, 9.0, 9.0});

    EXPECT_EQ(sorted(grid.query({0.5, 0.5, 2.0, 1.0})),
              (std::vector<int>{1, 2}));
    grid.update(3, {2.5, 0.0, 3.5, 1.0});
    EXPECT_EQ(sorted(grid.query({2.0, 0.0, 4.0, 1.0})),
              (std::vector<int>{2, 3}));
    EXPECT_TRUE(grid.remove(1));
    EXPECT_FALSE(grid.remove(1));
}

TEST(BvhTest, BuildsAndQueriesStaticEntries) {
    Bvh<int, double> bvh{1};
    bvh.build({{1, {0.0, 0.0, 1.0, 1.0}},
               {2, {3.0, 0.0, 4.0, 1.0}},
               {3, {8.0, 8.0, 9.0, 9.0}}});

    EXPECT_EQ(sorted(bvh.query({0.5, 0.5, 3.5, 1.0})),
              (std::vector<int>{1, 2}));
    EXPECT_TRUE(bvh.query({5.0, 5.0, 6.0, 6.0}).empty());
}

TEST(RTreeTest, InsertsQueriesUpdatesAndRemovesEntries) {
    RTree<int, double> tree{2};
    tree.insert(1, {0.0, 0.0, 1.0, 1.0});
    tree.insert(2, {3.0, 0.0, 4.0, 1.0});
    tree.insert(3, {8.0, 8.0, 9.0, 9.0});

    EXPECT_EQ(sorted(tree.query({0.5, 0.5, 3.5, 1.0})),
              (std::vector<int>{1, 2}));
    tree.update(3, {2.5, 0.0, 3.5, 1.0});
    EXPECT_EQ(sorted(tree.query({2.0, 0.0, 4.0, 1.0})),
              (std::vector<int>{2, 3}));
    EXPECT_TRUE(tree.remove(2));
    EXPECT_EQ(sorted(tree.query({2.0, 0.0, 4.0, 1.0})),
              (std::vector<int>{3}));
    EXPECT_FALSE(tree.remove(2));
}

} // namespace common::unittest
