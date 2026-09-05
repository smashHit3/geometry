#include "common/rtree/RTree.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

namespace common::rtree::unittest {
namespace {

std::vector<int> sorted(std::vector<int> ids) {
    std::sort(ids.begin(), ids.end());
    return ids;
}

} // namespace

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

} // namespace common::rtree::unittest
