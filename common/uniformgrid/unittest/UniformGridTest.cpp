#include "common/uniformgrid/UniformGrid.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

namespace common::uniformgrid::unittest {
namespace {

std::vector<int> sorted(std::vector<int> ids) {
    std::sort(ids.begin(), ids.end());
    return ids;
}

} // namespace

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

} // namespace common::uniformgrid::unittest
