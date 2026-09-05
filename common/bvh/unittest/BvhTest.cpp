#include "common/bvh/Bvh.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

namespace common::bvh::unittest {
namespace {

std::vector<int> sorted(std::vector<int> ids) {
    std::sort(ids.begin(), ids.end());
    return ids;
}

} // namespace

TEST(BvhTest, BuildsAndQueriesStaticEntries) {
    Bvh<int, double> bvh{1};
    bvh.build({{1, {0.0, 0.0, 1.0, 1.0}},
               {2, {3.0, 0.0, 4.0, 1.0}},
               {3, {8.0, 8.0, 9.0, 9.0}}});

    EXPECT_EQ(sorted(bvh.query({0.5, 0.5, 3.5, 1.0})),
              (std::vector<int>{1, 2}));
    EXPECT_TRUE(bvh.query({5.0, 5.0, 6.0, 6.0}).empty());
}

} // namespace common::bvh::unittest
