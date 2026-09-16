#include "trapezoid/OctilinearTrapezoidMap.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <vector>

namespace trapezoid::unittest {
namespace {

using Map = OctilinearTrapezoidMap<long long>;
using Segment = Map::Segment;
using Point = Map::Point;

Segment seg(long long x1, long long y1, long long x2, long long y2) {
    return {Point(x1, y1), Point(x2, y2)};
}

// Axis-aligned rectangle boundary (CCW).
std::vector<Segment> rectangle(long long x1, long long y1, long long x2,
                               long long y2) {
    return {seg(x1, y1, x2, y1), seg(x2, y1, x2, y2), seg(x2, y2, x1, y2),
            seg(x1, y2, x1, y1)};
}

// Diamond (45-degree square) boundary.
std::vector<Segment> diamond(long long cx, long long cy, long long radius) {
    return {seg(cx, cy - radius, cx + radius, cy),
            seg(cx + radius, cy, cx, cy + radius),
            seg(cx, cy + radius, cx - radius, cy),
            seg(cx - radius, cy, cx, cy - radius)};
}

long long cellArea(const Map::Cell& cell) {
    const long long width = cell.xmax - cell.xmin;
    const long long height_sum =
        (cell.yUpperLeft - cell.yLowerLeft) +
        (cell.yUpperRight - cell.yLowerRight);
    return width * height_sum; // doubled area
}

void expectValidDcel(const Map& map) {
    ASSERT_EQ(map.faces().size(), map.cells().size() + 1);
    ASSERT_EQ(map.exteriorFace(), map.cells().size());
    ASSERT_FALSE(map.vertices().empty());
    ASSERT_FALSE(map.halfEdges().empty());

    for (std::size_t i = 0; i < map.halfEdges().size(); ++i) {
        const Map::HalfEdge& edge = map.halfEdges()[i];
        ASSERT_LT(edge.origin, map.vertices().size());
        ASSERT_LT(edge.twin, map.halfEdges().size());
        ASSERT_LT(edge.next, map.halfEdges().size());
        ASSERT_LT(edge.previous, map.halfEdges().size());
        ASSERT_LT(edge.face, map.faces().size());
        EXPECT_EQ(map.halfEdges()[edge.twin].twin, i);
        EXPECT_EQ(map.halfEdges()[edge.next].previous, i);
        EXPECT_EQ(map.halfEdges()[edge.previous].next, i);
    }

    for (std::size_t face = 0; face < map.faces().size(); ++face) {
        const std::size_t first = map.faces()[face].incidentHalfEdge;
        ASSERT_LT(first, map.halfEdges().size());
        std::size_t edge = first;
        std::size_t steps = 0;
        do {
            ASSERT_EQ(map.halfEdges()[edge].face, face);
            edge = map.halfEdges()[edge].next;
            ASSERT_LE(++steps, map.halfEdges().size());
        } while (edge != first);
    }
}

} // namespace

TEST(OctilinearTrapezoidMapTest, DecomposesSingleRectangle) {
    Map map(rectangle(0, 0, 10, 8));

    expectValidDcel(map);
    EXPECT_EQ(map.cells().size(), 5u);
    EXPECT_EQ(map.faceCount(), 2u);

    const std::int64_t inner = map.locateFace(5, 4);
    EXPECT_EQ(inner, map.locateFace(1, 1));
    const std::int64_t outer = map.locateFace(5, -1);
    EXPECT_EQ(outer, map.locateFace(-1, 4));
    EXPECT_EQ(outer, map.locateFace(11, 4));
    EXPECT_NE(inner, outer);
    EXPECT_EQ(map.locate(-100, -100), -1);

    std::size_t input_half_edges = 0;
    for (const Map::HalfEdge& edge : map.halfEdges()) {
        input_half_edges += edge.isInputEdge;
    }
    EXPECT_EQ(input_half_edges, 8u);
}

TEST(OctilinearTrapezoidMapTest, DecomposesDiamond) {
    Map map(diamond(2, 2, 2)); // (2,0)-(4,2)-(2,4)-(0,2)

    expectValidDcel(map);
    EXPECT_EQ(map.faceCount(), 2u);
    const std::int64_t inside = map.locateFace(2, 2);
    EXPECT_EQ(inside, map.locateFace(2, 1));
    EXPECT_EQ(inside, map.locateFace(3, 2));
    EXPECT_EQ(inside, map.locateFace(1, 2));

    const std::int64_t outside = map.locateFace(4, 3);
    EXPECT_EQ(outside, map.locateFace(0, 0));
    EXPECT_EQ(outside, map.locateFace(2, 5));
    EXPECT_NE(inside, outside);
}

TEST(OctilinearTrapezoidMapTest, DecomposesHouseShapedPentagon) {
    // Rectangle body plus a 45-degree roof.
    const std::vector<Segment> house = {
        seg(0, 0, 4, 0), seg(4, 0, 4, 2), seg(4, 2, 2, 4),
        seg(2, 4, 0, 2), seg(0, 2, 0, 0)};
    Map map(house);

    expectValidDcel(map);
    const std::int64_t inside = map.locateFace(2, 1);
    EXPECT_EQ(inside, map.locateFace(3, 1));
    EXPECT_EQ(inside, map.locateFace(1, 1));
    EXPECT_EQ(inside, map.locateFace(2, 3)); // roof interior

    const std::int64_t outside = map.locateFace(0, 3);
    EXPECT_EQ(outside, map.locateFace(4, 3));
    EXPECT_EQ(outside, map.locateFace(2, 5));
    EXPECT_NE(inside, outside);
}

TEST(OctilinearTrapezoidMapTest, PartitionsWithoutGaps) {
    Map map(rectangle(0, 0, 10, 8));
    long long doubled_area = 0;
    for (const auto& cell : map.cells()) {
        doubled_area += cellArea(cell);
    }
    // Auto margin 1 -> 12 x 10 bounding box, doubled area 240.
    EXPECT_EQ(doubled_area, 240);

    Map diamond_map(diamond(2, 2, 2));
    const auto bounds = diamond_map.bounds();
    long long diamond_box_doubled =
        2 * (bounds.xmax - bounds.xmin) * (bounds.ymax - bounds.ymin);
    long long diamond_doubled = 0;
    for (const auto& cell : diamond_map.cells()) {
        diamond_doubled += cellArea(cell);
    }
    EXPECT_EQ(diamond_doubled, diamond_box_doubled);
}

TEST(OctilinearTrapezoidMapTest, RejectsInvalidInputs) {
    const std::vector<Segment> zero_length = {seg(0, 0, 0, 0)};
    EXPECT_THROW(Map{zero_length}, std::invalid_argument);

    // Unsupported slope (1/2).
    const std::vector<Segment> bad_slope = {seg(0, 0, 2, 1)};
    EXPECT_THROW(Map{bad_slope}, std::invalid_argument);

    // Two diagonals crossing in an X.
    const std::vector<Segment> crossing_diagonals = {
        seg(0, 0, 4, 4), seg(0, 4, 4, 0)};
    EXPECT_THROW(Map{crossing_diagonals}, std::invalid_argument);

    // Diagonal crossing a horizontal through both interiors.
    const std::vector<Segment> crossing_horizontal = {
        seg(0, 2, 4, 2), seg(1, 1, 3, 3)};
    EXPECT_THROW(Map{crossing_horizontal}, std::invalid_argument);

    // Diagonal crossing a vertical through both interiors.
    const std::vector<Segment> crossing_vertical = {
        seg(0, 0, 4, 4), seg(2, 0, 2, 3)};
    EXPECT_THROW(Map{crossing_vertical}, std::invalid_argument);

    // Collinear overlap.
    const std::vector<Segment> overlap = {seg(0, 0, 3, 3), seg(2, 2, 5, 5)};
    EXPECT_THROW(Map{overlap}, std::invalid_argument);
}

TEST(OctilinearTrapezoidMapTest, AllowsEndpointContacts) {
    // A horizontal segment ends exactly on the interior of a 45-degree edge.
    const std::vector<Segment> t_junction = {
        seg(0, 0, 4, 4), seg(0, 2, 2, 2)};
    EXPECT_NO_THROW(Map{t_junction});

    // A chevron made of two diagonals meeting at a vertex.
    const std::vector<Segment> chevron = {seg(0, 0, 2, 2), seg(2, 2, 4, 0)};
    EXPECT_NO_THROW(Map{chevron});
}

TEST(OctilinearTrapezoidMapTest, MatchesBruteForceOnUniformRectangles) {
    std::mt19937 rng(987654321u);
    for (int trial = 0; trial < 20; ++trial) {
        const int rows = 3 + rng() % 5;
        const int cols = 3 + rng() % 5;
        constexpr long long pitch = 6;
        std::vector<Segment> segments;
        std::vector<std::pair<long long, long long>> mins;
        std::vector<std::pair<long long, long long>> maxs;
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                const long long x1 = c * pitch + 1;
                const long long y1 = r * pitch + 1;
                const long long x2 = x1 + 2 + rng() % 2;
                const long long y2 = y1 + 2 + rng() % 2;
                auto rect = rectangle(x1, y1, x2, y2);
                segments.insert(segments.end(), rect.begin(), rect.end());
                mins.emplace_back(x1, y1);
                maxs.emplace_back(x2, y2);
            }
        }
        Map map(segments);
        const auto bounds = map.bounds();
        std::vector<std::int64_t> owner_by_face(map.faceCount(), -2);

        for (long long x = bounds.xmin; x <= bounds.xmax; ++x) {
            for (long long y = bounds.ymin; y <= bounds.ymax; ++y) {
                bool on_boundary = false;
                std::int64_t owner = -1;
                for (std::size_t k = 0; k < mins.size(); ++k) {
                    const bool inside =
                        x >= mins[k].first && x <= maxs[k].first &&
                        y >= mins[k].second && y <= maxs[k].second;
                    const bool strict =
                        x > mins[k].first && x < maxs[k].first &&
                        y > mins[k].second && y < maxs[k].second;
                    if (strict) {
                        owner = static_cast<std::int64_t>(k);
                    }
                    if (inside && !strict) {
                        on_boundary = true;
                    }
                }
                if (on_boundary) {
                    continue;
                }
                const std::int64_t cell_id = map.locate(x, y);
                ASSERT_GE(cell_id, 0) << "gap at (" << x << "," << y << ")";
                const std::size_t face = map.cells()[cell_id].face;
                std::int64_t& slot = owner_by_face[face];
                if (slot == -2) {
                    slot = owner;
                } else {
                    EXPECT_EQ(slot, owner)
                        << "face merges distinct regions";
                }
            }
        }
    }
}

TEST(OctilinearTrapezoidMapTest, MatchesBruteForceOnUniformDiamonds) {
    std::mt19937 rng(135792468u);
    for (int trial = 0; trial < 20; ++trial) {
        const int rows = 2 + rng() % 4;
        const int cols = 2 + rng() % 4;
        constexpr long long pitch = 8;
        constexpr long long radius = 3;
        std::vector<Segment> segments;
        std::vector<std::pair<long long, long long>> centers;
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                const long long cx = c * pitch + pitch / 2;
                const long long cy = r * pitch + pitch / 2;
                auto d = diamond(cx, cy, radius);
                segments.insert(segments.end(), d.begin(), d.end());
                centers.emplace_back(cx, cy);
            }
        }
        Map map(segments);
        const auto bounds = map.bounds();
        std::vector<std::int64_t> owner_by_face(map.faceCount(), -2);

        for (long long x = bounds.xmin; x <= bounds.xmax; ++x) {
            for (long long y = bounds.ymin; y <= bounds.ymax; ++y) {
                std::int64_t owner = -1;
                bool on_boundary = false;
                for (std::size_t k = 0; k < centers.size(); ++k) {
                    const long long l1 =
                        std::llabs(x - centers[k].first) +
                        std::llabs(y - centers[k].second);
                    if (l1 < radius) {
                        owner = static_cast<std::int64_t>(k);
                    } else if (l1 == radius) {
                        on_boundary = true;
                    }
                }
                if (on_boundary) {
                    continue;
                }
                const std::int64_t cell_id = map.locate(x, y);
                ASSERT_GE(cell_id, 0) << "gap at (" << x << "," << y << ")";
                const std::size_t face = map.cells()[cell_id].face;
                std::int64_t& slot = owner_by_face[face];
                if (slot == -2) {
                    slot = owner;
                } else {
                    EXPECT_EQ(slot, owner)
                        << "face merges distinct regions";
                }
            }
        }
    }
}

TEST(OctilinearTrapezoidMapTest, MatchesBruteForceOnMixedLayout) {
    std::mt19937 rng(246813579u);
    for (int trial = 0; trial < 15; ++trial) {
        const int rows = 3 + rng() % 4;
        const int cols = 3 + rng() % 4;
        constexpr long long pitch = 8;
        std::vector<Segment> segments;
        std::vector<std::array<long long, 5>> shapes; // cx,cy or x1,y1,x2,y2
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                if ((r + c + trial) % 2 == 0) {
                    auto rect =
                        rectangle(c * pitch + 1, r * pitch + 1,
                                  c * pitch + 6, r * pitch + 6);
                    segments.insert(segments.end(), rect.begin(), rect.end());
                    shapes.push_back({0, c * pitch + 1, r * pitch + 1,
                                      c * pitch + 6, r * pitch + 6});
                } else {
                    const long long cx = c * pitch + 4;
                    const long long cy = r * pitch + 4;
                    auto d = diamond(cx, cy, 2);
                    segments.insert(segments.end(), d.begin(), d.end());
                    shapes.push_back({1, cx, cy, 2, 0});
                }
            }
        }
        Map map(segments);
        const auto bounds = map.bounds();
        std::vector<std::int64_t> owner_by_face(map.faceCount(), -2);

        for (long long x = bounds.xmin; x <= bounds.xmax; ++x) {
            for (long long y = bounds.ymin; y <= bounds.ymax; ++y) {
                std::int64_t owner = -1;
                bool on_boundary = false;
                for (std::size_t k = 0; k < shapes.size(); ++k) {
                    bool inside = false;
                    bool boundary = false;
                    if (shapes[k][0] == 0) {
                        inside = x > shapes[k][1] && x < shapes[k][3] &&
                                 y > shapes[k][2] && y < shapes[k][4];
                        boundary = x >= shapes[k][1] && x <= shapes[k][3] &&
                                   y >= shapes[k][2] && y <= shapes[k][4] &&
                                   !inside;
                    } else {
                        const long long l1 =
                            std::llabs(x - shapes[k][1]) +
                            std::llabs(y - shapes[k][2]);
                        inside = l1 < shapes[k][3];
                        boundary = l1 == shapes[k][3];
                    }
                    if (inside) {
                        owner = static_cast<std::int64_t>(k);
                    }
                    if (boundary) {
                        on_boundary = true;
                    }
                }
                if (on_boundary) {
                    continue;
                }
                const std::int64_t cell_id = map.locate(x, y);
                ASSERT_GE(cell_id, 0) << "gap at (" << x << "," << y << ")";
                const std::size_t face = map.cells()[cell_id].face;
                std::int64_t& slot = owner_by_face[face];
                if (slot == -2) {
                    slot = owner;
                } else {
                    EXPECT_EQ(slot, owner) << "face merges distinct regions";
                }
            }
        }
    }
}

} // namespace trapezoid::unittest
