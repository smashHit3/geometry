# Kattis Point in Polygon

- Problem: <https://open.kattis.com/problems/pointinpolygon>
- Algorithm: persistent sweep-line trapezoidal map
- Preprocessing: `O(n log n)` time and `O(n log n)` space
- Query: `O(log n)` time

## Approach

Draw a vertical line through every distinct polygon-vertex x-coordinate. Between
two consecutive lines, no polygon edges start, end, or cross, so the vertical
order of all intersecting non-vertical edges is fixed. Those ordered edges
partition the slab into trapezoids.

Sweep the vertex x-coordinates from left to right while maintaining the
intersecting non-vertical edges in height order. At an event x-coordinate,
edges ending there are removed in their left-side order and edges starting
there are inserted in their right-side order.

The active set is a persistent AVL tree. Each insertion or deletion creates
only `O(log n)` new nodes, and the root after each event represents the next
slab. Thus all slab search structures are constructed in `O(n log n)` total
time instead of copying every active edge into every slab.

A query binary-searches its slab and walks that slab's AVL root, accumulating
the number of edges below the point from subtree sizes. Starting below the
polygon, each crossed boundary alternates between outside and inside, so the
crossing parity gives `in` or `out`.

All comparisons use exact `__int128` cross multiplication. Polygon vertices
and vertical edges are indexed separately, allowing boundary queries to return
`on` without a linear scan.

The persistent versions share unchanged subtrees, resulting in `O(n log n)`
worst-case storage and deterministic bounds.

## Build and run

```sh
cmake -S . -B build
cmake --build build --target kattis_pointinpolygon_trapezoidal_map
./build/problems/kattis_pointinpolygon_trapezoidal_map < input.txt
```
