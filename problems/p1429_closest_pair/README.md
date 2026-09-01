# P1429 Planar Closest Pair

- Problem: <https://www.luogu.com.cn/problem/P1429>
This directory contains four independent C++17 implementations:

| File | Approach | Expected / typical time | Space |
| --- | --- | --- | --- |
| `main.cpp` | Divide and conquer: split by x, merge by y, and inspect the cross-half strip. | `O(n log n)` worst case | `O(n)` |
| `bvh.cpp` | Static median-split bounding-volume hierarchy over point bounding boxes; query each point with nearest-neighbor branch-and-bound. | `O(n log n)` build and typically near `O(n log n)` total querying | `O(n)` |
| `kd_tree.cpp` | Static median-split K-D tree with a bounding rectangle at every subtree; query each point with nearest-neighbor branch-and-bound. | `O(n log n)` build and typically near `O(n log n)` total querying | `O(n)` |
| `r_tree.cpp` | Static STR-style bulk-loaded R-tree with MBRs; query each point with best-first nearest-neighbor branch-and-bound. | `O(n log n)` bulk loading and typically near `O(n log n)` total querying | `O(n)` |

The three nearest-neighbor tree variants exclude the queried point by index, so coincident
points are handled correctly. Their pruning is exact, but spatial tree nearest-neighbor
queries can have adversarial `O(n)` cost each (and therefore `O(n^2)` total); use
`main.cpp` for P1429's guaranteed `O(n log n)` bound.

Build from the repository root:

```sh
cmake -S . -B build
cmake --build build --target \
  luogu_p1429_closest_pair \
  luogu_p1429_closest_pair_bvh \
  luogu_p1429_closest_pair_kd_tree \
  luogu_p1429_closest_pair_r_tree
```

Run with standard input:

```sh
./build/problems/luogu_p1429_closest_pair < input.txt
```

Replace the executable name with any target from the table to run another implementation.
