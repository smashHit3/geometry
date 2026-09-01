# P4148 Simple Problem

- Problem: <https://www.luogu.com.cn/problem/P4148>
- Algorithm: maintain inserted points in a dynamically rebuilt two-dimensional K-D tree.
  Each node stores its bounding rectangle, subtree size, and subtree sum. A subtree is
  rebuilt when one child exceeds three quarters of its size; rectangle queries prune
  disjoint boxes and aggregate wholly contained boxes.
- Time: amortized `O(log² n)` insertion; typical `O(sqrt(n) + log n)` rectangle query
  (K-D-tree worst case is `O(n)`).
- Space: `O(n)`

Build from the repository root:

```sh
cmake -S . -B build
cmake --build build --target luogu_p4148_simple_problem
```

Run with standard input:

```sh
./build/problems/luogu_p4148_simple_problem < input.txt
```
