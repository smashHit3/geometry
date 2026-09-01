# P2093 JZPFAR

- Problem: <https://www.luogu.com.cn/problem/P2093>
- Algorithm: build a static two-dimensional K-D tree with bounding rectangles.  For each
  query, visit subtrees in decreasing upper-bound distance order and retain the best `k`
  candidates in a heap. Equal distances use smaller input ids as farther.
- Time: build `O(n log n)`; query is `O(log n + k)` on typical balanced spatial data
  (with K-D-tree worst case `O(n)`).
- Space: `O(n + k)`

Build from the repository root:

```sh
cmake -S . -B build
cmake --build build --target luogu_p2093_jzpfar
```

Run with standard input:

```sh
./build/problems/luogu_p2093_jzpfar < input.txt
```
