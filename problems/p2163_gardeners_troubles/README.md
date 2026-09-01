# P2163 [SHOI2007] Gardener's Troubles

- Problem: <https://www.luogu.com.cn/problem/P2163>
- Algorithm: reduce every rectangle count to four inclusive dominance counts using
  inclusion-exclusion. Sort those events by x-coordinate and sweep the points while a
  Fenwick tree counts y-coordinates.
- Time: `O((n + m) log(n + m))`
- Space: `O(n + m)`

Build from the repository root:

```sh
cmake -S . -B build
cmake --build build --target luogu_p2163_gardeners_troubles
```

Run with standard input:

```sh
./build/problems/luogu_p2163_gardeners_troubles < input.txt
```
