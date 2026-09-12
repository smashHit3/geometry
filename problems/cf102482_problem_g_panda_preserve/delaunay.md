# Panda Preserve: Divide-and-Conquer Delaunay Solution

## 1. Problem reduction

For every point `x` inside the polygon, only its nearest polygon vertex
matters. The function

```text
f(x) = min distance(x, polygon vertex)
```

is represented by the Voronoi diagram of the polygon vertices.

Within one Voronoi cell, the nearest site is fixed. The maximum distance in
the polygon therefore occurs at:

1. A Voronoi-cell vertex inside the polygon, or
2. An intersection between a Voronoi-cell edge and the polygon boundary.

This implementation obtains Voronoi cells from the dual Delaunay
triangulation.

## 2. Exact Delaunay predicates

Input coordinates are integers. The triangulation uses `__int128_t` for:

- Orientation tests.
- In-circle determinants.

This avoids topology errors caused by floating-point predicates. Floating
point is introduced only after the Delaunay adjacency graph has been built.

The orientation predicate determines which side of a directed edge contains a
point. The in-circle predicate decides whether an edge violates the Delaunay
condition during merging.

The implementation keeps these decisions exact:

```cpp
using WideInteger = __int128_t;

bool leftOf(int point, QuadEdge* edge) const
{
    return cross(points_[point], points_[edge->origin],
                 points_[edge->destination()]) > 0;
}

bool inCircle(int first, int second, int third, int fourth) const
{
    return incircleDeterminant(
               points_[first], points_[second],
               points_[third], points_[fourth]) > 0;
}
```

The algorithm never makes a triangulation decision from an approximate
floating-point value.

## 3. Quad-edge representation

Each geometric edge is represented by four linked `QuadEdge` records:

- Two directed primal edges.
- Two directed dual edges.

Important operations are:

- `reversed()`: the same primal edge in the opposite direction.
- `destination()`: the destination vertex.
- `leftNext()`: the next edge around the face on the left.
- `originPrevious()`: the previous edge around the origin.
- `splice()`: joins or separates edge rings.
- `connect()`: creates an edge between two existing subdivisions.
- `removeEdge()`: removes an invalid Delaunay edge.

Edges are allocated from owned storage and marked as removed rather than
freed. This keeps all topology pointers stable throughout recursive merging.

The compact record derives all navigation from `rotated` and `next`:

```cpp
struct QuadEdge {
    int origin = -1;
    QuadEdge* rotated = nullptr;
    QuadEdge* next = nullptr;
    bool removed = false;

    QuadEdge* reversed() const { return rotated->rotated; }
    int destination() const { return reversed()->origin; }
    QuadEdge* leftNext() const
    {
        return rotated->reversed()->next->rotated;
    }
    QuadEdge* originPrevious() const
    {
        return rotated->next->rotated;
    }
};
```

`splice()` is the fundamental topology operation:

```cpp
static void splice(QuadEdge* first, QuadEdge* second)
{
    std::swap(first->next->rotated->next,
              second->next->rotated->next);
    std::swap(first->next, second->next);
}
```

Swapping the primal rings and the corresponding dual rings keeps the
subdivision consistent. `connect()` uses two splices to insert a new diagonal:

```cpp
QuadEdge* connect(QuadEdge* first, QuadEdge* second)
{
    QuadEdge* edge =
        makeEdge(first->destination(), second->origin);
    splice(edge, first->leftNext());
    splice(edge->reversed(), second);
    return edge;
}
```

## 4. Divide-and-conquer triangulation

### 4.1 Sorting

All sites are sorted lexicographically by `(x, y)`.

### 4.2 Base cases

For two points:

1. Create one edge.
2. Return its two directions as the outer hull references.

For three points:

1. Connect the first pair and second pair.
2. Splice them around the middle vertex.
3. If the points are non-collinear, add the closing edge.
4. Return the left and right outer hull edges according to orientation.

The collinear case remains an open chain.

The two-point base case illustrates the meaning of the returned pair:

```cpp
if (right - left == 1) {
    QuadEdge* edge = makeEdge(order[left], order[right]);
    return {edge, edge->reversed()};
}
```

The first value is the left outer reference and the second is the right outer
reference. These are the handles needed by the parent merge.

### 4.3 Recursive split

`buildRange()` divides the sorted range into two halves and recursively builds
their Delaunay triangulations.

Each recursive result returns:

- Its left outer hull edge.
- Its right outer hull edge.

```cpp
const int middle = (left + right) / 2;
auto [left_outer, left_inner] =
    buildRange(order, left, middle);
auto [right_inner, right_outer] =
    buildRange(order, middle + 1, right);
```

### 4.4 Lower common tangent

The merge starts by walking the inner hull edges until they form the lower
common tangent of the two triangulations.

That tangent is connected with a new base edge.

```cpp
while (true) {
    if (leftOf(right_inner->origin, left_inner)) {
        left_inner = left_inner->leftNext();
    } else if (rightOf(left_inner->origin, right_inner)) {
        right_inner = right_inner->reversed()->next;
    } else {
        break;
    }
}

QuadEdge* base = connect(right_inner->reversed(), left_inner);
```

Each iteration moves one inner hull handle toward the true tangent. When
neither endpoint lies on the forbidden side, `base` can safely join the two
triangulations.

### 4.5 Delaunay merge

The base edge advances upward:

1. Select the next candidate edge from the left triangulation.
2. Remove left edges whose opposite point lies inside the relevant
   circumcircle.
3. Perform the symmetric operation on the right.
4. If neither candidate is valid, merging is complete.
5. Otherwise, use an in-circle test to choose the left or right candidate.
6. Connect the selected candidate to form the next base edge.

Every removed edge violates the local Delaunay condition. When no valid
candidate remains, the two halves form one Delaunay triangulation.

The left-side deletion loop is representative:

```cpp
QuadEdge* left_candidate = base->reversed()->next;
if (valid(left_candidate)) {
    while (inCircle(
        base->destination(), base->origin,
        left_candidate->destination(),
        left_candidate->next->destination())) {
        QuadEdge* next = left_candidate->next;
        removeEdge(left_candidate);
        left_candidate = next;
    }
}
```

If the next point lies inside the candidate circumcircle, the candidate edge
cannot belong to a Delaunay triangulation and is deleted. The right side is
cleaned symmetrically. The remaining candidates are compared to choose the
next base:

```cpp
if (!left_valid ||
    (right_valid &&
     inCircle(left_candidate->destination(),
              left_candidate->origin,
              right_candidate->origin,
              right_candidate->destination()))) {
    base = connect(right_candidate, base->reversed());
} else {
    base = connect(base->reversed(),
                   left_candidate->reversed());
}
```

### 4.6 Extracting neighbors

After recursion completes, all non-removed primal edges are scanned.
For every edge `(u, v)`:

```text
v is added to neighbors[u]
u is added to neighbors[v]
```

Sets remove duplicate adjacency entries.

```cpp
for (QuadEdge* edge : primal_edges_) {
    if (edge->removed) {
        continue;
    }
    const int first = edge->origin;
    const int second = edge->destination();
    neighbor_sets[first].insert(second);
    neighbor_sets[second].insert(first);
}
```

## 5. Constructing Voronoi cells

The Voronoi cell of site `p` is the intersection of half-planes

```text
distance(x, p) <= distance(x, q)
```

for every Delaunay neighbor `q`.

After expanding and cancelling quadratic terms:

```text
2 * (q - p) dot x <= |q|^2 - |p|^2
```

`buildVoronoiCell()` starts with a padded bounding rectangle and clips it
against each neighbor half-plane.

`clipPolygon()` performs one Sutherland-Hodgman clipping step:

1. Evaluate both endpoints of every polygon edge against the half-plane.
2. Preserve an inside endpoint.
3. If the edge crosses the boundary line, append the intersection.

Only Delaunay neighbors are required because they define all non-redundant
boundaries of a Voronoi cell.

The equation is translated directly into code:

```cpp
const Point& site = sites[site_index];
for (const int neighbor_index : neighbors) {
    const Point& neighbor = sites[neighbor_index];
    const Point normal = neighbor - site;
    const Real offset =
        (dot(neighbor, neighbor) - dot(site, site)) / 2;
    cell = clipPolygon(cell, normal, offset);
}
```

For one half-plane, Sutherland-Hodgman clipping retains inside vertices and
adds a line intersection whenever inside/outside state changes:

```cpp
const bool current_inside = current_value <= EPSILON;
const bool next_inside = next_value <= EPSILON;

if (current_inside) {
    clipped.push_back(current);
}
if (current_inside != next_inside) {
    const Real ratio =
        current_value / (current_value - next_value);
    clipped.push_back(current + (next - current) * ratio);
}
```

## 6. Intersecting cells with the park

For each site's completed Voronoi cell:

1. Test every cell vertex with `pointInPolygon()`.
2. Intersect every cell edge with every park boundary edge.
3. For every valid candidate, compute its squared distance to the cell's
   site.
4. Keep the largest squared distance globally.

Using squared distances avoids unnecessary square roots. One square root is
taken when printing the final answer.

The outer loop associates each cell with its nearest site:

```cpp
for (int site_index = 0; site_index < vertex_count; ++site_index) {
    const std::vector<Point> cell = buildVoronoiCell(
        site_index, polygon, neighbors[site_index],
        minimum_x - padding, minimum_y - padding,
        maximum_x + padding, maximum_y + padding);

    for (const Point& vertex : cell) {
        if (pointInPolygon(vertex, polygon)) {
            answer_squared = std::max(
                answer_squared,
                squaredDistance(vertex, polygon[site_index]));
        }
    }

    // Every cell edge is also checked against every park edge.
}
```

At a point belonging to this cell, `polygon[site_index]` is a nearest input
vertex, so the measured distance is exactly the objective function.

## 7. Final flow

```text
read integer polygon vertices
        |
        v
sort sites by coordinate
        |
        v
divide into halves recursively
        |
        v
merge with lower tangents and in-circle tests
        |
        v
extract Delaunay neighbor graph
        |
        v
for each site:
    clip a bounding rectangle by neighbor bisectors
    test cell vertices inside the park
    intersect cell edges with the park boundary
        |
        v
take the largest squared candidate distance
        |
        v
print its square root
```

The central part of `main()` mirrors that flow:

```cpp
const std::vector<std::vector<int>> neighbors =
    delaunayNeighbors(integer_vertices);

for (int site_index = 0; site_index < vertex_count; ++site_index) {
    const std::vector<Point> cell = buildVoronoiCell(
        site_index, polygon, neighbors[site_index],
        minimum_x - padding, minimum_y - padding,
        maximum_x + padding, maximum_y + padding);

    // Check cell vertices inside the park.
    // Check intersections of cell edges and park edges.
}

std::cout << std::fixed << std::setprecision(10)
          << std::sqrt(answer_squared) << '\n';
```

## 8. Correctness outline

1. The divide-and-conquer merge preserves the Delaunay empty-circumcircle
   condition, so the resulting adjacency graph is a Delaunay triangulation.
2. The dual Delaunay neighbors define every boundary half-plane of each
   Voronoi cell.
3. Half-plane clipping therefore constructs the correct cell for every site.
4. Within a cell, squared distance to its site is convex.
5. A convex function over each clipped polygonal portion reaches its maximum
   at a boundary vertex.
6. The implementation examines all relevant Voronoi vertices and all
   Voronoi-cell/polygon intersections, so it finds the global optimum.

## 9. Complexity

Let `n` be the number of polygon vertices.

- Divide-and-conquer Delaunay triangulation: `O(n log n)`.
- Delaunay storage: `O(n)`.
- Voronoi-cell construction and polygon-boundary checks: `O(n^2)` overall in
  this implementation.
- Memory: `O(n)`, excluding temporary cell polygons.

Overall:

```text
Time:  O(n^2)
Space: O(n)
```
