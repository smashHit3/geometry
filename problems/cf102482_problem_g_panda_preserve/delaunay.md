# Panda Preserve: Divide-and-Conquer Delaunay Solution

## 1. Problem reduction

Let the polygon be `P` and let its vertices be the sites

```text
S = {s0, s1, ..., s(n - 1)}.
```

For every point `x` inside the polygon, only its nearest site matters. Define

```text
f(x) = min(i) |x - si|.
```

The answer is

```text
R = max(x in P) f(x).
```

Equivalently, `R` is the radius of the largest empty circle centered in the
park, where empty means that no polygon vertex lies in its interior. The
circle itself does not have to remain inside the polygon.

The Voronoi cell of site `si` is

```text
Vi = {x : |x - si| <= |x - sj| for every j}.
```

Within `Vi`, the objective is simply

```text
f(x)^2 = |x - si|^2.
```

This is a convex quadratic. Over a polygonal region `P intersect Vi`, it
attains a maximum at a boundary vertex: an interior point can be moved away
from `si` until a boundary is reached, and on each straight boundary segment
the function is a one-dimensional convex quadratic whose maximum is at an
endpoint.

Therefore the only candidates are:

1. A Voronoi-cell vertex inside the polygon.
2. An intersection between a Voronoi-cell edge and the polygon boundary.
3. A polygon vertex inside a cell; these have value zero because every polygon
   vertex is itself a site, so they do not need to be tested explicitly.

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

For points `o`, `a`, and `b`, orientation is the signed doubled area

```text
orient(o, a, b) = (a - o) cross (b - o).
```

Its sign is positive for a counterclockwise turn, negative for a clockwise
turn, and zero for collinear points.

The in-circle test first translates all coordinates by the query point `d`.
For counterclockwise `a`, `b`, `c`, the determinant

```text
| ax  ay  ax^2 + ay^2 |
| bx  by  bx^2 + by^2 |
| cx  cy  cx^2 + cy^2 |
```

is positive exactly when `d` lies strictly inside their circumcircle. The
merge calls the predicate with an orientation consistent with this sign
convention.

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
floating-point value. `long double` is used only later, when bisector
intersections and Euclidean distances must be represented.

## 3. Quad-edge representation

Each geometric edge is represented by four linked `QuadEdge` records:

- `e`: one directed primal edge.
- `e->rotated`: the corresponding directed dual edge.
- `e->reversed()`: the primal edge in the opposite direction.
- `e->rotated->reversed()`: the opposite dual direction.

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

Every directed edge participates in an origin ring: repeatedly following
`next` walks counterclockwise through edges with the same origin. The dual
rings encode face adjacency. From only `rotated` and `next`, standard
quad-edge identities derive all other navigation.

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

The new edge runs from `first->destination()` to `second->origin`, joining the
left face of `first` to the face containing `second`.

Deletion is also expressed entirely through splices:

```cpp
splice(edge, edge->originPrevious());
splice(edge->reversed(), edge->reversed()->originPrevious());
```

This detaches both directed sides from their origin rings. The records remain
allocated, but `removed` excludes the edge when neighbors are extracted.

## 4. Divide-and-conquer triangulation

### 4.1 Sorting

All sites are sorted lexicographically by `(x, y)`. Recursive ranges are
therefore separated from left to right, which is the geometric precondition
for finding their common tangents during a merge.

`buildRange(left, right)` returns two directed hull handles:

- an edge incident to the leftmost site and exposed on the lower hull;
- an edge incident to the rightmost site and exposed on the lower hull.

The parent recursion uses these handles as starting points for its tangent
walk. The complete triangulation stays in the shared quad-edge structure, so
the function does not need to return all of its edges.

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

For a non-collinear triple, `connect(second, first)` closes the triangle. The
orientation determines which directed versions of the three edges are the
correct outer handles. No in-circle test is needed because three
non-collinear points have a unique triangulation.

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

The names `left_inner` and `right_inner` refer to the hull handles facing the
other half. `left_outer` and `right_outer` remain candidates for the outer
hull of the merged result.

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

After creating `base`, the code updates an outer handle if the tangent reaches
the corresponding extreme site:

```cpp
if (left_inner->origin == left_outer->origin) {
    left_outer = base->reversed();
}
if (right_inner->origin == right_outer->origin) {
    right_outer = base;
}
```

Those adjusted handles are returned after the upward merge finishes.

### 4.5 Delaunay merge

Starting from the lower tangent, the base edge advances upward and stitches
the two triangulations together:

1. Select the next candidate edge from the left triangulation.
2. Remove left edges whose opposite point lies inside the relevant
   circumcircle.
3. Perform the symmetric operation on the right.
4. If neither candidate is valid, merging is complete.
5. Otherwise, use an in-circle test to choose the left or right candidate.
6. Connect the selected candidate to form the next base edge.

Every removed edge violates the local Delaunay condition. When no valid
candidate remains, the two halves form one Delaunay triangulation.

For the current `base`, a candidate is geometrically usable only when its
destination lies to the right of `base`. Candidates on the other side cannot
form the next cross-edge without leaving the unmerged gap:

```cpp
const auto valid = [&](QuadEdge* edge) {
    return rightOf(edge->destination(), base);
};
```

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

If only one side has a valid candidate, that side is forced. If both are
valid, the in-circle predicate selects the diagonal satisfying the local
Delaunay condition. Each new connection raises `base`; when neither candidate
is valid, `base` has reached the upper common tangent and the merge is
complete.

Although a merge loop can delete edges, each edge is created once and removed
at most once. This amortization is what keeps all merge work at one recursion
level linear.

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

The algorithm only needs this undirected neighbor graph after triangulation.
It does not explicitly construct dual Voronoi vertices or retain face objects.

## 5. Constructing Voronoi cells

The Voronoi cell of site `p` is the intersection of half-planes

```text
distance(x, p) <= distance(x, q)
```

for every Delaunay neighbor `q`.

After expanding and cancelling quadratic terms:

```text
|x - p|^2 <= |x - q|^2

|x|^2 - 2p dot x + |p|^2
    <= |x|^2 - 2q dot x + |q|^2

(q - p) dot x <= (|q|^2 - |p|^2) / 2.
```

`buildVoronoiCell()` starts with a padded bounding rectangle and clips it
against each neighbor half-plane.

It is sufficient to use only Delaunay neighbors. Every Voronoi edge separating
the cells of `p` and `q` is dual to a Delaunay edge `(p, q)`, and all
non-neighbor inequalities are redundant once the neighbor inequalities have
been applied.

The initial rectangle is the polygon bounding box expanded in every direction
by

```text
padding = max(box width, box height) + 1.
```

Every point of the park is strictly inside this rectangle. The code only uses
the resulting cell inside the park or where its edges meet the park boundary,
so replacing an unbounded Voronoi cell by its intersection with this rectangle
cannot remove a relevant candidate.

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

The interpolation ratio follows from linearity. If

```text
v(t) = value(current + t(next - current)),
```

then `v(t) = current_value + t(next_value - current_value)`. Setting `v(t) = 0`
gives

```text
t = current_value / (current_value - next_value).
```

Because the cell starts convex and every clipping operation intersects it with
a half-plane, it remains a convex polygon throughout.

## 6. Intersecting cells with the park

For each site's completed Voronoi cell:

1. Test every cell vertex with `pointInPolygon()`.
2. Intersect every cell edge with every park boundary edge.
3. For every valid candidate, compute its squared distance to the cell's
   site.
4. Keep the largest squared distance globally.

Using squared distances avoids unnecessary square roots. One square root is
taken when printing the final answer.

`pointInPolygon()` uses the standard horizontal-ray parity test. A separate
`pointOnSegment()` check makes the park boundary inclusive, as required by
the maximization domain.

For two non-parallel segments,

```text
a + tu = b + sv,
```

cross products give

```text
t = cross(b - a, v) / cross(u, v)
s = cross(b - a, u) / cross(u, v).
```

The intersection is retained when both parameters lie in `[0, 1]`, with
`EPSILON` tolerance. Parallel or collinear pairs do not need a separate
overlap candidate: overlap endpoints are already cell or polygon vertices,
and the convex distance function reaches its maximum at one of those
endpoints.

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

Notice that the park itself may be non-convex. The implementation never clips
the Voronoi cell by the whole park using a convex-polygon algorithm. Instead,
it tests cell vertices for containment and independently intersects every
cell edge with every park edge. This enumerates the boundary vertices of every
possibly disconnected component of `Vi intersect P`.

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

The proof is organized into four lemmas.

### Lemma 1: The extracted graph contains all Delaunay neighbors

The base cases produce valid Delaunay subdivisions. During a merge, the lower
common tangent is a valid first cross-edge. The candidate-deletion loops remove
exactly edges that violate the empty-circumcircle condition, and the candidate
selection chooses the locally Delaunay next cross-edge. When the upper tangent
is reached, the merged subdivision is Delaunay. By induction over recursive
ranges, every surviving edge belongs to a Delaunay triangulation, and scanning
those edges extracts its complete neighbor graph.

### Lemma 2: Half-plane clipping constructs each Voronoi cell

A point belongs to site `p`'s cell exactly when it lies on `p`'s side of every
perpendicular bisector. Voronoi boundaries correspond exactly to Delaunay
neighbors, so inequalities from other sites are redundant. Starting from a
rectangle containing the whole park and successively intersecting all neighbor
half-planes therefore constructs the part of `p`'s Voronoi cell relevant to
the park.

### Lemma 3: The implementation evaluates all possible maxima

Inside a Voronoi cell, squared nearest-site distance equals squared distance to
that cell's site and is convex. Its maximum over each polygonal component of
the cell/park intersection occurs at a boundary vertex. Such a vertex is
either a Voronoi-cell vertex inside the park, which `pointInPolygon()` accepts,
or a cell-edge/park-edge intersection, which the nested segment loops find.
Polygon vertices contribute only zero.

### Theorem

By Lemma 2, every point of the park is considered in the cell of one of its
nearest sites. By Lemma 3, the algorithm evaluates a candidate attaining the
maximum in every cell/park intersection. Taking the maximum over all sites
therefore yields `max(x in P) min(i) |x - si|`, the required receiver radius.

## 9. Complexity

Let `n` be the number of polygon vertices.

- Sorting costs `O(n log n)`. At each recursion level, tangent walks,
  connections, and amortized deletions cost `O(n)`, so Delaunay construction
  costs `O(n log n)`.
- A planar Delaunay triangulation has `O(n)` edges, hence the sum of all site
  degrees is `O(n)`.
- For a site of degree `d`, its clipped cell has `O(d)` edges. Building it by
  repeated half-plane clipping costs `O(d^2)` in the straightforward
  implementation. Summed over all sites, this is at most `O(n^2)`.
- Testing all cell vertices with `pointInPolygon()` and all cell edges against
  all `n` park edges costs `O(nd)` for one site. Since the degree sum is
  `O(n)`, the total is `O(n^2)`.
- The triangulation, neighbor graph, and all cells considered one at a time
  require `O(n)` memory.

Overall:

```text
Time:  O(n^2)
Space: O(n)
```
