# Panda Preserve: Fortune Sweep Solution

## 1. Problem reduction

For a point `x` in the park, define

```text
f(x) = min distance(x, polygon vertex)
```

The required receiver range is

```text
max f(x), over all x inside or on the polygon.
```

The nearest polygon vertex is constant within each Voronoi cell. Inside one
cell, `f(x)` is the distance to that cell's site. Since squared distance is a
convex function, its maximum over a clipped Voronoi edge occurs at an extreme
point of the part lying inside the park.

Consequently, an optimum occurs at one of:

1. A Voronoi vertex inside the polygon.
2. An extreme intersection between a Voronoi edge and the polygon boundary.

The implementation first constructs the Voronoi diagram with Fortune's
algorithm and then finds these candidates with sweep-line clipping.

## 2. Numeric representation

All geometric calculations use `long double` through the `Real` alias.

Two scale-dependent tolerances are used by the Fortune sweep:

- `length_epsilon_` for coordinate and distance comparisons.
- `area_epsilon_` for orientation and determinant comparisons.

Input vertices are distinct integer points, but Voronoi vertices and
intersections are generally non-integral.

## 3. Fortune Voronoi construction

### 3.1 Events

The event queue processes events from larger `y` to smaller `y`:

- **Site event:** a new parabola enters the beach line.
- **Circle event:** an existing beach-line arc disappears.

At equal heights, site events are processed before circle events. Site events
are ordered from left to right.

Each circle event stores:

- The event position.
- The disappearing arc.
- A validity flag.

The corresponding code keeps both event kinds in one priority queue:

```cpp
enum class EventKind { Site, Circle };

struct Event {
    EventKind kind;
    Real y;
    Real x;
    std::size_t site = 0;
    Arc* arc = nullptr;
    bool valid = true;
    std::uint64_t sequence = 0;
};
```

`site` is used only by a site event, while `arc` identifies the arc removed by
a circle event. `sequence` gives deterministic ordering when all geometric
keys are equal.

When an arc's neighborhood changes, its pending circle event is invalidated
instead of being removed from the priority queue. Invalid events are ignored
when popped.

### 3.2 Beach line

The beach line is represented in two simultaneous forms:

- A doubly linked list using `previous` and `next`, which provides neighboring
  arcs in geometric order.
- An AVL tree using `left`, `right`, `parent`, and `height`, which locates the
  arc above a new site in logarithmic time.

`findArcAbove()` walks the AVL tree. At each node it evaluates the breakpoint
between the current arc and its next arc at the current directrix.

AVL maintenance consists of:

- `insertBetween()`
- `removeFromTree()`
- `rotateLeft()`
- `rotateRight()`
- `rebalanceFrom()`

The same `Arc` object participates in both representations:

```cpp
struct Arc {
    std::size_t site;

    // Beach-line order.
    Arc* previous = nullptr;
    Arc* next = nullptr;

    // AVL search structure.
    Arc* left = nullptr;
    Arc* right = nullptr;
    Arc* parent = nullptr;
    int height = 1;

    EdgeRecord* edgeToNext = nullptr;
    Event* circle = nullptr;
};
```

The search itself is a normal binary-tree descent, except that the key is a
moving parabola breakpoint:

```cpp
Arc* findArcAbove(Real x, Real directrix) const
{
    Arc* current = root_;
    Arc* candidate = nullptr;
    while (current != nullptr) {
        if (current->next == nullptr ||
            x <= breakpointX(current->site, current->next->site,
                             directrix)) {
            candidate = current;
            current = current->left;
        } else {
            current = current->right;
        }
    }
    return candidate == nullptr ? last_ : candidate;
}
```

At a fixed directrix the breakpoint order agrees with the linked-list order,
so the AVL tree can search by this dynamically evaluated key without rebuilding
the tree after every event.

### 3.3 Site event

For a normal site event:

1. Locate the beach-line arc above the new site.
2. Invalidate that arc's pending circle event.
3. Replace the old arc with three arcs:

   ```text
   old-site, new-site, old-site
   ```

4. Create the Voronoi edge between the old and new sites.
5. Update both the linked list and AVL tree.
6. Test the two newly formed triples for circle events.

Sites with the same `y` coordinate need special handling. Their parabolas are
degenerate at the current directrix, so splitting an arc into three would
create a zero-width duplicate arc. Instead, the new arc is inserted directly
after the existing equal-height arc:

```text
split -> inserted -> next
```

The edge ownership invariant is:

```text
arc->edgeToNext joins arc->site and arc->next->site
```

Therefore:

- `split->edgeToNext` becomes the `split/inserted` edge.
- `inserted->edgeToNext` becomes the `inserted/next` edge.

The equal-height branch implements this invariant directly:

```cpp
Arc* inserted = makeArc(event->site);
EdgeRecord* edge = edgeFor(split->site, event->site);

inserted->previous = split;
inserted->next = next;
inserted->edgeToNext =
    next == nullptr ? nullptr : edgeFor(inserted->site, next->site);

split->next = inserted;
split->edgeToNext = edge;
if (next != nullptr) {
    next->previous = inserted;
} else {
    last_ = inserted;
}
insertBetween(split, next, inserted);
```

For different heights, the old arc is replaced by three nodes. Both
breakpoints initially trace the same bisector:

```cpp
Arc* left = makeArc(split->site);
Arc* middle = makeArc(event->site);
Arc* right = makeArc(split->site);

left->next = middle;
left->edgeToNext = split_edge;
middle->previous = left;
middle->next = right;
middle->edgeToNext = split_edge;
right->previous = middle;
right->next = next;
right->edgeToNext = split->edgeToNext;
```

### 3.4 Circle event

When three consecutive arcs define a valid empty circle:

1. Compute their circumcenter.
2. Schedule the bottom of the circle as the event position.
3. When the event is processed, remove the middle arc.
4. Add the circumcenter to the two ending Voronoi edges.
5. Start the new Voronoi edge between the remaining neighboring sites.
6. Reschedule circle events for the changed neighbors.

Near-immediate events are accepted even when they lie extremely close to the
current directrix. Rejecting them with a large epsilon leaves expired arcs in
the beach line and makes later breakpoint searches invalid.

The circle-event update is the topological inverse of a site split:

```cpp
Arc* first = middle->previous;
Arc* third = middle->next;

addVertex(first->edgeToNext, center, third->site);
addVertex(middle->edgeToNext, center, first->site);
EdgeRecord* new_edge = edgeFor(first->site, third->site);
addVertex(new_edge, center, middle->site);

first->next = third;
first->edgeToNext = new_edge;
third->previous = first;
removeFromTree(middle);
```

The two old breakpoints terminate at `center`; the newly adjacent outer arcs
start a new Voronoi edge at the same point.

### 3.5 Edge completion

`EdgeRecord` stores:

- The two generating sites.
- Finite Voronoi vertices discovered by circle events.
- The number of completed beach-line breakpoints.

`clipEdge()` converts each record into a segment inside the polygon's bounding
box:

- Two completed breakpoints produce a finite segment.
- One completed breakpoint produces an unbounded ray.
- Zero completed breakpoints are discarded because, for this non-degenerate
  polygon input, they represent transient beach-line adjacencies rather than
  complete Voronoi edges.

The resulting `VoronoiEdge` retains its generating site indices, its bounded
segment, and its real circle-event vertices.

The distinction between a segment and a ray appears explicitly in
`clipEdge()`:

```cpp
if (edge.completedBreakpoints >= 2U) {
    // Two finite Voronoi vertices: clip the segment between them.
    clipParameterRange(edge.vertices[0].point,
                       edge.vertices[1].point - edge.vertices[0].point,
                       0.0, 1.0, start, end);
} else if (edge.completedBreakpoints == 1U) {
    // One finite vertex: orient and clip the unbounded ray.
    clipParameterRange(vertex.point, ray_direction, 0.0,
                       std::numeric_limits<Real>::infinity(),
                       start, end);
} else {
    return std::nullopt;
}
```

## 4. `O(n log n)` polygon clipping

Enumerating every polygon/Voronoi intersection can take quadratic time. The
implementation instead uses a sweep inspired by active-edge polygon clippers.

### 4.1 Why only extreme intersections matter

Parameterize a Voronoi edge as

```text
p(t) = p0 + t * direction.
```

For either generating site `s`,

```text
distance(p(t), s)^2
```

is a convex quadratic in `t`. Even if a non-convex polygon intersects the edge
in many disjoint intervals, the maximum over all retained intervals occurs at
the globally smallest or largest retained parameter.

It is therefore sufficient to find:

- The first polygon crossing from the left endpoint.
- The first polygon crossing from the right endpoint.
- Any genuine Voronoi endpoint that lies inside the polygon.

There is no need to enumerate intermediate crossings.

### 4.2 Preparing sweep segments

`makeSweepSegments()` converts polygon edges and Voronoi edges into
`SweepSegment` objects.

A fixed rotation gives the segments a generic sweep direction:

```cpp
Point rotateForClipping(const Point& point)
{
    constexpr Real cosine = 0.9923888851137123702L;
    constexpr Real sine = 0.1231434151823108666L;
    return {
        point.x() * cosine - point.y() * sine,
        point.x() * sine + point.y() * cosine,
    };
}
```

This is the standard two-dimensional rotation

```text
x' = x cos(theta) - y sin(theta)
y' = x sin(theta) + y cos(theta)
```

The angle is fixed and has no geometric significance. It avoids vertical input
segments, which would otherwise need special cases in `yAt()` and the active
segment comparator. Because the transformation is a rigid rotation, it
preserves:

- Distances, so the final answer is unchanged.
- Orientation and incidence.
- Segment intersections.
- Whether a point lies inside the polygon.

After rotation, every segment is normalized so that `start.x() < end.x()`.

Polygon edges are classified using their original direction:

- `LowerBoundary`
- `UpperBoundary`

Voronoi segments also store one generating site so candidate distances can be
evaluated in constant time.

### 4.3 Active-edge ordering

`SweepSegmentLess` orders active segments by their vertical positions over
their common `x` interval. Polygon edges do not cross each other, and Voronoi
edges do not cross each other except at shared endpoints. The only relevant
order changes are polygon/Voronoi crossings.

Each active segment can evaluate its height at the current horizontal
position:

```cpp
Real yAt(Real x) const
{
    return start.y() +
           (end.y() - start.y()) *
               (x - start.x()) / (end.x() - start.x());
}
```

The comparator evaluates two segments at the beginning and end of their
overlapping x-range. The stable `id` tie-break handles shared endpoints.

### 4.4 Extreme-intersection sweeps

`findExtremeIntersections()` processes:

- Segment insertion events.
- Segment removal events.
- Polygon/Voronoi crossing events.

Only adjacent active segments can be the next pair to cross. When an adjacent
polygon/Voronoi pair intersects:

1. Evaluate the distance from the intersection to the Voronoi edge's site.
2. Schedule removal of that Voronoi segment at the crossing.

Removing the Voronoi segment ensures that only its first crossing in the
current direction is processed. The sweep runs twice:

1. Original coordinates, finding one extreme.
2. Reflected `x` coordinates, finding the opposite extreme.

Events whose `x` coordinates differ by at most `EPSILON` are treated as
simultaneous. Removal events run before insertion and crossing events to keep
the active ordering valid at shared endpoints.

The key operation is `schedule_crossing`:

```cpp
const auto schedule_crossing =
    [&](const SweepSegment& first, const SweepSegment& second, Real) {
        if ((first.kind == SweepSegmentKind::Voronoi) ==
            (second.kind == SweepSegmentKind::Voronoi)) {
            return;
        }
        const auto intersection = segmentIntersection(first, second);
        if (!intersection.has_value()) {
            return;
        }

        const SweepSegment& voronoi =
            first.kind == SweepSegmentKind::Voronoi ? first : second;
        answer_squared = std::max(
            answer_squared,
            squaredDistance(*intersection, sites[voronoi.site]));

        events.push({intersection->x(), SweepEventKind::Crossing,
                     voronoi.id, sequence++});
    };
```

A crossing event names the Voronoi segment. When processed, that segment is
removed from the active set, so later crossings on the same directed scan are
never generated.

### 4.5 Voronoi-vertex point location

`findInteriorVoronoiVertices()` performs a third offline sweep:

1. Insert and remove polygon edges at their endpoint events.
2. Query every genuine circle-event vertex.
3. Find the first active polygon edge above the query.
4. Use its `LowerBoundary`/`UpperBoundary` classification to determine whether
   the query lies inside the polygon.

Artificial endpoints introduced when unbounded Voronoi rays are clipped to
the bounding box are deliberately excluded from these queries.

For a query point, `lower_bound` finds the first boundary edge above it:

```cpp
const SweepSegment query{
    segments.size() + event.sequence,
    point, point, SweepSegmentKind::Voronoi, site};
const auto above = active.lower_bound(query);

if (above != active.end() &&
    above->kind == SweepSegmentKind::UpperBoundary) {
    answer_squared = std::max(
        answer_squared, squaredDistance(point, sites[site]));
}
```

The polygon orientation and edge direction determine whether the region below
that edge is interior. This replaces a separate `O(n)` point-in-polygon test
for every Voronoi vertex.

## 5. Final flow

```text
read polygon vertices
        |
        v
run Fortune site/circle sweep
        |
        v
clip Voronoi edges to polygon bounding box
        |
        v
rotate polygon and Voronoi geometry
        |
        +--> point-location sweep for real Voronoi vertices
        |
        +--> left-to-right first-crossing sweep
        |
        +--> right-to-left first-crossing sweep
        |
        v
take the largest squared candidate distance
        |
        v
print its square root
```

The coordinating functions are intentionally small. `clipPolygon()` combines
the three independent candidate searches:

```cpp
Real answer_squared = findInteriorVoronoiVertices(
    segments, rotated_sites, vertices, polygon.size());
answer_squared = std::max(
    answer_squared,
    findExtremeIntersections(segments, rotated_sites, false));
answer_squared = std::max(
    answer_squared,
    findExtremeIntersections(segments, rotated_sites, true));
```

`main()` then shows the complete high-level solution:

```cpp
int main()
{
    int vertex_count;
    std::cin >> vertex_count;

    std::vector<Point> polygon(vertex_count);
    for (Point& point : polygon) {
        std::cin >> point.x() >> point.y();
    }

    const std::vector<VoronoiEdge> voronoi_edges =
        buildFortuneEdges(polygon);
    const Real maximum_squared_radius =
        clipPolygon(voronoi_edges, polygon);

    std::cout << std::setprecision(12)
              << std::sqrt(maximum_squared_radius) << '\n';
}
```

## 6. Correctness outline

1. Fortune's event sweep records every adjacency change of the beach line, so
   its completed breakpoints form the Voronoi diagram.
2. A site event creates exactly the two breakpoints surrounding the new arc.
3. A valid circle event removes exactly one disappearing arc and joins its two
   neighbors at the correct circumcenter.
4. On a Voronoi edge, squared distance to either generating site is a convex
   quadratic in the edge parameter.
5. Therefore only the first retained point from each direction can maximize
   the distance on that edge.
6. The two directional sweeps find those points, while the point-location
   sweep includes all genuine Voronoi vertices inside the park.
7. These are all possible maxima, so the largest recorded distance is the
   required radius.

## 7. Complexity

Let `n` be the number of polygon vertices.

- Fortune construction: `O(n log n)`.
- Number of Voronoi edges and vertices: `O(n)`.
- Sweep-line clipping and point location: `O(n log n)`.
- Memory: `O(n)`.

Overall:

```text
Time:  O(n log n)
Space: O(n)
```
