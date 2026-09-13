# Panda Preserve: Fortune Sweep Solution

## 1. Problem reduction

Let the polygon vertices, which are also the possible transmitter sites, be

```text
S = {s0, s1, ..., s(n - 1)}.
```

For a point `x` in the park `P`, define

```text
f(x) = min(i) |x - si|.
```

The required receiver range is therefore

```text
R = max(x in P) f(x).
```

This is the radius of the largest empty circle whose center is constrained to
lie in the park and whose interior contains no polygon vertex. The circle may
cross the polygon boundary; only its center has to belong to the park.

The Voronoi cell of site `si` is

```text
Vi = {x : |x - si| <= |x - sj| for every j}.
```

For `x` in `Vi`, the nearest site is `si`, hence

```text
f(x)^2 = |x - si|^2.
```

The function on the right is convex. It has no strict local maximum in the
interior of a two-dimensional region: moving from an interior point in the
direction away from `si` increases the distance until the boundary is
reached. Thus a maximum over `P intersect Vi` lies on its boundary.

That boundary consists of pieces of:

- the park boundary, and
- Voronoi edges.

On a straight boundary piece parameterized by `p(t) = a + td`,

```text
|p(t) - si|^2
    = |d|^2 t^2 + 2 d dot (a - si) t + |a - si|^2,
```

which is a convex quadratic. Its maximum over any closed interval occurs at
an endpoint. The endpoints of all relevant pieces are either Voronoi vertices
inside the park or intersections of Voronoi edges with the park boundary.
Polygon vertices themselves need not be added separately: each one is a site,
so its nearest-site distance is zero.

Consequently, an optimum occurs at one of:

1. A Voronoi vertex inside the polygon.
2. An endpoint of a connected portion of a Voronoi edge inside the polygon,
   equivalently a Voronoi-edge/polygon-boundary intersection.

The implementation first constructs the Voronoi diagram with Fortune's
algorithm and then finds these candidates with sweep-line clipping.

## 2. Numeric representation

All geometric calculations use `long double` through the `Real` alias.

Two scale-dependent tolerances are used by the Fortune sweep:

- `length_epsilon_ = scale * 1e-11` for coordinate and distance comparisons.
- `area_epsilon_ = scale^2 * 1e-12` for orientation and determinant
  comparisons.

Here `scale` is at least `1` and at least the absolute value of every input
coordinate. The different dimensions matter: coordinate differences scale
linearly, while cross products scale quadratically.

Input vertices are distinct integer points, but Voronoi vertices and
intersections are generally non-integral.

The later polygon-clipping sweeps use the fixed `EPSILON = 1e-9`. Distances
are accumulated as squared distances and only the final answer is square
rooted. Besides avoiding repeated square roots, this keeps comparisons
monotone and slightly more stable.

## 3. Fortune Voronoi construction

Fortune's algorithm moves a horizontal sweep line, called the **directrix**,
from top to bottom. The processed sites lie above the directrix. Their
parabolas separate points already known to be closer to a processed site from
points not reached by the sweep. The lower envelope of those parabolas is the
**beach line**.

Each maximal parabola portion on the envelope is an arc. A breakpoint between
two neighboring arcs traces a Voronoi edge because every point on it is
equidistant from the two corresponding sites.

### 3.1 Events

The event queue processes events from larger `y` to smaller `y`:

- **Site event at `(sx, sy)`:** the directrix reaches a new site and its
  parabola enters the beach line.
- **Circle event at `(cx, cy - r)`:** three neighboring arcs meet at the
  circumcenter `(cx, cy)` and the middle arc disappears when the directrix
  reaches the bottom of their circumcircle.

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

This lazy invalidation is important because `std::priority_queue` cannot erase
an arbitrary element efficiently. The test

```cpp
middle->circle == event
```

also ensures that an old event cannot remove an arc after a newer event has
been scheduled for it.

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

#### Computing a breakpoint

For a focus `(px, py)` and directrix `y = l`, its parabola satisfies

```text
(x - px)^2 + (y - py)^2 = (y - l)^2.
```

Solving for `y` gives

```text
y = ((x - px)^2 + py^2 - l^2) / (2(py - l)).
```

Equating this expression for two neighboring sites produces the quadratic
solved by `breakpointX()`. There may be two algebraic roots, but only one is
the breakpoint on the lower envelope. The implementation selects the lower
root when the left focus is higher and the upper root otherwise.

The quadratic formula is evaluated with

```text
q = -0.5 * (b + sign(b) * sqrt(discriminant))
x1 = q / a
x2 = c / q
```

to avoid catastrophic cancellation when `b` and the square root have similar
magnitudes. Equal focus heights and a focus lying directly on the directrix
are handled separately because the general formula degenerates there.

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

The original `split` object is removed from the AVL tree but remains safely
allocated in `arcs_`. Stable allocation is required because stale queue events
may still contain its pointer; lazy invalidation makes those pointers harmless.

After either insertion form, only triples touching the modified neighborhood
can have gained or lost a circle event. Rechecking those local arcs is enough;
all other triples are unchanged.

### 3.4 Circle event

For consecutive sites `a`, `b`, and `c`, a disappearing middle arc requires
the correct turn for a downward sweep. `scheduleCircle()` therefore rejects
non-clockwise and near-collinear triples. It computes the circumcenter from
cross products, then schedules

```text
event_y = center.y - radius.
```

Events above the current directrix are stale or geometrically impossible and
are rejected.

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

`addVertex()` increments `completedBreakpoints` even if the geometric point is
already stored. This distinction handles cocircular degeneracies: several
topological breakpoint completions can coincide at one geometric Voronoi
vertex.

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

The clipping box is exactly the polygon's axis-aligned bounding box. This is
sufficient because every point of the polygon lies in that box, so no
discarded part of a Voronoi edge can intersect the park.

For sites `p` and `q`, a direction vector along their perpendicular bisector
is

```text
d = (py - qy, qx - px).
```

For a ray with one finite Voronoi vertex, there are two possible signs of
`d`. The `opposite` site stored with the circle event identifies the side
occupied by the Delaunay triangle. The ray must point away from that site, so
the implementation flips `d` when

```text
dot(opposite - p, d) > 0.
```

`clipParameterRange()` is a slab intersection. It starts with the allowed
parameter interval (`[0, 1]` for a segment or `[0, +infinity)` for a ray) and
intersects it with the parameter intervals imposed by the box's `x` and `y`
slabs. An empty interval means the edge never reaches the box.

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

Suppose the part of one Voronoi edge inside the non-convex polygon is a union
of parameter intervals:

```text
[a1, b1] union [a2, b2] union ... union [ak, bk].
```

Convexity implies that the maximum over their union occurs at the globally
smallest or globally largest retained parameter, namely `a1` or `bk` after
sorting the intervals. It is therefore sufficient to find:

- The first polygon crossing from the left endpoint.
- The first polygon crossing from the right endpoint.
- Any genuine Voronoi endpoint that lies inside the polygon.

There is no need to enumerate intermediate crossings.

If the bounded segment starts inside the polygon, its endpoint is a genuine
Voronoi vertex and is covered by the point-location sweep. Otherwise, the
first retained point is precisely the first boundary crossing. The symmetric
argument applies at the other end.

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

For the counterclockwise polygon order used by the problem, an edge directed
left-to-right has polygon interior above it and is a `LowerBoundary`. An edge
whose original direction is right-to-left is normalized by swapping its
endpoints; its interior lies below it and it receives the `UpperBoundary`
label. This label is later enough to answer point-location queries without
counting all ray crossings.

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

This avoids storing a mutable global sweep coordinate inside the comparator,
which would violate the ordering requirements of `std::set`. For segment
pairs that cannot cross while simultaneously active, their order at the ends
of the overlap determines a consistent ordering throughout the overlap.
Polygon/Voronoi pairs are removed at their first crossing, before their order
would reverse.

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

At insertion, only the predecessor and successor can be the first segment
crossed. At removal, those two neighbors become adjacent and are checked
together. Thus each set update schedules only `O(1)` intersection tests.

A crossing event is represented as removal of the Voronoi segment rather than
an order swap because this pass needs only that segment's first crossing.
Stale crossing or endpoint-removal events are harmless: `active.find()` fails
after the segment has already been removed.

Reflection for the second pass maps `x` to `-x`, then swaps each normalized
segment's endpoints. The same left-to-right implementation consequently sees
the original geometry from right to left.

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

Why does the first boundary above the query determine containment? On a
vertical line through a point not on the boundary, polygon crossings alternate
between entering and leaving the polygon. Above the highest interior interval
there is no polygon. Moving downward across a boundary labeled
`UpperBoundary` enters an interior interval, while crossing a `LowerBoundary`
leaves one. Therefore a query lies inside exactly when the first edge above it
has the `UpperBoundary` label.

Boundary points are valid centers as well. The event ordering and epsilon
tie-breaking retain them consistently; evaluating their distance again does
not affect the maximum.

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

The argument can be separated into four lemmas.

### Lemma 1: Fortune's sweep constructs the Voronoi edges

A site event inserts exactly the beach-line arc belonging to the new site and
creates the two breakpoints between it and the split arc. A valid circle event
occurs exactly when a middle arc shrinks to zero; it terminates its two
incident breakpoints at the three sites' circumcenter and starts the
breakpoint between the newly adjacent outer arcs. Lazy invalidation prevents
events for obsolete triples from changing the beach line. Hence all and only
Voronoi adjacencies and vertices are recorded.

### Lemma 2: The candidate set contains an optimum

Inside a Voronoi cell, the objective squared is distance squared to one fixed
site. This convex function has a maximum on the boundary of the cell portion
inside the park. Its restriction to every straight boundary piece is a convex
quadratic, whose maximum is at an endpoint. Such endpoints are Voronoi
vertices inside the park or Voronoi-edge/park-boundary intersections.

### Lemma 3: The sweeps examine every necessary candidate

The point-location sweep includes every genuine Voronoi vertex in the park.
For a Voronoi edge, convexity means only the two extreme points of all
inside-park portions can maximize the objective. If an extreme is a Voronoi
endpoint, point location includes it; otherwise it is the first polygon
crossing when scanning from that end. The forward and reflected sweeps find
exactly these first crossings.

### Theorem

By Lemma 2, some optimum belongs to the candidate set. By Lemma 3, the
algorithm evaluates every candidate needed to include that optimum. Each
distance is measured to a generating site of the containing Voronoi edge or
vertex, so it equals the nearest-site distance there. Therefore the largest
recorded distance is exactly the required receiver radius.

## 7. Complexity

Let `n` be the number of polygon vertices.

- Fortune processes `n` site events and `O(n)` valid circle events. Each queue
  or AVL operation costs `O(log n)`, for `O(n log n)` time.
- A planar Voronoi diagram has `O(n)` edges and vertices. Bounding-box
  clipping takes constant work per edge.
- Each clipping sweep has `O(n)` segment endpoint/removal events. Every event
  performs `O(1)` set or queue operations, each costing `O(log n)`, for
  `O(n log n)` time per sweep.
- The point-location sweep sorts `O(n)` events and performs `O(n)` balanced
  tree operations, also `O(n log n)`.
- The event queues, beach line, Voronoi records, and active sets all use
  `O(n)` memory.

Overall:

```text
Time:  O(n log n)
Space: O(n)
```
