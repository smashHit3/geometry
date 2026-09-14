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
algorithm and then checks its vertices and its intersections with the polygon
boundary.

![A rectangular park whose four corner sites produce a central Voronoi vertex,
four boundary intersections, and an empty circle through every site](fortune-overview.svg)

The red corners are both polygon vertices and Voronoi sites. Purple dashed
segments separate the four shaded Voronoi cells. Gold points are
edge/boundary intersections, while the blue center is an interior Voronoi
vertex. In this exact example, the blue circle is centered at the answer and
passes through all four sites.

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

The later polygon clipping uses the fixed `EPSILON = 1e-9`. Distances are
accumulated as squared distances and only the final answer is square rooted.
Besides avoiding repeated square roots, this keeps comparisons monotone and
slightly more stable.

## 3. Fortune Voronoi construction

Fortune's algorithm moves a horizontal sweep line, called the **directrix**,
from top to bottom. The processed sites lie above the directrix. Their
parabolas separate points already known to be closer to a processed site from
points not reached by the sweep. The lower envelope of those parabolas is the
**beach line**.

Each maximal parabola portion on the envelope is an arc. A breakpoint between
two neighboring arcs traces a Voronoi edge because every point on it is
equidistant from the two corresponding sites.

![A site event splitting a beach-line arc and a circle event removing an
arc](fortune-events.svg)

The left panel shows the `old, new, old` beach-line replacement at a site
event. The right panel shows three sites becoming cocircular: the middle arc
shrinks to zero, its two breakpoints meet at a Voronoi vertex, and the outer
arcs become adjacent.

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

![The same beach-line arcs stored in geometric linked-list order and in an
AVL search tree](fortune-beach-line.svg)

The linked list makes local event updates explicit, while the AVL tree finds
the arc above a site's `x` coordinate in logarithmic time. Both structures
contain pointers to the same stable `Arc` objects.

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

![Finite, ray, and incomplete Voronoi edge records clipped against the
polygon bounding box](fortune-edge-completion.svg)

`completedBreakpoints` distinguishes a finite segment from a ray. The box
turns either form into the bounded segment consumed by polygon clipping.

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

## 4. Polygon clipping

A planar Voronoi diagram has only `O(n)` edges and vertices. With
`n <= 2000`, the simplest reliable clipping stage is fast enough:

1. Test every genuine Voronoi vertex with ray-casting point-in-polygon.
2. Intersect every clipped Voronoi edge with every polygon edge.
3. Evaluate the distance at every retained candidate.

Both operations take `O(n^2)` time. This avoids maintaining a geometric active
set whose ordering changes at polygon/Voronoi crossings.

`pointInPolygon()` treats boundary points as inside. For other points, it
casts a horizontal ray and toggles containment at every polygon edge crossing:

```cpp
if ((first.y() > point.y()) != (second.y() > point.y())) {
    const Real crossing_x =
        first.x() +
        (second.x() - first.x()) *
            (point.y() - first.y()) /
            (second.y() - first.y());
    if (crossing_x > point.x()) {
        inside = !inside;
    }
}
```

For every Voronoi edge, `clipPolygon()` first checks its real circle-event
vertices, excluding artificial bounding-box endpoints. It then tests the
bounded edge segment against all polygon edges:

```cpp
for (const VoronoiEdge& edge : voronoi_edges) {
    for (const Point& vertex : edge.vertices) {
        if (pointInPolygon(vertex, polygon)) {
            updateAnswer(vertex, edge.firstSite);
        }
    }
    for (const Segment& boundary : polygon_edges) {
        if (const auto intersection =
                segmentIntersection(edge.segment, boundary)) {
            updateAnswer(*intersection, edge.firstSite);
        }
    }
}
```

Every point on a Voronoi edge is nearest to either generating site, so storing
one of those site indices is sufficient for evaluating the objective.

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
test real Voronoi vertices for polygon containment
        |
        v
intersect every Voronoi edge with every polygon edge
        |
        v
take the largest squared candidate distance
        |
        v
print its square root
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

The argument can be separated into three lemmas.

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

### Lemma 3: Polygon clipping examines every necessary candidate

Every genuine Voronoi vertex is tested for polygon containment. Every
Voronoi-edge/polygon-boundary intersection is found by testing its pair of
segments. Therefore every candidate described by Lemma 2 is evaluated.

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
- Testing `O(n)` Voronoi vertices against an `n`-edge polygon takes
  `O(n^2)` time.
- Intersecting `O(n)` Voronoi edges with all `n` polygon edges also takes
  `O(n^2)` time.
- The event queue, beach line, and Voronoi records use `O(n)` memory.

Overall:

```text
Time:  O(n^2)
Space: O(n)
```
