#ifndef TRAPEZOID_TRAPEZOID_H
#define TRAPEZOID_TRAPEZOID_H

#include "common/BasePoint.h"
#include "common/BasePointUtil.h"
#include "common/CommonUtil.h"
#include "trapezoid/detail/TrapezoidTypes.h"

#include <algorithm>
#include <deque>
#include <limits>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

namespace trapezoid {

// TrapezoidalMap implements the randomized incremental construction of a
// trapezoidal decomposition for a set of non-crossing segments, following
// Seidel's algorithm (1991).  After build(), locate() answers point-location
// queries in expected O(log n) time.
template <typename T = double>
class TrapezoidalMap {
public:
    using Point = common::BasePoint<T>;
    using Seg = detail::Segment<T>;
    using Trap = detail::Trapezoid<T>;
    using Node = detail::DagNode<T>;

    // Build the trapezoidal map from a list of segment endpoint pairs.
    // Each pair (p, q) defines one segment; endpoints are normalised
    // internally so that left.x <= right.x.
    void build(std::vector<std::pair<Point, Point>> segments) {
        // 1. Normalise segments and collect bounding box.
        std::vector<Seg*> segment_ptrs;
        segment_ptrs.reserve(segments.size());
        for (std::size_t i = 0; i < segments.size(); ++i) {
            auto& [a, b] = segments[i];
            // Ensure left.x <= right.x (tie-break on y).
            if (a.x() > b.x() || (common::almostEqual(a.x(), b.x()) && a.y() > b.y())) {
                std::swap(a, b);
            }
            m_segments.push_back(Seg{a, b, i});
            segment_ptrs.push_back(&m_segments.back());
        }

        // 2. Random shuffle (key to expected O(n log n)).
        std::mt19937 rng(static_cast<std::mt19937::result_type>(42));
        std::shuffle(segment_ptrs.begin(), segment_ptrs.end(), rng);

        // 3. Compute bounding box and initialise.
        initializeBoundingBox(segment_ptrs);

        // 4. Insert segments one by one.
        for (auto* s : segment_ptrs) {
            insertSegment(s);
        }

        // 5. Collect active (non-removed) trapezoids.
        for (auto& t : m_trapezoids) {
            if (!t.removed) {
                m_active.push_back(&t);
            }
        }
    }

    // Point-location query: return the trapezoid containing q, or nullptr
    // if q lies outside the bounding box.
    const Trap* locate(Point q) const {
        if (!m_root) {
            return nullptr;
        }
        const auto* node = m_root;
        while (node->type != detail::DagNodeType::LEAF) {
            if (node->type == detail::DagNodeType::X) {
                // Compare x-coordinates; on exact tie use y.
                if (common::almostLess(q.x(), node->point->x())) {
                    node = node->left;
                } else if (common::almostGreater(q.x(), node->point->x())) {
                    node = node->right;
                } else {
                    node = (q.y() < node->point->y()) ? node->left : node->right;
                }
            } else { // Y-node
                if (pointAbove(q, *node->segment)) {
                    node = node->left;
                } else {
                    node = node->right;
                }
            }
        }
        return node->trapezoid;
    }

    // All active trapezoids after build().
    const std::vector<Trap*>& trapezoids() const { return m_active; }

    // --- Diagnostics (mainly for testing) ---
    std::size_t dagDepth() const {
        return dagDepthHelper(m_root);
    }

private:
    // Object pools — deque guarantees pointer stability.
    std::deque<Seg> m_segments;
    std::deque<Point> m_points;
    std::deque<Trap> m_trapezoids;
    std::deque<Node> m_dag_nodes;

    Node* m_root = nullptr;
    std::vector<Trap*> m_active;

    // Sentinel segments for the bounding-box top and bottom edges.
    Seg* m_top_sentinel = nullptr;
    Seg* m_bottom_sentinel = nullptr;

    // --- Geometric predicates ---

    // Returns true if point q lies strictly above segment s.
    static bool pointAbove(const Point& q, const Seg& s) {
        // cross > 0  =>  q is to the left of s (p->q direction),
        // which in our left-to-right convention means "above".
        const auto cr = common::cross(s.right - s.left, q - s.left);
        return cr > 0;
    }

    // --- Object allocation ---

    Seg* allocSeg(Seg s) {
        m_segments.push_back(std::move(s));
        return &m_segments.back();
    }

    Point* allocPoint(Point p) {
        m_points.push_back(std::move(p));
        return &m_points.back();
    }

    Trap* allocTrap() {
        m_trapezoids.push_back(Trap{});
        return &m_trapezoids.back();
    }

    Node* allocNode() {
        m_dag_nodes.push_back(Node{});
        return &m_dag_nodes.back();
    }

    Node* makeLeaf(Trap* t) {
        auto* n = allocNode();
        n->type = detail::DagNodeType::LEAF;
        n->trapezoid = t;
        t->dag_node = n;
        return n;
    }

    Node* makeXNode(Point* p, Node* left, Node* right) {
        auto* n = allocNode();
        n->type = detail::DagNodeType::X;
        n->point = p;
        n->left = left;
        n->right = right;
        return n;
    }

    Node* makeYNode(Seg* s, Node* left, Node* right) {
        auto* n = allocNode();
        n->type = detail::DagNodeType::Y;
        n->segment = s;
        n->left = left;
        n->right = right;
        return n;
    }

    // --- Bounding box initialisation ---

    void initializeBoundingBox(const std::vector<Seg*>& segs) {
        T min_x = std::numeric_limits<T>::max();
        T max_x = std::numeric_limits<T>::lowest();
        T min_y = std::numeric_limits<T>::max();
        T max_y = std::numeric_limits<T>::lowest();
        for (const auto* s : segs) {
            min_x = std::min({min_x, s->left.x(), s->right.x()});
            max_x = std::max({max_x, s->left.x(), s->right.x()});
            min_y = std::min({min_y, s->left.y(), s->right.y()});
            max_y = std::max({max_y, s->left.y(), s->right.y()});
        }
        // Add margin so all segments are strictly inside.
        const auto margin = T{1};
        min_x -= margin;
        max_x += margin;
        min_y -= margin;
        max_y += margin;

        // Four bounding-box corner points.
        auto* bl = allocPoint(Point{min_x, min_y}); // bottom-left
        auto* br = allocPoint(Point{max_x, min_y}); // bottom-right
        auto* tl = allocPoint(Point{min_x, max_y}); // top-left
        auto* tr = allocPoint(Point{max_x, max_y}); // top-right

        // Sentinel segments for top and bottom edges.
        m_top_sentinel = allocSeg(Seg{*tl, *tr, static_cast<std::size_t>(-1)});
        m_bottom_sentinel = allocSeg(Seg{*bl, *br, static_cast<std::size_t>(-2)});

        // Single trapezoid = entire bounding box.
        auto* t = allocTrap();
        t->top = m_top_sentinel;
        t->bottom = m_bottom_sentinel;
        t->leftp = bl;
        t->rightp = br;
        m_root = makeLeaf(t);
    }

    // --- Segment insertion ---

    void insertSegment(Seg* s) {
        // Locate the trapezoid containing s->left.
        auto* start = locatePointInternal(*s->left);

        // Collect all trapezoids crossed by s.
        std::vector<Trap*> crossed;
        crossed.push_back(start);
        auto* curr = start;
        while (common::almostGreater(s->right.x(), curr->rightp->x())) {
            // Determine whether s exits curr through its upper-right or
            // lower-right neighbour, by comparing the segment's y value
            // at x = curr->rightp->x with rightp itself.
            const auto seg_y = s->evalY(curr->rightp->x());
            if (common::almostGreater(seg_y, curr->rightp->y())) {
                curr = curr->upper_right;
            } else {
                curr = curr->lower_right;
            }
            crossed.push_back(curr);
        }

        if (crossed.size() == 1) {
            insertInSingleTrapezoid(s, crossed[0]);
        } else {
            insertInMultipleTrapezoids(s, crossed);
        }
    }

    Trap* locatePointInternal(const Point& q) {
        auto* node = m_root;
        while (node->type != detail::DagNodeType::LEAF) {
            if (node->type == detail::DagNodeType::X) {
                if (common::almostLess(q.x(), node->point->x())) {
                    node = node->left;
                } else if (common::almostGreater(q.x(), node->point->x())) {
                    node = node->right;
                } else {
                    node = (q.y() < node->point->y()) ? node->left : node->right;
                }
            } else {
                if (pointAbove(q, *node->segment)) {
                    node = node->left;
                } else {
                    node = node->right;
                }
            }
        }
        return node->trapezoid;
    }

    // Case 1: segment s lies entirely inside a single trapezoid delta.
    //
    //   Before:          After:
    //   +--------+       +--+------+--+
    //   |        |       |L |  A   |R |
    //   |  delta |  -->  |  +------+  |
    //   |        |       |  |  B   |  |
    //   +--------+       +--+------+--+
    //                     s.left   s.right
    void insertInSingleTrapezoid(Seg* s, Trap* delta) {
        auto* A = allocTrap(); // above s
        auto* B = allocTrap(); // below s
        auto* L = allocTrap(); // left
        auto* R = allocTrap(); // right

        A->top = delta->top;
        A->bottom = s;
        A->leftp = s->left.left;  // TODO: use allocPoint
        A->rightp = s->left.right;

        // We need actual point objects for leftp/rightp.
        auto* sleft = allocPoint(s->left.left);
        auto* sright = allocPoint(s->left.right);
        A->leftp = sleft;
        A->rightp = sright;

        B->top = s;
        B->bottom = delta->bottom;
        B->leftp = sleft;
        B->rightp = sright;

        L->top = delta->top;
        L->bottom = delta->bottom;
        L->leftp = delta->leftp;
        L->rightp = sleft;

        R->top = delta->top;
        R->bottom = delta->bottom;
        R->leftp = sright;
        R->rightp = delta->rightp;

        // Neighbour links.
        // L inherits delta's left neighbours.
        L->upper_left = delta->upper_left;
        L->lower_left = delta->lower_left;
        L->upper_right = A;
        L->lower_right = B;
        patchLeftNeighbours(delta, L, L);

        // R inherits delta's right neighbours.
        R->upper_right = delta->upper_right;
        R->lower_right = delta->lower_right;
        R->upper_left = A;
        R->lower_left = B;
        patchRightNeighbours(delta, R, R);

        // A and B are between L and R.
        A->upper_left = L;
        A->lower_left = nullptr;
        A->upper_right = R;
        A->lower_right = nullptr;

        B->upper_left = nullptr;
        B->lower_left = L;
        B->upper_right = nullptr;
        B->lower_right = R;

        // Mark old trapezoid as removed.
        delta->removed = true;

        // Update DAG: replace delta's leaf with a sub-DAG.
        auto* leaf = delta->dag_node;
        auto* y_node = makeYNode(s, makeLeaf(A), makeLeaf(B));
        auto* right_x = makeXNode(sright, y_node, makeLeaf(R));
        leaf->type = detail::DagNodeType::X;
        leaf->point = sleft;
        leaf->left = makeLeaf(L);
        leaf->right = right_x;
        leaf->trapezoid = nullptr;
    }

    // Case 2: segment s crosses multiple trapezoids crossed[0..k-1].
    void insertInMultipleTrapezoids(Seg* s, std::vector<Trap*>& crossed) {
        const auto k = crossed.size();

        // Allocate the left/right endpoint points once.
        auto* sleft = allocPoint(s->left.left);
        auto* sright = allocPoint(s->left.right);

        // We build upper/lower trapezoid chains, merging adjacent pieces
        // that share the same top/bottom.
        //
        // For each crossed trapezoid we create an "upper" piece (above s)
        // and a "lower" piece (below s).  The first and last also get an
        // extra L / R trapezoid.

        // --- First trapezoid (crossed[0]) ---
        auto* delta0 = crossed[0];

        auto* L = allocTrap();
        L->top = delta0->top;
        L->bottom = delta0->bottom;
        L->leftp = delta0->leftp;
        L->rightp = sleft;
        L->upper_left = delta0->upper_left;
        L->lower_left = delta0->lower_left;
        patchLeftNeighbours(delta0, L, L);

        // Upper piece of first trapezoid.
        auto* upper_curr = allocTrap();
        upper_curr->top = delta0->top;
        upper_curr->bottom = s;
        upper_curr->leftp = sleft;
        upper_curr->rightp = delta0->rightp;
        upper_curr->upper_left = L;
        upper_curr->lower_left = nullptr;

        // Lower piece of first trapezoid.
        auto* lower_curr = allocTrap();
        lower_curr->top = s;
        lower_curr->bottom = delta0->bottom;
        lower_curr->leftp = sleft;
        lower_curr->rightp = delta0->rightp;
        lower_curr->upper_left = nullptr;
        lower_curr->lower_left = L;

        L->upper_right = upper_curr;
        L->lower_right = lower_curr;

        // Replace delta0's DAG leaf.
        replaceLeafWithYNode(delta0, s, upper_curr, lower_curr);
        delta0->removed = true;

        // --- Middle trapezoids (crossed[1] .. crossed[k-2]) ---
        for (std::size_t i = 1; i < k - 1; ++i) {
            auto* delta = crossed[i];

            // Check if we can merge: same top (for upper) / same bottom (for lower).
            bool merge_upper = (delta->top == upper_curr->top);
            bool merge_lower = (delta->bottom == lower_curr->bottom);

            if (merge_upper) {
                // Extend upper_curr to the right.
                upper_curr->rightp = delta->rightp;
            } else {
                // Finalise upper_curr, start a new upper piece.
                auto* new_upper = allocTrap();
                new_upper->top = delta->top;
                new_upper->bottom = s;
                new_upper->leftp = delta->leftp;
                new_upper->rightp = delta->rightp;
                // Link: upper_curr's right neighbour is new_upper.
                upper_curr->upper_right = delta->upper_right; // carry old upper-right
                new_upper->upper_left = upper_curr->upper_right;
                new_upper->lower_left = nullptr;
                upper_curr = new_upper;
            }

            if (merge_lower) {
                lower_curr->rightp = delta->rightp;
            } else {
                auto* new_lower = allocTrap();
                new_lower->top = s;
                new_lower->bottom = delta->bottom;
                new_lower->leftp = delta->leftp;
                new_lower->rightp = delta->rightp;
                lower_curr->lower_right = delta->lower_right;
                new_lower->lower_left = lower_curr->lower_right;
                new_lower->upper_left = nullptr;
                lower_curr = new_lower;
            }

            // If neither merged, we created new trapezoids for this delta.
            if (!merge_upper && !merge_lower) {
                replaceLeafWithYNode(delta, s, upper_curr, lower_curr);
            } else if (merge_upper && merge_lower) {
                // Both merged — delta's leaf is now stale, just remove.
                delta->removed = true;
                // The DAG leaf for delta needs to point to the merged
                // upper/lower trapezoids.
                replaceLeafWithYNode(delta, s, upper_curr, lower_curr);
            } else if (merge_upper) {
                replaceLeafWithYNode(delta, s, upper_curr, lower_curr);
            } else {
                replaceLeafWithYNode(delta, s, upper_curr, lower_curr);
            }

            delta->removed = true;
        }

        // --- Last trapezoid (crossed[k-1]) ---
        auto* deltaN = crossed[k - 1];

        auto* R = allocTrap();
        R->top = deltaN->top;
        R->bottom = deltaN->bottom;
        R->leftp = sright;
        R->rightp = deltaN->rightp;
        R->upper_right = deltaN->upper_right;
        R->lower_right = deltaN->lower_right;
        patchRightNeighbours(deltaN, R, R);

        // Extend / finalise upper and lower chains.
        bool merge_upper = (deltaN->top == upper_curr->top);
        bool merge_lower = (deltaN->bottom == lower_curr->bottom);

        if (merge_upper) {
            upper_curr->rightp = sright;
        } else {
            auto* new_upper = allocTrap();
            new_upper->top = deltaN->top;
            new_upper->bottom = s;
            new_upper->leftp = deltaN->leftp;
            new_upper->rightp = sright;
            new_upper->upper_left = upper_curr;
            new_upper->lower_left = nullptr;
            upper_curr->upper_right = new_upper;
            upper_curr = new_upper;
        }

        if (merge_lower) {
            lower_curr->rightp = sright;
        } else {
            auto* new_lower = allocTrap();
            new_lower->top = s;
            new_lower->bottom = deltaN->bottom;
            new_lower->leftp = deltaN->leftp;
            new_lower->rightp = sright;
            new_lower->upper_left = nullptr;
            new_lower->lower_left = lower_curr;
            lower_curr->lower_right = new_lower;
            lower_curr = new_lower;
        }

        // Connect upper/lower to R.
        upper_curr->upper_right = R;
        lower_curr->lower_right = R;
        R->upper_left = upper_curr;
        R->lower_left = lower_curr;

        replaceLeafWithYNode(deltaN, s, upper_curr, lower_curr);
        deltaN->removed = true;
    }

    // --- Neighbour patching helpers ---

    // In all left-neighbours of old, replace pointers to old with
    // upper_new (for upper links) and lower_new (for lower links).
    static void patchLeftNeighbours(Trap* old, Trap* upper_new, Trap* lower_new) {
        if (old->upper_left) {
            if (old->upper_left->upper_right == old)
                old->upper_left->upper_right = upper_new;
            if (old->upper_left->lower_right == old)
                old->upper_left->lower_right = lower_new;
        }
        if (old->lower_left) {
            if (old->lower_left->upper_right == old)
                old->lower_left->upper_right = upper_new;
            if (old->lower_left->lower_right == old)
                old->lower_left->lower_right = lower_new;
        }
    }

    // In all right-neighbours of old, replace pointers to old with
    // upper_new / lower_new.
    static void patchRightNeighbours(Trap* old, Trap* upper_new, Trap* lower_new) {
        if (old->upper_right) {
            if (old->upper_right->upper_left == old)
                old->upper_right->upper_left = upper_new;
            if (old->upper_right->lower_left == old)
                old->upper_right->lower_left = lower_new;
        }
        if (old->lower_right) {
            if (old->lower_right->upper_left == old)
                old->lower_right->upper_left = upper_new;
            if (old->lower_right->lower_left == old)
                old->lower_right->lower_left = lower_new;
        }
    }

    // Replace a trapezoid's DAG leaf with a Y-node (segment) whose
    // left child is the upper trapezoid and right child is the lower.
    void replaceLeafWithYNode(Trap* delta, Seg* s, Trap* upper, Trap* lower) {
        auto* leaf = delta->dag_node;
        auto* upper_leaf = makeLeaf(upper);
        auto* lower_leaf = makeLeaf(lower);

        // Convert the leaf node in-place to a Y-node.
        leaf->type = detail::DagNodeType::Y;
        leaf->segment = s;
        leaf->left = upper_leaf;
        leaf->right = lower_leaf;
        leaf->trapezoid = nullptr;
    }

    // --- Diagnostics ---
    std::size_t dagDepthHelper(const Node* node) const {
        if (!node || node->type == detail::DagNodeType::LEAF) {
            return 0;
        }
        return 1 + std::max(dagDepthHelper(node->left), dagDepthHelper(node->right));
    }
};

} // namespace trapezoid

#endif // TRAPEZOID_TRAPEZOID_H
