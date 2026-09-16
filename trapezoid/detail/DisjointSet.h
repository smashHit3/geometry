#ifndef TRAPEZOID_DETAIL_DISJOINTSET_H
#define TRAPEZOID_DETAIL_DISJOINTSET_H

#include <cstddef>
#include <numeric>
#include <vector>

namespace trapezoid::detail {

// Union-find with path compression and union by size.
class DisjointSet {
public:
    explicit DisjointSet(std::size_t count = 0)
        : m_parent(count), m_size(count, 1) {
        std::iota(m_parent.begin(), m_parent.end(), std::size_t{0});
    }

    std::size_t add() {
        const std::size_t id = m_parent.size();
        m_parent.push_back(id);
        m_size.push_back(1);
        return id;
    }

    std::size_t find(std::size_t element) {
        std::size_t root = element;
        while (m_parent[root] != root) {
            root = m_parent[root];
        }
        while (m_parent[element] != element) {
            const std::size_t next = m_parent[element];
            m_parent[element] = root;
            element = next;
        }
        return root;
    }

    void unite(std::size_t first, std::size_t second) {
        first = find(first);
        second = find(second);
        if (first == second) {
            return;
        }
        if (m_size[first] < m_size[second]) {
            std::swap(first, second);
        }
        m_parent[second] = first;
        m_size[first] += m_size[second];
    }

private:
    std::vector<std::size_t> m_parent;
    std::vector<std::size_t> m_size;
};

} // namespace trapezoid::detail

#endif // TRAPEZOID_DETAIL_DISJOINTSET_H
