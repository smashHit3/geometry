#ifndef COMMON_BVH_H
#define COMMON_BVH_H

#include "common/Aabb.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace common {

template <typename id_type, typename coordinate_type>
class Bvh {
public:
    using Bounds = Aabb<coordinate_type>;

    struct Entry {
        id_type id;
        Bounds bounds;
    };

    explicit Bvh(std::size_t leaf_size = 4)
        : m_leaf_size(leaf_size) {
        if (leaf_size == 0) {
            throw std::invalid_argument("BVH leaf size must be positive");
        }
    }

    void build(std::vector<Entry> entries) {
        m_entries = std::move(entries);
        m_nodes.clear();
        if (!m_entries.empty()) {
            buildNode(0, m_entries.size());
        }
    }

    std::vector<id_type> query(const Bounds& bounds) const {
        std::vector<id_type> result;
        if (!m_nodes.empty()) {
            queryNode(0, bounds, result);
        }
        return result;
    }

private:
    // Internal BVH node. Leaf nodes reference a contiguous range in m_entries;
    // internal nodes reference their two child nodes in m_nodes.
    struct Node {
        Bounds bounds;       // Merged bounds of every entry in this subtree.
        std::size_t first;   // First entry index; meaningful for leaf nodes.
        std::size_t count;   // Number of entries; meaningful for leaf nodes.
        std::size_t left = 0;   // Left child node index; meaningful for internal nodes.
        std::size_t right = 0;  // Right child node index; meaningful for internal nodes.
        bool leaf;              // True when this node directly stores entry indices.
    };

    std::size_t buildNode(std::size_t first, std::size_t count) {
        const auto node_index = m_nodes.size();
        m_nodes.push_back({boundsFor(first, count), first, count, 0, 0, true});
        if (count <= m_leaf_size) {
            return node_index;
        }

        const auto bounds = m_nodes[node_index].bounds;
        const bool split_x = bounds.maxX() - bounds.minX() >=
                             bounds.maxY() - bounds.minY();
        const auto middle = first + count / 2;
        std::nth_element(
            m_entries.begin() + static_cast<std::ptrdiff_t>(first),
            m_entries.begin() + static_cast<std::ptrdiff_t>(middle),
            m_entries.begin() + static_cast<std::ptrdiff_t>(first + count),
            [split_x](const Entry& left, const Entry& right) {
                const auto left_center = split_x ? left.bounds.minX() + left.bounds.maxX()
                                                 : left.bounds.minY() + left.bounds.maxY();
                const auto right_center = split_x ? right.bounds.minX() + right.bounds.maxX()
                                                  : right.bounds.minY() + right.bounds.maxY();
                return left_center < right_center;
            });
        const auto left = buildNode(first, count / 2);
        const auto right = buildNode(middle, count - count / 2);
        m_nodes[node_index].left = left;
        m_nodes[node_index].right = right;
        m_nodes[node_index].leaf = false;
        return node_index;
    }

    Bounds boundsFor(std::size_t first, std::size_t count) const {
        auto bounds = m_entries[first].bounds;
        for (auto index = first + 1; index < first + count; ++index) {
            bounds = bounds.merged(m_entries[index].bounds);
        }
        return bounds;
    }

    void queryNode(std::size_t node_index, const Bounds& bounds,
                   std::vector<id_type>& result) const {
        const auto& node = m_nodes[node_index];
        if (!node.bounds.intersects(bounds)) {
            return;
        }
        if (node.leaf) {
            for (auto index = node.first; index < node.first + node.count; ++index) {
                if (m_entries[index].bounds.intersects(bounds)) {
                    result.push_back(m_entries[index].id);
                }
            }
            return;
        }
        queryNode(node.left, bounds, result);
        queryNode(node.right, bounds, result);
    }

    std::size_t m_leaf_size;
    std::vector<Entry> m_entries;
    std::vector<Node> m_nodes;
};

} // namespace common

#endif // COMMON_BVH_H
