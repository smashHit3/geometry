#ifndef COMMON_RTREE_RTREE_H
#define COMMON_RTREE_RTREE_H

#include "common/geometry/Aabb.h"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace common::rtree {

template <typename id_type, typename coordinate_type>
class RTree {
public:
    using Bounds = common::geometry::Aabb<coordinate_type>;

    struct Entry {
        id_type id;
        Bounds bounds;
    };

    explicit RTree(std::size_t max_entries = 8)
        : m_max_entries(max_entries) {
        if (max_entries < 2) {
            throw std::invalid_argument("R-tree node capacity must be at least two");
        }
    }

    void insert(const id_type& id, const Bounds& bounds) {
        if (containsId(id)) {
            throw std::invalid_argument("R-tree already contains this id");
        }
        Entry entry{id, bounds};
        if (!m_root) {
            m_root = std::make_unique<Node>(true);
            m_root->entries.push_back(std::move(entry));
            refreshBounds(*m_root);
            return;
        }

        auto sibling = insertInto(*m_root, std::move(entry));
        if (sibling) {
            auto root = std::make_unique<Node>(false);
            root->children.push_back(std::move(m_root));
            root->children.push_back(std::move(sibling));
            refreshBounds(*root);
            m_root = std::move(root);
        }
    }

    bool remove(const id_type& id) {
        std::vector<Entry> entries;
        collectEntries(entries);
        const auto found = std::find_if(entries.begin(), entries.end(),
                                        [&id](const Entry& entry) {
                                            return entry.id == id;
                                        });
        if (found == entries.end()) {
            return false;
        }
        entries.erase(found);
        rebuild(entries);
        return true;
    }

    void update(const id_type& id, const Bounds& bounds) {
        if (!remove(id)) {
            throw std::out_of_range("R-tree does not contain this id");
        }
        insert(id, bounds);
    }

    std::vector<id_type> query(const Bounds& bounds) const {
        std::vector<id_type> result;
        if (m_root) {
            queryNode(*m_root, bounds, result);
        }
        return result;
    }

private:
    struct Node {
        explicit Node(bool is_leaf)
            : leaf(is_leaf) {}

        bool leaf;
        Bounds bounds{coordinate_type{}, coordinate_type{}, coordinate_type{},
                      coordinate_type{}};
        std::vector<Entry> entries;
        std::vector<std::unique_ptr<Node>> children;
    };

    bool containsId(const id_type& id) const {
        std::vector<Entry> entries;
        collectEntries(entries);
        return std::any_of(entries.begin(), entries.end(), [&id](const Entry& entry) {
            return entry.id == id;
        });
    }

    std::unique_ptr<Node> insertInto(Node& node, Entry entry) {
        if (node.leaf) {
            node.entries.push_back(std::move(entry));
        } else {
            auto* child = chooseChild(node, entry.bounds);
            auto sibling = insertInto(*child, std::move(entry));
            if (sibling) {
                node.children.push_back(std::move(sibling));
            }
        }
        refreshBounds(node);
        if (size(node) <= m_max_entries) {
            return nullptr;
        }
        return split(node);
    }

    Node* chooseChild(Node& node, const Bounds& bounds) const {
        return (*std::min_element(
                    node.children.begin(), node.children.end(),
                    [&bounds](const std::unique_ptr<Node>& left,
                              const std::unique_ptr<Node>& right) {
                        const auto left_growth =
                            left->bounds.merged(bounds).area() - left->bounds.area();
                        const auto right_growth =
                            right->bounds.merged(bounds).area() - right->bounds.area();
                        if (left_growth != right_growth) {
                            return left_growth < right_growth;
                        }
                        return left->bounds.area() < right->bounds.area();
                    }))
            .get();
    }

    std::unique_ptr<Node> split(Node& node) {
        auto sibling = std::make_unique<Node>(node.leaf);
        const auto split_index = size(node) / 2;
        if (node.leaf) {
            std::sort(node.entries.begin(), node.entries.end(),
                      [](const Entry& left, const Entry& right) {
                          return left.bounds.minX() + left.bounds.maxX() <
                                 right.bounds.minX() + right.bounds.maxX();
                      });
            sibling->entries.insert(
                sibling->entries.end(),
                std::make_move_iterator(node.entries.begin() +
                                        static_cast<std::ptrdiff_t>(split_index)),
                std::make_move_iterator(node.entries.end()));
            node.entries.erase(node.entries.begin() +
                                   static_cast<std::ptrdiff_t>(split_index),
                               node.entries.end());
        } else {
            std::sort(node.children.begin(), node.children.end(),
                      [](const std::unique_ptr<Node>& left,
                         const std::unique_ptr<Node>& right) {
                          return left->bounds.minX() + left->bounds.maxX() <
                                 right->bounds.minX() + right->bounds.maxX();
                      });
            sibling->children.insert(
                sibling->children.end(),
                std::make_move_iterator(node.children.begin() +
                                        static_cast<std::ptrdiff_t>(split_index)),
                std::make_move_iterator(node.children.end()));
            node.children.erase(node.children.begin() +
                                     static_cast<std::ptrdiff_t>(split_index),
                                node.children.end());
        }
        refreshBounds(node);
        refreshBounds(*sibling);
        return sibling;
    }

    static std::size_t size(const Node& node) {
        return node.leaf ? node.entries.size() : node.children.size();
    }

    static void refreshBounds(Node& node) {
        if (node.leaf) {
            node.bounds = node.entries.front().bounds;
            for (const auto& entry : node.entries) {
                node.bounds = node.bounds.merged(entry.bounds);
            }
            return;
        }
        node.bounds = node.children.front()->bounds;
        for (const auto& child : node.children) {
            node.bounds = node.bounds.merged(child->bounds);
        }
    }

    static void queryNode(const Node& node, const Bounds& bounds,
                          std::vector<id_type>& result) {
        if (!node.bounds.intersects(bounds)) {
            return;
        }
        if (node.leaf) {
            for (const auto& entry : node.entries) {
                if (entry.bounds.intersects(bounds)) {
                    result.push_back(entry.id);
                }
            }
            return;
        }
        for (const auto& child : node.children) {
            queryNode(*child, bounds, result);
        }
    }

    void collectEntries(std::vector<Entry>& entries) const {
        if (m_root) {
            collectNodeEntries(*m_root, entries);
        }
    }

    static void collectNodeEntries(const Node& node, std::vector<Entry>& entries) {
        if (node.leaf) {
            entries.insert(entries.end(), node.entries.begin(), node.entries.end());
            return;
        }
        for (const auto& child : node.children) {
            collectNodeEntries(*child, entries);
        }
    }

    void rebuild(const std::vector<Entry>& entries) {
        m_root.reset();
        for (const auto& entry : entries) {
            insert(entry.id, entry.bounds);
        }
    }

    std::size_t m_max_entries;
    std::unique_ptr<Node> m_root;
};

} // namespace common::rtree

#endif // COMMON_RTREE_RTREE_H
