#ifndef COMMON_UNIFORMGRID_UNIFORMGRID_H
#define COMMON_UNIFORMGRID_UNIFORMGRID_H

#include "common/geometry/Aabb.h"

#include <cmath>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace common::uniformgrid {

template <typename id_type, typename coordinate_type,
          typename hash_type = std::hash<id_type>>
class UniformGrid {
public:
    using Bounds = common::geometry::Aabb<coordinate_type>;

    explicit UniformGrid(coordinate_type cell_size)
        : m_cell_size(cell_size) {
        if (cell_size <= 0) {
            throw std::invalid_argument("Grid cell size must be positive");
        }
    }

    void insert(const id_type& id, const Bounds& bounds) {
        if (m_bounds_by_id.find(id) != m_bounds_by_id.end()) {
            throw std::invalid_argument("Grid already contains this id");
        }
        m_bounds_by_id.emplace(id, bounds);
        forEachCell(bounds, [&](const Cell& cell) { m_cells[cell].insert(id); });
    }

    bool remove(const id_type& id) {
        const auto found = m_bounds_by_id.find(id);
        if (found == m_bounds_by_id.end()) {
            return false;
        }
        forEachCell(found->second, [&](const Cell& cell) {
            const auto cell_found = m_cells.find(cell);
            cell_found->second.erase(id);
            if (cell_found->second.empty()) {
                m_cells.erase(cell_found);
            }
        });
        m_bounds_by_id.erase(found);
        return true;
    }

    void update(const id_type& id, const Bounds& bounds) {
        if (!remove(id)) {
            throw std::out_of_range("Grid does not contain this id");
        }
        insert(id, bounds);
    }

    std::vector<id_type> query(const Bounds& bounds) const {
        std::unordered_set<id_type, hash_type> result;
        forEachCell(bounds, [&](const Cell& cell) {
            const auto found = m_cells.find(cell);
            if (found == m_cells.end()) {
                return;
            }
            for (const auto& id : found->second) {
                if (m_bounds_by_id.at(id).intersects(bounds)) {
                    result.insert(id);
                }
            }
        });
        return {result.begin(), result.end()};
    }

private:
    struct Cell {
        long long x;
        long long y;

        bool operator==(const Cell& other) const {
            return x == other.x && y == other.y;
        }
    };

    struct CellHash {
        std::size_t operator()(const Cell& cell) const {
            return std::hash<long long>{}(cell.x) ^
                   (std::hash<long long>{}(cell.y) << 1U);
        }
    };

    template <typename function_type>
    void forEachCell(const Bounds& bounds, function_type&& function) const {
        const auto min_x = cellCoordinate(bounds.minX());
        const auto max_x = cellCoordinate(bounds.maxX());
        const auto min_y = cellCoordinate(bounds.minY());
        const auto max_y = cellCoordinate(bounds.maxY());
        for (auto x = min_x; x <= max_x; ++x) {
            for (auto y = min_y; y <= max_y; ++y) {
                function(Cell{x, y});
            }
        }
    }

    long long cellCoordinate(coordinate_type value) const {
        return static_cast<long long>(
            std::floor(static_cast<long double>(value) / m_cell_size));
    }

    coordinate_type m_cell_size;
    std::unordered_map<Cell, std::unordered_set<id_type, hash_type>, CellHash> m_cells;
    std::unordered_map<id_type, Bounds, hash_type> m_bounds_by_id;
};

} // namespace common::uniformgrid

#endif // COMMON_UNIFORMGRID_UNIFORMGRID_H
