#pragma once

#include <unordered_map>
#include <vector>
#include "geometry.hpp"
#include "distance.hpp"

// 构建物面节点的二叉树，用以查询空间节点到物面节点的最小距离
void _set_wall_tree(GridBTree<int, Point<double>>& global_wall_tree,
                    const std::vector<Node>& wall);

// 根据节点到物面的距离给节点分等级：
//   - 距离 <= D    → 1
//   - 距离 <= 5*D  → 2
//   - 其他         → 3
void _set_nd2wall_lvl_mp(std::unordered_map<int, int>& _mp,
                         const std::vector<Node>& all_nds,
                         const GridBTree<int, Point<double>>& global_wall_tree,
                         double D);
