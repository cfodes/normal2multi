#pragma once

#include <metis.h>
#include <vector>
#include <unordered_set>
#include "geometry.hpp"
#include "distance.hpp"

// 网格块类定义
class mesh_block 
{
public:
    idx_t block_id;                         // 分区号，也就是该块所属的分区编号
    double block_D;                         // 分区插值系统所用的D, 对应对应于DRRBF里的psi = 1-r/D
    std::vector<Node> m_nodes;              // 所有属于该分区的成员节点
    std::unordered_set<int> m_nodes_id;     //用以查找内部节点id的set 
    std::vector<element> internal_elements; // 内部网格单元：仅属于该分块的单元
    std::vector<Node> internal_nodes;       // 内部节点：只出现在本分块内的节点
    std::vector<Node> boundary_nodes;       // 边界节点：位于多个分块之间共享的节点
    std::unordered_set<int> bndry_nds_id;   //用以查找边界节点id的set 
    GridBTree<int, Point<double>> block_tree;   //根据内部节点构造的二叉树
    Point<double> aabb_min;         // block_tree的包围盒最小值，外层树构建所需
    Point<double> aabb_max;         // block_tree的包围盒最小值，外层树构建所需
    bool has_aabb = false;          // 用于判断block_tree是否为空

    // 默认构造函数
    mesh_block() = default;

    // 计算网格内部节点的待变形量 = 已知变形量 - 已变形的量（针对wall_nodes没有更新的情况）
    void set_internal_nodes_df(const std::vector<Node>& wall_nodes_0, const std::vector<Node>& wall_nodes, std::unordered_map<int, int>& id2node);

    // 计算网格内部节点的待变形量 = 已知变形量 - 已变形的量（针对wall_nodes已经更新的情况
    void set_blk_nodes_df(const std::vector<Node>& wall_prev, const std::unordered_map<int, int>& id2node);

    // 内部节点建立后，构建内部节点的二叉树，同时缓存AABB
    void set_block_tree();
};

// 根据新的物面坐标重置blocks的节点坐标
void Reset_blocks(std::vector<mesh_block>& blocks, const std::vector<Node>& wall_nodes, const std::unordered_map<int, int>& id2node);

// 构建一个基于_unique_boundary_nodes的二叉树，这是由于每个分区的二叉树并不包含这些点
GridBTree<int, Point<double>> set_unique_boundary_nodes_tree(const std::vector<Node>& unique_boundary_nodes);  