#include "block.hpp"
#include "distance.hpp"
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace std;

void mesh_block::set_internal_nodes_df(const std::vector<Node>& wall_nodes_0, const std::vector<Node>& wall_nodes, std::unordered_map<int, int>& id2node)
//计算网格内部节点的待变形量 = 已知变形量 - 已变形的量
//传入参数：wall_nodes的历史数据 wall_nodes_0
//          wall_nodes 
//          wall_nodes的id和node查找集id2node
//该函数针对wall_nodes没有更新的情况
{
        Eigen::Vector3d df_temp{ 0.0,0.0,0.0 };
        double df_max = 0.0;
        for (auto& inode : internal_nodes)
        {
            auto it = id2node.find(inode.id);   //查找该内部节点在物面节点集的索引
            if (it != id2node.end())
            {
                df_temp = wall_nodes_0[it->second].df - wall_nodes[it->second].df;  //df_new = df_0-df_1
                inode.df = df_temp;
                df_max = std::max(df_max, df_temp.norm());
                df_temp.setZero();
            }
        }
        for (auto& inode : boundary_nodes)
        {
            auto it = id2node.find(inode.id);   //查找该内部节点在物面节点集的索引
            if (it != id2node.end())
            {
                df_temp = wall_nodes_0[it->second].df - wall_nodes[it->second].df;  //df_new = df_0-df_1
                inode.df = df_temp;
                df_max = std::max(df_max, df_temp.norm());
                df_temp.setZero();
            }
        }
        block_D = 5 * df_max;
        std::cout << "=====================" << endl;
        std::cout << "block " << block_id << " block_D is determined as " << block_D << endl;

    }

void mesh_block::set_blk_nodes_df(const std::vector<Node>& wall_prev, const std::unordered_map<int, int>& id2node)
//计算网格内部节点的待变形量 = 已知变形量 - 已变形的量
//传入参数：wall_prev, 物面节点的坐标和待变形量都已更新
//          wall_prev的id和node查找集id2node
//该函数针对wall_nodes已经更新的情况，适用于multi_partition
{
        Eigen::Vector3d df_temp{ 0.0,0.0,0.0 };
        double df_max = 0.0;
        for (auto& inode : internal_nodes)
        {
            auto it = id2node.find(inode.id);   //查找该内部节点在物面节点集的索引
            if (it != id2node.end())
            {
                df_temp = wall_prev[it->second].df;  //df_new = df_0-df_1
                inode.df = df_temp;
                df_max = std::max(df_max, df_temp.norm());
                df_temp.setZero();
            }
        }
        for (auto& inode : boundary_nodes)
        {
            auto it = id2node.find(inode.id);   //查找该内部节点在物面节点集的索引
            if (it != id2node.end())
            {
                df_temp = wall_prev[it->second].df;  //df_new = df_0-df_1
                inode.df = df_temp;
                df_max = std::max(df_max, df_temp.norm());
                df_temp.setZero();
            }
        }
        block_D = 5 * df_max;
        //std::cout << "=====================" << endl;
        //std::cout << "block " << block_id << " block_D is determined as " << block_D << endl;

}

void mesh_block::set_block_tree()  //内部节点建立后，构建内部节点的二叉树，同时缓存AABB
{
        block_tree.setGridSize(internal_nodes.size());
        for (int i = 0; i < internal_nodes.size(); ++i)
        {
            block_tree.setGrid(i, internal_nodes[i].id, internal_nodes[i].point);  //构建二叉树，参数分别为索引、点的id和点
        }
        block_tree.construct();

        // 从block_tree的根节点读出AABB并缓存（无需额外扫描每个点）
        Point<double> bmin, bmax;
        if (block_tree.root_bound(bmin, bmax))
        {
            aabb_min = bmin;
            aabb_max = bmax;
            has_aabb = true;
        }
        else
        {
            has_aabb = false;    //空树不参与外部树构建
        }
}

void Reset_blocks(std::vector<mesh_block>& blocks, const std::vector<Node>& wall_nodes, const std::unordered_map<int, int>& id2node)
//根据新的变形结果重新设值blocks内部每个节点的坐标
//传入参数： 分块网格数组 blocks
//          物面节点数组 wall_nodes (需要使用变形后的节点构建)
//           wall_nodes的id和node查找集  id2node
{
    //遍历每个分块，更新内部节点和边界节点
    for (auto& blk : blocks)
    {
        //更新内部节点
        for (auto& node : blk.internal_nodes)
        {
            auto it = id2node.find(node.id);
            if (it != id2node.end())
            {
                node.point = wall_nodes[it->second].point;
            }
        }

        //更新边界节点
        for (auto& node : blk.boundary_nodes)
        {
            auto it = id2node.find(node.id);
            if (it != id2node.end())
            {
                node.point = wall_nodes[it->second].point;
            }
        }
    }
    
}

GridBTree<int, Point<double>> set_unique_boundary_nodes_tree(const std::vector<Node>& unique_boundary_nodes)
//该函数构建一个基于_unique_boundary_nodes的二叉树，这是由于每个分区的二叉树并不包含这些点
//传入参数： _unique_boundary_points
{
    GridBTree<int, Point<double>> tree_temp;
    tree_temp.setGridSize(unique_boundary_nodes.size());
    for (int i = 0; i < unique_boundary_nodes.size(); ++i)
    {
        tree_temp.setGrid(i, unique_boundary_nodes[i].id, unique_boundary_nodes[i].point);  //构建二叉树，参数分别为索引、点的id和点
    }
    tree_temp.construct();
    return tree_temp;
}