#pragma once
#include <vector>
#include <metis.h>
#include "geometry.hpp"

class mesh_block; //前置声明

class metis_block
    //封装metis对网格分区的功能
    //根据外部输入的数据建立metis可用的数据结构
    //完成网格分区并公开该数据
{
public:
    idx_t m_k;   //分区数
    idx_t ne;    //单元数
    idx_t nn;    //节点数
    std::vector<idx_t> m_eptr;   //metis的数据结构:m_eptr[i]表示第i个单元在m_eind[]中起始下标
    std::vector<idx_t> m_eind;   //m_eptr[i+1]表示第i个单元的终止下标+1, 该eind的节点编号采用外部的节点编号
    std::vector<idx_t> m_eind_temp;   //该m_eind_temp实现从0开始编号的节点编号，与m_eind的节点编号一一映射

    std::vector<idx_t> m_epart;  //metis的分区结果:m_epart[e] 表示第 e 个单元的分区号 (0 ~ k-1)
    std::vector<idx_t> m_npart;  //m_npart[n] 表示第 n 个节点的分区号 (0 ~ k-1)

    //构造函数，默认分区数k为2（传入边界进行构造）
    explicit metis_block(const boundary& b, const std::vector<Node>& wall_nodes, idx_t k = 2);

    //构造函数，默认分区数k为2（传入单元进行构造）
    explicit metis_block(const std::vector<element>& b, const std::vector<Node>& wall_nodes, idx_t k = 2);

    //使用metis函数进行分区，默认边共享ncommon=2
    void partitionMeshDual(idx_t ncommon = 2);

    //根据m_k/m_epart/m_npart构建网格块数组(传入边界)
    void generate_mesh_blocks(const boundary& b, const std::vector<Node>& wall_nodes, std::vector<mesh_block>& blocks, std::vector<Node>& unique_boundary_nodes);

    //根据m_k/m_epart/m_npart构建网格块数组(传入单元)
    void generate_mesh_blocks(const std::vector<element>& b,const std::vector<Node>& wall_nodes, std::vector<mesh_block>& blocks, std::vector<Node>& unique_boundary_nodes);
};