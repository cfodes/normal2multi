#include "metis_block.hpp"
#include "block.hpp"
#include <unordered_set>
#include <iostream>
#include <cassert>

using namespace std;

// ====================== 构造函数 ======================
// 构造函数，默认分区数k为2（传入边界进行构造）
explicit metis_block::metis_block(const boundary &b, const std::vector<Node> &wall_nodes, idx_t k = 2)
{
    // 参数设置
    m_k = k;
    ne = static_cast<idx_t>(b.bound_elements.size());
    nn = wall_nodes.size();

    // 构建外部节点编号和由0开始的节点编号的一一映射
    std::unordered_map<idx_t, idx_t> node_id_2_0based_id;
    for (idx_t i = 0; i < nn; ++i)
    {
        node_id_2_0based_id[wall_nodes[i].id] = i;
    }

    // 构造m_eptr/m_eind

    // 计算m_eptr元素
    m_eptr.resize(ne + 1);
    m_eptr[0] = 0;
    for (idx_t i = 0; i < ne; i++)
    {
        m_eptr[i + 1] = m_eptr[i] + b.bound_elements[i].node_id.size();
    }
    // 构建m_eind数组：
    m_eind.resize(m_eptr[ne]); // m_eptr最后一个元素对应m_eind的大小
    m_eind_temp.resize(m_eptr[ne]);
    idx_t m_eind_id = 0;
    for (idx_t i = 0; i < ne; ++i)
    {
        for (auto node : b.bound_elements[i].node_id) // 不同类型单元的节点数不一样
        {
            // m_eind保存外部存储的节点编号
            m_eind[m_eind_id] = static_cast<idx_t>(node); // bound_elements[i].node_id是int类型，这里做了一个安全的类型转换

            // m_eind保存0based的节点编号
            auto it = node_id_2_0based_id.find(node);
            if (it == node_id_2_0based_id.end()) // 检查，若该节点不存在于wall_node里，程序出错
            {
                std::cerr << "Error: Node ID " << node << " not found in wall_nodes!" << std::endl;
                exit(EXIT_FAILURE);
            }
            // 节点存在的情况：
            m_eind_temp[m_eind_id] = it->second;

            m_eind_id++;
        }
    }
}

// 构造函数，默认分区数k为2（传入单元进行构造）
explicit metis_block::metis_block(const std::vector<element> &b, const std::vector<Node> &wall_nodes, idx_t k = 2)
{
    // 参数设置
    m_k = k;
    ne = static_cast<idx_t>(b.size());
    nn = wall_nodes.size();

    // 构建外部节点编号和由0开始的节点编号的一一映射
    std::unordered_map<idx_t, idx_t> node_id_2_0based_id;
    for (idx_t i = 0; i < nn; ++i)
    {
        node_id_2_0based_id[wall_nodes[i].id] = i;
    }

    // 构造m_eptr/m_eind

    // 计算m_eptr元素
    m_eptr.resize(ne + 1);
    m_eptr[0] = 0;
    for (idx_t i = 0; i < ne; i++)
    {
        m_eptr[i + 1] = m_eptr[i] + b[i].node_id.size();
    }
    // 构建m_eind数组：
    m_eind.resize(m_eptr[ne]); // m_eptr最后一个元素对应m_eind的大小
    m_eind_temp.resize(m_eptr[ne]);
    idx_t m_eind_id = 0;
    for (idx_t i = 0; i < ne; ++i)
    {
        for (auto node : b[i].node_id) // 不同类型单元的节点数不一样
        {
            // m_eind保存外部存储的节点编号
            m_eind[m_eind_id] = static_cast<idx_t>(node); // bound_elements[i].node_id是int类型，这里做了一个安全的类型转换

            // m_eind保存0based的节点编号
            auto it = node_id_2_0based_id.find(node);
            if (it == node_id_2_0based_id.end()) // 检查，若该节点不存在于wall_node里，程序出错
            {
                std::cerr << "Error: Node ID " << node << " not found in wall_nodes!" << std::endl;
                exit(EXIT_FAILURE);
            }
            // 节点存在的情况：
            m_eind_temp[m_eind_id] = it->second;

            m_eind_id++;
        }
    }
}

// ====================== METIS 分区函数 ======================
// 使用metis函数进行分区，默认边共享ncommon=2
void metis_block::partitionMeshDual(idx_t ncommon = 2)
{
    assert(ne > 0); // 不允许传入空的单元
    m_epart.resize(ne);
    m_npart.resize(nn);

    // objval：METIS 返回的分区目标值（如切分边数或通信量等）
    idx_t objval = 0;
    int status = METIS_PartMeshDual(
        &ne,                // [输入] 单元数
        &nn,                // [输入] 节点数
        m_eptr.data(),      // [输入] eptr
        m_eind_temp.data(), // [输入] eind，需要0based的节点编号
        nullptr,            // [输入] vwgt (单元权重)，此处不使用
        nullptr,            // [输入] vsize (单元大小)，此处不使用
        &ncommon,           // [输入] ncommon
        &m_k,               // [输入] nparts
        nullptr,            // [输入] tpwgts (各分区目标权重)，此处不使用
        nullptr,            // [输入] options (高级选项)，此处不使用
        &objval,            // [输出] 分区的目标值
        m_epart.data(),     // [输出] 每个单元的分区号
        m_npart.data()      // [输出] 每个节点的分区号
    );

    if (status != METIS_OK)
    {
        std::cerr << "METIS error: " << status << std::endl;
    }
}

// ====================== 根据分区结果生成网格块 ======================
// 从boundary构建网格块数组
void metis_block::generate_mesh_blocks(const boundary &b, const std::vector<Node> &wall_nodes,
                                       std::vector<mesh_block> &blocks, std::vector<Node> &unique_boundary_nodes)
// 根据m_k/m_epart/m_npart构建网格块数组
// 同时将网格块之间共享的边界节点唯一地保存在_unique_boundary_nodes
// 传入参数：boundary
//          wall_nodes
//          blocks
//          unique_boundary_nodes
//!!!boundary里边界单元数组每个单元索引和m_eptr/m_epart的索引是一一对应的
{
    idx_t part_id = 0;  // 用于循环里临时记录分区号
    int inode_id = 0;   // 用于循环里临时记录节点id
    blocks.resize(m_k); // 根据分区数确定网格块数组的大小
    for (idx_t i = 0; i < m_k; ++i)
    {
        blocks[i].block_id = i; // 设定每个网格块的分区号
    }
    // 遍历wall，根据METIS的m_epart将单元分配到对应的块中
    for (size_t e = 0; e < b.bound_elements.size(); e++)
    {
        part_id = m_epart[e];
        blocks[part_id].internal_elements.push_back(b.bound_elements[e]);
        // b.bound_elements的索引和m_epart是一一对应的
        // m_epart[e]是索引为e的单元的分区号
    }

    // 首先建立一个从节点编号到其相邻单元分区集合的映射
    // key 为节点编号，value 为该节点所在的所有单元的分区集合
    unordered_map<int, unordered_set<idx_t>> node_partition_map;
    for (size_t e = 0; e < b.bound_elements.size(); e++)
    {
        part_id = m_epart[e];
        for (auto node : b.bound_elements[e].node_id)
        {
            node_partition_map[node].insert(part_id); // 在map的每个节点插入所属单元的分区编号，unordered_set确保这些编号是唯一存在的
        }
    }
    // 为了方便快速查找每个节点的坐标，利用 wall_nodes 建立一个映射：
    // key 为节点编号，value 为对应的 Point<double> 对象（包含坐标）
    unordered_map<int, Node> wall_node_map;
    for (const auto &p : wall_nodes)
    {
        wall_node_map[p.id] = p;
    }
    // 接下来遍历映射，对于每个节点，如果它只出现于一个分区，则视为内部节点；
    // 如果它出现在多个分区（即共享节点），则从 wall_node_map 中获取该节点的完整信息
    for (auto &entry : node_partition_map)
    {
        inode_id = entry.first;      // entry是一个pair<node_id,node_part_id>
        if (entry.second.size() > 1) // 该点是共享的（出现在多个分区）
        {
            if (wall_node_map.find(inode_id) != wall_node_map.end()) // 检查 wall_node_map 中是否存在该节点
            {
                const Node &node = wall_node_map[inode_id];
                // 将该节点插入对应分区的 boundary_nodes 里（可能属于多个分区）
                for (auto block_id : entry.second)
                {
                    blocks[block_id].boundary_nodes.push_back(node);
                    blocks[block_id].m_nodes.push_back(node); // 边界节点插入该分区的成员节点中
                }
                // 同时将共享节点添加到全局 unique_boundary_nodes 中
                // 首先检查该节点是否已经存在（这里依据节点 id 判断）
                bool found = false;
                for (const auto &p : unique_boundary_nodes)
                {
                    if (p.id == node.id)
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    unique_boundary_nodes.push_back(node);
                }
            }
            else // 物面没有该节点，说明程序出现错误
            {
                std::cerr << "Error: cannot find corresponding node " << inode_id
                          << " in wall_nodes. Partition node and wall nodes do not match!"
                          << std::endl;
                exit(EXIT_FAILURE); // 或者 throw std::runtime_error("..."); 根据需求选择中断方式
            }
        }
        else // 该点只属于一个分区，将其插入对应分区的内部节点里
        {
            part_id = *(entry.second.begin());
            if (wall_node_map.find(inode_id) != wall_node_map.end())
            {
                blocks[part_id].internal_nodes.push_back(wall_node_map[inode_id]);
                blocks[part_id].m_nodes.push_back(wall_node_map[inode_id]); // 内部节点插入该分区的成员节点中
            }

            else
            {
                std::cerr << "Error: cannot find corresponding node " << inode_id
                          << " in wall_nodes. Partition node and wall nodes do not match!"
                          << std::endl;
                exit(EXIT_FAILURE); // 或者 throw std::runtime_error("..."); 根据需求选择中断方式
            }
        }
    }
}

// 从边界单元elements构建网格块数组
void metis_block::generate_mesh_blocks(const std::vector<element> &b, const vector<Node> &wall_nodes, vector<mesh_block> &blocks, vector<Node> &unique_boundary_nodes)
// 根据m_k/m_epart/m_npart构建网格块数组
// 同时将网格块之间共享的边界节点唯一地保存在_unique_boundary_nodes
// 传入参数：boundary
//          wall_nodes
//          blocks
//          unique_boundary_nodes
//!!!boundary里边界单元数组每个单元索引和m_eptr/m_epart的索引是一一对应的
{
    idx_t part_id = 0;  // 用于循环里临时记录分区号
    int inode_id = 0;   // 用于循环里临时记录节点id
    blocks.resize(m_k); // 根据分区数确定网格块数组的大小
    for (idx_t i = 0; i < m_k; ++i)
    {
        blocks[i].block_id = i; // 设定每个网格块的分区号
    }
    // 遍历wall，根据METIS的m_epart将单元分配到对应的块中
    for (size_t e = 0; e < b.size(); e++)
    {
        part_id = m_epart[e];
        blocks[part_id].internal_elements.push_back(b[e]);
        // b.bound_elements的索引和m_epart是一一对应的
        // m_epart[e]是索引为e的单元的分区号
    }

    // 首先建立一个从节点编号到其相邻单元分区集合的映射
    // key 为节点编号，value 为该节点所在的所有单元的分区集合
    unordered_map<int, unordered_set<idx_t>> node_partition_map;
    for (size_t e = 0; e < b.size(); e++)
    {
        part_id = m_epart[e];
        for (auto node : b[e].node_id)
        {
            node_partition_map[node].insert(part_id); // 在map的每个节点插入所属单元的分区编号，unordered_set确保这些编号是唯一存在的
        }
    }
    // 为了方便快速查找每个节点的坐标，利用 wall_nodes 建立一个映射：
    // key 为节点编号，value 为对应的 Point<double> 对象（包含坐标）
    unordered_map<int, Node> wall_node_map;
    for (const auto &p : wall_nodes)
    {
        wall_node_map[p.id] = p;
    }
    // 接下来遍历映射，对于每个节点，如果它只出现于一个分区，则视为内部节点；
    // 如果它出现在多个分区（即共享节点），则从 wall_node_map 中获取该节点的完整信息
    for (auto &entry : node_partition_map)
    {
        inode_id = entry.first;      // entry是一个pair<node_id,node_part_id>
        if (entry.second.size() > 1) // 该点是共享的（出现在多个分区）
        {
            if (wall_node_map.find(inode_id) != wall_node_map.end()) // 检查 wall_node_map 中是否存在该节点
            {
                const Node &node = wall_node_map[inode_id];
                // 将该节点插入对应分区的 boundary_nodes 里（可能属于多个分区）
                for (auto block_id : entry.second)
                {
                    blocks[block_id].boundary_nodes.push_back(node);
                    blocks[block_id].m_nodes.push_back(node);      // 边界节点插入该分区成员节点中
                    blocks[block_id].m_nodes_id.insert(node.id);   // 将边界节点的id存入成员节点id查找集
                    blocks[block_id].bndry_nds_id.insert(node.id); // 将边界节点的id存入边界节点id查找集
                }
                // 同时将共享节点添加到全局 unique_boundary_nodes 中
                // 首先检查该节点是否已经存在（这里依据节点 id 判断）
                bool found = false;
                for (const auto &p : unique_boundary_nodes)
                {
                    if (p.id == node.id)
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    unique_boundary_nodes.push_back(node);
                }
            }
            else // 物面没有该节点，说明程序出现错误
            {
                std::cerr << "Error: cannot find corresponding node " << inode_id
                          << " in wall_nodes. Partition node and wall nodes do not match!"
                          << std::endl;
                exit(EXIT_FAILURE); // 或者 throw std::runtime_error("..."); 根据需求选择中断方式
            }
        }
        else // 该点只属于一个分区，将其插入对应分区的内部节点里
        {
            part_id = *(entry.second.begin());
            if (wall_node_map.find(inode_id) != wall_node_map.end())
            {
                blocks[part_id].internal_nodes.push_back(wall_node_map[inode_id]);
                blocks[part_id].m_nodes.push_back(wall_node_map[inode_id]);    // 内部节点插入该分区成员节点中
                blocks[part_id].m_nodes_id.insert(wall_node_map[inode_id].id); // 将内部节点的id存入成员节点id查找集
            }

            else
            {
                std::cerr << "Error: cannot find corresponding node " << inode_id
                          << " in wall_nodes. Partition node and wall nodes do not match!"
                          << std::endl;
                exit(EXIT_FAILURE); // 或者 throw std::runtime_error("..."); 根据需求选择中断方式
            }
        }
    }
}