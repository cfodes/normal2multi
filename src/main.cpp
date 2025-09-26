#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <unordered_set>  //std::unorder_set
#include <unordered_map>     //std::unordered_map
#include<string>
#include<sstream>
#include <list>
#include <stdio.h>
#include <iomanip>  //fixed and setprecision
#include<cmath>
#include<cassert>  //assert

#include <metis.h>
//#include "ADT.hpp"
#include "./distance.hpp"
#include "geometry.hpp"

#include <Eigen/Core>          
#include <Eigen/SVD>
#include <Eigen/LU>

#include <filesystem>
namespace fs = std::filesystem;
using namespace std;




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

    //构造函数，默认分区数k为2
    //传入参数：boundary
    //         _wall_nodes
    //根据传入参数构造m_eptr/m_eind
    explicit metis_block(const boundary& b, const vector<Node>& wall_nodes, idx_t k = 2)
    {
        //参数设置
        m_k = k;
        ne = static_cast<idx_t>(b.bound_elements.size());
        nn = wall_nodes.size();

        //构建外部节点编号和由0开始的节点编号的一一映射
        std::unordered_map<idx_t, idx_t> node_id_2_0based_id;
        for (idx_t i = 0; i < nn; ++i)
        {
            node_id_2_0based_id[wall_nodes[i].id] = i;
        }

        //构造m_eptr/m_eind

        //计算m_eptr元素
        m_eptr.resize(ne + 1);
        m_eptr[0] = 0;
        for (idx_t i = 0; i < ne; i++)
        {
            m_eptr[i + 1] = m_eptr[i] + b.bound_elements[i].node_id.size();
        }
        //构建m_eind数组：
        m_eind.resize(m_eptr[ne]);        //m_eptr最后一个元素对应m_eind的大小
        m_eind_temp.resize(m_eptr[ne]);  
        idx_t m_eind_id = 0;
        for (idx_t i = 0; i < ne; ++i)
        {
            for (auto node : b.bound_elements[i].node_id)  //不同类型单元的节点数不一样
            {
                //m_eind保存外部存储的节点编号
                m_eind[m_eind_id] = static_cast<idx_t>(node);  //bound_elements[i].node_id是int类型，这里做了一个安全的类型转换

                //m_eind保存0based的节点编号
                auto it = node_id_2_0based_id.find(node);
                if (it == node_id_2_0based_id.end()) //检查，若该节点不存在于wall_node里，程序出错
                {
                    std::cerr << "Error: Node ID " << node << " not found in wall_nodes!" << std::endl;
                    exit(EXIT_FAILURE);
                }
                //节点存在的情况：
                m_eind_temp[m_eind_id] = it->second;

                m_eind_id++;
            }
        }
    }

    explicit metis_block(const std::vector<element>& b, const vector<Node>& wall_nodes, idx_t k = 2)
    {
        //参数设置
        m_k = k;
        ne = static_cast<idx_t>(b.size());
        nn = wall_nodes.size();

        //构建外部节点编号和由0开始的节点编号的一一映射
        std::unordered_map<idx_t, idx_t> node_id_2_0based_id;
        for (idx_t i = 0; i < nn; ++i)
        {
            node_id_2_0based_id[wall_nodes[i].id] = i;
        }

        //构造m_eptr/m_eind

        //计算m_eptr元素
        m_eptr.resize(ne + 1);
        m_eptr[0] = 0;
        for (idx_t i = 0; i < ne; i++)
        {
            m_eptr[i + 1] = m_eptr[i] + b[i].node_id.size();
        }
        //构建m_eind数组：
        m_eind.resize(m_eptr[ne]);        //m_eptr最后一个元素对应m_eind的大小
        m_eind_temp.resize(m_eptr[ne]);
        idx_t m_eind_id = 0;
        for (idx_t i = 0; i < ne; ++i)
        {
            for (auto node : b[i].node_id)  //不同类型单元的节点数不一样
            {
                //m_eind保存外部存储的节点编号
                m_eind[m_eind_id] = static_cast<idx_t>(node);  //bound_elements[i].node_id是int类型，这里做了一个安全的类型转换

                //m_eind保存0based的节点编号
                auto it = node_id_2_0based_id.find(node);
                if (it == node_id_2_0based_id.end()) //检查，若该节点不存在于wall_node里，程序出错
                {
                    std::cerr << "Error: Node ID " << node << " not found in wall_nodes!" << std::endl;
                    exit(EXIT_FAILURE);
                }
                //节点存在的情况：
                m_eind_temp[m_eind_id] = it->second;

                m_eind_id++;
            }
        }
    }

    void partitionMeshDual(idx_t ncommon = 2)  //使用metis函数进行分区，默认边共享ncommon=2
    {
        assert(ne > 0);  //不允许传入空的单元
        m_epart.resize(ne);
        m_npart.resize(nn);

        // objval：METIS 返回的分区目标值（如切分边数或通信量等）
        idx_t objval = 0;
        int status = METIS_PartMeshDual(
            &ne,           // [输入] 单元数
            &nn,           // [输入] 节点数
            m_eptr.data(),        // [输入] eptr
            m_eind_temp.data(),   // [输入] eind，需要0based的节点编号
            nullptr,         // [输入] vwgt (单元权重)，此处不使用
            nullptr,         // [输入] vsize (单元大小)，此处不使用
            &ncommon,        // [输入] ncommon
            &m_k,            // [输入] nparts
            nullptr,         // [输入] tpwgts (各分区目标权重)，此处不使用
            nullptr,         // [输入] options (高级选项)，此处不使用
            &objval,         // [输出] 分区的目标值
            m_epart.data(),  // [输出] 每个单元的分区号
            m_npart.data()   // [输出] 每个节点的分区号
        );
        

        if (status != METIS_OK) {
            std::cerr << "METIS error: " << status << std::endl;
        }

    }

    
    void generate_mesh_blocks(const boundary& b, const vector<Node>& wall_nodes, vector<mesh_block>& blocks, vector<Node>& unique_boundary_nodes)
    //根据m_k/m_epart/m_npart构建网格块数组
    //同时将网格块之间共享的边界节点唯一地保存在_unique_boundary_nodes
    //传入参数：boundary    
    //         wall_nodes
    //         blocks
    //         unique_boundary_nodes
    //!!!boundary里边界单元数组每个单元索引和m_eptr/m_epart的索引是一一对应的
    {
        idx_t part_id = 0;    //用于循环里临时记录分区号
        int inode_id = 0;     //用于循环里临时记录节点id
        blocks.resize(m_k);    //根据分区数确定网格块数组的大小
        for (idx_t i = 0; i < m_k; ++i)
        {
            blocks[i].block_id = i;   //设定每个网格块的分区号
        }
        //遍历wall，根据METIS的m_epart将单元分配到对应的块中
        for (size_t e = 0; e < b.bound_elements.size(); e++)
        {
            part_id = m_epart[e];
            blocks[part_id].internal_elements.push_back(b.bound_elements[e]);
            //b.bound_elements的索引和m_epart是一一对应的
            //m_epart[e]是索引为e的单元的分区号
        }

        // 首先建立一个从节点编号到其相邻单元分区集合的映射
        // key 为节点编号，value 为该节点所在的所有单元的分区集合
        unordered_map<int, unordered_set<idx_t>> node_partition_map;
        for (size_t e = 0; e < b.bound_elements.size(); e++)
        {
            part_id = m_epart[e];
            for (auto node : b.bound_elements[e].node_id)
            {
                node_partition_map[node].insert(part_id);   //在map的每个节点插入所属单元的分区编号，unordered_set确保这些编号是唯一存在的
            }
        }
        // 为了方便快速查找每个节点的坐标，利用 wall_nodes 建立一个映射：
        // key 为节点编号，value 为对应的 Point<double> 对象（包含坐标）
        unordered_map<int, Node> wall_node_map;
        for (const auto& p : wall_nodes)
        {
            wall_node_map[p.id] = p;
        }
        // 接下来遍历映射，对于每个节点，如果它只出现于一个分区，则视为内部节点；
        // 如果它出现在多个分区（即共享节点），则从 wall_node_map 中获取该节点的完整信息
        for (auto& entry : node_partition_map)
        {
            inode_id = entry.first;   //entry是一个pair<node_id,node_part_id>
            if (entry.second.size() > 1)  // 该点是共享的（出现在多个分区）
            {
                if (wall_node_map.find(inode_id) != wall_node_map.end()) // 检查 wall_node_map 中是否存在该节点
                {
                    const Node& node = wall_node_map[inode_id];
                    // 将该节点插入对应分区的 boundary_nodes 里（可能属于多个分区）
                    for (auto block_id : entry.second)
                    {
                        blocks[block_id].boundary_nodes.push_back(node);
                        blocks[block_id].m_nodes.push_back(node);       //边界节点插入该分区的成员节点中
                    }
                    // 同时将共享节点添加到全局 unique_boundary_nodes 中
                    // 首先检查该节点是否已经存在（这里依据节点 id 判断）
                    bool found = false;
                    for (const auto& p : unique_boundary_nodes)
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
                else   //物面没有该节点，说明程序出现错误
                {
                    std::cerr << "Error: cannot find corresponding node " << inode_id
                        << " in wall_nodes. Partition node and wall nodes do not match!"
                        << std::endl;
                    exit(EXIT_FAILURE);  // 或者 throw std::runtime_error("..."); 根据需求选择中断方式
                }
            }
            else  //该点只属于一个分区，将其插入对应分区的内部节点里
            {
                part_id = *(entry.second.begin());
                if (wall_node_map.find(inode_id) != wall_node_map.end())
                {
                    blocks[part_id].internal_nodes.push_back(wall_node_map[inode_id]);
                    blocks[part_id].m_nodes.push_back(wall_node_map[inode_id]);      //内部节点插入该分区的成员节点中
                }
                    
                else
                {
                    std::cerr << "Error: cannot find corresponding node " << inode_id
                        << " in wall_nodes. Partition node and wall nodes do not match!"
                        << std::endl;
                    exit(EXIT_FAILURE);  // 或者 throw std::runtime_error("..."); 根据需求选择中断方式
                }
            }

        }
    }

    void generate_mesh_blocks(const std::vector<element>& b,const vector<Node>& wall_nodes, vector<mesh_block>& blocks, vector<Node>& unique_boundary_nodes)
        //根据m_k/m_epart/m_npart构建网格块数组
        //同时将网格块之间共享的边界节点唯一地保存在_unique_boundary_nodes
        //传入参数：boundary    
        //         wall_nodes
        //         blocks
        //         unique_boundary_nodes
        //!!!boundary里边界单元数组每个单元索引和m_eptr/m_epart的索引是一一对应的
    {
        idx_t part_id = 0;    //用于循环里临时记录分区号
        int inode_id = 0;     //用于循环里临时记录节点id
        blocks.resize(m_k);    //根据分区数确定网格块数组的大小
        for (idx_t i = 0; i < m_k; ++i)
        {
            blocks[i].block_id = i;   //设定每个网格块的分区号
        }
        //遍历wall，根据METIS的m_epart将单元分配到对应的块中
        for (size_t e = 0; e < b.size(); e++)
        {
            part_id = m_epart[e];
            blocks[part_id].internal_elements.push_back(b[e]);
            //b.bound_elements的索引和m_epart是一一对应的
            //m_epart[e]是索引为e的单元的分区号
        }

        // 首先建立一个从节点编号到其相邻单元分区集合的映射
        // key 为节点编号，value 为该节点所在的所有单元的分区集合
        unordered_map<int, unordered_set<idx_t>> node_partition_map;
        for (size_t e = 0; e < b.size(); e++)
        {
            part_id = m_epart[e];
            for (auto node : b[e].node_id)
            {
                node_partition_map[node].insert(part_id);   //在map的每个节点插入所属单元的分区编号，unordered_set确保这些编号是唯一存在的
            }
        }
        // 为了方便快速查找每个节点的坐标，利用 wall_nodes 建立一个映射：
        // key 为节点编号，value 为对应的 Point<double> 对象（包含坐标）
        unordered_map<int, Node> wall_node_map;
        for (const auto& p : wall_nodes)
        {
            wall_node_map[p.id] = p;
        }
        // 接下来遍历映射，对于每个节点，如果它只出现于一个分区，则视为内部节点；
        // 如果它出现在多个分区（即共享节点），则从 wall_node_map 中获取该节点的完整信息
        for (auto& entry : node_partition_map)
        {
            inode_id = entry.first;   //entry是一个pair<node_id,node_part_id>
            if (entry.second.size() > 1)  // 该点是共享的（出现在多个分区）
            {
                if (wall_node_map.find(inode_id) != wall_node_map.end()) // 检查 wall_node_map 中是否存在该节点
                {
                    const Node& node = wall_node_map[inode_id];
                    // 将该节点插入对应分区的 boundary_nodes 里（可能属于多个分区）
                    for (auto block_id : entry.second)
                    {
                        blocks[block_id].boundary_nodes.push_back(node);
                        blocks[block_id].m_nodes.push_back(node);     //边界节点插入该分区成员节点中
                        blocks[block_id].m_nodes_id.insert(node.id);    //将边界节点的id存入成员节点id查找集
                        blocks[block_id].bndry_nds_id.insert(node.id);    //将边界节点的id存入边界节点id查找集
                    }
                    // 同时将共享节点添加到全局 unique_boundary_nodes 中
                    // 首先检查该节点是否已经存在（这里依据节点 id 判断）
                    bool found = false;
                    for (const auto& p : unique_boundary_nodes)
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
                else   //物面没有该节点，说明程序出现错误
                {
                    std::cerr << "Error: cannot find corresponding node " << inode_id
                        << " in wall_nodes. Partition node and wall nodes do not match!"
                        << std::endl;
                    exit(EXIT_FAILURE);  // 或者 throw std::runtime_error("..."); 根据需求选择中断方式
                }
            }
            else  //该点只属于一个分区，将其插入对应分区的内部节点里
            {
                part_id = *(entry.second.begin());
                if (wall_node_map.find(inode_id) != wall_node_map.end())
                {
                    blocks[part_id].internal_nodes.push_back(wall_node_map[inode_id]);
                    blocks[part_id].m_nodes.push_back(wall_node_map[inode_id]);       //内部节点插入该分区成员节点中
                    blocks[part_id].m_nodes_id.insert(wall_node_map[inode_id].id);    //将内部节点的id存入成员节点id查找集
                }
                    
                else
                {
                    std::cerr << "Error: cannot find corresponding node " << inode_id
                        << " in wall_nodes. Partition node and wall nodes do not match!"
                        << std::endl;
                    exit(EXIT_FAILURE);  // 或者 throw std::runtime_error("..."); 根据需求选择中断方式
                }
            }

        }
    }
};






void _calculate_deformed_coordinates(vector<Node>& every_node)   //计算节点变形后的坐标
    //传入参数：节点数组 every_nodes 
    //该函数会引发传入参数的xyz改变
{
    std::cout << "=====================" << endl;
    std::cout << "Starting to calculate deformed coordinates ... " << std::endl;
    for (int i = 0; i < every_node.size(); ++i)
    {
        //计算每个维度的变形后的值
        every_node[i].point.x += every_node[i].df[0];
        every_node[i].point.y += every_node[i].df[1];
        every_node[i].point.z += every_node[i].df[2];
    }
}

void _write_wall(vector<Node>& wall_nodes,string filename = "wall.dat")
{
    ofstream ofs;
    ofs.open(filename, ios::out | ios::trunc);
    if (!ofs.is_open())
    {
        std::cout << "文件打开失败！" << endl;
        return;
    }
    else
    {
        for (int i = 0; i < wall_nodes.size(); ++i)
        {
            ofs << wall_nodes[i].point.x << " " << wall_nodes[i].point.y << wall_nodes[i].point.z << std::endl;
        }
    }
    ofs.close();
}



//测试网格分区的结果
void _writeBlocksToFile(const std::vector<mesh_block>& blocks, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }

    // 设置固定小数格式，保留6位小数
    outfile << std::fixed << std::setprecision(6);

    // 遍历所有分块
    for (const auto& block : blocks) {
        const idx_t current_block_id = block.block_id;

        // 写入内部节点（每个节点带分块ID）
        for (const auto& node : block.m_nodes) {
            outfile << node.point.x << " "
                << node.point.y << " "
                << node.point.z << " "
                << current_block_id << std::endl;
        }
    }

    outfile.close();
}

void _Reset_blocks(vector<mesh_block>& blocks, const vector<Node>& wall_nodes, const unordered_map<int, int>& id2node)
//根据新的变形结果重新设值blocks内部每个节点的坐标
//传入参数： 分块网格数组 blocks
//          物面节点数组 wall_nodes (需要使用变形后的节点构建)
{
    ////构建wall_nodes的id2node查找集  其中node是节点在wall_nodes的索引！！！
    //id2node.reserve(wall_nodes.size());
    //for (int i = 0; i < wall_nodes.size(); ++i)
    //{
    //    id2node[wall_nodes[i].id] = i;
    //}

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

// multi_partition 类：负责多级 METIS 分区及分组 RBF 算法
class multi_partition {
public:
    // 构造：传入每层分区数，例如 {4,3,2,2}
    explicit multi_partition(const std::vector<idx_t>& parts_per_level)
        : levels(parts_per_level.size()), parts_per_level(parts_per_level) {}

    //根据外部传入的_nd2wall_lvl_mp设置m_nd2wall_lvl_mp
    void set_m_nd2wall_lvl_mp(const std::unordered_map<int, int>& nd2wall_lvl)
    {
        m_nd2wall_lvl_mp = nd2wall_lvl;
    }

    // 执行多级分区：
    //  - elems: 全局边界单元列表
    //  - wall_nodes: 全局节点列表
    // 分区结果保存在 blocks_per_level[level] 中，每个元素为 mesh_block
    // 唯一边界节点保存在 unique_boundary_per_level[level][block_index]
    void divide_wall(const std::vector<element>& elems,
        const std::vector<Node>& wall_nodes) {
        // 预分配各层容器
        blocks_per_level.clear();
        unique_boundary_per_level.clear();
        blocks_per_level.reserve(levels);
        unique_boundary_per_level.reserve(levels);
        unique_bndry2blkid_per_lvl.resize(levels);

        // 根分区：一个 mesh_block 包含所有元素与节点
        std::vector<mesh_block> current;
        current.reserve(1);
        mesh_block root;
        root.internal_elements = elems;
        root.m_nodes = wall_nodes;
        current.push_back(std::move(root));

        // 复用容器，减少分配
        std::vector<mesh_block> next_blocks;
        std::vector<std::vector<Node>> level_unique_boundary;

        // 迭代每个分区层级
        for (size_t lvl = 0; lvl < levels; ++lvl) {
            next_blocks.clear();
            level_unique_boundary.clear();
            // 粗略估计下一级块数并预分配
            size_t est = current.size() * parts_per_level[lvl];
            next_blocks.reserve(est);
            level_unique_boundary.reserve(current.size());

            // 对上一级的每个块分区
            for (auto& blk : current) {
                // 构造 metis_block
                metis_block mblk(blk.internal_elements, blk.m_nodes, parts_per_level[lvl]);
                mblk.partitionMeshDual();
                // 调用 generate_mesh_blocks，获取子 mesh_block 和它们的 unique boundary
                std::vector<mesh_block> children;
                std::vector<Node> unique_bndry;
                children.reserve(parts_per_level[lvl]);
                mblk.generate_mesh_blocks(blk.internal_elements,
                    blk.m_nodes,
                    children,
                    unique_bndry);


                // 除了第一次整体的分区以外，后续的分区在 generate_mesh_blocks 之后，
                // 其分出的新子分区有可能失去一部分边界节点。
                // 原因：generate_mesh_blocks 接收的只是当前大分区的信息，
                //       它并不知道该大分区和其他大分区之间共享的边界节点。
                // 因此这里需要额外的补充操作：
                //
                // 遍历父分区 (blk) 的 boundary_nodes，检查哪些节点也属于子分区 (iblk) 的 m_nodes：
                //   - 如果是父分区的边界节点，同时存在于子分区的成员节点中，
                //   - 且该节点还没被标记为子分区的边界节点，
                //   - 则将其补充进子分区的 boundary_nodes 和 bndry_nds_id。
                //
                // 这样，每个子分区的边界节点集合会被补齐，包含：
                //   1. 子分区与兄弟子分区共享的节点
                //   2. 子分区与父分区外部分区共享的节点
                //
                // 效果：保证了跨层级的 unique_bndry 节点不会丢失，
                //      后续构建 unique_bndry2blkid_per_lvl 时，映射表就是完整的。
                if (lvl != 0)
                {
                    for (auto& iblk : children)
                    {
                        for (auto& ind : blk.boundary_nodes)    //blk是上一次分区的结果，每个大分区的边界节点包含了其和其他大分区共享的边界节点
                        {
                            if (iblk.m_nodes_id.find(ind.id) != iblk.m_nodes_id.end())   //ind是大分区的边界节点，如果他属于小分区的成员节点
                            {
                                if (iblk.bndry_nds_id.find(ind.id) == iblk.bndry_nds_id.end())   //且不属于小分区的边界节点，则插入该节点
                                {
                                    iblk.boundary_nodes.push_back(ind);   //插入边界节点
                                    iblk.bndry_nds_id.insert(ind.id);     //插入边界节点的id
                                }
                            }
                        }
                    }
                }
                

                // 移动子分区及其唯一边界
                next_blocks.insert(
                    next_blocks.end(),
                    std::make_move_iterator(children.begin()),
                    std::make_move_iterator(children.end()));
                level_unique_boundary.push_back(std::move(unique_bndry));
            }

            // 保存本层结果
            blocks_per_level.push_back(std::move(next_blocks));
            set_blk_id(blocks_per_level[lvl]);    //重新调整分区编号
            unique_boundary_per_level.push_back(std::move(level_unique_boundary));

            build_unique_bndry2blkid_for_lvl(lvl);    //构建unique_bndry_nodes的id和分区编号的查找集

            // 为下一层准备：交换以复用内存，避免拷贝
            current.swap(blocks_per_level.back());
        }

    }

    // 分组 RBF 算法主入口：
    // tol_steps: 每层的误差容限数组，长度 == levels
    // wall_nodes: 全局物面节点数组，包含最终的已知变形结果
    // all_nodes_coords: 全局所有待变形节点坐标
    void multi_partition_rbf_algorithm(
        const std::vector<double>& tol_steps,
        std::vector<Node>& wall_nodes,
        std::vector<Node>& all_nodes_coords)
    {
        assert(tol_steps.size() == levels);    //检查误差设置的级数与分区的级数是否一致
        set_wall_id2nodes(wall_id2nodes, wall_nodes);    //构建物面的索引和id的查找集
        std::vector<Node> wall_accurate = set_accurate_wall(wall_nodes);   //根据物面节点初始位置和待变形量计算准确的变形结果
        std::vector<Node> wall_prev = wall_nodes;        //当前物面结果
        std::vector<Node> wall_0 = wall_nodes;           //上一步物面结果
        block_rbf_per_level.clear(); 
        block_rbf_per_level.resize(levels);

        for (size_t lvl = 0; lvl < levels; ++lvl) {
            double tol = tol_steps[lvl];

            // 步骤 0: 设置当前步下的wall_nodes(wall_prev)，wall_prev的节点坐标和待变形量需要更新
            if (lvl != 0)    //lvl==0的时候是整个物面系统进行第一次变形，此时物面节点就应该保持初始状态，lvl!=0开始更新物面节点
            {
                auto t_pre_start = clock();

                update_wall_prev(wall_prev, wall_accurate, all_nodes_coords);             //根据更新后的all_node_coord和准确的物面变形结果更新wall_prev的节点坐标和待变形量
                _Reset_blocks(blocks_per_level[lvl], wall_prev, wall_id2nodes);    //根据更新后的wall_prev重新设置分区内部的节点坐标, 变形量由后续的set_internal_nodes_df处理
                update_unique_bndry_nodes(wall_prev, lvl);                        //更新已经用过的unique_bndry的位置，这些点已经准确变形，无需再变
                set_blocks_df_and_tree(wall_prev, lvl);                            //根据更新后的wall_prev重新设置分区内部的节点待变形量，并且计算了分区的block_D!!!
                set_unique_bndry_tree(lvl, unique_bndry_tree, all_nodes_coords);   //构建unique_bndry的二叉树

                auto t_pre_end = clock();
                std::cout << "[lvl = " << lvl << "] Preprocessing time: " << double(t_pre_end - t_pre_start) << " ms" << std::endl;
            }

            // 步骤 1: 构建 RBF 插值系统（第0层只有1个系统，从第1层开始每个块一个）
            auto t_rbf_start = clock();

            if (lvl == 0) {
                build_single_rbf_system(wall_nodes, tol, lvl);
            }
            else{    //构建分区RBF系统
                build_multiple_rbf_systems(wall_prev, wall_nodes, tol, lvl);
            }

            auto t_rbf_end = clock();
            std::cout << "[lvl = " << lvl << "] rbf coeff calculating time: " << double(t_rbf_end - t_rbf_start) << " ms" << std::endl;


            
            // 步骤 2: 执行节点变形计算
            auto t_df_start = clock();

            if (lvl == 0) {
                apply_simple_rbf_deformation(all_nodes_coords, lvl);
            }
            else {
                apply_drrbf_deformation(all_nodes_coords, lvl);
            }

            auto t_df_end = clock();
            std::cout << "[lvl = " << lvl << "] df calculating time: " << double(t_df_end - t_df_start) << " ms" << std::endl;


            auto t_deform_start = clock();
            _calculate_deformed_coordinates(all_nodes_coords);   //计算变形后的节点坐标
            auto t_deform_end = clock();
            std::cout << "[lvl = " << lvl << "] deforming nodes time: " << double(t_deform_end - t_deform_start) << " ms" << std::endl;
            std::cout << std::endl;

            //test_temp(all_nodes_coords, lvl);
            //test_Wall(all_nodes_coords, lvl);
       
        }
    }

    //测试当前级数的结果
    void test_temp(const std::vector<Node>& all_coords, int level)
    {
        _node_coords = all_coords;
        std::string writefile_name = "test_su2_deformed_lvl" + std::to_string(level) + "_naca0012.su2";
        writefile(writefile_name);
    }

    void test_Wall(const std::vector<Node>& all_coords, int level)
    {
        vector<Node> wall_temp;
        wall_temp.resize(0);
        wall_temp = Set_wall_nodes(_every_boundary, all_coords);
        std::string writefile_name = "test_wall_lvl" + std::to_string(level) + ".dat";
        _write_wall(wall_temp,writefile_name);
    }

    
    void set_blk_id(std::vector<mesh_block>& blk)
        //由于generate_mesh_blocks只能根据其接受的分区设置新分出的小分区的id
        //此函数将局部的id转换成全局的id
        //由于只对一个物面进行分区，故按顺序分配id即可
        //blk为blocks_per_level[lvl]的mesh_block的数组
    {
        for (int i = 0; i < blk.size(); ++i)
        {
            blk[i].block_id = i;
        }
    }

    // 获取某一级的分区块
    std::vector<mesh_block>& get_level_blocks(size_t level) {
        return blocks_per_level.at(level);
    }
    // 获取某一级的分区块（只读）
    const std::vector<mesh_block>& get_level_blocks(size_t level) const{
        return blocks_per_level.at(level);
    }


    // 获取某一级中第block_index块的唯一边界节点列表
    std::vector<Node>& get_unique_boundary(size_t level, size_t block_index) {
        return unique_boundary_per_level.at(level).at(block_index);
    }

    // 获取某一级中第block_index块的唯一边界节点列表 （只读）
    const std::vector<Node>& get_unique_boundary (size_t level, size_t block_index) const{
        return unique_boundary_per_level.at(level).at(block_index);
    }

    // 为 multi_partition 添加查询总层数的接口
    size_t get_total_levels() const {
        return levels;
    }

    // 查询接口：获取某一级别某节点所在的所有分区编号
    const std::vector<int>& query_unique_bndry2blkid(size_t lvl, int nd_id) const
        //输入是当前层级的lvl，对应的查找集放在lvl - 1的位置
        //仔细检查unique_bndry2blkid_per_lvl建立的逻辑即可查清
    {
        auto it = unique_bndry2blkid_per_lvl.at(lvl - 1).find(nd_id);   //查找节点编号对应的分区编号
        if (it == unique_bndry2blkid_per_lvl[lvl - 1].end())    //如果没有找到传入的节点
        {
            //输出警告信息
            std::cerr << "[Error] unique_bndry node " << nd_id
                << " has no corresponding block at level " << lvl
                << ". Please check partitioning.\n";
            std::abort();  //中止程序
        }
        return it->second;   //返回分区编号数组
    }



private:
    size_t levels;    // 分区层数
    std::unordered_map<int, int> wall_id2nodes;   //物面的索引和id查找集
    std::vector<idx_t> parts_per_level;  // 每层的分区数
    // blocks_per_level[lvl] = vector<mesh_block> 当前层所有的 mesh_block
    std::vector<std::vector<mesh_block>> blocks_per_level;
    // unique_boundary_per_level[lvl][i] = 第lvl层第i个块的唯一边界节点
    std::vector<std::vector<std::vector<Node>>> unique_boundary_per_level;
    std::vector<std::vector<deform_calculate>> block_rbf_per_level;   // 每层每块对应的 RBF 插值系统

    std::unordered_set<int> global_unique_bndry_set;     //已经是unique_bndry在all_noord_coords的索引的查找集
    std::vector<int> global_unique_bndry_id;             //unique_bndry在all_noord_coords的索引
    GridBTree<int, Point<double>> unique_bndry_tree;    //当前变形级数前所有的使用过的Unique_bndry_nodes
    std::unordered_map<int, int> m_nd2wall_lvl_mp;      //节点等级映射（成员）
    std::vector<std::unordered_map<int, std::vector<int>>> unique_bndry2blkid_per_lvl;  //每个层级的unique_bndry_nodes对应的分区编号

    //构建物面的索引和id查找集
    void set_wall_id2nodes(std::unordered_map<int, int>& wall_id2nd, const std::vector<Node>& wall_nodes)
    {
        wall_id2nd.reserve(wall_nodes.size());
        for (int i = 0; i < wall_nodes.size(); ++i)
        {
            wall_id2nd[wall_nodes[i].id] = i;
        }
    }

    //计算准确的物面变形后的物面节点
    std::vector<Node> set_accurate_wall(const std::vector<Node>& wall_0)
        //wall_0：外部传入的_wall_nodes，需有节点原始位置和待变形量
    {
        std::vector<Node> wall_temp;
        wall_temp.resize(wall_0.size());
        int point_size = wall_0[0].point.size();      //计算点的维度
        for (int i = 0; i < wall_0.size(); ++i)       //遍历物面数组
        {
            for (int k = 0; k < point_size; ++k)
            {
                wall_temp[i].point[k] = wall_0[i].point[k] + wall_0[i].df[k];   //变形后的节点位置为原来的位置加上设定的位移量
                wall_temp[i].id = wall_0[i].id;           //保留id编号
            }
        }
        return wall_temp;
    }

    //根据更新后的all_nodes_coords修改wall_prev的节点坐标和待变形量
    //all_node_coord 是按照节点id排序的全局节点数组，每次循环的RBF变形后，其节点坐标和df都发生改变
    //wall_prev是当前的物面节点数组
    //wall_0是上一步的物面节点，用以提供上一步的待变形量
    void update_wall_prev(std::vector<Node>& wall_prev, const std::vector<Node>& wall_accurate, const std::vector<Node>& all_node_coord)
    {
        int df_size = wall_prev[0].df.size();                  //计算df的维度
        for (auto& inode : wall_prev)
        {
            assert(inode.id == all_node_coord[inode.id].id);   //确保在all_node_coord里索引的节点和当前的物面节点一致，节点编号是唯一的
            inode = all_node_coord[inode.id];
            auto it = wall_id2nodes.find(inode.id);
            assert(it != wall_id2nodes.end());                //确保没有在wall里查找wall以外的节点
            for (int i = 0; i < df_size; ++i)
            {
                inode.df[i] = wall_accurate[it->second].point[i] - inode.point[i];    //新的待变形量为准确的节点位置减去已变形的节点位置：s*-sk
            }
        }
    }

    //根据更新后的wall_prev更新unique_boundary_nodes内部的节点坐标和待变形量
    //wall_prev: 更新后的物面节点
    //lvl: 级数
    void update_unique_bndry_nodes(const std::vector<Node>& wall_prev, int lvl)
    {
        for (auto& iboundary : unique_boundary_per_level[lvl])    //iboundary 每个分区的共享边界点数组
        {
            for (auto& inode : iboundary)      //inode 共享边界点
            {
                auto it = wall_id2nodes.find(inode.id);
                assert(it != wall_id2nodes.end());               //确保没有在wall里查找wall以外的节点
                inode = wall_prev[it->second];                   //将wall_prev对应的节点的坐标和df传给unique_bndry_nodes
            }
        }
    }

    //计算当前级数下的blocks的节点的df并且构建相应的二叉树
    //wall_prev: 更新后的物面节点
    //lvl: 级数
    void set_blocks_df_and_tree(const std::vector<Node>& wall_prev, int lvl)
    {
        for (auto& iblock : blocks_per_level[lvl])    //iblock：每个小分区
        {
            iblock.set_blk_nodes_df(wall_prev, wall_id2nodes);      //根据已经更新的wall_prev设定blocks的节点待变形量
            iblock.set_block_tree();
        }
    }

    // 第0层构建单个RBF系统
    // wall_final: 最終的物面变形结果
    // tol: 误差容限
    void build_single_rbf_system(const std::vector<Node>& wall_final,
        double tol,
        size_t lvl)
    {
        deform_calculate single_rbf;
        single_rbf.Greedy_algorithm(unique_boundary_per_level[lvl][0], tol);
        block_rbf_per_level[lvl].resize(1);
        block_rbf_per_level[lvl][0] = std::move(single_rbf);
    }

    // 第1层及以上构建多个RBF系统
    void build_multiple_rbf_systems(const std::vector<Node>& wall_prev,
        const std::vector<Node>& wall_final,
        double tol,
        size_t lvl)
    {
        auto& blocks = blocks_per_level[lvl];
        auto& unique_bndry = unique_boundary_per_level[lvl];
        auto& rbfs = block_rbf_per_level[lvl];
        rbfs.resize(blocks.size());

        

        if (lvl < levels - 1)   //非最后一步使用共享边界点构建插值系统
        {
            set_block_rbf(rbfs, unique_bndry, blocks);
        }
        else  //最后一步使用内部点构建插值系统
        {
            set_block_rbf(rbfs, blocks);
        }

        for (size_t i = 0; i < blocks.size(); ++i)
        {
            rbfs[i].Greedy_algorithm(tol);    //使用待支撑点计算系数
        }
    }

    void set_block_rbf(vector<deform_calculate>& block_rbf,  vector<vector<Node>> unique_bndry, const vector<mesh_block>& blocks)
    //使用unique_bndry的点构造每个分区的插值系统
    //插值系统的支撑点由unique_bnry的点构成，它包括：
    //1.该网格块的unique_bndry点
    //2.该网格块的block_D形成的圆覆盖的网格块的点
    {
        double d_temp = 0.0;     //block_D的临时变量
        double d_r2omega_temp = 0.0;    //到网格快的距离的临时变量
        int key_temp = 0;

        for (int i = 0; i < block_rbf.size(); ++i)
        {
            unordered_set<int> nodes_id_set;     //加入插值系统的待支撑点(external points)的查找集，用于确保点是唯一的
            d_temp = blocks[i].block_D;          //设置当前block[i]的覆盖范围
            for (int j = 0; j < blocks.size(); ++j)    //计算其他区间的unique_bndry点到block[i]的距离
            {
                if (j == i)  //该网格块本身，直接插入其内部的unique_bndry
                {
                    for (auto& inode : unique_bndry[i])
                    {
                        block_rbf[i].external_suppoints.push_back(inode);  //向RBF系统插入该网格快的unique_bndry点
                        nodes_id_set.insert(inode.id);   //将已插入的点放入查找集里
                    }
                } 
                else   //非该网格块，需检查每个网格块的边界点和内部点是否被该网格块的block_D形成的圆覆盖
                {
                    for (auto& inode : blocks[j].internal_nodes)  //检查内部点
                    {
                        blocks[i].block_tree.search(d_r2omega_temp, key_temp, inode.point); //查找当前点到blocks[i]的最小距离
                        if (d_r2omega_temp <= d_temp)   //这个点被该网格块的block_D形成的圆覆盖，则插入
                        {
                            if (nodes_id_set.find(inode.id) == nodes_id_set.end())  //检查该点是否在查找集内
                            {
                                block_rbf[i].external_suppoints.push_back(inode);
                                block_rbf[i].external_suppoints.back().df.setZero();  //非内部点的df都为0
                                nodes_id_set.insert(inode.id);    //将已插入的点放入查找集里
                            }
                        }
                    }
                    for (auto& inode : blocks[j].boundary_nodes)  //检查边界点
                    {
                        blocks[i].block_tree.search(d_r2omega_temp, key_temp, inode.point); //查找当前点到blocks[i]的最小距离
                        if (d_r2omega_temp <= d_temp)   //这个点被该网格块的block_D形成的圆覆盖，则插入
                        {
                            if (nodes_id_set.find(inode.id) == nodes_id_set.end())  //检查该点是否在查找集内
                            {
                                block_rbf[i].external_suppoints.push_back(inode);
                                block_rbf[i].external_suppoints.back().df.setZero();  //非内部点的df都为0
                                nodes_id_set.insert(inode.id);    //将已插入的点放入查找集里
                            }
                        }
                    }
                }
                //for (int k = 0; k < unique_bndry[j].size(); ++k)    //循环每个分区的unique_bndry
                //{
                //    auto& inode = unique_bndry[j][k];    //inode是当前点
                //    if (j == i)    //该网格块本身，直接插入其内部的unique_bndry
                //    {
                //        if (nodes_id_set.find(inode.id) == nodes_id_set.end())  //该点不在查找集内，可插入
                //        {
                //            block_rbf[i].external_suppoints.push_back(inode);
                //            nodes_id_set.insert(inode.id);    //将已插入的点放入查找集里
                //        }
                //    }
                //    else
                //    {
                //        blocks[i].block_tree.search(d_r2omega_temp, key_temp, inode.point);   //查找当前点到blocks[i]的最小距离
                //        if (d_r2omega_temp <= d_temp)   //这个点被该网格块的block_D形成的圆覆盖，则插入
                //        {
                //            if (nodes_id_set.find(inode.id) == nodes_id_set.end())  //检查该点是否在查找集内
                //            {
                //                block_rbf[i].external_suppoints.push_back(inode);
                //                block_rbf[i].external_suppoints.back().df.setZero();  //非内部点的df都为0
                //                nodes_id_set.insert(inode.id);    //将已插入的点放入查找集里
                //            }
                //        }
                //    }
                //}
            }
        }
    }
    void set_block_rbf(vector<deform_calculate>& block_rbf, const vector<mesh_block>& blocks)
     //使用网格块内部的点构造每个分区的插值系统（分级的最后一步）
     //插值系统的支撑点由网格块内部的点构成，它包括：
     //1.该网格块的内部点
     //2.该网格块的block_D形成的圆覆盖的网格块的点
    {
        double d_temp = 0.0;     //block_D的临时变量
        double d_r2omega_temp = 0.0;    //到网格快的距离的临时变量
        int key_temp = 0;

        for (int i = 0; i < block_rbf.size(); ++i)
        {
            unordered_set<int> nodes_id_set;     //加入插值系统的待支撑点(external points)的查找集，用于确保点是唯一的
            d_temp = blocks[i].block_D;          //设置当前block[i]的覆盖范围
            for (int j = 0; j < blocks.size(); ++j)    //计算其他区间的unique_bndry点到block[i]的距离
            {
                for (int k = 0; k < blocks[j].internal_nodes.size(); ++k)  //循环每个分区的内部点
                {
                    auto& inode = blocks[j].internal_nodes[k];    //inode是当前点
                    if (j == i)    //该网格块本身，直接插入其内部的点
                    {
                        if (nodes_id_set.find(inode.id) == nodes_id_set.end())  //该点不在查找集内，可插入
                        {
                            block_rbf[i].external_suppoints.push_back(inode);
                            nodes_id_set.insert(inode.id);    //将已插入的点放入查找集里
                        }
                    }
                    else
                    {
                        blocks[i].block_tree.search(d_r2omega_temp, key_temp, inode.point);  //查找当前区间的内部点到变形区间的最小距离
                        if (d_r2omega_temp <= d_temp)  //这个点被该网格块的block_D形成的圆覆盖，则插入
                        {
                            if (nodes_id_set.find(inode.id) == nodes_id_set.end())
                            {
                                block_rbf[i].external_suppoints.push_back(inode);
                                block_rbf[i].external_suppoints.back().df.setZero();   //非内部点的df都为0
                                nodes_id_set.insert(inode.id);    //将已插入的点放入查找集里
                            }
                        }
                    }
                }
                for (int k = 0; k < blocks[j].boundary_nodes.size(); ++k)  //循环每个分区的边界点
                {
                    auto& inode = blocks[j].boundary_nodes[k];  //inode是当前点
                    blocks[i].block_tree.search(d_r2omega_temp, key_temp, inode.point);  //查找当前区间的边界点到变形区间的最小距离
                    if (d_r2omega_temp <= d_temp)    //这个点被该网格块的block_D形成的圆覆盖，则插入
                    {
                        if(nodes_id_set.find(inode.id) == nodes_id_set.end())   //查找该点是否在查找集
                        {
                            block_rbf[i].external_suppoints.push_back(inode);
                            block_rbf[i].external_suppoints.back().df.setZero();   //非内部点的df都为0
                            nodes_id_set.insert(inode.id);    //将已插入的点放入查找集里
                        }
                    }
                }
            }
        }
    }


    // 简单 RBF 插值变形
    void apply_simple_rbf_deformation(std::vector<Node>& coords,
        size_t lvl)
    //第一级变形，只需使用简单的RBF变形计算即可
    {
        for (auto& node : coords) {
            block_rbf_per_level[lvl][0].calculate__deform_RBF(node);    //第一次变形使用简单的RBF计算
        }
    }

    // 分组 DRRBF 插值变形
    void apply_drrbf_deformation(std::vector<Node>& all_noords, 
        size_t lvl)
    //第二级以上的变形，使用DDRBF变形
    //需计算每个点到动边界的距离和到静止边界的距离
    //静止距离的计算应该包括每个分区内部的点以及整体的unqiue_boundry点
    {
        const auto& blocks = blocks_per_level[lvl];
        auto& rbfs = block_rbf_per_level[lvl];

        for (auto& inode : all_noords)
        {
            // --- 新增等级筛选 ---
            auto it = m_nd2wall_lvl_mp.find(inode.id);
            if (it != m_nd2wall_lvl_mp.end())
            {
                int node_lvl = it->second;
                if (node_lvl > 3 - static_cast<int>(lvl))
                    continue; // 等级大于 3-lvl，不处理
            }
            int i_id = 0;     //动边界的分区编号
            int key_temp = 0;
            double i_d_r2omega1 = 0;   //到动边界的最小距离
            double i_d_r2omega2 = 1.0e10;   //到其他静止边界分区的最小距离
            double i_d_r2omega2_temp = 0;   //临时距离
            

            //查找到动边界和静边界的距离
            find_moving_and_static_bndry(inode, blocks, lvl, i_id, i_d_r2omega1, i_d_r2omega2);



            if (global_unique_bndry_set.find(inode.id) == global_unique_bndry_set.end())  //unique_bndry点无需变形
            {
                rbfs[i_id].calculate_deform_RRBF(inode, i_d_r2omega1, i_d_r2omega2, blocks[i_id].block_D);  //_unique_boundary_nodes以外的点使用分组变形
            }
            else
            {
                inode.df.setZero();   //已经用过的unique_bndry的节点变形量为0
            }

        }
    }


    void set_unique_bndry_tree(size_t lvl, GridBTree<int, Point<double>>& unique_bndry_tree_temp, std::vector<Node> all_noords)
        //将之前的unique_bndry点构建成一个二叉树用以查找到unique_bndry点的最小距离
        //lvl：当前变形的级数
        //unique_bndry_tree_temp: 用unique_bndry点构建的二叉树
    {
        collect_unique_bndry_nodes(lvl, global_unique_bndry_set, global_unique_bndry_id);   //收集上一级的unique_bndry节点
        unique_bndry_tree_temp.setGridSize(global_unique_bndry_id.size());    //重新调整二叉树的大小
        for (int i = 0; i < global_unique_bndry_id.size(); ++i)    //遍历global_unique_bndry_id的所有索引
        {
            auto& id = global_unique_bndry_id[i];
            assert(id == all_noords[id].id);   //确保点的id和它的索引是一致的
            unique_bndry_tree_temp.setGrid(i, id, all_noords[id].point);  //将点放入二叉树内
        }
        unique_bndry_tree_temp.construct();
    }

    void collect_unique_bndry_nodes(size_t lvl, std::unordered_set<int>& bndry_set, std::vector<int>& bndry_id)
        //将之前所有分区存在的unique_bndry节点所对应的索引存放到bndry_id
        //bndry_set确保点是唯一的
        //只需要不断地将上一级的存入即可
        //该函数确保了drrbf查找到静边界距离的时候，不会漏掉unique_bndry的点
        //lvl：变形当前的级数，其上一级是lvl-1
    {
        for (const auto& blk_bndry : unique_boundary_per_level[lvl - 1])
        {
            for (const auto& inode : blk_bndry)
            {
                if (bndry_set.insert(inode.id).second)  //如果节点未被重复，则插入
                {
                    bndry_id.push_back(inode.id);
                }
            }
        }
    }

    void build_unique_bndry2blkid_for_lvl(size_t lvl)
        //建立每个lvl对应的unique_bndry_nodes和blk_id的查找集
        //输入: lvl: 层级
    {
        auto& mp = unique_bndry2blkid_per_lvl[lvl];  //取当前层级的unique_bndry_nodes和blk_id的查找集
        mp.clear();

        //遍历本层的每个分区
        for (const auto& blk : blocks_per_level[lvl])
        {
            //遍历该分区的每个unique_bndry_id
            for (int nd_id : blk.bndry_nds_id)
            {
                mp[nd_id].push_back(blk.block_id);   //插入节点id和分区id的对应关系
            }
        }
    }

    void find_moving_and_static_bndry(
        const Node& ind,
        const std::vector<mesh_block>& blocks,
        int lvl,
        int& i_id,
        double& i_d_r2omega1,
        double& i_d_r2omega2)const
        //查找距离当前空间节点最近的动边界和静止边界
        //算法思路：首先查找到unique_bndry_nodes集合的最小距离d1和对应节点nd1
        //         查找空间节点到nd1所属的分区编号所对应的所有分区的距离di
        //         将d1和di里最小的距离作为动边界，确定对应分区编号，第二小的距离作为静止边界
        //输入：空间节点 ind
        //     分区数组 blocks
        //     层级 lvl
        //     动边界编号 i_id
        //     到动边界距离 i_d_r2omega1
        //     到静边界距离 i_d_r2omega2
    {
        int nearest_id2unique = -1;
        double min_dist2unique = 0.0;

        //查找最近的unique_bndry_nodes
        unique_bndry_tree.search(min_dist2unique, nearest_id2unique, ind.point);
        
        //查找分区集合
        const auto& candidate_blks = query_unique_bndry2blkid(lvl, nearest_id2unique);

        //计算到这些分区的距离
        std::vector<std::pair<double, int>> dist_blk_pairs;    //距离和分区编号构成的配对数组
        dist_blk_pairs.reserve(candidate_blks.size());
        
        int key_temp = 0;
        for (int blk_id : candidate_blks)
        {
            double d_temp = 0.0;
            blocks[blk_id].block_tree.search(d_temp, key_temp, ind.point);
            dist_blk_pairs.emplace_back(d_temp, blk_id);
        }
        
        
        //排序
        std::sort(dist_blk_pairs.begin(), dist_blk_pairs.end(), [](auto& a, auto& b) {return a.first < b.first; });

        //输出结果
        if (!dist_blk_pairs.empty()) {
            i_d_r2omega1 = dist_blk_pairs[0].first;   // 动边界距离
            i_id = dist_blk_pairs[0].second;  // 动边界分区编号

            // 默认静边界距离 = 最近 unique_bndry 距离
            i_d_r2omega2 = min_dist2unique;

            // 如果存在第二小的候选分区距离，则取两者较小值
            if (dist_blk_pairs.size() > 1) {
                i_d_r2omega2 = std::min(i_d_r2omega2, dist_blk_pairs[1].first);
            }
        }
        else {
            std::cerr << "[Error] no candidate blocks found for unique_bndry "
                << nearest_id2unique << " at level " << lvl << "\n";
            std::abort();
        }
    }


};





class TEST {
public:
    // 输出 multi_partition 的所有 block 到指定目录
        // 输出 multi_partition 的所有 block 及唯一边界节点 到指定目录
    static void writeMultiLevelBlocksToDirectory(const multi_partition& mp,
        const std::string& directory_path) {
        // 创建目录（若不存在）
        if (!fs::exists(directory_path)) {
            fs::create_directories(directory_path);
        }

        size_t levels = mp.get_total_levels();
        // 1. 输出每层的 blocks
        for (size_t lvl = 0; lvl < levels; ++lvl) {
            const auto& blocks = mp.get_level_blocks(lvl);
            std::string filename = directory_path + "/level_" + std::to_string(lvl) + ".txt";
            _writeBlocksToFile(blocks, filename);
        }

        // 2. 从 level 0 开始输出每层的 unique boundary nodes
        for (size_t lvl = 0; lvl < levels; ++lvl) {
            std::string filename = directory_path + "/unique_level_" + std::to_string(lvl) + ".txt";
            std::ofstream outfile(filename);
            if (!outfile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << " for writing unique boundary nodes." << std::endl;
                continue;
            }
            outfile << std::fixed << std::setprecision(6);
            size_t num_blocks = mp.get_level_blocks(lvl).size();
            // 遍历本层每个块的唯一边界节点
            for (size_t b = 0; b < num_blocks; ++b) {
                const auto& ub = mp.get_unique_boundary(lvl, b);
                for (const auto& node : ub) {
                    outfile << node.point.x << " "
                        << node.point.y << " "
                        << node.point.z << " "
                        << b << std::endl;
                }
            }
            outfile.close();
        }
    }
};


void _set_wall_tree(GridBTree<int, Point<double>>& global_wall_tree ,const std::vector<Node>& wall)
//构建物面节点的二叉树，用以查询空间节点到物面节点的最小距离
//输入：global_wall_tree：全局的物面二叉树
//      wall：全局的物面节点数组
{
    global_wall_tree.setGridSize(wall.size());    //设置全局的物面二叉树的大小
    for (int i = 0; i < wall.size(); ++i)
    {
        global_wall_tree.setGrid(i, wall[i].id, wall[i].point);  //构建二叉树，参数分别为索引、点的id和点
    }
    global_wall_tree.construct();
}

void _set_nd2wall_lvl_mp(std::unordered_map<int, int>& _mp,const std::vector<Node>& all_nds,const GridBTree<int, Point<double>>& global_wall_tree, double D)
//根据节点到物面的距离给节点分等级，其中距离物面D的设为1，距离物面5D的设为2，其他设为3
//输入：_mp：空间节点的id和对应等级的映射
//      all_nds：空间节点坐标
//      global_wall_tree：物面节点的二叉树
//      D：D是五倍最大的位移变形量
{
    int key_temp = 0;       //最近邻节点的id
    double min_d = 0.0;     //最近距离
    for (const auto& ind : all_nds)
    {
        //查找ind到物面的最近距离
        global_wall_tree.search(min_d, key_temp, ind.point);
        int lvl = 3;   //默认等级为3
        if (min_d <= 5.0 * D)
        {
            lvl = 2;
        }
        if (min_d <= D)
        {
            lvl = 1;
        }

        _mp[ind.id] = lvl;  //把节点 id -> 等级存入 map
    }
}

GridBTree<int, Point<double>> _wall_tree;   //初始化物面节点的二叉树
std::unordered_map<int, int> _nd2wall_lvl_mp;  //空间节点的id和对应等级的映射

int main() {
    
    ////O型结果
    //string readfile_name = "naca0012.su2";
    //string writefile_name = "su2_deformed_naca0012.su2";

    //C型结果
    //string readfile_name = "naca0012_Cmesh_test.su2";
    //string writefile_name = "naca0012_Cmesh_test_result.su2";

    string readfile_name = "../data/naca0012_Cmesh.su2";
    string writefile_name = "../output/naca0012_Cmesh_result.su2";

    //double greedy_tol_0 = 1e-14;  //设置贪心算法误差
    //double greedy_tol_1 = 1e-8;  //设置贪心算法误差
    //idx_t metis_k = 16;  //metis分区数
    //std::vector<deform_calculate> block_RBF;      //分区的RBF系统数组

    readfile(readfile_name);  //读取网格文件
    _wall_nodes = Set_wall_nodes(_every_boundary, _node_coords);   //设置物面节点数组
    _select_R(_node_coords, _R);           //选出网格点间近似的最大距离
    calculat_wall_deformation(_wall_nodes,_D);                        //计算物面节点位移
   
    _set_wall_tree(_wall_tree, _wall_nodes);      //构建物面节点的二叉树
    _set_nd2wall_lvl_mp(_nd2wall_lvl_mp, _node_coords, _wall_tree, _D);   //构建空间节点id和对应等级的映射



    ////***************全局RBF算法开始**********************
    //auto t_global_start = clock();

    //auto t_global_rbf_start = clock();
    //_RBF.Greedy_algorithm(_wall_nodes, 1e-7);   //和分组RBF相同的误差设置
    //auto t_global_rbf_end = clock();

    //auto t_global_df_start = clock();
    //_RBF.calculate_every_node_deform(_node_coords);
    //auto t_global_df_end = clock();

    //    std::cout << "_RBF coeff calculating time: " << fixed << setprecision(2) << double(t_global_rbf_end - t_global_rbf_start) << " ms" << std::endl;
    //std::cout << "_RBF df calculating time: " << fixed << setprecision(2) << double(t_global_df_end - t_global_df_start) << " ms" << std::endl;

    //auto t_global_deform_start = clock();
    //_calculate_deformed_coordinates(_node_coords);
    //auto t_global_deform_end = clock();

    //auto t_global_end = clock();

    //std::cout << "_RBF deforming nodes time: " << fixed << setprecision(2) << double(t_global_deform_end - t_global_deform_start) << " ms" << std::endl;
    //std::cout << "_RBF total time: " << fixed << setprecision(2) << double(t_global_end - t_global_start) << " ms" << std::endl;
    ////***************全局RBF算法结束**********************
     
     
     
   
    
    
    //测试class multi_partition
    std::vector<idx_t> v_temp{9,2,2};
    std::vector<double> steps_tol{ 1e-14,1e-14,1e-7};
    multi_partition mp(v_temp);
    //mp.set_m_nd2wall_lvl_mp(_nd2wall_lvl_mp);
    clock_t blk_divide_start, blk_divide_end;
    clock_t blk_rbf_start, blk_rbf_end;

    blk_divide_start = clock();
    mp.divide_wall(_every_boundary[0].bound_elements, _wall_nodes);
    blk_divide_end = clock();

    blk_rbf_start = clock();
    mp.multi_partition_rbf_algorithm(steps_tol, _wall_nodes, _node_coords);
    blk_rbf_end = clock();

    std::cout << "divide wall use time: " << fixed << setprecision(2) << double(blk_divide_end - blk_divide_start) << std::endl;
    std::cout << "block rbf coeff calculating and deforming mesh use time: " << fixed << setprecision(2) << double(blk_rbf_end - blk_rbf_start) << std::endl;
    
    ////测试壁面数据
    //_wall_nodes.resize(0);
    //_wall_nodes = Set_wall_nodes(_every_boundary, _node_coords);
    //
    //_write_wall(_wall_nodes);
    //
    //
    writefile(writefile_name);

   


    return 0;
}