#include "deformation.hpp"
#include <iostream>
#include <unordered_set>
#include <cmath>
#define PI acos(-1)
using namespace std;

// 两点距离
double _distance(const Point<double>& P1, const Point<double>& P2) {
    return std::sqrt((P1.x - P2.x) * (P1.x - P2.x)
                   + (P1.y - P2.y) * (P1.y - P2.y)
                   + (P1.z - P2.z) * (P1.z - P2.z));
}

// 物面节点提取函数
//确定物面节点的数组，该函数确保节点不重复
//!!!!该函数目前只支持一个wall的情况
vector<Node> Set_wall_nodes(const vector<boundary>& every_boundary, const vector<Node>& node_coords)
{
    int inode_id = 0;
    unordered_set<int> wall_nodes_set;
    vector<Node> local_wall_nodes;
    
    for (int i = 0; i < every_boundary.size(); ++i)
    {
        if (every_boundary[i].boundary_tag == "wall")  //查找边界数组里壁面的数组
        {
            std::cout << "=====================" << endl;
            std::cout << "Setting wall's nodes array ..." << std::endl;
            for (int j = 0; j < every_boundary[i].bound_elements.size(); ++j)  //遍历壁面数组的每个单元
            {
                for (int k = 0; k < every_boundary[i].bound_elements[j].node_id.size(); ++k)  //遍历壁面数组一个单元的所有节点id
                {
                    inode_id = every_boundary[i].bound_elements[j].node_id[k];
                    if (wall_nodes_set.find(inode_id) == wall_nodes_set.end())  //如果set里没有该节点编号，则会返回set.end()，所以该条件判断set里有没有该节点编号
                    {
                        wall_nodes_set.insert(inode_id);    //set插入新的id
                        local_wall_nodes.push_back(node_coords[inode_id]);   //wall_nodes插入对应节点
                    }
                }
                
            }
        }
    }
    std::cout << "Setting wall's nodes array successfully" << std::endl;
    return local_wall_nodes;
}


// 计算支撑距离R，写入State
// 计算 R（包围盒尺度）
void select_R(const std::vector<Node>& nodes, double& R)
 {
    // 计算包围盒对角线长度作为 R
    // Point<double> max_pt(-std::numeric_limits<double>::infinity(),
    //                      -std::numeric_limits<double>::infinity(),
    //                      -std::numeric_limits<double>::infinity());
    // Point<double> min_pt( std::numeric_limits<double>::infinity(),
    //                       std::numeric_limits<double>::infinity(),
    //                       std::numeric_limits<double>::infinity());
    // for (const auto& node : nodes) {
    //     max_pt.x = std::max(max_pt.x, node.point.x);
    //     max_pt.y = std::max(max_pt.y, node.point.y);
    //     max_pt.z = std::max(max_pt.z, node.point.z);
    //     min_pt.x = std::min(min_pt.x, node.point.x);
    //     min_pt.y = std::min(min_pt.y, node.point.y);
    //     min_pt.z = std::min(min_pt.z, node.point.z);
    // }
    // double dx = max_pt.x - min_pt.x;
    // double dy = max_pt.y - min_pt.y;
    // double dz = max_pt.z - min_pt.z;
    // R = std::sqrt(dx * dx + dy * dy + dz * dz);

    // 用户设置（网格较大的时候，根据变形范围设置）
    R = 10;

    std::cout << "R is determined as " << R << std::endl;
}

//物面变形函数 计算的是变形量
Eigen::Vector3d deforming_fuc(const Point<double>& wall_nodes)
{
    Eigen::Vector3d v1{ 0,0,0 };
    v1[0] = 0;   //x方向的变形量
    v1[1] = 0.1 * sin(2 * PI * wall_nodes.x);//y方向的变形量
    v1[2] = 0;   //z方向的变形量
    return v1;
}

//物面变形量计算，将D写入State.D
void calculat_wall_deformation(vector<Node>& wall_nodes, double& D)   //该函数设置物面变形，同时计算物面最大变形量，_D是五倍最大的位移变形量
{
    double D_max = 0;
    double D_temp = 0;
    for (int i = 0; i < wall_nodes.size(); ++i)
    {
        wall_nodes[i].df = deforming_fuc(wall_nodes[i].point);
        D_temp = wall_nodes[i].df.norm();
        if (D_temp > D_max)
        {
            D_max = D_temp;
        }
    }
    D = 5 * D_max;
    std::cout << "=====================" << endl;
    std::cout << "points' deformations on wall have benn set" << endl;
    std::cout << "D is determined as " << D << endl;
}

// 计算节点变形后的坐标
void calculate_deformed_coordinates(vector<Node>& every_node)   //计算节点变形后的坐标
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
    std::cout << "Calculating deformed coordinates successfully" << std::endl;
}   