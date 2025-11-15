#pragma once

#include <vector>
#include <unordered_map>
#include <Eigen/Core>
#include "distance.hpp"
#include "geometry.hpp"

struct State
{
    //读取文件相关包括： 网格维度、单元数、节点数
    int Dimension = 0;   //网格维度
    int Nelements = 0;   //单元数
    int NPoints = 0;    //节点数
    int wall_id = 0;    //物面边界在边界数组的索引

    //网格数据相关包括： 所有节点坐标、所有单元信息、所有边界数组、物面节点数组
    std::vector<Node> node_coords;  //保存所有节点的坐标的数组，在node_coords的索引就对应着具体编号
    std::vector<element> every_elements;     //保存所有单元的数组，每个元素保存有其节点编号
    std::vector<boundary> every_boundary;  //保存所有边界单元的数组，每个元素是一个边界
    std::vector<Node> wall_nodes;        //物面节点数组

    //算法参数
    double R = 30.0;      //支撑距离 phi = phi(||r-ri||/R)
    double R_squared = R * R;  // R的平方，用于优化distance在基函数的sqrt计算
    double invR = 1 / R;    //invR，用于优化除法效率
    double D = 0.0;       //限制距离 psi = 1-r/D, D是最大的位移变形量的五倍
    double alpha = 10.0;  //DRRBF参数, 用于将静止的物面网格排出变形区域
    double beta = 1.0;   //DRRBF参数
    double invBeta = 1 / beta;  // 用于优化除法效率

    //查询结构
    GridBTree<int, Point<double>> wall_tree;   //根据物面节点构造的二叉树
    std::unordered_map<int, int> wall_id2node; //物面节点id和node的查找集
};
