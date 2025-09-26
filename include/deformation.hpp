#pragma once

#include <vector>
#include <Eigen/Core>
#include "geometry.hpp"

// 两点距离
double _distance(const Point<double>& P1, const Point<double>& P2);

// 物面节点提取函数
std::vector<Node> Set_wall_nodes(const std::vector<boundary>& every_boundary, 
    const std::vector<Node>& node_coords);

// 计算支撑距离R，写入State
void select_R(const std::vector<Node>& node_coords, double& R);

// 物面变形函数，计算的是变形量（可替换）
Eigen::Vector3d deforming_fuc(const Point<double>& wall_nodes);

// 计算物面变形并更新 State.D
void calculat_wall_deformation(std::vector<Node>& wall_nodes, double& D);

// 计算节点变形后的坐标
void calculate_deformed_coordinates(std::vector<Node>& every_nodes);