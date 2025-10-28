#pragma once

#include <vector>
#include <Eigen/Core>
#include "geometry.hpp"

// 两点距离(传点类版本)

inline double _distance(const Point<double>& P1, const Point<double>& P2)
{
    return std::sqrt((P1.x - P2.x) * (P1.x - P2.x)
                     + (P1.y - P2.y) * (P1.y - P2.y)
                     + (P1.z - P2.z) * (P1.z - P2.z));
}

// 两点距离(传点坐标版本，优化效率用)
inline double _distance(double xi, double yi, double zi,
    double xj, double yj, double zj) noexcept
{
    const double dx = xi - xj;
    const double dy = yi - yj;
    const double dz = zi - zj;

    // 返回两点欧式距离
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// 物面节点提取函数
std::vector<Node> Set_wall_nodes(const std::vector<boundary>& every_boundary, 
    const std::vector<Node>& node_coords, int& wall_id);

// 计算支撑距离R，写入State
void select_R(const std::vector<Node>& node_coords, double& R);

// 物面变形函数，计算的是变形量（可替换）
Eigen::Vector3d deforming_fuc(const Point<double>& wall_nodes);

// 计算物面变形并更新 State.D
void calculat_wall_deformation(std::vector<Node>& wall_nodes, double& D);

// 计算节点变形后的坐标
void calculate_deformed_coordinates(std::vector<Node>& every_nodes);