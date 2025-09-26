#include "rbf.hpp"
#include "State.hpp"
#include "deformation.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <unordered_set>
using namespace std;

//====================== RBFInterpolator 类实现======================

// 构造函数，初始化A、b、coeff的大小为3
RBFInterpolator::RBFInterpolator()
{
    A.resize(3);     // 初始化A为3个元素
    b.resize(3);     // 初始化b为3个元素
    coeff.resize(3); // 初始化coeff为3个元素
}

// 计算物面变形与误差的函数
void RBFInterpolator::calculate_df(const vector<Eigen::VectorXd> &coeff, const vector<Node> &suppoints,
                                   vector<Node> &P_df, vector<double> &interp_tol,
                                   int &max_tol_id, const vector<Node> &wall_nodes, const State &S)
// 计算物面变形与误差的函数
// 传入参数：插值系数 coeff
//          支撑点数组 suppoints
//          物面节点数组 wall_nodes
//          最大插值误差的索引 max_tol_id
//          变形量数组 P_df
//          插值误差数组 interp_tol
//          参数结构体 S
// 该函数会引发传入参数P_df和interp_tol的改变!!!
{
    Eigen::Vector3d df_ij{0, 0, 0};
    double tol_ij = 0.0;
    double max_tol = 0.0;
    max_tol_id = -1;
    for (int i = 0; i < wall_nodes.size(); ++i)
    {
        for (int j = 0; j < suppoints.size(); ++j)
        {
            for (int k = 0; k < coeff.size(); ++k) // coeff是个size为3的vector
            {
                df_ij[k] += coeff[k][j] * rbf_func_Wendland(_distance(wall_nodes[i].point, suppoints[j].point), S.R); // 第i个空间节点在第j个支撑点上第k个空间维度上的变形值
            }
        }
        P_df[i].df = df_ij;                         // 计算得到的变形量
        tol_ij = (wall_nodes[i].df - df_ij).norm(); // 计算误差，df是目标变形量
        interp_tol[i] = tol_ij;
        if (tol_ij > max_tol)
        {
            max_tol = tol_ij;
            max_tol_id = i;
        }
        df_ij.setZero();
    }
}

// 使用贪心算法选择支撑点并且计算，该算法适用于rbf系统还没有设置external_suppoints的情况
void RBFInterpolator::Greedy_algorithm(double tol, const State &S) // 使用贪心算法选择支撑点并且计算，该算法适用于rbf系统已设置external_suppoints的情况
                                                                   // 传入参数：误差 tol
                                                                   //          参数结构体 S
{
    // std::cout << "=====================" << endl;
    // std::cout << "Starting Greedy algorithm ... " << std::endl;

    interp_tol.resize(external_suppoints.size());
    P_df = external_suppoints;
    int n = 0; // 支撑点数，也是插值系统的维度
    double eta_ij = 0;
    int max_tol_id = 0;
    double max_tol_i = 0;
    for (int i = 0; i < external_suppoints.size(); ++i)
    {
        if (i == 0)
        {
            suppoints.push_back(external_suppoints[i]); // 第一步任选一个点即可

            for (int k = 0; k < b.size(); ++k)
            {
                A[k].resize(1, 1);
                b[k].resize(1);
                eta_ij = _distance(suppoints[i].point, suppoints[i].point);
                A[k](i, i) = rbf_func_Wendland(eta_ij, S.R); // 对角线元素
                b[k](i) = suppoints[i].df[k];
                coeff[k] = A[k].llt().solve(b[k]); // Cholesky分解插值系数
            }

            calculate_df(coeff, suppoints, P_df, interp_tol, max_tol_id, external_suppoints, S);
            // 计算变形量，插值误差以及最大误差所对应的max_tol_id，这个对应的是在wall_nodes里的索引
            // std::cout << "Step: " << i + 1 << std::endl;
            // std::cout << "max tol: " << std::fixed << std::setprecision(13) << interp_tol[max_tol_id] << std::endl;
        }
        else
        {
            suppoints.push_back(external_suppoints[max_tol_id]);
            n = suppoints.size();
            for (int k = 0; k < b.size(); ++k)
            {
                A[k].conservativeResize(n, n);
                b[k].conservativeResize(n);
                for (int j = 0; j < n - 1; ++j)
                {
                    eta_ij = _distance(suppoints[n - 1].point, suppoints[j].point);
                    A[k](n - 1, j) = rbf_func_Wendland(eta_ij, S.R); // 插入矩阵新的一行
                    A[k](j, n - 1) = rbf_func_Wendland(eta_ij, S.R); // 插入矩阵新的一列，矩阵是对称正定的

                    eta_ij = _distance(suppoints[n - 1].point, suppoints[n - 1].point);
                    A[k](n - 1, n - 1) = rbf_func_Wendland(eta_ij, S.R); // 对角线元素
                    b[k](n - 1) = suppoints[n - 1].df[k];                // b向量最后一个元素
                    coeff[k] = A[k].llt().solve(b[k]);                   // Cholesky分解插值系数
                }
            }
            calculate_df(coeff, suppoints, P_df, interp_tol, max_tol_id, external_suppoints, S);
            // 计算变形量，插值误差以及最大误差所对应的max_tol_i，这个对应的是在wall_nodes里的索引
            if (i % 99 == 0) // 打印信息
            {
                // 计算变形量，插值误差以及最大误差所对应的max_tol_id，这个对应的是在wall_nodes里的索引
                // std::cout << "Step: " << i + 1 << std::endl;
                // std::cout << "max tol: " << std::fixed << std::setprecision(13) << interp_tol[max_tol_id] << std::endl;
            }
            if (interp_tol[max_tol_id] < tol)
            {
                // std::cout << "Step: " << i + 1 << std::endl;
                // std::cout << "max tol: " << std::fixed << std::setprecision(13) << interp_tol[max_tol_id] << std::endl;
                // std::cout << "The setting tolerance has been reached" << std::endl;
                break;
            }
            if (i == external_suppoints.size() - 1)
            {
                // std::cout << "Step: " << i + 1 << std::endl;
                // std::cout << "max tol: " << std::fixed << std::setprecision(13) << interp_tol[max_tol_id] << std::endl;
                std::cout << "every boundary nodes have been chosen (this will happen especially when the nodes are selected in boundary nodes shared by blocks)" << std::endl;
            }
        }
    }
}

// 使用贪心算法选择支撑点并且计算，该算法适用于rbf系统已设置external_suppoints的情况
void RBFInterpolator::Greedy_algorithm(const std::vector<Node> &wall_nodes, double tol, const State &S)
// 传入参数：物面节点 _wall_nodes
//          误差 tol
{
    // std::cout << "=====================" << endl;
    // std::cout << "Starting Greedy algorithm ... " << std::endl;

    interp_tol.resize(wall_nodes.size());
    P_df = wall_nodes;
    int n = 0; // 支撑点数，也是插值系统的维度
    double eta_ij = 0;
    int max_tol_id = 0;
    double max_tol_i = 0;
    for (int i = 0; i < wall_nodes.size(); ++i)
    {
        if (i == 0)
        {
            suppoints.push_back(wall_nodes[i]); // 第一步任选一个点即可

            for (int k = 0; k < b.size(); ++k)
            {
                A[k].resize(1, 1);
                b[k].resize(1);
                eta_ij = _distance(suppoints[i].point, suppoints[i].point);
                A[k](i, i) = rbf_func_Wendland(eta_ij, S.R); // 对角线元素
                b[k](i) = suppoints[i].df[k];
                coeff[k] = A[k].llt().solve(b[k]); // Cholesky分解插值系数
            }

            calculate_df(coeff, suppoints, P_df, interp_tol, max_tol_id, wall_nodes, S);
            // 计算变形量，插值误差以及最大误差所对应的max_tol_id，这个对应的是在wall_nodes里的索引
            // std::cout << "Step: " << i + 1 << std::endl;
            // std::cout << "max tol: " << std::fixed << std::setprecision(13) << interp_tol[max_tol_id] << std::endl;
        }
        else
        {
            suppoints.push_back(wall_nodes[max_tol_id]);
            n = suppoints.size();
            for (int k = 0; k < b.size(); ++k)
            {
                A[k].conservativeResize(n, n);
                b[k].conservativeResize(n);
                for (int j = 0; j < n - 1; ++j)
                {
                    eta_ij = _distance(suppoints[n - 1].point, suppoints[j].point);
                    A[k](n - 1, j) = rbf_func_Wendland(eta_ij, S.R); // 插入矩阵新的一行
                    A[k](j, n - 1) = rbf_func_Wendland(eta_ij, S.R); // 插入矩阵新的一列，矩阵是对称正定的

                    eta_ij = _distance(suppoints[n - 1].point, suppoints[n - 1].point);
                    A[k](n - 1, n - 1) = rbf_func_Wendland(eta_ij, S.R); // 对角线元素
                    b[k](n - 1) = suppoints[n - 1].df[k];                // b向量最后一个元素
                    coeff[k] = A[k].llt().solve(b[k]);                   // Cholesky分解插值系数
                }
            }
            calculate_df(coeff, suppoints, P_df, interp_tol, max_tol_id, wall_nodes, S);
            // 计算变形量，插值误差以及最大误差所对应的max_tol_i，这个对应的是在wall_nodes里的索引
            if (i % 99 == 0) // 打印信息
            {
                // 计算变形量，插值误差以及最大误差所对应的max_tol_id，这个对应的是在wall_nodes里的索引
                // std::cout << "Step: " << i + 1 << std::endl;
                // std::cout << "max tol: " << std::fixed << std::setprecision(13) << interp_tol[max_tol_id] << std::endl;
            }
            if (interp_tol[max_tol_id] < tol)
            {
                // std::cout << "Step: " << i + 1 << std::endl;
                // std::cout << "max tol: " << std::fixed << std::setprecision(13) << interp_tol[max_tol_id] << std::endl;
                // std::cout << "The setting tolerance has been reached" << std::endl;
                break;
            }
            if (i == wall_nodes.size() - 1)
            {
                // std::cout << "Step: " << i + 1 << std::endl;
                // std::cout << "max tol: " << std::fixed << std::setprecision(13) << interp_tol[max_tol_id] << std::endl;
                // std::cout << "every boundary nodes have been chosen (this will happen especially when the nodes are selected in boundary nodes shared by blocks)" << std::endl;
            }
        }
    }
}

// 贪心算法误差输出至文件查看
void RBFInterpolator::write_greedy_tol(const std::string &filename)
{
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }

    // 设置固定小数格式，保留12位小数
    outfile << std::fixed << std::setprecision(12);

    // Tecplot 标准头部
    outfile << "TITLE = \"" << filename << "\"\n";
    outfile << "VARIABLES = \"step\" \"tol\"\n";

    // 输出误差曲线数据
    for (int i = 0; i < max_tol_step.size(); ++i)
    {
        outfile << i + 1 << " " << max_tol_step[i] << "\n";
    }

    outfile.close();
}

//====================== DeformCalculator ======================
void DeformCalculator::calculate_deform(std::vector<Node> &every_nodes, const State &S) const
// 计算所有节点的变形（纯RBF）
// 传入参数：所有节点 every_nodes
//          参数结构体 S
//          S包括数据: R
{
    for (auto &node : every_nodes)
    {
        calculate_deform_RBF(node, S);
    }
}

void DeformCalculator::calculate_deform_RBF(Node &inode, const State &S) const
// 计算单个节点的变形（纯RBF）
// 传入参数：单个节点 inode
//          参数结构体 S
//          S包括数据: R
{
    Eigen::Vector3d df_ij{0, 0, 0};
    for (int j = 0; j < rbf.suppoints.size(); ++j)
    {
        for (int k = 0; k < rbf.coeff.size(); ++k) // coeff是个size为3的vector
        {
            df_ij[k] += rbf.coeff[k][j] * rbf_func_Wendland(_distance(inode.point, rbf.suppoints[j].point), S.R); // 第i个空间节点在第j个支撑点上第k个空间维度上的变形值
        }
    }
    inode.df = df_ij; // 计算得到的变形量
}

void DeformCalculator::calculate_deform_DRRBF(Node &inode, double d_r2omega1, double d_r2omega2, const State &S) const
// 计算单个节点的变形（RRBF）
// 传入参数：单个节点 inode
//          待插值点到动边界的距离 d_r2omega1
//          待插值点到静止边界的距离 d_r2omega2
//          参数结构体 S
//          S包括数据: alpha, beta, D, R
{
    Eigen::Vector3d df_ij{0, 0, 0};
    for (int j = 0; j < rbf.suppoints.size(); ++j)
    {
        for (int k = 0; k < rbf.coeff.size(); ++k) // coeff是个size为3的vector
        {
            if (d_r2omega1 > S.D)
                continue; // 超过限制半径就直接为0
            df_ij[k] += rbf.coeff[k][j] * psi(d_r2omega1, d_r2omega2, S) * rbf_func_Wendland(_distance(inode.point, rbf.suppoints[j].point), S.R);
        }
    }
    inode.df = df_ij; // 计算得到的变形量
}

void set_block_rbf(std::vector<DeformCalculator> &block_rbf, const std::vector<mesh_block> &blocks)
// =====================================================
// 生成每个网格块的 RBF 插值系统候选支撑点集 external_suppoints
// - block_rbf[i] 对应 blocks[i]
// - external_suppoints 包括：
//   * 该块的所有内部点
//   * 在 block_D 覆盖范围内的其它块内部点和边界点
// =====================================================
{
    // 临时变量
    double d_temp = 0.0;
    double d_r2omgea_temp = 0.0;
    int key_temp = 0;

    // 遍历每个网格块
    for (int i = 0; i < blocks.size(); ++i)
    {
        const auto &blk = blocks[i];                           // 当前块
        RBFInterpolator &rbf = block_rbf[i].get_rbf_mutable(); // 当前块的 RBF 插值系统

        std::unordered_set<int> nodes_id_set; // 用于存储支撑点的唯一 ID
        d_temp = blk.block_D;                 // 当前块的支撑距离

        // 遍历所有分区
        for (int j = 0; j < blocks.size(); ++j)
        {
            for (const auto &nd : blocks[j].internal_nodes) // 遍历分区j的所有内部节点
            {
                if (j == i)
                {
                    // 自己分区的内部节点，直接加入
                    rbf.external_suppoints.push_back(nd);
                    nodes_id_set.insert(nd.id);
                }
                else
                {
                    // 其他分区的内部节点
                    if (nodes_id_set.find(nd.id) == nodes_id_set.end()) // 如果该节点ID不在集合中
                    {
                        blk.block_tree.search(d_r2omgea_temp, key_temp, nd.point); // 查询该节点到当前块的最小距离
                        if (d_r2omgea_temp <= d_temp)                              // 如果距离小于等于支撑距离
                        {
                            Node nd_copy = nd;                         // 复制节点
                            nd_copy.df.setZero();                      // 将变形量设为零(其他分区为静止边界)
                            rbf.external_suppoints.push_back(nd_copy); // 加入支撑点集
                            nodes_id_set.insert(nd.id);                // 将节点ID加入集合
                        }
                    }
                }
            }

            for (const auto &nd : blocks[j].boundary_nodes) // 遍历分区j的所有边界节点
            {
                if (nodes_id_set.find(nd.id) == nodes_id_set.end()) // 如果该节点ID不在集合中
                {
                    blk.block_tree.search(d_r2omgea_temp, key_temp, nd.point); // 查询该节点到当前块的最小距离
                    if (d_r2omgea_temp <= d_temp)                              // 如果距离小于等于支撑距离
                    {
                        Node nd_copy = nd; // 复制节点
                        nd_copy.df.setZero(); // 将变形量设为零
                        rbf.external_suppoints.push_back(nd_copy); // 加入支撑点集
                        nodes_id_set.insert(nd.id);                // 将节点ID加入集合
                    }
                }
            }
        }
    }
}