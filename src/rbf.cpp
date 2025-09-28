#include "rbf.hpp"
#include "block.hpp"
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

// ------------------------------------------------------------
// 重载 1：使用 unique_bndry 作为每块的“内部核心候选”
// external_suppoints = 本块 unique_bndry
//                    + 以 block_D 覆盖到的其他块 internal_nodes / boundary_nodes（df 置 0）
// ------------------------------------------------------------
void set_block_rbf(std::vector<DeformCalculator>& block_rbf,
                   const std::vector<std::vector<Node>>& unique_bndry,
                   const std::vector<mesh_block>& blocks)
{
    // 基本一致性检查
    if (block_rbf.size() != blocks.size() || unique_bndry.size() != blocks.size()) {
        std::cerr << "[set_block_rbf(unique_bndry)] size mismatch: rbfs="
                  << block_rbf.size() << ", unique_bndry=" << unique_bndry.size()
                  << ", blocks=" << blocks.size() << std::endl;
        std::abort();
    }

    double d_temp = 0.0;        // 当前块的 block_D
    double d_near = 0.0;        // 最近距离
    int    key_tmp = 0;         // 树查询返回的 id 占位

    for (int i = 0; i < static_cast<int>(blocks.size()); ++i)
    {
        const auto& blk_i = blocks[i];              //获取当前的网格块
        auto& rbf = block_rbf[i].get_rbf_mutable();  //获取当前网格块对应的RBF对象

        // 清空旧的候选集，避免多次调用累积
        rbf.external_suppoints.clear();

        std::unordered_set<int> used_ids;   // 加入插值系统的待支撑点(external points)的查找集，用于确保点是唯一的（去重）
        d_temp = blk_i.block_D;             // 设置当前block[i]的覆盖范围

        // 1) 本块的 unique_bndry 直接加入（保留其 df）
        for (const auto& nd : unique_bndry[i]) 
        {
            if (used_ids.insert(nd.id).second)    //将已插入的点放入查找集里，若点已存在则该条件为否
            {
                rbf.external_suppoints.push_back(nd);  //向RBF系统插入该网格块的unique_bndry点
            }
        }

        // 2) 其它块被 block_D 覆盖到的 internal / boundary
        for (int j = 0; j < static_cast<int>(blocks.size()); ++j) 
        {
            if (j == i) continue; // 自己的 unique 已经在上一步加入

            // 2.1 internal
            for (const auto& nd : blocks[j].internal_nodes) //检查内部点
            {
                blk_i.block_tree.search(d_near, key_tmp, nd.point);  //查找当前点到blocks[i]的最小距离
                if (d_near <= d_temp && used_ids.insert(nd.id).second)        //这个点被该网格块的block_D形成的圆覆盖，且不在查找集内则插入
                {
                    Node nd_copy = nd;
                    nd_copy.df.setZero();                   // 他块视作静止
                    rbf.external_suppoints.push_back(nd_copy);
                }
            }

            // 2.2 boundary
            for (const auto& nd : blocks[j].boundary_nodes) //检查边界点
            {
                blk_i.block_tree.search(d_near, key_tmp, nd.point);  //查找当前点到blocks[i]的最小距离
                if (d_near <= d_temp && used_ids.insert(nd.id).second)     //这个点被该网格块的block_D形成的圆覆盖，且不在查找集内则插入
                {
                    Node nd_copy = nd;
                    nd_copy.df.setZero();    // 他块视作静止
                    rbf.external_suppoints.push_back(nd_copy);
                }
            }
        }
    }
}

// ------------------------------------------------------------
// 重载 2：仅依据各块 own internal + block_D 扫入其它块的 internal/boundary
// external_suppoints = 本块 internal（保留 df）
//                    + 以 block_D 覆盖到的所有块 boundary（含本块）df=0
//                    + 以 block_D 覆盖到的其它块 internal df=0
// ------------------------------------------------------------
void set_block_rbf(std::vector<DeformCalculator>& block_rbf,
                   const std::vector<mesh_block>& blocks)
{
    if (block_rbf.size() != blocks.size()) {
        std::cerr << "[set_block_rbf(blocks)] size mismatch: rbfs="
                  << block_rbf.size() << ", blocks=" << blocks.size() << std::endl;
        std::abort();
    }

    double d_temp = 0.0;  // block_D
    double d_near = 0.0;  // 最近距离
    int    key_tmp = 0;   // 树查询返回的 id 占位

    for (int i = 0; i < static_cast<int>(blocks.size()); ++i)
    {
        const auto& blk_i = blocks[i];            //获取当前的网格块
        auto& rbf = block_rbf[i].get_rbf_mutable(); //获取当前网格块对应的RBF对象

        // 清空旧的候选集，避免多次调用累积
        rbf.external_suppoints.clear();

        std::unordered_set<int> used_ids;
        d_temp = blk_i.block_D;

        // 1) 自己块的 internal 直接加入（保留 df）
        for (const auto& nd : blocks[i].internal_nodes) 
        {
            if (used_ids.insert(nd.id).second)    //将已插入的点放入查找集里，若点已存在则该条件为否
            {
                rbf.external_suppoints.push_back(nd);
            }
        }

        // 2) 所有块的 boundary，如果被覆盖则加入（包含本块 boundary）
        for (int j = 0; j < static_cast<int>(blocks.size()); ++j) 
        {
            for (const auto& nd : blocks[j].boundary_nodes)  //检查边界点
            {
                blk_i.block_tree.search(d_near, key_tmp, nd.point);   //查找当前点到blocks[i]的最小距离
                if (d_near <= d_temp && used_ids.insert(nd.id).second)     //这个点被该网格块的block_D形成的圆覆盖，且不在查找集内则插入
                {
                    Node nd_copy = nd;
                    nd_copy.df.setZero();                 // 边界默认静止
                    rbf.external_suppoints.push_back(nd_copy);
                }
            }
        }

        // 3) 其它块的 internal，如果被覆盖则加入（df=0）
        for (int j = 0; j < static_cast<int>(blocks.size()); ++j) 
        {
            if (j == i) continue;
            for (const auto& nd : blocks[j].internal_nodes)  //检查内部点
            {
                blk_i.block_tree.search(d_near, key_tmp, nd.point);   //查找当前点到blocks[i]的最小距离
                if (d_near <= d_temp && used_ids.insert(nd.id).second)   //这个点被该网格块的block_D形成的圆覆盖，且不在查找集内则插入
                {
                    Node nd_copy = nd;
                    nd_copy.df.setZero();
                    rbf.external_suppoints.push_back(nd_copy);
                }
            }
        }
    }
}
