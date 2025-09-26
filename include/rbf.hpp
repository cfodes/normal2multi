#pragma once

#include "State.hpp"
#include "block.hpp"
#include <vector>
#include <Eigen/Dense>



// rbf基函数定义：Wendland's C2 function
inline double rbf_func_Wendland(double eta, double R)
{
    eta = eta / R;    //传入的eta = ri-rj, 为了更好的插值精度计算：eta = (ri-rj)/R
    if (eta >= 0 && eta <= 1)
    {
        return pow((1 - eta), 4) * (4 * eta + 1);
    }
    else
    {
        return 0;
    }
}

// DDRBF的psi函数定义
inline double psi(double d_r2omega1,double d_r2omega2, const State& S)  //计算RRBF：I(r)=~phi*(omgea_i*phi_i)  哑标代表求和
//传入参数：d_r2omega1   点到运动物面点的最小距离
//         d_r2omega2   点到静止物面点的最小距离
//         参数结构体 S
//         S包括数据: alpha, beta, D
{
    double psi1 = 0; double psi2 = 0;
    double psi2_temp = 0;
    if (d_r2omega1 >= 0 && d_r2omega1 <= 1)
    {
        psi1 = 1 - d_r2omega1 / S.D;
    }
    else if (d_r2omega1 > 0)   //实际上为了减少运算，这段可以省略，在待插值点循环的时候直接判断即可
    {
        psi1 = 0;
    }
    psi2_temp = (S.alpha * d_r2omega2 - d_r2omega1) / S.beta;
    psi2 = std::min(1.0, std::max(0.0, psi2_temp));
    return std::min(psi1, psi2);
}

//RBF插值系统类，该类完成插值系数计算
class RBFInterpolator 
{
public:
    std::vector<Node> suppoints;   //选择的支撑点
    std::vector<Node> external_suppoints;   //外部传入的待选取的支撑点
    std::vector<double> interp_tol;   //插值误差数组，索引和wall_nodes对应
    std::vector<Node> P_df;           //变形量矩阵，需要和物面节点对应，该数组变形量会变化
    std::vector<double> max_tol_step; //每一步的最大误差
    std::vector<Eigen::MatrixXd> A;      //系数矩阵，xyz三个方向有三个
    std::vector<Eigen::VectorXd> b;      //已知位移向量，xyz三个方向有三个
    std::vector<Eigen::VectorXd> coeff;   //插值系数，xyz三个方向有三个

    // 构造函数，确保A, b, coeff的size为3
    RBFInterpolator(); 

    //计算物面变形与误差的函数
    void calculate_df(const std::vector<Eigen::VectorXd>& coeff, const std::vector<Node>& suppoints,
        std::vector<Node>& P_df, std::vector<double>& interp_tol,
         int& max_tol_id, const std::vector<Node>& wall_nodes, const State& S);
    
    //使用贪心算法选择支撑点并且计算，该算法适用于rbf系统已设置external_suppoints的情况
    void Greedy_algorithm(double tol, const State& S);

    //使用贪心算法选择支撑点并且计算，该算法适用于rbf系统还没有设置external_suppoints的情况
    void Greedy_algorithm(const std::vector<Node>& wall_nodes, double tol, const State& S);

    //贪心算法误差输出至文件查看
    void write_greedy_tol(const std::string& filename);
};

// 利用插值系数计算变形的工具类
class DeformCalculator
{
private:
    RBFInterpolator& rbf;

public:
    explicit DeformCalculator(RBFInterpolator& r): rbf(r) {}

    // 获取rbf对象的引用
    RBFInterpolator& get_rbf_mutable() { return rbf; }

    // 只读的rbf对象
    const RBFInterpolator& get_rbf() const { return rbf; }

    // 计算所有节点的变形（纯RBF）
    void calculate_deform(std::vector<Node>& every_nodes, const State& S) const;

    // 计算单个节点的变形（纯RBF）
    void calculate_deform_RBF(Node& inode, const State& S) const;  

    // 计算单个节点的变形（RRBF）
    void calculate_deform_DRRBF(Node& inode, 
                                double d_r2omega1, double d_r2omega2, 
                                const State& S) const;
    
};

//根据网格块组生成每个网格块的rbf插值系统的支撑点集
void set_block_rbf(std::vector<DeformCalculator>& block_rbf, const std::vector<mesh_block>& blocks);