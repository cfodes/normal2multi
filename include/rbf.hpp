#pragma once

#include "State.hpp"
#include "block.hpp"
#include <vector>
#include <Eigen/Dense>



// rbf基函数定义：Wendland's C2 function
inline double rbf_func_Wendland(double dist, double R, double invR)
//传入的dist为||ri-rj||
// R是支撑半径
//invR为 1/R，由于R是支撑半径，提前已知，该值用于优化效率
//Wendland's C2 function的定义为： phi(eta)=(1-eta)^4(4eta+1) , 0<=eta<1 
//                                          0 , eta>=1
//令t=1-eta, phi = t^4(5-4t)
{
    if (dist >= R) return 0.0;
    const double eta = dist * invR;
    const double t = 1.0 - eta;
    const double t2 = t * t;
    return(5.0 - 4.0 * t) * (t2 * t2);
}





//// DDRBF的psi函数定义
//// 效率优化版本
//inline double psi(double d_r2omega1, double d_r2omega2, double D, const State& S)
//{
//    double psi1 = 0.0, psi2 = 0.0;
//
//    // 用 D 作为半径阈值
//    if (d_r2omega1 >= 0 && d_r2omega1 <= D)
//        psi1 = 1.0 - d_r2omega1 / D;
//    else if (d_r2omega1 > D)  // 超过限制半径就直接为0
//        psi1 = 0.0;
//
//    const double psi2_temp = (S.alpha * d_r2omega2 - d_r2omega1) / S.beta;
//    psi2 = std::min(1.0, std::max(0.0, psi2_temp));
//
//    return std::min(psi1, psi2);
//}

// DDRBF的psi函数定义
// 效率优化版本，尽量避免除法和多余乘法的运算
inline double psi(double d1, double d2, double D, const State& S)
{
    if (d1 >= D) return 0.0;   //超过限制半径直接为0
    const double invD = 1.0 / D;
    
    // ~psi =1-d1/D
    const double psi1 = 1 - d1 * invD;

    // ~psi2 = min(1,max(0,(alpha*d2-d1)/beta))
    const double tmp = (S.alpha * d2 - d1) * S.invBeta;
    const double psi2 = std::min(1.0, std::max(0.0, tmp));

    // 返回较小者
    return std::min(psi1, psi2);
}


//RBF插值系统类，该类完成插值系数计算
class RBFInterpolator 
{
public:


    // public 增加
    void set_debug_trace(bool on, const std::string& tag = {}) {
        trace_ = on;
        trace_tag_ = tag;
    }


    std::vector<Node> suppoints;   //选择的支撑点
    std::vector<Node> external_suppoints;   //外部传入的待选取的支撑点
    std::vector<double> interp_tol;   //插值误差数组，索引和wall_nodes对应
    std::vector<Node> P_df;           //变形量矩阵，需要和物面节点对应，该数组变形量会变化
    std::vector<double> max_tol_step; //每一步的最大误差
    std::vector<Eigen::MatrixXd> A;      //系数矩阵，xyz三个方向有三个
    std::vector<Eigen::VectorXd> b;      //已知位移向量，xyz三个方向有三个
    std::vector<Eigen::VectorXd> coeff;   //插值系数，xyz三个方向有三个

    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> supp_pos;   //连续缓冲（行主序）用以优化效率 Ns*3:(x y z)
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> coeff_mat;   //Ns*3:(c0 c1 c2)


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

    // 直接用传入的所有候选点构建插值系统（不贪心）
    // reg 为可选对角正则（如 1e-12），默认 0
    void BuildAll(const std::vector<Node>& candidates, const State& S, double reg = 0.0);

    // 使用已准备好的 external_suppoints 构建（不贪心）
    void BuildAllFromExternal(const State& S, double reg = 0.0);

    //贪心算法误差输出至文件查看
    void write_greedy_tol(const std::string& filename);

    //从现有的suppoints/coeff 构建紧凑缓冲，用以优化效率
    void rebuild_compact_buffers_from_current();
private:

    // RBFInterpolator 类里（private 增加）
    bool trace_ = false;
    std::string trace_tag_;

    static inline Eigen::VectorXd robust_symmetric_solve(Eigen::MatrixXd& A,
        const Eigen::VectorXd& b,
        double* out_added_reg = nullptr)
    {
        // 先尝试 LLT
        {
            Eigen::LLT<Eigen::MatrixXd> llt(A);
            if (llt.info() == Eigen::Success) {
                if (out_added_reg) *out_added_reg = 0.0;
                return llt.solve(b);
            }
        }
        // 自适应极小对角正则再试 LLT（nugget）
        const double diag_mean = A.diagonal().cwiseAbs().mean();
        const double eps_base = std::max(1.0, diag_mean);
        double       lambda = 1e-12 * eps_base; // 先加一个很小的
        for (int trial = 0; trial < 2; ++trial) {
            Eigen::MatrixXd A2 = A;
            A2.diagonal().array() += lambda;
            Eigen::LLT<Eigen::MatrixXd> llt2(A2);
            if (llt2.info() == Eigen::Success) {
                if (out_added_reg) *out_added_reg = lambda;
                return llt2.solve(b);
            }
            lambda *= 1e3; // 再强一点
        }
        // 兜底：LDLT（带主元）——对称不定/半正定更稳健
        if (out_added_reg) *out_added_reg = -1.0; // 标记走了 LDLT
        Eigen::LDLT<Eigen::MatrixXd> ldlt;
        ldlt.compute(A);
        return ldlt.solve(b);
    }

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
    void calculate_deform_DRRBF(Node &inode, double d_r2omega1, double d_r2omega2, double D, const State &S) const;
    
};

// 重载1：使用每块的 unique_bndry 作为“本块的核心候选”，再用 block_D 扫入其它块内/边界点
void set_block_rbf(std::vector<DeformCalculator>& block_rbf,
                   const std::vector<std::vector<Node>>& unique_bndry,
                   const std::vector<mesh_block>& blocks);

// 重载2：仅依据各块的 internal / boundary + block_D扫入的其他块内/边界点作为“本块的核心候选”
void set_block_rbf(std::vector<DeformCalculator>& block_rbf, const std::vector<mesh_block>& blocks);

