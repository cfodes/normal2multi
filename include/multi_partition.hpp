#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <metis.h>
#include <cassert>
#include <utility>
#include <cstddef>
#include "geometry.hpp"
#include "distance.hpp"
#include "block.hpp"          // mesh_block / Reset_blocks
#include "rbf.hpp"            // RBFInterpolator / DeformCalculator / free set_block_rbf(...)
#include "metis_block.hpp"
#include "State.hpp"

// Holds per-level timing statistics for profiling
struct LevelTiming {
    double preprocess_ms = 0.0;
    double build_ms = 0.0;
    double compute_ms = 0.0;
    double update_ms = 0.0;
    double distance_ms = 0.0;
    double apply_ms = 0.0;
    double max_block_D = 0.0;
};

// Extra diagnostics collected in test mode
struct BlockTestInfo {
    int block_id = -1;
    std::size_t internal_points = 0;
    std::size_t candidate_points = 0;
    std::size_t support_points = 0;
    double block_D = 0.0;
};

// 多级 METIS 分区 + 分组 RBF 主控类
class multi_partition {
public:
    // 传入每层分区数，例如 {4,3,2}
    explicit multi_partition(const std::vector<idx_t>& parts_per_level);

    // 设置节点等级映射，用于加速网格DRRBF变形
    void set_m_nd2wall_lvl_mp(const std::unordered_map<int, int>& nd2wall_lvl);

    // 对边界单元 elems 和 wall_nodes 进行多级分区
    void divide_wall(const std::vector<element>& elems,
                     const std::vector<Node>& wall_nodes);

    // 分级 RBF 主流程
    // tol_steps.size() == levels_
    void multi_partition_rbf_algorithm(const std::vector<double>& tol_steps,
                                       const State& S,
                                       std::vector<Node>& wall_nodes,
                                       std::vector<Node>& all_nodes_coords);

    // 查询
    std::vector<mesh_block>&       get_level_blocks(size_t level);
    const std::vector<mesh_block>& get_level_blocks(size_t level) const;

    std::vector<Node>&             get_unique_boundary(size_t level, size_t block_index);
    const std::vector<Node>&       get_unique_boundary(size_t level, size_t block_index) const;

    size_t get_total_levels() const { return levels_; }

    // 获取在第 lvl 层，一个 unique 边界节点 id 所在的所有分区编号
    const std::vector<int>& query_unique_bndry2blkid(size_t lvl, int nd_id) const;

    // Accessors for the latest reports (timing + optional block diagnostics)
    const std::vector<LevelTiming>& last_timing_report() const { return timing_report_; }
    const std::vector<std::vector<BlockTestInfo>>& block_info_report() const { return block_info_report_; }
    // Enable or disable collection of block-level test info
    void set_collect_test_info(bool flag) { collect_test_info_ = flag; }
    // Control whether mid-level partitions build RBF systems via greedy selection
    void set_greedy_intermediate(bool flag) { use_greedy_intermediate_ = flag; }

private:
    // ===== 成员数据 =====
    size_t levels_;                               // 层数
    std::vector<idx_t> parts_per_level_;          // 每层分区数

    // 每层分区块
    std::vector<std::vector<mesh_block>> blocks_per_level_;
    // 每层（按父块）收集的“该父块的 children 之间共享的 unique 边界节点”
    std::vector<std::vector<std::vector<Node>>> unique_boundary_per_level_;

    // 每层每块：RBFInterpolator 的实际存储（保证生命周期）
    std::vector<std::vector<RBFInterpolator>> rbf_systems_per_level_;
    // 每层每块：DeformCalculator（引用上面的 RBFInterpolator）
    std::vector<std::vector<DeformCalculator>> block_rbf_per_level_;

    // 已经进入“静边界集合”的 unique 边界点
    std::unordered_set<int> global_unique_bndry_set_;
    std::vector<int>        global_unique_bndry_id_;
    GridBTree<int, Point<double>> unique_bndry_tree_;

    // wall id <-> index
    std::unordered_map<int, int> wall_id2nodes_;
    // 节点等级映射，用于加速网格DRRBF变形
    std::unordered_map<int, int> m_nd2wall_lvl_mp_;

    // 前置构建：unique 边界到分区编号的映射（按层存）
    std::vector<std::unordered_map<int, std::vector<int>>> unique_bndry2blkid_per_lvl_;

    // 最近一次运行的计时结果
    std::vector<LevelTiming> timing_report_;
    std::vector<std::vector<BlockTestInfo>> block_info_report_; // populated only when test info is requested
    bool collect_test_info_ = false;                            // toggled by callers that need detailed diagnostics
    bool use_greedy_intermediate_ = false;                      // 是否在非末级使用贪心法构建块内RBF

private:
    // ===== 内部工具 =====
    static void set_blk_id(std::vector<mesh_block>& blks);

    // 构建 wall 节点 id -> 索引
    static void set_wall_id2nodes(std::unordered_map<int, int>& wall_id2nd,
                                  const std::vector<Node>& wall_nodes);

    // 基于 wall 初始坐标 + df 得到理论“准确墙面”位置（s*）
    static std::vector<Node> set_accurate_wall(const std::vector<Node>& wall_0);

    // 用 all_nodes_coords（按 id 排序）更新 wall_prev 的坐标，并重置待变形量为 s* - s_k
    void update_wall_prev(std::vector<Node>& wall_prev,
                          const std::vector<Node>& wall_accurate,
                          const std::vector<Node>& all_nodes_coords);

    // 将第 lvl 层的 unique 边界点更新为 wall_prev 中的最新位置/df
    void update_unique_bndry_nodes(const std::vector<Node>& wall_prev, int lvl);

    // 用 wall_prev 更新本层 blocks 的节点待变形量，并重建块的 kd-tree（并计算 block_D）
    double set_blocks_df_and_tree(const std::vector<Node>& wall_prev, int lvl);

    // 第 0 层：建立全局 RBF 系统（用本层收集到的 unique 边界点作为候选集）
    void build_single_rbf_system(double tol, const State& S, size_t lvl);

    // 第 1 层及以上：建立每块 RBF 系统
    void build_multiple_rbf_systems(double tol, const State& S, size_t lvl);

    // 第 0 层：纯 RBF
    void apply_simple_rbf_deformation(std::vector<Node>& coords,
                                      const State& S, size_t lvl);

    // 第 ≥1 层：DRRBF（考虑动/静边界距离）
    void apply_drrbf_deformation(std::vector<Node>& coords,
                                 const State& S, size_t lvl,
                                 LevelTiming& timing_row);

    // 组建“历史 unique 边界” kd-tree（供静边界距离查询）
    void set_unique_bndry_tree(size_t lvl,
                               GridBTree<int, Point<double>>& tree,
                               const std::vector<Node>& all_nodes);

    // 收集前一层的 unique 边界节点 id（去重）
    static void collect_unique_bndry_nodes(
        size_t lvl,
        const std::vector<std::vector<std::vector<Node>>>& unique_boundary_per_level,
        std::unordered_set<int>& bndry_set,
        std::vector<int>& bndry_id);

    // 建立“本层 unique 边界 id → 覆盖的分区 id 列表”的映射
    void build_unique_bndry2blkid_for_lvl(size_t lvl);

    // 查询某个点到“动边界”和“静边界”的距离
    void find_moving_and_static_bndry(const Node& ind,
                                  const std::vector<mesh_block>& blocks,
                                  int& i_id,
                                  double& i_d_r2omega1,
                                  double& i_d_r2omega2) const;


};
