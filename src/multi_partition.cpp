#include "multi_partition.hpp"
#include "deformation.hpp"
#include "rbf.hpp"
#include <algorithm>
#include <iostream>
#include <chrono>
using namespace std::chrono;

// ===== 构造 =====
multi_partition::multi_partition(const std::vector<idx_t> &parts_per_level)
    : levels_(parts_per_level.size()), parts_per_level_(parts_per_level) {}

// 设置节点等级映射，用于加速网格DRRBF变形
void multi_partition::set_m_nd2wall_lvl_mp(
    const std::unordered_map<int, int> &nd2wall_lvl) {
  m_nd2wall_lvl_mp_ = nd2wall_lvl;
}

// 对边界单元 elems 和 wall_nodes 进行多级分区
// elems: 物面单元数组
// wall_nodes: 物面节点数组
void multi_partition::divide_wall(const std::vector<element> &elems,
                                  const std::vector<Node> &wall_nodes) {
  // 预分配各层容器
  blocks_per_level_.clear();
  unique_boundary_per_level_.clear();
  rbf_systems_per_level_.clear();
  block_rbf_per_level_.clear();

  blocks_per_level_.reserve(levels_);
  unique_boundary_per_level_.reserve(levels_);
  rbf_systems_per_level_.resize(levels_);
  block_rbf_per_level_.resize(levels_);
  unique_bndry2blkid_per_lvl_.resize(levels_);

  // 根分区：一个 mesh_block 包含所有元素与节点
  std::vector<mesh_block> current;
  current.reserve(1);
  {
    mesh_block root;
    root.internal_elements = elems;
    root.m_nodes = wall_nodes;
    current.push_back(std::move(root));
  }
  // 复用容器，减少分配
  std::vector<mesh_block> next_blocks;
  std::vector<std::vector<Node>> level_unique_boundary;

  // 逐层分区
  for (size_t lvl = 0; lvl < levels_; ++lvl) {
    next_blocks.clear();
    level_unique_boundary.clear();
    next_blocks.reserve(current.size() * parts_per_level_[lvl]);
    level_unique_boundary.reserve(current.size());

    for (auto &blk : current) {
      // 构造 metis_block
      metis_block mblk(blk.internal_elements, blk.m_nodes,
                       parts_per_level_[lvl]);
      mblk.partitionMeshDual();

      // 调用 generate_mesh_blocks，获取子 mesh_block 和它们的 unique boundary
      std::vector<mesh_block> children;
      std::vector<Node> unique_bndry; // 该父块下，children 之间的共享边界
      children.reserve(parts_per_level_[lvl]);

      mblk.generate_mesh_blocks(blk.internal_elements, blk.m_nodes, children,
                                unique_bndry);

      // 除了第一次整体的分区以外，后续的分区在 generate_mesh_blocks
      // 之后， 其分出的新子分区有可能失去一部分边界节点。
      // 原因：generate_mesh_blocks 接收的只是当前大分区的信息，
      //       它并不知道该大分区和其他大分区之间共享的边界节点。
      // 因此这里需要额外的补充操作：
      //
      // 遍历父分区 (blk) 的 boundary_nodes，检查哪些节点也属于子分区
      // (iblk) 的 m_nodes：
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
      if (lvl != 0) {
        for (auto &cblk : children) {
          for (
              const auto &nd :
              blk.boundary_nodes) // blk是上一次分区的结果，每个大分区的边界节点包含了其和其他大分区共享的边界节点
          {
            if (cblk.m_nodes_id.find(nd.id) !=
                    cblk.m_nodes_id
                        .end() && // ind是大分区的边界节点，如果他属于小分区的成员节点
                cblk.bndry_nds_id.find(nd.id) ==
                    cblk.bndry_nds_id
                        .end()) // 且不属于小分区的边界节点，则插入该节点
            {
              cblk.boundary_nodes.push_back(nd); // 插入边界节点
              cblk.bndry_nds_id.insert(nd.id);   // 插入边界节点的id
            }
          }
        }
      }

      // 移动子分区及其唯一边界
      next_blocks.insert(next_blocks.end(),
                         std::make_move_iterator(children.begin()),
                         std::make_move_iterator(children.end()));
      level_unique_boundary.push_back(std::move(unique_bndry));
    }

    // 保存本层结果
    blocks_per_level_.push_back(std::move(next_blocks));
    set_blk_id(blocks_per_level_[lvl]); // 重新调整分区编号
    unique_boundary_per_level_.push_back(std::move(level_unique_boundary));

    // 构建unique_bndry_nodes的id和分区编号的查找集
    build_unique_bndry2blkid_for_lvl(lvl);

    // 为下一层准备：交换以复用内存，避免拷贝
    current.swap(blocks_per_level_.back());
  }
}

// 分组 RBF 算法主入口：
// tol_steps: 每层的误差容限数组，长度 == levels
// wall_nodes: 全局物面节点数组，包含最终的已知变形结果
// all_nodes_coords: 全局所有待变形节点坐标
void multi_partition::multi_partition_rbf_algorithm(
    const std::vector<double> &tol_steps, const State &S,
    std::vector<Node> &wall_nodes, std::vector<Node> &all_nodes_coords) {
  // 检查误差设置的级数与分区的级数是否一致
  assert(tol_steps.size() == levels_ && "tol_steps size must match levels");

  set_wall_id2nodes(wall_id2nodes_, wall_nodes); // 构建物面的索引和id的查找集
  std::vector<Node> wall_accurate = set_accurate_wall(
      wall_nodes); // 根据物面节点初始位置和待变形量计算准确的变形结果
  std::vector<Node> wall_prev =
      wall_nodes; // 当前物面结果，第 0 层保持原始 wall

  std::vector<LevelTiming> timing_report(levels_);
  if (collect_test_info_) {
    // Pre-size per-level block diagnostic container when test mode is on
    block_info_report_.assign(levels_, {});
  } else {
    block_info_report_.clear();
  }

  for (size_t lvl = 0; lvl < levels_; ++lvl) {
    const double tol = tol_steps[lvl];

    //每一步迭代的预处理
    const auto t_pre_start = high_resolution_clock::now();  //预处理计时开始
    if (lvl != 0) {
      // 根据更新后的all_node_coord和准确的物面变形结果更新wall_prev的节点坐标和待变形量
      update_wall_prev(wall_prev, wall_accurate, all_nodes_coords);
      // 根据更新后的wall_prev重新设置分区内部的节点坐标,
      // 变形量由后续的set_internal_nodes_df处理
      Reset_blocks(blocks_per_level_[lvl], wall_prev, wall_id2nodes_);
      // 更新已经用过的unique_bndry的位置，这些点已经准确变形，无需再变
      update_unique_bndry_nodes(wall_prev, static_cast<int>(lvl));
      // 根据更新后的wall_prev重新设置分区内部的节点待变形量，并且计算了分区的block_D!!!
      const double max_block_D = set_blocks_df_and_tree(wall_prev, static_cast<int>(lvl));
      timing_report[lvl].max_block_D = max_block_D;
      // 建立“历史 unique 边界”的 kd-tree（静边界距离要用）
      set_unique_bndry_tree(lvl, unique_bndry_tree_, all_nodes_coords);

      //预处理计时结束
      const auto t_pre_end = high_resolution_clock::now();
      const double ms = duration<double, std::milli>(t_pre_end - t_pre_start).count();
      timing_report[lvl].preprocess_ms = ms;
      std::cout << "[lvl=" << lvl << "] preprocess time: " << ms << " ms\n";
    }

    // 步骤 1: 构建 RBF 插值系统（第0层只有1个系统，从第1层开始每个块一个）
    const auto t_build_start = high_resolution_clock::now(); //插值系统计时开始
    if (lvl == 0) {
      build_single_rbf_system(tol, S, lvl);
    } else { // 构建分区RBF系统
      build_multiple_rbf_systems(tol, S, lvl);
    }
    //插值系统计时结束
    const auto t_build_end = high_resolution_clock::now();
    const double build_ms = duration<double, std::milli>(t_build_end - t_build_start).count();
    timing_report[lvl].build_ms = build_ms;
    std::cout << "[lvl=" << lvl << "] build RBF system time: "
        << build_ms
        << " ms\n";

    // 步骤 2: 执行节点变形量计算
    const auto t_df_start = high_resolution_clock::now();  //变形量计算计时开始
    if (lvl == 0) {
      apply_simple_rbf_deformation(all_nodes_coords, S, lvl);
    } else // 第 ≥1 层：DRRBF（考虑动/静边界距离）
    {
      apply_drrbf_deformation(all_nodes_coords, S, lvl, timing_report[lvl]);
    }
    //变形量计算计时结束
    const auto t_df_end = high_resolution_clock::now();
    const double compute_ms = duration<double, std::milli>(t_df_end - t_df_start).count();
    timing_report[lvl].compute_ms = compute_ms;
    std::cout << "[lvl=" << lvl << "] compute df time: "
        << compute_ms
        << " ms\n";

    // 3) 计算变形后的节点坐标 s_{k+1} = s_k + df
    const auto t_update_start = high_resolution_clock::now();  //计算节点坐标计时开始
    calculate_deformed_coordinates(all_nodes_coords);
    //计算节点坐标计时结束
    const auto t_update_end = high_resolution_clock::now();
    const double update_ms = duration<double, std::milli>(t_update_end - t_update_start).count();
    timing_report[lvl].update_ms = update_ms;
    std::cout << "[lvl=" << lvl << "] update coordinates time: "
        << update_ms
        << " ms\n";
  }

  timing_report_.resize(timing_report.size());
  for (size_t lvl = 0; lvl < timing_report.size(); ++lvl) {
    timing_report_[lvl].preprocess_ms = timing_report[lvl].preprocess_ms;
    timing_report_[lvl].build_ms = timing_report[lvl].build_ms;
    timing_report_[lvl].compute_ms = timing_report[lvl].compute_ms;
    timing_report_[lvl].update_ms = timing_report[lvl].update_ms;
    timing_report_[lvl].distance_ms = timing_report[lvl].distance_ms;
    timing_report_[lvl].apply_ms = timing_report[lvl].apply_ms;
    timing_report_[lvl].max_block_D = timing_report[lvl].max_block_D;
  }
}

// ===== 查询 =====
// 获取某层的所有分区
std::vector<mesh_block> &multi_partition::get_level_blocks(size_t level) {
  return blocks_per_level_.at(level);
}

// 获取某层的所有分区（只读）
const std::vector<mesh_block> &
multi_partition::get_level_blocks(size_t level) const {
  return blocks_per_level_.at(level);
}

// 获取某层某分区的唯一边界节点
std::vector<Node> &multi_partition::get_unique_boundary(size_t level,
                                                        size_t block_index) {
  return unique_boundary_per_level_.at(level).at(block_index);
}

// 获取某层某分区的唯一边界节点（只读）
const std::vector<Node> &
multi_partition::get_unique_boundary(size_t level, size_t block_index) const {
  return unique_boundary_per_level_.at(level).at(block_index);
}

// 获取unique_bndry节点id对应的分区编号列表
const std::vector<int> &
multi_partition::query_unique_bndry2blkid(size_t lvl, int nd_id) const 
{
  assert(lvl > 0 && lvl <= levels_ && "lvl must be in [1, levels_]");
  const auto &mp = unique_bndry2blkid_per_lvl_.at(lvl - 1);
  auto it = mp.find(nd_id);
  if (it == mp.end()) {
    std::cerr << "[Error] unique_bndry node " << nd_id
              << " has no corresponding block at level " << lvl << '\n';
    std::abort();
  }
  return it->second;
}

// ===== 内部工具 =====
// 构建物面的索引和id查找集
void multi_partition::set_wall_id2nodes(std::unordered_map<int, int>& wall_id2nd,
                       const std::vector<Node>& wall_nodes) 
{
  wall_id2nd.clear();  
  wall_id2nd.reserve(wall_nodes.size());
  for (size_t i = 0; i < wall_nodes.size(); ++i) 
  {
    wall_id2nd[wall_nodes[i].id] = static_cast<int>(i);
  }
}

void multi_partition::set_blk_id(std::vector<mesh_block> &blk)
// 由于generate_mesh_blocks只能根据其接受的分区设置新分出的小分区的id
// 此函数将局部的id转换成全局的id
// 由于只对一个物面进行分区，故按顺序分配id即可
// blk为blocks_per_level[lvl]的mesh_block的数组
{
  for (int i = 0; i < blk.size(); ++i) 
  {
    blk[i].block_id = i;
  }
}

std::vector<Node> multi_partition::set_accurate_wall(const std::vector<Node>& wall_0)
// 计算准确的物面变形后的物面节点
// wall_0：外部传入的_wall_nodes，需有节点原始位置和待变形量
{
    std::vector<Node> wall_temp(wall_0);
    const int dim = static_cast<int>(wall_0[0].point.size());    //获取点的维度
    for (int i = 0; i < static_cast<int>(wall_0.size()); ++i)    //遍历物面数组
    {
        for (int k = 0; k < dim; ++k)
            wall_temp[i].point[k] = wall_0[i].point[k] + wall_0[i].df[k];   //变形后的节点位置为原来的位置加上设定的位移量
        wall_temp[i].id = wall_0[i].id;   //保留id编号
    }
    return wall_temp;
}

void multi_partition::update_wall_prev(std::vector<Node>& wall_prev,
                                       const std::vector<Node>& wall_accurate,
                                       const std::vector<Node>& all_nodes_coords)
//根据更新后的all_nodes_coords修改wall_prev的节点坐标和待变形量
//all_node_coord 是按照节点id排序的全局节点数组，每次循环的RBF变形后，其节点坐标和df都发生改变
//wall_prev是当前的物面节点数组
//wall_0是上一步的物面节点，用以提供上一步的待变形量
{
    const int dim = static_cast<int>(wall_prev[0].df.size());   //获取df的维度
    for (auto& inode : wall_prev) 
    {
        // all_nodes_coords 下标等于 id（你的工程就是这样组织的）
        assert(inode.id == all_nodes_coords[inode.id].id);   //确保在all_node_coord里索引的节点和当前的物面节点一致，节点编号是唯一的
        inode = all_nodes_coords[inode.id];
        auto it = wall_id2nodes_.find(inode.id);   //确保没有在wall里查找wall以外的节点
        assert(it != wall_id2nodes_.end());

        for (int k = 0; k < dim; ++k) {
            inode.df[k] = wall_accurate[it->second].point[k] - inode.point[k];   //新的待变形量为准确的节点位置减去已变形的节点位置：s*-sk
        }
    }
}

void multi_partition::update_unique_bndry_nodes(const std::vector<Node> &wall_prev, int lvl)
// 根据更新后的wall_prev更新unique_boundary_nodes内部的节点坐标和待变形量
// wall_prev: 更新后的物面节点
// lvl: 级数
{
  for (auto &by_parent : unique_boundary_per_level_.at(lvl))  //by_parent 每个分区的共享边界点数组
  {
    for (auto &inode : by_parent)   //inode 共享边界点
    {
      auto it = wall_id2nodes_.find(inode.id);
      assert(it != wall_id2nodes_.end());     //确保没有在wall里查找wall以外的节点
      inode = wall_prev[it->second];          //将wall_prev对应的节点的坐标和df传给unique_bndry_nodes
    } 
  }
}

double multi_partition::set_blocks_df_and_tree(const std::vector<Node> &wall_prev,
                                               int lvl)
// 计算当前级数下的blocks的节点的df并且构建相应的二叉树
// wall_prev: 更新后的物面节点
// lvl: 级数
{
  auto &blks = blocks_per_level_.at(lvl);    //blks：当前层级的小分区数组
  double max_block_D = 0.0;
  if (collect_test_info_) {
    // Prepare the current level's diagnostic storage
    if (block_info_report_.size() <= static_cast<std::size_t>(lvl)) {
      block_info_report_.resize(levels_);
    }
    auto &lvl_report = block_info_report_[lvl];
    lvl_report.clear();
    lvl_report.reserve(blks.size());
    for (auto &b : blks) 
    {
      b.set_blk_nodes_df(wall_prev, wall_id2nodes_);     //根据已经更新的wall_prev设定blocks的节点待变形量
      b.set_block_tree();
      max_block_D = std::max(max_block_D, b.block_D);
      lvl_report.push_back(BlockTestInfo{b.block_id,
                                         b.internal_nodes.size(),
                                         0,
                                         0,
                                         b.block_D});
    }
  } else {
    for (auto &b : blks) 
    {
      b.set_blk_nodes_df(wall_prev, wall_id2nodes_);     //根据已经更新的wall_prev设定blocks的节点待变形量
      b.set_block_tree();
      max_block_D = std::max(max_block_D, b.block_D);
    }
  }
  return max_block_D;
}

// ===== 第 0 层：建立单个 RBF 系统 =====
void multi_partition::build_single_rbf_system(double tol, const State &S,
                                              size_t lvl)
// wall_final: 最終的物面变形结果
// tol: 误差容限
{
  // 我们用第 0 层收集到的 unique 边界点作为候选集（位于unique_boundary_per_level_[lvl][0]）
  rbf_systems_per_level_[lvl].resize(1);
  block_rbf_per_level_[lvl].clear();
  block_rbf_per_level_[lvl].emplace_back(rbf_systems_per_level_[lvl][0]);

  //// 直接调用“未预设 external_suppoints”的 Greedy 版本
  //rbf_systems_per_level_[lvl][0].Greedy_algorithm(
  //    unique_boundary_per_level_.at(lvl).at(0), tol, S);

  // ---- 直接用所有候选点一次性构建（不使用贪心）----
  constexpr double kReg = 1e-12; // 矩阵迭代求解误差
  rbf_systems_per_level_[lvl][0].BuildAll(
      unique_boundary_per_level_.at(lvl).at(0), S, kReg);
}

// ===== 第 ≥1 层：每块一个 RBF 系统 =====
void multi_partition::build_multiple_rbf_systems(double tol, const State& S, size_t lvl)
{
    auto& blks = blocks_per_level_.at(lvl);    //获取当前层的分区数组
    const int nb = static_cast<int>(blks.size());                      //当前层的分区数

    // 先准备 RBF 存储和 Calculator
    rbf_systems_per_level_[lvl].assign(nb, RBFInterpolator{});
    block_rbf_per_level_[lvl].clear();
    block_rbf_per_level_[lvl].reserve(nb);
    for (int i = 0; i < nb; ++i) 
    {
        block_rbf_per_level_[lvl].emplace_back(rbf_systems_per_level_[lvl][i]);
    }

    // ★ 根据是不是最后一层，调用不同的 set_block_rbf
    if (lvl < levels_ - 1) {
        // 非最后一层：用 unique_bndry + block_D
        set_block_rbf(block_rbf_per_level_[lvl], unique_boundary_per_level_[lvl], blks);
    } else {
        // 最后一层：用 internal_nodes + block_D
        set_block_rbf(block_rbf_per_level_[lvl], blks);
    }

    if (collect_test_info_ && block_info_report_.size() > static_cast<std::size_t>(lvl)) {
        // Record candidate pool size per block for diagnostics
        auto &lvl_report = block_info_report_[lvl];
        for (int i = 0; i < nb && i < static_cast<int>(lvl_report.size()); ++i) {
            lvl_report[i].candidate_points = rbf_systems_per_level_[lvl][i].external_suppoints.size();
        }
    }

    //// 每块做一次 Greedy_algorithm（使用 external_suppoints）
    //for (int i = 0; i < nb; ++i) {
    //    rbf_systems_per_level_[lvl][i].Greedy_algorithm(tol, S);
    //}

    // ---- 每块直接一次性构建（使用已填充的 external_suppoints）----
    constexpr double kReg = 1e-12; // 矩阵迭代求解误差
    for (int i = 0; i < nb; ++i) {
        rbf_systems_per_level_[lvl][i].BuildAllFromExternal(S, kReg);
    }

    if (collect_test_info_ && block_info_report_.size() > static_cast<std::size_t>(lvl)) {
        // Record resulting support count per block for diagnostics
        auto &lvl_report = block_info_report_[lvl];
        for (int i = 0; i < nb && i < static_cast<int>(lvl_report.size()); ++i) {
            lvl_report[i].support_points = rbf_systems_per_level_[lvl][i].suppoints.size();
        }
    }
}

// ===== 第 0 层：纯 RBF 变形 =====
void multi_partition::apply_simple_rbf_deformation(std::vector<Node>& coords,
                                                   const State& S, size_t lvl)
{
    // 第 0 层只有一个全局系统
    block_rbf_per_level_[lvl][0].calculate_deform(coords, S);
}

// ===== 第 ≥1 层：DRRBF 变形 =====
void multi_partition::apply_drrbf_deformation(std::vector<Node>& coords,
                                              const State& S, size_t lvl,
                                              LevelTiming& timing_row)
{
    const auto& blks = blocks_per_level_.at(lvl);   // 当前层的所有分区
    auto& calculators = block_rbf_per_level_.at(lvl);  // 当前层的所有 RBF 计算器

    struct DistanceRecord {
        int moving_blk_id = -1;
        double d_r2omega1 = 0.0;
        double d_r2omega2 = 0.0;
        bool skip = false;
    };

    std::vector<DistanceRecord> distance_records(coords.size());

    const auto t_dist_start = high_resolution_clock::now();
    for (size_t idx = 0; idx < coords.size(); ++idx) {
        const auto& inode = coords[idx];
        DistanceRecord record;

        if (global_unique_bndry_set_.find(inode.id) != global_unique_bndry_set_.end()) {
            record.skip = true;
            distance_records[idx] = record;
            continue;
        }

        find_moving_and_static_bndry(
            inode,
            blks,
            record.moving_blk_id,
            record.d_r2omega1,
            record.d_r2omega2);

        distance_records[idx] = record;
    }
    const auto t_dist_end = high_resolution_clock::now();
    const double dist_ms = duration<double, std::milli>(t_dist_end - t_dist_start).count();
    timing_row.distance_ms = dist_ms;
    std::cout << "[lvl=" << lvl << "] distance query time: "
              << dist_ms
              << " ms\n";

    const auto t_apply_start = high_resolution_clock::now();
    for (size_t idx = 0; idx < coords.size(); ++idx) {
        auto& inode = coords[idx];
        const auto& record = distance_records[idx];

        if (record.skip || record.moving_blk_id < 0) {
            inode.df.setZero();
            continue;
        }

        calculators[record.moving_blk_id].calculate_deform_DRRBF(
            inode,
            record.d_r2omega1,
            record.d_r2omega2,
            blks[record.moving_blk_id].block_D,
            S);
    }
    const auto t_apply_end = high_resolution_clock::now();
    const double apply_ms = duration<double, std::milli>(t_apply_end - t_apply_start).count();
    timing_row.apply_ms = apply_ms;
    std::cout << "[lvl=" << lvl << "] DRRBF apply time: "
              << apply_ms
              << " ms\n";
}



// ===== 历史 unique 边界 kd-tree =====
void multi_partition::set_unique_bndry_tree(size_t lvl,
                                            GridBTree<int, Point<double>>& tree,
                                            const std::vector<Node>& all_nodes)
//将之前的unique_bndry点构建成一个二叉树用以查找到unique_bndry点的最小距离
//lvl：当前变形的级数
//unique_bndry_tree_temp: 用unique_bndry点构建的二叉树
{
    // 收集“上一层”的 unique 边界 id
    collect_unique_bndry_nodes(lvl, unique_boundary_per_level_,
                               global_unique_bndry_set_,
                               global_unique_bndry_id_);
    tree.setGridSize(static_cast<int>(global_unique_bndry_id_.size()));    //重新调整二叉树的大小
    for (int i = 0; i < static_cast<int>(global_unique_bndry_id_.size()); ++i)  //遍历global_unique_bndry_id的所有索引
    {
        int id = global_unique_bndry_id_[i];
        assert(id == all_nodes[id].id);    //确保点的id和它的索引是一致的
        tree.setGrid(i, id, all_nodes[id].point);   //将点放入二叉树内
    }
    tree.construct();
}

void multi_partition::collect_unique_bndry_nodes(
    size_t lvl,
    const std::vector<std::vector<std::vector<Node>>>& unique_boundary_per_level,
    std::unordered_set<int>& bndry_set,
    std::vector<int>& bndry_id)
//将之前所有分区存在的unique_bndry节点所对应的索引存放到bndry_id
//bndry_set确保点是唯一的
//只需要不断地将上一级的存入即可
//该函数确保了drrbf查找到静边界距离的时候，不会漏掉unique_bndry的点
//lvl：变形当前的级数，其上一级是lvl-1
{
    // 把 lvl-1 层所有 parent 的 unique 边界节点 id 汇总（去重）
    const auto& per_parent = unique_boundary_per_level.at(lvl - 1);
    for (const auto& one_parent_list : per_parent) 
    {
        for (const auto& nd : one_parent_list) 
        {
            if (bndry_set.insert(nd.id).second)      //如果节点未被重复，则插入
            {
                bndry_id.push_back(nd.id);    
            }
        }
    }
}

// ===== 建立“unique 边界 id → 分区 id 列表”的映射 =====
void multi_partition::build_unique_bndry2blkid_for_lvl(size_t lvl)
//建立每个lvl对应的unique_bndry_nodes和blk_id的查找集
//输入: lvl: 层级
{
    auto& mp = unique_bndry2blkid_per_lvl_[lvl];    //取当前层级的unique_bndry_nodes和blk_id的查找集
    mp.clear();
    for (const auto& blk : blocks_per_level_.at(lvl))
    {
        for (int nd_id : blk.bndry_nds_id) 
        {
            mp[nd_id].push_back(blk.block_id);    //插入节点id和分区id的对应关系
        }
    }
}



void multi_partition::find_moving_and_static_bndry(
    const Node& ind,
    const std::vector<mesh_block>& blocks,
    int& i_id,
    double& i_d_r2omega1,
    double& i_d_r2omega2) const
// 查找到空间节点ind最近的动边界和静止边界
// ind: 待查询的空间节点
// blocks: 当前层的所有分区
// i_id: 输出，距离最近的动边界所在分区的 id
// i_d_r2omega1: 输出，距离最近的动边界的距离
// i_d_r2omega2: 输出，距离最近的静边界的距离
{
    // 遍历所有 block 的内部节点树，记录最小和次小的距离
    assert(!blocks.empty() && "find_moving_and_static_bndry: empty blocks");  // 确保分区非空
    double min_dist = 1.0e20;
    double second_min_dist = 1.0e20;
    int nearest_blk_id = -1;
    int key_temp = 0;

    for (const auto& blk : blocks) {
        double d_temp = 0.0;
        blk.block_tree.search(d_temp, key_temp, ind.point);

        if (d_temp < min_dist) {
            // 更新最小距离，同时把原来的最小距离下放为次小
            second_min_dist = min_dist;
            min_dist = d_temp;
            nearest_blk_id = blk.block_id;
        } else if (d_temp < second_min_dist) {
            // 更新次小距离
            second_min_dist = d_temp;
        }
    }

    // 动边界：离节点最近的 block
    assert(nearest_blk_id >= 0 && "nearest block not found");
    i_id = nearest_blk_id;
    i_d_r2omega1 = min_dist;

    // 静边界候选1：次近 block
    double candidate_static = second_min_dist;

    // 静边界候选2：global unique_bndry
    double d_unique = 1.0e20;
    if (!global_unique_bndry_id_.empty()) {
        unique_bndry_tree_.search(d_unique, key_temp, ind.point);
    }

    // 静边界取二者较小
    i_d_r2omega2 = std::min(candidate_static, d_unique);
}
