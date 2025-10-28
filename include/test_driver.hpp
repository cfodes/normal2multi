#pragma once

#include <filesystem>
#include <metis.h>
#include <string>
#include <unordered_map>
#include <vector>
#include "State.hpp"
#include "multi_partition.hpp"


class TestDriver {
public:
    // 构造函数，初始化输入输出文件名、分区配置和误差容限
    TestDriver(const std::string& input_file,
               const std::string& output_file,
               const std::vector<idx_t>& parts_per_level,
               const std::vector<double>& tol_steps);

    // 主执行流程
    void run();

private:
    // === 输入输出 ===
    std::string input_file_;   // 输入网格文件名
    std::string output_file_;  // 输出网格文件名

    // === 算法配置 ===
    std::vector<idx_t> parts_per_level_;  // 每层的分区数
    std::vector<double> tol_steps_;       // 每层的误差容限

    // === 状态数据 ===
    State S_;                      // 状态参数
    std::unordered_map<int, int> nd2wall_lvl_;    // 节点id到物面等级的映射

    // === 核心流程 ===
    void read_mesh();     // 读取网格文件
    void preprocess();    // 预处理，构建物面树和节点等级映射
    void run_multi_partition();   // 运行多级分区RBF算法
    void write_mesh();    // 写出变形后的网格文件
    void write_reports(); // 输出计时与信息到文件

    // === 统计信息 ===
    struct BlockSummary {
        int block_id = -1;   //分区编号
        size_t candidate_points = 0;   // 该block的候选点数量
        size_t support_points = 0;     // 该block的支撑点数量
        double block_D = 0.0;  //分区的D
    };

    std::vector<multi_partition::LevelTiming> level_timings_;
    std::vector<multi_partition::LevelStatistics> level_statistics_;
    std::vector<std::vector<BlockSummary>> block_summaries_;
    double initial_D_ = 0.0;

    void write_time_report(const std::filesystem::path& time_path) const;
    void write_info_report(const std::filesystem::path& info_path) const;
};
