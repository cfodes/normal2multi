#pragma once

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

    // Accessors used by batch runners after a call to run()
    const std::vector<LevelTiming>& get_last_timing_report() const { return last_timing_report_; }
    const std::vector<std::vector<BlockTestInfo>>& get_last_block_info_report() const { return last_block_info_report_; }
    // Switch on/off collection of extra test diagnostics
    void set_collect_test_info(bool flag) { collect_test_info_ = flag; }

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

    std::vector<LevelTiming> last_timing_report_;
    std::vector<std::vector<BlockTestInfo>> last_block_info_report_;
    bool collect_test_info_ = false; // keeps track of whether the next run() should gather test diagnostics
 
    // === 核心流程 ===
    void read_mesh();     // 读取网格文件
    void preprocess();    // 预处理，构建物面树和节点等级映射
    void run_multi_partition();   // 运行多级分区RBF算法
    void write_mesh();    // 写出变形后的网格文件
};
