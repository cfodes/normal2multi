#pragma once
#include "State.hpp"
#include "rbf.hpp"
#include "deformation.hpp"
#include <metis.h>
#include <string>
#include <vector>

class RBFTest
{
public:
    // 完整运行流程：读入 → 变形 → 输出
    static void run_global_test(const std::string& input_file,
                                const std::string& output_file,
                                double tol,
                                const State& S);
};

struct PartitionBatchCase
{
    std::string run_name;              // 个案名称，用于生成文件名
    std::vector<idx_t> parts;          // 分组设置
    std::vector<double> tolerances;    // 误差容限
};

class MultiPartitionBatch
{
public:
    MultiPartitionBatch(std::string test_name,
                        std::string input_file,
                        std::string output_root = "output");

    // 添加一条测试配置
    void add_case(const std::string& run_name,
                  const std::vector<idx_t>& parts,
                  const std::vector<double>& tolerances);

    // 清空已添加的配置
    void clear_cases();

    // 调整输出根目录（默认 output）
    void set_output_root(const std::string& output_root);

    // 运行所有配置，将结果写入 output_root/test_name/run_name.su2
    void run_all() const;

    const std::string& test_name() const { return test_name_; }
    const std::string& input_file() const { return input_file_; }
    const std::vector<PartitionBatchCase>& cases() const { return cases_; }

private:
    std::string test_name_;
    std::string input_file_;
    std::string output_root_;
    std::vector<PartitionBatchCase> cases_;
};
