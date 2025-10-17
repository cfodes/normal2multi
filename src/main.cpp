#include "Test.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <exception>
#include <streambuf>

namespace {

// 简单的空设备输出缓冲，用于静默模式屏蔽 cout
struct NullBuffer : public std::streambuf {
    int overflow(int c) override { return traits_type::not_eof(c); }
};

// 支持的日志模式
enum class LogMode { Console, Silent };

static LogMode parse_log_mode(int argc, char** argv) {
    // 支持 --log console | --log silent | -q (silent)
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-q" || arg == "--quiet") return LogMode::Silent;
        if (arg == "--log") {
            if (i + 1 < argc) {
                std::string val = argv[++i];
                if (val == "console") return LogMode::Console;
                if (val == "silent") return LogMode::Silent;
            }
        } else if (arg.rfind("--log=", 0) == 0) {
            std::string val = arg.substr(6);
            if (val == "console") return LogMode::Console;
            if (val == "silent") return LogMode::Silent;
        }
    }
    return LogMode::Console; // 默认直接输出到控制台
}

struct GlobalExampleConfig {
    std::string input_file;
    std::string output_file;
    double tol;
};

MultiPartitionBatch make_buildall_example()
{
    MultiPartitionBatch batch("onera_m6_buildall",
                              "data/onera_m6_2M4.su2");

    batch.add_case("case_40_8_2",
                   {40, 8, 2},
                   {1e-3, 1e-3, 1e-3},
                   false);

    batch.add_case("case_20_16_2",
                   {20, 16, 2},
                   {1e-3, 1e-3, 1e-3},
                   false);

    batch.add_case("case_10_32_2",
                   {10, 32, 2},
                   {1e-3, 1e-3, 1e-3},
                   false);

    batch.add_case("case_2_4_5_2",
                   {2, 4, 5, 2},
                   {1e-3, 1e-3, 1e-3, 1e-3},
                   false);

    return batch;
}

MultiPartitionBatch make_greedy_example()
{
    MultiPartitionBatch batch("onera_m6_greedy",
                              "data/onera_m6_2M4.su2");

    batch.add_case("greedy_40_8_2",
                   {40, 8, 2},
                   {1e-4, 1e-4, 1e-3},
                   true);

    batch.add_case("greedy_20_16_2",
                   {20, 16, 2},
                   {1e-4, 1e-4, 1e-3},
                   true);

    batch.add_case("greedy_10_32_2",
                   {10, 32, 2},
                   {1e-4, 1e-4, 1e-3},
                   true);

    batch.add_case("greedy_2_4_5_2",
                   {2, 4, 5, 2},
                   {1e-4, 1e-4, 1e-4, 1e-3},
                   true);

    return batch;
}

GlobalExampleConfig make_global_example()
{
    return GlobalExampleConfig{
        // 可执行文件位于 build/，data 与 output 与 build 同级
        "../data/onera_m6_2M4.su2",
        "../output/onera_m6_global_result.su2",
        1e-3
    };
}

void run_batch_example(MultiPartitionBatch batch, bool with_test_info)
{
    if (with_test_info) {
        batch.run_all_with_test_info();
    } else {
        batch.run_all();
    }
}

void run_global_example(const GlobalExampleConfig& cfg)
{
    RBFTest::run_global_test(cfg.input_file,
                             cfg.output_file,
                             cfg.tol,
                             State{});
}

} // namespace

int main(int argc, char** argv)
{
    // 设置日志模式（控制窗口输出或静默）
    LogMode log_mode = parse_log_mode(argc, argv);
    NullBuffer null_buf;
    std::streambuf* original_cout_buf = nullptr;
    if (log_mode == LogMode::Silent) {
        original_cout_buf = std::cout.rdbuf(&null_buf);
    }

    // === 可切换的测试案例 ===
    const bool run_buildall_example = true;   // 分组RBF传统方法示例
    const bool run_greedy_example = true;     // 分组RBF贪心算法示例
    const bool run_global_example_case = true; // 全局 RBF 测试示例

    if (run_buildall_example) {
        run_batch_example(make_buildall_example(), false);
    }

    if (run_greedy_example) {
        run_batch_example(make_greedy_example(), true);
    }

    if (run_global_example_case) {
        run_global_example(make_global_example());
    }

    // 恢复 cout（静默模式）
    if (original_cout_buf) {
        std::cout.rdbuf(original_cout_buf);
    }

    return 0;
}
