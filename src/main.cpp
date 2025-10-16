#include "Test.hpp"
#include "StreamCapture.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <exception>

namespace {

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
        "data/onera_m6_2M4.su2",
        "output/onera_m6_global_result.su2",
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

int main()
{
    StreamCapture capture(std::cout);

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

    capture.stop();
    capture.dump_to_stream(std::cout);

    const std::string log_file = "run.log";
    const std::vector<std::string> log_locations = {
        "../output/" + log_file,
        "output/" + log_file
    };

    bool saved = false;
    for (const auto& path : log_locations) {
        try {
            capture.save_to_file(path);
            saved = true;
            break;
        } catch (const std::exception&) {
            // try next candidate
        }
    }

    if (!saved) {
        std::cerr << "Failed to write log file to output directory.\n";
    }

    return 0;
}
