#include "test_driver.hpp"
#include "Test.hpp"
#include <iostream>
#include <metis.h>
#include <string>
#include <vector>

struct PartitionCaseConfig {
    std::string input_file;
    std::string output_file;
    std::vector<idx_t> parts;
    std::vector<double> tol;
};

struct GlobalCaseConfig {
    std::string input_file;
    std::string output_file;
    double tol;
};

int main()
{
    std::vector<PartitionCaseConfig> partition_cases = {
        {
            "../data/onera_m6_2M4.su2",
            "../output/onera_m6_2M4_80_16_g_eps6.su2",
            {80, 16, 2},
            {1e-6, 1e-6, 1e-5}
        },
        {
            "../data/onera_m6_2M4.su2",
            "../output/onera_m6_2M4_40_32_g_eps6.su2",
            {40, 32, 2},
            {1e-6, 1e-6, 1e-5}
        },
        {
            "../data/onera_m6_2M4.su2",
            "../output/onera_m6_2M4_40_16_g_eps6.su2",
            {40, 16, 2},
            {1e-6, 1e-6, 1e-5}
        },
        {
            "../data/onera_m6_2M4.su2",
            "../output/onera_m6_2M4_20_32_g_eps6.su2",
            {20, 32, 2},
            {1e-6, 1e-6, 1e-5}
        },
        {
            "../data/onera_m6_2M4.su2",
            "../output/onera_m6_2M4_20_16_g_eps6.su2",
            {20, 16, 2},
            {1e-6, 1e-6, 1e-5}
        }
        // 可以在此处继续添加更多分组算例
    };

    std::vector<GlobalCaseConfig> global_cases = {
        {
            "../data/onera_m6_2M4.su2",
            "../output/onera_m6_2M4_global_result.su2",
            1e-5
        }
        // 可以在此处继续添加更多全局RBF算例
    };

    for (const auto& cfg : partition_cases) {
        std::cout << "Running partitioned case with input: " << cfg.input_file
                  << " -> output: " << cfg.output_file << '\n';
        TestDriver driver(cfg.input_file, cfg.output_file, cfg.parts, cfg.tol);
        driver.run();
    }

    State default_state;
    for (const auto& cfg : global_cases) {
        std::cout << "Running global RBF case with input: " << cfg.input_file
                  << " -> output: " << cfg.output_file << '\n';
        RBFTest::run_global_test(cfg.input_file, cfg.output_file, cfg.tol, default_state);
    }

    return 0;
}
