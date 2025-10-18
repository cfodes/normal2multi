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
            "../data/onera_m6_294K.su2",
            "../output/onera_m6_294K_result.su2",
            {30, 4, 2},
            {1e-14, 1e-14, 1e-7}
        }
        // 可以在此处继续添加更多分组算例
    };

    std::vector<GlobalCaseConfig> global_cases = {
        {
            "../data/onera_m6_294K.su2",
            "../output/onera_m6_294K_global_result.su2",
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
