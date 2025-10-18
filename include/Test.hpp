#pragma once
#include "State.hpp"
#include "rbf.hpp"
#include "deformation.hpp"
#include <string>

class RBFTest
{
public:
    // 完整运行流程：读入 → 变形 → 输出
    static void run_global_test(const std::string& input_file,
                                const std::string& output_file,
                                double tol,
                                const State& S);
};
