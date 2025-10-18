#include "Test.hpp"
#include "meshio.hpp"
#include "geometry.hpp"
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <system_error>
using namespace std::chrono;

void RBFTest::run_global_test(const std::string& input_file,
                              const std::string& output_file,
                              double tol,
                              const State& S)
{
    std::cout << "\n========== Global RBF Test ==========\n";

    // 1. 读取网格
    State d_S = S;
    readfile(input_file, d_S);

    // 2. 提取物面节点
    int wall_id = 0;
    std::vector<Node> wall_nodes = Set_wall_nodes(d_S.every_boundary, d_S.node_coords, wall_id);
    std::cout << "Read mesh from: " << input_file
              << " | wall nodes: " << wall_nodes.size() << std::endl;

    // 3. 计算物面变形
    select_R(d_S.node_coords, d_S.R);
    calculat_wall_deformation(wall_nodes, d_S.D);
    std::cout << "Wall deformation set. D = " << d_S.D << std::endl;

    // 4. 构建并运行 RBF
    RBFInterpolator rbf;
    const auto t0 = steady_clock::now();
    rbf.Greedy_algorithm(wall_nodes, tol, d_S);
    const auto t1 = steady_clock::now();
    std::cout << "RBF system built in "
        << duration<double, std::milli>(t1 - t0).count() << " ms\n";

    // 5. 计算所有节点的变形
    DeformCalculator deform_calc(rbf);
    deform_calc.calculate_deform(d_S.node_coords, d_S);
    calculate_deformed_coordinates(d_S.node_coords);

    // 6. 输出文件
    writefile(output_file, d_S);
    std::cout << "Deformed mesh written to: " << output_file << std::endl;
    std::cout << "=====================================\n";

    namespace fs = std::filesystem;
    fs::path out_path(output_file);
    fs::path out_dir = out_path.parent_path();
    if (out_dir.empty()) {
        out_dir = fs::path{"."};
    }

    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "Failed to create output directory " << out_dir << ": " << ec.message() << '\n';
    }

    const std::string base_name = out_path.stem().string();
    const fs::path info_file = out_dir / (base_name + "_info.dat");

    std::ofstream ofs(info_file);
    if (!ofs.is_open()) {
        std::cerr << "Failed to write global info to " << info_file << '\n';
        return;
    }

    ofs << std::fixed << std::setprecision(6);
    ofs << "wall_nodes " << wall_nodes.size() << '\n';
    ofs << "support_points " << rbf.suppoints.size() << '\n';
}
