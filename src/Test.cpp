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

    // 4. 构建并运行 RBF（计时：build）
    RBFInterpolator rbf;
    const auto t_build_s = steady_clock::now();
    rbf.Greedy_algorithm(wall_nodes, tol, d_S);
    const auto t_build_e = steady_clock::now();
    const double build_ms = duration<double, std::milli>(t_build_e - t_build_s).count();
    std::cout << "RBF system built in " << build_ms << " ms\n";

    // 5. 计算所有节点的变形（计时：compute）
    DeformCalculator deform_calc(rbf);
    const auto t_deform_s = steady_clock::now();
    deform_calc.calculate_deform(d_S.node_coords, d_S);
    const auto t_deform_e = steady_clock::now();
    const double deform_ms = duration<double, std::milli>(t_deform_e - t_deform_s).count();

    // 6. 更新坐标（计时：update）
    const auto t_update_s = steady_clock::now();
    calculate_deformed_coordinates(d_S.node_coords);
    const auto t_update_e = steady_clock::now();
    const double update_ms = duration<double, std::milli>(t_update_e - t_update_s).count();

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
    const fs::path time_file = out_dir / (base_name + "_time.dat");
    const fs::path info_file = out_dir / (base_name + "_info.dat");

    // 写 timing.dat，格式参考 partition_cases 的 *_time.dat
    {
        std::ofstream t_ofs(time_file);
        if (!t_ofs.is_open()) {
            std::cerr << "Failed to write global timing to " << time_file << '\n';
        } else {
            t_ofs << std::fixed << std::setprecision(6);
            t_ofs << "# level preprocess_ms build_rbf_ms distance_search_ms drrbf_deform_ms compute_df_ms update_coords_ms\n";
            // 全局测试仅有一行（视为 level 0），无 preprocess/距离搜索/DRRBF 分离，置 0
            t_ofs << 0 << ' ' << 0.0 << ' ' << build_ms << ' ' << 0.0 << ' ' << 0.0
                 << ' ' << deform_ms << ' ' << update_ms << '\n';
        }
    }

    std::ofstream ofs(info_file);
    if (!ofs.is_open()) {
        std::cerr << "Failed to write global info to " << info_file << '\n';
        return;
    }

    ofs << std::fixed << std::setprecision(6);
    ofs << "wall_nodes " << wall_nodes.size() << '\n';
    ofs << "support_points " << rbf.suppoints.size() << '\n';
}
