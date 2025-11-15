#include "test_driver.hpp"
#include "deformation.hpp"
#include "meshio.hpp"  // readfile / writefile
#include "multi_partition.hpp"
#include "preprocess_utils.hpp"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <system_error>

TestDriver::TestDriver(const std::string& input_file,
                       const std::string& output_file,
                       const std::vector<idx_t>& parts_per_level,
                       const std::vector<double>& tol_steps)
    : input_file_(input_file),
      output_file_(output_file),
      parts_per_level_(parts_per_level),
      tol_steps_(tol_steps)
{}

void TestDriver::run() 
// 主执行流程
{
    level_timings_.clear();
    level_statistics_.clear();
    block_summaries_.clear();
    nd2wall_lvl_.clear();
    S_ = State{};

    read_mesh();  // 读取网格文件,s设置物面节点
    preprocess();  // 预处理，构建物面树和节点等级映射
    run_multi_partition();   // 运行多级分区RBF算法
    write_mesh();  // 写出变形后的网格文件
    write_reports(); // 输出统计信息
}

void TestDriver::read_mesh() 
{
    readfile(input_file_, S_);  // 填充 S_ 的数据
    // 设置物面节点
    S_.wall_nodes = Set_wall_nodes(S_.every_boundary, S_.node_coords, S_.wall_id);
    std::cout << "the number of wall nodes is "<<S_.wall_nodes.size() << std::endl;
}

void TestDriver::preprocess() 
{
    // 选择 R
    select_R(S_.node_coords, S_.R);
    // 计算 D
    calculat_wall_deformation(S_.wall_nodes, S_.D);
    initial_D_ = S_.D;

    // 构建 wall_tree + 节点等级映射
    GridBTree<int, Point<double>> wall_tree_tmp;
    _set_wall_tree(wall_tree_tmp, S_.wall_nodes);
    _set_nd2wall_lvl_mp(nd2wall_lvl_, S_.node_coords, wall_tree_tmp, S_.D);
}

void TestDriver::run_multi_partition()
{
    multi_partition mp(parts_per_level_);
    mp.set_m_nd2wall_lvl_mp(nd2wall_lvl_);


    mp.divide_wall(S_.every_boundary[S_.wall_id].bound_elements, S_.wall_nodes);
    mp.multi_partition_rbf_algorithm(tol_steps_, S_, S_.wall_nodes, S_.node_coords);

    level_timings_ = mp.get_level_timings();
    level_statistics_ = mp.get_level_statistics();

    //将数据传给block_summaries
    block_summaries_.assign(level_timings_.size(), {});
    for (size_t lvl = 0; lvl < block_summaries_.size(); ++lvl) {
        const auto stats = mp.get_block_stats(lvl,S_.wall_nodes);
        auto& summaries = block_summaries_[lvl];
        summaries.reserve(stats.size());
        for (const auto& st : stats) {
            BlockSummary bs;
            bs.block_id = st.block_id;
            bs.candidate_points = st.candidate_points;
            bs.support_points = st.support_points;
            bs.block_D = st.block_D;
            bs.internal_nodes = st.internal_points;
            bs.boundary_nodes = st.boundary_points;
            summaries.push_back(std::move(bs));
        }
    }

}

void TestDriver::write_mesh() 
{
    writefile(output_file_, S_);  // 输出 S_ 的数据
}

void TestDriver::write_reports()
{
    namespace fs = std::filesystem;
    fs::path out_path(output_file_);
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

    write_time_report(time_file);
    write_info_report(info_file);
}

void TestDriver::write_time_report(const std::filesystem::path& time_path) const
{
    std::ofstream ofs(time_path);
    if (!ofs.is_open()) {
        std::cerr << "Failed to write timing report to " << time_path << '\n';
        return;
    }

    ofs << std::fixed << std::setprecision(6);
    ofs << "# level preprocess_ms build_rbf_ms distance_search_ms drrbf_deform_ms compute_df_ms update_coords_ms\n";
    for (size_t lvl = 0; lvl < level_timings_.size(); ++lvl) {
        const auto& t = level_timings_[lvl];
        ofs << lvl << ' ' << t.preprocess_ms << ' ' << t.build_rbf_ms << ' '
            << t.distance_search_ms << ' ' << t.drrbf_deform_ms << ' '
            << t.compute_df_ms << ' ' << t.update_coords_ms << '\n';
    }
}

void TestDriver::write_info_report(const std::filesystem::path& info_path) const
{
    std::ofstream ofs(info_path);
    if (!ofs.is_open()) {
        std::cerr << "Failed to write info report to " << info_path << '\n';
        return;
    }

    ofs << std::fixed << std::setprecision(6);

    for (size_t lvl = 0; lvl < level_statistics_.size(); ++lvl) {
        ofs << "level " << lvl << '\n';
        ofs << "blocks\n";
        ofs << "# block_id candidate_p support_p internal_nodes boundary_nodes block_D\n"; // 新增列头

        if (lvl < block_summaries_.size()) {
            const auto& summaries = block_summaries_[lvl];
            for (const auto& blk : summaries) {
                double block_D_value = blk.block_D;
                if (lvl == 0) {
                    block_D_value = initial_D_;
                }

                ofs << "  "
                    << blk.block_id << ' '
                    << blk.candidate_points << ' '
                    << blk.support_points << ' '
                    << blk.internal_nodes << ' '   // 输出 internal_nodes
                    << blk.boundary_nodes << ' '   // 输出 boundary_nodes
                    << block_D_value << '\n';
            }
        }
        ofs << '\n';
    }
}
