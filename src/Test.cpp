#include "Test.hpp"
#include "meshio.hpp"
#include "geometry.hpp"
#include "test_driver.hpp"
#include "XlsxWriter.hpp"
#include <filesystem>
#include <iostream>
#include <ctime>
#include <chrono>
#include <stdexcept>
#include <vector>
using namespace std::chrono;

void RBFTest::run_global_test(const std::string& input_file,
                              const std::string& output_file,
                              double tol,
                              const State& S)
{
    std::cout << "\n========== Global RBF Test ==========\n";

    // 1. 读取网格（支持从 build/ 运行，优先解析 ../data/...）
    State d_S = S;
    namespace fs = std::filesystem;
    fs::path in_path(input_file);
    if (!fs::exists(in_path)) {
        const std::vector<fs::path> candidates = {
            fs::path("..") / in_path,
            fs::path("./") / in_path
        };
        for (const auto& p : candidates) {
            if (fs::exists(p)) { in_path = p; break; }
        }
    }
    readfile(in_path.string(), d_S);

    // 2. 提取物面节点
    int wall_id = 0;
    std::vector<Node> wall_nodes = Set_wall_nodes(d_S.every_boundary, d_S.node_coords, wall_id);
    std::cout << "Read mesh from: " << input_file
              << " | wall nodes: " << wall_nodes.size() << std::endl;

    // 3. 计算物面变形
    select_R(d_S.node_coords, d_S.R);
    calculat_wall_deformation(wall_nodes, d_S.D);
    std::cout << "Wall deformation set. D = " << d_S.D << std::endl;

    const fs::path output_path(output_file);
    fs::path output_parent = output_path.parent_path();
    if (output_parent.empty()) {
        output_parent = fs::path{"."};
    }
    const fs::path case_dir = output_parent / output_path.stem();
    std::error_code ec;
    fs::create_directories(case_dir, ec);
    if (ec) {
        throw std::runtime_error("Failed to create global test output directory '" +
                                 case_dir.string() + "': " + ec.message());
    }
    const fs::path case_output_path = case_dir / output_path.filename();

    // 4. 构建并运行 RBF（记录用时）
    RBFInterpolator rbf;
    const auto t_build_start = steady_clock::now();
    rbf.Greedy_algorithm(wall_nodes, tol, d_S);
    const auto t_build_end = steady_clock::now();
    const double build_ms = duration<double, std::milli>(t_build_end - t_build_start).count();
    std::cout << "RBF system built in " << build_ms << " ms\n";

    // 5. 计算所有节点的变形（记录用时）
    DeformCalculator deform_calc(rbf);
    const auto t_deform_start = steady_clock::now();
    deform_calc.calculate_deform(d_S.node_coords, d_S);
    const auto t_deform_end = steady_clock::now();
    const double deform_ms = duration<double, std::milli>(t_deform_end - t_deform_start).count();

    // 6. 更新坐标（记录用时）
    const auto t_update_start = steady_clock::now();
    calculate_deformed_coordinates(d_S.node_coords);
    const auto t_update_end = steady_clock::now();
    const double update_ms = duration<double, std::milli>(t_update_end - t_update_start).count();

    // 6. 输出文件
    writefile(case_output_path.string(), d_S);
    std::cout << "Deformed mesh written to: " << case_output_path.string() << std::endl;
    std::cout << "=====================================\n";

    // === 写出全局 timing 信息 ===
    std::vector<std::vector<XlsxCell>> timing_rows{
        {XlsxCell::String("phase"), XlsxCell::String("ms")},
        {XlsxCell::String("build_rbf"), XlsxCell::Number(build_ms)},
        {XlsxCell::String("compute_deform"), XlsxCell::Number(deform_ms)},
        {XlsxCell::String("update_coordinates"), XlsxCell::Number(update_ms)}
    };

    XlsxWriter timing_writer;
    timing_writer.add_sheet("timing", timing_rows);
    const fs::path timing_xlsx = case_dir / (output_path.stem().string() + "_global_timings.xlsx");
    timing_writer.save(timing_xlsx.string());

    // === 写出 test info ===
    const std::size_t candidate_points = wall_nodes.size();
    const std::size_t support_points = rbf.suppoints.size();
    const std::size_t deformed_points = d_S.node_coords.size();

    std::vector<std::vector<XlsxCell>> summary_rows{
        {XlsxCell::String("metric"), XlsxCell::String("value")},
        {XlsxCell::String("candidate_points"), XlsxCell::Number(static_cast<double>(candidate_points))},
        {XlsxCell::String("support_points"), XlsxCell::Number(static_cast<double>(support_points))},
        {XlsxCell::String("deformed_nodes"), XlsxCell::Number(static_cast<double>(deformed_points))}
    };

    std::vector<std::vector<XlsxCell>> greedy_rows;
    greedy_rows.push_back({XlsxCell::String("step"), XlsxCell::String("max_interp_error")});
    for (std::size_t i = 0; i < rbf.max_tol_step.size(); ++i) {
        greedy_rows.push_back({XlsxCell::Number(static_cast<double>(i + 1)),
                               XlsxCell::Number(rbf.max_tol_step[i])});
    }

    XlsxWriter info_writer;
    info_writer.add_sheet("summary", summary_rows);
    info_writer.add_sheet("greedy_steps", greedy_rows);
    const fs::path info_xlsx = case_dir / (output_path.stem().string() + "_global_test_info.xlsx");
    info_writer.save(info_xlsx.string());
}

MultiPartitionBatch::MultiPartitionBatch(std::string test_name,
                                         std::string input_file,
                                         std::string output_root)
    : test_name_(std::move(test_name)),
      input_file_(std::move(input_file)),
      output_root_(output_root.empty() ? "output" : std::move(output_root))
{}

void MultiPartitionBatch::add_case(const std::string& run_name,
                                   const std::vector<idx_t>& parts,
                                   const std::vector<double>& tolerances,
                                   bool use_greedy_nonfinal)
{
    if (run_name.empty()) {
        throw std::invalid_argument("run_name must not be empty");
    }
    if (parts.empty()) {
        throw std::invalid_argument("parts configuration must not be empty");
    }
    if (tolerances.size() != parts.size()) {
        throw std::invalid_argument("parts and tolerances must have the same length");
    }
    cases_.push_back(PartitionBatchCase{run_name, parts, tolerances, use_greedy_nonfinal});
}

void MultiPartitionBatch::clear_cases()
{
    cases_.clear();
}

void MultiPartitionBatch::set_output_root(const std::string& output_root)
{
    output_root_ = output_root.empty() ? "output" : output_root;
}

void MultiPartitionBatch::run_all() const
{
    run_all_impl(false);
}

void MultiPartitionBatch::run_all_with_test_info() const
{
    run_all_impl(true);
}

void MultiPartitionBatch::run_all_impl(bool generate_test_info) const
{
    if (cases_.empty()) {
        std::cout << "[MultiPartitionBatch] No cases to run for test '" << test_name_ << "'.\n";
        return;
    }

    namespace fs = std::filesystem;
    const std::vector<fs::path> candidates = {
        fs::path("..") / output_root_ / test_name_,
        fs::path(output_root_) / test_name_
    };

    fs::path base_dir;
    std::error_code ec;
    for (const auto& candidate : candidates) {
        if (candidate.empty()) {
            continue;
        }
        ec.clear();
        fs::create_directories(candidate, ec);
        if (!ec || fs::exists(candidate)) {
            base_dir = candidate;
            break;
        }
    }

    if (base_dir.empty()) {
        throw std::runtime_error("Failed to create output directory for test '" + test_name_ + "'.");
    }

    std::vector<std::pair<std::string, std::vector<LevelTiming>>> reports;
    reports.reserve(cases_.size());

    struct BlockReport {
        std::string run_name;
        std::vector<idx_t> parts;
        std::vector<std::vector<BlockTestInfo>> data;
    };
    std::vector<BlockReport> block_reports;
    if (generate_test_info) {
        block_reports.reserve(cases_.size());
    }

    fs::path resolved_input(input_file_);
    if (!fs::exists(resolved_input)) {
        const std::vector<fs::path> candidates = {
            fs::path("..") / input_file_,
            fs::path("../") / input_file_,
            fs::path("./") / input_file_
        };
        bool found = false;
        for (const auto& candidate : candidates) {
            if (fs::exists(candidate)) {
                resolved_input = candidate;
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("Input file not found: " + input_file_);
        }
    }

    for (const auto& c : cases_) {
        fs::path file_name = c.run_name;
        if (file_name.extension().empty()) {
            file_name += ".su2";
        }
        const std::string file_stem = file_name.stem().string();
        const fs::path case_dir = base_dir / file_stem;
        ec.clear();
        fs::create_directories(case_dir, ec);
        if (ec) {
            throw std::runtime_error("Failed to create case output directory '" +
                                     case_dir.string() + "': " + ec.message());
        }
        const fs::path output_path = case_dir / file_name.filename();

        std::cout << "[MultiPartitionBatch] >>> Start case '" << c.run_name
                  << "' -> " << output_path.string() << std::endl;

        TestDriver driver(resolved_input.string(), output_path.string(), c.parts, c.tolerances);
        driver.set_collect_test_info(generate_test_info);
        driver.set_use_greedy(c.use_greedy_nonfinal);
        driver.run();
        std::cout << "[MultiPartitionBatch] <<< Finished case '" << c.run_name << "'"
                  << std::endl;
        {
            XlsxWriter case_writer;
            std::vector<std::vector<XlsxCell>> case_rows;
            case_rows.push_back({XlsxCell::String("lvl"),
                                 XlsxCell::String("preprocess_ms"),
                                 XlsxCell::String("build_ms"),
                                 XlsxCell::String("compute_ms"),
                                 XlsxCell::String("update_ms"),
                                 XlsxCell::String("distance_ms"),
                                 XlsxCell::String("apply_ms")});

            const auto& timings = driver.get_last_timing_report();
            for (size_t lvl = 0; lvl < timings.size(); ++lvl) {
                case_rows.push_back({XlsxCell::Number(static_cast<double>(lvl)),
                                     XlsxCell::Number(timings[lvl].preprocess_ms),
                                     XlsxCell::Number(timings[lvl].build_ms),
                                     XlsxCell::Number(timings[lvl].compute_ms),
                                     XlsxCell::Number(timings[lvl].update_ms),
                                     XlsxCell::Number(timings[lvl].distance_ms),
                                     XlsxCell::Number(timings[lvl].apply_ms)});
            }

            case_writer.add_sheet("timing", case_rows);
            const fs::path case_timing_path = case_dir / (file_stem + "_timings.xlsx");
            case_writer.save(case_timing_path.string());
            std::cout << "[MultiPartitionBatch] Case timing workbook written to: "
                      << case_timing_path.string() << std::endl;
        }

        if (generate_test_info) {
            XlsxWriter case_block_writer;
            std::vector<std::vector<XlsxCell>> rows;
            rows.push_back({XlsxCell::String("lvl"),
                            XlsxCell::String("part_count"),
                            XlsxCell::String("block_id"),
                            XlsxCell::String("internal_points"),
                            XlsxCell::String("candidate_points"),
                            XlsxCell::String("support_points"),
                            XlsxCell::String("block_D")});

            const auto& block_info = driver.get_last_block_info_report();
            for (size_t lvl = 0; lvl < block_info.size(); ++lvl) {
                const auto& blocks = block_info[lvl];
                const double part_count = (lvl < c.parts.size())
                    ? static_cast<double>(c.parts[lvl])
                    : 0.0;
                if (blocks.empty()) {
                    rows.push_back({XlsxCell::Number(static_cast<double>(lvl)),
                                    XlsxCell::Number(part_count),
                                    XlsxCell::String("N/A"),
                                    XlsxCell::String("N/A"),
                                    XlsxCell::String("N/A"),
                                    XlsxCell::String("N/A"),
                                    XlsxCell::String("N/A")});
                    continue;
                }
                for (const auto& entry : blocks) {
                    rows.push_back({XlsxCell::Number(static_cast<double>(lvl)),
                                    XlsxCell::Number(part_count),
                                    XlsxCell::Number(static_cast<double>(entry.block_id)),
                                    XlsxCell::Number(static_cast<double>(entry.internal_points)),
                                    XlsxCell::Number(static_cast<double>(entry.candidate_points)),
                                    XlsxCell::Number(static_cast<double>(entry.support_points)),
                                    XlsxCell::Number(entry.block_D)});
                }
            }

            case_block_writer.add_sheet("blocks", rows);
            const fs::path case_info_path = case_dir / (file_stem + "_test_info.xlsx");
            case_block_writer.save(case_info_path.string());
            std::cout << "[MultiPartitionBatch] Case test-info workbook written to: "
                      << case_info_path.string() << std::endl;
        }

        reports.emplace_back(c.run_name, driver.get_last_timing_report());
        if (generate_test_info) {
            block_reports.push_back(BlockReport{c.run_name, c.parts, driver.get_last_block_info_report()});
        }
    }

    if (!reports.empty()) {
        XlsxWriter writer;
        for (const auto& report : reports) {
            std::vector<std::vector<XlsxCell>> rows;
            rows.push_back({XlsxCell::String("lvl"),
                            XlsxCell::String("preprocess_ms"),
                            XlsxCell::String("build_ms"),
                            XlsxCell::String("compute_ms"),
                            XlsxCell::String("update_ms"),
                            XlsxCell::String("distance_ms"),
                            XlsxCell::String("apply_ms")});

            const auto& timings = report.second;
            for (size_t lvl = 0; lvl < timings.size(); ++lvl) {
                rows.push_back({XlsxCell::Number(static_cast<double>(lvl)),
                                XlsxCell::Number(timings[lvl].preprocess_ms),
                                XlsxCell::Number(timings[lvl].build_ms),
                                XlsxCell::Number(timings[lvl].compute_ms),
                                XlsxCell::Number(timings[lvl].update_ms),
                                XlsxCell::Number(timings[lvl].distance_ms),
                                XlsxCell::Number(timings[lvl].apply_ms)});
            }

            writer.add_sheet(report.first, rows);
        }

        const fs::path workbook_path = base_dir / (test_name_ + "_timings.xlsx");
        writer.save(workbook_path.string());
        std::cout << "[MultiPartitionBatch] Timing workbook written to: "
                  << workbook_path.string() << std::endl;
    }

    if (generate_test_info && !block_reports.empty()) {
        XlsxWriter block_writer;
        for (const auto& report : block_reports) {
            std::vector<std::vector<XlsxCell>> rows;
            rows.push_back({XlsxCell::String("lvl"),
                            XlsxCell::String("part_count"),
                            XlsxCell::String("block_id"),
                            XlsxCell::String("internal_points"),
                            XlsxCell::String("candidate_points"),
                            XlsxCell::String("support_points"),
                            XlsxCell::String("block_D")});

            for (size_t lvl = 0; lvl < report.data.size(); ++lvl) {
                const auto& blocks = report.data[lvl];
                const double part_count = (lvl < report.parts.size())
                    ? static_cast<double>(report.parts[lvl])
                    : 0.0;
                if (blocks.empty()) {
                    rows.push_back({XlsxCell::Number(static_cast<double>(lvl)),
                                    XlsxCell::Number(part_count),
                                    XlsxCell::String("N/A"),
                                    XlsxCell::String("N/A"),
                                    XlsxCell::String("N/A"),
                                    XlsxCell::String("N/A"),
                                    XlsxCell::String("N/A")});
                    continue;
                }
                for (const auto& entry : blocks) {
                    rows.push_back({XlsxCell::Number(static_cast<double>(lvl)),
                                    XlsxCell::Number(part_count),
                                    XlsxCell::Number(static_cast<double>(entry.block_id)),
                                    XlsxCell::Number(static_cast<double>(entry.internal_points)),
                                    XlsxCell::Number(static_cast<double>(entry.candidate_points)),
                                    XlsxCell::Number(static_cast<double>(entry.support_points)),
                                    XlsxCell::Number(entry.block_D)});
                }
            }

            block_writer.add_sheet(report.run_name, rows);
        }

        const fs::path workbook_path = base_dir / (test_name_ + "_test_info.xlsx");
        block_writer.save(workbook_path.string());
        std::cout << "[MultiPartitionBatch] test info workbook written to: "
                  << workbook_path.string() << std::endl;
    }
}
