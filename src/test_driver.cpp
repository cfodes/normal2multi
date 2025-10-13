#include "test_driver.hpp"
#include "deformation.hpp"
#include "meshio.hpp"  // readfile / writefile
#include "multi_partition.hpp"
#include "preprocess_utils.hpp"
#include <ctime>

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
    read_mesh();  // 读取网格文件,s设置物面节点
    preprocess();  // 预处理，构建物面树和节点等级映射
    run_multi_partition();   // 运行多级分区RBF算法
    write_mesh();  // 写出变形后的网格文件
}

#include <iostream>
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
    last_timing_report_ = mp.last_timing_report();
    last_blockD_report_ = mp.block_d_report();

}

void TestDriver::write_mesh() 
{
    writefile(output_file_, S_);  // 输出 S_ 的数据
}
