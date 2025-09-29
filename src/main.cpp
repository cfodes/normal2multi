#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <unordered_set>  //std::unorder_set
#include <unordered_map>     //std::unordered_map
#include<string>
#include<sstream>
#include <list>
#include <stdio.h>
#include <iomanip>  //fixed and setprecision
#include<cmath>
#include<cassert>  //assert

#include <metis.h>
//#include "ADT.hpp"
#include "./distance.hpp"
#include "geometry.hpp"

#include <Eigen/Core>          
#include <Eigen/SVD>
#include <Eigen/LU>

#include <filesystem>
namespace fs = std::filesystem;
using namespace std;



int main() {
    
    ////O型结果
    //string readfile_name = "naca0012.su2";
    //string writefile_name = "su2_deformed_naca0012.su2";

    //C型结果
    //string readfile_name = "naca0012_Cmesh_test.su2";
    //string writefile_name = "naca0012_Cmesh_test_result.su2";

    string readfile_name = "../data/naca0012_Cmesh.su2";
    string writefile_name = "../output/naca0012_Cmesh_result.su2";

    //double greedy_tol_0 = 1e-14;  //设置贪心算法误差
    //double greedy_tol_1 = 1e-8;  //设置贪心算法误差
    //idx_t metis_k = 16;  //metis分区数
    //std::vector<deform_calculate> block_RBF;      //分区的RBF系统数组

    readfile(readfile_name);  //读取网格文件
    _wall_nodes = Set_wall_nodes(_every_boundary, _node_coords);   //设置物面节点数组
    _select_R(_node_coords, _R);           //选出网格点间近似的最大距离
    calculat_wall_deformation(_wall_nodes,_D);                        //计算物面节点位移
   
    _set_wall_tree(_wall_tree, _wall_nodes);      //构建物面节点的二叉树
    _set_nd2wall_lvl_mp(_nd2wall_lvl_mp, _node_coords, _wall_tree, _D);   //构建空间节点id和对应等级的映射



    ////***************全局RBF算法开始**********************
    //auto t_global_start = clock();

    //auto t_global_rbf_start = clock();
    //_RBF.Greedy_algorithm(_wall_nodes, 1e-7);   //和分组RBF相同的误差设置
    //auto t_global_rbf_end = clock();

    //auto t_global_df_start = clock();
    //_RBF.calculate_every_node_deform(_node_coords);
    //auto t_global_df_end = clock();

    //    std::cout << "_RBF coeff calculating time: " << fixed << setprecision(2) << double(t_global_rbf_end - t_global_rbf_start) << " ms" << std::endl;
    //std::cout << "_RBF df calculating time: " << fixed << setprecision(2) << double(t_global_df_end - t_global_df_start) << " ms" << std::endl;

    //auto t_global_deform_start = clock();
    //_calculate_deformed_coordinates(_node_coords);
    //auto t_global_deform_end = clock();

    //auto t_global_end = clock();

    //std::cout << "_RBF deforming nodes time: " << fixed << setprecision(2) << double(t_global_deform_end - t_global_deform_start) << " ms" << std::endl;
    //std::cout << "_RBF total time: " << fixed << setprecision(2) << double(t_global_end - t_global_start) << " ms" << std::endl;
    ////***************全局RBF算法结束**********************
     
     
     
   
    
    
    //测试class multi_partition
    std::vector<idx_t> v_temp{9,2,2};
    std::vector<double> steps_tol{ 1e-14,1e-14,1e-7};
    multi_partition mp(v_temp);
    //mp.set_m_nd2wall_lvl_mp(_nd2wall_lvl_mp);
    clock_t blk_divide_start, blk_divide_end;
    clock_t blk_rbf_start, blk_rbf_end;

    blk_divide_start = clock();
    mp.divide_wall(_every_boundary[0].bound_elements, _wall_nodes);
    blk_divide_end = clock();

    blk_rbf_start = clock();
    mp.multi_partition_rbf_algorithm(steps_tol, _wall_nodes, _node_coords);
    blk_rbf_end = clock();

    std::cout << "divide wall use time: " << fixed << setprecision(2) << double(blk_divide_end - blk_divide_start) << std::endl;
    std::cout << "block rbf coeff calculating and deforming mesh use time: " << fixed << setprecision(2) << double(blk_rbf_end - blk_rbf_start) << std::endl;
    
    ////测试壁面数据
    //_wall_nodes.resize(0);
    //_wall_nodes = Set_wall_nodes(_every_boundary, _node_coords);
    //
    //_write_wall(_wall_nodes);
    //
    //
    writefile(writefile_name);

   


    return 0;
}