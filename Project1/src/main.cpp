#include "test_driver.hpp"
#include <vector>
#include <string>

int main() 
{
    

    //// 算例Naca0012 C型网格
    //// 输入文件
    //std::string input_file = "D:/oujc/works/mesh_normal2multi/Project1/data/naca0012_Cmesh_test.su2";
    //// 输出文件
    //std::string output_file = "D:/oujc/works/mesh_normal2multi/Project1/output/naca0012_Cmesh_result.su2";
    //// 设置分区数和误差容限
    //std::vector<idx_t> parts{11,3,2};
    //std::vector<double> tol{1e-14,1e-14,1e-7};

    // 算例ONERA M6结构网格
    // 输入文件
    std::string input_file = "D:/oujc/works/mesh_normal2multi/Project1/data/onera_m6_test.su2";
    // 输出文件
    std::string output_file = "D:/oujc/works/mesh_normal2multi/Project1/output/onera_m6_result.su2";
    // 设置分区数和误差容限
    std::vector<idx_t> parts{ 16,4,2 };
    std::vector<double> tol{ 1e-14,1e-14,1e-7 };

    TestDriver driver(input_file,
                      output_file,
                      parts, tol);
    driver.run();
}
