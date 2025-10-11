#include "test_driver.hpp"
#include "Test.hpp"
#include <vector>
#include <string>

int main() 
{
    

    //// 算例Naca0012 C型网格 分组RBF算法
    //// 输入文件
    //std::string input_file = "D:/oujc/works/mesh_normal2multi/Project1/data/naca0012_Cmesh_test.su2";
    //// 输出文件
    //std::string output_file = "D:/oujc/works/mesh_normal2multi/Project1/output/naca0012_Cmesh_result.su2";
    //// 设置分区数和误差容限
    //std::vector<idx_t> parts{11,3,2};
    //std::vector<double> tol{1e-14,1e-14,1e-7};

    //// 算例ONERA M6混合网格 分组RBF算法
    //// 输入文件
    //std::string input_file = "../data/onera_m6_2M4.su2";
    //// 输出文件
    //std::string output_file = "../output/onera_m6_2M4_result.su2";
    //// 设置分区数和误差容限
    //std::vector<idx_t> parts{ 199,24,2 };
    //std::vector<double> tol{ 1e-14,1e-14,1e-7 };

    // 算例ONERA M6结构网格 分组RBF算法
    // 输入文件
    std::string input_file = "../data/onera_m6_294K.su2";
    // 输出文件
    std::string output_file = "../output/onera_m6_294K_result.su2";
    // 设置分区数和误差容限
    std::vector<idx_t> parts{ 30,4,2 };
    std::vector<double> tol{ 1e-14,1e-14,1e-7 };

    TestDriver driver(input_file,
                      output_file,
                      parts, tol);
    driver.run();

    //// 算例ONERA M6纯RBF方法
    //State _S;
    //RBFTest::run_global_test("../data/onera_m6_294K.su2",
    //    "../output/onera_m6_294K_result.su2",
    //    1e-5, _S);

    return 0;
}
