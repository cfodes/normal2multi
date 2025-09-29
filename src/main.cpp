#include "test_driver.hpp"
#include <vector>
#include <string>

int main() 
{
    

    // 算例Naca0012 C型网格
    // 输入文件
    std::string input_file = "../data/naca0012_Cmesh.su2";
    // 输出文件
    std::string output_file = "../output/naca0012_Cmesh_result.su2";
    // 设置分区数和误差容限
    std::vector<idx_t> parts{11,3,2};
    std::vector<double> tol{1e-14,1e-14,1e-7};

    TestDriver driver(input_file,
                      output_file,
                      parts, tol);
    driver.run();
}
