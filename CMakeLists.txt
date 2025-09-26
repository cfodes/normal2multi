cmake_minimum_required(VERSION 3.10)
project(MyProject LANGUAGES CXX)

# C++ 标准和构建类型
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)

# 让 CMake 导出 compile_commands.json，供 clangd 使用
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 查找 Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 查找 Boost
find_package(Boost REQUIRED)

# 查找 Metis
find_path(METIS_INCLUDE_DIR metis.h)
find_library(METIS_LIB metis)

if(NOT METIS_INCLUDE_DIR OR NOT METIS_LIB)
    message(FATAL_ERROR "METIS not found. Please install libmetis-dev")
endif()

# 源文件：把 src 下所有 cpp 文件都包含进来
file(GLOB SRC_FILES src/*.cpp)

# 构建可执行文件
add_executable(MyProject ${SRC_FILES})

# 添加头文件搜索路径
target_include_directories(MyProject PRIVATE
    ${Boost_INCLUDE_DIRS}
    ${METIS_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/include   # 你自己的头文件
)

# 链接库
target_link_libraries(MyProject PRIVATE Eigen3::Eigen ${METIS_LIB})
