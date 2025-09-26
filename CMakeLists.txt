cmake_minimum_required(VERSION 3.10)
project(MyProject LANGUAGES CXX)

# C++ 标准和构建类型
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)

# Homebrew 安装路径
set(HOMEBREW_PREFIX "/opt/homebrew/opt")

# 查找 Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Boost 和 Metis 路径
set(BOOST_INCLUDE_DIR "${HOMEBREW_PREFIX}/boost/include")
set(METIS_INCLUDE_DIR "${HOMEBREW_PREFIX}/metis/include")
set(METIS_LIB        "${HOMEBREW_PREFIX}/metis/lib/libmetis.dylib")

# 源文件：把 src 下所有 cpp 文件都包含进来
file(GLOB SRC_FILES src/*.cpp)

# 构建可执行文件
add_executable(MyProject ${SRC_FILES})

# 添加头文件搜索路径
target_include_directories(MyProject PRIVATE
    ${BOOST_INCLUDE_DIR}
    ${METIS_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/include   # 你自己的头文件
    ${PROJECT_SOURCE_DIR}/src       # 如果以后有内部头文件
)

# 链接库
target_link_libraries(MyProject PRIVATE Eigen3::Eigen ${METIS_LIB})
