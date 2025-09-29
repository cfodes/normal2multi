cmake_minimum_required(VERSION 3.10)
project(MyProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)

# Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Metis (Linux)
find_path(METIS_INCLUDE_DIR metis.h)
find_library(METIS_LIB metis)

if (NOT METIS_INCLUDE_DIR OR NOT METIS_LIB)
    message(FATAL_ERROR "METIS not found. Please install libmetis-dev")
endif()

# 源文件
file(GLOB SRC_FILES src/*.cpp)

add_executable(MyProject ${SRC_FILES})

target_include_directories(MyProject PRIVATE
    ${METIS_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(MyProject PRIVATE Eigen3::Eigen ${METIS_LIB})
