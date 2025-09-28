cmake_minimum_required(VERSION 3.10)
project(Rasterizer)

find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
get_filename_component(EIGEN3_ROOT_DIR ${EIGEN3_INCLUDE_DIR} DIRECTORY)
include_directories(${EIGEN3_ROOT_DIR})

link_directories(${Eigen_INCLUDE_DIRS})
set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/local/include)

add_executable(Rasterizer main.cpp rasterizer.hpp rasterizer.cpp Triangle.hpp Triangle.cpp)
target_link_libraries(Rasterizer ${OpenCV_LIBRARIES})
