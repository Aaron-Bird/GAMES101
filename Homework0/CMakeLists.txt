cmake_minimum_required (VERSION 2.8.11)
project (Transformation)

find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

get_filename_component(EIGEN3_ROOT_DIR ${EIGEN3_INCLUDE_DIR} DIRECTORY)
include_directories(${EIGEN3_ROOT_DIR})

add_executable (Transformation main.cpp)