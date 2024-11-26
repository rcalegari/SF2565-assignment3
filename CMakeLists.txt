cmake_minimum_required(VERSION 3.29)
project(ass3)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Eigen3
find_package(Eigen3 REQUIRED PATHS /usr/local/include/eigen3)

add_executable(ass3 main.cpp
        utilities/point.h
        utilities/curve.h
        utilities/grid.h)

# Link Eigen3
target_link_libraries(ass3 PRIVATE Eigen3::Eigen)
