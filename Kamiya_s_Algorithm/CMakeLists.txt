cmake_minimum_required(VERSION 3.5)

project(Kamiya_s_Algorithm)

find_package(Ceres REQUIRED)

# Kamiya_s_Algorithm
add_executable(testKamiya testKamiya.cpp)
target_link_libraries(testKamiya Ceres::ceres)


