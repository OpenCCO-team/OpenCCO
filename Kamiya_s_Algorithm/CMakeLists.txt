cmake_minimum_required(VERSION 3.5)

project(Kamiya_s_Algorithm)

find_package(Ceres REQUIRED)

# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/../ext/" )
FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


# Kamiya_s_Algorithm
add_executable(testKamiya testKamiya.cpp)
target_link_libraries(testKamiya Ceres::ceres ${DGTAL_LIBRARIES})


