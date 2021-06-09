cmake_minimum_required(VERSION 3.5)

project(Schreiner_s_Algorithm)

find_package(Ceres REQUIRED)

# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/../ext/" )
FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


# Schreiner_s_Algorithm
add_executable(testSchreiner testSchreiner.cpp)
target_link_libraries(testSchreiner ${DGTAL_LIBRARIES})


