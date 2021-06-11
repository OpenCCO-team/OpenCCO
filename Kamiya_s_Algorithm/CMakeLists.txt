cmake_minimum_required(VERSION 3.5)

project(Kamiya_s_Algorithm)

find_package(Ceres REQUIRED)


SET(SOURCE_HEADER_FILES
   ${PROJECT_SOURCE_DIR}/../helpers/geomhelpers.h
)


# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/../ext/" "${PROJECT_SOURCE_DIR}/../helpers/" )
FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


# Kamiya_s_Algorithm
add_executable(testKamiya testKamiya.cpp ${SOURCE_HEADER_FILES})
target_link_libraries(testKamiya Ceres::ceres ${DGTAL_LIBRARIES})


