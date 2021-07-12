cmake_minimum_required(VERSION 3.5)
project(Schreiner_s_Algorithm)
set (CMAKE_CXX_STANDARD 11)


SET(SOURCE_HEADER_FILES
  ${PROJECT_SOURCE_DIR}/../base/CoronaryArteryTree.h
    ${PROJECT_SOURCE_DIR}/../helpers/geomhelpers.h
)

find_package(Ceres REQUIRED)


# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/../ext/" "${PROJECT_SOURCE_DIR}/../helpers/" "${PROJECT_SOURCE_DIR}/../base/" )
FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


# Schreiner_s_Algorithm
add_executable(testSchreiner testSchreiner.cpp ../base/CoronaryArteryTree.cpp ${SOURCE_HEADER_FILES} )
target_link_libraries(testSchreiner Ceres::ceres  ${DGTAL_LIBRARIES})


