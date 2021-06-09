cmake_minimum_required(VERSION 3.5)
project(Schreiner_s_Algorithm)
set (CMAKE_CXX_STANDARD 11)



# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/../ext/" "${PROJECT_SOURCE_DIR}/../helpers/" "${PROJECT_SOURCE_DIR}/src/" )
FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


# Schreiner_s_Algorithm
add_executable(testSchreiner testSchreiner.cpp src/CoronaryArteryTree)
target_link_libraries(testSchreiner  ${DGTAL_LIBRARIES})


