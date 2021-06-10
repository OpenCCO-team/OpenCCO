cmake_minimum_required(VERSION 3.5)
project(CCO_Algorithm)
set (CMAKE_CXX_STANDARD 11)


SET(SOURCE_HEADER_FILES
  ${PROJECT_SOURCE_DIR}/../base/CoronaryArteryTree.h
    ${PROJECT_SOURCE_DIR}/../helpers/geomhelpers.h
)



# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/../ext/" "${PROJECT_SOURCE_DIR}/../helpers/" "${PROJECT_SOURCE_DIR}/../base/" )
FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


# Schreiner_s_Algorithm
add_executable(testCCO testCCO.cpp ../base/CoronaryArteryTree.cpp ../base/CoronaryArteryTree.h ../helpers/geomhelpers.h)
target_link_libraries(testCCO  ${DGTAL_LIBRARIES})


