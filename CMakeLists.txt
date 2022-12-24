cmake_minimum_required(VERSION 3.5)
project(CCO_Implementation)
set (CMAKE_CXX_STANDARD 11)


SET(SOURCE_BASE_HEADER_FILES
  ${PROJECT_SOURCE_DIR}/helpers/GeomHelpers.h
  ${PROJECT_SOURCE_DIR}/helpers/ConstructionHelpers.h
  ${PROJECT_SOURCE_DIR}/base/CoronaryArteryTree.h
  ${PROJECT_SOURCE_DIR}/helpers/XmlHelpers.h
)
SET(SOURCE_BASE_FILES
  ${PROJECT_SOURCE_DIR}/base/CoronaryArteryTree.ih
)

find_package(Ceres REQUIRED)



# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/ext/" "${PROJECT_SOURCE_DIR}/helpers/" "${PROJECT_SOURCE_DIR}/base/" )
FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS} )

add_subdirectory(tests)
add_subdirectory(bin)

