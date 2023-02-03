cmake_minimum_required(VERSION 3.15.1)

project(GraphAnalysis LANGUAGES Swift)

if(CMAKE_SYSTEM_NAME STREQUAL Windows OR CMAKE_SYSTEM_NAME STREQUAL Darwin)
  option(BUILD_SHARED_LIBS "Build shared libraries by default" YES)
endif()

add_executable(xml2graph SrcXml2graph/main.swift)
add_executable(graph2statBifRad  SrcGraph2statBifRad/main.swift)



