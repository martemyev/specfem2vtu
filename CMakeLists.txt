cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(specfem2vtu CXX)

set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Build type" FORCE)
set(MY_CXX_FLAGS "-std=c++0x -Wall")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")

file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/*.cpp") # .cpp files
file(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/*.hpp") # .hpp files

# Triangle software
set(Triangle_DIR Triangle)
include_directories(${Triangle_DIR})
add_subdirectory(${Triangle_DIR})
#message("Triangle_LIB = ${Triangle_LIB}")

#find_library(Triangle_LIB triangle PATHS ${Triangle_DIR} PATH_SUFFIXES build)
#if(Triangle_LIB)
#  message(STATUS "triangle lib: ${Triangle_LIB}")
#else()
#  message(FATAL_ERROR "triangle library was not found")
#endif()

add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC triangle)

