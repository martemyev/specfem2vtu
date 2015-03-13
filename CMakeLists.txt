cmake_minimum_required(VERSION 2.8)

project(specfem2vtu)

set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Build type" FORCE)
set(MY_CXX_FLAGS "-std=c++0x -Wall")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")

file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/*.cpp") # .cpp files
file(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/*.hpp") # .hpp files

add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})


