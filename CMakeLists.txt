cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(specfem2vtu CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build: Debug | Release" FORCE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)

set(MY_CXX_FLAGS "-std=c++0x -Wall")

if(BUILD_TYPE MATCHES DEBUG)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${MY_CXX_FLAGS}")
elseif(BUILD_TYPE MATCHES RELEASE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MY_CXX_FLAGS}")
endif()

file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/*.cpp") # .cpp files
file(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/*.hpp") # .hpp files

# Triangle software
set(Triangle_DIR Triangle)
include_directories(${Triangle_DIR})
add_subdirectory(${Triangle_DIR})

message("\n------------ specfem2vtu configuration information ------------")
message("build type     " ${BUILD_TYPE})
message("compiler       " ${CMAKE_CXX_COMPILER})
if(BUILD_TYPE MATCHES DEBUG)
  message("compiler flags " ${CMAKE_CXX_FLAGS_DEBUG})
elseif(BUILD_TYPE MATCHES RELEASE)
  message("compiler flags " ${CMAKE_CXX_FLAGS_RELEASE})
endif()

add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC triangle)

