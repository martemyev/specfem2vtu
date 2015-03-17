message("Entering Triangle subdirectory")

enable_language(C)
set(MY_C_FLAGS "-Wall -DTRILIBRARY")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")

if(BUILD_TYPE MATCHES DEBUG)
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${MY_C_FLAGS}")
elseif(BUILD_TYPE MATCHES RELEASE)
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${MY_C_FLAGS}")
endif()

message("\n------------ Triangle configuration information ------------")
message("build type     " ${BUILD_TYPE})
message("compiler       " ${CMAKE_C_COMPILER})
if(BUILD_TYPE MATCHES DEBUG)
  message("compiler flags " ${CMAKE_C_FLAGS_DEBUG})
elseif(BUILD_TYPE MATCHES RELEASE)
  message("compiler flags " ${CMAKE_C_FLAGS_RELEASE})
endif()

add_library(triangle triangle.c triangle.h)
target_include_directories(triangle PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

message("Exiting Triangle subdirectory")
