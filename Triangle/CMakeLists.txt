message("Entering Triangle subdirectory")

enable_language(C)
set(MY_C_FLAGS "-O3 -Wall -DTRILIBRARY") # -DANSI_DECLARATORS")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")
message(STATUS "Triangle: compilation flags ${CMAKE_C_FLAGS}")

#set(Triangle_LIB triangle)
#add_library(${Triangle_LIB} triangle.c triangle.h)
add_library(triangle triangle.c triangle.h)
target_include_directories(triangle PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

message("Exiting Triangle subdirectory")
