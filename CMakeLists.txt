cmake_minimum_required (VERSION 2.8)

project (gvcfmerge)

include_directories("${PROJECT_BINARY_DIR}")

add_custom_command(
    OUTPUT  ${PROJECT_SOURCE_DIR}/htslib-1.5/libhts.a
    WORKING_DIR ${PROJECT_SOURCE_DIR}/htslib-1.5/
    COMMAND make
)

include_directories("${PROJECT_BINARY_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}/htslib-1.5/")

add_executable(gvcfmerge gvcfmerge.cpp)
target_link_libraries(gvcfmerge ${PROJECT_SOURCE_DIR}/htslib-1.5/libhts.a)
