cmake_minimum_required (VERSION 2.8)

project (gvcfmerge)

add_custom_command(
	OUTPUT ${PROJECT_SOURCE_DIR}/htslib-1.5/libhts.a
	COMMAND make -j4 libhts.a
    	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/htslib-1.5/
)

add_custom_target(
    libhts.a
    DEPENDS ${PROJECT_SOURCE_DIR}/htslib-1.5/libhts.a
)

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_SOURCE_DIR}/htslib-1.5/")

file(GLOB SRC "${PROJECT_SOURCE_DIR}/src/*")
add_executable(gvcfmerge ${SRC})

target_link_libraries(gvcfmerge ${PROJECT_SOURCE_DIR}/htslib-1.5/libhts.a)
add_dependencies(gvcfmerge libhts.a)


