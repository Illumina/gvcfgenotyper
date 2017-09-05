cmake_minimum_required (VERSION 2.8)

project (gvcfmerge C CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")

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

file(GLOB SRC "${PROJECT_SOURCE_DIR}/src/*.cpp" "${PROJECT_SOURCE_DIR}/src/*.c")
add_executable(gvcfmerge ${SRC})

target_link_libraries(gvcfmerge ${PROJECT_SOURCE_DIR}/htslib-1.5/libhts.a)
add_dependencies(gvcfmerge libhts.a)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( gvcfmerge ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )


