cmake_minimum_required (VERSION 2.8)

project (gvcfmerge C CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")

add_custom_command(
	OUTPUT  ${PROJECT_SOURCE_DIR}/external/htslib-1.5/libhts.a
	COMMAND make -j4 libhts.a
    	WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}/external/htslib-1.5/
)

add_custom_target(
    libhts.a
    DEPENDS  ${PROJECT_SOURCE_DIR}/external/htslib-1.5/libhts.a
)

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/external/htslib-1.5/")

file(GLOB SRC "${PROJECT_SOURCE_DIR}/src/lib/*.cpp" "${PROJECT_SOURCE_DIR}/src/lib/*.c" "${PROJECT_SOURCE_DIR}/src/lib/tools/*.cpp")
add_executable(gvcfmerge src/gvcfmerge.cpp ${SRC})

target_link_libraries(gvcfmerge  ${PROJECT_SOURCE_DIR}/external/htslib-1.5/libhts.a)
add_dependencies(gvcfmerge libhts.a)


## google test stuff
add_subdirectory(${PROJECT_SOURCE_DIR}/external/googletest-release-1.8.0)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

file(GLOB TEST_SRC "${PROJECT_SOURCE_DIR}/src/test/*.cpp")
add_executable(test_gvcfmerge src/test_gvcfmerge.cpp ${SRC} ${TEST_SRC}) 
target_link_libraries(test_gvcfmerge  ${PROJECT_SOURCE_DIR}/external/htslib-1.5/libhts.a)
add_dependencies(test_gvcfmerge libhts.a)
target_link_libraries(test_gvcfmerge gtest gtest_main)


find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( gvcfmerge ${ZLIB_LIBRARIES} )
    target_link_libraries( test_gvcfmerge ${ZLIB_LIBRARIES} )    
endif( ZLIB_FOUND )

file (COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR})

