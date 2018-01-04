cmake_minimum_required (VERSION 2.8)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project (gvcfgenotyper C CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")

set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS_RELEASE "-O2")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -Werror")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
        message(FATAL_ERROR "GCC version must be at least 4.7")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

else()
    message(WARNING "You are using an unsupported compiler.")
endif()

add_custom_command(
	OUTPUT  ${PROJECT_SOURCE_DIR}/external/htslib-1.6/libhts.a
	COMMAND make -j4 libhts.a
    	WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}/external/htslib-1.6/
)

add_custom_target(
    libhts.a
    DEPENDS  ${PROJECT_SOURCE_DIR}/external/htslib-1.6/libhts.a
)

include_directories("${PROJECT_SOURCE_DIR}/src/cpp/lib/" "${PROJECT_SOURCE_DIR}/src/c/" "${PROJECT_SOURCE_DIR}/src/cpp/test/")
include_directories("${PROJECT_SOURCE_DIR}/external/htslib-1.6/")

file(GLOB SRC "${PROJECT_SOURCE_DIR}/src/cpp/lib/*.cpp" "${PROJECT_SOURCE_DIR}/src/c/*.c")
add_executable(gvcfgenotyper src/cpp/gvcfgenotyper.cpp ${SRC} )

target_link_libraries(gvcfgenotyper  ${PROJECT_SOURCE_DIR}/external/htslib-1.6/libhts.a)
add_dependencies(gvcfgenotyper libhts.a)


## google test stuff
add_subdirectory(${PROJECT_SOURCE_DIR}/external/googletest-release-1.8.0)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

file(GLOB TEST_SRC "${PROJECT_SOURCE_DIR}/src/test/*.cpp")
add_executable(test_gvcfgenotyper src/cpp/test/test_gvcfgenotyper.cpp ${SRC} ${TEST_SRC})
target_link_libraries(test_gvcfgenotyper  ${PROJECT_SOURCE_DIR}/external/htslib-1.6/libhts.a)
add_dependencies(test_gvcfgenotyper libhts.a)
target_link_libraries(test_gvcfgenotyper gtest gtest_main)


find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( gvcfgenotyper ${ZLIB_LIBRARIES} )
    target_link_libraries( test_gvcfgenotyper ${ZLIB_LIBRARIES} )    
endif( ZLIB_FOUND )

file (COPY ${PROJECT_SOURCE_DIR}/src DESTINATION ${PROJECT_BINARY_DIR})
file (COPY ${PROJECT_SOURCE_DIR}/src/bash DESTINATION ${PROJECT_BINARY_DIR})

enable_testing()

add_test(UNIT_TESTS test_gvcfgenotyper)
add_test(SMOKE_TESTS bash -e bash/run_smoke_tests.sh)
