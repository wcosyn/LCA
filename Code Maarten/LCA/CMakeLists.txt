cmake_minimum_required(VERSION 2.8)
project(NuclearMomentumDistributions)

set (NuclearMomentumDistributions_VERSION_MAJOR 1)
set (NuclearMomentumDistributions_VERSION_MINOR 0)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

MESSAGE(" -- [info] c++ ${CMAKE_CXX_COMPILER}")
MESSAGE(" -- [info] c ${CMAKE_C_COMPILER}")
MESSAGE(" -- [info] Type ${CMAKE_BUILD_TYPE}")

find_package(GSL REQUIRED)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -fopenmp")

set (LIBS ${LIBS} ${GSL_LIBRARIES})

include_directories("${CMAKE_SOURCE_DIR}/include")
message("-- [info] CXX FLAGS ARE ${CMAKE_CXX_FLAGS}")
add_subdirectory(${CMAKE_SOURCE_DIR}/src)




