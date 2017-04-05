cmake_minimum_required(VERSION 2.8)
project(NuclearMomentumDistributions)

set (NuclearMomentumDistributions_VERSION_MAJOR 1)
set (NuclearMomentumDistributions_VERSION_MINOR 0)
set (nmd_version "${NuclearMomentumDistributions_VERSION_MAJOR}.${NuclearMomentumDistributions_VERSION_MINOR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

MESSAGE(" -- [info] c++ ${CMAKE_CXX_COMPILER}")
MESSAGE(" -- [info] c ${CMAKE_C_COMPILER}")
MESSAGE(" -- [info] Type ${CMAKE_BUILD_TYPE}")

find_package(GSL REQUIRED)
find_package(Doxygen)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wunused -std=c++11 -fopenmp")
set(LIBS "${LIBS}${GSL_LIBRARIES}" nmd)

# make sure our executables will find dynimcal library:
set( CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/lib")



include_directories("${CMAKE_SOURCE_DIR}/include")
message("-- [info] CXX FLAGS ARE ${CMAKE_CXX_FLAGS}")
add_subdirectory(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(${CMAKE_SOURCE_DIR}/test)
add_subdirectory(${CMAKE_SOURCE_DIR}/main)


## Doxygen configuration
# configure_file replaces variables referenced with @VAR@ with their values defined in CMAKE
# custom target will make a doc target that we can use as
# $> make doc
##
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.in ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target( doc
        COMMAND "${DOXYGEN_EXECUTABLE}" "${CMAKE_BINARY_DIR}/Doxyfile"
        SOURCES "${CMAKE_BINARY_DIR}/Doxyfile"
        COMMENT " -- [info] Building DOXYGEN info" VERBATIM
        )
endif()

