CMAKE_MINIMUM_REQUIRED(VERSION 3.6.0 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
project(NuclearMomentumDistributions)

set (NuclearMomentumDistributions_VERSION_MAJOR 1)
set (NuclearMomentumDistributions_VERSION_MINOR 0)
set (nmd_version "${NuclearMomentumDistributions_VERSION_MAJOR}.${NuclearMomentumDistributions_VERSION_MINOR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)


#set defaults
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT STATIC_BUILD)
  set(STATIC_BUILD no)
endif()


MESSAGE(" -- [info] c++ ${CMAKE_CXX_COMPILER}")
MESSAGE(" -- [info] c ${CMAKE_C_COMPILER}")
MESSAGE(" -- [info] Type ${CMAKE_BUILD_TYPE}")

find_package(GSL REQUIRED)
find_package(Doxygen)
find_package(OpenMP REQUIRED)

set(RELEASE_FLAGS -Wall -Wextra -Wunused)
set(DEBUG_FLAGS ${RELEASE_FLAGS})
set(PROFILE_FLAGS ${RELEASE_FLAGS})

if (CMAKE_COMPILER_IS_GNUCXX)
MESSAGE("--- Found GNU compiler collection")
  set(RELEASE_FLAGS ${RELEASE_FLAGS} -march=native -mtune=native -pipe -fopenmp)
  # set(DEBUG_FLAGS ${DEBUG_FLAGS} -g )
  set(PROFILE_FLAGS ${PROFILE_FLAGS} -g -O3)
endif (CMAKE_COMPILER_IS_GNUCXX)

# if(${CMAKE_BUILD_TYPE} MATCHES "Profile")
#   MESSAGE ("----profiling")
#   SET( CMAKE_EXE_LINKER_FLAGS  ${CMAKE_EXE_LINKER_FLAGS} -pg)
# endif(${CMAKE_BUILD_TYPE} MATCHES "Profile")


# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wunused -std=c++11 -fopenmp")
set(LIBS "${LIBS}${GSL_LIBRARIES}")

# make sure our executables will find dynamical library:
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
set(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH 1)


set(CMAKE_EXPORT_COMPILE_COMMANDS ON) #to create compile_commands.json file used in VS code

set(CMAKE_VERBOSE_MAKEFILE ON)


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

