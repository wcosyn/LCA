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
find_package(Doxygen)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wunused -O0 -std=c++11 -fopenmp")

set (LIBS ${LIBS} ${GSL_LIBRARIES})

set( SRC 
	${CMAKE_SOURCE_DIR}/src/correlation_functions.cpp
	${CMAKE_SOURCE_DIR}/src/density_angle.cpp
	${CMAKE_SOURCE_DIR}/src/density_ob3.cpp
	${CMAKE_SOURCE_DIR}/src/density_ob_integrand3.cpp
	${CMAKE_SOURCE_DIR}/src/density_ob_integrand_cf.cpp
	${CMAKE_SOURCE_DIR}/src/density_ob_integrand_pow.cpp
	${CMAKE_SOURCE_DIR}/src/density_rel.cpp
	${CMAKE_SOURCE_DIR}/src/density_rel_integrand2.cpp
	${CMAKE_SOURCE_DIR}/src/density_td.cpp
	${CMAKE_SOURCE_DIR}/src/kinenergy_ob.cpp
	${CMAKE_SOURCE_DIR}/src/newcoef.cpp
	${CMAKE_SOURCE_DIR}/src/norm_ob.cpp
	${CMAKE_SOURCE_DIR}/src/norm_tb.cpp
	${CMAKE_SOURCE_DIR}/src/nucleusall.cpp
	${CMAKE_SOURCE_DIR}/src/nucleus.cpp
	${CMAKE_SOURCE_DIR}/src/nucleusnn.cpp
	${CMAKE_SOURCE_DIR}/src/nucleusnp.cpp
	${CMAKE_SOURCE_DIR}/src/nucleuspp.cpp
	${CMAKE_SOURCE_DIR}/src/operator_virtual.cpp
	${CMAKE_SOURCE_DIR}/src/operator_virtual_ob.cpp
	${CMAKE_SOURCE_DIR}/src/paircoef.cpp
	${CMAKE_SOURCE_DIR}/src/pair.cpp
	${CMAKE_SOURCE_DIR}/src/recmosh.cpp
	${CMAKE_SOURCE_DIR}/src/rms_ob.cpp
	${CMAKE_SOURCE_DIR}/src/shell.cpp
	${CMAKE_SOURCE_DIR}/src/speedy.cpp
	${CMAKE_SOURCE_DIR}/src/threebodycoef.cpp
	${CMAKE_SOURCE_DIR}/src/threej.cpp
	${CMAKE_SOURCE_DIR}/src/tripletcoef.cpp
	${CMAKE_SOURCE_DIR}/src/triplet.cpp
	${CMAKE_SOURCE_DIR}/src/wavefunctionp.cpp
	${CMAKE_SOURCE_DIR}/src/wsexpansion.cpp
	${CMAKE_SOURCE_DIR}/src/wsnucleus.cpp
	${CMAKE_SOURCE_DIR}/src/wsnucleusnp.cpp
	${CMAKE_SOURCE_DIR}/src/wsnucleuspp.cpp
	${CMAKE_SOURCE_DIR}/src/wspair.cpp
	${CMAKE_SOURCE_DIR}/src/wswf.cpp )

include_directories("${CMAKE_SOURCE_DIR}/include")
message("-- [info] CXX FLAGS ARE ${CMAKE_CXX_FLAGS}")
add_subdirectory(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(${CMAKE_SOURCE_DIR}/test)


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

