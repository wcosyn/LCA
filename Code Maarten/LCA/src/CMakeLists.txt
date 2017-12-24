set( SRC
	correlation_functions.cpp
	density_angle.cpp
	density_ob3.cpp
	density_ob_integrand3.cpp
	density_ob_integrand_cf.cpp
	density_ob_integrand_pow.cpp
	density_rel.cpp
	density_rel_integrand2.cpp
	density_td.cpp
	isopaircoef.cpp
	kinenergy_ob.cpp
	newcoef.cpp
	norm_ob.cpp
	norm_tb.cpp
	nucleusall.cpp
	nucleus.cpp
	nucleusnn.cpp
	nucleusnp.cpp
	nucleuspp.cpp
	nucleus_iso.cpp
	operator_virtual.cpp
	operator_virtual_ob.cpp
	paircoef.cpp
	pair.cpp
	recmosh.cpp
	rms_ob.cpp
	shell.cpp
	speedy.cpp
	threebodycoef.cpp
	threej.cpp
	tripletcoef.cpp
	triplet.cpp
	wavefunctionp.cpp
	wsexpansion.cpp
	wsnucleus.cpp
	wsnucleusnp.cpp
	wsnucleuspp.cpp
	wspair.cpp
	wswf.cpp )

add_executable(mainwswf mainwswf.cpp)
target_link_libraries(mainwswf ${LIBS} nmd)
target_compile_options(mainwswf PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_FLAGS}>")
target_compile_options(mainwswf PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_FLAGS}>")
target_compile_options(mainwswf PUBLIC "$<$<CONFIG:PROFILE>:${PROFILE_FLAGS}>")

add_executable(mainop mainop.cpp)
target_link_libraries(mainop ${LIBS} nmd)
target_compile_options(mainop PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_FLAGS}>")
target_compile_options(mainop PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_FLAGS}>")
target_compile_options(mainop PUBLIC "$<$<CONFIG:PROFILE>:${PROFILE_FLAGS}>")

add_library(nmd SHARED ${SRC})
target_include_directories(nmd PUBLIC "${CMAKE_SOURCE_DIR}/include" ${GSL_INCLUDE_DIRS})
target_link_libraries(nmd "${GSL_LIBRARIES}" ${OpenMP_CXX_FLAGS})
target_compile_options(nmd PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_FLAGS}>")
target_compile_options(nmd PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_FLAGS}>")
target_compile_options(nmd PUBLIC "$<$<CONFIG:PROFILE>:${PROFILE_FLAGS}>")
set_target_properties(nmd 
    PROPERTIES
        VERSION   "${nmd_version}"
        SOVERSION "${nmd_version}"
        SUFFIX    ".so"
)

install(TARGETS mainwswf mainop nmd
RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)



if(STATIC_BUILD)
		add_library(nmdST STATIC ${SRC})
		target_compile_options(nmdST PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_FLAGS}>")
		target_compile_options(nmdST PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_FLAGS}>")
		target_compile_options(nmdST PUBLIC "$<$<CONFIG:PROFILE>:${PROFILE_FLAGS}>")
		target_include_directories(nmdST PUBLIC ${GSL_INCLUDE_DIRS})
		target_link_libraries(nmdST PUBLIC ${GSL_LIBRARIES} -Wl,-no-whole-archive ${OpenMP_CXX_FLAGS})
		set_target_properties(nmdST PROPERTIES OUTPUT_NAME nmd)
		install(TARGETS nmdST DESTINATION "${CMAKE_SOURCE_DIR}/lib")


		add_executable(mainwswfST mainwswf.cpp)
		target_link_libraries(mainwswfST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)
		target_compile_options(mainwswfST PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_FLAGS}>")
		target_compile_options(mainwswfST PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_FLAGS}>")
		target_compile_options(mainwswfST PUBLIC "$<$<CONFIG:PROFILE>:${PROFILE_FLAGS}>")
		
		add_executable(mainopST mainop.cpp)
		target_link_libraries(mainopST PUBLIC "-static" nmdST -Wl,-no-whole-archive)
		target_compile_options(mainopST PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_FLAGS}>")
		target_compile_options(mainopST PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_FLAGS}>")
		target_compile_options(mainopST PUBLIC "$<$<CONFIG:PROFILE>:${PROFILE_FLAGS}>")

		install(TARGETS mainwswfST mainopST DESTINATION ${CMAKE_SOURCE_DIR}/bin)

endif(STATIC_BUILD)
	  