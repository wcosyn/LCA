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
	kinenergy_ob.cpp
	newcoef.cpp
	norm_ob.cpp
	norm_tb.cpp
	nucleusall.cpp
	nucleus.cpp
	nucleusnn.cpp
	nucleusnp.cpp
	nucleuspp.cpp
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

add_executable(mainwswf mainwswf.cpp ${SRC})
target_link_libraries(mainwswf ${LIBS})

add_executable(mainop mainop.cpp ${SRC})
target_link_libraries(mainop ${LIBS})

add_library(nmd SHARED ${SRC})
target_link_libraries(nmd "${GSL_LIBRARIES}")
set_target_properties(nmd 
    PROPERTIES
        VERSION   "${nmd_version}"
        SOVERSION "${nmd_version}"
        SUFFIX    ".so"
)

install(TARGETS mainwswf mainop nmd
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
        LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)
