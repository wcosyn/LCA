
add_executable(obwd obwd.cpp)
target_link_libraries(obwd ${LIBS} nmd)

add_executable(normob normob.cpp)
target_link_libraries(normob ${LIBS} nmd)

add_executable(normob_iso normob_iso.cpp)
target_link_libraries(normob_iso ${LIBS} nmd)

add_executable(obmd obmd.cpp)
target_link_libraries(obmd ${LIBS} nmd)

add_executable(rms rms.cpp)
target_link_libraries(rms ${LIBS} nmd)

add_executable(rms_iso rms_iso.cpp)
target_link_libraries(rms_iso ${LIBS} nmd)

add_executable(rms_calc rms_calc.cpp)
target_link_libraries(rms_calc ${LIBS} nmd)

add_executable(normtb normtb.cpp)
target_link_libraries(normtb ${LIBS} nmd)

add_executable(tbmd_rel tbmd_rel.cpp)
target_link_libraries(tbmd_rel ${LIBS} nmd)

#add_executable(fitcode fitcode.cpp)
#target_link_libraries(fitcode ${LIBS} nmd)

add_executable(fitcode_NZ fitcode_NZ.cpp)
target_link_libraries(fitcode_NZ ${LIBS} nmd)

install(TARGETS obwd rms rms_iso rms_calc normob normob_iso obmd normtb tbmd_rel fitcode_NZ
	RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

	if(STATIC_BUILD)
	add_executable(normobST normob.cpp)
	target_link_libraries(normobST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(normob_isoST normob_iso.cpp)
	target_link_libraries(normob_isoST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(rmsST rms.cpp)
	target_link_libraries(rmsST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(rms_isoST rms_iso.cpp)
	target_link_libraries(rms_isoST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(rms_calcST rms_calc.cpp)
	target_link_libraries(rms_calcST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(obmdST obmd.cpp)
	target_link_libraries(obmdST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(normtbST normtb.cpp)
	target_link_libraries(normtbST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(tbmd_relST tbmd_rel.cpp)
	target_link_libraries(tbmd_relST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(fitcodeST fitcode.cpp)
	target_link_libraries(fitcodeST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	install(TARGETS rmsST rms_isoST rms_calcST normobST normob_isoST obmdST normtbST tbmd_relST fitcodeST
	RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

endif(STATIC_BUILD)

