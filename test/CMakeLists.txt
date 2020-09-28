

add_executable(threej_test threej_test.cpp)
target_link_libraries(threej_test ${LIBS} nmd)

add_executable(nucleus_test nucleus_test.cpp)
target_link_libraries(nucleus_test ${LIBS} nmd)

add_executable(ho_test hotest.cpp)
target_link_libraries(ho_test ${LIBS} nmd)

add_executable(newcoef_test newcoef_test.cpp)
target_link_libraries(newcoef_test ${LIBS} nmd)

add_executable(hoft_test hoft_test.cpp)
target_link_libraries(hoft_test ${LIBS} nmd)

install(TARGETS threej_test nucleus_test ho_test newcoef_test hoft_test
	RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

if(STATIC_BUILD)
	add_executable(threej_testST threej_test.cpp)
	target_link_libraries(threej_testST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(nucleus_testST nucleus_test.cpp)
	target_link_libraries(nucleus_testST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(ho_testST hotest.cpp)
	target_link_libraries(ho_testST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(newcoef_testST newcoef_test.cpp)
	target_link_libraries(newcoef_testST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	add_executable(hoft_testST hoft_test.cpp)
	target_link_libraries(hoft_testST PUBLIC "-static" ${LIBS} nmdST -Wl,-no-whole-archive)

	install(TARGETS threej_testST nucleus_testST ho_testST newcoef_testST hoft_testST
	RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

endif(STATIC_BUILD)
	