

add_executable(threej_test threej_test.cpp)
target_link_libraries(threej_test ${LIBS})

add_executable(nucleus_test nucleus_test.cpp)
target_link_libraries(nucleus_test ${LIBS})

add_executable(normob_test normob_test.cpp)
target_link_libraries(normob_test ${LIBS})

add_executable(ho_test hotest.cpp)
target_link_libraries(ho_test ${LIBS})

add_executable(newcoef_test newcoef_test.cpp)
target_link_libraries(newcoef_test ${LIBS})

add_executable(hoft_test hoft_test.cpp)
target_link_libraries(hoft_test ${LIBS})


install(TARGETS threej_test nucleus_test normob_test ho_test newcoef_test hoft_test
	RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

