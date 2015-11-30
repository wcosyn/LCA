

add_executable(threej_test threej_test.cpp ${CMAKE_SOURCE_DIR}/src/threej.cpp)
target_link_libraries(threej_test ${LIBS})

add_executable(nucleus_test nucleus_test.cpp)
target_link_libraries(nucleus_test ${LIBS})

install(TARGETS threej_test nucleus_test
	RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

