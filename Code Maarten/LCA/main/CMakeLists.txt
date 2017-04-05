

add_executable(normob normob.cpp)
target_link_libraries(normob ${LIBS})

add_executable(obmd obmd.cpp)
target_link_libraries(obmd ${LIBS})

add_executable(normtb normtb.cpp)
target_link_libraries(normtb ${LIBS})

add_executable(tbmd_rel tbmd_rel.cpp)
target_link_libraries(tbmd_rel ${LIBS})

install(TARGETS normob obmd normtb tbmd_rel
	RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

