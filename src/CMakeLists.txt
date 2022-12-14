set(MODULE_NAME computational_module)

file(GLOB_RECURSE SOURCES *.h *.cpp *.hpp)

set(TARGET_LIBS Eigen)

add_library(${MODULE_NAME} STATIC ${SOURCES})

set_target_properties(${MODULE_NAME} PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(${MODULE_NAME} PUBLIC ${TARGET_LIBS})
target_include_directories(${MODULE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})