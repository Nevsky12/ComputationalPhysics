set(MODULE_NAME closure)

add_library(Eigen INTERFACE)
target_include_directories(Eigen INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(${LIB_NAME} PROPERTIES LINKER_LANGUAGE CXX)