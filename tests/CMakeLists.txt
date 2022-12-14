file(GLOB_RECURSE TEST_SOURCES test*.cpp)
set(LINK_LIBS ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(googletest)

enable_testing()

message("Testing started!")

foreach (TEST_FILE IN LISTS TEST_SOURCES)

    get_filename_component(TEST_FILE_NAME ${TEST_FILE} NAME_WE)
    set(EXECUTABLE_NAME run_${TEST_FILE_NAME})
    add_executable(${EXECUTABLE_NAME} ${TEST_FILE})
    target_link_libraries(${EXECUTABLE_NAME} Eigen gtest gtest_main computational_module)
    add_test(NAME ${EXECUTABLE_NAME} COMMAND ${EXECUTABLE_NAME})

endforeach ()

message("Testing ended!")