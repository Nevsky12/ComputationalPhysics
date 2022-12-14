cmake_minimum_required(VERSION 3.23)
project(computationalPhysics)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a -pthread")
set(CMAKE_CXX_STANDARD 20)
message(${CMAKE_CXX_FLAGS})

add_subdirectory(src)
add_subdirectory(third_party)
include_directories(third_party)

enable_testing()
add_subdirectory(tests)
