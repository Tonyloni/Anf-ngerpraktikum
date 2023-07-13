cmake_minimum_required(VERSION 3.13)
project(test VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE RELEASE)


add_library(my_lib STATIC argtable3.c graph_ch.cpp graph.cpp graph_ch2.cpp adressable_pq.cpp minheap.cpp dijkstra.cpp)

add_executable(contract contract.cpp)
add_executable(graph main.cpp)
add_executable(tests tests.cpp)

target_link_libraries(contract PRIVATE my_lib)
target_link_libraries(graph PRIVATE my_lib)
target_link_libraries(tests PRIVATE my_lib)




