cmake_minimum_required(VERSION 3.0)
project(test)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread -mavx2 -O3")
add_compile_options(-Wunknown-pragmas)

add_executable(test src/main.cpp)
target_link_libraries(test pbc
        "/lib/x86_64-linux-gnu/libgmpxx.a"
        "/lib/x86_64-linux-gnu/libgmp.a")