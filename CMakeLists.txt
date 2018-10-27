cmake_minimum_required(VERSION 3.12)
project(cvtree)
#set(CMAKE_GENERATOR "Unix")
#message("generator is set to ${CMAKE_GENERATOR}")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc -fopenmp")
set(arch_hint "x64")
include_directories(${INCLUDE_DIR})
add_executable(cvtree
     parallel.cpp)


