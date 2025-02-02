cmake_minimum_required(VERSION 3.12)

project(functors LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(functors main.cpp)

target_link_libraries(functors Kokkos::kokkos)
