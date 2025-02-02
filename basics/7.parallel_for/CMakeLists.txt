cmake_minimum_required(VERSION 3.12)

project(parallel_for LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(parallel_for main.cpp)

target_link_libraries(parallel_for Kokkos::kokkos)
