cmake_minimum_required(VERSION 3.12)

project(matmul_begin LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(matmul_begin main.cpp)

target_link_libraries(matmul_begin Kokkos::kokkos)
