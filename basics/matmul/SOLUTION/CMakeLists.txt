cmake_minimum_required(VERSION 3.12)

project(matmul_solution LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(matmul_solution main.cpp)

target_link_libraries(matmul_solution Kokkos::kokkos)
