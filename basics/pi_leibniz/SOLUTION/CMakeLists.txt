cmake_minimum_required(VERSION 3.12)

project(pi_leibniz_solution LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(pi_leibniz_solution main.cpp)

target_link_libraries(pi_leibniz_solution Kokkos::kokkos)
