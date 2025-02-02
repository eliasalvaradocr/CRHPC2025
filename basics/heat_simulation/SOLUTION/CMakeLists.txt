cmake_minimum_required(VERSION 3.12)

project(heat_simulation_kokkos LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(heat_simulation_kokkos main.cpp)

target_link_libraries(heat_simulation_kokkos Kokkos::kokkos)
