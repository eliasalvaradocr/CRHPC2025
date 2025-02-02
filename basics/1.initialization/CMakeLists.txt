cmake_minimum_required(VERSION 3.12)

project(initialize_kokkos LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(initialization main.cpp)

target_link_libraries(initialization Kokkos::kokkos)
