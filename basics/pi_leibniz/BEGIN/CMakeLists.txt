cmake_minimum_required(VERSION 3.12)

project(pi_leibniz_begin LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(pi_leibniz_begin main.cpp)

target_link_libraries(pi_leibniz_begin Kokkos::kokkos)
