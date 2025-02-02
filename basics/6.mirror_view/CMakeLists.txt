cmake_minimum_required(VERSION 3.12)

project(mirror_view LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(mirror_view main.cpp)

target_link_libraries(mirror_view Kokkos::kokkos)
