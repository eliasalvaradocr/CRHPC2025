cmake_minimum_required(VERSION 3.12)

project(views LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(views main.cpp)

target_link_libraries(views Kokkos::kokkos)
