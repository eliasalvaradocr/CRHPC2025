cmake_minimum_required(VERSION 3.12)

project(lambda LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(lambda main.cpp)

target_link_libraries(lambda Kokkos::kokkos)
