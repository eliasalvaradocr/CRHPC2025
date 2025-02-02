cmake_minimum_required(VERSION 3.12)

project(lambda_correct LANGUAGES CXX)

find_package(Kokkos REQUIRED)

add_executable(lambda_correct main.cpp)

target_link_libraries(lambda_correct Kokkos::kokkos)
