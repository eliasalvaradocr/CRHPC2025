cmake_minimum_required(VERSION 3.12)

project(KMeans_solution LANGUAGES CXX)

find_package(Kokkos)

add_executable(kmeans_solution main.cpp)
target_link_libraries(kmeans_solution Kokkos::kokkos)
