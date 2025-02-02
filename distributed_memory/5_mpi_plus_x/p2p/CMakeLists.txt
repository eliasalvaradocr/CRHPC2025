cmake_minimum_required(VERSION 3.10)

project(P2PKokkos)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find MPI package
find_package(MPI REQUIRED)
find_package(Kokkos REQUIRED)

# Include directories (MPI and VTK)
include_directories(${MPI_INCLUDE_PATH})

# Add executable
add_executable(p2p p2p.cpp)

# Link MPI libraries
target_link_libraries(p2p ${MPI_LIBRARIES})
target_link_libraries(p2p Kokkos::kokkos) # Link Kokkos