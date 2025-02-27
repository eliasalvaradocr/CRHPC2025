cmake_minimum_required(VERSION 3.10)

project(findpi)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find VTK package (adjust version or components as needed)
find_package(VTK REQUIRED)
find_package(MPI REQUIRED)
find_package(Kokkos REQUIRED)

# Include directories (MPI and VTK)
include_directories(${VTK_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_DIRS})

# Add executable
add_executable(findpi findpi.cpp)

# Link MPI libraries
target_link_libraries(findpi ${VTK_LIBRARIES})
target_link_libraries(findpi ${MPI_LIBRARIES})
target_link_libraries(findpi Kokkos::kokkos)