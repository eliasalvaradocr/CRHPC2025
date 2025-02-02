cmake_minimum_required(VERSION 3.10)

project(wave_equation)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find VTK package (adjust version or components as needed)
find_package(VTK REQUIRED)

find_package(Kokkos REQUIRED)

# Include directories (MPI and VTK)
include_directories(${VTK_INCLUDE_DIRS})

# Add executable
add_executable(wave_equation wave_equation.cpp)

# Link MPI libraries
target_link_libraries(wave_equation ${VTK_LIBRARIES})
target_link_libraries(wave_equation Kokkos::kokkos)