# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Project name
project(VTK_FileWritingExamples)

# Find the VTK package (make sure VTK is installed and findable)
find_package(VTK REQUIRED)

# Include VTK's CMake module for proper configuration
include_directories(${VTK_INCLUDE_DIRS})

# Add the source files for all examples
add_executable(sphere sphere.cpp)
add_executable(structured_grid structured_grid.cpp)
add_executable(line line.cpp)

# Link VTK libraries to each executable
target_link_libraries(sphere PRIVATE ${VTK_LIBRARIES})
target_link_libraries(structured_grid PRIVATE ${VTK_LIBRARIES})
target_link_libraries(line PRIVATE ${VTK_LIBRARIES})

# Optionally set C++ standard to C++11 or higher
set_target_properties(sphere structured_grid line PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES)
