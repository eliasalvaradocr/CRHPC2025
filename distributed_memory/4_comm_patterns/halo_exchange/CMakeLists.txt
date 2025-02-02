cmake_minimum_required(VERSION 3.10)

project(MPI_Halo_Exchange)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find MPI package
find_package(MPI REQUIRED)

# Find VTK package (adjust version or components as needed)
find_package(VTK REQUIRED)

# Include directories (MPI and VTK)
include_directories(${MPI_INCLUDE_PATH} ${VTK_INCLUDE_DIRS})

# Add executable
add_executable(mpi_halo_exchange write_halo.h halo_exchange.cpp)

# Link MPI libraries
target_link_libraries(mpi_halo_exchange ${MPI_LIBRARIES})
target_link_libraries(mpi_halo_exchange ${VTK_LIBRARIES})

# Optional: Set the compiler flags for optimization or debugging
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Optional: Set the MPI include directories if MPI_INCLUDE_PATH isn't working
# include_directories(/path/to/mpi/include)
# link_directories(/path/to/mpi/lib)
