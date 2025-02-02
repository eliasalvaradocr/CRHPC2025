/* @file main.cpp
 * @brief This program simulates a simple heat distribution in a 2D grid
 *        using the Kokkos library.
 */

#include <Kokkos_Core.hpp>
#include <fstream>
#include <iostream>

/*
 * @brief This function initializes the heat distribution in a 2D grid.
 * @param[in] grid The 2D grid where heat distribution will be initialized.
 * @param[in] nx Number of rows in the grid.
 * @param[in] ny Number of columns in the grid.
 */
void initialize_grid(Kokkos::View<double **> grid, const int nx, const int ny) {

  Kokkos::parallel_for(
      "initialize_grid", nx * ny, KOKKOS_LAMBDA(const int idx) {
        int i = idx / ny;
        int j = idx % ny;
        if ((i == 0 || i == nx - 1) || (j == 0 || j == ny - 1)) {
          grid(i, j) = 100.0; // Boundary condition: set to 100
        } else {
          grid(i, j) = 0.0; // Initial temperature inside the grid is 0
        }
      });
}

/*
 * @brief This function performs one iteration of heat diffusion on a 2D grid.
 * @param[in] current The current state of the heat distribution in the grid.
 * @param[out] next The next state of the heat distribution after diffusion.
 * @param[in] nx Number of rows in the grid.
 * @param[in] ny Number of columns in the grid.
 */
void diffuse_heat(Kokkos::View<const double **> current,
                  Kokkos::View<double **> next, const int nx, const int ny) {
  constexpr double grid_spacing = 0.5;
  constexpr double delta_t = 0.05;
  constexpr double factor = delta_t / (grid_spacing * grid_spacing);
  const auto grid_policy = Kokkos::MDRangePolicy({1, 1}, {nx - 1, ny - 1});
  Kokkos::parallel_for(
      "Diffuse Heat", grid_policy, KOKKOS_LAMBDA(const int i, const int j) {
        next(i, j) =
            current(i, j) + factor * (current(i + 1, j) + current(i - 1, j) +
                                      current(i, j - 1) + current(i, j + 1) -
                                      4 * current(i, j));
      });
  // Apply boundary conditions: keep the edges at 100 degrees
  Kokkos::parallel_for(
      "Apply Boundary Conditions", Kokkos::RangePolicy<>(0, nx),
      KOKKOS_LAMBDA(const int i) {
        next(i, 0) = 100.0;
        next(i, ny - 1) = 100.0;
      });
  Kokkos::parallel_for(
      "Apply Boundary Conditions", Kokkos::RangePolicy<>(0, ny),
      KOKKOS_LAMBDA(const int j) {
        next(0, j) = 100.0;
        next(nx - 1, j) = 100.0;
      });
}

/*
 * @brief This function outputs a vtk file for visualization.
 *
 * @param current The current state of the heat distribution.
 * @param nx Number of grid points in x-direction.
 * @param ny Number of grid points in y-direction.
 * @param step Current time step
 */
void output_vtk(Kokkos::View<const double **> current, const int nx,
                const int ny, const int step) {
  std::ofstream file;
  std::string filename = "heat_" + std::to_string(step) + ".vtk";
  auto host_current =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), current);
  file.open(filename);
  if (file.is_open()) {
    file << "# vtk DataFile Version 3.0\n";
    file << "Heat distribution\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "ORIGIN 0 0 0\n";
    file << "SPACING 1 1 1\n";
    file << "POINT_DATA " << nx * ny << "\n";
    file << "SCALARS temperature float\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny;
         ++j) { // Loop over y direction first for VTK file format
      for (int i = 0; i < nx; ++i) {        // Then loop over x direction
        file << host_current(j, i) << "\n"; // Output temperature values
      }
    }
    file.close();
  } else {
    std::cerr << "Unable to open file for writing\n";
  }
}

// Main function to simulate heat distribution over time
int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  std::cout << "******************\n";
  std::cout << "Heat Simulation\n";
  std::cout << "******************\n";
  int nx = 100;    // Number of grid points in x direction
  int ny = 100;    // Number of grid points in y direction
  int steps = 500; // Total number of time steps
  int output_freq = 1000;
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0]
              << " <nx> <ny> <steps> <output frequency>\n";
    exit(-1);
  } else {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    steps = atoi(argv[3]);
    output_freq = atoi(argv[4]);
  }
  {
    Kokkos::View<double **> heat_map("heat_map", nx, ny);
    Kokkos::View<double **> new_heat_map("new_heat_map", nx, ny);
    // Initialize the heat distribution in the grid
    std::cout << "Initializing grid...\n";
    initialize_grid(heat_map, nx, ny);

    std::cout << "Starting simulation...\n";
    // Perform the heat simulation for the specified number of steps
    auto timer = Kokkos::Timer();
    double total_step_time = 0.;
    double total_io_time = 0.;
    for (int step = 0; step < steps; ++step) {
      timer.reset();
      diffuse_heat(heat_map, new_heat_map, nx, ny);
      total_step_time += timer.seconds();
      // Swap the views to use the updated heat distribution in the next
      // iteration
      Kokkos::kokkos_swap(heat_map, new_heat_map);
      if (step % output_freq == 0) {
        std::cout << "Writing step " << step << '\n';
        timer.reset();
        output_vtk(heat_map, nx, ny, step);
        total_io_time += timer.seconds();
      }
    }
    auto average_step_time = total_step_time / steps;
    auto average_io_time = total_io_time / steps;
    std::cout << "Simulation complete...\n";
    std::cout << "Total step time: " << total_step_time << '\n';
    std::cout << "Average step time: " << average_step_time << '\n';
    std::cout << "Total I/O time: " << total_io_time << '\n';
    std::cout << "Average I/O time: " << average_io_time << '\n';
  }
  Kokkos::finalize();
  return 0;
}
