#include <Kokkos_Core.hpp>
#include <cmath>
#include <iostream>

#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

// Simulation parameters
const int Nx = 100;       // Number of grid points in x
const int Ny = 100;       // Number of grid points in y
const int Nt = 1000;      // Number of time steps
const double c = 1.0;     // Wave speed
const double dx = 0.01;   // Spatial step in x
const double dy = 0.01;   // Spatial step in y
const double dt = 0.005;  // Time step
const double dt2 = dt * dt;

// Function to write simulation data to a VTK structured grid file
void write_to_vtk(const Kokkos::View<double**>& data, int Nx, int Ny, double dx, double dy, int timestep) {
    // Create a vtkStructuredGrid object
    auto grid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Set the dimensions of the grid
    grid->SetDimensions(Nx, Ny, 1);

    // Create points for the structured grid
    auto points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            points->InsertNextPoint(i * dx, j * dy, 0.0);
        }
    }
    grid->SetPoints(points);

    // Create a vtkDoubleArray to hold the scalar field (wave amplitude)
    auto waveAmplitude = vtkSmartPointer<vtkDoubleArray>::New();
    waveAmplitude->SetName("WaveAmplitude");
    waveAmplitude->SetNumberOfComponents(1);
    waveAmplitude->SetNumberOfTuples(Nx * Ny);

    // Copy the simulation data from Kokkos to the vtkDoubleArray
    //auto data_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), data);
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            waveAmplitude->SetTuple1(i * Ny + j, data(i, j));
        }
    }

    // Add the scalar field to the structured grid
    grid->GetPointData()->SetScalars(waveAmplitude);

    // Write the structured grid to a VTK file
    auto writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();

    // Generate file name
    std::ostringstream filename;
    filename << "wave_" << std::setw(4) << std::setfill('0') << timestep << ".vts";

    writer->SetFileName(filename.str().c_str());
    writer->SetInputData(grid);
    writer->Write();

    std::cout << "VTK file written: " << filename.str() << std::endl;
}

// Add kokkos timer

// Wave equation solver using Kokkos
void wave_equation_solver() {
    using View2D = Kokkos::View<double**, Kokkos::Device<Kokkos::DefaultExecutionSpace,Kokkos::SharedSpace>>;

    // Allocate views for the wave field (current, next, and previous time steps)
    View2D u_prev("u_prev", Nx, Ny);
    View2D u_curr("u_curr", Nx, Ny);
    View2D u_next("u_next", Nx, Ny);

    // Initialize the wave field (e.g., Gaussian pulse at the center)
    Kokkos::parallel_for("Init", Nx * Ny, KOKKOS_LAMBDA(int idx) {
        int i = idx / Ny;
        int j = idx % Ny;
        double x = i * dx - 0.5;
        double y = j * dy - 0.5;
        u_curr(i, j) = exp(-100 * (x * x + y * y)); // Gaussian pulse
        u_prev(i, j) = u_curr(i, j);
    });

    // Time-stepping loop
    for (int t = 0; t < Nt; ++t) {
        // Compute the next time step
        Kokkos::parallel_for("Update", Nx * Ny, KOKKOS_LAMBDA(int idx) {
            int i = idx / Ny;
            int j = idx % Ny;

            if (i > 0 && i < Nx - 1 && j > 0 && j < Ny - 1) {
                double d2u_dx2 = (u_curr(i + 1, j) - 2 * u_curr(i, j) + u_curr(i - 1, j)) / (dx * dx);
                double d2u_dy2 = (u_curr(i, j + 1) - 2 * u_curr(i, j) + u_curr(i, j - 1)) / (dy * dy);
                u_next(i, j) = 2 * u_curr(i, j) - u_prev(i, j) + dt2 * c * c * (d2u_dx2 + d2u_dy2);
            }
        });

        write_to_vtk(u_curr, Nx, Ny, dx, dy, t);
        // Swap views for the next iteration
        Kokkos::deep_copy(u_prev, u_curr);
        Kokkos::deep_copy(u_curr, u_next);
    }

    // Output the final wave field (e.g., for visualization)
    /*auto u_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), u_curr);
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            std::cout << u_host(i, j) << " ";
        }
        std::cout << "\n";
    }*/
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);

    wave_equation_solver();

    Kokkos::finalize();
    return 0;
}
