#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h> //
#include <vtkPolyLine.h> //
#include <vtkXMLPolyDataWriter.h>

#include <cstdlib>
#include <ctime>

#define PI 3.141592653589793238462643

int main(int argv, char* argc[]) {
  Kokkos::initialize(argv, argc);
  {
    Kokkos::Random_XorShift64_Pool<> random_pool(/*seed=*/12345);

    // Parameters
    const int numParticles = 1000;   // Number of particles
    const int numSteps = 1000;       // Steps per particle
    const double R = 0.001;         // Step size
    const int n_write = 2;         // write frequency

    // 1. Kokkos View to Store All Particles' Random Walks
    Kokkos::View<double**, Kokkos::Device<Kokkos::DefaultExecutionSpace,Kokkos::SharedSpace>> positions("positions", numParticles, 3);

    // Initialize the starting point for all particles
    Kokkos::parallel_for("initialize_particles", numParticles, KOKKOS_LAMBDA(const int i) {
        positions(i, 0) = 0.0; // x
        positions(i, 1) = 0.0; // y
        positions(i, 2) = 0.0; // z
    });

    for (int i = 0; i < numSteps; i++)
    {

      std::cout << "Step: " << i <<  std::endl; 

      if (i % n_write == 0)
      {
        // VTK Routines for plotting
        vtkNew<vtkPoints> points;

        /*
        positions_gpu
        positions_cpu

        Kokkos::deep_copy(positions_cpu, positions_gpu);
        */
        
        // Iterate over particles to create VTK structures
        for (int j = 0; j < numParticles; ++j) {
          points->InsertNextPoint(
              positions(j, 0), 
              positions(j, 1), 
              positions(j, 2)
          );
        }
        
        std::ostringstream filename;
        filename << "random_walk_" << i << ".vtp";

        vtkNew<vtkPolyData> polyData;
        polyData->SetPoints(points);
        vtkNew<vtkXMLPolyDataWriter> writer;
        writer->SetFileName(filename.str().c_str());
        writer->SetInputData(polyData);
        writer->Write();
      }

      // Generate random walks for all particles
      Kokkos::parallel_for("random_walk_particles", numParticles, KOKKOS_LAMBDA(const int i) {
        auto generator = random_pool.get_state();
        
        double theta = PI*generator.drand(0., 1.);
        double phi = 2*PI*generator.drand(0., 1.);

        double dx = R*sin(theta)*cos(phi);
        double dy = R*sin(theta)*sin(phi);
        double dz = R*cos(theta);

        // do not forget to release the state of the engine
        random_pool.free_state(generator);

        // Compute next position
        positions(i, 0) += dx;
        positions(i, 1) += dy;
        positions(i, 2) += dz;
      });
    }
  }
  Kokkos::finalize();
  return 0;
}
