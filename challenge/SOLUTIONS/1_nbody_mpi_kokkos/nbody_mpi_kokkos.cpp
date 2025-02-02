#include <Kokkos_Random.hpp>
#include <mpi.h>

#include <cstdlib>
#include <ctime>

#include "write_particles.h"

#define G 1
#define M 1000
#define PI 3.141592653589793238462643

void calculate_forces(View& particles_i, View& particles_j, int Np) {
  //using TeamPolicy = Kokkos::TeamPolicy<>;
  //using MemberType = TeamPolicy::member_type;

  //TeamPolicy policy(Np, Kokkos::AUTO()); // One team per particle

  //Kokkos::parallel_for("calculate_forces", policy, KOKKOS_LAMBDA(const MemberType& team) {
  //  const int i = team.league_rank(); // Each team handles one particle
  Kokkos::parallel_for("calculate_forces", Np, KOKKOS_LAMBDA(const int i) {
    double fx = 0.0, fy = 0.0, fz = 0.0;

    for (int i = 0; i < Np; i++)
    {
      if (i != j) {
        double dx = particles_i(i,1) - particles_j(j,1);
        double dy = particles_i(i,2) - particles_j(j,2);
        double dz = particles_i(i,3) - particles_j(j,3);
        
        double r = sqrt(dx*dx + dy*dy + dz*dz) + 1e-6;
        double Fmag = G * (particles_i(i,0) * particles_j(j,0)) / (r * r);

        fx += -Fmag * dx / r;
        fy += -Fmag * dy / r;
        fz += -Fmag * dz / r;
      }
    }

    /*
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, Np), [&](const int j, double& fx_sum, double& fy_sum, double& fz_sum) {
      if (i != j) {
        double dx = particles_i(i,1) - particles_j(j,1);
        double dy = particles_i(i,2) - particles_j(j,2);
        double dz = particles_i(i,3) - particles_j(j,3);
        
        double r = sqrt(dx*dx + dy*dy + dz*dz) + 1e-6;
        double Fmag = G * (particles_i(i,0) * particles_j(j,0)) / (r * r);

        fx_sum += -Fmag * dx / r;
        fy_sum += -Fmag * dy / r;
        fz_sum += -Fmag * dz / r;
      }
    }, fx, fy, fz);*/

    // Add central mass force
    double dx = particles_i(i,1);
    double dy = particles_i(i,2);
    double dz = particles_i(i,3);

    double r = sqrt(dx*dx + dy*dy + dz*dz) + 1e-6;
    double Fmag = G * (M * particles_i(i,0)) / (r * r);

    fx += -Fmag * dx / r;
    fy += -Fmag * dy / r;
    fz += -Fmag * dz / r;

    // Store results (outside reduction)
    particles_i(i, 7) = fx;
    particles_i(i, 8) = fy;
    particles_i(i, 9) = fz;
  });
}


void gather_forces_ranks(View& rank_particles, View& send_particles, View& recv_particles, int Np, int source, int dest, int size, MPI_Comm cart_comm, MPI_Request requests[2]) {

  // zero out forces
  Kokkos::parallel_for("zero_out_forces", Np, KOKKOS_LAMBDA(const int i) {
    rank_particles(i, 7) = 0.0;
    rank_particles(i, 8) = 0.0;
    rank_particles(i, 9) = 0.0;
  });

  // calculate force on own rank
  calculate_forces(rank_particles, rank_particles, Np);

  // Send our own particles
  Kokkos::deep_copy(send_particles, rank_particles);
  Kokkos::fence();
  
  for (int i = 0; i < size - 1; i++) {
    //std::cout << "Neigh Left: " << source << " Neigh Right: " << dest << std::endl; 
    MPI_Irecv(recv_particles.data(), Np*11, MPI_DOUBLE, source, 1, cart_comm, &requests[1]);
    MPI_Isend(send_particles.data(), Np*11, MPI_DOUBLE, dest, 1, cart_comm, &requests[0]);
    MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

    calculate_forces(rank_particles, recv_particles, Np);

    Kokkos::deep_copy(send_particles, recv_particles);
    Kokkos::fence();
  }
}

int main(int argc, char* argv[]){
  
  // Setup MPI environment
  MPI_Init(&argc, &argv);
  int rank, size, err, period, coords;

  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  err = MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Comm cart_comm;

  period = 1;    
  err = MPI_Cart_create(MPI_COMM_WORLD, 1, &size, &period, 0, &cart_comm);
  err = MPI_Cart_coords(cart_comm, rank, 1, &coords);

  int nbr_left, nbr_right;
  err = MPI_Cart_shift(cart_comm, 0, 1, &nbr_left, &nbr_right); 
  MPI_Request requests[2];
  // Setup MPI environment 

  // Initialize VTK controller 
  vtkNew<vtkMPIController> mpicontr;
  mpicontr->Initialize(&argc, &argv, 1); 
  // Initialize VTK controller

  // Initialize Kokkos
  Kokkos::initialize(argc, argv);
  {
    Kokkos::Random_XorShift64_Pool<> random_pool(/*seed=*/122394+rank);

    // Parse Command-Line arguments 
    int Np = atoi(argv[1]);       // Number of particles in simulation
    int Nt = atoi(argv[2]);       // Number of timesteps for simulation
    int N_write = atoi(argv[3]);  // Write frequency for output
    double dt = atof(argv[4]);    // timestep for simulation
    // Parse Command-Line arguments 

    // Initialize data structure to hold particles 
    View particles("particles", Np, 11);            // each particle (mass,x,y,z,vx,vy,vz,ax,ay,az)
    View send_particles("send_particles", Np, 11);  // each particle (mass,x,y,z,vx,vy,vz,ax,ay,az)
    View recv_particles("recv_particles", Np, 11);  // each particle (mass,x,y,z,vx,vy,vz,ax,ay,az)
    // Initialize data structure to hold particles 

    // Initialize the starting point for all particles
    Kokkos::parallel_for("initialize_particles", Np, KOKKOS_LAMBDA(const int i) {
      auto generator = random_pool.get_state();

      double r = generator.drand(10., 50.);
      double theta = 2*PI*generator.drand(0., 1.);
      double z = 10*generator.drand(0., 1.);

      particles(i, 0) = generator.drand(0.5, 1.5);

      particles(i, 1) = r*cos(theta); // x
      particles(i, 2) = r*sin(theta); // y
      particles(i, 3) = z; // z

      double eps = generator.drand(-0.01, 0.01);
      double vt = 1.5*sqrt(G*M/r)*(1+eps);

      // calculate tangential velocity (r_hat x k_hat)

      particles(i, 4) =  vt*particles(i, 2)/r;
      particles(i, 5) = -vt*particles(i, 1)/r;
      particles(i, 6) = 0.0;

      particles(i, 10) = (double) rank;

      random_pool.free_state(generator);
    });
    // Initialize the starting point for all particles

    // Simulation Loop
    for (int step = 0; step < Nt; step++) {
      // Visualization
      if (step % N_write == 0) {
        if (rank == 0) {
          std::cout << "Step: " << step << std::endl;
        }
        std::ostringstream filename;
        filename << "output_rank_" << step;
        write_vtk_polydata(filename.str(), particles, Np, rank, size, mpicontr);
      }
      // Visualization

      // Step update
      gather_forces_ranks(particles, send_particles, recv_particles, Np, nbr_left, nbr_right, size, cart_comm, requests);

      Kokkos::parallel_for("kick_leap_1", Np, KOKKOS_LAMBDA(const int i) {
        particles(i,4) += 0.5*dt*particles(i,7) / particles(i,0);
        particles(i,5) += 0.5*dt*particles(i,8) / particles(i,0);
        particles(i,6) += 0.5*dt*particles(i,9) / particles(i,0);

        particles(i,1) += particles(i,4)*dt;
        particles(i,2) += particles(i,5)*dt;
        particles(i,3) += particles(i,6)*dt;
      }); 

      gather_forces_ranks(particles, send_particles, recv_particles, Np, nbr_left, nbr_right, size, cart_comm, requests);

      Kokkos::parallel_for("kick_leap_2", Np, KOKKOS_LAMBDA(const int i) {
        particles(i,4) += 0.5*dt*particles(i,7) / particles(i,0);
        particles(i,5) += 0.5*dt*particles(i,8) / particles(i,0);
        particles(i,6) += 0.5*dt*particles(i,9) / particles(i,0);
      });
      // Step update
    }
    // Simulation Loop
  } 
  Kokkos::finalize();

  mpicontr->Finalize(1);

  MPI_Finalize();
  return 0;
}


