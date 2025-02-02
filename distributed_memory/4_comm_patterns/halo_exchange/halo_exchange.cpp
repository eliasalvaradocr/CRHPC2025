#include <mpi.h>
#include <iostream>
#include <sstream>
#include "write_halo.h"

int main(int argc, char *argv[]) {
    // Initialize MPI communicator
    MPI_Init(&argc, &argv);

    // Initialize VTK MPI handler
    vtkNew<vtkMPIController> mpicontr;
    mpicontr->Initialize(&argc, &argv, 1); 
    // 

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define the 3D grid dimensions
    int dims[3] = {2, 2, 2};    // 2x2x2 Cartesian grid
    int periods[3] = {1, 1, 1}; // Non-periodic boundaries
    int coords[3];              // Rank's coordinates in the grid
    
    // Define cartesian communicator for cartesian topology
    MPI_Comm cart_comm;

    // Create Cartesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, coords);
    // Neighbors in each direction

    int nbr_left, nbr_right, nbr_up, nbr_down, nbr_front, nbr_back;
    MPI_Cart_shift(cart_comm, 0, 1, &nbr_left, &nbr_right);  // x-direction
    MPI_Cart_shift(cart_comm, 1, 1, &nbr_down, &nbr_up);     // y-direction
    MPI_Cart_shift(cart_comm, 2, 1, &nbr_back, &nbr_front);  // z-direction

    const int N = 4;                        // Local size of the grid in each dimension

    // Calculate local extents
    //int N {10}, N{10}, N{10};        //local dimensions of the process's grid
    int physdims[3] = {N+1, N+1, N+1};
    int global_extent[6] = {0, dims[0]*N, 
                            0, dims[1]*N,
                            0, dims[2]*N};

    int local_extent[6] = {coords[0]*N, (coords[0]+1)*N,
                           coords[1]*N, (coords[1]+1)*N,
                           coords[2]*N, (coords[2]+1)*N};

    

    // Local 3D grid for each process
    int *grid = new int[N * N * N];         // Allocate 3D grid dynamicalN
    int *grid_x = new int[N * N * N];       // Allocate 3D grid for exchange
    int *grid_y = new int[N * N * N];       // Allocate 3D grid for exchange
    int *grid_z = new int[N * N * N];       // Allocate 3D grid for exchange


    for (int i = 0; i < N * N * N; ++i) {
        grid[i] = rank;                     // Initialize local grid with the rank
        grid_x[i] = rank;                     // Initialize local grid with the rank
        grid_y[i] = rank;                     // Initialize local grid with the rank
        grid_z[i] = rank;                     // Initialize local grid with the rank
    }

    std::ostringstream filename_1;
    filename_1 << "output_rank_1";
    
    write_vtk(filename_1.str(), grid, N, rank, size, 
      coords, local_extent, global_extent, physdims, mpicontr);

    // Buffers for halo exchange
    int *send_x = new int[N * N];
    int *recv_x = new int[N * N];
    int *send_y = new int[N * N];
    int *recv_y = new int[N * N];
    int *send_z = new int[N * N];
    int *recv_z = new int[N * N];

    // Prepare send buffers (example: sending first layer of each dimension)


    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < N; j++)
      {
        send_x[j*N+i] = grid[j*N*N + i*N];
        send_y[j*N+i] = grid[j*N*N+i];   
        send_z[j*N+i] = grid[j*N+i];   

      }
    }

 
     // Exchange in x-direction
    MPI_Sendrecv(send_x, N * N, MPI_INT, nbr_left, 0,
                 recv_x, N * N, MPI_INT, nbr_right, 0,
                 cart_comm, MPI_STATUS_IGNORE);

    // Exchange in y-direction
    MPI_Sendrecv(send_y, N * N, MPI_INT, nbr_down, 1,
                 recv_y, N * N, MPI_INT, nbr_up, 1,
                 cart_comm, MPI_STATUS_IGNORE);

    // Exchange in z-direction
    MPI_Sendrecv(send_z, N * N, MPI_INT, nbr_back, 2,
                 recv_z, N * N, MPI_INT, nbr_front, 2,
                 cart_comm, MPI_STATUS_IGNORE);

    // Overwrite received data
    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < N; j++)
      {
        grid_x[j*N*N + i*N + (N-1)] = recv_x[j*N+i];
        grid_y[j*N*N + (N-1)*N + i] = recv_y[j*N+i];   
        grid_z[(N-1)*N*N + j*N + i] = recv_z[j*N+i];
      }
    }

    std::cout << "Rank " << rank << " completed exchanges." << std::endl;

    // Print results

    std::ostringstream filename_2_x;
    filename_2_x << "output_rank_2_x";
    
    std::ostringstream filename_2_y;
    filename_2_y << "output_rank_2_y";
    
    std::ostringstream filename_2_z;
    filename_2_z << "output_rank_2_z";

    write_vtk(filename_2_x.str(), grid_x, N, rank, size, 
      coords, local_extent, global_extent, physdims, mpicontr);

    write_vtk(filename_2_y.str(), grid_y, N, rank, size, 
      coords, local_extent, global_extent, physdims, mpicontr);
    
    write_vtk(filename_2_z.str(), grid_z, N, rank, size, 
      coords, local_extent, global_extent, physdims, mpicontr);

    // Clean up dynamicalN allocated memory
    delete[] grid;
    delete[] grid_x;
    delete[] grid_y;
    delete[] grid_z;
    delete[] send_x;
    delete[] recv_x;
    delete[] send_y;
    delete[] recv_y;
    delete[] send_z;
    delete[] recv_z;

    // Finalize VTK processes
    mpicontr->Finalize(1);

    // Finalize MPI processes
    MPI_Finalize();
    return 0;
}
