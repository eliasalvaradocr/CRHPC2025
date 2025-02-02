#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv); // Initialize the MPI environment

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of the process
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the number of processes

    int data;

    if (rank == 0) {
        // Root process initializes the data
        data = 10;
        std::cout << "Process " << rank << " is broadcasting data: " << data << std::endl;
    }

    // Broadcast the value of `data` from process 0 to all processes
    MPI_Bcast(&data, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // All processes (including root) print the received data
    std::cout << "Process " << rank << " received data: " << data << std::endl;

    MPI_Finalize(); // Finalize the MPI environment
    return 0;
}
