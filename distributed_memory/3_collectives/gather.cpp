#include <mpi.h>
#include <iostream>


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv); // Initialize the MPI environment

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of the process
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the number of processes

    // Each process generates some data to send
    int send_data = rank * 10;

    int *recv_data;  // Allocate space for the received data

    // Root process will receive data from all processes
    if (rank == 0) {
        recv_data = new int[size];
    }

    // Gather data from all processes to the root process (rank 0)
    MPI_Gather(&send_data, 1, MPI_INT, recv_data, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Display results on the root process
    if (rank == 0) {
        std::cout << "Root process received the following data:" << std::endl;
        for (int i = 0; i < size; ++i) {
            std::cout << "recv_data[" << i << "] = " << recv_data[i] << std::endl;
        }
    }

    MPI_Finalize(); // Finalize the MPI environment
    return 0;
}
