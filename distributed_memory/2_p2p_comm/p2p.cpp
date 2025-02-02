#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv); // Initialize the MPI environment

    int rank, size, err;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of the process
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the number of processes

    if (size != 2) {
        if (rank == 0) {
            std::cout << "This program requires only two processes." << std::endl;
        }
        MPI_Finalize();
        return 1; // 
    }

    if (rank == 0) {
        // Process 0 sends a message
        int data = 42; // Data to send
        err = MPI_Send(&data, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        std::cout << "Process " << rank << " sent data " << data << " to process 1" << std::endl;
    } else if (rank == 1) {
        // Process 1 receives a message
        int received_data;
        MPI_Status status; // To check the details of the received message
        err = MPI_Recv(&received_data, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        std::cout << "Process " << rank << " received data " << received_data << " from process 0" << std::endl;
    }

    MPI_Finalize(); // Finalize the MPI environment
    return 0;
}
