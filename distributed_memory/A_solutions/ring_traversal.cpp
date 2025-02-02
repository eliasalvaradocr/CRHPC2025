#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialize the value each process sends
    int send_value = rank;
    int recv_value = -1;

    // Determine neighbors in the ring
    int next = (rank + 1) % size;
    int prev = (rank - 1 + size) % size;

    for (int step = 0; step < size; step++) {
        // Send and receive values in the ring
        MPI_Sendrecv(&send_value, 1, MPI_INT, next, 0,
                     &recv_value, 1, MPI_INT, prev, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Update the value to be sent in the next step
        send_value = recv_value;

        // Print the current state of each process
        MPI_Barrier(MPI_COMM_WORLD);
        printf("Step %d: Process %d has value %d\n", step, rank, send_value);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // Final state: Ensure each process has its original value
    printf("Process %d ends with value %d\n", rank, send_value);

    MPI_Finalize();
    return 0;
}
