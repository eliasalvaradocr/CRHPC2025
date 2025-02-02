#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int send_data = rank;
    int recv_data;
    int next = (rank + 1) % size;
    int prev = (rank - 1 + size) % size;

    // Ring exchange
    MPI_Sendrecv(&send_data, 1, MPI_INT, next, 0,
                 &recv_data, 1, MPI_INT, prev, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    std::cout << "Rank " << rank << " sent " << send_data
              << ", received " << recv_data << std::endl;

    MPI_Finalize();
    return 0;
}
