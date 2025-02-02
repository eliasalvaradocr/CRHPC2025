#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int send_data = rank * 2;
    int *recv_data = new int[size];

    MPI_Allgather(&send_data, 1, MPI_INT, recv_data, 1, MPI_INT, MPI_COMM_WORLD);

    std::cout << "Process " << rank << " received: ";
    for (int i = 0; i < size; i++) {
        if (i < size - 1){std::cout << recv_data[i] << " ";} else {std::cout << recv_data[i] << " " <<std::endl;}
    }
    MPI_Finalize();
    return 0;
}
