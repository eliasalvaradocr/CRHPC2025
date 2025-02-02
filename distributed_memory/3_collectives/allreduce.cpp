#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int send_data = rank + 1;
    int recv_data;

    MPI_Allreduce(&send_data, &recv_data, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);

    std::cout << "Process " << rank << " computed product: " << recv_data << std::endl;

    MPI_Finalize();
    return 0;
}
