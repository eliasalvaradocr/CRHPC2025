#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int *send_data;
    int recv_data;

    if (rank == 0) {
        send_data = new int[size];
        for (int i = 0; i < size; ++i) {
            send_data[i] = i * 10;
        }
    }
    // Documentation
    // int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //           void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
    //           MPI_Comm comm)
    
    MPI_Scatter(send_data, 1, MPI_INT, &recv_data, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::cout << "Process " << rank << " received data: " << recv_data << std::endl;

    MPI_Finalize();
    return 0;
}
