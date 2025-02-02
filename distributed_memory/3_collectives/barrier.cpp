#include <mpi.h>
#include <iostream>
#include <chrono>
#include <thread>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::this_thread::sleep_for(std::chrono::seconds(2));
    }
    else { 
        std::cout << "Process " << rank << " keeps going" << std::endl;
    }
    std::cout << "Process " << rank << " reached barrier" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << "Process " << rank << " passed the barrier" << std::endl;

    MPI_Finalize();
    return 0;
}
