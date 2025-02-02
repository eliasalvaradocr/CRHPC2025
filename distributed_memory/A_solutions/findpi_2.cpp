#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int total_points = ;
    int points_per_process = total_points / size;
    int local_count = 0;

    // Seed for random number generation
    srand(rank + 42);

    // Monte Carlo simulation
    for (int i = 0; i < points_per_process; i++) {
        double x = (double)rand() / RAND_MAX * 2.0 - 1.0; // [-1, 1]
        double y = (double)rand() / RAND_MAX * 2.0 - 1.0; // [-1, 1]
        if (x * x + y * y <= 1.0) {
            local_count++;
        }
    }

    // Reduce to get the total count of points inside the circle
    int global_count = 0;
    MPI_Reduce(&local_count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double pi_estimate = 4.0 * global_count / total_points;
        std::cout << "Estimated value of π: " << pi_estimate << std::endl;
    }

    MPI_Finalize();
    return 0;
}
