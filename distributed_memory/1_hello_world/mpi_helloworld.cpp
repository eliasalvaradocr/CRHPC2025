#include <iostream>
#include <mpi.h>

int main(int argc, char*argv[])
{
    int error;
    error = MPI_Init(&argc, &argv); // error-checking

    int rank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of the process


    char buffer[MPI_MAX_PROCESSOR_NAME];
    int resultlen = 0;

    MPI_Get_processor_name(buffer, &resultlen);

    std::cout << "Hello World! from: " << buffer << " rank: " << rank << std::endl;
    
    error = MPI_Finalize();
    
    return 0;
}