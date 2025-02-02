#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[])
{
    MPI_Init(&argv, &argc);

    // Get our rank and world size
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // we are working with an (Nrows, Ncols) matrix
    int Nrows = 6; 
    int Ncols = 6; 
    int rows_procs = Nrows / size;
    int extra_rows = Nrows % size;

    int local_rows = rows_procs + (rank < extra_rows ? 1 : 0); 

    // posistion of rown in global matrix
    int offset = rank * rows_procs + (rank < extra_rows ? rank : extra_rows);

    //std::cout << "Rank: " << rank << " Offset: " << offset << std::endl;

    // allocate local matrix A and x vector and b vector result
    int *local_A = new int[local_rows*Ncols];
    int *x = new int[Ncols];
    int *local_b = new int[local_rows];

    // Initialize x vector on root and broadcast to all processes
    if (rank == 0) {
        for (int i = 0; i < Ncols; i++) {
            x[i] = i + 1;
        }
    }    

    MPI_Bcast(x, Ncols, MPI_INT, 0, MPI_COMM_WORLD);

    // fill matrix and vectors with data
    for (int j = 0; j < local_rows; j++) {
        local_b[j] = 0;
        for (int i = 0; i < Ncols; i++) {
            local_A[j*Ncols+i] = rank;
        } 
    }

    // Uncomment to check local matrix A
    //printf("Rank %d local_A:\n", rank);
    //for (int j = 0; j < local_rows; j++) {
    //    for (int i = 0; i < Ncols; i++) {
    //        printf("%d ", local_A[i]);
    //    }
    //    printf("\n");
    //}

    // Uncomment to check vector x
    //printf("Rank %d x:\n", rank);
    //for (int i = 0; i < Ncols; i++) {
    //    printf("%d ", x[i]);
    //}

    // compute partial A*x
    for (int j = 0; j < local_rows; j++) {
        for (int i = 0; i < Ncols; i++){
            local_b[j] += local_A[j*Ncols+i]*x[i];
        }
    }

    // null variables for other procs
    int *recv_counts = NULL;
    int *displs = NULL;
    int *global_b = NULL; 

    // define variables for Vector Gather communication

    if (rank == 0){
        global_b = new int[Nrows];
        recv_counts = new int[size];
        displs = new int[size];
    
        for (int i = 0; i < size; i++) {
            recv_counts[i] = rows_procs + (i < extra_rows ? 1 : 0);
            displs[i] = i * rows_procs + (i < extra_rows ? i : extra_rows);
        }
    }

    MPI_Gatherv(local_b, local_rows, MPI_INT,
                global_b, recv_counts, displs, MPI_INT,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Global result:\n");
        for (int i = 0; i < Nrows; i++) {
            printf("%d ", global_b[i]);
        }
        printf("\n");

        // Cleanup
        delete[] global_b;
        delete[] recv_counts;
        delete[] displs;
    }


    // Cleanup local memory
    delete[] local_A;
    delete[] x;
    delete[] local_b;

    MPI_Finalize();
    return 0;
}