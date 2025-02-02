#include "mpi.h"
#include <math.h>
#include <cstdlib>

#define PI25DT 3.141592653589793238462643

int main(int argc, char* argv[])
{
  int n, rank, size;
  double mypi, pi, h, sum, x;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  n = atoi(argv[1]);

  // we need to have the number of intervals in all our ranks
  // MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  h   = 1.0 / (double) n; // width of subdivision
  sum = 0.0;              // initialize sum to zero
  double start = (n/size) * rank; // 0, 1, 2, ... 
  double stop = (n/size) * (rank + 1); // 1, 2, 3, ... 
  
  for (int i = start; i <= stop; i += 1) {
    x = h * ((double) i - 0.5);  // midpoint calculation
    sum += 4.0 / (1.0 + x*x);
  }

  mypi = h * sum;

  // reduce and communicate result to rank 
  MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (rank == 0)
      printf("pi is approximately %.16f, Error is %.16f\n",
      pi, fabs(pi - PI25DT));

  MPI_Finalize();
  return 0;
}
