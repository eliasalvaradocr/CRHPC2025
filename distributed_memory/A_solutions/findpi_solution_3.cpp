#include "mpi.h"
#include <math.h>
#include <cstdlib>

#define PI25DT 3.141592653589793238462643

int main(int argc, char* argv[])
{
  int N, rank, size, i;
  double mypi, pi, h, sum, x;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  N = atoi(argv[1]);

  int rect_procs =  N / size;
  int extra_rects = N % size;

  int local_rects = rect_procs;  //  + (rank < extra_rects ? 1 : 0);

  if (rank < extra_rects ){
    local_rects += 1
  }

  int start = rank * rect_procs + (rank < extra_rects ? rank : extra_rects);

  h   = 1.0 / (double) N; // width of subdivision
  sum = 0.0;              // initialize sum to zero
  
  for (i = start; i <= start+local_rects; i++) {
    x = h * ((double)i - 0.5);  // midpoint calculation
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
