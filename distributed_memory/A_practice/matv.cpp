#include <iostream>
#include <chrono>

int main(int argc, char **argv){

  int Nrows; 
  int Ncols;

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
            << " <Nrows> <Ncols>\n";
    exit(-1);
  } else {
    Nrows = atoi(argv[1]);
    Ncols = atoi(argv[2]);
  }

  int *A = new int[Nrows*Ncols];
  int *x = new int[Ncols];
  int *b = new int[Ncols];

  for (int i = 0; i < Ncols; i++) {
      x[i] = i + 1;
  }    

  for (int j = 0; j < Nrows; j++) {
      b[j] = 0;
      for (int i = 0; i < Ncols; i++) {
          A[j*Ncols+i] = 1;
      } 
  }

  // compute A*x
  std::chrono::duration<double> total_time;
  auto start{std::chrono::steady_clock::now()};
  for (int j = 0; j < Nrows; j++) {
      for (int i = 0; i < Ncols; i++){
          b[j] += A[j*Ncols+i]*x[i];
      }
  }
  total_time += std::chrono::steady_clock::now() - start;

  printf("Total time: %.8f seconds\n", total_time.count());    
  //printf("Result:\n");
  //for (int i = 0; i < Nrows; i++) {
  //    printf("%d ", b[i]);
  //}
  //printf("\n");
  return 0; 
}