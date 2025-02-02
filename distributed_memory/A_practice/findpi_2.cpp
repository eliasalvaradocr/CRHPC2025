#include <chrono>
#include <iostream>
#include <cstdio>
#include <math.h>

#define PI25DT 3.141592653589793238462643


int main(int argc, char** argv) {
  int N;

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
            << " <N> \n";
    exit(-1);
  } else {
    N = std::atoi(argv[1]);
  }

  int local_count = 0;
    
  // Seed for random number generation
  srand(42);

  std::chrono::duration<double> total_time;
  auto start{std::chrono::steady_clock::now()};
  for (int i = 0; i < N; i++){
    double x = (double)rand() / RAND_MAX * 2.0 - 1.0; // [-1, 1]
    double y = (double)rand() / RAND_MAX * 2.0 - 1.0; // [-1, 1]
    if (x * x + y * y <= 1.0) {
        local_count++;
    }
  }
    double pi = 4.0 * local_count / N;
    
    total_time += std::chrono::steady_clock::now() - start;
    
    printf("Total time: %.5f seconds\n", total_time.count());    
    printf("pi is approximately %.16f, Error is %.16f\n",
      pi, fabs(pi - PI25DT));
}