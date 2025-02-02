#include <iostream>

#include <Kokkos_Core.hpp>

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  int total_iterations = -1;
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <number of iterations>\n";
    std::exit(EXIT_FAILURE);
  } else {
    total_iterations = std::atoi(argv[1]);
  }

  double sum = 0.;
  Kokkos::parallel_reduce(
      total_iterations,
      KOKKOS_LAMBDA(const int i, double &partial_sum) {
        int sign = i % 2 == 0 ? 1 : -1;
        partial_sum += sign / (2. * i + 1.);
      },
      sum);
  std::cout << 4 * sum << '\n';
  Kokkos::finalize();
  return 0;
}
