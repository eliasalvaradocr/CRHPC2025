#include <iostream>

#include <Kokkos_Core.hpp>

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  const int N = 10;
  double *data;
  data = new double[N];

  for (int i = 0; i < N; ++i) {
    data[i] = i + 1.0;
  }

  std::cout << "Initial data: ";
  for (int i = 0; i < N; ++i) {
    std::cout << data[i] << ' ';
  }
  std::cout << '\n';

  Kokkos::parallel_for("SquareData", N,
                       [=](const int i) { data[i] = data[i] * data[i]; });

  std::cout << "Squared data: ";
  for (int i = 0; i < N; ++i) {
    std::cout << data[i] << ' ';
  }
  std::cout << '\n';

  Kokkos::finalize();
  delete[] data;
}
