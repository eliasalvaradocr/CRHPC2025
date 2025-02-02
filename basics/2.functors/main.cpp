#include <iostream>

#include <Kokkos_Core.hpp>

struct SquareFunctor {
  double *data;

  SquareFunctor(double *data) : data(data) {}

  void operator()(const int i) const { data[i] = data[i] * data[i]; }
};

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  const int N = 10;
  double data[N];

  for (int i = 0; i < N; ++i) {
    data[i] = i + 1.0;
  }

  std::cout << "Initial data: ";
  for (int i = 0; i < N; ++i) {
    std::cout << data[i] << ' ';
  }
  std::cout << '\n';

  Kokkos::parallel_for("SquareData", N, SquareFunctor(data));

  std::cout << "Squared data: ";
  for (int i = 0; i < N; ++i) {
    std::cout << data[i] << ' ';
  }
  std::cout << '\n';

  Kokkos::finalize();
}
