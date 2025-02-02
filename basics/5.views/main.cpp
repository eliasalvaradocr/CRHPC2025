#include <iostream>

#include <Kokkos_Core.hpp>

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  const int N = 10;
  {
    Kokkos::View<double *> data("data", N);

    Kokkos::parallel_for(
        "InitializeData", N, KOKKOS_LAMBDA(const int i) { data(i) = i + 1.0; });

    std::cout << "Initial data: ";
    Kokkos::parallel_for(
        "PrintData", N, KOKKOS_LAMBDA(const int i) {
          Kokkos::printf("i(%d,%.2f) ", i, data(i));
        });
    std::cout << '\n';

    Kokkos::parallel_for(
        "SquareData", N,
        KOKKOS_LAMBDA(const int i) { data(i) = data(i) * data(i); });

    std::cout << "Squared data: ";
    Kokkos::parallel_for(
        "PrintData", N, KOKKOS_LAMBDA(const int i) {
          Kokkos::printf("f(%d,%.2f) ", i, data(i));
        });

    std::cout << '\n';
  }

  Kokkos::finalize();
}
