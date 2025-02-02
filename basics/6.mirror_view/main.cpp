#include <fstream>
#include <iostream>

#include <Kokkos_Core.hpp>

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  const int N = 10;
  {
    Kokkos::View<double *> data("data", N);

    Kokkos::parallel_for(
        "InitializeData", N, KOKKOS_LAMBDA(const int i) { data(i) = i + 1.0; });

    Kokkos::parallel_for(
        "SquareData", N,
        KOKKOS_LAMBDA(const int i) { data(i) = data(i) * data(i); });

    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), data);

    std::ofstream outfile("output.txt");

    for (int i = 0; i < N; ++i) {
      outfile << host_view(i);
      if (i < N - 1) {
        outfile << " ";
      }
    }

    outfile.close();
  }

  Kokkos::finalize();
}
