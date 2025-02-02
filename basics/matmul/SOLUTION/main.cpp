#include <Kokkos_Core.hpp>
#include <cstdlib>
#include <iostream>

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  const int width = 2;
  const int height = 2;
  int *matA = (int *)malloc(sizeof(int) * width * height);
  int *matB = (int *)malloc(sizeof(int) * width * height);
  int *matC = (int *)malloc(sizeof(int) * width * height);

  Kokkos::parallel_for(
      width * height, KOKKOS_LAMBDA(int i) {
        matA[i] = i;
        matB[i] = 10 - i;
        matC[i] = 0;
      });

  Kokkos::parallel_for(
      width, KOKKOS_LAMBDA(int i) {
        for (int j = 0; j < height; ++j) {
          for (int k = 0; k < height; ++k) {
            matC[i * height + j] += matA[i * height + k] * matB[k * height + j];
          }
        }
      });

  for (int i = 0; i < width; ++i) {
    for (int j = 0; j < height; ++j) {
      std::cout << matC[i * height + j] << ' ';
    }
    std::cout << '\n';
  }

  free(matA);
  free(matB);
  free(matC);
  Kokkos::finalize();
}
