#include <iostream>

#include <Kokkos_Core.hpp>

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  std::cout << "******************\n";
  std::cout << "Hello, world!\n";
  for (int i = 0; i < argc; ++i) {
    std::cout << i << ": " << argv[i] << '\n';
  }
  Kokkos::finalize();
}
