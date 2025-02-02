#include <iostream>

int main(int argc, char **argv) {
  int total_iterations = -1;
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <number of iterations>\n";
    std::exit(EXIT_FAILURE);
  } else {
    total_iterations = std::atoi(argv[1]);
  }

  double sum = 0.;
  for (int i = 0; i < total_iterations; ++i) {
    int sign = i % 2 == 0 ? 1 : -1;
    sum += sign / (2. * i + 1.);
  }
  std::cout << 4 * sum << '\n';
  return 0;
}
