#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include <Kokkos_Core.hpp>

template <typename T> std::vector<T> read_vector(std::ifstream &file) {
  std::string line;
  std::getline(file, line);
  std::istringstream iss(line);
  T value;
  std::vector<T> vec;

  while (iss >> value) {
    vec.push_back(value);
  }
  return vec;
}

template <typename T>
void write_vector(const Kokkos::View<const T *> &vec, const int size,
                  const std::string &filename) {
  auto vec_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vec);
  std::ofstream outfile(filename);

  for (int i = 0; i < size; ++i) {
    outfile << vec_host(i);
    if (i < size - 1) {
      outfile << " ";
    }
  }
  outfile << '\n';
  outfile.close();
}

template <typename T> void print_vector(const std::vector<T> &vec) {
  for (int i = 0; i < vec.size(); ++i) {
    std::cout << vec[i] << ' ';
  }
  std::cout << '\n';
}

void spmv_csr(const Kokkos::View<const double *> &values,
              const Kokkos::View<const int *> &col_indices,
              const Kokkos::View<const int *> &row_ptr,
              const Kokkos::View<const double *> &x,
              const Kokkos::View<double *> &y) {
  const int n = row_ptr.size() - 1;
  Kokkos::parallel_for(
      "Spmv_csr", n, KOKKOS_LAMBDA(const int i) {
        double sum = 0.;
        for (int j = row_ptr(i); j < row_ptr(i + 1); ++j) {
          sum += values(j) * x(col_indices(j));
        }
        y(i) = sum;
      });
}

template <typename T> Kokkos::View<T *> copy_vector(const std::vector<T> &vec) {
  auto exec_policy =
      Kokkos::RangePolicy(Kokkos::DefaultHostExecutionSpace(), 0, vec.size());
  Kokkos::View<T *, Kokkos::HostSpace> host_view("Host View", vec.size());
  Kokkos::parallel_for(
      "Host View Copy", exec_policy,
      KOKKOS_LAMBDA(const int i) { host_view(i) = vec[i]; });
  auto dev_view = Kokkos::create_mirror_view_and_copy(
      Kokkos::DefaultExecutionSpace(), host_view);
  return dev_view;
}

int main(int argc, char *argv[]) {
  std::string filename;
  std::cout << "*************************************\n";
  std::cout << "Sparse Matrix x Vector Multiplication\n";
  std::cout << "*************************************\n";
  Kokkos::initialize(argc, argv);
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <filename>\n";
    exit(EXIT_FAILURE);
  } else {
    filename = argv[1];
  }
  std::ifstream file(filename);
  int num_rows;
  file >> num_rows;
  const auto max_size = std::numeric_limits<std::streamsize>::max();
  file.ignore(max_size, '\n');

  std::vector<double> values = read_vector<double>(file);
  std::vector<int> col_indices = read_vector<int>(file);
  std::vector<int> row_ptr = read_vector<int>(file);
  std::vector<double> x = read_vector<double>(file);
  file.close();

  {
    auto values_view = copy_vector(values);
    auto col_indices_view = copy_vector(col_indices);
    auto row_ptr_view = copy_vector(row_ptr);
    auto x_view = copy_vector(x);
    Kokkos::View<double *> y("Result", num_rows);
    spmv_csr(values_view, col_indices_view, row_ptr_view, x_view, y);
    write_vector<double>(y, num_rows, "result_kokkos.txt");
  }

  Kokkos::finalize();

  return 0;
}
