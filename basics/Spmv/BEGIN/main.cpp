#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

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
void write_vector(const std::vector<T> &vec, const std::string &filename) {
  std::ofstream outfile(filename);

  for (int i = 0; i < vec.size(); ++i) {
    outfile << vec[i];
    if (i < vec.size() - 1) {
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

void spmv_csr(const std::vector<double> &values,
              const std::vector<int> &col_indices,
              const std::vector<int> &row_ptr, const std::vector<double> &x,
              std::vector<double> &y) {
  const int n = row_ptr.size() - 1;
  for (int i = 0; i < n; ++i) {
    double sum = 0.;
    for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
      sum += values[j] * x[col_indices[j]];
    }
    y[i] = sum;
  }
}

int main(int argc, char *argv[]) {
  std::string filename;
  std::cout << "*************************************\n";
  std::cout << "Sparse Matrix x Vector Multiplication\n";
  std::cout << "*************************************\n";
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

  std::vector<double> y(num_rows);

  spmv_csr(values, col_indices, row_ptr, x, y);
  write_vector(y, "result_serial.txt");

  return 0;
}
