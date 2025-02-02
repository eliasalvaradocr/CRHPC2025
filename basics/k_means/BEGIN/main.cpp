#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

std::vector<std::pair<double, double>> read_csv(const std::string &filename) {
  std::vector<std::pair<double, double>> points;
  std::ifstream file(filename);
  std::string line;

  while (std::getline(file, line)) {
    double x, y;
    if (sscanf(line.c_str(), "%lf,%lf", &x, &y) == 2) {
      points.emplace_back(x, y);
    }
  }

  return points;
}

constexpr double distance_squared(const double x1, const double y1,
                                  const double x2, const double y2) {
  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

void write_csv(const std::string &filename_base,
               const std::vector<std::pair<double, double>> points,
               const std::vector<int> &clusters, const int iteration) {
  const std::string filename =
      filename_base + std::to_string(iteration) + ".csv";

  std::ofstream file(filename);
  if (file.is_open()) {
    file << "x,y,cluster\n";
    for (int i = 0; i < points.size(); ++i) {
      file << points[i].first << "," << points[i].second << "," << clusters[i]
           << "\n";
    }
    file.close();
  } else {
    std::cerr << "Unable to open file for writing\n";
  }
}

void k_means(const std::vector<std::pair<double, double>> &points, const int k,
             const int max_iters, std::vector<int> clusters) {
  const size_t num_points = points.size();

  std::vector<std::pair<double, double>> centroids(k);
  std::cout << "Initializing centroids...\n";
  for (int centr = 0; centr < k; ++centr) {
    centroids[centr].first = points[rand() % num_points].first;
    centroids[centr].second = points[rand() % num_points].second;
  }

  std::vector<std::pair<double, double>> new_centroids(k);
  std::vector<int> counts(k);

  std::cout << "Starting iterations...\n";
  for (int iter = 0; iter < max_iters; ++iter) {
    for (int pt_idx = 0; pt_idx < num_points; ++pt_idx) {
      double min_distance = std::numeric_limits<double>::max();
      int best_cluster = 0;

      for (int clus = 0; clus < k; ++clus) {
        const double dist =
            distance_squared(points[pt_idx].first, points[pt_idx].second,
                             centroids[clus].first, centroids[clus].second);
        if (dist < min_distance) {
          min_distance = dist;
          best_cluster = clus;
        }
      }
      clusters[pt_idx] = best_cluster;
    }

    for (int clus = 0; clus < k; ++clus) {
      new_centroids[clus].first = 0.;
      new_centroids[clus].second = 0.;
      counts[clus] = 0;
    }

    for (int pt_idx = 0; pt_idx < num_points; ++pt_idx) {
      int cluster = clusters[pt_idx];
      new_centroids[cluster].first += points[pt_idx].first;
      new_centroids[cluster].second += points[pt_idx].second;
      counts[cluster] += 1;
    }

    for (int clus = 0; clus < k; ++clus) {
      if (counts[clus] > 0) {
        centroids[clus].first = new_centroids[clus].first / counts[clus];
        centroids[clus].second = new_centroids[clus].second / counts[clus];
      }
    }
    // std::cout << "Writing file...\n";
    write_csv("clusters", points, clusters, iter);
  }
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input_csv> <k> <max_iterations>\n";
    exit(EXIT_FAILURE);
  }
  srand(time(NULL));
  std::cout << "******************\n";
  std::cout << "K-means clustering\n";
  std::cout << "******************\n";

  const std::string input_csv = argv[1];
  const int k = std::stoi(argv[2]);
  const int max_iters = std::stoi(argv[3]);

  std::cout << "Reading file...\n";
  auto points = read_csv(input_csv);

  std::cout << "Clustering data, " << points.size() << "...\n";
  std::vector<int> clusters(points.size());
  k_means(points, k, max_iters, clusters);
}
