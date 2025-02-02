#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

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
  file.close();
  return points;
}

KOKKOS_INLINE_FUNCTION
constexpr double distance_squared(const double x1, const double y1,
                                  const double x2, const double y2) {
  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

void write_csv(const std::string &filename_base,
               const Kokkos::View<const double *[2]> &points,
               const int total_points,
               const Kokkos::View<const int *> &clusters, const int k,
               const int iteration) {
  const std::string filename =
      filename_base + std::to_string(iteration) + ".csv";

  std::ofstream file(filename);
  auto points_host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), points);
  auto clusters_host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), clusters);
  if (file.is_open()) {
    file << "x,y,cluster\n";
    for (int i = 0; i < total_points; ++i) {
      file << points_host(i, 0) << "," << points_host(i, 1) << ","
           << clusters_host(i) << "\n";
    }
    file.close();
  } else {
    std::cerr << "Unable to open file for writing\n";
  }
}

void k_means(const Kokkos::View<const double *[2]> &points,
             const int num_points, const int k, const int max_iters,
             Kokkos::View<int *> &clusters) {

  Kokkos::View<double *[2]> centroids("centroids", k);
  std::cout << "Initializing centroids...\n";
  int centr_stride = num_points / k;
  Kokkos::parallel_for(
      "Initialize centroids", k, KOKKOS_LAMBDA(const int centr) {
        centroids(centr, 0) = points(centr_stride * centr, 0);
        centroids(centr, 1) = points(centr_stride * centr, 1);
      });

  Kokkos::View<double *[2]> new_centroids("new_centroids", k);
  Kokkos::View<int *> counts("counts", k);

  std::cout << "Starting iterations...\n";

  for (int iter = 0; iter < max_iters; ++iter) {
    Kokkos::parallel_for(
        "Iterate points", num_points, KOKKOS_LAMBDA(const int pt_idx) {
          // double min_distance = std::numeric_limits<double>::max();
          double min_distance = 1.7e308;
          int best_cluster = 0;

          for (int clus = 0; clus < k; ++clus) {
            const double dist =
                distance_squared(points(pt_idx, 0), points(pt_idx, 1),
                                 centroids(clus, 0), centroids(clus, 1));
            if (dist < min_distance) {
              min_distance = dist;
              best_cluster = clus;
            }
          }
          clusters(pt_idx) = best_cluster;
        });
    Kokkos::parallel_for(
        "Reset centroids", k, KOKKOS_LAMBDA(const int clus) {
          new_centroids(clus, 0) = 0.;
          new_centroids(clus, 1) = 0.;
          counts(clus) = 0;
        });
    Kokkos::parallel_for(
        "Assign centroids", num_points, KOKKOS_LAMBDA(const int pt_idx) {
          const int cluster = clusters(pt_idx);
          Kokkos::atomic_add(&new_centroids(cluster, 0), points(pt_idx, 0));
          Kokkos::atomic_add(&new_centroids(cluster, 1), points(pt_idx, 1));
          Kokkos::atomic_add(&counts(cluster), 1);
        });
    Kokkos::parallel_for(
        "Finalize centroids", k, KOKKOS_LAMBDA(const int clus) {
          if (counts(clus) > 0) {
            centroids(clus, 0) = new_centroids(clus, 0) / counts(clus);
            centroids(clus, 1) = new_centroids(clus, 1) / counts(clus);
          }
        });
    write_csv("clusters", points, num_points, clusters, k, iter);
  }
}

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input_csv> <k> <max_iterations>\n";
    Kokkos::finalize();
    exit(EXIT_FAILURE);
  }
  srand(time(NULL));
  std::cout << "******************\n";
  std::cout << "K-means clustering\n";
  std::cout << "******************\n";

  const std::string input_csv = argv[1];
  const int k = std::stoi(argv[2]);
  const int max_iters = std::stoi(argv[3]);
  {
    std::cout << "Reading file...\n";
    auto points = read_csv(input_csv);
    const int total_points = points.size();

    Kokkos::View<double *[2]> points_view("Points", total_points);
    auto points_view_host = Kokkos::create_mirror_view(points_view);

    for (int pt_idx = 0; pt_idx < total_points; ++pt_idx) {
      points_view_host(pt_idx, 0) = points[pt_idx].first;
      points_view_host(pt_idx, 1) = points[pt_idx].second;
    }

    Kokkos::deep_copy(points_view, points_view_host);

    std::cout << "Clustering data...\n";
    Kokkos::View<int *> clusters("clusters", total_points);
    k_means(points_view, total_points, k, max_iters, clusters);
  }

  Kokkos::finalize();
}
