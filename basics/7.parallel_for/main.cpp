#include <Kokkos_Core_fwd.hpp>
#include <iostream>

#include <Kokkos_Core.hpp>

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  const int N = 10;

  const auto cpu_exec_policy =
      Kokkos::RangePolicy(Kokkos::DefaultHostExecutionSpace(), 0, N);

  Kokkos::View<double *, Kokkos::HostSpace> data_host("Host Data", N);
  Kokkos::parallel_for(
      "CPU execution", cpu_exec_policy,
      KOKKOS_LAMBDA(const int i) { data_host(i) = i; });

  Kokkos::View<double *, Kokkos::DefaultExecutionSpace> data_device(
      "Device Data", N);

  const auto dev_exec_policy =
      Kokkos::RangePolicy(Kokkos::DefaultExecutionSpace(), 0, N);

  Kokkos::parallel_for(
      "GPU execution", dev_exec_policy,
      KOKKOS_LAMBDA(const int i) { data_device(i) = N - i; });

  Kokkos::finalize();
}
