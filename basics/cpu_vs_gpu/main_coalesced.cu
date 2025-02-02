#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

float *alloc_matrix_1D(const int n, const int m) {
  float *matrix = (float *)malloc(sizeof(float) * n * m);
  return matrix;
}

void print_matrix_1D(const int n, const int m, const float *matrix) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      printf("%.2f\t", matrix[i * m + j]);
    }
    printf("\n");
  }
}

__global__ void transpose_matrix_coalesced(const int n, const int m,
                                           const float *origin, float *result) {

  __shared__ float tile[16][16 + 1];
  const int i = blockIdx.y * 16 + threadIdx.y;
  const int j = blockIdx.x * 16 + threadIdx.x;

  if (i < n && j < m) {
    tile[threadIdx.y][threadIdx.x] = origin[j * m + i];
  }

  __syncthreads();

  if (i < n && j < m) {
    result[j * n + i] = tile[threadIdx.x][threadIdx.y];
  }
}

float *run_transpose_coalesced(const int n, const int m,
                               const float *host_origin) {
  float *dev_origin, *dev_result, *host_result;
  host_result = alloc_matrix_1D(m, n);
  cudaMalloc(&dev_origin, sizeof(float) * n * m);
  cudaMalloc(&dev_result, sizeof(float) * m * n);

  cudaMemcpy(dev_origin, host_origin, n * m * sizeof(float),
             cudaMemcpyHostToDevice);

  dim3 blockSize(16, 16);
  dim3 gridSize((n + blockSize.x - 1) / blockSize.x,
                (m + blockSize.y - 1) / blockSize.y);

  transpose_matrix_coalesced<<<gridSize, blockSize>>>(n, m, dev_origin,
                                                      dev_result);
  cudaDeviceSynchronize();
  cudaMemcpy(host_result, dev_result, m * n * sizeof(float),
             cudaMemcpyDeviceToHost);
  cudaFree(dev_origin);
  cudaFree(dev_result);
  return host_result;
}

void init_matrix_1D(const int n, const int m, const int base, float *matrix) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      matrix[i * m + j] = base * i * j + i + 2 * j;
    }
  }
}

void free_matrix_1D(float *matrix) { free(matrix); }

int main(int argc, char **argv) {
  printf("Hello, world!\n");
  int n = atol(argv[1]);
  int m = atol(argv[2]);
  float *matA = alloc_matrix_1D(n, m);
  float *matB;
  init_matrix_1D(n, m, 2, matA);
  matB = run_transpose_coalesced(n, m, matA);
  free_matrix_1D(matA);
  free_matrix_1D(matB);
  return 0;
}
