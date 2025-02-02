#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

void transpose_matrix_block(const int n, const int m, float **origin,
                            float **result) {
  for (int i = 0; i < n; i += BLOCK_SIZE) {
    for (int j = 0; j < n; j += BLOCK_SIZE) {
      for (int ii = i; ii < i + BLOCK_SIZE && ii < n; ++ii) {
        for (int jj = j; jj < j + BLOCK_SIZE && jj < m; ++jj) {
          result[jj][ii] = origin[ii][jj];
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  printf("Hello, world!\n");
  int n = atol(argv[1]);
  int m = atol(argv[2]);
  float **matA = alloc_matrix(n, m);
  float **matB = alloc_matrix(m, n);
  init_matrix(n, m, 2, matA);
  transpose_matrix_block(n, m, matA, matB);
  free_matrix(n, m, matA);
  free_matrix(m, n, matB);
  return 0;
}
