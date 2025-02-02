#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

float **alloc_matrix(const int n, const int m) {
  float **matrix = malloc(sizeof(float *) * n);
  for (int i = 0; i < n; ++i) {
    matrix[i] = malloc(sizeof(float) * m);
  }
  return matrix;
}

void print_matrix(const int n, const int m, float **matrix) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      printf("%.2f\t", matrix[i][j]);
    }
    printf("\n");
  }
}

void free_matrix(const int n, const int m, float **matrix) {
  for (int i = 0; i < n; ++i) {
    free(matrix[i]);
  }
  free(matrix);
}

void init_matrix(const int n, const int m, int base, float **matrix) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      matrix[i][j] = base * i * j + i + 2 * j;
      ++base;
    }
  }
}
