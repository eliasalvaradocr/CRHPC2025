#ifndef MATRIX_H
#define MATRIX_H
#define BLOCK_SIZE 32

#ifdef __cplusplus
extern "C" {
#endif
float **alloc_matrix(const int n, const int m);
void print_matrix(const int n, const int m, float **matrix);
void free_matrix(const int n, const int m, float **matrix);
void init_matrix(const int n, const int m, int base, float **matrix);
#ifdef __cplusplus
}
#endif

#endif
