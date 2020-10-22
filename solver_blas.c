/*
 * Tema 2 ASC
 * 2019 Spring
 * Catalin Olaru / Vlad Spoiala
 */
#include "utils.h"
#include "cblas.h"

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {
	double *B_trans;
	double *A_trans;
	int i, j;
	double *result_first;
	double *result_second;
	double *result_final;

	result_first = (double*) malloc(N * N * sizeof(double));
	result_second = (double*) malloc(N * N * sizeof(double));
	result_final = (double*) malloc(N * N * sizeof(double));
	B_trans = (double*) malloc(N * N * sizeof(double));
	A_trans = (double*) malloc(N * N * sizeof(double));

	for (i = 0; i <  N; i ++) {
		for (j = 0; j < N; j++)  {
			A_trans[i * N + j] = A[j * N + i];
			B_trans[i * N + j] = B[j * N + i];
		}
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					N, N, N, 1.0, A_trans, N, B, N, 0.0, result_first, N);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					N, N, N, 1.0, B_trans, N, A, N, 0.0, result_second, N);

	for (i = 0; i <  N; i ++) {
		for (j = 0; j < N; j++)  {
			result_final[i * N + j] = result_first[i * N + j] + result_second[i * N + j];
		}
	}

	for (i = 0; i <  N; i ++) {
		for (j = 0; j < N; j++) {
			if (i > j)
			result_final[i * N + j] = 0.0;
		}
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					N, N, N, 1.0, result_final, N, result_final, N, 0.0, result_first, N);

	printf("BLAS SOLVER\n");
	free(A_trans);
	free(B_trans);
	free(result_second);
	free(result_final);

	return result_first;
}
