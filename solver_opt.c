/*
 * Tema 2 ASC
 * 2019 Spring
 * Catalin Olaru / Vlad Spoiala
 */
#include "utils.h"
#include <string.h>
/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	double *B_trans;
	double *A_trans;
	int i, k, j;
	double *result_first;
	double *result_second;
	double *result_final;
	double sum1, sum2, sum3;

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

	for (i = 0; i <  N; i ++) {
		for (k = 0; k < N; k++)  {
			sum1 = A_trans[i * N + k];
			sum2 = B_trans[i * N + k];
			for (j = 0; j < N ; j ++) {
				result_first[i * N + j] += sum1 * B[k * N + j];
				result_second[i * N + j] += sum2 * A[k * N + j];
			}
		}
	}
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

	memset(result_first, 0, N * N * sizeof(double));

	for (i = 0; i <  N; i ++) {
		for (k = 0; k < N; k++) {
			sum3 = result_final[i * N + k];
			for(j = 0; j < N; j++) {
				result_first[i * N + j] += sum3 * result_final[k * N + j];
			}
		}
	}
	
	printf("OPT SOLVER\n");
	free(A_trans);
	free(B_trans);
	free(result_second);
	free(result_final);

	return result_first;
}
