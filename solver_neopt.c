/*
 * Tema 2 ASC
 * 2019 Spring
 * Catalin Olaru / Vlad Spoiala
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	double *B_trans;
	double *A_trans;
	int i, k, j;
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

	for (i = 0; i <  N; i ++) {
		for (j = 0; j < N; j++)  {
			result_first[i * N + j] = 0.0;
			result_second[i * N + j] = 0.0;
			for (k = 0; k < N ; k ++) {
				result_first[i * N + j] += A_trans[i * N + k] * B[k * N + j];
				result_second[i * N + j] += B_trans[i * N + k] * A[k * N + j];
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

	for (i = 0; i <  N; i ++) {
		for (j = 0; j < N; j++) {
			result_first[i * N + j] = 0.0;
			for(k = 0; k < N; k++) {
				result_first[i * N + j] += result_final[i * N + k] * result_final[k * N + j];
			}
		}
	}
	
	printf("NEOPT SOLVER\n");
	free(A_trans);
	free(B_trans);
	free(result_second);
	free(result_final);

	return result_first;
}
