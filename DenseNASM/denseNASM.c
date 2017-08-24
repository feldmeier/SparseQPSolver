#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../LAPACK/LAPACKE/include/lapacke.h"
#include "../LAPACK/CBLAS/include/cblas.h"

#include <time.h>
#include <stdlib.h>

#define PARSER_DEBUG 0

#include "parser.h"

#ifndef INF
#define INF 10E18
#endif

#ifndef EPS
#define EPS 10E-8
#endif

#define SOLVER_DEBUG 1

void copyVector(double* v, double* t, int elements) {
	for(int i = 0; i < elements; i++) {
		t[i] = v[i];
	}
}

void selectRow(double* m, int idx, double* target, int rows, int cols) {
	for(int i = 0; i < cols; i++) {
		target[i]  = m[i*rows+idx];
	}
}

void selectColumn(double* m, int idx, double* target, int rows, int cols) {
	for(int i = 0; i < rows; i++) {
		target[i] = m[idx*cols+i];
	}
}


void override(double* v, int N, double value) {
	for(int i = 0; i < N; i++)
		v[i] = value;
}

void ones(double* v, int N) {
	override(v, N, 1);
}

void zeros(double* v, int N) {
	override(v, N, 0);
}

//can be avoided by multiplication with 1-vector, still chose that way cause of datatypes (booleans as doubles are just nonsense)
int nnz(char* v, int N) {
	int res = 0;
	for(int i = 0; i < N; i++)
		res += v[i] != 0;
	return res;
}

int nnz_d(double* v, int N) {
	int res = 0;
	for(int i = 0; i < N; i++) {

		res += v[i] >= EPS || v[i] <= -EPS;
	}
	return res;
}
void mat2vec(double** mat, double* v, int rows, int cols, char byRow) {
	if(byRow) {
		int idx = 0;
		for(int i = 0; i < rows; i++) {
			copyVector(mat[i], &(v[idx]), cols);
			idx += cols;
		}		
	} else {
		int idx = 0;
		for(int i = 0; i < cols; i++) {
			for(int j = 0; j < rows; j++) {
				v[idx++] = mat[j][i];
			}
		}
	}
}

double min(double a, double b){
	return a < b ? a : b;
}


void blas_copy(double* x, double* y, int N) {
	int incx = 1;
	int incy = 1;
	dcopy_(&N, x, &incx, y, &incy);
}

void blas_axpy(double* x, double* y, double* tar, double alpha, int N) {
	blas_copy(y, tar, N);
	int incx = 1;
	int incy = 1;
	daxpy_(&N, &alpha, x, &incx, tar, &incy);
}

double blas_dot(double* x, double* y, int N) {
	int incx=1;
	int incy=1;
	//return ddot_(&N, x, &incx, y, &incy);
	double res = 0;
	for(int i = 0; i < N; i++)
		res += x[i] * y[i];
	return res;
}

void print_vector(void* vec, int N, char type) {
	double* dv = (double*)vec;
	char* cv = (char*)vec;
	for(int i = 0; i < N; i++) {
		if(type == 'd')
			printf("%lf ", dv[i]);
		if(type == 'c')
			printf("%d ", cv[i]);
	}
	printf("\n");
}

void print_matrix(double** mat, int N_ROWS, int N_COLS) {
	for(int i = 0; i < N_ROWS; i++) {
		print_vector(mat[i], N_COLS, 'd');
	}
}

//does not work probably
//matrix A: M rows, k cols
//matris B: K rows, N Cols
//matrix C: K rows, k cols
void matmatprod(double* A_arr, double* B_arr, double* target_arr, char TRANA, char TRANB, int M, int K, int N) {
	double alpha = 1;
	double beta = 0;
	int lda = N;
	int ldb = K;
	int ldc = N;
	dgemm_(&TRANA, &TRANB, &M, &N, &K, &alpha, A_arr, &lda, B_arr, &ldb, &beta, target_arr, &ldc);
}

void vecmatprod(double* A_arr, double* x, double* target,  int N_ROWS, int N_COLS) {
	char TRANS = 'T';

	double alpha = 1;
	int lda = N_ROWS;
	int incx = 1;
	double beta = 0;
	int incy = 1;
	dgemv_(&TRANS, &N_COLS, &N_ROWS, &alpha, A_arr, &lda, x, &incx, &beta, target, &incy);
}


void matvecprod(double* A_arr, double* x, double* target,  int N_ROWS, int N_COLS) {
	char TRANS = 'N';

	double alpha = 1;
	int lda = N_ROWS;
	int incx = 1;
	double beta = 0;
	int incy = 1;
	dgemv_(&TRANS, &N_ROWS, &N_COLS, &alpha, A_arr, &lda, x, &incx, &beta, target, &incy);
}

void matvecprodadd(double* A_arr, double* x, double* y, double* target, int N_ROWS, int N_COLS) {
	blas_copy(y, target, N_ROWS);
	char TRANS = 'N';
	double alpha = 1;
	int lda = N_ROWS;
	int incx = 1;
	double beta = 1;
	int incy = 1;
	dgemv_(&TRANS, &N_ROWS, &N_COLS, &alpha, A_arr, &lda, x, &incx, &beta, target, &incy);
}

void solve(double* A_arr, double * b, double* target, char tran, int N_ROWS, int N_COLS) {
	blas_copy(b, target, N_ROWS);
	
	// LAPACKE stuff
	char TRANS = tran ? 'T' : 'N';	//transposed? 'N' or 'T'
	int INFO = 3;
	int LDA = N_ROWS;
	int LDB = N_ROWS;
	int NRHS = 1;
	int IPIV[N_ROWS];

	// compute LR factorization of A
	dgetrf_(&N_ROWS, &N_COLS, A_arr, &LDA, IPIV, &INFO);

	if(INFO) {
		printf("Error in solve (LR) %d\n", INFO);
		exit(1);
	} else {
		// solve Ax=b
		dgetrs_(&TRANS, &N_ROWS, &NRHS, A_arr, &LDA, IPIV, target, &LDB, &INFO);
		if(INFO) {
			printf("Error in solve %d\n", INFO);
			exit(1);
		}
	}
}

void invertMatrix(double* A_arr, double* workspace, int N) {
	
    int IPIV[N];
    int LWORK = N*N;
    int INFO = 3;

    dgetrf_(&N, &N, A_arr, &N, IPIV, &INFO);
    if(INFO) {
    	printf("error in invert (LR): %d\n", INFO);
    	exit(1);
    }
    dgetri_(&N, A_arr, &N, IPIV, workspace, &LWORK, &INFO); //, 
	if(INFO) {
    	printf("error in invert: %d\n", INFO);
    	exit(1);
    }
}

void computeActiveSet(double* A, double* x, double* b, double* activeA, int M, int N, char	* W) {
	double* temp = malloc(M* sizeof(double));
	matvecprod(A, x, temp, M, N);
	for(int i = 0; i < M; i++) {
		W[i] = (temp[i] - b[i]) <= EPS && (temp[i] - b[i]) >= -EPS;
		//W[i] =  temp[i] == b[i];
	}

	printf("Active Set: ");
	print_vector(W, M, 'c');
	printf("\n\n");
	free(temp);
}

void computeActiveATrans(double* A, char* W, double* activeA, int N, int M) {
	int active = nnz(W, M);
	int idx = 0;
	int noActive = 0;
	for(int i = 0; i < M; i++) {
		if(W[i] && noActive < N) {
			for(int j = 0; j < N; j++) {
				activeA[idx++] = A[j*M+i];
			}
			noActive++;
		}
	}
}

void computeInactiveATrans(double* A, char* W, double* activeA, int N, int M) {
	int active = nnz(W, M);
	int idx = 0;
	for(int i = 0; i < M; i++) {
		if(!W[i]) {
			for(int j = 0; j < N; j++) {
				activeA[idx++] = A[j*M+i];
			}
		}
	}
}

void extendActiveA(double* activeA, double* nonactiveA, int l, int N) {
	int idx = 0;
	for(int i = (N-l)*N; i < N*N; i++){
		activeA[i] = (double)rand()/RAND_MAX;
		//activeA[i] = nonactiveA[idx++];
	}
}

void computeReducedHessian(double* G, double* Z, double* M, int N, int L) {
	if(SOLVER_DEBUG) {
		printf("Analysis of current siuation: \n");
		printf("G: ");
		print_vector(G, N*N, 'd');
		printf("Z: ");
		print_vector(Z, L*N, 'd');
	}
	//compute M = Z' * G * Z
	//first, compute Temp = Z' * G
	double* temp = malloc(N*L*sizeof(double));

	char TRANA = 'T';
	char TRANB = 'N';
	double alpha = 1;
	double beta = 0;
	int lda = N;
	int ldb = N;
	int ldc = L;
	dgemm_(&TRANA, &TRANB, &L, &N, &N, &alpha, Z, &lda, G, &ldb, &beta, temp, &ldc);

	if(SOLVER_DEBUG) {
		printf("temp: ");
		print_vector(temp, N*L, 'd');
	}

	//now compute temp * Z
	matmatprod(temp, Z, M, 'N', 'N', L, N, L);
	if(SOLVER_DEBUG) {
		printf("M: ");
		print_vector(M, L*L, 'd');
	}
	free(temp);
}

void computeU(double* M, double* Z, double * g, double* u, int L, int N) {
	//solve for u: M * u = -Z' *g
	double* temp = malloc(L*sizeof(double));
	//printf("L: %d, N: %d\n", L, N);
	//printf("Z: ");
	//print_vector(Z, L*N, 'd');
	//printf("g: ");
	//print_vector(g, N, 'd');

	char TRANA = 'T';
	char TRANB = 'N';
	double alpha = -1;
	double beta = 0;
	int lda = N;
	int ldb = N;
	int ldc = L;
	int ONE = 1;
	dgemm_(&TRANA, &TRANB, &L, &ONE, &N, &alpha, Z, &lda, g, &ldb, &beta, temp, &ldc);
	solve(M, temp, u, 0, L, L);
	free(temp);
}

void computeSearchDirection(double* A, double* G, double* g, char* W, double* activeA, double* p, int N, int M) {
	int active = nnz(W,M);
	printf("%d constraints active.\n", active);
	double* nonactiveA = malloc(N*(M-active)*sizeof(double));

	computeInactiveATrans(A, W, nonactiveA, N, M);
	if(SOLVER_DEBUG) {
		printf("inactive A: ");
		print_vector(nonactiveA, N*(M-active), 'd');
	}
	int l = N - active;
	if(l == 0) { 
		//vertex solution, V is vacuous, => p = 0
		zeros(p, N);
	} else {
		printf("Extend Active A by %d columns\n", l);
		extendActiveA(activeA, nonactiveA, l, N);
		if(SOLVER_DEBUG) {
			printf("Extended Active A: ");
			print_vector(activeA, N*N, 'd');
		}
		//inverse transpose of active A
		double* workspace = malloc(N*N*sizeof(double));
		invertMatrix(activeA, workspace, N);
		if(SOLVER_DEBUG) {
			printf("Inverted extended active A: ");
			print_vector(activeA, N*N, 'd');
		}

		free(workspace);
		//instead of transposing the inverse and taking columns, we just take the last l rows -> become columns of Z
		double* Z = malloc(N*l*sizeof(double));
		double* row = malloc(N*sizeof(double));
		for(int i = 0; i < l; i++) {
			selectRow(activeA, i+(N-l), row, N, N);
			copyVector(row, &(Z[i*N]), N);
		}
		free(row);
		if(SOLVER_DEBUG) {
			printf("Z: ");
			print_vector(Z, N*l, 'd');
		}
		double* u = malloc(l * sizeof(double));
		double* M_arr = malloc(l*l*sizeof(double));
		//compute M = Z' * G * Z
		computeReducedHessian(G, Z, M_arr, N, l);
		//if M pos def:
		//solve for u: M * u = -Z' *g
		computeU(M_arr, Z, g, u, l, N);
		free(M_arr);
		matvecprod(Z, u, p, N, l); //p = Z*u
		free(u);
		free(Z);
	}
	free(nonactiveA);
}

int findConstraintToRemoveFromActiveSet(double* activeA, double* g, char* W, char* Inequality, int N, int M){
	double* lambda = malloc(N*sizeof(double));
	solve(activeA, g, lambda,0, N, N);
	
	if(SOLVER_DEBUG){
		printf("g: ");
		print_vector(g, N, 'd');
		printf("lambda: ");
		print_vector(lambda, N, 'd');
	}

	int lambdaIdx = -1; int minLambdaConstraintIdx = -1;
	int minLambdaIdx = -1;
	double minLambdaVal = INF;
	for(int i = 0; i < M; i++) {
		if(Inequality[i] && W[i]) {
			if(lambda[++lambdaIdx] < 0) {
				if(lambda[lambdaIdx] < minLambdaVal) {
					minLambdaVal = lambda[lambdaIdx];
					minLambdaIdx = lambdaIdx;
					minLambdaConstraintIdx = i;
				}
			}
		}
	}
	free(lambda);
	return minLambdaConstraintIdx;
}

double computeSteplength(double* A_arr, double* b, double* x, double* p, char* W, int N, int M) {
	double* temp = malloc(M*sizeof(double));
	double* temp2 = malloc(M*sizeof(double));
	matvecprod(A_arr, p, temp,  M, N);

	double* row = malloc(N*sizeof(double));
	int minIdx = -1;
	double minValue = INF;
	for(int i = 0; i < M; i++) {
		if(!W[i] && temp[i] < 0) {
			selectRow(A_arr, i, row, M, N);
			double upper = (b[i] - blas_dot(row, x, N));
			double lower = blas_dot(row, p, N);
			temp2[i] = upper / lower;
			if(temp2[i] < minValue) {
				minValue = temp2[i];
				minIdx = i;
			}
		}
	}
	if (minValue < 1) {
		W[minIdx] = 1;
		printf("add constraint %d to active set.\n", minIdx);
	}
	free(row);
	free(temp);
	free(temp2);
	return min(1, minValue);
}


void nullspaceActiveSetMethod(double* G_arr, double* d, double* A_arr, double* b, char* Inequality, double* x, int N, int M){
	char* W = malloc(M * sizeof(char));
	double* p = malloc(N * sizeof(double));
	
	double* activeA = malloc(N*N*sizeof(double));

	double alpha = 0;
	
	//compute active set
	computeActiveSet(A_arr, x, b, activeA, M, N, W);

	int iter = 1;
	char stop = 0;
	double* g = malloc(N*sizeof(double));
	while(!stop && iter < 100) {
		printf("\nIteration %d start.\n", iter);
		
		matvecprodadd(G_arr, x, d, g, N, N);
		computeActiveATrans(A_arr, W, activeA, N, M);
		//exit(1);
		if(SOLVER_DEBUG) {
			printf("Active A: ");
			print_vector(activeA, N*N, 'd');
		}
		computeSearchDirection(A_arr, G_arr, g, W, activeA, p, N, M);
		if(SOLVER_DEBUG) {
			printf("Direction: ");
			print_vector(p, N, 'd');
		}
		if (nnz_d(p, N) == 0) {//if p = 0
			int minLambdaConstraintIdx = findConstraintToRemoveFromActiveSet(activeA, g, W, Inequality, N, M);
			if(minLambdaConstraintIdx < 0) {
				printf("done!\n\n");
				stop=1;
			} else {
				printf("remove constrait %d from active set.\n", minLambdaConstraintIdx);
				W[minLambdaConstraintIdx] = 0;
			}
		} else {
			double alpha = computeSteplength(A_arr, b, x, p, W, N, M);
			printf("steplength: %lf\n", alpha);
			blas_axpy(p, x, x, alpha, N);
			if(SOLVER_DEBUG) {
				printf("new x: ");
				print_vector(x, N, 'd');
			}
		}	
		iter++;
	}
	free(g);
	free(activeA);
	free(p);
	free(W);
}

int main(void) {
	srand(time(NULL)); 
	

	FILE* f = fopen("testcase01.qplib", "r");
	struct SOLVER_INSTANCE instance = loadInstance(f);
	fclose(f);

	double* G_arr = malloc(instance.N * instance.N * sizeof(double));
	double* A_arr = malloc(instance.N * instance.M * sizeof(double));
	mat2vec(instance.Hessian, G_arr, instance.N, instance.N, 0);
	mat2vec(instance.A, A_arr, instance.M, instance.N, 0);

	char* Inequality = malloc(instance.M * sizeof(char));
	int i;
	for(i=0; i < instance.M; i++) {
		Inequality[i] = 1;
	}

	nullspaceActiveSetMethod(G_arr, instance.d, A_arr, instance.b, Inequality, instance.x0, instance.N, instance.M);

	printf("Result:\n");
	print_vector(instance.x0, instance.N, 'd');


	free(G_arr);
	free(A_arr);
	free(Inequality);

	//free stuff from instance

}