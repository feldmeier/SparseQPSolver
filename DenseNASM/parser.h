#ifndef PARSER_H
#define PARSER_H
#endif


#ifndef PARSER_DEBUG
	#define PARSER_DEBUG 0
#endif

#ifndef INF
#define INF 10E18
#endif

#ifndef EPS
#define EPS 10E-8
#endif

struct SOLVER_INSTANCE {
	int N;
	int M;
	double** Hessian;
	double* d;
	double** A;
	double* b;
	double* x0;
};

struct SOLVER_INSTANCE loadInstance(FILE* f);