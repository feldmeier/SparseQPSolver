#include <stdio.h>
#include <stdlib.h>

#include "parser.h"

struct QPLIB_INSTANCE {
	char name[20];
	int n;
	int m;
	int entries;
	double** hessian;
	double* d;
	double** coefficients;
	double* LHS;
	double* RHS;
	double* lower;
	double* upper;
	double objectiveConstant;
	double* primalStart;
};



void readuntilendofline(FILE* fp) {
	char c = 'a';
	while(c != '\n') {
		fscanf(fp, "%c", &c);
	}
}

void ignoreLine(FILE* fp) {
	char buf[255];
	fscanf(fp, "%s", buf);
}

struct SOLVER_INSTANCE loadInstance(FILE* f) {
	char buf[255];

	struct QPLIB_INSTANCE instance;
	fscanf(f, "%s", instance.name);	//read instance id
	ignoreLine(f);	//ignore 3 digit problem type code
	ignoreLine(f); //ignore minimize/maximize
	fscanf(f, "%d", &instance.n);	//read number of variables
	readuntilendofline(f);	
	fscanf(f, "%d", &instance.m); //read number of constraints
	readuntilendofline(f);
	fscanf(f, "%d", &instance.entries); //read number of nonzero entries in hessian
	readuntilendofline(f);


	printf("Reading instance %s from file. Instance has:\n%d variables\n%d constraints\n%d nonzero entries in hessian\n\nAllocating memory...", instance.name, instance.n, instance.m, instance.entries);

	instance.d = malloc(instance.n * sizeof(double));
	instance.hessian = malloc(instance.n * sizeof(double*));
	int i;
	for(i = 0; i < instance.n; i++) {
		instance.hessian[i] = malloc(instance.n * sizeof(double));
	}
	instance.RHS = malloc(instance.m * sizeof(double));
	instance.LHS = malloc(instance.m * sizeof(double));
	instance.coefficients = malloc(instance.m * sizeof(double*));
	for(i = 0; i < instance.m; i++) {
		instance.coefficients[i] = malloc(instance.n * sizeof(double));
	}

	instance.upper = malloc(instance.n * sizeof(double));
	instance.lower = malloc(instance.n * sizeof(double));

	instance.primalStart = malloc(instance.n * sizeof(double));

	printf("initializing values..\n");

	int j;
	for(i = 0; i < instance.n; i++) {
		for(j = 0; j < instance.n; j++) {
			instance.hessian[i][j] = 0.0;
		}
		instance.lower[i]= 0.0;
		instance.upper[i]= 0.0;
	}
	for(i = 0; i < instance.m; i++) {
		for(j = 0; j < instance.n; j++) {
			instance.coefficients[i][j] = 0.0;
		}
		instance.RHS[i]= 0.0;
		instance.LHS[i]= 0.0;
	}

	

	printf("Reading hessian values from file..\n");

	for(i = 0; i < instance.entries; i++) {
		int x1, x2;
		double val;
		fscanf(f, "%d %d %lf", &x1, &x2, &val);
		instance.hessian[x1-1][x2-1] = val;
		if(PARSER_DEBUG) printf("%d %d %lf\n", x1, x2, val);
	}

	printf("Reading linear coefficients from file.. \n");

	double defaultVal;
	fscanf(f, "%lf", &defaultVal); //default value for linear coefficients
	readuntilendofline(f); 
	printf("Default value for linear coefficients: %lf.\n", defaultVal);
	int lines;
	fscanf(f, "%d", &lines); //read number of non-default linear coefficients
	printf("%d non-default linear coefficients.\n", lines);
	readuntilendofline(f); 

	printf("Initializing default values..\n");
	for(i = 0; i < instance.n; i++) {
		instance.d[i] = 0.0;
	}

	for(i = 0; i < lines; i++) {
		int x1;
		double val;
		fscanf(f, "%d %lf", &x1, &val);
		instance.d[x1-1] = val;
		if(PARSER_DEBUG) printf("%d %lf\n", x1, val);
	}

	
	fscanf(f, "%lf", &instance.objectiveConstant);
	readuntilendofline(f); 
	printf("Objective Constant: %lf\n", instance.objectiveConstant);

	printf("Reading linear constraints from file..\n");
	fscanf(f, "%d", &lines);
	readuntilendofline(f); 
	printf("%d linear terms in all constraints\n", lines);

	for(i = 0; i < lines; i++) {
		int con, var;
		double val;
		fscanf(f, "%d %d %lf", &con , &var, &val);
		if(PARSER_DEBUG) printf("%d %d %lf\n", con, var, val);
		instance.coefficients[con-1][var-1] = val;
	}

	//value for inf
	double infty;
	fscanf(f, "%lf", &infty);
	readuntilendofline(f);

	//read in LHS values
	double defaultLHS;
	int nLHS;
	fscanf(f, "%lf", &defaultLHS);
	readuntilendofline(f);
	fscanf(f, "%d", &nLHS);
	readuntilendofline(f);
	for(i=0; i < instance.m; i++) {
		instance.LHS[i] = defaultLHS;
	}
	for(i = 0; i < nLHS; i++) {
		int x;
		double val;
		fscanf(f, "%d %lf", &x, &val);
		instance.LHS[x-1] = val;
		if(PARSER_DEBUG) printf("%d %lf\n", x, val);
	}

	//read in RHS values
	double defaultRHS;
	int nRHS;
	fscanf(f, "%lf", &defaultRHS);
	readuntilendofline(f);
	fscanf(f, "%d", &nRHS);
	readuntilendofline(f);
	for(i=0; i < instance.m; i++) {
		instance.RHS[i] = defaultRHS;
	}
	for(i = 0; i < nRHS; i++) {
		int x;
		double val;
		fscanf(f, "%d %lf", &x, &val);
		instance.RHS[x-1] = val;
		if(PARSER_DEBUG) printf("%d %lf\n", x, val);
	}


	//read in lower bounds
	double defaultLower;
	int nLower;
	fscanf(f, "%lf", &defaultLower);
	readuntilendofline(f);
	fscanf(f, "%d", &nLower);
	printf("%d lower bounds:\n", nLower);
	readuntilendofline(f);
	for(i=0; i < instance.n; i++) {
		instance.lower[i] = defaultLower;
	}
	for(i = 0; i < nLower; i++) {
		int x;
		double val;
		fscanf(f, "%d %lf", &x, &val);
		instance.lower[x-1] = val;
		if(PARSER_DEBUG) printf("%d %lf\n", x, val);
	}

	//read in upper bounds
	double defaultUpper;
	int nUpper;
	fscanf(f, "%lf", &defaultUpper);
	readuntilendofline(f);
	fscanf(f, "%d", &nUpper);
	readuntilendofline(f);
	for(i=0; i < instance.n; i++) {
		instance.upper[i] = defaultUpper;
	}
	for(i = 0; i < nUpper; i++) {
		int x;
		double val;
		fscanf(f, "%d %lf", &x, &val);
		instance.upper[x-1] = val;
		if(PARSER_DEBUG) printf("%d %lf\n", x, val);
	}


	//read primal Start point
	double defaultPrimal;
	int nPrimal;
	fscanf(f, "%lf", &defaultPrimal);
	readuntilendofline(f);
	fscanf(f, "%d", &nPrimal);
	readuntilendofline(f);
	for(i=0; i < instance.n; i++) {
		instance.primalStart[i] = defaultPrimal;
	}
	for(i = 0; i < nPrimal; i++) {
		int x;
		double val;
		fscanf(f, "%d %lf", &x, &val);
		instance.primalStart[x-1] = val;
		if(PARSER_DEBUG) printf("%d %lf\n", x, val);
	}
	

	//convert to solver instance
	struct SOLVER_INSTANCE solverInstance;
	//constraints: RHS, LHS -> 2* M, upper, lower -> 2*N
	solverInstance.M = 2*instance.m + 2* instance.n;
	solverInstance.N = instance.n;
	solverInstance.Hessian = instance.hessian;
	solverInstance.d = instance.d;
	solverInstance.x0 = instance.primalStart;


	solverInstance.A = malloc(solverInstance.M * sizeof(double*));
	solverInstance.b = malloc(solverInstance.M * sizeof(double));
	for(i = 0; i < solverInstance.M; i++) {
		solverInstance.A[i] = malloc(solverInstance.N * sizeof(double));
		int j;
		for(j = 0; j < solverInstance.N; j++) {
			solverInstance.A[i][j] = 0.0;
		}
		solverInstance.b[i] = 0.0;
	}

	//constraints and LHS a_iTx >= bi
	printf("Part 1");
	int con = 0;
	for(i = 0; i < instance.m; i++) {
		solverInstance.b[con] = instance.LHS[i];
		if(solverInstance.b[con] > -INF) {
			//printf("Add constraint ");
			int j;
			for(j = 0; j < solverInstance.N; j++) {
				solverInstance.A[con][j] = instance.coefficients[i][j];
				//printf("%lfx%d ", solverInstance.A[con][j], j+1);
			}
			
			//printf(">= %lf\n", solverInstance.b[con]);
			con++;
		}
	} 

	//constraints and RHS
	printf("Part 2");
	for(i = 0; i < instance.m; i++) {
		solverInstance.b[con] = -instance.RHS[i];
		if(solverInstance.b[con] > -INF) {
			//printf("Add constraint ");
			int j;
			for(j = 0; j < solverInstance.N; j++) {
				solverInstance.A[con][j] = -instance.coefficients[i][j];
				//printf("%lfx%d ", solverInstance.A[con][j], j+1);
			}
			
			//printf(">= %lf\n", solverInstance.b[con]);
			con++;
		}
	} 

	//lower bounds
	printf("Part 3");
	for(i = 0; i < instance.n; i++) {
		solverInstance.b[con] = instance.lower[i];
		if(solverInstance.b[con] > -INF) {
			solverInstance.A[con][i] = 1;
			//printf("Add constraint x%d >= %lf\n", i+1, solverInstance.b[con]);
			con++;
		}
	}

	//upper bounds
	printf("Part 4");
	for(i = 0; i < instance.n; i++) {
		solverInstance.b[con] = -instance.upper[i];
		if(solverInstance.b[con] > -INF) {
			solverInstance.A[con][i] = -1;
			//printf("Add constraint x%d >= %lf\n", i+1, solverInstance.b[con]);
			con++;
		}
	}
	solverInstance.M = con;

	//TODO: dual start


	//free all stuff of qplib instance
	free(instance.RHS);
	free(instance.LHS);
	for(i = 0; i < instance.m; i++) {
		free(instance.coefficients[i]);
	}
	free(instance.coefficients);

	free(instance.upper);
	free(instance.lower); 


	printf("Instance reading finished.\n");
	return solverInstance;
}
