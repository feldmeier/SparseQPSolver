/*
 * This is the example for the dual active set method
 */

HMatrix getHessian() {
	int n = 3;
	HMatrix hessian;
	int* hstart = new int[n+1];
	int* hindex = new int[n];
	double* hvalue = new double[n];
	int* nonbasic = new int[n];
	for(int i = 0; i < n; i++) {
		hstart[i] = i;
		hindex[i] = i;
		hvalue[i] = 1;
		nonbasic[i] = 1;
	}
	hstart[n] = n;
	hessian.setup(n, n, hstart, hindex, hvalue, nonbasic);
	return hessian;
}

HFactor getDiagonalHessianFactor(int n) {
	HFactor factor;
	int* hstart = new int[n+1];
	int* hindex = new int[n];
	double* hvalue = new double[n];
	int* basicidx = new int[n];
	for(int i = 0; i < n; i++) {
		hstart[i] = i;
		hindex[i] = i;
		hvalue[i] = 1;
		basicidx[i] = i;
	}
	hstart[n] = n;
	factor.setup(n, n, hstart, hindex, hvalue, basicidx);
	factor.build();
	return factor;
}

HVector getExampleFeasibleStartVector() {
	HVector x;
	double* x_vec = new double[3];
	int* x_base = new int[3];
	x_vec[0] = 2;
	x_vec[1] = 0;
	x_vec[2] = 0;
	x_base[0] = 0;
	x_base[1] = 1;
	x_base[2] = 2;
	x = createVector(x_vec, x_base, 3, 3);
	return x;
}

HVector getExampleCostVector() {
	HVector c;
	double* c_vec = new double[3];
	int* c_base = new int[3];
	c_vec[0] = 0;
	c_vec[1] = -5;
	c_vec[2] = 0;
	c_base[0] = 0;
	c_base[1] = 1;
	c_base[2] = 2;
	c = createVector(c_vec, c_base, 3, 3);
	return c;
}

HFactor getExampleConstraintMatrix() {
	HFactor A;
	int* astart = new int[4];
	int* aindex = new int[6];
	double* avalue = new double[6];
	int* basicidx = new int[3];
	astart[0] = 0;
	astart[1] = 2;
	astart[2] = 4;
	astart[3] = 6;

	aindex[0] = 0;
	aindex[1] = 1;
	aindex[2] = 0;
	aindex[3] = 1;
	aindex[4] = 1;
	aindex[5] = 2;

	avalue[0] = -4;
	avalue[1] = -3;
	avalue[2] = 2;
	avalue[3] = 1;
	avalue[4] = -2;
	avalue[5] = 1;

	
	basicidx[0] = 0;
	basicidx[1] = 1;
	basicidx[2] = 2;
	A.setup(3, 3, astart, aindex, avalue, basicidx);
	A.build();
	return A;
}


HVector getExampleRHSVector() {
	double* b_vec = new double[3];
	int* b_base = new int[3];
	b_vec[0] = -8;
	b_vec[1] = 2;
	b_vec[2] = 0;
	b_base[0] = 0;
	b_base[1] = 1;
	b_base[2] = 2;
	HVector b = createVector(b_vec, b_base, 3, 3);
	return b;
}

bool* getExampleInequality() {
	bool* ineq = new bool[3];
	ineq[0] = true;
	ineq[1] = true;
	ineq[2] = true;
	return ineq;
}
