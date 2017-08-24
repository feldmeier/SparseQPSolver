HMatrix getHessian() {
	int n = 2;
	HMatrix hessian;
	int* hstart = new int[n+1];
	int* hindex = new int[n];
	double* hvalue = new double[n];
	int* nonbasic = new int[n];
	for(int i = 0; i < n; i++) {
		hstart[i] = i;
		hindex[i] = i;
		hvalue[i] = 2;
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
		hvalue[i] = 2;
		basicidx[i] = i;
	}
	hstart[n] = n;
	factor.setup(n, n, hstart, hindex, hvalue, basicidx);
	factor.build();
	return factor;
}

HVector getExampleFeasibleStartVector() {
	HVector x;
	double* x_vec = new double[2];
	int* x_base = new int[2];
	x_vec[0] = 2;
	//x_vec[1] = 0;
	x_base[0] = 0;
	//x_base[1] = 1;
	x = createVector(x_vec, x_base, 1, 2);
	return x;
}

HVector getExampleCostVector() {
	HVector c;
	double* c_vec = new double[2];
	int* c_base = new int[2];
	c_vec[0] = -2;
	c_vec[1] = -5;
	c_base[0] = 0;
	c_base[1] = 1;
	c = createVector(c_vec, c_base, 2, 2);
	return c;
}

HFactor getExampleConstraintMatrix() {
	HFactor A;
	int* astart = new int[6];
	int* aindex = new int[10];
	double* avalue = new double[10];
	int* basicidx = new int[2];
	astart[0] = 0;
	astart[1] = 2;
	astart[2] = 4;
	astart[3] = 6;
	astart[4] = 8;
	astart[5] = 10;
	aindex[0] = 0;
	aindex[1] = 1;
	aindex[2] = 0;
	aindex[3] = 1;
	aindex[4] = 0;
	aindex[5] = 1;
	aindex[6] = 0;
	aindex[7] = 1;
	aindex[8] = 0;
	aindex[9] = 1;
	avalue[0] = 1;
	avalue[1] = -2;
	avalue[2] = -1;
	avalue[3] = -2;
	avalue[4] = -1;
	avalue[5] = 2;
	avalue[6] = 1;
	avalue[7] = 0;
	avalue[8] = 0;
	avalue[9] = 1;
	basicidx[0] = 2;
	basicidx[1] = 4;
	A.setup(5, 2, astart, aindex, avalue, basicidx);
	A.build();
	return A;
}


HVector getExampleRHSVector() {
	double* b_vec = new double[5];
	int* b_base = new int[5];
	b_vec[0] = -2;
	b_vec[1] = -6;
	b_vec[2] = -2;
	b_vec[3] = 0;
	b_vec[4] = 0;
	b_base[0] = 0;
	b_base[1] = 1;
	b_base[2] = 2;
	b_base[3] = 3;
	b_base[4] = 4;
	HVector b = createVector(b_vec, b_base, 5, 5);
	return b;
}

bool* getExampleInequality() {
	bool* ineq = new bool[5];
	ineq[0] = true;
	ineq[1] = true;
	ineq[2] = true;
	ineq[3] = true;
	ineq[4] = true;
	return ineq;
}
