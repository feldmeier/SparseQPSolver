HMatrix getHessian() {
	int n = 6;
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
	double* x_vec = new double[6];
	int* x_base = new int[6];
	x_vec[0] = 1;
	x_vec[1] = 0;
	x_vec[2] = 1;
	x_vec[3] = 0;
	x_vec[4] = 0;
	x_vec[5] = 1;
	x_base[0] = 0;
	x_base[1] = 1;
	x_base[2] = 2;
	x_base[3] = 3;
	x_base[4] = 4;
	x_base[5] = 5;
	x = createVector(x_vec, x_base, 6, 6);
	return x;
}

HVector getExampleCostVector() {
	HVector c;
	double* c_vec = new double[6];
	int* c_base = new int[6];
	c_vec[0] = 0;
	c_vec[1] = 0;
	c_vec[2] = 0;
	c_vec[3] = 1;
	c_vec[4] = 0;
	c_vec[5] = 0;
	c_base[0] = 0;
	c_base[1] = 1;
	c_base[2] = 2;
	c_base[3] = 3;
	c_base[4] = 4;
	c_base[5] = 5;
	c = createVector(c_vec, c_base, 6, 6);
	return c;
}

HFactor getExampleConstraintMatrix() {
	//try doing this non-sparse and see if it works
	HFactor A;
	int* astart = new int[9];
	int* aindex = new int[28];
	double* avalue = new double[28];
	int* basicidx = new int[6];
	for(int i = 0; i < 28; i++) {
		avalue[i] = 1;
	}
	astart[0] = 0;
	aindex[0] = 0;
	aindex[1] = 3;
	aindex[2] = 4;
	aindex[3] = 5;
	astart[1] = 4;
	aindex[4] = 0;
	aindex[5] = 1;
	aindex[6] = 3;
	aindex[7] = 4;
	astart[2] = 8;
	aindex[8] = 0;
	aindex[9] = 1;
	aindex[10] = 5;
	astart[3] = 11;
	aindex[11] = 1;
	aindex[12] = 2;
	aindex[13] = 4;
	aindex[14] = 5;
	astart[4] = 15;
	aindex[15] = 0;
	aindex[16] = 2;
	aindex[17] = 5;
	astart[5] = 18;
	aindex[18] = 0;
	aindex[19] = 2;
	astart[6] = 20;
	aindex[20] = 0;
	aindex[21] = 1;
	aindex[22] = 2;
	aindex[23] = 3;
	aindex[24] = 5;
	astart[7] = 25;
	aindex[25] = 1;
	aindex[26] = 3;
	aindex[27] = 5;
	astart[8] = 28;
	
	basicidx[0] = 0;
	basicidx[1] = 1;
	basicidx[2] = 2;
	basicidx[3] = 3;
	basicidx[4] = 4;
	basicidx[5] = 5;
	A.setup(8, 6, astart, aindex, avalue, basicidx);
	A.build();
	return A;
}


HVector getExampleRHSVector() {
	double* b_vec = new double[8];
	int* b_base = new int[8];
	b_vec[0] = 2;
	b_vec[1] = 1;
	b_vec[2] = 2;
	b_vec[3] = 2;
	b_vec[4] = 3;
	b_vec[5] = 2;
	b_vec[6] = 3;
	b_vec[7] = 1;
	b_base[0] = 0;
	b_base[1] = 1;
	b_base[2] = 2;
	b_base[3] = 3;
	b_base[4] = 4;
	b_base[5] = 5;
	b_base[6] = 6;
	b_base[7] = 7;
	HVector b = createVector(b_vec, b_base, 8, 8);
	return b;
}

bool* getExampleInequality() {
	bool* ineq = new bool[8];
	ineq[0] = true;
	ineq[1] = false; //ineq[1] = true;
	ineq[2] = true;
	ineq[3] = true;
	ineq[4] = true;
	ineq[5] = false; //ineq[5] = true;
	ineq[6] = true;
	ineq[7] = false; //ineq[7] = true;
	return ineq;
}
