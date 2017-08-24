#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#define JAJH_dev 2
#include "HModel.h"
#include "HMatrix.h"

#include "Routines.h"

#include "loadExample1.c"


void updateBasis(HFactor& BasisFactor, int in, int out, int p) {
	//in: index within ALL constraints
	//out: index within the basic constraints (e.g.: basis{2,4}, 4 leaves basis -> out = 1
	HVector column;
	column.setup(BasisFactor.numRow);
	HVector row_ep;
	row_ep.setup(BasisFactor.numRow);

	//column entering the basis
	HMatrix ConstraintMatrix;
	factor_to_matrix(BasisFactor, ConstraintMatrix);
	ConstraintMatrix.collect_aj(column, in, 1.0);
	BasisFactor.ftran(column, 1);

	//p = 0 or p = 1?
	//row p of inverse of B
	row_ep = createUnitVector(p, BasisFactor.numRow);
	row_ep.packFlag = true;
	BasisFactor.btran(row_ep, 1);
	int hint = 0;
	BasisFactor.update(&column, &row_ep, &out, &hint);
	if(hint) {
		cout << "Better to reinvert matrix" << endl;
	} else {
		cout << "Update successfull" << endl;
	}
}


/*
 * activeIdx active constraints in basis
 * nonactiveIdx nonactiveconstraints in basis
 */
void computeNullspaceMatrix(int* nonactiveIdx, int nonactiveCount, HFactor& basis, HMatrix& Z) {
	int rows = basis.numRow;
	vector<HVector> zcols;
	vector<int> perm = basis.getPermute();
	for(int i = 0; i < nonactiveCount; i++) {
		int columnToExtract = perm[nonactiveIdx[i]];
		HVector unit = createUnitVector(columnToExtract, rows);
		basis.btran(unit, 1.0);
		HVector zcol; zcol.setup(rows); zcol.copy(&unit);
		//TODO verify that bug is gone now
		int count = 0;
		for(int j = 0; j < rows; j++){
			//zcol.index[j] = perm[unit.index[j]]; //this should work in theory but doesn't
			if(abs(unit.array[j]) > HSOL_CONST_TINY) {
				zcol.index[count++] = j;
			}
		}
		zcols.push_back(zcol);
	}
	vectors_to_matrix(zcols, Z, nonactiveCount, basis.numRow, true);
}

void computeReducedHessian(HMatrix& Z, HMatrix& Hessian, HMatrix& result) {
	compute_matT_mat_prod(Z, Hessian, result);
	compute_mat_mat_prod(result, Z, result);
}

void computeSearchDirection(HMatrix Z, HMatrix Hessian, HVector gradient, HVector& result) {
	//compute Z'*G*Z * u = - Z'*gk
	HMatrix M;
	computeReducedHessian(Z, Hessian, M);

	HVector rhs;
	rhs.setup(Z.numCol);

	HVector temp;
	temp.setup(Z.numCol);
	compute_matT_vec_prod(Z, gradient, temp);
	rhs.saxpy(-1, &temp);
	HVector u;
	u.setup(M.numRow);
	u.copy(&rhs);
	result.setup(Z.numRow);
	cout << "Reduced Hessian matrix" << endl;
	M.rp_mtx();
	if(Z.numCol > 1) {
		HFactor Mfac = matrix_to_factor(M);
		vector<int> permute = Mfac.getPermute();
		//TODO verify correctness in general case
		Mfac.ftran(rhs, 1.0);
		for(int i = 0; i < rhs.count; i++) {
			u.index[i] = i;
			u.array[i] = rhs.array[permute[i]];
		}
		Z.compute_matB_vec(&u.array[0], &u.index[0], &result);
	} else {
		double* val = new double[1];
		int* base = new int[1];
		base[0] = 0;
		val[0] = rhs.array[0] / (M.getAvalue()[0]);
		HVector u = createVector(val, base, 1, 1);
		Z.compute_matB_vec(&u.array[0], &u.index[0], &result);
		delete[] val;
		delete[] base;
	}
}

void computeLagrangeMultipliers(HFactor& A, HVector& gradient, HVector& lambda) {
	//TODO might lead to errors depending on permute?
	HVector temp;
	temp.setup(A.numRow);
	temp.copy(&gradient);
	vector<int> permute = A.getPermute();
	A.ftran(temp, 1.0);
	for(int i = 0; i < /*A.numRow*/ temp.count; i++) {
		//lambda.index[i] = i;
		//lambda.array[i] = temp.array[permute[i]];
		//lambda.index[i] = i;
		//lambda.array[i] = temp.array[i];
	}
	lambda.count = temp.count; //A.numRow;
}

double computeSteplength(HMatrix constraints, HVector rhs, HVector p, HVector x, vector<int> nonactiveConstraints, int* minIdx) {
	double alpha = 1;
	*minIdx = -1;

	HVector aitpk; aitpk.setup(constraints.numCol);
	compute_matT_vec_prod(constraints, p, aitpk);
	HVector aitxk; aitxk.setup(constraints.numCol);
	compute_matT_vec_prod(constraints, x, aitxk);
	for(int i = 0; i < nonactiveConstraints.size(); i++) {
		int constraintId = nonactiveConstraints[i];
		//check if a_i*p < 0
		assert(aitpk.index[constraintId] == constraintId);
		if(aitpk.array[constraintId] < 0) {
			assert(rhs.index[constraintId] == constraintId);
			double bi = rhs.array[constraintId];
			double a = (bi - aitxk.array[constraintId]) / aitpk.array[constraintId];
			if(a < alpha) {
				alpha = a;
				*minIdx = constraintId;
			}
		}
	}

	return alpha;
}


void computeGradient(HMatrix G, HVector c, HVector x, HVector& result) {
	G.compute_matB_vec(&x.packValue[0], &x.packIndex[0], &result);
	result.saxpy(1, &c);
}

void makestep(HVector& x, double alpha, HVector& s) {
	x.tight();
	x.saxpy(alpha, &s);
	x.packFlag = true;
	x.pack();
}



void primal_nullspace_active_set_method(HMatrix Hessian, HVector x, HVector c, HFactor BasisFactor, HVector b, bool* ineq){
	HMatrix ConstraintMatrix;
	factor_to_matrix(BasisFactor, ConstraintMatrix);
	HMatrix Z;
	HVector gradient;
	gradient.setup(Hessian.numCol);
	HVector p;
	HVector lambda;
	lambda.setup(BasisFactor.numRow);
	vector<int> nonactiveConstraints;
	for(int i = 0; i < ConstraintMatrix.numCol; i++) {
		nonactiveConstraints.push_back(i);
	}
	vector<int> activeConstraints;

	bool terminate = false;
	int iter = 1;

	//Initialize basisFactor
	int i;
	int* basicId = new int[BasisFactor.numRow];
	bool* isZcolumn = new bool[BasisFactor.numRow];
	for(i = 0; i < activeConstraints.size(); i++) {
		basicId[i] = activeConstraints[i];
		isZcolumn[i] = false;
	}
	for(int j = 0; i < BasisFactor.numRow; j++) {
		isZcolumn[i] = true;
		basicId[i++] = nonactiveConstraints[j];

	}
	BasisFactor.setup(BasisFactor.numCol, BasisFactor.numRow, (int*)BasisFactor.getAstart(), (int*)BasisFactor.getAindex(), (double*)BasisFactor.getAvalue(), basicId);
	BasisFactor.build();


	while(!terminate) {
		std::cout << "Iteration " << iter << " start" << endl;
		iter++;

		int l = BasisFactor.numRow - activeConstraints.size();
		computeGradient(Hessian, c, x, gradient);
		if(l > 0) {
			int* nonactiveIdx = new int[l];
			int j = 0;
			for(int i = 0; i < BasisFactor.numRow; i++) {
				if(isZcolumn[i])
				nonactiveIdx[j++] = i;
			}
			computeNullspaceMatrix(nonactiveIdx, l, BasisFactor, Z);
			cout << "Nullspace matrix" << endl;
			Z.rp_mtx();
			computeSearchDirection(Z, Hessian, gradient, p);
		} else {
			p.setup(BasisFactor.numRow);
		}
		if(p.isnullvector()) {
			computeLagrangeMultipliers(BasisFactor, gradient, lambda);
			int minIdx = -1;
			int minVal = 0;
			//find inequality constraint with negative Lagrange multiplier
			for(int i = 0; i < lambda.count; i++) {
				if(ineq[lambda.index[i]]) {
					if(lambda.array[i] < minVal) {
						minVal = lambda.array[i];
						minIdx = lambda.index[i];
					}
				}
			}
			if(minIdx < 0) {
				terminate = true;
			} else {
				//remove constraint minIdx from working set
				//no need to update BasisFactor, the column will simply become part of V
				isZcolumn[minIdx] = true;
				nonactiveConstraints.push_back(activeConstraints[minIdx]);
				activeConstraints.erase(std::remove(activeConstraints.begin(), activeConstraints.end(), activeConstraints[minIdx]), activeConstraints.end());
			}
		} else {
			//compute steplength alpha
			int minIdx;
			double alpha = computeSteplength(ConstraintMatrix, b, p, x, nonactiveConstraints, &minIdx);
			makestep(x, alpha, p);
			if(alpha < 1) {
				//TODO clear up code, very dirty
				//add one of the constraints that limits alpha to the working set (i.e., update BasisFactor)
				vector<int> perm = BasisFactor.getPermute();
				int j = -1;
				//check if the column minIdx is already in the basis
				for(int i = 0; i < BasisFactor.numRow; i++){
					if(basicId[i] == minIdx) {
						assert(isZcolumn[i]);
						j = perm[i];
						isZcolumn[j] = false;
						//no basis update necessary
						break;
					}
				}
				if(j == -1) {
					//if not, find first Z column to leave basis
					for(int i = 0; i < BasisFactor.numRow; i++){
						if(isZcolumn[i]) {
							j = i;
							isZcolumn[j] = false;
							updateBasis(BasisFactor, minIdx, j, j);
							basicId[j] = minIdx;
							break;
						}
					}
				}
				activeConstraints.push_back(minIdx);
				nonactiveConstraints.erase(std::remove(nonactiveConstraints.begin(), nonactiveConstraints.end(), minIdx), nonactiveConstraints.end());
			}
		}
	}
	reportVector(x, "solution");
	std::cout << "Final active set: ";
	for(int i = 0; i < activeConstraints.size(); i++) {
		std::cout << activeConstraints[i] << " ";
	}
}


int main(int argc, char **argv) {
	HMatrix Hessian = getHessian();
	HVector x = getExampleFeasibleStartVector();
	HVector c = getExampleCostVector();
	HFactor BasisFactor = getExampleConstraintMatrix();

	HVector b = getExampleRHSVector();
	bool* ineq = getExampleInequality();


	primal_nullspace_active_set_method(Hessian, x, c, BasisFactor, b, ineq);


	// cleanup memory
	delete[] ineq;
}
