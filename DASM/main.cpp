#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#define JAJH_dev 2
#include "HModel.h"
#include "HMatrix.h"
#include "Routines.h"

#include "loadExample3.c"


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

void compute_inverse_mat(HMatrix& mat, HMatrix& res) {
	//TODO verify that this is indeed the inverse and not the inverse transpose
	HFactor factor = matrix_to_factor(mat);
	int* idx = new int[mat.numCol];
	for(int i = 0; i < mat.numCol; i++) {
		idx[i] = i;
	}
	computeNullspaceMatrix(idx, mat.numCol, factor, res);
	delete[] idx;
}

void dual_compute_slack(HMatrix constraints, HVector x, HVector b, HVector& slack) {
	slack.setup(constraints.numCol);
	compute_matT_vec_prod(constraints, x, slack);
	saxpy(slack, -1, b);
	reportVector(slack, "Slack");
}

int dual_choose_most_violated_constraint(HVector slack, vector<int> nonactiveSet) {
	int minIdx = -1;
	double minVal = 10E20;
	for(int i = 0; i < nonactiveSet.size(); i++) {
		if(slack.array[nonactiveSet[i]] < minVal) {
			minIdx = nonactiveSet[i];
			minVal = slack.array[nonactiveSet[i]];
		}
	}
	cout << "Violated constraint " << minIdx << " chosen." << endl;
	return minIdx;
}

int dual_choose_first_violated_constraint(HVector slack, vector<int> nonactiveSet) {
	for(int i = 0; i < nonactiveSet.size(); i++) {
		if(slack.array[nonactiveSet[i]] < 0) {
			cout << "Violated constraint " << nonactiveSet[i] << " chosen." << endl;
			return nonactiveSet[i];
		}
	}
	return -1;
}

void dual_compute_N_star(HMatrix& N, HMatrix& Hinv, HMatrix& res, int ndim) {
	if(ndim > 0) {
		//computes N* = (N^T G⁻¹ N)⁻¹N^T G⁻¹
		HMatrix temp1;
		compute_matT_mat_prod(N, Hinv, temp1);
		HMatrix temp2;
		compute_mat_mat_prod(temp1, N, temp2);
		HMatrix temp3;
		compute_inverse_mat(temp2, temp3);
		compute_matT_mat_prod(N, Hinv, temp1);
		compute_mat_mat_prod(temp3, temp1, res);
		cout << "inverse of N" << endl;
		N.rp_mtx();
	}
}

void dual_compute_reduced_inverse_hessian(HMatrix& Hinv, HMatrix& I, HMatrix& N, HMatrix& N_star, HMatrix& res, int ndim) {
	if(ndim > 0) {
		HMatrix temp1;
		compute_mat_mat_prod(N, N_star, temp1);
		HMatrix temp2;
		compute_mat_mat_sub(I, temp1, temp2);
		compute_mat_mat_prod(Hinv, temp2, res);
	} else {
		res.setup_lgBs(Hinv.numCol, Hinv.numRow, (int*)Hinv.getAstart(), (int*)Hinv.getAindex(), (double*)Hinv.getAvalue());
	}
	cout << "reduced inverse hessian" << endl;
	res.rp_mtx();
}

void dual_compute_N(vector<int> activeSet, HMatrix& constraints, HMatrix& N, int ndim) {
	if(ndim > 0) {
		vector<HVector> columns;
		for(int i = 0; i < activeSet.size(); i++) {
			HVector column; column.setup(constraints.numRow);
			constraints.collect_aj(column, activeSet[i], 1);
			columns.push_back(column);
		}
		vectors_to_matrix(columns, N, activeSet.size(), constraints.numRow, true);
		cout << "N" << endl;
		N.rp_mtx();
	}
}

double dual_compute_primal_steplength(HVector& slack, HMatrix& constraints, HVector& primal_search_direction,int q) {
	double alpha_primal = - slack.array[q] / constraints.compute_dot(primal_search_direction, q);
	cout << "Primal step length: " << alpha_primal << endl;
	return alpha_primal;

}

void dual_compute_dual_derivative(HMatrix& NStar, HVector& normal_q, HVector& r, int ndim) {
	if(ndim > 0) {
		NStar.compute_matB_vec(&normal_q.array[0], &normal_q.index[0], &r);
		reportVector(r, "directional derivative of dual variables along primal search direction");
	}
}

double dual_compute_dual_steplength(HVector& r, HVector& u, vector<int> activeSet, bool* ineq, int* limiting_constraint_idx, int* limiting_constraint_pos) {
	double alpha_dual = HSOL_CONST_INF;
	if(activeSet.size() > 0) {
		for(int i = 0; i < activeSet.size(); i++) {
			if(r.array[i] > 0 && ineq[activeSet[i]]) {
				double lim = u.array[i] / r.array[i];
				if(lim < alpha_dual) {
					alpha_dual = lim;
					*limiting_constraint_idx = activeSet[i];
					*limiting_constraint_pos = i;
				}
			}
		}
	}

	cout << "Dual step length: " << alpha_dual << " limited by constraint " << *limiting_constraint_idx << endl;
	return alpha_dual;
}

void dual_variables_update(HVector& u, HVector& r, double alpha) {
	append_element_to_vector(&r, r.count, -1);
	saxpy(u, -alpha, r);
	reportVector(u, "updated dual variables");
}


void dual_active_set_method(HMatrix H, HVector c, HMatrix constraints, HVector b, bool* ineq, HMatrix I) {
	//Initialise active set
	vector<int> activeSet;
	vector<int> nonactiveSet;
	for(int i = 0; i < constraints.numCol; i++) {
		nonactiveSet.push_back(i);
	}

	//primal variables
	HVector x;	x.setup(H.numCol);
	//dual variables
	HVector u; u.setup(constraints.numCol);

	HMatrix Hinv;
	compute_inverse_mat(H, Hinv);
	HMatrix N;
	HMatrix NStar;
	HMatrix Minv;
	int ndim = 0;

	//compute start point
	HVector tempx; tempx.setup(H.numCol);
	Hinv.compute_matB_vec(&c.array[0], &c.index[0], &tempx);
	saxpy(x, -1, tempx);

	int iteration = 1;
	bool terminate = false;
	bool skip_augment_u = false;
	while(!terminate) {
		cout << "Iteration " << iteration << " start." << endl;
		dual_compute_N(activeSet, constraints, N, ndim);
		dual_compute_N_star(N, Hinv, NStar, ndim);
		dual_compute_reduced_inverse_hessian(Hinv, I, N, NStar, Minv, ndim);
		HVector slack;
		dual_compute_slack(constraints, x, b, slack);
		int q = dual_choose_first_violated_constraint(slack, nonactiveSet);
		if(q != -1) {
			if(skip_augment_u) {
				skip_augment_u = false;
			} else {
				append_element_to_vector(&u, u.count, 0);
			}
			HVector normal_q; normal_q.setup(H.numCol);
			constraints.collect_aj(normal_q, q, 1);
			saxpy(normal_q, 0, normal_q);
			reportVector(normal_q, "nq");
			HVector primal_search_direction; primal_search_direction.setup(H.numCol);
			Minv.compute_matB_vec(&normal_q.array[0], &normal_q.index[0], &primal_search_direction);
			reportVector(primal_search_direction, "primal search direction");
			HVector r; r.setup(ndim);
			dual_compute_dual_derivative(NStar, normal_q, r, ndim);
			int limiting_constraint_idx = -1;
			int limiting_constraint_pos = -1;
			double alpha_dual = dual_compute_dual_steplength(r, u, activeSet, ineq, &limiting_constraint_idx, &limiting_constraint_pos);
			if(!primal_search_direction.isnullvector()) {
				double alpha_primal = dual_compute_primal_steplength(slack, constraints, primal_search_direction, q);
				double alpha = min(alpha_primal, alpha_dual);
				saxpy(x, alpha, primal_search_direction);
				reportVector(x, "new iterate");
				dual_variables_update(u, r, alpha);
				if(alpha == alpha_primal) {
					activeSet.push_back(q);
					nonactiveSet.erase(std::remove(nonactiveSet.begin(), nonactiveSet.end(), q), nonactiveSet.end());
					ndim++;
				} else {
					nonactiveSet.push_back(limiting_constraint_idx);
					activeSet.erase(std::remove(activeSet.begin(), activeSet.end(), limiting_constraint_idx), activeSet.end());
					ndim--;
					//delete u component
					delete_element_from_vector(&u, u.count, limiting_constraint_pos);
					skip_augment_u = true;
				}
			} else {
				if(alpha_dual == HSOL_CONST_INF) {
					cout << "Dual problem unbounded, primal problem infeasible." << endl;
					terminate = true;
				} else {
					dual_variables_update(u, r, alpha_dual);
					nonactiveSet.push_back(limiting_constraint_idx);
					activeSet.erase(std::remove(activeSet.begin(), activeSet.end(), limiting_constraint_idx), activeSet.end());
					ndim--;
					//delete u component
					delete_element_from_vector(&u, u.count, limiting_constraint_pos);
					skip_augment_u = true;
				}
			}
		} else {
			terminate = true;
			cout << "Optimal solution found." << endl;
			reportVector(x, "Optimal primal solution");
			reportVector(u, "Optimal dual solution");
		}
		iteration++;
	}
}





int main(int argc, char **argv) {
	HMatrix Hessian = getHessian();
	HVector x = getExampleFeasibleStartVector();
	HVector c = getExampleCostVector();
	HFactor BasisFactor = getExampleConstraintMatrix();
	HMatrix constraints;
	factor_to_matrix(BasisFactor, constraints);

	HVector b = getExampleRHSVector();
	bool* ineq = getExampleInequality();

	//setup identity matrix
	HMatrix I;
	int* id_start = new int[Hessian.numCol+1];
	int* id_index = new int[Hessian.numCol];
	double* id_value = new double[Hessian.numCol];
	for(int i = 0; i < Hessian.numCol; i++) {
		id_start[i] = i;
		id_index[i] = i;
		id_value[i] = 1;
	}
	id_start[Hessian.numCol] = Hessian.numCol;
	I.setup_lgBs(Hessian.numCol, Hessian.numCol, id_start, id_index, id_value);

	dual_active_set_method(Hessian, c, constraints, b, ineq, I);



	// cleanup memory
	delete[] ineq;
}
