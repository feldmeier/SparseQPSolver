#ifndef ROUTINES_H
#define ROUTINES_H

#include "HModel.h"
#include "HMatrix.h"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

void reportVector(HVector v, const char* name);
void delete_element_from_vector(HVector* u, int current_length, int deleteIdx);
HVector createUnitVector(int idx, int length);
HVector createVector(double* vec, int* base, int nonzeros, int length);
void saxpy(HVector& x, double p, HVector& y);
void append_element_to_vector(HVector* u, int current_length, double value);

HFactor matrix_to_factor(HMatrix& mat);
void factor_to_matrix(HFactor& fac, HMatrix& mat);
void vectors_to_matrix(vector<HVector> vectors, HMatrix& mat, int numCol, int numRow, bool columnwise);

void compute_matT_vec_prod(HMatrix A, HVector v, HVector& r);
void compute_mat_mat_prod(HMatrix& A, HMatrix& B, HMatrix& R);
void compute_matT_mat_prod(HMatrix& A, HMatrix& B, HMatrix& R);
void compute_mat_mat_sub(HMatrix& mat1, HMatrix& mat2, HMatrix& res);






#endif