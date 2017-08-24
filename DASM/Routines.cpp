#include "Routines.h"


/* VECTOR ROUTINES */

//outputs vector v
void reportVector(HVector v, const char* name) {
	std::cout << "Vector " << name <<  "\n";
	std::cout << "Idx\tVal" <<  "\n";
	for(int i = 0; i < v.count; i++) {
		if(v.array[v.index[i]] != HSOL_CONST_ZERO && v.array[v.index[i]] != 0)
			std::cout << v.index[i] << "\t" << v.array[v.index[i]] << "\n";
	}
}

//deletes the element at position deleIdx from vector u
void delete_element_from_vector(HVector* u, int current_length, int deleteIdx) {
	HVector new_u;
	int* base = new int[current_length - 1];
	double* val = new double[current_length - 1];
	int idx = 0;
	for(int i = 0; i < current_length; i++) {
		if(i == deleteIdx) {
			continue;
		}
		base[idx] = idx;
		val[idx] = u->array[i];
		idx++;
	}
	new_u.setup(current_length - 1, base, val);
	reportVector(new_u, "diminished dual variables");
	*u = new_u;
}

HVector createUnitVector(int idx, int length) {
	double* vec = new double[1];
	int* base = new int[1];
	vec[0] = 1;
	base[0] = idx;
	return createVector(vec, base, 1, length);
}

HVector createVector(double* vec, int* base, int nonzeros, int length) {
	HVector result;
	result.setup(length);
	for(int i = 0; i < length; i++) {
		result.index[i] = i;
		result.array[i] = 0;
	}

	for(int i = 0; i < nonzeros; i++) {
		result.array[base[i]] = vec[i];
	}
	result.count = length;
	result.packFlag = true;
	result.pack();
	return result;
}

//computes x = x + p * y
void saxpy(HVector& x, double p, HVector& y) {
	//TODO: Implement as efficient sparse operation
	int* index = new int[x.size];
	double* array = new double[x.size];
	for(int i = 0; i < x.size; i++) {
		index[i] = i;
		array[i] = 0;
	}
	for(int i = 0; i < x.count; i++) {
		array[x.index[i]] = x.array[x.index[i]];
	}
	for(int i = 0; i < y.count; i++) {
		array[y.index[i]] += p * y.array[y.index[i]];
	}
	x.setup(x.size, index, array);
}

// appends an element at the end of a vector
void append_element_to_vector(HVector* u, int current_length, double value) {
	HVector new_u;
	int* base = new int[current_length + 1];
	double* val = new double[current_length + 1];
	for(int i = 0; i < current_length; i++) {
		base[i] = i;
		val[i] = u->array[i];
	}
	base[current_length] = current_length;
	val[current_length] = value;
	new_u.setup(current_length + 1, base, val);
	reportVector(new_u, "augmented dual variables");
	*u = new_u;
}

/* CONVERSION ROUTINES */

// convertes HMatrix to HFactor
HFactor matrix_to_factor(HMatrix& mat) {
	HFactor res;
	int* basicIdx = new int[mat.numRow];
	for(int i = 0; i < mat.numRow; i++){
		basicIdx[i] = i;
	}
	int* start = (int*)mat.getAstart();
	int* index = (int*)mat.getAindex();
	double* value = (double*)mat.getAvalue();

	res.setup(mat.numCol, mat.numRow, start, index, value, basicIdx);
	res.build();
	return res;
}

// convertes HFactor to HMatrix
void factor_to_matrix(HFactor& fac, HMatrix& mat) {
	int* start = (int*)fac.getAstart();
	int* index = (int*)fac.getAindex();
	double* value = (double*)fac.getAvalue();
	mat.setup_lgBs(fac.numCol, fac.numRow, start, index, value);
}

// convertes a number of vectors to HMatrix (indicate if rowvector or columnvector)
void vectors_to_matrix(vector<HVector> vectors, HMatrix& mat, int numCol, int numRow, bool columnwise) {
	//count elements of matrix
	int elementCount = 0;
	for(int i = 0; i < vectors.size(); i++) {
		vectors[i].packFlag = true;
		vectors[i].pack();
		elementCount += vectors[i].count;
	}
	//build matrix
	int* matstart = new int[numCol + 1];
	int* matindex = new int[elementCount];
	double* matvalue = new double[elementCount];

	if(columnwise) { // each vector is a column of the matrix
		int entryCounter = 0;
		for(int i = 0; i < numCol; i++) {
			matstart[i] = entryCounter;
			HVector column = vectors[i];
			for(int j = 0; j < column.count; j++) {
				matindex[entryCounter] = column.packIndex[j];
				matvalue[entryCounter] = column.packValue[j];
				entryCounter++;
			}
		}
		matstart[numCol] = entryCounter;
	}
	else { //each vector is a row of the matrix
		int entryCounter = 0;
		for(int i = 0; i < numCol; i++) {
			matstart[i] = entryCounter;
			//in all vectors see if they have an entry for index i
			for(int j = 0; j < numRow; j++) {
				HVector row = vectors[j];
				for(int k = 0; k < row.count; k++) {
					if(row.packIndex[k] == i) {
						matindex[entryCounter] = j;
						matvalue[entryCounter] = row.packValue[k];
						entryCounter++;
						break;
					}
				}
			}
		}
		matstart[numCol] = entryCounter;
	}
	mat.setup_lgBs(numCol, numRow, matstart, matindex, matvalue);
}


/* MATRIX ROUTINES */
// computes r = A^T * v
void compute_matT_vec_prod(HMatrix A, HVector v, HVector& r) {
	int* resbase = new int[A.numCol];
	double* resvalue = new double[A.numCol];

	int* astart = (int*)A.getAstart();
	int* aindex = (int*)A.getAindex();
	double* avalue = (double*)A.getAvalue();

	for(int i = 0; i < A.numCol; i++) {
		//compute column i of A times vector v
		double val = 0;
		for(int j = astart[i]; j < astart[i+1]; j++) {
			int index = aindex[j];
			//check if v has an entry at index
			for(int k = 0; k < v.count; k++) {
				if(v.index[k] == index) {
					val += avalue[j] * v.array[v.index[k]];
				}
			}
		}
		resbase[i] = i;
		resvalue[i] = val;
	}
	r.setup(A.numCol, resbase, resvalue);
}

// computes R = A * B
void compute_mat_mat_prod(HMatrix& A, HMatrix& B, HMatrix& R) {
	assert(A.numCol = B.numRow);
	vector<HVector> rows;
	int* arstart = (int*)A.getARstart();
	int* arindex = (int*)A.getARindex();
	double* arvalue = (double*)A.getARvalue();
	for(int i = 0; i < A.numRow; i++) {
		double* rowvalue = &arvalue[arstart[i]];
		int* rowbase = &arindex[arstart[i]];
		int rowlength = arstart[i+1] - arstart[i];
		HVector row = createVector(rowvalue, rowbase, rowlength, A.numCol);
		B.compute_vec_mat(row, row);
		rows.push_back(row);
	}
	vectors_to_matrix(rows, R, B.numCol, A.numRow, false);
}

// computes R = A^T * B
void compute_matT_mat_prod(HMatrix& A, HMatrix& B, HMatrix& R) {
	assert(A.numRow == B.numRow);
	vector<HVector> rows;
	int* astart = (int*)A.getAstart();
	int* aindex = (int*)A.getAindex();
	double* avalue = (double*)A.getAvalue();
	for(int i = 0; i < A.numCol; i++) {
		double* rowvalue = &avalue[astart[i]];
		int* rowbase = &aindex[astart[i]];
		int rowlength = astart[i+1] - astart[i];
		HVector row = createVector(rowvalue, rowbase, rowlength, A.numRow);
		B.compute_vec_mat(row, row);
		rows.push_back(row);
	}
	vectors_to_matrix(rows, R, B.numCol, A.numCol, false);
}

// computes res = mat1 - mat2
void compute_mat_mat_sub(HMatrix& mat1, HMatrix& mat2, HMatrix& res) {
	int* astart = new int[mat1.numCol + 1];
	int* aindex = new int[mat1.numCol * mat1.numRow];
	double* avalue = new double[mat1.numCol * mat1.numRow];

	int* m1start = (int*)mat1.getAstart();
	int* m1index = (int*)mat1.getAindex();
	double* m1value = (double*)mat1.getAvalue();
	int* m2start = (int*)mat2.getAstart();
	int* m2index = (int*)mat2.getAindex();
	double* m2value = (double*)mat2.getAvalue();

	//initialize result matrix: set all to zero
	int idx = 0;
	for(int i = 0; i < mat1.numCol; i++) {
		astart[i] = idx;
		for(int j = 0; j < mat1.numRow; j++) {
			aindex[idx] = j;
			avalue[idx] = 0;
			idx++;
		}
	}
	astart[mat1.numCol] = idx;

	//add mat1
	for(int i = 0; i < mat1.numCol; i++) {
		int colstart = m1start[i];
		int colend = m1start[i+1];
		for(int j = 0; j < colend - colstart; j++) {
			//
			int index = m1index[colstart + j];
			double val = m1value[colstart + j];
			if(abs(val) > HSOL_CONST_TINY)
				avalue[i * mat1.numRow + index] = val;
		}
	}
	//sub mat2
	for(int i = 0; i < mat2.numCol; i++) {
		int colstart = m2start[i];
		int colend = m2start[i+1];
		for(int j = 0; j < colend - colstart; j++) {
			//
			int index = m2index[colstart + j];
			double val = m2value[colstart + j];
			if(abs(val) > HSOL_CONST_TINY)
				avalue[i * mat2.numRow + index] -= val;
		}
	}
	res.setup_lgBs(mat1.numCol, mat1.numRow, astart, aindex, avalue);
}