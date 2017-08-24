#include "HVector.h"
#include "HConst.h"

#include <cmath>
#include <cassert>

bool HVector::isnullvector() {
	for(int i = 0; i < count; i++) {
		if(abs(array[index[i]]) > HSOL_CONST_ZERO) {
			return false;
		}
	}
	return true;
}

bool HVector::ispositive() {
	for(int i = 0; i < count; i++) {
		if(array[index[i]] <= 0) {
			return false;
		}
	}
	return true;
}

void HVector::setup(int size_) {
	next = 0;
    size = size_;
    count = 0;
    index.resize(size);
    array.assign(size, 0);
    cwork.assign(size + 6400, 0); // MAX invert
    iwork.assign(size * 4, 0);

    packCount = 0;
    packIndex.resize(size);
    packValue.resize(size);
}

void HVector::setup(int size, int* base, double* value) {
	setup(size);
	for(int i = 0; i < size; i++) {
		index[i] = base[i];
		array[i] = value[i];
	}
	count = size;
	packFlag = true;
	pack();
}

void HVector::clear() {
    int clearVector_inDense = count < 0 || count > size * 0.3;
    if (clearVector_inDense) {
        array.assign(size, 0);
    } else {
        for (int i = 0; i < count; i++)
            array[index[i]] = 0;
    }
    packFlag = false;
    count = 0;
    pseudoTick = 0;
    fakeTick = 0;
    next = 0;
}

void HVector::tight() {
    int totalCount = 0;
    for (int i = 0; i < count; i++) {
        const int my_index = index[i];
        const double value = array[my_index];
        if (fabs(value) > HSOL_CONST_TINY) {
            index[totalCount++] = my_index;
        } else {
            array[my_index] = 0;
        }
    }
    count = totalCount;
}

void HVector::pack() {
    if (packFlag) {
        packFlag = false;
        packCount = 0;

        for (int i = 0; i < count; i++) {
            const int ipack = index[i];
            packIndex[packCount] = ipack;
            packValue[packCount] = array[ipack];
            packCount++;
        }
    }
}

void HVector::copy(const HVector *from) {
    clear();
    fakeTick = from->fakeTick;
    pseudoTick = from->pseudoTick;
    const int fromCount = count = from->count;
    const int *fromIndex = &from->index[0];
    const double *fromArray = &from->array[0];
    for (int i = 0; i < fromCount; i++) {
        const int iFrom = fromIndex[i];
        const double xFrom = fromArray[iFrom];
        index[i] = iFrom;
        array[iFrom] = xFrom;
    }
}

double HVector::norm2() {
    const int workCount = count;
    const int *workIndex = &index[0];
    const double *workArray = &array[0];

    double result = 0;
    for (int i = 0; i < workCount; i++) {
        double value = workArray[workIndex[i]];
        result += value * value;
    }
    return result;
}

void HVector::saxpy(const double pivotX, const HVector *pivot) {
    int workCount = count;
    int *workIndex = &index[0];
    double *workArray = &array[0];

    const int pivotCount = pivot->count;
    const int *pivotIndex = &pivot->index[0];
    const double *pivotArray = &pivot->array[0];

    for (int k = 0; k < pivotCount; k++) {
        const int iRow = pivotIndex[k];
        const double x0 = workArray[iRow];
        const double x1 = x0 + pivotX * pivotArray[iRow];
        if (x0 == 0)
            workIndex[workCount++] = iRow;
        workArray[iRow] = (fabs(x1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : x1;
    }
    count = workCount;
}

