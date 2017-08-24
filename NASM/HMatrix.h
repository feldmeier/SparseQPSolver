#ifndef HMATRIX_H_
#define HMATRIX_H_

#include "HVector.h"

#include <vector>
using namespace std;

class HMatrix {
public:
	const int *getAstart() const {
		return &Astart[0];
	}
	const int *getAindex() const {
	    return &Aindex[0];
	}
	const double *getAvalue() const {
	    return &Avalue[0];
	}
	const int *getARstart() const {
	    return &ARstart[0];
	}
	const int *getARindex() const {
	    return &ARindex[0];
	}
	const double *getARvalue() const {
	    return &ARvalue[0];
	}
    void setup(int numCol, int numRow, const int *Astart, const int *Aindex,
            const double *Avalue, const int *nonbasicFlag);
    void setup_lgBs(int numCol, int numRow, const int *Astart, const int *Aindex,
            const double *Avalue);
    bool setup_ok(const int *nonbasicFlag);
    void price_by_col(HVector& row_ap, HVector& row_ep) const;
    void price_by_row(HVector& row_ap, HVector& row_ep) const;
    void update(int columnIn, int columnOut);
    double compute_dot(HVector& vector, int iCol) const;
    void collect_aj(HVector& vector, int iCol, double multi) const;

    void compute_vecT_matB(const double *vec, const int *base, HVector *res);
    void compute_matB_vec(const double *vec, const int *base, HVector *res);
    void compute_mat_mat(HMatrix& mat, HMatrix& res);
    void compute_vec_mat(HVector& vec, HVector& res);
    void rp_mtx();
    int numCol;
    int numRow;
private:
    vector<int> Astart;
    vector<int> Aindex;
    vector<double> Avalue;

    vector<int> ARstart;
    vector<int> AR_Nend;
    vector<int> ARindex;
    vector<double> ARvalue;

};

#endif /* HMATRIX_H_ */
