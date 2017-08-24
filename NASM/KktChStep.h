#ifndef KKTCHSTEP_H_
#define KKTCHSTEP_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "HConst.h"
#include <cmath>
#include <stack>


#include "KktCheck.h"



class KktChStep {
	//model: full matrix in AR (row-wise) and working copy(column-wise)
	
	int RnumCol;
    int RnumRow;
public:
	vector<int> ARstart;
    vector<int> ARindex;
    vector<double> ARvalue;
    
private:
    //the 4 vectors below always of full length 
    vector<double> RcolCost;
    vector<double> RcolLower;
    vector<double> RcolUpper;
	//vector<double> Rb;
    vector<double> RrowLower;
    vector<double> RrowUpper;

	vector<bool> flagCol;
	vector<bool> flagRow;

	//testing 
	void printA();
	void printAR();
		
	public:
	//data for actual check
	int numCol;
    int numRow;
	vector<int> Astart;
    vector<int> Aindex;
    vector<double> Avalue;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper;
	vector<double> rowLower;
	vector<double> rowUpper;
	int print; 
	
	//solution
	vector<double> colValue;
	vector<double> colDual;  
	vector<double> rowDual;  

	//stack<vector<double> > bs;
	stack<vector<pair<int, double> > > rLowers;
	stack<vector<pair<int, double> > > rUppers;
	stack<vector<pair<int, double> > > cLowers;
	stack<vector<pair<int, double> > > cUppers;
	stack<vector<pair<int, double> > > costs;
	//stack<double> M;
	
	void passSolution(const vector<double>& colVal, const vector<double>& colDu, const vector<double>& rDu);
	//full matrix
	void setMatrixAR(int nCol, int nRow, const vector<int>& ARstart_, const  vector<int>& ARindex_, const  vector<double>& ARvalue_);
	void setBoundsCostRHS(const  vector<double>& colUpper_, const  vector<double>& colLower_,const vector<double>& cost, const vector<double>& rowLower_, const vector<double>& rowUpper_);
	void addChange(int type, int row, int col, double valC, double dualC, double dualR);
	void setFlags(vector<bool>& r, vector<bool>& c);
	void makeKKTCheck();
	void resizeProblemMatrix(KktCheck& checker);
	void addCost(int col, double value);
	
};
#endif /* KKTCHSTEP_H_ */
