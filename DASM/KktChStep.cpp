
#include "KktChStep.h"
#include <utility>
using namespace std;

void KktChStep::passSolution(const vector<double>& colVal, const vector<double>& colDu, const vector<double>& rDu) {

	colValue.resize(RnumCol);
	colDual.resize(RnumCol);
	rowDual.resize(RnumRow);
	
	numRow = rDu.size();
	numCol = colDu.size();
	//arrays to keep track of indices
	vector<int> rIndex(RnumRow, -1);
	vector<int> cIndex(RnumCol, -1);
	int nR = 0;
	int nC = 0;
	
	for (int i=0;i<RnumRow;i++)
		if (flagRow[i]) {
			rIndex[i] = nR;
			nR++;
			}
			
	for (int i=0;i<RnumCol;i++)
		if (flagCol[i]) {
			cIndex[i] = nC;
			nC++;
		}
		
	//set corresponding parts of solution vectors:
	int j=0;
	vector<int> eqIndexOfReduced(numCol, -1);
	vector<int> eqIndexOfReduROW(numRow, -1);
	for (int i=0;i<RnumCol;i++) 
		if (cIndex[i]>-1) {
			eqIndexOfReduced[j] = i;
			j++;
		}
	j=0;
	for (int i=0;i<RnumRow;i++) 
		if (rIndex[i]>-1) {
			eqIndexOfReduROW[j] = i;
			j++;
		}
	
	for (int i=0;i<numCol;i++) {
		colValue[eqIndexOfReduced[i]] = colVal[i];
		colDual[eqIndexOfReduced[i]] = colDu[i];
	}
	
	for (int i=0;i<numRow;i++) 
		rowDual[eqIndexOfReduROW[i]] = rDu[i];
	
}

//full matrix
void KktChStep::setMatrixAR(int nCol, int nRow, const vector<int>& ARstart_, const  vector<int>& ARindex_, const  vector<double>& ARvalue_) {
	RnumCol = nCol;
	RnumRow = nRow;
	ARstart = ARstart_;
	ARindex = ARindex_;
	ARvalue = ARvalue_;
}

void KktChStep::setBoundsCostRHS(const  vector<double>& colUpper_, const  vector<double>& colLower_,const vector<double>& cost, const vector<double>& rowLower_, const vector<double>& rowUpper_) {
	RcolLower = colLower_;
	RcolUpper = colUpper_;
	RrowLower = rowLower_;
	RrowUpper = rowUpper_;
	RcolCost = cost;
}

void KktChStep::addCost(int col, double value) {
	RcolCost[col] = value;
}


void KktChStep::addChange(int type, int row, int col, double valC, double dualC, double dualR) {
//when updating fill new values for b, c, bounds in Rb RcolCost RcolUpper, RcolLower
	vector<pair<int, double>> upd;
	
	switch(type) {
	case 171: //new bounds from doubleton equation, retrieve old ones
		upd =  rLowers.top();
		rLowers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RrowLower[ind] = get<1>(upd[i]);
		}
		upd =  rUppers.top();
		rUppers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RrowUpper[ind] = get<1>(upd[i]);
		}
		break;
	case 172:{ //matrix transformation from doubleton equation, retrieve old value
		//here col is x, the value of which we are changing back.
		//valC is the old value of the x entry in the matrix
		int i = ARstart[row];
		while (i<ARstart[row+1] && ARindex[i]!= col) i++;
		ARvalue[i] = valC;
		break;
	}
	case 173:{ //matrix transformation from doubleton equation, retrieve old value
		//here col is x, the value of which we are changing back.
		//valC is the old value of the x entry in the matrix
		int i = ARstart[row];
		while (i<ARstart[row+1] && ARindex[i]!= col) i++;
		ARvalue[i] = valC;
		ARindex[i] = (int) dualC;
		break;
	}
	case 0:  //empty row
		flagRow[row] = true;
		break;
	case 1:  //row singleton
		flagRow[row] = true;
		colValue[col] = valC;
		colDual[col] = dualC;
		rowDual[row] = dualR;
		if (valC != 0) {
			upd =  cLowers.top();
			cLowers.pop();
			for (int i=0;i<upd.size();i++) {
				int ind = get<0>(upd[i]);
				RcolLower[ind] = get<1>(upd[i]);
			}
			upd =  cUppers.top();
			cUppers.pop();
			for (int i=0;i<upd.size();i++) {
				int ind = get<0>(upd[i]);
				RcolUpper[ind] = get<1>(upd[i]);
			}
		}
		break;
	case 2: //each variable at forcing row: rowDual is cost here
		colValue[col] = valC;
		colDual[col] = dualC;
		flagCol[col] = true;
		RcolCost[col] = dualR;
		break;
	case 21: //
		rowDual[row] = dualR;
		break;
	case 22: //
		upd =  rLowers.top();
		rLowers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RrowLower[ind] = get<1>(upd[i]);
		}
		upd =  rUppers.top();
		rUppers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RrowUpper[ind] = get<1>(upd[i]);
		}
		break;
	case 3: //the row that is forcing 
		rowDual[row] = dualR;
		flagRow[row] = true;
		if (valC != 0) {
			upd =  rLowers.top();
			rLowers.pop();
			for (int i=0;i<upd.size();i++) {
				int ind = get<0>(upd[i]);
				RrowLower[ind] = get<1>(upd[i]);
			}
			upd =  rUppers.top();
			rUppers.pop();
			for (int i=0;i<upd.size();i++) {
				int ind = get<0>(upd[i]);
				RrowUpper[ind] = get<1>(upd[i]);
			}
		}
		break;
	case 4: //implied free column singleton (also from duplicate row)
		flagRow[row] = true;
		flagCol[col] = true;
		colValue[col] = valC;
		colDual[col] = dualC;
		rowDual[row] = dualR;
		upd =  costs.top();
		costs.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RcolCost[ind] = get<1>(upd[i]);
		}
		break;
	case 5: //doubleton eq with singleton col
		flagRow[row] = true;
		flagCol[col] = true;
		colValue[col] = valC;
		colDual[col] = dualC;
		rowDual[row] = dualR;
		upd =  cLowers.top();
		cLowers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RcolLower[ind] = get<1>(upd[i]);
		}
		upd =  cUppers.top();
		cUppers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RcolUpper[ind] = get<1>(upd[i]);
		}
		upd =  costs.top();
		costs.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RcolCost[ind] = get<1>(upd[i]);
		}
		break;
	case 17: {//doubleton equation
		flagRow[row] = true;
		flagCol[col] = true;
		colValue[col] = valC;
		colDual[col] = dualC;
		rowDual[row] = dualR;
		upd =  cLowers.top();
		cLowers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RcolLower[ind] = get<1>(upd[i]);
		}
		upd =  cUppers.top();
		cUppers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RcolUpper[ind] = get<1>(upd[i]);
		}
		upd =  costs.top();
		costs.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RcolCost[ind] = get<1>(upd[i]);
		}
		break;
	}
	case 6: //empty column, dominated column or weakly dominated
		colValue[col] = valC;
		colDual[col] = dualC;
		flagCol[col] = true;
		if (valC != 0) {
			upd =  rLowers.top();
			rLowers.pop();
			for (int i=0;i<upd.size();i++) {
				int ind = get<0>(upd[i]);
				RrowLower[ind] = get<1>(upd[i]);
			}
			upd =  rUppers.top();
			rUppers.pop();
			for (int i=0;i<upd.size();i++) {
				int ind = get<0>(upd[i]);
				RrowUpper[ind] = get<1>(upd[i]);
			}
		}
		break;
	case 7: //fixed variable
		colValue[col] = valC;
		colDual[col] = dualC;
		flagCol[col] = true;
		if (valC != 0) {
			upd =  rLowers.top();
			rLowers.pop();
			for (int i=0;i<upd.size();i++) {
				int ind = get<0>(upd[i]);
				RrowLower[ind] = get<1>(upd[i]);
			}
			upd =  rUppers.top();
			rUppers.pop();
			for (int i=0;i<upd.size();i++) {
				int ind = get<0>(upd[i]);
				RrowUpper[ind] = get<1>(upd[i]);
			}
		}
		break;
	case 11: //empty row from duplucate rows
		upd =  rLowers.top();
		rLowers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RrowLower[ind] = get<1>(upd[i]);
		}
		upd =  rUppers.top();
		rUppers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RrowUpper[ind] = get<1>(upd[i]);
		}
		flagRow[row] = true;
		rowDual[row] = dualR;
		break;
	case 12: //doubleton eq from dupliocate rows;
		flagRow[row] = true;
		flagCol[col] = true;
		colValue[col] = valC;
		colDual[col] = dualC;
		rowDual[row] = dualR;
				upd =  cLowers.top();
		cLowers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RcolLower[ind] = get<1>(upd[i]);
		}
		upd =  cUppers.top();
		cUppers.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RcolUpper[ind] = get<1>(upd[i]);
		}
		upd =  costs.top();
		costs.pop();
		for (int i=0;i<upd.size();i++) {
			int ind = get<0>(upd[i]);
			RcolCost[ind] = get<1>(upd[i]);
		}
		break;
	case 121: //
		colDual[col] = dualC;
		break;/*
	case 14: //two duplicate columns by one 
		colValue[col] = valC;
		colDual[col] = dualC;
		RcolLower = cLowers.top(); cLowers.pop();
		RcolUpper = cUppers.top(); cUppers.pop();
		break;
	case 15: //sing variable on initEq
		flagRow[row] = true;
		rowDual[row] = dualR;
		break;*/
	}
}


void KktChStep::setFlags(vector<bool>& r, vector<bool>& c) {
	flagRow = r;
	flagCol = c;
}

 
void KktChStep::resizeProblemMatrix(KktCheck& checker) {
	//full matrix data in AR copy, A copy made to be passed to checker
	int i, j, k;
	int nz = 0;
	int nR = 0;
	int nC = 0;
	
	//nzRow is the current # on nonzeros, not the general one! 
	vector<int> nzRow(RnumRow, 0);
	for (i=0;i<RnumRow;i++) 
		if (flagRow[i])
			for (k=ARstart[i]; k<ARstart[i+1]; k++) 
				if (flagCol[ARindex[k]]) 
					nzRow[i]++;
	
	numRow = RnumRow;
	numCol = RnumCol;
	//arrays to keep track of indices
	vector<int> rIndex(numRow, -1);
	vector<int> cIndex(numCol, -1);
		 
	for (i=0;i<numRow;i++)
		if (flagRow[i]) {
			nz += nzRow[i];
			rIndex[i] = nR;
			nR++;
			}
			
	for (i=0;i<numCol;i++)
		if (flagCol[i]) {
			cIndex[i] = nC;
			nC++;
		}

	//set indices for checker
	checker.setIndexVectors(rIndex, cIndex);
		
	//counts		
	int numRowEq = numRow;
	int numColEq = numCol;
	numRow = nR;
	numCol = nC;
	
	//matrix
    vector<int> iwork(numCol, 0);
    Astart.resize(numCol + 1, 0);
    Aindex.resize(nz);
    Avalue.resize(nz);
    
        
    for (i = 0;i<numRowEq; i++) 
    	if (flagRow[i])
    	    for (int k = ARstart[i]; k < ARstart[i+1]; k++) {
    	    	j = ARindex[k];
    	    	if (flagCol[j])
        			iwork[cIndex[j]]++;
        		}     
    for (i = 1; i <= numCol; i++)
        Astart[i] = Astart[i - 1] + iwork[i - 1];
    for (i = 0; i < numCol; i++)
        iwork[i] = Astart[i];
    for (i = 0; i < numRowEq; i++) {
    	if (flagRow[i]) {
			int iRow = rIndex[i];
		    for (k = ARstart[i]; k < ARstart[i + 1]; k++) {
		        j = ARindex[k];
		        if (flagCol[j]) {
		        	int iCol = cIndex[j];
				    int iPut = iwork[iCol]++;
				    Aindex[iPut] = iRow;
				    Avalue[iPut] = ARvalue[k];
				}
		    }
		}
    }
	
	//resize cost, RHS, bounds    
    colCost.resize(numCol);
    colLower.resize(numCol);
    colUpper.resize(numCol);
    
    k=0;
    for (i=0;i<numColEq;i++) 
    	if (flagCol[i]) {
    		colCost[k]  = RcolCost[i];
    		colLower[k] = RcolLower[i];
    		colUpper[k] = RcolUpper[i];
    		k++;
	    }
    
    rowLower.resize(numRow);
    rowUpper.resize(numRow);

    k=0;
    for (i=0;i<numRowEq;i++) 
    	if (flagRow[i]) {
    		rowUpper[k] = RrowUpper[i];
    		rowLower[k] = RrowLower[i];

    		//b[k] = Rb[i];
    		k++;
	    }
	
	
		
}

void KktChStep::makeKKTCheck() {
	KktCheck checker;
	resizeProblemMatrix(checker);
	//KktCheck checker;
	checker.setMatrix(Astart, Aindex, Avalue);
	checker.setBounds(colUpper, colLower);
	checker.setNumbersCostRHS(numCol, numRow, rowLower, rowUpper, colCost);
	if (print) 	
		checker.print = print;
	
	//resize and pass solutions
	vector<double> cV;
    vector<double> cD;
    vector<double> rD;
    
    cV.resize(numCol);
    cD.resize(numCol);    
    
    int k=0;
    for (int i=0;i<RnumCol;i++) 
    	if (flagCol[i]) {
    		cV[k]  = colValue[i];
    		cD[k] = colDual[i];
    		k++;
	    }
    
    rD.resize(numRow);
    k=0;
    for (int i=0;i<RnumRow;i++) 
    	if (flagRow[i]) {
    		rD[k] = rowDual[i];
    		k++;
	    }
	
	
	
	checker.passSolution(cV, cD, rD);
	checker.checkKKT();
	}

void KktChStep::printA() {
	char buff [4];
	cout<<"\n-----cost-----\n";

	for (int i=0;i<numCol;i++) { 
		sprintf(buff, "%2.1g ", colCost[i]);
		cout<<std::setw(5)<<buff; 
	}
	cout<<endl;
	cout<<"------A-|-b-----\n";
	for (int i=0;i<numRow;i++) {
		for (int j=0;j<numCol;j++) {

			int ind = Astart[j];
			while (Aindex[ind]!=i && ind<Astart[j+1]) 
				ind++;
			//if a_ij is nonzero print
			if (Aindex[ind]==i && ind<Astart[j+1])
			{	
				sprintf(buff, "%2.1g ", Avalue[ind]);
				cout<<std::setw(5)<<buff; 
				}
			else cout<<std::setw(5)<<"   ";
			
		}
		cout<<"  |   "<<std::setw(5)<<RrowLower[i]<<" < < "<<RrowUpper[i]<<endl;
	}
	cout<<"------l------\n";
	for (int i=0;i<numCol;i++) { 
		if (colLower[i]>-HSOL_CONST_INF)
			sprintf(buff, "%2.1g ", colLower[i]);
		else 
			sprintf(buff, "-inf");
		cout<<setw(5)<<buff; 
	}
	cout<<endl;
	cout<<"------u------\n";
	for (int i=0;i<numCol;i++) { 
		if (colUpper[i]<HSOL_CONST_INF)
			sprintf(buff, "%2.1g ", colUpper[i]);
		else 
			sprintf(buff, "inf");
		cout<<setw(5)<<buff; 
	}
	cout<<endl;

}

void KktChStep::printAR() {
	char buff [4];
	cout<<"\n-----cost-----\n";

	for (int i=0;i<numCol;i++) { 
		sprintf(buff, "%2.1g ", colCost[i]);
		cout<<std::setw(5)<<buff; 
	}
	cout<<endl;
	cout<<"------AR-|-b-----\n";
	for (int i=0;i<RnumRow;i++) {
		for (int j=0;j<RnumCol;j++) {

			int ind = ARstart[i];
			while (ARindex[ind]!=j && ind<ARstart[i+1]) 
				ind++;
			//if a_ij is nonzero print
			if (ARindex[ind]==j && ind<ARstart[i+1])
			{	
				sprintf(buff, "%2.1g ", ARvalue[ind]);
				cout<<std::setw(5)<<buff; 
				}
			else cout<<std::setw(5)<<"   ";
			
		}
		
		cout<<"  |   "<<std::setw(5)<<RrowLower[i]<<" < < "<<RrowUpper[i]<<endl;
		
	}
	
	cout<<endl;
}
	
