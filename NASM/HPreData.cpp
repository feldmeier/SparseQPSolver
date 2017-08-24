/*
 * HPreData.cpp
 *
 *  Created on: 8 Jan 2017
 *      Author: ivet
 */

#include "HPreData.h"

HPreData::HPreData() {

}


double HPreData::getRowValue(int i) {
	double sum=0;
	for (int k = ARstart[i]; k < ARstart[i+1]; k++)
	 if (flagRow[ARindex[k]])
		 sum+=ARvalue[k]*valuePrimal[ARindex[k]];
	return sum;
}


void HPreData::printSolution() {
	char buff [10];
    cout<<endl<<"Col value: ";
    for (int i=0;i<numColOriginal;i++) {
    	sprintf(buff, "%2.2f ", valuePrimal[i]);
    	cout<<setw(5)<<buff;
    	if ((i%30) == 0)
    		cout<<std::flush;
    	}

    cout<<endl<<endl;
}

double HPreData::getaij(int i, int j) {
	int k=ARstart[i];
    while (j != ARindex[k] && k<=ARstart[i+1])
            k++;
    if (k==ARstart[i+1]) {
    	//cout<<"Error: No such element in A: row "<<i<<", column "<<j<<endl;
    	//exit(0);
    	}
    return ARvalue[k];
}

bool HPreData::isZeroA(int i, int j) {
	int k=ARstart[i];
    while (j != ARindex[k] && k<ARstart[i+1])
            k++;
    if (k==ARstart[i+1]) {
    	return true;
    	}
    return false;
}


void HPreData::makeARCopy() {
	// Make a AR copy
    int i,j,k;
    vector<int> iwork(numRow, 0);
    ARstart.resize(numRow + 1, 0);
    int AcountX = Aindex.size();
    ARindex.resize(AcountX);
    ARvalue.resize(AcountX);
    for (int k = 0; k < AcountX; k++)
        iwork[Aindex[k]]++;
    for (i = 1; i <= numRow; i++)
        ARstart[i] = ARstart[i - 1] + iwork[i - 1];
    for (i = 0; i < numRow; i++)
        iwork[i] = ARstart[i];
    for (int iCol = 0; iCol < numCol; iCol++) {
        for (k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
            int iRow = Aindex[k];
            int iPut = iwork[iRow]++;
            ARindex[iPut] = iCol;
            ARvalue[iPut] = Avalue[k];
        }
    }
}

void HPreData::makeACopy() {
    // Make a A copy

    int i,j,k;
    vector<int> iwork(numColOriginal, 0);
    Astart.assign(numColOriginal + 1, 0);
    int AcountX = ARindex.size();
    Aindex.resize(AcountX);
    Avalue.resize(AcountX);
    for (int k = 0; k < AcountX; k++)
    	if (ARindex[k] < numColOriginal)
    		iwork[ARindex[k]]++;
    for (i = 1; i <= numColOriginal; i++)
        Astart[i] = Astart[i - 1] + iwork[i - 1];
    for (i = 0; i < numColOriginal; i++)
        iwork[i] = Astart[i];
    for (int iRow = 0; iRow < numRowOriginal; iRow++) {
        for (k = ARstart[iRow]; k < ARstart[iRow + 1]; k++) {
            int iColumn = ARindex[k];
            if (iColumn != numColOriginal) {
				int iPut = iwork[iColumn]++;
				Aindex[iPut] = iRow;
				Avalue[iPut] = ARvalue[k];
            }
        }
    }

    Aend.resize(numColOriginal + 1, 0);
    for (i = 0; i < numColOriginal; i++)
        Aend[i] = Astart[i + 1];
}

void HPreData::print(int k) {
	cout<<"N="<<numCol<<",  M="<<numRow<<",  NZ= "<<Astart[numCol]<<'\n';
	cout<<"\n-----in-------\n";

	char buff [4];
	cout<<"\n-----cost-----\n";

	if (k==0) {
		for (int i=0;i<colCost.size();i++) {
			sprintf(buff, "%2.1g ", colCost[i]);
			cout<<std::setw(5)<<buff;
		}
	}

	if (k==1) {
		for (int i=0;i<colCostAtEl.size();i++) {
			sprintf(buff, "%2.1g ", colCostAtEl[i]);
			cout<<std::setw(5)<<buff;
		}
	}

	if (k==2) {
		for (int i=0;i<colCostAtEl.size();i++) {
			sprintf(buff, "%2.1g ", colCostAtEl[i]);
			cout<<std::setw(5)<<buff;
		}
	}
	cout<<endl;
	cout<<"------A-|-b-----\n";
	int rows= numRow;
	if (k)
		 rows= numRowOriginal;

	for (int i=0;i<rows;i++) {
		if (flagRow[i] ) {
			for (int j=0;j<Astart.size()-1;j++) {

				int ind = Astart[j];
				while (Aindex[ind]!=i && ind<Aend[j])
					ind++;

				if (!flagCol[j])
					continue;

				//if a_ij is nonzero print
				if (Aindex[ind]==i && ind<Aend[j])
				{
					sprintf(buff, "%2.1g ", Avalue[ind]);
					cout<<std::setw(5)<<buff;
					}
				else cout<<std::setw(5)<<"   ";

			}
			cout<<"  |   "<<std::setw(5)<<rowLower[i]<<" < < "<<rowUpper[i]<<endl;
		}
	}
	cout<<"------l------\n";
	for (int i=0;i<colLower.size();i++) {
		if (colLower[i]>-HSOL_CONST_INF)
			sprintf(buff, "%2.1g ", colLower[i]);
		else
			sprintf(buff, "-inf");
		cout<<setw(5)<<buff;
	}
	cout<<endl;
	cout<<"------u------\n";
	for (int i=0;i<colUpper.size();i++) {
		if (colUpper[i]<HSOL_CONST_INF)
			sprintf(buff, "%2.1g ", colUpper[i]);
		else
			sprintf(buff, "inf");
		cout<<setw(5)<<buff;
	}
	cout<<endl;
}

void HPreData::printAR(int i) {
	int rows=numRow, cols=numCol;
	if (i) {
		rows = numRowOriginal;
		cols = numColOriginal;
	}

	cout<<"\n-----cost-----\n";

	char buff [4];
	for (int i=0;i<colCost.size();i++) {
		sprintf(buff, "%2.1g ", colCost[i]);
		cout<<std::setw(5)<<buff;
	}
	cout<<endl;
	cout<<"------AR-|-b-----\n";
	for (int i=0;i<rows;i++) {
		for (int j=0;j<cols;j++) {

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
		cout<<"  |   "<<std::setw(5)<<rowLower[i]<<" < < "<<rowUpper[i]<<endl;
	}
	cout<<"------l------\n";
	for (int i=0;i<cols;i++) {
		if (colLower[i]>-HSOL_CONST_INF)
			sprintf(buff, "%2.1g ", colLower[i]);
		else
			sprintf(buff, "-inf");
		cout<<setw(5)<<buff;
	}
	cout<<endl;
	cout<<"------u------\n";
	for (int i=0;i<cols;i++) {
		if (colUpper[i]<HSOL_CONST_INF)
			sprintf(buff, "%2.1g ", colUpper[i]);
		else
			sprintf(buff, "inf");
		cout<<setw(5)<<buff;
	}
	cout<<endl;
}

void HPreData::writeNewFormat(string fileName) {

}

