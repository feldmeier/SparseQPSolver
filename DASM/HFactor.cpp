#include "HFactor.h"
#include "HConst.h"
#include "HTimer.h"

#include "HMatrix.h"

#include <cmath>
#include <cassert>
#include <iostream>
#include <stdexcept>
using namespace std;

vector<int> HFactor::getPermute() {
	return permute;
}

void HFactor::copyFrom(const HFactor *from) {
    numRow = from->numRow;
    numCol = from->numCol;
    Astart = from->Astart;
    Aindex = from->Aindex;
    Avalue = from->Avalue;
    baseIndex = from->baseIndex;
    updateMethod = from->updateMethod;

    nwork = from->nwork;
    iwork = from->iwork;
    dwork = from->dwork;

    LpivotLookup = from->LpivotLookup;
    LpivotIndex = from->LpivotIndex;
    Lstart = from->Lstart;
    Lindex = from->Lindex;
    Lvalue = from->Lvalue;
    LRstart = from->LRstart;
    LRindex = from->LRindex;
    LRvalue = from->LRvalue;

    UpivotLookup = from->UpivotLookup;
    UpivotIndex = from->UpivotIndex;
    UpivotValue = from->UpivotValue;
    UmeritX = from->UmeritX;
    UtotalX = from->UtotalX;
    Ustart = from->Ustart;
    Ulastp = from->Ulastp;
    Uindex = from->Uindex;
    Uvalue = from->Uvalue;
    URstart = from->URstart;
    URlastp = from->URlastp;
    URspace = from->URspace;
    URindex = from->URindex;
    URvalue = from->URvalue;

    PFpivotValue = from->PFpivotValue;
    PFpivotIndex = from->PFpivotIndex;
    PFstart = from->PFstart;
    PFindex = from->PFindex;
    PFvalue = from->PFvalue;
}

void solveMatrixT(const int Xstart, const int Xend, const int Ystart,
        const int Yend, const int *Tindex, const double *Tvalue,
        const double Tpivot, int *RHScount, int *RHSindex, double *RHSarray) {

    // Collect by X
    double pivotX = 0;
    for (int k = Xstart; k < Xend; k++)
        pivotX += Tvalue[k] * RHSarray[Tindex[k]];

    // Scatter by Y
    if (fabs(pivotX) > HSOL_CONST_TINY) {
        int workCount = *RHScount;

        pivotX /= Tpivot;
        for (int k = Ystart; k < Yend; k++) {
            const int index = Tindex[k];
            const double value0 = RHSarray[index];
            const double value1 = value0 - pivotX * Tvalue[k];
            if (value0 == 0)
                RHSindex[workCount++] = index;
            RHSarray[index] =
                    (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
        }

        *RHScount = workCount;
    }
}

void solveHyper(const int Hsize, const int *Hlookup, const int *HpivotIndex,
        const double *HpivotValue, const int *Hstart, const int *Hend,
        const int *Hindex, const double *Hvalue, HVector *rhs) {

    int RHS_Tick = rhs->pseudoTick;
    int RHScount = rhs->count;
    int *RHSindex = &rhs->index[0];
    double *RHSarray = &rhs->array[0];

    // Take count

    // Build list
    char *listMark = &rhs->cwork[0];
    int *listIndex = &rhs->iwork[0];
    int *listStack = &rhs->iwork[Hsize];
    int listCount = 0;

    int countPivot = 0;
    int countEntry = 0;

    for (int i = 0; i < RHScount; i++) {
        // Skip touched index
        int iTrans = Hlookup[RHSindex[i]]; // XXX: this contains a bug iTran
        if (listMark[iTrans]) // XXX bug here
            continue;

        int Hi = iTrans; // H matrix pivot index
        int Hk = Hstart[Hi]; // H matrix non zero position
        int nStack = -1; // Usage of the stack (-1 not used)

        listMark[Hi] = 1; // Mark this as touched

        for (;;) {
            if (Hk < Hend[Hi]) {
                int Hi_sub = Hlookup[Hindex[Hk++]];
                if (listMark[Hi_sub] == 0) { // Go to a child
                    listMark[Hi_sub] = 1; // Mark as touched
                    listStack[++nStack] = Hi; // Store current into stack
                    listStack[++nStack] = Hk;
                    Hi = Hi_sub; // Replace current with child
                    Hk = Hstart[Hi];
                    if (Hi >= Hsize) {
                        countPivot++;
                        countEntry += Hend[Hi] - Hstart[Hi];
                        RHS_Tick += 10;
                    }
                }
            } else {
                listIndex[listCount++] = Hi;
                if (nStack == -1) // Quit on empty stack
                    break;
                Hk = listStack[nStack--]; // Back to last in stack
                Hi = listStack[nStack--];
            }
        }
    }

    rhs->fakeTick += countPivot * 20 + countEntry * 10;

    // Solve with list
    if (HpivotValue == 0) {
        RHScount = 0;
        for (int iList = listCount - 1; iList >= 0; iList--) {
            int i = listIndex[iList];
            listMark[i] = 0;
            int pivotRow = HpivotIndex[i];
            double pivotX = RHSarray[pivotRow];
            if (fabs(pivotX) > HSOL_CONST_TINY) {
                RHSindex[RHScount++] = pivotRow;
                const int start = Hstart[i];
                const int end = Hend[i];
                RHS_Tick += end - start;
                for (int k = start; k < end; k++)
                    RHSarray[Hindex[k]] -= pivotX * Hvalue[k];
            } else
                RHSarray[pivotRow] = 0;
        }
        rhs->pseudoTick = RHS_Tick;
        rhs->count = RHScount;
    } else {
        RHScount = 0;
        for (int iList = listCount - 1; iList >= 0; iList--) {
            int i = listIndex[iList];
            listMark[i] = 0;
            int pivotRow = HpivotIndex[i];
            double pivotX = RHSarray[pivotRow];
            if (fabs(pivotX) > HSOL_CONST_TINY) {
                pivotX /= HpivotValue[i];
                RHSarray[pivotRow] = pivotX;
                RHSindex[RHScount++] = pivotRow;
                const int start = Hstart[i];
                const int end = Hend[i];
                RHS_Tick += end - start;
                for (int k = start; k < end; k++)
                    RHSarray[Hindex[k]] -= pivotX * Hvalue[k];
            } else
                RHSarray[pivotRow] = 0;
        }
        rhs->pseudoTick = RHS_Tick;
        rhs->count = RHScount;
    }
}

void HFactor::setup(int numCol_, int numRow_, int *Astart_, int *Aindex_,
		    double *Avalue_, int *baseIndex_,
		    int updateMethod_) {
  //baseIndex_ used to be passed but isn't used
    // Copy Problem size and (pointer to) coefficient matrix
    numRow = numRow_;
    numCol = numCol_;
    Astart = Astart_;
    Aindex = Aindex_;
    Avalue = Avalue_;
    baseIndex = baseIndex_;
    updateMethod = updateMethod_;

    // Allocate for working buffer
    iwork.reserve(numRow * 2);
    dwork.assign(numRow, 0);

    // Find Basis matrix limit size
    int BlimitX = 0;
    iwork.assign(numRow + 1, 0);
    for (int i = 0; i < numCol; i++)
        iwork[Astart[i + 1] - Astart[i]]++;
    for (int i = numRow, counted = 0; i >= 0 && counted < numRow; i--)
        BlimitX += i * iwork[i], counted += iwork[i];
    BlimitX += numRow;

    // Allocate spaces for basis matrix, L, U factor and Update buffer
    Bstart.resize(numRow + 1, 0);
    Bindex.resize(BlimitX);
    Bvalue.resize(BlimitX);

    // Allocate spaces for pivot records
    permute.resize(numRow);

    // Allocate spaces for Markowitz matrices
    MCstart.resize(numRow);
    MCcountA.resize(numRow);
    MCcountN.resize(numRow);
    MCspace.resize(numRow);
    MCminpivot.resize(numRow);
    MCindex.resize(BlimitX * 2);
    MCvalue.resize(BlimitX * 2);

    MRstart.resize(numRow);
    MRcount.resize(numRow);
    MRspace.resize(numRow);
    MRcountb4.resize(numRow);
    MRindex.resize(BlimitX * 2);

    McolumnMark.assign(numRow, 0);
    McolumnIndex.resize(numRow);
    McolumnArray.assign(numRow, 0);

    // Allocate spaces for count-link-list
    clinkFirst.assign(numRow + 1, -1);
    clinkNext.resize(numRow);
    clinkLast.resize(numRow);

    rlinkFirst.assign(numRow + 1, -1);
    rlinkNext.resize(numRow);
    rlinkLast.resize(numRow);

    // Allocate spaces for L factor
    LpivotLookup.resize(numRow);
    LpivotIndex.reserve(numRow);
    Lstart.reserve(numRow + 1);
    Lindex.reserve(BlimitX * 3);
    Lvalue.reserve(BlimitX * 3);

    LRstart.reserve(numRow + 1);
    LRindex.reserve(BlimitX * 3);
    LRvalue.reserve(BlimitX * 3);

    // Allocate spaces for U factor
    UpivotLookup.resize(numRow);
    UpivotIndex.reserve(numRow + 1000);
    UpivotValue.reserve(numRow + 1000);

    Ustart.reserve(numRow + 1000 + 1);
    Ulastp.reserve(numRow + 1000);
    Uindex.reserve(BlimitX * 3);
    Uvalue.reserve(BlimitX * 3);

    URstart.reserve(numRow + 1000 + 1);
    URlastp.reserve(numRow + 1000);
    URspace.reserve(numRow + 1000);
    URindex.reserve(BlimitX * 3);
    URvalue.reserve(BlimitX * 3);

    // Allocate spaces for Update buffer
    PFpivotValue.reserve(1000);
    PFpivotIndex.reserve(1000);
    PFstart.reserve(2000 + 1);
    PFindex.reserve(BlimitX * 4);
    PFvalue.reserve(BlimitX * 4);
}

void HFactor::change(int updateMethod_) {
    updateMethod = updateMethod_;
}

void HFactor::build() {
    pseudoTick = 0;

    HTimer timer;
    timer.reset();

    // Build the L, U factor
    fakeTick = 0;
    realTick = timer.getTick();
    buildSimple();
    buildKernel();
    buildFinish();

    realTick = timer.getTick() - realTick;

    // I can't believe I made this bug!!!
    pseudoTick *= 8;

//    printf("realTick = %10.0f, COMPARE = %10.4f  COMPAREX = %10.4f\n", realTick,
//            (realTick) / (fakeTick), pseudoTick / fakeTick);
}

void HFactor::ftran(HVector& vector, double hist_dsty) const {
    ftranL(vector, hist_dsty);
    ftranU(vector, hist_dsty);
}

void HFactor::btran(HVector& vector, double hist_dsty) const {
    btranU(vector, hist_dsty);
    btranL(vector, hist_dsty);
}

void HFactor::update(HVector *aq, HVector *ep, int *iRow, int *hint) {
    // Special case
    if (aq->next) {
        updateCFT(aq, ep, iRow, hint);
        return;
    }

    if (updateMethod == UPDATE_METHOD_FT)
        updateFT(aq, ep, *iRow, hint);
    if (updateMethod == UPDATE_METHOD_PF)
        updatePF(aq, ep, *iRow, hint);
    if (updateMethod == UPDATE_METHOD_MPF)
        updateMPF(aq, ep, *iRow, hint);
    if (updateMethod == UPDATE_METHOD_APF)
        updateAPF(aq, ep, *iRow, hint);
}

void HFactor::buildSimple() {
    /**
     * 0. Clear L and U factor
     */
    Lstart.clear();
    Lstart.push_back(0);
    Lindex.clear();
    Lvalue.clear();

    UpivotIndex.clear();
    UpivotValue.clear();
    Ustart.clear();
    Ustart.push_back(0);
    Uindex.clear();
    Uvalue.clear();

    /**
     * 1. Prepare basis matrix and deal with unit columns
     */

    int BcountX = 0;
    fill_n(&MRcountb4[0], numRow, 0);
    nwork = 0;
    for (int iCol = 0; iCol < numRow; iCol++) {
        int iMat = baseIndex[iCol];
        int iRow = -1;

        if (iMat >= numCol) {
            // 1.1 Logical column
            iRow = iMat - numCol;
        } else {
            // 1.2 Structural column
            int start = Astart[iMat];
            int count = Astart[iMat + 1] - start;
            if (count == 1 && Avalue[start] == 1) {
                iRow = Aindex[start];
            } else {
                for (int k = start; k < start + count; k++) {
                    MRcountb4[Aindex[k]]++;
                    Bindex[BcountX] = Aindex[k];
                    Bvalue[BcountX++] = Avalue[k];
                }
                iwork[nwork++] = iCol;
            }
        }

        if (iRow >= 0) {
            // 1.3 Record unit column
            permute[iCol] = iRow;
            Lstart.push_back(Lindex.size());
            UpivotIndex.push_back(iRow);
            UpivotValue.push_back(1);
            Ustart.push_back(Uindex.size());
            MRcountb4[iRow] = -numRow;
        }
        Bstart[iCol + 1] = BcountX;
    }
    BtotalX = numRow - nwork + BcountX;
    pseudoTick += (numRow - nwork) * 2;
    pseudoTick += BcountX * 1.5;

//    count1 = 0;
    // Comments: for pds-20, dfl001: 60 / 80
    // Comments: when system is large: enlarge
    // Comments: when system is small: decrease
    fakeTick += BcountX * 60 + (numRow - nwork) * 80;

    /**
     * 2. Search for and deal with singletons
     */
    double t2_search = 0;
    double t2_storeL = Lindex.size();
    double t2_storeU = Uindex.size();
    double t2_storep = nwork;
    while (nwork > 0) {
        int nworkLast = nwork;
        nwork = 0;
        for (int i = 0; i < nworkLast; i++) {
            const int iCol = iwork[i];
            const int start = Bstart[iCol];
            const int end = Bstart[iCol + 1];
            int pivot_k = -1;
            int found_row_singleton = 0;
            int count = 0;

            // 2.1 Search for singleton
            t2_search += end - start;
            pseudoTick += (end - start) * 0.5;
            for (int k = start; k < end; k++) {
                const int iRow = Bindex[k];
                if (MRcountb4[iRow] == 1) {
                    pivot_k = k;
                    found_row_singleton = 1;
                    break;
                }
                if (MRcountb4[iRow] > 1) {
                    pivot_k = k;
                    count++;
                }
            }

            if (found_row_singleton) {
                // 2.2 Deal with row singleton
                const double pivotX = 1 / Bvalue[pivot_k];
                for (int section = 0; section < 2; section++) {
                    int p0 = section == 0 ? start : pivot_k + 1;
                    int p1 = section == 0 ? pivot_k : end;
                    for (int k = p0; k < p1; k++) {
                        int iRow = Bindex[k];
                        if (MRcountb4[iRow] > 0) {
                            Lindex.push_back(iRow);
                            Lvalue.push_back(Bvalue[k] * pivotX);
                        } else {
                            Uindex.push_back(iRow);
                            Uvalue.push_back(Bvalue[k]);
                        }
                        MRcountb4[iRow]--;
                    }

                }
                int iRow = Bindex[pivot_k];
                MRcountb4[iRow] = 0;
                permute[iCol] = iRow;
                Lstart.push_back(Lindex.size());

                UpivotIndex.push_back(iRow);
                UpivotValue.push_back(Bvalue[pivot_k]);
                Ustart.push_back(Uindex.size());
            } else if (count == 1) {
                // 2.3 Deal with column singleton
                for (int k = start; k < pivot_k; k++) {
                    Uindex.push_back(Bindex[k]);
                    Uvalue.push_back(Bvalue[k]);
                }
                for (int k = pivot_k + 1; k < end; k++) {
                    Uindex.push_back(Bindex[k]);
                    Uvalue.push_back(Bvalue[k]);
                }

                int iRow = Bindex[pivot_k];
                MRcountb4[iRow] = 0;
                permute[iCol] = iRow;
                Lstart.push_back(Lindex.size());

                UpivotIndex.push_back(iRow);
                UpivotValue.push_back(Bvalue[pivot_k]);
                Ustart.push_back(Uindex.size());
            } else {
                iwork[nwork++] = iCol;
            }
        }

        // No singleton found in the last pass
        pseudoTick += (nworkLast - nwork) * 2; // For pivot
        if (nworkLast == nwork)
            break;
    }
    t2_storeL = Lindex.size() - t2_storeL;
    t2_storeU = Uindex.size() - t2_storeU;
    t2_storep = t2_storep - nwork;

    fakeTick += t2_search * 20 + (t2_storep + t2_storeL + t2_storeU) * 80;
    pseudoTick += (Lindex.size() + Uindex.size()) * 1.5; // For factor

    /**
     * 3. Prepare the kernel parts
     */
    // 3.1 Prepare row links, row matrix spaces
    rlinkFirst.assign(numRow + 1, -1);
    MRcount.assign(numRow, 0);
    int MRcountX = 0;
    for (int iRow = 0; iRow < numRow; iRow++) {
        int count = MRcountb4[iRow];
        if (count > 0) {
            MRstart[iRow] = MRcountX;
            MRspace[iRow] = count * 2;
            MRcountX += count * 2;
            rlinkAdd(iRow, count);
        }
    }
    MRindex.resize(MRcountX);

    // 3.2 Prepare column links, kernel matrix
    clinkFirst.assign(numRow + 1, -1);
    MCindex.clear();
    MCvalue.clear();
    MCcountA.assign(numRow, 0);
    MCcountN.assign(numRow, 0);
    int MCcountX = 0;
    for (int i = 0; i < nwork; i++) {
        int iCol = iwork[i];
        MCstart[iCol] = MCcountX;
        MCspace[iCol] = (Bstart[iCol + 1] - Bstart[iCol]) * 2;
        MCcountX += MCspace[iCol];
        MCindex.resize(MCcountX);
        MCvalue.resize(MCcountX);

        for (int k = Bstart[iCol]; k < Bstart[iCol + 1]; k++) {
            const int iRow = Bindex[k];
            const double value = Bvalue[k];
            if (MRcountb4[iRow] > 0) {
                colInsert(iCol, iRow, value);
                rowInsert(iCol, iRow);
            } else {
                colStoreN(iCol, iRow, value);
            }
        }
        colFixMax(iCol);
        clinkAdd(iCol, MCcountA[iCol]);
    }

    pseudoTick += nwork * 2; // DLL
    pseudoTick += MCcountX; // MATRIX

    fakeTick += (numRow + nwork + MCcountX) * 40 + MRcountX * 20;
}

void HFactor::buildKernel() {
    // Deal with the kernel part by 'n-work' pivoting

    double fake_search = 0;
    double fake_fill = 0;
    double fake_eliminate = 0;

    while (nwork-- > 0) {
        /**
         * 1. Search for the pivot
         */
        int jColPivot = -1;
        int iRowPivot = -1;

        // 1.1. Setup search merits
        int searchLimit = min(nwork, 8);
        int searchCount = 0;
        double meritLimit = 1.0 * numRow * numRow;
        double meritPivot = meritLimit;

        // 1.2. Search for local singletons
        bool foundPivot = false;
        if (!foundPivot && clinkFirst[1] != -1) {
            jColPivot = clinkFirst[1];
            iRowPivot = MCindex[MCstart[jColPivot]];
            foundPivot = true;
        }
        if (!foundPivot && rlinkFirst[1] != -1) {
            iRowPivot = rlinkFirst[1];
            jColPivot = MRindex[MRstart[iRowPivot]];
            foundPivot = true;
        }

        // 1.3. Major search loop
        for (int count = 2; !foundPivot && count <= numRow; count++) {
            // 1.3.1 Search for columns
            for (int j = clinkFirst[count]; j != -1; j = clinkNext[j]) {
                double minpivot = MCminpivot[j];
                int start = MCstart[j];
                int end = start + MCcountA[j];
                pseudoTick += MCcountA[j] * 0.3;
                for (int k = start; k < end; k++) {
                    if (fabs(MCvalue[k]) >= minpivot) {
                        int i = MCindex[k];
                        int rowCount = MRcount[i];
                        double meritLocal = 1.0 * (count - 1) * (rowCount - 1);
                        if (meritPivot > meritLocal) {
                            meritPivot = meritLocal;
                            jColPivot = j;
                            iRowPivot = i;
                            foundPivot = foundPivot || (rowCount < count);
                        }
                    }
                }

                if (searchCount++ >= searchLimit && meritPivot < meritLimit)
                    foundPivot = true;
                if (foundPivot)
                    break;

                fake_search += count;
            }

            // 1.3.2 Search for rows
            for (int i = rlinkFirst[count]; i != -1; i = rlinkNext[i]) {
                int start = MRstart[i];
                int end = start + MRcount[i];
                pseudoTick += MRcount[i] * 0.3;
                for (int k = start; k < end; k++) {
                    int j = MRindex[k];
                    int columnCount = MCcountA[j];
                    double meritLocal = 1.0 * (count - 1) * (columnCount - 1);
                    if (meritLocal < meritPivot) {
                        int ifind = MCstart[j];
                        while (MCindex[ifind] != i)
                            ifind++;
                        pseudoTick += columnCount * 0.05;
                        if (fabs(MCvalue[ifind]) >= MCminpivot[j]) {
                            meritPivot = meritLocal;
                            jColPivot = j;
                            iRowPivot = i;
                            foundPivot = foundPivot || (columnCount <= count);
                        }
                    }
                }
                if (searchCount++ >= searchLimit && meritPivot < meritLimit)
                    foundPivot = true;
                if (foundPivot)
                    break;
            }

            fake_search += count;
        }

        // 1.4. If we found nothing: tell singular
        if (!foundPivot)
            throw runtime_error("singular-basis-matrix");

        /**
         * 2. Elimination other elements by the pivot
         */

        // XXX NEED TO UPDATE THIS
        pseudoTick += MCcountA[jColPivot] * 2 + MCcountN[jColPivot];

        // 2.1. Delete the pivot
        double pivotX = colDelete(jColPivot, iRowPivot);
        rowDelete(jColPivot, iRowPivot);
        clinkDel(jColPivot);
        rlinkDel(iRowPivot);
        permute[jColPivot] = iRowPivot;

        // 2.2. Store active pivot column to L
        int start_A = MCstart[jColPivot];
        int end_A = start_A + MCcountA[jColPivot];
        int McolumnCount = 0;
        for (int k = start_A; k < end_A; k++) {
            const int iRow = MCindex[k];
            const double value = MCvalue[k] / pivotX;
            McolumnIndex[McolumnCount++] = iRow;
            McolumnArray[iRow] = value;
            McolumnMark[iRow] = 1;
            Lindex.push_back(iRow);
            Lvalue.push_back(value);
            MRcountb4[iRow] = MRcount[iRow];
            rowDelete(jColPivot, iRow);
        }
        Lstart.push_back(Lindex.size());
        fake_fill += 2 * MCcountA[jColPivot];

        // 2.3. Store non active pivot column to U
        int end_N = start_A + MCspace[jColPivot];
        int start_N = end_N - MCcountN[jColPivot];
        for (int i = start_N; i < end_N; i++) {
            Uindex.push_back(MCindex[i]);
            Uvalue.push_back(MCvalue[i]);
        }
        UpivotIndex.push_back(iRowPivot);
        UpivotValue.push_back(pivotX);
        Ustart.push_back(Uindex.size());
        fake_fill += end_N - start_N;

        // 2.4. Loop over pivot row to eliminate other column
        const int row_start = MRstart[iRowPivot];
        const int row_end = row_start + MRcount[iRowPivot];
        for (int row_k = row_start; row_k < row_end; row_k++) {

            // 2.4.1. My pointer
            int iCol = MRindex[row_k];
            const int my_count = MCcountA[iCol];
            const int my_start = MCstart[iCol];
            const int my_end = my_start + my_count - 1;
            double my_pivot = colDelete(iCol, iRowPivot);
            colStoreN(iCol, iRowPivot, my_pivot);

            // 2.4.2. Elimination on the overlapping part
            int nFillin = McolumnCount;
            int nCancel = 0;
            pseudoTick += McolumnCount;
            for (int my_k = my_start; my_k < my_end; my_k++) {
                int iRow = MCindex[my_k];
                double value = MCvalue[my_k];
                if (McolumnMark[iRow]) {
                    McolumnMark[iRow] = 0;
                    nFillin--;
                    value -= my_pivot * McolumnArray[iRow];
                    if (fabs(value) < HSOL_CONST_TINY) {
                        value = 0;
                        nCancel++;
                    }
                    MCvalue[my_k] = value;
                }
            }
            fake_eliminate += McolumnCount;
            fake_eliminate += nFillin * 2;

            // 2.4.3. Remove cancellation gaps
            if (nCancel > 0) {
                int new_end = my_start;
                for (int my_k = my_start; my_k < my_end; my_k++) {
                    if (MCvalue[my_k] != 0) {
                        MCindex[new_end] = MCindex[my_k];
                        MCvalue[new_end++] = MCvalue[my_k];
                    } else {
                        rowDelete(iCol, MCindex[my_k]);
                    }
                }
                MCcountA[iCol] = new_end - my_start;
            }

            // 2.4.4. Insert fill-in
            if (nFillin > 0) {
                pseudoTick += nFillin * 5;
                // 2.4.4.1 Check column size
                if (MCcountA[iCol] + MCcountN[iCol] + nFillin > MCspace[iCol]) {
                    // p1&2=active, p3&4=non active, p5=new p1, p7=new p3
                    int p1 = MCstart[iCol];
                    int p2 = p1 + MCcountA[iCol];
                    int p3 = p1 + MCspace[iCol] - MCcountN[iCol];
                    int p4 = p1 + MCspace[iCol];
                    MCspace[iCol] += max(MCspace[iCol], nFillin);
                    int p5 = MCstart[iCol] = MCindex.size();
                    int p7 = p5 + MCspace[iCol] - MCcountN[iCol];
                    MCindex.resize(p5 + MCspace[iCol]);
                    MCvalue.resize(p5 + MCspace[iCol]);
                    copy(&MCindex[p1], &MCindex[p2], &MCindex[p5]);
                    copy(&MCvalue[p1], &MCvalue[p2], &MCvalue[p5]);
                    copy(&MCindex[p3], &MCindex[p4], &MCindex[p7]);
                    copy(&MCvalue[p3], &MCvalue[p4], &MCvalue[p7]);
                }

                // 2.4.4.2 Fill into column copy
                for (int i = 0; i < McolumnCount; i++) {
                    int iRow = McolumnIndex[i];
                    if (McolumnMark[iRow])
                        colInsert(iCol, iRow, -my_pivot * McolumnArray[iRow]);
                }

                // 2.4.4.3 Fill into the row copy
                for (int i = 0; i < McolumnCount; i++) {
                    int iRow = McolumnIndex[i];
                    if (McolumnMark[iRow]) {
                        // Expand row space
                        if (MRcount[iRow] == MRspace[iRow]) {
                            int p1 = MRstart[iRow];
                            int p2 = p1 + MRcount[iRow];
                            int p3 = MRstart[iRow] = MRindex.size();
                            MRspace[iRow] *= 2;
                            MRindex.resize(p3 + MRspace[iRow]);
                            copy(&MRindex[p1], &MRindex[p2], &MRindex[p3]);
                        }
                        rowInsert(iCol, iRow);
                    }
                }
            }

            // 2.4.5. Reset pivot column mark
            for (int i = 0; i < McolumnCount; i++)
                McolumnMark[McolumnIndex[i]] = 1;

            // 2.4.6. Fix max value and link list
            colFixMax(iCol);
            if (my_count != MCcountA[iCol]) {
                clinkDel(iCol);
                clinkAdd(iCol, MCcountA[iCol]);
            }
        }

        // 2.5. Clear pivot column buffer
        for (int i = 0; i < McolumnCount; i++)
            McolumnMark[McolumnIndex[i]] = 0;

        // 2.6. Correct row links for the remain active part
        for (int i = start_A; i < end_A; i++) {
            int iRow = MCindex[i];
            if (MRcountb4[iRow] != MRcount[iRow]) {
                rlinkDel(iRow);
                rlinkAdd(iRow, MRcount[iRow]);
            }
        }
    }
    fakeTick += fake_search * 20 + fake_fill * 160 + fake_eliminate * 80;
}

void HFactor::buildFinish() {
    // The look up table
    for (int i = 0; i < numRow; i++)
        UpivotLookup[UpivotIndex[i]] = i;
    LpivotIndex = UpivotIndex;
    LpivotLookup = UpivotLookup;

    // LR space
    int LcountX = Lindex.size();
    LRindex.resize(LcountX);
    LRvalue.resize(LcountX);

    // LR pointer
    iwork.assign(numRow, 0);
    for (int k = 0; k < LcountX; k++)
        iwork[LpivotLookup[Lindex[k]]]++;

    LRstart.assign(numRow + 1, 0);
    for (int i = 1; i <= numRow; i++)
        LRstart[i] = LRstart[i - 1] + iwork[i - 1];

    // LR elements
    iwork.assign(&LRstart[0], &LRstart[numRow]);
    for (int i = 0; i < numRow; i++) {
        const int index = LpivotIndex[i];
        for (int k = Lstart[i]; k < Lstart[i + 1]; k++) {
            int iRow = LpivotLookup[Lindex[k]];
            int iPut = iwork[iRow]++;
            LRindex[iPut] = index;
            LRvalue[iPut] = Lvalue[k];
        }
    }

    // U pointer
    Ulastp.assign(&Ustart[1], &Ustart[numRow + 1]);
    Ustart.resize(numRow);

    // UR space
    int UcountX = Uindex.size();
    int URstuffX = updateMethod == UPDATE_METHOD_FT ? 5 : 0;
    int URcountX = UcountX + URstuffX * numRow;
    URindex.resize(URcountX);
    URvalue.resize(URcountX);

    // UR pointer
    URstart.assign(numRow + 1, 0);
    URlastp.assign(numRow, 0);
    URspace.assign(numRow, URstuffX);
    for (int k = 0; k < UcountX; k++)
        URlastp[UpivotLookup[Uindex[k]]]++;
    for (int i = 1; i <= numRow; i++)
        URstart[i] = URstart[i - 1] + URlastp[i - 1] + URstuffX;
    URstart.resize(numRow);

    // UR element
    URlastp = URstart;
    for (int i = 0; i < numRow; i++) {
        const int index = UpivotIndex[i];
        for (int k = Ustart[i]; k < Ulastp[i]; k++) {
            int iRow = UpivotLookup[Uindex[k]];
            int iPut = URlastp[iRow]++;
            URindex[iPut] = index;
            URvalue[iPut] = Uvalue[k];
        }
    }

    // Re-factor merit
    UmeritX = numRow + (LcountX + UcountX) * 1.5;
    UtotalX = UcountX;
    if (updateMethod == UPDATE_METHOD_PF)
        UmeritX = numRow + UcountX * 4;
    if (updateMethod == UPDATE_METHOD_MPF)
        UmeritX = numRow + UcountX * 3;

    // Clear update buffer
    PFpivotValue.clear();
    PFpivotIndex.clear();
    PFstart.clear();
    PFstart.push_back(0);
    PFindex.clear();
    PFvalue.clear();

    // Finally, permute the base index
    iwork.assign(baseIndex, baseIndex + numRow);
    for (int i = 0; i < numRow; i++)
        baseIndex[permute[i]] = iwork[i];

    pseudoTick += LcountX + UcountX + numRow;
    FtotalX = LcountX + UcountX + numRow;

    fakeTick += numRow * 80 + (LcountX + UcountX) * 60;
}


void HFactor::report() {
    int* nonbasic = new int[numRow];
    for(int i = 0; i < numRow; i++) nonbasic[i] = 1;

    HMatrix lmat;
    lmat.setup(numRow, numRow, &this->Lstart[0],  &this->Lindex[0], &this->Lvalue[0], nonbasic);
    lmat.rp_mtx();

    HMatrix umat;
    umat.setup(numRow, numRow, &this->Ustart[0], &this->Uindex[0], &this->Uvalue[0], nonbasic);
    umat.rp_mtx();

    const int *Ustart = &this->Ustart[0];
    const int *Uend = &this->Ulastp[0];
    const int *Uindex = &this->Uindex[0];
    const double *Uvalue = &this->Uvalue[0];
}

void HFactor::matvecprod(HVector& vector, HVector& result) {
    report();
    const int *Lstart = &this->Lstart[0];
    const int *Lindex = &this->Lindex[0];
    const double *Lvalue = &this->Lvalue[0];
}


void HFactor::ftranL(HVector& rhs, double hist_dsty) const {
    const double hyperFTRANL = 0.15;
    const double hyperCANCEL = 0.05;

    if (updateMethod == UPDATE_METHOD_APF) {
        rhs.tight();
        rhs.pack();
        ftranAPF(rhs);
        rhs.tight();
    }

    double curr_dsty = 1.0 * rhs.count / numRow;
    if (curr_dsty > hyperCANCEL || hist_dsty > hyperFTRANL) {
        // Alias to RHS
        int RHScount = 0;
        int *RHSindex = &rhs.index[0];
        double *RHSarray = &rhs.array[0];

        // Alias to factor L
        const int *Lstart = &this->Lstart[0];
        const int *Lindex = &this->Lindex[0];
        const double *Lvalue = &this->Lvalue[0];

        // Transform
        for (int i = 0; i < numRow; i++) {
            int pivotRow = LpivotIndex[i];
            const double pivotX = RHSarray[pivotRow];
            if (fabs(pivotX) > HSOL_CONST_TINY) {
                RHSindex[RHScount++] = pivotRow;
                const int start = Lstart[i];
                const int end = Lstart[i + 1];
                for (int k = start; k < end; k++)
                    RHSarray[Lindex[k]] -= pivotX * Lvalue[k];
            } else
                RHSarray[pivotRow] = 0;
        }

        // Save the count
        rhs.count = RHScount;
    } else {
        solveHyper(numRow, &LpivotLookup[0], &LpivotIndex[0], 0, &Lstart[0],
                &Lstart[1], &Lindex[0], &Lvalue[0], &rhs);
    }
}

void HFactor::btranL(HVector& rhs, double hist_dsty) const {
    const double hyperBTRANL = 0.10;
    const double hyperCANCEL = 0.05;

    double curr_dsty = 1.0 * rhs.count / numRow;
    if (curr_dsty > hyperCANCEL || hist_dsty > hyperBTRANL) {
        // Alias to RHS
        int RHScount = 0;
        int *RHSindex = &rhs.index[0];
        double *RHSarray = &rhs.array[0];

        // Alias to factor L
        const int *LRstart = &this->LRstart[0];
        const int *LRindex = &this->LRindex[0];
        const double *LRvalue = &this->LRvalue[0];

        // Transform
        for (int i = numRow - 1; i >= 0; i--) {
            int pivotRow = LpivotIndex[i];
            const double pivotX = RHSarray[pivotRow];
            if (fabs(pivotX) > HSOL_CONST_TINY) {
                RHSindex[RHScount++] = pivotRow;
                RHSarray[pivotRow] = pivotX;
                const int start = LRstart[i];
                const int end = LRstart[i + 1];
                for (int k = start; k < end; k++)
                    RHSarray[LRindex[k]] -= pivotX * LRvalue[k];
            } else
                RHSarray[pivotRow] = 0;
        }

        // Save the count
        rhs.count = RHScount;
    } else {
        solveHyper(numRow, &LpivotLookup[0], &LpivotIndex[0], 0, &LRstart[0],
                &LRstart[1], &LRindex[0], &LRvalue[0], &rhs);
    }

    if (updateMethod == UPDATE_METHOD_APF) {
        btranAPF(rhs);
        rhs.tight();
        rhs.pack();
    }
}

void HFactor::ftranU(HVector& rhs, double hist_dsty) const {
    // The update part
    if (updateMethod == UPDATE_METHOD_FT) {
        ftranFT(rhs);
        rhs.tight();
        rhs.packFlag = true;
        rhs.pack();
    }
    if (updateMethod == UPDATE_METHOD_MPF) {
        ftranMPF(rhs);
        rhs.tight();
        rhs.pack();
    }

    // The regular part
    const double hyperFTRANU = 0.10;
    const double hyperCANCEL = 0.05;

    double curr_dsty = 1.0 * rhs.count / numRow;
    if (curr_dsty > hyperCANCEL || hist_dsty > hyperFTRANU) {
        // Alias to non constant
        int RHS_Tick = rhs.pseudoTick;
        int RHScount = 0;
        int *RHSindex = &rhs.index[0];
        double *RHSarray = &rhs.array[0];

        // Alias to the factor
        const int *Ustart = &this->Ustart[0];
        const int *Uend = &this->Ulastp[0];
        const int *Uindex = &this->Uindex[0];
        const double *Uvalue = &this->Uvalue[0];

        // Transform
        double RHS_start = RHS_Tick;
        int UpivotCount = UpivotIndex.size();
        for (int iLogic = UpivotCount - 1; iLogic >= 0; iLogic--) {
            // Skip void
            if (UpivotIndex[iLogic] == -1)
                continue;

            // Normal part
            const int pivotRow = UpivotIndex[iLogic];
            double pivotX = RHSarray[pivotRow];
            if (fabs(pivotX) > HSOL_CONST_TINY) {
                pivotX /= UpivotValue[iLogic];
                RHSindex[RHScount++] = pivotRow;
                RHSarray[pivotRow] = pivotX;
                const int start = Ustart[iLogic];
                const int end = Uend[iLogic];
                if (iLogic >= numRow)
                    RHS_Tick += (end - start);
                for (int k = start; k < end; k++)
                    RHSarray[Uindex[k]] -= pivotX * Uvalue[k];
            } else
                RHSarray[pivotRow] = 0;
        }

        // Save the count
        rhs.count = RHScount;
        rhs.pseudoTick = RHS_Tick + (UpivotCount - numRow);

        rhs.fakeTick += (RHS_Tick - RHS_start) * 15
                + (UpivotCount - numRow) * 10;

    } else {
        solveHyper(numRow, &UpivotLookup[0], &UpivotIndex[0], &UpivotValue[0],
                &Ustart[0], &Ulastp[0], &Uindex[0], &Uvalue[0], &rhs);
    }

    if (updateMethod == UPDATE_METHOD_PF) {
        ftranPF(rhs);
        rhs.tight();
        rhs.pack();
    }

}

void HFactor::btranU(HVector& rhs, double hist_dsty) const {
    if (updateMethod == UPDATE_METHOD_PF) {
        btranPF(rhs);
    }

    // The regular part
    const double hyperBTRANU = 0.15;
    const double hyperCANCEL = 0.05;

    double curr_dsty = 1.0 * rhs.count / numRow;
    if (curr_dsty > hyperCANCEL || hist_dsty > hyperBTRANU) {
        // Alias to non constant
        int RHS_Tick = rhs.pseudoTick;
        int RHScount = 0;
        int *RHSindex = &rhs.index[0];
        double *RHSarray = &rhs.array[0];

        // Alias to the factor
        const int *URstart = &this->URstart[0];
        const int *URend = &this->URlastp[0];
        const int *URindex = &this->URindex[0];
        const double *URvalue = &this->URvalue[0];

        // Transform
        double RHS_start = RHS_Tick;

        int UpivotCount = UpivotIndex.size();
        for (int iLogic = 0; iLogic < UpivotCount; iLogic++) {
            // Skip void
            if (UpivotIndex[iLogic] == -1)
                continue;

            // Normal part
            const int pivotRow = UpivotIndex[iLogic];
            double pivotX = RHSarray[pivotRow];
            if (fabs(pivotX) > HSOL_CONST_TINY) {
                pivotX /= UpivotValue[iLogic];
                RHSindex[RHScount++] = pivotRow;
                RHSarray[pivotRow] = pivotX;
                const int start = URstart[iLogic];
                const int end = URend[iLogic];
                if (iLogic >= numRow)
                    RHS_Tick += end - start;
                for (int k = start; k < end; k++)
                    RHSarray[URindex[k]] -= pivotX * URvalue[k];
            } else
                RHSarray[pivotRow] = 0;
        }

        // Save the count
        rhs.count = RHScount;
        rhs.pseudoTick = RHS_Tick + (UpivotCount - numRow);

        rhs.fakeTick += (RHS_Tick - RHS_start) * 15
                + (UpivotCount - numRow) * 10;

    } else {
        solveHyper(numRow, &UpivotLookup[0], &UpivotIndex[0], &UpivotValue[0],
                &URstart[0], &URlastp[0], &URindex[0], &URvalue[0], &rhs);
    }

    // The update part
    if (updateMethod == UPDATE_METHOD_FT) {
        rhs.tight();
        rhs.pack();
        btranFT(rhs);
        rhs.tight();
    }
    if (updateMethod == UPDATE_METHOD_MPF) {
        rhs.tight();
        rhs.pack();
        btranMPF(rhs);
        rhs.tight();
    }

}

void HFactor::ftranFT(HVector& vector) const {
    // Alias to PF buffer
    const int PFpivotCount = PFpivotIndex.size();
    const int *PFpivotIndex = &this->PFpivotIndex[0];
    const int *PFstart = &this->PFstart[0];
    const int *PFindex = &this->PFindex[0];
    const double *PFvalue = &this->PFvalue[0];

    // Alias to non constant
    int RHS_tick = vector.pseudoTick;
    int RHScount = vector.count;
    int *RHSindex = &vector.index[0];
    double *RHSarray = &vector.array[0];

    // Forwardly apply row ETA
    for (int i = 0; i < PFpivotCount; i++) {
        int iRow = PFpivotIndex[i];
        double value0 = RHSarray[iRow];
        double value1 = value0;
        const int start = PFstart[i];
        const int end = PFstart[i + 1];
        RHS_tick += end - start;
        for (int k = start; k < end; k++)
            value1 -= RHSarray[PFindex[k]] * PFvalue[k];
        // This would skip the situation where they are both zeros
        if (value0 || value1) {
            if (value0 == 0)
                RHSindex[RHScount++] = iRow;
            RHSarray[iRow] =
                    (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
        }
    }

    // Save count back
    vector.count = RHScount;
    vector.pseudoTick = RHS_tick + PFpivotCount * 0.5;

    vector.fakeTick += PFpivotCount * 20 + PFstart[PFpivotCount] * 5;
    if (PFstart[PFpivotCount] / (PFpivotCount + 1) < 5)
        vector.fakeTick += PFstart[PFpivotCount] * 5;

}

void HFactor::btranFT(HVector& vector) const {
    // Alias to PF buffer
    const int PFpivotCount = PFpivotIndex.size();
    const int *PFpivotIndex = &this->PFpivotIndex[0];
    const int *PFstart = &this->PFstart[0];
    const int *PFindex = &this->PFindex[0];
    const double *PFvalue = &this->PFvalue[0];

    // Alias to non constant
    int RHS_tick = vector.pseudoTick;
    int RHScount = vector.count;
    int *RHSindex = &vector.index[0];
    double *RHSarray = &vector.array[0];

    // Backwardly apply row ETA
    double RHS_start = RHS_tick;
    for (int i = PFpivotCount - 1; i >= 0; i--) {
        int pivotRow = PFpivotIndex[i];
        double pivotX = RHSarray[pivotRow];
        if (pivotX) {
            const int start = PFstart[i];
            const int end = PFstart[i + 1];
            RHS_tick += end - start;
            for (int k = start; k < end; k++) {
                int iRow = PFindex[k];
                double value0 = RHSarray[iRow];
                double value1 = value0 - pivotX * PFvalue[k];
                if (value0 == 0)
                    RHSindex[RHScount++] = iRow;
                RHSarray[iRow] =
                        (fabs(value1) < HSOL_CONST_TINY) ?
                                HSOL_CONST_ZERO : value1;
            }
        }
    }

    vector.fakeTick += (RHS_tick - RHS_start) * 15 + PFpivotCount * 10;

    // Save count back
    vector.pseudoTick = RHS_tick + PFpivotCount * 0.3;
    vector.count = RHScount;
}

void HFactor::ftranPF(HVector& vector) const {
    // Alias to PF buffer
    const int PFpivotCount = PFpivotIndex.size();
    const int *PFpivotIndex = &this->PFpivotIndex[0];
    const double *PFpivotValue = &this->PFpivotValue[0];
    const int *PFstart = &this->PFstart[0];
    const int *PFindex = &this->PFindex[0];
    const double *PFvalue = &this->PFvalue[0];

    // Alias to non constant
    int RHScount = vector.count;
    int *RHSindex = &vector.index[0];
    double *RHSarray = &vector.array[0];

    // Forwardly
    for (int i = 0; i < PFpivotCount; i++) {
        int pivotRow = PFpivotIndex[i];
        double pivotX = RHSarray[pivotRow];
        if (fabs(pivotX) > HSOL_CONST_TINY) {
            pivotX /= PFpivotValue[i];
            RHSarray[pivotRow] = pivotX;
            for (int k = PFstart[i]; k < PFstart[i + 1]; k++) {
                const int index = PFindex[k];
                const double value0 = RHSarray[index];
                const double value1 = value0 - pivotX * PFvalue[k];
                if (value0 == 0)
                    RHSindex[RHScount++] = index;
                RHSarray[index] =
                        (fabs(value1) < HSOL_CONST_TINY) ?
                                HSOL_CONST_ZERO : value1;
            }
        }
    }

    // Save count
    vector.count = RHScount;
}

void HFactor::btranPF(HVector& vector) const {
    // Alias to PF buffer
    const int PFpivotCount = PFpivotIndex.size();
    const int *PFpivotIndex = &this->PFpivotIndex[0];
    const double *PFpivotValue = &this->PFpivotValue[0];
    const int *PFstart = &this->PFstart[0];
    const int *PFindex = &this->PFindex[0];
    const double *PFvalue = &this->PFvalue[0];

    // Alias to non constant
    int RHScount = vector.count;
    int *RHSindex = &vector.index[0];
    double *RHSarray = &vector.array[0];

    // Backwardly
    for (int i = PFpivotCount - 1; i >= 0; i--) {
        int pivotRow = PFpivotIndex[i];
        double pivotX = RHSarray[pivotRow];
        for (int k = PFstart[i]; k < PFstart[i + 1]; k++)
            pivotX -= PFvalue[k] * RHSarray[PFindex[k]];
        pivotX /= PFpivotValue[i];

        if (RHSarray[pivotRow] == 0)
            RHSindex[RHScount++] = pivotRow;
        RHSarray[pivotRow] = (fabs(pivotX) < HSOL_CONST_TINY) ? 1e-100 : pivotX;
    }

    // Save count
    vector.count = RHScount;
}

void HFactor::ftranMPF(HVector& vector) const {
    // Alias to non constant
    int RHScount = vector.count;
    int *RHSindex = &vector.index[0];
    double *RHSarray = &vector.array[0];

    // Forwardly
    int PFpivotCount = PFpivotValue.size();
    for (int i = 0; i < PFpivotCount; i++) {
        solveMatrixT(PFstart[i * 2 + 1], PFstart[i * 2 + 2], PFstart[i * 2],
                PFstart[i * 2 + 1], &PFindex[0], &PFvalue[0], PFpivotValue[i],
                &RHScount, RHSindex, RHSarray);
    }

    // Remove cancellation
    vector.count = RHScount;
}

void HFactor::btranMPF(HVector& vector) const {
    // Alias to non constant
    int RHScount = vector.count;
    int *RHSindex = &vector.index[0];
    double *RHSarray = &vector.array[0];

    // Backwardly
    for (int i = PFpivotValue.size() - 1; i >= 0; i--) {
        solveMatrixT(PFstart[i * 2], PFstart[i * 2 + 1], PFstart[i * 2 + 1],
                PFstart[i * 2 + 2], &PFindex[0], &PFvalue[0], PFpivotValue[i],
                &RHScount, RHSindex, RHSarray);
    }

    // Remove cancellation
    vector.count = RHScount;
}

void HFactor::ftranAPF(HVector& vector) const {
    // Alias to non constant
    int RHScount = vector.count;
    int *RHSindex = &vector.index[0];
    double *RHSarray = &vector.array[0];

    // Backwardly
    int PFpivotCount = PFpivotValue.size();
    for (int i = PFpivotCount - 1; i >= 0; i--) {
        solveMatrixT(PFstart[i * 2 + 1], PFstart[i * 2 + 2], PFstart[i * 2],
                PFstart[i * 2 + 1], &PFindex[0], &PFvalue[0], PFpivotValue[i],
                &RHScount, RHSindex, RHSarray);
    }

    // Remove cancellation
    vector.count = RHScount;
}

void HFactor::btranAPF(HVector& vector) const {
    // Alias to non constant
    int RHScount = vector.count;
    int *RHSindex = &vector.index[0];
    double *RHSarray = &vector.array[0];

    // Forwardly
    int PFpivotCount = PFpivotValue.size();
    for (int i = 0; i < PFpivotCount; i++) {
        solveMatrixT(PFstart[i * 2], PFstart[i * 2 + 1], PFstart[i * 2 + 1],
                PFstart[i * 2 + 2], &PFindex[0], &PFvalue[0], PFpivotValue[i],
                &RHScount, RHSindex, RHSarray);
    }
    vector.count = RHScount;
}

void HFactor::updateCFT(HVector *aq, HVector *ep, int *iRow, int *hint) {

    /*
     * In the major update loop, the prefix
     *
     * c(p) = current working pivot
     * p(p) = previous pivot  (0 =< pp < cp)
     */

    int numUpdate = 0;
    for (HVector *vec = aq; vec != 0; vec = vec->next)
        numUpdate++;

    HVector* aqWork[numUpdate];
    HVector* epWork[numUpdate];

    for (int i = 0; i < numUpdate; i++) {
        aqWork[i] = aq;
        epWork[i] = ep;
        aq = aq->next;
        ep = ep->next;
    }

    // Pivot related buffers
    int PFnp0 = PFpivotIndex.size();
    int pLogic[numUpdate];
    double pValue[numUpdate];
    double pAlpha[numUpdate];
    for (int cp = 0; cp < numUpdate; cp++) {
        int cRow = iRow[cp];
        int iLogic = UpivotLookup[cRow];
        pLogic[cp] = iLogic;
        pValue[cp] = UpivotValue[iLogic];
        pAlpha[cp] = aqWork[cp]->array[cRow];
    }

    // Temporary U pointers
    int Tstart[numUpdate + 1];
    double Tpivot[numUpdate];
    Tstart[0] = Uindex.size();

    // Logically sorted previous row_ep
    vector<pair<int, int> > sorted_pp;

    // Major update loop
    for (int cp = 0; cp < numUpdate; cp++) {
        // 1. Expand partial FTRAN result to buffer
        iwork.clear();
        for (int i = 0; i < aqWork[cp]->packCount; i++) {
            int index = aqWork[cp]->packIndex[i];
            double value = aqWork[cp]->packValue[i];
            iwork.push_back(index);
            dwork[index] = value;
        }

        // 2. Update partial FTRAN result by recent FT matrix
        for (int pp = 0; pp < cp; pp++) {
            int pRow = iRow[pp];
            double value = dwork[pRow];
            int PFpp = pp + PFnp0;
            for (int i = PFstart[PFpp]; i < PFstart[PFpp + 1]; i++)
                value -= dwork[PFindex[i]] * PFvalue[i];
            iwork.push_back(pRow); // OK to duplicate
            dwork[pRow] = value;
        }

        // 3. Store the partial FTRAN result to matirx U
        double ppaq = dwork[iRow[cp]]; // pivot of the partial aq
        dwork[iRow[cp]] = 0;
        int UcountX = Tstart[cp];
        int UstartX = UcountX;
        for (unsigned i = 0; i < iwork.size(); i++) {
            int index = iwork[i];
            double value = dwork[index];
            dwork[index] = 0; // This effectively removes all duplication
            if (fabs(value) > HSOL_CONST_TINY) {
                Uindex.push_back(index);
                Uvalue.push_back(value);
            }
        }
        UcountX = Uindex.size();
        Tstart[cp + 1] = UcountX;
        Tpivot[cp] = pValue[cp] * pAlpha[cp];

        // 4. Expand partial BTRAN result to buffer
        iwork.clear();
        for (int i = 0; i < epWork[cp]->packCount; i++) {
            int index = epWork[cp]->packIndex[i];
            double value = epWork[cp]->packValue[i];
            iwork.push_back(index);
            dwork[index] = value;
        }

        // 5. Delete logical later rows (in logical order)
        for (int isort = 0; isort < cp; isort++) {
            int pp = sorted_pp[isort].second;
            int pRow = iRow[pp];
            double multiplier = -pValue[pp] * dwork[pRow];
            if (fabs(dwork[pRow]) > HSOL_CONST_TINY) {
                for (int i = 0; i < epWork[pp]->packCount; i++) {
                    int index = epWork[pp]->packIndex[i];
                    double value = epWork[pp]->packValue[i];
                    iwork.push_back(index);
                    dwork[index] += value * multiplier;
                }
            }
            dwork[pRow] = 0; // Force to be 0
        }

        // 6. Update partial BTRAN result by recent U columns
        for (int pp = 0; pp < cp; pp++) {
            int kpivot = iRow[pp];
            double value = dwork[kpivot];
            for (int k = Tstart[pp]; k < Tstart[pp + 1]; k++)
                value -= dwork[Uindex[k]] * Uvalue[k];
            value /= Tpivot[pp];
            iwork.push_back(kpivot);
            dwork[kpivot] = value; // Again OK to duplicate
        }

        // 6.x compute current alpha
        double thex = 0;
        for (int k = UstartX; k < UcountX; k++) {
            int index = Uindex[k];
            double value = Uvalue[k];
            thex += dwork[index] * value;
        }
        Tpivot[cp] = ppaq + thex * pValue[cp];

        // 7. Store BTRAN result to FT elimination, update logic helper
        dwork[iRow[cp]] = 0;
        double pivotX = -pValue[cp];
        for (unsigned i = 0; i < iwork.size(); i++) {
            int index = iwork[i];
            double value = dwork[index];
            dwork[index] = 0;
            if (fabs(value) > HSOL_CONST_TINY) {
                PFindex.push_back(index);
                PFvalue.push_back(value * pivotX);
            }
        }
        PFpivotIndex.push_back(iRow[cp]);
        FtotalX += PFindex.size() - PFstart.back() + 1;
        UtotalX += PFindex.size() - PFstart.back();
        PFstart.push_back(PFindex.size());

        // 8. Update the sorted ep
        sorted_pp.push_back(make_pair(pLogic[cp], cp));
        sort(sorted_pp.begin(), sorted_pp.end());
    }

    // Now modify the U matrix
    for (int cp = 0; cp < numUpdate; cp++) {
        // 1. Delete pivotal row from U
        int cIndex = iRow[cp];
        int cLogic = pLogic[cp];
        FtotalX -= URlastp[cLogic] - URstart[cLogic];
        UtotalX -= URlastp[cLogic] - URstart[cLogic];
        for (int k = URstart[cLogic]; k < URlastp[cLogic]; k++) {
            // Find the pivotal position
            int iLogic = UpivotLookup[URindex[k]];
            int iFind = Ustart[iLogic];
            int iLast = --Ulastp[iLogic];
            for (; iFind <= iLast; iFind++)
                if (Uindex[iFind] == cIndex)
                    break;
            // Put last to find, and delete last
            Uindex[iFind] = Uindex[iLast];
            Uvalue[iFind] = Uvalue[iLast];
        }

        // 2. Delete pivotal column from UR
        FtotalX -= Ulastp[cLogic] - Ustart[cLogic];
        UtotalX -= Ulastp[cLogic] - Ustart[cLogic];
        for (int k = Ustart[cLogic]; k < Ulastp[cLogic]; k++) {
            // Find the pivotal position
            int iLogic = UpivotLookup[Uindex[k]];
            int iFind = URstart[iLogic];
            int iLast = --URlastp[iLogic];
            for (; iFind <= iLast; iFind++)
                if (URindex[iFind] == cIndex)
                    break;
            // Put last to find, and delete last
            URspace[iLogic]++;
            URindex[iFind] = URindex[iLast];
            URvalue[iFind] = URvalue[iLast];
        }

        // 3. Insert the (stored) partial FTRAN to the row matrix
        int UstartX = Tstart[cp];
        int UendX = Tstart[cp + 1];
        FtotalX += UendX - UstartX + 1;
        UtotalX += UendX - UstartX;
        // Store column as UR elements
        for (int k = UstartX; k < UendX; k++) {
            // Which ETA file
            int iLogic = UpivotLookup[Uindex[k]];

            // Move row to the end if necessary
            if (URspace[iLogic] == 0) {
                // Make pointers
                int row_start = URstart[iLogic];
                int row_count = URlastp[iLogic] - row_start;
                int new_start = URindex.size();
                int new_space = row_count * 1.1 + 5;

                // Check matrix UR
                URindex.resize(new_start + new_space);
                URvalue.resize(new_start + new_space);

                // Move elements
                int iFrom = row_start;
                int iEnd = row_start + row_count;
                int iTo = new_start;
                copy(&URindex[iFrom], &URindex[iEnd], &URindex[iTo]);
                copy(&URvalue[iFrom], &URvalue[iEnd], &URvalue[iTo]);

                // Save new pointers
                URstart[iLogic] = new_start;
                URlastp[iLogic] = new_start + row_count;
                URspace[iLogic] = new_space - row_count;
            }

            // Put into the next available space
            URspace[iLogic]--;
            int iPut = URlastp[iLogic]++;
            URindex[iPut] = cIndex;
            URvalue[iPut] = Uvalue[k];
        }

        // 4. Save pointers
        Ustart.push_back(UstartX);
        Ulastp.push_back(UendX);

        URstart.push_back(URstart[cLogic]);
        URlastp.push_back(URstart[cLogic]);
        URspace.push_back(URspace[cLogic] + URlastp[cLogic] - URstart[cLogic]);

        UpivotLookup[cIndex] = UpivotIndex.size();
        UpivotIndex[cLogic] = -1;
        UpivotIndex.push_back(cIndex);
        UpivotValue.push_back(Tpivot[cp]);
    }

//    // See if we want refactor
//    if (UtotalX > UmeritX && PFpivotIndex.size() > 100)
//        *hint = 1;

}

void HFactor::updateFT(HVector *aq, HVector *ep, int iRow, int *hint) {
    // Store pivot
    int pLogic = UpivotLookup[iRow];
    double pivot = UpivotValue[pLogic];
    double alpha = aq->array[iRow];
    UpivotIndex[pLogic] = -1;

    // Delete pivotal row from U
    FtotalX -= URlastp[pLogic] - URstart[pLogic];
    for (int k = URstart[pLogic]; k < URlastp[pLogic]; k++) {
        // Find the pivotal position
        int iLogic = UpivotLookup[URindex[k]];
        int iFind = Ustart[iLogic];
        int iLast = --Ulastp[iLogic];
        for (; iFind <= iLast; iFind++)
            if (Uindex[iFind] == iRow)
                break;
        // Put last to find, and delete last
        Uindex[iFind] = Uindex[iLast];
        Uvalue[iFind] = Uvalue[iLast];
    }

    // Delete pivotal column from UR
    FtotalX -= Ulastp[pLogic] - Ustart[pLogic];
    for (int k = Ustart[pLogic]; k < Ulastp[pLogic]; k++) {
        // Find the pivotal position
        int iLogic = UpivotLookup[Uindex[k]];
        int iFind = URstart[iLogic];
        int iLast = --URlastp[iLogic];
        for (; iFind <= iLast; iFind++)
            if (URindex[iFind] == iRow)
                break;
        // Put last to find, and delete last
        URspace[iLogic]++;
        URindex[iFind] = URindex[iLast];
        URvalue[iFind] = URvalue[iLast];
    }

    // Store column to U
    Ustart.push_back(Uindex.size());
    for (int i = 0; i < aq->packCount; i++)
        if (aq->packIndex[i] != iRow) {
            Uindex.push_back(aq->packIndex[i]);
            Uvalue.push_back(aq->packValue[i]);
        }
    Ulastp.push_back(Uindex.size());
    int UstartX = Ustart.back();
    int UendX = Ulastp.back();
    UtotalX += UendX - UstartX + 1;
    FtotalX += UendX - UstartX + 1;

    // Store column as UR elements
    for (int k = UstartX; k < UendX; k++) {
        // Which ETA file
        int iLogic = UpivotLookup[Uindex[k]];

        // Move row to the end if necessary
        if (URspace[iLogic] == 0) {
            // Make pointers
            int row_start = URstart[iLogic];
            int row_count = URlastp[iLogic] - row_start;
            int new_start = URindex.size();
            int new_space = row_count * 1.1 + 5;

            // Check matrix UR
            URindex.resize(new_start + new_space);
            URvalue.resize(new_start + new_space);

            // Move elements
            int iFrom = row_start;
            int iEnd = row_start + row_count;
            int iTo = new_start;
            copy(&URindex[iFrom], &URindex[iEnd], &URindex[iTo]);
            copy(&URvalue[iFrom], &URvalue[iEnd], &URvalue[iTo]);

            // Save new pointers
            URstart[iLogic] = new_start;
            URlastp[iLogic] = new_start + row_count;
            URspace[iLogic] = new_space - row_count;
        }

        // Put into the next available space
        URspace[iLogic]--;
        int iPut = URlastp[iLogic]++;
        URindex[iPut] = iRow;
        URvalue[iPut] = Uvalue[k];
    }

    // Store UR pointers
    URstart.push_back(URstart[pLogic]);
    URlastp.push_back(URstart[pLogic]);
    URspace.push_back(URspace[pLogic] + URlastp[pLogic] - URstart[pLogic]);

    // Update pivot count
    UpivotLookup[iRow] = UpivotIndex.size();
    UpivotIndex.push_back(iRow);
    UpivotValue.push_back(pivot * alpha);

    // Store row_ep as R matrix
    for (int i = 0; i < ep->packCount; i++) {
        if (ep->packIndex[i] != iRow) {
            PFindex.push_back(ep->packIndex[i]);
            PFvalue.push_back(-ep->packValue[i] * pivot);
        }
    }
    UtotalX += PFindex.size() - PFstart.back();
    FtotalX += PFindex.size() - PFstart.back() + 1;

    // Store R matrix pivot
    PFpivotIndex.push_back(iRow);
    PFstart.push_back(PFindex.size());

    // Update total countX
    UtotalX -= Ulastp[pLogic] - Ustart[pLogic];
    UtotalX -= URlastp[pLogic] - URstart[pLogic];

//    // See if we want refactor
//    if (UtotalX > UmeritX && PFpivotIndex.size() > 100)
//        *hint = 1;
}

void HFactor::updatePF(HVector *aq, HVector *ep, int iRow, int *hint) {
    // Check space
    const int columnCount = aq->packCount;
    const int *columnIndex = &aq->packIndex[0];
    const double *columnArray = &aq->packValue[0];

    // Copy the pivotal column
    for (int i = 0; i < columnCount; i++) {
        int index = columnIndex[i];
        double value = columnArray[i];
        if (index != iRow) {
            PFindex.push_back(index);
            PFvalue.push_back(value);
        }
    }

    // Save pivot
    PFpivotIndex.push_back(iRow);
    PFpivotValue.push_back(aq->array[iRow]);
    PFstart.push_back(PFindex.size());

    // Check refactor
    UtotalX += aq->packCount;
    if (UtotalX > UmeritX)
        *hint = 1;
    FtotalX += aq->packCount;
}

void HFactor::updateMPF(HVector *aq, HVector *ep, int iRow, int *hint) {
    int PFcountX0 = PFindex.size();

    // Store elements
    for (int i = 0; i < aq->packCount; i++) {
        PFindex.push_back(aq->packIndex[i]);
        PFvalue.push_back(aq->packValue[i]);
    }
    int pLogic = UpivotLookup[iRow];
    int UstartX = Ustart[pLogic];
    int UendX = Ustart[pLogic + 1];
    for (int k = UstartX; k < UendX; k++) {
        PFindex.push_back(Uindex[k]);
        PFvalue.push_back(-Uvalue[k]);
    }
    PFindex.push_back(iRow);
    PFvalue.push_back(-UpivotValue[pLogic]);
    PFstart.push_back(PFindex.size());

    for (int i = 0; i < ep->packCount; i++) {
        PFindex.push_back(ep->packIndex[i]);
        PFvalue.push_back(ep->packValue[i]);
    }
    PFstart.push_back(PFindex.size());

    // Store pivot
    PFpivotValue.push_back(aq->array[iRow]);

    // Refactor or not
    UtotalX += aq->packCount + ep->packCount;
    FtotalX += PFindex.size() - PFcountX0;
    if (UtotalX > UmeritX)
        *hint = 1;
}

void HFactor::updateAPF(HVector *aq, HVector *ep, int iRow, int *hint) {
    int PFcountX0 = PFindex.size();

    // Store elements
    for (int i = 0; i < aq->packCount; i++) {
        PFindex.push_back(aq->packIndex[i]);
        PFvalue.push_back(aq->packValue[i]);
    }

    int columnOut = baseIndex[iRow];
    if (columnOut >= numCol) {
        PFindex.push_back(columnOut - numCol);
        PFvalue.push_back(-1);
    } else {
        for (int k = Astart[columnOut]; k < Astart[columnOut + 1]; k++) {
            PFindex.push_back(Aindex[k]);
            PFvalue.push_back(-Avalue[k]);
        }
    }
    PFstart.push_back(PFindex.size());

    for (int i = 0; i < ep->packCount; i++) {
        PFindex.push_back(ep->packIndex[i]);
        PFvalue.push_back(ep->packValue[i]);
    }
    PFstart.push_back(PFindex.size());
    FtotalX += PFindex.size() - PFcountX0;

    // Store pivot
    PFpivotValue.push_back(aq->array[iRow]);

}

