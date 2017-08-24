/*
 * HCrash.cpp
 *
 * Created on: 20 Oct 2016
 * Author: Julian
 */

#include "HCrash.h"
//#include "HModel.h"
#include "HMatrix.h"
//TODO Why doesn't including HDual.h acquire Presolve_Mode?
#include "HDual.h"
//#include "HFactor.h"
//#include "HVector.h"
//#include "HDualRow.h"
//#include "HDualRHS.h"

#include <set>
#include <string>
#include <vector>
#include <cassert>
using namespace std;

void HCrash::crash(HModel *ptr_model, int Crash_Mode) {
  model = ptr_model;
  if (model->getNumRow() == 0) return;
  model->timer.reset();
  numRow = model->getNumRow();
  numCol = model->getNumCol();
  numTot = model->getNumTot();
  
  bool HaveImpliedBounds = model->impliedBoundsPresolve;
  //bool HaveImpliedBounds = Presolve_Mode == Presolve_Mode_On;
  //	printf("HaveImpliedBounds = %d\n", HaveImpliedBounds);
#ifdef JAJH_dev
  if (HaveImpliedBounds) {crsh_ck_an_impl_bd();}
#endif
  //      Indicate that we are not currently using the implied bounds from presolve
  model->usingImpliedBoundsPresolve = false;
  if (HaveImpliedBounds && Crash_Mode == Crash_Mode_Dev) {
    //		Overwrite the true bounds by the implied bounds
#ifdef JAJH_dev
    printf("Over-writing bounds with implied bounds\n");
#else
    printf("Over-writing bounds with implied bounds\n");
#endif
    model->copy_impliedBoundsToModelBounds();
  }
  
  if (Crash_Mode == Crash_Mode_Bixby || Crash_Mode == Crash_Mode_Dev) {
    bixby(ptr_model, Crash_Mode);
  } else {
    ltssf(ptr_model, Crash_Mode);
  }
  model->totalTime += model->timer.getTime();
}

void HCrash::bixby(HModel *ptr_model, int Crash_Mode) {
  model = ptr_model;
  const int *Astart = &model->Astart[0];
  const int *Aindex = &model->Aindex[0];
  const double *Avalue = &model->Avalue[0];;
  int numEl = Astart[numCol];
  
  bixby_no_nz_c_co = Crash_Mode == Crash_Mode_Dev;
  bixby_no_nz_c_co = false;
  
  bool perform_crash = bixby_iz_da(ptr_model);
  if (!perform_crash) return;
  
  //	bixby_rp_mrt(ptr_model);
  
  //	These multipliers are in Step 2(a) and Step 2(b) of the paper: default values 0.99 and 0.01
  bixby_mu_a = 0.99;
  bixby_mu_b = 0.01;
  double rp_v;
  
#ifdef JAJH_dev
  printf("\nBixby Crash");
  if (bixby_no_nz_c_co)
    printf(": No basic columns with nonzero costs\n");
  else
    printf(": Any basic columns regardless of cost\n");
#endif
  for (int ps_n = 0; ps_n < numCol; ps_n++) {
    //		In each pass:
    //		Consider column c_n
    int c_n = bixby_mrt_ix[ps_n];
    double c_mx_abs_v = crsh_mtx_c_mx_abs_v[c_n];
    //		Find the max |entry| over all rows with a count of zero
    double aa = 0;
    int r_o_mx_aa;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      if (bixby_r_k[r_n] == 0) {
	double lc_aa = abs(Avalue[el_n]);
	if (lc_aa > aa) {
	  aa = lc_aa;
	  r_o_mx_aa = r_n;
	}
      }
    }
#ifdef JAJH_dev
    if (Rp_Bixby_Ps) printf("Pass %3d: c_n = %3d; MxEn = %10g; aa = %10g; aa/MxEn = %10g",
			    ps_n, c_n, c_mx_abs_v, aa, aa / c_mx_abs_v);
#endif
    //		Scale aa by the max |entry| in the column since CPLEX assumes the
    //		matrix is scaled so the max entry in each column is 1
    aa /= c_mx_abs_v;
    bool nx_ps = false;
    if (aa >= bixby_mu_a) {
      //			Column pv_c_n becomes basic in row pv_r_n
      int pv_c_n = c_n;
      int pv_r_n = r_o_mx_aa;
      //			printf(" ** Type a: c_n = %d; pv_c_n = %d ** ", c_n, pv_c_n);
      bixby_pv_in_r[pv_r_n] = 1;
      bixby_vr_in_r[pv_r_n] = pv_c_n;
      bixby_pseudo_pv_v[pv_r_n] = aa;
      for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
	//				int r_n = Aindex[el_n];
	//				printf("\n Row %3d: value %g", r_n, Avalue[el_n]);
	bixby_r_k[Aindex[el_n]] += 1;
      }
      bixby_n_cdd_r -= 1;
#ifdef JAJH_dev
      if (Rp_Bixby_Ps) printf(": pv_r = %3d; n_cdd_r = %3d\n", pv_r_n, bixby_n_cdd_r);
#endif
      nx_ps = true;
    } else {
      //			Find out if there is some row l for which |entry| > bixby_mu_b * v_l
      for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
	int r_n = Aindex[el_n];
	//				If this value in the column would give an unacceptable multiplier then continue to next pass
	nx_ps = abs(Avalue[el_n])
	  > bixby_mu_b * bixby_pseudo_pv_v[r_n] * c_mx_abs_v;
	if (nx_ps) {
	  rp_v = abs(Avalue[el_n])/(bixby_pseudo_pv_v[r_n] * c_mx_abs_v);
	  break;
	}
      }
#ifdef JAJH_dev
      if (nx_ps) if (Rp_Bixby_Ps) printf(": Unacceptable multiplier of %g > %g\n", rp_v, bixby_mu_b);
#endif
    }
    //		Some value in the column would give an unacceptable multiplier so continue to next pass
    if (nx_ps) continue;
    //		Find out whether there is an entry in a row with no pivot
    aa = 0;
    r_o_mx_aa = no_ix;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      if (bixby_pv_in_r[r_n] == 0) {
	double lc_aa = abs(Avalue[el_n]);
	if (lc_aa > aa) {
	  aa = lc_aa;
	  r_o_mx_aa = r_n;
	}
      }
    }
    //		If there is no entry in a row with no pivot then continue to next pass
    if (r_o_mx_aa == no_ix) {
      if (Rp_Bixby_Ps) printf(": No entry in a row with no pivot\n");
      continue;}
    //		Scale aa by the max |entry| in the column since CPLEX assumes the
    //		matrix is scaled so the max entry in each column is 1
    aa /= c_mx_abs_v;
    //		Column pv_c_n becomes basic in row pv_r_n
    int pv_c_n = c_n;
    int pv_r_n = r_o_mx_aa;
    //		printf(" ** Type b: c_n = %d; pv_c_n = %d ** ", c_n, pv_c_n);
    bixby_pv_in_r[pv_r_n] = 1;
    bixby_vr_in_r[pv_r_n] = pv_c_n;
    bixby_pseudo_pv_v[pv_r_n] = aa;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      //			int r_n = Aindex[el_n];
      //			printf("\n Row %3d: value %g", r_n, Avalue[el_n]);
      bixby_r_k[Aindex[el_n]] += 1;
    }
    bixby_n_cdd_r -= 1;
#ifdef JAJH_dev
    if (Rp_Bixby_Ps) printf(": pv_r = %3d; n_cdd_r = %3d\n", pv_r_n, bixby_n_cdd_r);
#endif
    if (bixby_n_cdd_r == 0) break;
  }
  for (int r_n = 0; r_n < numRow; r_n++) {
    //		printf("Bixby vr in row %3d is %3d: Qy %d\n", r_n, bixby_vr_in_r[r_n], numCol+r_n);
    if (bixby_vr_in_r[r_n] == no_ix) continue;
    if (bixby_vr_in_r[r_n] == numCol+r_n) continue;
    int cz_r_n = r_n;
    int cz_c_n = bixby_vr_in_r[r_n];
    int columnIn = cz_c_n;
    int rowOut = cz_r_n;
    int columnOut = numCol + r_n;
    int sourceOut = 1;
    //			Update the basic/nonbasic variable info and the row-wise copy of the matrix
    model->updatePivots(columnIn, rowOut, sourceOut);
    if (model->mlFg_haveMatrixRowWise)
      model->updateMatrix(columnIn, columnOut);
    int vr_ty = crsh_r_ty[cz_r_n];
    crsh_vr_ty_rm_n_r[vr_ty] += 1;
    vr_ty = crsh_c_ty[cz_c_n];
    crsh_vr_ty_add_n_c[vr_ty] += 1;
  }
  
  crsh_an_r_c_st_af(ptr_model);
}

void HCrash::bixby_rp_mrt(HModel *ptr_model) {
  model = ptr_model; 
  const double *colCost = model->getcolCost();
  const double *colLower = model->getcolLower();
  const double *colUpper = model->getcolUpper();
  double mx_co_v = -HSOL_CONST_INF;
  for (int c_n = 0; c_n < numCol; c_n++)
    mx_co_v = max(abs(colCost[c_n]), mx_co_v);
  double co_v_mu = 1;
  if (mx_co_v > 0) co_v_mu = 1e3 * mx_co_v;
  int c_n = bixby_mrt_ix[0];
  double mrt_v = bixby_mrt_v[0];
  double prev_mrt_v0 = -HSOL_CONST_INF;
  double prev_mrt_v = -HSOL_CONST_INF;
  bool rp_c;
  bool rp_al_c = false;
  int n_mrt_v = 0;
  if (mx_co_v > 0) co_v_mu = 1e3 * mx_co_v;
  printf("\nAnalysis of sorted Bixby merits\n");
  for (int ps_n = 0; ps_n < numCol; ps_n++) {
    double mrt_v = bixby_mrt_v[ps_n];
    int c_n = bixby_mrt_ix[ps_n];
    double mrt_v0 = mrt_v - colCost[c_n] / co_v_mu;
    double c_lb = colLower[c_n];
    double c_ub = colUpper[c_n];
    if ((ps_n == 0) || (ps_n == numCol-1))
      rp_c = true;
    else if (
	     (crsh_c_ty[c_n] != crsh_c_ty[bixby_mrt_ix[ps_n-1]]) ||
	     (crsh_c_ty[c_n] != crsh_c_ty[bixby_mrt_ix[ps_n+1]])) {
      rp_c = true;
      prev_mrt_v = -HSOL_CONST_INF;
      prev_mrt_v0 = -HSOL_CONST_INF;
    }
    else if (rp_al_c) 
      rp_c = true;
    else
      rp_c = mrt_v0 > prev_mrt_v0;
    prev_mrt_v0 = mrt_v0;
    if (mrt_v > prev_mrt_v) {
      n_mrt_v += 1;
      prev_mrt_v = mrt_v;
    }
    if (rp_c) printf("%5d: Col %5d, Type = %1d; MrtV = %10.4g; MrtV0 = %10.4g; [%10.4g,%10.4g]\n", ps_n, c_n, crsh_c_ty[c_n], mrt_v, mrt_v0, c_lb, c_ub);
  }
  printf("\n%6d different Bixby merits\n", n_mrt_v);
}

bool HCrash::bixby_iz_da(HModel *ptr_model) {
  model = ptr_model;
  const int *Astart = &model->Astart[0];
  const int *Aindex = &model->Aindex[0];
  const double *Avalue = &model->Avalue[0];;
  int numEl = Astart[numCol];
  
  const double *colCost = model->getcolCost();
  const double *colLower = model->getcolLower();
  const double *colUpper = model->getcolUpper();
  
  //	const double *primalColLowerImplied = model->getprimalColLowerImplied();
  //	const double *primalColUpperImplied = model->getprimalColUpperImplied();
  //
  //	const double *dualColLowerImplied = model->getdualColLowerImplied();
  //	const double *dualColUpperImplied = model->getdualColUpperImplied();
  
  //Allocate the arrays required for crash
  crsh_mtx_c_mx_abs_v.resize(numCol);
  
  bixby_mrt_v.resize(numCol + 1);
  bixby_pseudo_pv_v.resize(numRow);
  bixby_mrt_ix.resize(numCol + 1);
  bixby_pv_in_r.resize(numRow);
  bixby_vr_in_r.resize(numRow);
  bixby_r_k.resize(numRow);
  //	bixby_ze_r_k.resize(numRow);
  
#ifdef JAJH_dev
  crsh_an_c_co(ptr_model);
#endif
  crsh_iz_vr_ty(ptr_model);
  
#ifdef JAJH_dev
  crsh_rp_r_c_st(0);
#endif
  
  //	Initialise the arrays required for the Bixby crash
  //	bixby_pseudo_pv_v: 	"v" in the paper: this is a pseudo pivot value for each row
  //	bixby_pv_in_r: 		"I" in the paper: this is a 0/1 flag to indicate whether there is a basic variable in each row
  //	bixby_vr_in_r:		"B" is the paper: this is the basic variable in a particular row
  //	bixby_r_k:			"r" in the paper: this is the number of entries in each row of the basis matrix
  bixby_n_cdd_r = numRow;
  for (int r_n = 0; r_n < numRow; r_n++) {
    bixby_pseudo_pv_v[r_n] = HSOL_CONST_INF;
    if (crsh_r_ty[r_n] == crsh_vr_ty_fx) {
      bixby_pv_in_r[r_n] = 0;
      bixby_vr_in_r[r_n] = no_ix;
      bixby_r_k[r_n] = 0;
      //			bixby_ze_r_k[r_n] = 1;
    } else {
      bixby_pv_in_r[r_n] = 1;
      bixby_vr_in_r[r_n] = numCol + r_n;
      bixby_r_k[r_n] = 1;
      //			bixby_ze_r_k[r_n] = 0;
      bixby_n_cdd_r -= 1;
    }
  }
  if (bixby_n_cdd_r == 0) return false;
  double mx_co_v = -HSOL_CONST_INF;
  for (int c_n = 0; c_n < numCol; c_n++) {
    //	Find largest |entry| in each column
    crsh_mtx_c_mx_abs_v[c_n] = 0.0;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      crsh_mtx_c_mx_abs_v[c_n] = max(abs(Avalue[el_n]),
				     crsh_mtx_c_mx_abs_v[c_n]);
    }
    //		printf("Column %3d has max |entry| = %g\n", c_n, crsh_mtx_c_mx_abs_v[c_n]);
    mx_co_v = max(abs(colCost[c_n]), mx_co_v);
  }
  double co_v_mu = 1;
  if (mx_co_v > 0) co_v_mu = 1e3 * mx_co_v;
  //	... and then updated with c_j/c_max, being colCost[c_n]/co_v_mu
  //	So, first compute the cost coefficient of maximum absolute value
  int os;
  int n_en;
  os = 0;
  //Free columns - impossible after presolve
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_fr) continue;
    if (bixby_no_nz_c_co && colCost[c_n] != 0.0) continue;
    n_en += 1;
    bixby_mrt_v[os+n_en] = colCost[c_n] / co_v_mu;
    bixby_mrt_ix[os+n_en] = c_n;
  }
  if (n_en >0) {
    build_maxheap(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    heapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os+1; en <= os+n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en+1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en+1];
    }
    os += 1+n_en;
  }
#ifdef JAJH_dev
  printf("%5d Free    cols (%1d)\n", n_en, crsh_vr_ty_fr);
  //	printf("%5d Free    cols (%1d): from %5d to %5d\n", n_en, crsh_vr_ty_fr, os-n_en, os-1);
#endif
  
  //1-sided columns
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_1_sd) continue;
    if (bixby_no_nz_c_co && colCost[c_n] != 0.0) continue;
    n_en += 1;
    if (colUpper[c_n] >= HSOL_CONST_INF) {
      bixby_mrt_v[os+n_en] =  colLower[c_n] + colCost[c_n] / co_v_mu;
    } else {
      bixby_mrt_v[os+n_en] = -colUpper[c_n] + colCost[c_n] / co_v_mu;
    }
    bixby_mrt_ix[os+n_en] = c_n;
  }
  if (n_en >0) {
    build_maxheap(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    heapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os+1; en <= os+n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en+1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en+1];
    }
    os += 1+n_en;
  }
#ifdef JAJH_dev
  printf("%5d 1-sided cols (%1d)\n", n_en, crsh_vr_ty_1_sd);
  //	printf("%5d 1-sided cols (%1d): from %5d to %5d\n", n_en, crsh_vr_ty_1_sd, os-n_en, os-1);
#endif
  
  //2-sided columns
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_2_sd) continue;
    if (bixby_no_nz_c_co && colCost[c_n] != 0.0) continue;
    n_en += 1;
    bixby_mrt_v[os+n_en] =  colLower[c_n]-colUpper[c_n] + colCost[c_n] / co_v_mu;
    bixby_mrt_ix[os+n_en] = c_n;
  }
  if (n_en >0) {
    build_maxheap(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    heapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os+1; en <= os+n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en+1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en+1];
    }
    os += 1+n_en;
  }
#ifdef JAJH_dev
  printf("%5d 2-sided cols (%1d)\n", n_en, crsh_vr_ty_2_sd);
  //	printf("%5d 2-sided cols (%1d): from %5d to %5d\n", n_en, crsh_vr_ty_2_sd, os-n_en, os-1);
#endif
  
  //Fixed columns - impossible after presolve
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_fx) continue;
    if (bixby_no_nz_c_co && colCost[c_n] != 0.0) continue;
    n_en += 1;
    bixby_mrt_v[os+n_en] = HSOL_CONST_INF;
    bixby_mrt_ix[os+n_en] = c_n;
  }
  if (n_en >0) {
    build_maxheap(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    heapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os+1; en <= os+n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en+1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en+1];
    }
    os += 1+n_en;
  }
#ifdef JAJH_dev
  printf("%5d Fixed   cols (%1d)\n", n_en, crsh_vr_ty_fx);
  //	printf("%5d Fixed   cols (%1d): from %5d to %5d\n", n_en, crsh_vr_ty_fx, os-n_en, os-1);
  
  //	printf("Last Bixby index touched is %d of %d\n", os, numCol);
#endif
  return true;
}

void HCrash::crsh_iz_vr_ty(HModel *ptr_model) {
  model = ptr_model;
  const double *colLower = model->getcolLower();
  const double *colUpper = model->getcolUpper();
  const double *rowLower = model->getrowLower();
  const double *rowUpper = model->getrowUpper();
  //Allocate the arrays required for crash
  crsh_r_ty.resize(numRow);
  crsh_c_ty.resize(numCol);
  for (int r_n = 0; r_n < numRow; r_n++) {
    if (rowUpper[r_n] >= HSOL_CONST_INF) {
      if (rowLower[r_n] <= -HSOL_CONST_INF)
	crsh_r_ty[r_n] = crsh_vr_ty_fr;    // Free row
      else
	crsh_r_ty[r_n] = crsh_vr_ty_1_sd; // Lower-bounded (1-sided) row
    } else {
      if (rowLower[r_n] <= -HSOL_CONST_INF)
	crsh_r_ty[r_n] = crsh_vr_ty_1_sd;  // Upper-bonded (1-sided) row
      else {
	//Two-sided row - maybe fixed (equality)
	if (rowLower[r_n] != rowUpper[r_n])
	  crsh_r_ty[r_n] = crsh_vr_ty_2_sd;    // 2-sided row
	else
	  crsh_r_ty[r_n] = crsh_vr_ty_fx;    // fixed (equality) row
      }
    }
  }
  // Set up the column variable types for crash
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (colUpper[c_n] >= HSOL_CONST_INF) {
      if (colLower[c_n] <= -HSOL_CONST_INF)
	crsh_c_ty[c_n] = crsh_vr_ty_fr;    // Free column
      else
	crsh_c_ty[c_n] = crsh_vr_ty_1_sd; // Lower-bounded (1-sided) column
    } else {
      if (colLower[c_n] <= -HSOL_CONST_INF)
	crsh_c_ty[c_n] = crsh_vr_ty_1_sd;  // Upper-bonded (1-sided) column
      else {
	//Two-sided row - maybe fixed (equality)
	if (colLower[c_n] != colUpper[c_n])
	  crsh_c_ty[c_n] = crsh_vr_ty_2_sd;    // 2-sided column
	else
	  crsh_c_ty[c_n] = crsh_vr_ty_fx;    // fixed column
      }
    }
  }
  //Allocate the arrays to analyse crash
  crsh_vr_ty_og_n_r.resize(crsh_l_vr_ty + 1);
  crsh_vr_ty_rm_n_r.resize(crsh_l_vr_ty + 1);
  crsh_vr_ty_og_n_c.resize(crsh_l_vr_ty + 1);
  crsh_vr_ty_add_n_c.resize(crsh_l_vr_ty + 1);
  crsh_bs_vr_ty_n_r.resize(crsh_l_vr_ty + 1);
  crsh_bs_vr_ty_n_c.resize(crsh_l_vr_ty + 1);
  crsh_nonbc_vr_ty_n_r.resize(crsh_l_vr_ty + 1);
  crsh_nonbc_vr_ty_n_c.resize(crsh_l_vr_ty + 1);
  //Initialise the counts of numbers and changes of variable types - just for reporting
  for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_l_vr_ty + 1; vr_ty++) {
    crsh_vr_ty_og_n_r[vr_ty] = 0;
    crsh_vr_ty_og_n_c[vr_ty] = 0;
    crsh_vr_ty_rm_n_r[vr_ty] = 0;
    crsh_vr_ty_add_n_c[vr_ty] = 0;
    crsh_bs_vr_ty_n_r[vr_ty] = 0;
    crsh_bs_vr_ty_n_c[vr_ty] = 0;
    crsh_nonbc_vr_ty_n_r[vr_ty] = 0;
    crsh_nonbc_vr_ty_n_c[vr_ty] = 0;
  }
  for (int r_n = 0; r_n < numRow; r_n++) crsh_vr_ty_og_n_r[crsh_r_ty[r_n]] += 1;
  for (int c_n = 0; c_n < numCol; c_n++) crsh_vr_ty_og_n_c[crsh_c_ty[c_n]] += 1;
}

void HCrash::crsh_an_c_co(HModel *ptr_model) {
  model = ptr_model;
  const double *colCost = model->getcolCost();
  const double *colLower = model->getcolLower();
  const double *colUpper = model->getcolUpper();
  const double *rowLower = model->getrowLower();
  const double *rowUpper = model->getrowUpper();
  
  int n_ze_c_co = 0;
  int n_fs_c_co = 0;
  
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (colCost[c_n] == 0.0) {
      n_ze_c_co += 1;
      n_fs_c_co += 1;
      continue;
    }
    if (colUpper[c_n] >= HSOL_CONST_INF) {
      // Free column: nonzero cost cannot be feasible
      if (colLower[c_n] > -HSOL_CONST_INF) {
	// Lower-bounded (1-sided) column: non-negative cost is feasible
	if (colCost[c_n] >= 0.0)
	  n_fs_c_co += 1;
      }
    }
    else {
      if (colLower[c_n] <= -HSOL_CONST_INF) {
	// Upper-bonded (1-sided) column: non-positive cost is feasible
	if (colCost[c_n] <= 0.0)
	  n_fs_c_co += 1;
      } else {
	// Two-sided column: any cost is feasible
	n_fs_c_co += 1;
      }
    }
  }
  printf(" Model has %7d Ze costs (%3d%%)\n", n_ze_c_co,
	 (100 * n_ze_c_co) / numCol);
  printf(" Model has %7d Fs costs (%3d%%)\n", n_fs_c_co,
	 (100 * n_fs_c_co) / numCol);
}

void HCrash::crsh_an_r_c_st_af(HModel *ptr_model) {
  model = ptr_model;
  const int *Astart = &model->Astart[0];
  const int *Aindex = &model->Aindex[0];
  const double *Avalue = &model->Avalue[0];;
  for (int k = 0; k < numRow; k++) {
    int vr_n = model->getBaseIndex()[k];
    if (vr_n < numCol) {
      int c_n = vr_n;
      crsh_bs_vr_ty_n_c[crsh_c_ty[c_n]] += 1;
    } else {
      int r_n = vr_n - numCol;
      crsh_bs_vr_ty_n_r[crsh_r_ty[r_n]] += 1;
      
    }
  }
  
  for (int vr_n = 0; vr_n < numTot; vr_n++) {
    if (model->getNonbasicFlag()[vr_n] == 0)
      continue;
    if (vr_n < numCol) {
      int c_n = vr_n;
      crsh_nonbc_vr_ty_n_c[crsh_c_ty[c_n]] += 1;
    } else {
      int r_n = vr_n - numCol;
      crsh_nonbc_vr_ty_n_r[crsh_r_ty[r_n]] += 1;
      
    }
  }
  int bs_mtx_n_struc_el = 0;
  for (int r_n = 0; r_n < numRow; r_n++) {
    int vr_n = model->getBaseIndex()[r_n];
    if (vr_n < numCol) {
      int c_n_el = Astart[vr_n + 1] - Astart[vr_n];
      bs_mtx_n_struc_el += c_n_el;
    }
  }
  
#ifdef JAJH_dev
  crsh_rp_r_c_st(1);
  crsh_rp_r_c_st(2);
  printf(" Basis    matrix    contains%7d structural entries\n",
	 bs_mtx_n_struc_el);
  crsh_rp_r_c_st(3);
#endif
}

void HCrash::crsh_rp_r_c_st(int mode) {
  string TyNm;
#ifdef JAJH_dev
  int ck_su_n_c = 0;
  int ck_su_n_r = 0;
  int ck_su_n_bc_vr = 0;
  int ck_su_n_nonbc_vr = 0;
#endif
  for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_l_vr_ty + 1; vr_ty++) {
    if (vr_ty == crsh_vr_ty_fx)
      TyNm = "Fx ";
    else if (vr_ty == crsh_vr_ty_2_sd)
      TyNm = "2sd";
    else if (vr_ty == crsh_vr_ty_1_sd)
      TyNm = "1sd";
    else if (vr_ty == crsh_vr_ty_fr)
      TyNm = "Fr ";
    else
      printf("Unrecognised type %d\n", vr_ty);
    if (mode == 0) {
#ifdef JAJH_dev
      ck_su_n_r += crsh_vr_ty_og_n_r[vr_ty];
#endif
      if (crsh_vr_ty_og_n_r[vr_ty] > 0)
	printf(" Model has %7d %3s rows (%3d%%)\n",
	       crsh_vr_ty_og_n_r[vr_ty], TyNm.c_str(),
	       (100 * crsh_vr_ty_og_n_r[vr_ty]) / numRow);
    } else if (mode == 1) {
      if (crsh_vr_ty_og_n_r[vr_ty] > 0)
	printf(" Removed %7d of %7d %3s rows (%3d%%)\n",
	       crsh_vr_ty_rm_n_r[vr_ty], crsh_vr_ty_og_n_r[vr_ty],
	       TyNm.c_str(),
	       (100 * crsh_vr_ty_rm_n_r[vr_ty])
	       / crsh_vr_ty_og_n_r[vr_ty]);
    } else if (mode == 2) {
#ifdef JAJH_dev
      ck_su_n_bc_vr += crsh_bs_vr_ty_n_r[vr_ty];
      ck_su_n_bc_vr += crsh_bs_vr_ty_n_c[vr_ty];
      ck_su_n_nonbc_vr += crsh_nonbc_vr_ty_n_r[vr_ty];
      ck_su_n_nonbc_vr += crsh_nonbc_vr_ty_n_c[vr_ty];
#endif
      if (crsh_bs_vr_ty_n_r[vr_ty] > 0)
	printf(" Basic    variables contain %7d %3s rows (%3d%%)\n",
	       crsh_bs_vr_ty_n_r[vr_ty], TyNm.c_str(),
	       (100 * crsh_bs_vr_ty_n_r[vr_ty]) / numRow);
      if (crsh_bs_vr_ty_n_c[vr_ty] > 0)
	printf(" Basic    variables contain %7d %3s cols (%3d%%)\n",
	       crsh_bs_vr_ty_n_c[vr_ty], TyNm.c_str(),
	       (100 * crsh_bs_vr_ty_n_c[vr_ty]) / numRow);
    } else {
      if (crsh_nonbc_vr_ty_n_c[vr_ty] > 0)
	printf(" Nonbasic variables contain %7d %3s cols (%3d%%)\n",
	       crsh_nonbc_vr_ty_n_c[vr_ty], TyNm.c_str(),
	       (100 * crsh_nonbc_vr_ty_n_c[vr_ty]) / numCol);
      if (crsh_nonbc_vr_ty_n_r[vr_ty] > 0)
	printf(" Nonbasic variables contain %7d %3s rows (%3d%%)\n",
	       crsh_nonbc_vr_ty_n_r[vr_ty], TyNm.c_str(),
	       (100 * crsh_nonbc_vr_ty_n_r[vr_ty]) / numCol);
    }
  }
#ifdef JAJH_dev
  if (mode == 0) assert(ck_su_n_r == numRow);
  if (mode == 2) assert(ck_su_n_bc_vr == numRow);
  if (mode == 2) assert(ck_su_n_nonbc_vr == numCol);
#endif
  if (mode <= 1) {
    for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_l_vr_ty + 1; vr_ty++) {
      if (vr_ty == crsh_vr_ty_fx)
	TyNm = "Fx ";
      else if (vr_ty == crsh_vr_ty_2_sd)
	TyNm = "2sd";
      else if (vr_ty == crsh_vr_ty_1_sd)
	TyNm = "1sd";
      else if (vr_ty == crsh_vr_ty_fr)
	TyNm = "Fr ";
      else
	printf("Unrecognised type %d\n", vr_ty);
      if (mode == 0)
#ifdef JAJH_dev
	ck_su_n_c += crsh_vr_ty_og_n_c[vr_ty];
#endif
	if (crsh_vr_ty_og_n_c[vr_ty] > 0)
	  printf(" Model has %7d %3s cols (%3d%%)\n",
		 crsh_vr_ty_og_n_c[vr_ty], TyNm.c_str(),
		 (100 * crsh_vr_ty_og_n_c[vr_ty]) / numCol);
      
	else if (mode == 1) {
	  if (crsh_vr_ty_og_n_c[vr_ty] > 0)
	    printf(" Added   %7d of %7d %3s cols (%3d%%)\n",
		   crsh_vr_ty_add_n_c[vr_ty],
		   crsh_vr_ty_og_n_c[vr_ty], TyNm.c_str(),
		   (100 * crsh_vr_ty_add_n_c[vr_ty])
		   / crsh_vr_ty_og_n_c[vr_ty]);
	}
    }
#ifdef JAJH_dev
    if (mode == 0) assert(ck_su_n_c == numCol);
#endif
  }
  
}

void HCrash::crsh_ck_an_impl_bd() {
  const double *colLower = model->getcolLower();
  const double *colUpper = model->getcolUpper();
  const double *rowLower = model->getrowLower();
  const double *rowUpper = model->getrowUpper();
  
  const double *primalColLowerImplied = model->getprimalColLowerImplied();
  const double *primalColUpperImplied = model->getprimalColUpperImplied();
  const double *primalRowLowerImplied = model->getprimalRowLowerImplied();
  const double *primalRowUpperImplied = model->getprimalRowUpperImplied();
  
  const double *dualColLowerImplied = model->getdualColLowerImplied();
  const double *dualColUpperImplied = model->getdualColUpperImplied();
  const double *dualRowLowerImplied = model->getdualRowLowerImplied();
  const double *dualRowUpperImplied = model->getdualRowUpperImplied();
  
  //	if (numTot < 100) {
  //		printf("\nReporting bounds and corresponding implied bounds\n");
  //		for (int c_n = 0; c_n < numCol; c_n++) {
  //			printf(
  //					"Col %5d has primal bounds [%12.4e,%12.4e] and implied primal bounds [%12.4e,%12.4e]\n",
  //					c_n, colLower[c_n], colUpper[c_n],
  //					primalColLowerImplied[c_n], primalColUpperImplied[c_n]);
  //		}
  //		for (int c_n = 0; c_n < numCol; c_n++) {
  //			printf(
  //					"Col %5d has primal bounds [%12.4e,%12.4e] and implied dual bounds [%12.4e,%12.4e]\n",
  //					c_n, colLower[c_n], colUpper[c_n],
  //					dualColLowerImplied[c_n], dualColUpperImplied[c_n]);
  //		}
  //		for (int r_n = 0; r_n < numRow; r_n++) {
  //			printf(
  //					"Row %5d has primal bounds [%12.4e,%12.4e] and implied primal bounds [%12.4e,%12.4e]\n",
  //					r_n, rowLower[r_n], rowUpper[r_n],
  //					primalRowLowerImplied[r_n], primalRowUpperImplied[r_n]);
  //		}
  //		for (int r_n = 0; r_n < numRow; r_n++) {
  //			printf(
  //					"Row %5d has primal bounds [%12.4e,%12.4e] and implied dual bounds [%12.4e,%12.4e]\n",
  //					r_n, rowLower[r_n], rowUpper[r_n],
  //					dualRowLowerImplied[r_n], dualRowUpperImplied[r_n]);
  //		}
  //	}
  printf("\nAnalysing bounds and corresponding implied bounds\n");
  int numTighterPrimalColLower = 0;
  int numTighterPrimalColUpper = 0;
  int numSlackerPrimalColLower = 0;
  int numSlackerPrimalColUpper = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (primalColLowerImplied[c_n] > colLower[c_n]) {
      numTighterPrimalColLower += 1;
    }
    if (primalColLowerImplied[c_n] < colLower[c_n]) {
      numSlackerPrimalColLower += 1;
      printf(
	     "Col %5d has primal bounds [%12.4e,%12.4e] and slacker implied primal bounds  [%12.4e,%12.4e] Query lower\n",
	     c_n, colLower[c_n], colUpper[c_n],
	     primalColLowerImplied[c_n], primalColUpperImplied[c_n]);
    }
    if (primalColUpperImplied[c_n] < colUpper[c_n]) {
      numTighterPrimalColUpper += 1;
    }
    if (primalColUpperImplied[c_n] > colUpper[c_n]) {
      numSlackerPrimalColUpper += 1;
      printf(
	     "Col %5d has primal bounds [%12.4e,%12.4e] and slacker implied primal bounds  [%12.4e,%12.4e] Query upper\n",
	     c_n, colLower[c_n], colUpper[c_n],
	     primalColLowerImplied[c_n], primalColUpperImplied[c_n]);
    }
  }
  int numTighterDualColLower = 0;
  int numTighterDualColUpper = 0;
  int numSlackerDualColLower = 0;
  int numSlackerDualColUpper = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (colLower[c_n] > -HSOL_CONST_INF) {
      //			Col lower >  -inf so dualColUpperImplied < inf is tighter
      if (dualColUpperImplied[c_n] < HSOL_CONST_INF) {
	numTighterDualColUpper += 1;
      }
    } else {
      //			Col lower <= -inf so dualColUpperImplied should be at most zero
      if (dualColUpperImplied[c_n] > 0) {
	numSlackerDualColUpper += 1;
	printf(
	       "Col %5d has primal bounds [%12.4e,%12.4e] and slacker implied dual bounds  [%12.4e,%12.4e] Query upper\n",
	       c_n, colLower[c_n], colUpper[c_n],
	       dualColLowerImplied[c_n], dualColUpperImplied[c_n]);
      }
    }
    if (colUpper[c_n] < HSOL_CONST_INF) {
      //			Col upper <   inf so dualColLowerImplied > -inf is tighter
      if (dualColLowerImplied[c_n] > -HSOL_CONST_INF) {
	numTighterDualColLower += 1;
      }
    } else {
      //			Col upper >=  inf so dualRowLowerImplied should be at least zero
      if (dualColLowerImplied[c_n] < 0) {
	numSlackerDualColLower += 1;
	printf(
	       "Col %5d has primal bounds [%12.4e,%12.4e] and slacker implied dual bounds  [%12.4e,%12.4e] Query lower\n",
	       c_n, colLower[c_n], colUpper[c_n],
	       dualColLowerImplied[c_n], dualColUpperImplied[c_n]);
      }
    }
  }
  int numTighterPrimalRowLower = 0;
  int numTighterPrimalRowUpper = 0;
  int numSlackerPrimalRowLower = 0;
  int numSlackerPrimalRowUpper = 0;
  for (int r_n = 0; r_n < numRow; r_n++) {
    if (primalRowLowerImplied[r_n] > rowLower[r_n]) {
      numTighterPrimalRowLower += 1;
    }
    if (primalRowLowerImplied[r_n] < rowLower[r_n]) {
      numSlackerPrimalRowLower += 1;
      printf(
	     "row %5d has primal bounds [%12.4e,%12.4e] and slacker implied primal bounds  [%12.4e,%12.4e] Query lower\n",
	     r_n, rowLower[r_n], rowUpper[r_n],
	     primalRowLowerImplied[r_n], primalRowUpperImplied[r_n]);
    }
    if (primalRowUpperImplied[r_n] < rowUpper[r_n]) {
      numTighterPrimalRowUpper += 1;
    }
    if (primalRowUpperImplied[r_n] > rowUpper[r_n]) {
      numSlackerPrimalRowUpper += 1;
      printf(
	     "row %5d has primal bounds [%12.4e,%12.4e] and slacker implied primal bounds  [%12.4e,%12.4e] Query upper\n",
	     r_n, rowLower[r_n], rowUpper[r_n],
	     primalRowLowerImplied[r_n], primalRowUpperImplied[r_n]);
    }
  }
  int numTighterDualRowLower = 0;
  int numTighterDualRowUpper = 0;
  int numSlackerDualRowLower = 0;
  int numSlackerDualRowUpper = 0;
  for (int r_n = 0; r_n < numRow; r_n++) {
    if (rowLower[r_n] > -HSOL_CONST_INF) {
      //			Row lower >  -inf so dualRowUpperImplied < inf is tighter
      if (dualRowUpperImplied[r_n] < HSOL_CONST_INF) {
	numTighterDualRowUpper += 1;
      }
    } else {
      //			Row lower <= -inf so dualRowUpperImplied should be at most zero
      if (dualRowUpperImplied[r_n] > 0) {
	numSlackerDualRowUpper += 1;
	printf(
	       "Row %5d has primal bounds [%12.4e,%12.4e] and slacker implied dual bounds  [%12.4e,%12.4e] Query upper\n",
	       r_n, rowLower[r_n], rowUpper[r_n],
	       dualRowLowerImplied[r_n], dualRowUpperImplied[r_n]);
      }
    }
    if (rowUpper[r_n] < HSOL_CONST_INF) {
      //			Row upper <   inf so dualColLowerImplied > -inf is tighter
      if (dualRowLowerImplied[r_n] > -HSOL_CONST_INF) {
	numTighterDualRowLower += 1;
      }
    } else {
      //			Row upper >=  inf so dualColLowerImplied should be at least zero
      if (dualRowLowerImplied[r_n] < 0) {
	numSlackerDualRowLower += 1;
	printf(
	       "Row %5d has primal bounds [%12.4e,%12.4e] and slacker implied dual bounds  [%12.4e,%12.4e] Query lower\n",
	       r_n, rowLower[r_n], rowUpper[r_n],
	       dualRowLowerImplied[r_n], dualRowUpperImplied[r_n]);
      }
    }
  }
  int numTighter = numTighterPrimalColLower + numTighterPrimalColUpper
    + numTighterDualColLower + numTighterDualColUpper
    + numTighterPrimalRowLower + numTighterPrimalRowUpper
    + numTighterDualRowLower + numTighterDualRowUpper;
  if (numTighter > 0) {
    printf(
	   "Model has %5d/%5d (%3d%%) tighter implied primal column lower bounds\n",
	   numTighterPrimalColLower, numCol,
	   (100 * numTighterPrimalColLower) / numCol);
    printf(
	   "Model has %5d/%5d (%3d%%) tighter implied primal column upper bounds\n",
	   numTighterPrimalColUpper, numCol,
	   (100 * numTighterPrimalColUpper) / numCol);
    printf(
	   "Model has %5d/%5d (%3d%%) tighter implied dual   column lower bounds\n",
	   numTighterDualColLower, numCol,
	   (100 * numTighterDualColLower) / numCol);
    printf(
	   "Model has %5d/%5d (%3d%%) tighter implied dual   column upper bounds\n",
	   numTighterDualColUpper, numCol,
	   (100 * numTighterDualColUpper) / numCol);
    printf(
	   "Model has %5d/%5d (%3d%%) tighter implied primal row    lower bounds\n",
	   numTighterPrimalRowLower, numRow,
	   (100 * numTighterPrimalRowLower) / numRow);
    printf(
	   "Model has %5d/%5d (%3d%%) tighter implied primal row    upper bounds\n",
	   numTighterPrimalRowUpper, numRow,
	   (100 * numTighterPrimalRowUpper) / numRow);
    printf(
	   "Model has %5d/%5d (%3d%%) tighter implied dual   row    lower bounds\n",
	   numTighterDualRowLower, numRow,
	   (100 * numTighterDualRowLower) / numRow);
    printf(
	   "Model has %5d/%5d (%3d%%) tighter implied dual   row    upper bounds\n",
	   numTighterDualRowUpper, numRow,
	   (100 * numTighterDualRowUpper) / numRow);
  } else {
    printf("\nModel has no tighter implied bounds\n");
  }
  int numSlacker = numSlackerPrimalColLower + numSlackerPrimalColUpper
    + numSlackerDualColLower + numSlackerDualColUpper
    + numSlackerPrimalRowLower + numSlackerPrimalRowUpper
    + numSlackerDualRowLower + numSlackerDualRowUpper;
  if (numSlacker > 0) {
    printf(
	   "Model has %5d/%5d (%3d%%) slacker implied primal column lower bounds\n",
	   numSlackerPrimalColLower, numCol,
	   (100 * numSlackerPrimalColLower) / numCol);
    printf(
	   "Model has %5d/%5d (%3d%%) slacker implied primal column upper bounds\n",
	   numSlackerPrimalColUpper, numCol,
	   (100 * numSlackerPrimalColUpper) / numCol);
    printf(
	   "Model has %5d/%5d (%3d%%) slacker implied dual   column lower bounds\n",
	   numSlackerDualColLower, numCol,
	   (100 * numSlackerDualColLower) / numCol);
    printf(
	   "Model has %5d/%5d (%3d%%) slacker implied dual   column upper bounds\n",
	   numSlackerDualColUpper, numCol,
	   (100 * numSlackerDualColUpper) / numCol);
    printf(
	   "Model has %5d/%5d (%3d%%) slacker implied primal row    lower bounds\n",
	   numSlackerPrimalRowLower, numRow,
	   (100 * numSlackerPrimalRowLower) / numRow);
    printf(
	   "Model has %5d/%5d (%3d%%) slacker implied primal row    upper bounds\n",
	   numSlackerPrimalRowUpper, numRow,
	   (100 * numSlackerPrimalRowUpper) / numRow);
    printf(
	   "Model has %5d/%5d (%3d%%) slacker implied dual   row    lower bounds\n",
	   numSlackerDualRowLower, numRow,
	   (100 * numSlackerDualRowLower) / numRow);
    printf(
	   "Model has %5d/%5d (%3d%%) slacker implied dual   row    upper bounds\n",
	   numSlackerDualRowUpper, numRow,
	   (100 * numSlackerDualRowUpper) / numRow);
  } else {
    printf("\nModel has no slacker implied bounds\n\n");
  }
  
}

void HCrash::ltssf(HModel *ptr_model, int Crash_Mode) {
  model = ptr_model;
  if (Crash_Mode == Crash_Mode_LTSSF_k) {
    crsh_fn_cf_pri_v = 1;
    crsh_fn_cf_k = 10;
    alw_al_bs_cg = false;
    no_ck_pv = false;
  } else if (Crash_Mode == Crash_Mode_LTSF_k) {
    crsh_fn_cf_pri_v = 1;
    crsh_fn_cf_k = 10;
    alw_al_bs_cg = false;
    no_ck_pv = true;
  } else if (Crash_Mode == Crash_Mode_LTSF) {
    crsh_fn_cf_pri_v = 1;
    crsh_fn_cf_k = 10;
    alw_al_bs_cg = true;
    no_ck_pv = true;
  } else if (Crash_Mode == Crash_Mode_LTSSF_pri) {
    crsh_fn_cf_pri_v = 10;
    crsh_fn_cf_k = 1;
    alw_al_bs_cg = false;
    no_ck_pv = false;
  } else if (Crash_Mode == Crash_Mode_LTSF_pri) {
    crsh_fn_cf_pri_v = 10;
    crsh_fn_cf_k = 1;
    alw_al_bs_cg = false;
    no_ck_pv = true;
  } else {
    //		Dev version
    crsh_fn_cf_pri_v = 1;
    crsh_fn_cf_k = 10;
    alw_al_bs_cg = false;
    no_ck_pv = false;
  }
  
  mn_co_tie_bk = false;
  
  model->mlFg_Report();
  
  const int *Astart = &model->Astart[0];
  const int *Aindex = &model->Aindex[0];
  const double *Avalue = &model->Avalue[0];;
  
  numRow = model->getNumRow();
  numCol = model->getNumCol();
  numTot = model->getNumTot();
  int numEl = Astart[numCol];
  printf("In crash: LTSSF after Alias to the matrix\n");cout<<flush;
  
  //	for (int r_n = 0; r_n < numRow; r_n++)
  //		printf(
  //				"Row: %d; BaseIndex = %d; Work[L,V,U] = [%g,%g,%g]; Base[L,V,U] = [%g,%g,%g]\n",
  //				r_n, model->getBaseIndex()[r_n], model->getWorkLower()[r_n],
  //				model->getWorkValue()[r_n], model->getWorkUpper()[r_n],
  //				model->getBaseLower()[r_n], model->getBaseValue()[r_n],
  //				model->getBaseUpper()[r_n]);
  //	for (int vr_n = 0; vr_n < numTot; vr_n++)
  //		printf("Var: %d; NonbasicFlag = %d; NonbasicMove = %d; Work[L,V,U] = [%g,%g,%g]\n", vr_n,
  //				model->getNonbasicFlag()[vr_n],
  //				model->getNonbasicMove()[vr_n], model->getWorkLower()[vr_n],
  //				model->getWorkValue()[vr_n], model->getWorkUpper()[vr_n]);
  
  //	  const double *colCost = model->getcolCost();
  //    const double *colLower = model->getcolLower();
  //    const double *colUpper = model->getcolUpper();
  //    const double *rowLower = model->getrowLower();
  //    const double *rowUpper = model->getrowUpper();
  
  
  ltssf_iz_da(ptr_model);
#ifdef JAJH_dev
  printf("\nLTSSF Crash\n");
  
  printf(" crsh_fn_cf_pri_v = %d\n", crsh_fn_cf_pri_v);
  printf(" crsh_fn_cf_k = %d\n", crsh_fn_cf_k);
  printf(" Objective f = %2d*pri_v - %2d*k\n", crsh_fn_cf_pri_v,
	 crsh_fn_cf_k);
  string lg;
  //	if (mn_co_tie_bk) lg = "T"; else lg = "F";
  //	printf(" Break priority function ties on smallest cost = %s\n", lg);
  if (alw_al_bs_cg)
    printf(" Allow any basis change regardless of priority = T\n");
  else
    printf(" Allow any basis change regardless of priority = F\n");
  if (no_ck_pv)
    printf(" Don't do numerical check on pivot = T\n");
  else
    printf(" Don't do numerical check on pivot = F\n");
  
  printf("Max row priority is %d\n", mx_r_pri);
  printf("Max col priority is %d\n", mx_c_pri);
#endif
  if ((!alw_al_bs_cg) && (mx_r_pri + mx_c_pri <= crsh_mx_pri_v)) {
#ifdef JAJH_dev
    printf(
	   "Max row priority of %d + Max col priority of %d = %d <= %d: no value in performing LTSSF crash\n",
	   mx_r_pri, mx_c_pri, mx_r_pri + mx_c_pri, crsh_mx_pri_v);
#endif
    // Save the solved results
    return;
  }
  //	for (int r_n = 0; r_n < numRow; r_n++) {
  //		if (crsh_act_r[r_n] != 0)
  //			printf("Row %2d , Pri = %1d; K = %2d\n", r_n, crsh_r_ty_pri_v[crsh_r_ty[r_n]], crsh_r_k[r_n]);
  //	}
#ifdef JAJH_dev
  if (ltssf_ck_fq > 0) {
    printf("\nCHECKING LTSSF DATA NOW AND EVERY %d PASS(ES)!!\n\n",
	   ltssf_ck_fq);
    ltssf_ck_da();
  }
#endif
  if (Rp_TwoD_da)
    ltssf_rp_pri_k_da();
  
#ifdef JAJH_dev
  crsh_rp_r_c_st(0);
#endif
  ltssf_iterate(ptr_model);
#ifdef JAJH_dev
  printf(" %d/%d basis changes from %d passes\n", n_crsh_bs_cg, numRow,
	 n_crsh_ps);
  printf(
	 " Absolute tolerance (%6.4f): Rejected %7d pivots: min absolute pivot value = %6.4e\n",
	 tl_crsh_abs_pv_v, n_abs_pv_no_ok, mn_abs_pv_v);
  printf(
	 " Relative tolerance (%6.4f): Rejected %7d pivots: min relative pivot value = %6.4e\n",
	 tl_crsh_rlv_pv_v, n_rlv_pv_no_ok, mn_rlv_pv_v);
  crsh_an_r_c_st_af(ptr_model);
#endif
}

void HCrash::ltssf_iz_mode(int Crash_Mode) {
  //	TODO Would be good to be able to call this from HDual::setCrash
  crsh_fn_cf_pri_v = 1;
  crsh_fn_cf_k = 10;
  alw_al_bs_cg = false;
  no_ck_pv = false;
  
}

void HCrash::ltssf_iterate(HModel *ptr_model) {
  //	LTSSF Main loop
  n_crsh_ps = 0;
  n_crsh_bs_cg = 0;
  bool ltssf_stop = false;
  for (;;) {
    //		printf("\n LTSSF Crash pass %d\n", n_crsh_ps);
    //		ltssf_rp_r_pri();
    ltssf_cz_r();
    if (cz_r_n == no_ix)
      break;
    cz_r_pri_v = crsh_r_ty_pri_v[crsh_r_ty[cz_r_n]];
    ltssf_cz_c(ptr_model);
    bool bs_cg = cz_c_n != no_ix;
    if (bs_cg) {
      if (Rp_TwoD_da)
	printf("Pass %2d; cz_r = %2d; cz_c = %2d\n", n_crsh_ps, cz_r_n,
	       cz_c_n);
      //			A basis change has occurred
      n_crsh_bs_cg += 1;
      double abs_pv_v = abs(pv_v);
      double rlv_pv_v = abs_pv_v / crsh_mtx_c_mx_abs_v[cz_c_n];
      mn_abs_pv_v = min(abs_pv_v, mn_abs_pv_v);
      mn_rlv_pv_v = min(rlv_pv_v, mn_rlv_pv_v);
      int columnIn = cz_c_n;
      int rowOut = cz_r_n;
      int columnOut = numCol + cz_r_n;
      int sourceOut = 1;
      //			Update the basic/nonbasic variable info and the row-wise copy of the matrix
      model->updatePivots(columnIn, rowOut, sourceOut);
      if (model->mlFg_haveMatrixRowWise)
	model->updateMatrix(columnIn, columnOut);
      //		    Update the count of this type of removal and addition
      int vr_ty = crsh_r_ty[cz_r_n];
      crsh_vr_ty_rm_n_r[vr_ty] += 1;
      vr_ty = crsh_c_ty[cz_c_n];
      crsh_vr_ty_add_n_c[vr_ty] += 1;
    } else {
      if (Rp_TwoD_da) printf("Pass %2d; cz_r = %2d: No basis change\n", n_crsh_ps, cz_r_n);
    }
    if (Rp_TwoD_da) ltssf_rp_pri_k_da();
    ltssf_u_da(ptr_model);
    //for (int r_n = 0; r_n < numRow; r_n++) {
    //if (crsh_act_r[r_n] != 0)
    //printf("Row %2d , Pri = %1d; K = %2d\n", r_n, crsh_r_ty_pri_v[crsh_r_ty[r_n]], crsh_r_k[r_n]);
    //}
    //Check LTSSF data every ltssf_ck_fq passes (if ltssf_ck_fq>0)
    if ((ltssf_ck_fq > 0) && (n_crsh_ps % ltssf_ck_fq == 0)) ltssf_ck_da();
    if (Rp_TwoD_da) ltssf_rp_pri_k_da();
    int mn_r_k = numCol + 1;
    if (TwoD_hdr) {
      //			Determine whether the are still rows worth removing
      mx_r_pri = crsh_mn_pri_v - 1;
      for (int pri_v = crsh_mx_pri_v; pri_v > crsh_mn_pri_v; pri_v--) {
	if (crsh_r_pri_mn_r_k[pri_v] < numCol + 1) {
	  mx_r_pri = pri_v;
	  break;
	}
      }
      if ((!alw_al_bs_cg) && (mx_r_pri + mx_c_pri <= crsh_mx_pri_v)) {
	printf("Max active row priority of %d + Max original col priority of %d = %d <= %d: no value in performing further LTSSF crash\n",
	       mx_r_pri, mx_c_pri, mx_r_pri + mx_c_pri, crsh_mx_pri_v);
	ltssf_stop = true;
      }
    }
    
    n_crsh_ps += 1;
    if (ltssf_stop) break;
    //if (n_crsh_ps==6) break;
  }
}

void HCrash::ltssf_cz_r() {
  cz_r_n = no_ix;
  //		Choose a row from the active submatrix
  if (crsh_fn_cf_pri_v > crsh_fn_cf_k) {
    //		Search according to maximum priority function
    if (OneD_hdr) {
      int mx_r_pri_fn_v = -HSOL_CONST_I_INF;
      for (int pri_v = crsh_mx_pri_v; pri_v > crsh_mn_pri_v; pri_v--) {
	int r_n = crsh_r_pri_hdr[pri_v];
	if (r_n == no_lk)
	  continue;
	//     Note the limiting priority function value so that the search can be stopped if it is attained
	//     corresponds to a maximum priority
	int lm_r_pri_fn_v = crsh_fn_cf_pri_v * pri_v - crsh_fn_cf_k;
	do {
	  int k = crsh_r_k[r_n];
	  int r_pri_fn_v = crsh_fn_cf_pri_v * pri_v
	    - crsh_fn_cf_k * k;
	  if (r_pri_fn_v > mx_r_pri_fn_v) {
	    cz_r_n = r_n;
	    mx_r_pri_fn_v = r_pri_fn_v;
	    if (mx_r_pri_fn_v == lm_r_pri_fn_v)
	      break;
	  }
	  r_n = crsh_r_pri_lkf[r_n];
	} while (r_n != no_lk);
	if (cz_r_n != no_ix)
	  break;
      }
    }
    if (TwoD_hdr) {
      for (int pri_v = crsh_mx_pri_v; pri_v > crsh_mn_pri_v; pri_v--) {
	int r_k = crsh_r_pri_mn_r_k[pri_v];
	if (r_k > numCol)
	  continue;
	cz_r_n = crsh_r_pri_k_hdr[pri_v * (numCol + 1) + r_k];
	if (cz_r_n == no_ix) {
	  printf(
		 "ERROR: header for pri_v = %d and count = %d is empty for crsh_r_pri_mn_r_k[pri_v] = %d\n",
		 pri_v, r_k, crsh_r_pri_mn_r_k[pri_v]);
	}
	break;
      }
    }
  } else {
    //		Search according to number of active entries
    if (OneD_hdr) {
      int mx_r_pri_fn_v = -HSOL_CONST_I_INF;
      for (int k = 1; k < numCol + 1; k++) {
	int r_n = crsh_r_k_hdr[k];
	// printf("Priority %d: Header %d\n", k, r_n);
	if (r_n == no_ix)
	  continue;
	// Note the limiting priority function value so that the search can be stopped if it is attained
	// corresponds to a maximum priority
	int lm_r_pri_fn_v = crsh_fn_cf_pri_v * mx_r_pri_v
	  - crsh_fn_cf_k * k;
	do {
	  int pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
	  int r_pri_fn_v = crsh_fn_cf_pri_v * pri_v
	    - crsh_fn_cf_k * k;
	  if (r_pri_fn_v > mx_r_pri_fn_v) {
	    cz_r_n = r_n;
	    mx_r_pri_fn_v = r_pri_fn_v;
	    if (mx_r_pri_fn_v == lm_r_pri_fn_v)
	      break;
	  }
	  r_n = crsh_r_k_lkf[r_n];
	} while (r_n != no_lk);
	if (cz_r_n != no_ix)
	  break;
      }
    }
    if (TwoD_hdr) {
      int mn_r_k = numCol + 1;
      for (int pri_v = crsh_mx_pri_v; pri_v > crsh_mn_pri_v; pri_v--) {
	int r_k = crsh_r_pri_mn_r_k[pri_v];
	if (r_k < mn_r_k) {
	  cz_r_n = crsh_r_pri_k_hdr[pri_v * (numCol + 1) + r_k];
	  if (cz_r_n == no_ix) {
	    printf(
		   "ERROR: header for pri_v = %d and count = %d is empty for crsh_r_pri_mn_r_k[pri_v] = %d\n",
		   pri_v, r_k, crsh_r_pri_mn_r_k[pri_v]);
	  }
	  mn_r_k = r_k;
	  if (mn_r_k == 1)
	    break;
	}
      }
    }
  }
}

void HCrash::ltssf_cz_c(HModel *ptr_model) {
  model = ptr_model;
  //  matrix = model->getMatrix();
  const double *colCost = model->getcolCost();
  const double *colLower = model->getcolLower();
  const double *colUpper = model->getcolUpper();
  const double *rowLower = model->getrowLower();
  const double *rowUpper = model->getrowUpper();
  
  cz_c_n = no_ix;
  int su_r_c_pri_v_lm = crsh_mx_pri_v;
  //     Choose a column which has maximium priority function amongst those
  //     with entries in the selected row---making sure that the pivot is
  //     acceptable numerically.
  //
  //     Note
  //
  //     *  Numerical checking is switched off by setting the tolerances to
  //     .  zero.
  //
  //     *  Make sure that tl_crsh_rlv_pv_v < 1 otherwise test
  //     .  abs_c_v/crsh_mtx_c_mx_abs_v(c_n) .gt. tl_crsh_rlv_pv_v
  //     .  will never be satisfied
  
  if (alw_al_bs_cg)
    int su_r_c_pri_v_lm = -crsh_mx_pri_v;
  else
    int su_r_c_pri_v_lm = crsh_mx_pri_v;
  
  n_eqv_c = 0;
  pv_v = 0.0;
  double mn_co = HSOL_CONST_INF;
  int mx_c_pri_fn_v = -HSOL_CONST_I_INF;
  for (int el_n = CrshARstart[cz_r_n]; el_n < CrshARstart[cz_r_n + 1];
       el_n++) {
    int c_n = CrshARindex[el_n];
    if (crsh_act_c[c_n] == crsh_vr_st_no_act)
      continue;
    //     Don't allow the row to be replaced by a column whose priority to
    //     remain nonbasic is the same or greater.
    if (!alw_al_bs_cg
	&& (crsh_c_ty_pri_v[crsh_c_ty[c_n]] + cz_r_pri_v
	    <= su_r_c_pri_v_lm))
      continue;
    //			If column is worse than current best then break
    int c_pri_fn_v = crsh_fn_cf_pri_v * crsh_c_ty_pri_v[crsh_c_ty[c_n]]
      - crsh_fn_cf_k * crsh_c_k[c_n];
    if (c_pri_fn_v < mx_c_pri_fn_v)
      continue;
    //Pivot is OK if pivots are not being checked
    bool pv_ok = no_ck_pv;
    if (!no_ck_pv) {
      //     Check the matrix entry if it may be used as a pivot.
      double abs_c_v = abs(CrshARvalue[el_n]);
      bool abs_pv_v_ok = abs_c_v > tl_crsh_abs_pv_v;
      bool rlv_pv_v_ok = abs_c_v
	> tl_crsh_rlv_pv_v * crsh_mtx_c_mx_abs_v[c_n];
      if (!abs_pv_v_ok)
	n_abs_pv_no_ok += 1;
      if (!rlv_pv_v_ok)
	n_rlv_pv_no_ok += 1;
      pv_ok = abs_pv_v_ok && rlv_pv_v_ok;
    }
    //If the pivot is not OK then break
    if (!pv_ok)
      continue;
    //    Pivot is OK [or not being checked]
    if (c_pri_fn_v > mx_c_pri_fn_v) {
      //				Better column than best so far
      mx_c_pri_fn_v = c_pri_fn_v;
      cz_c_n = c_n;
      pv_v = CrshARvalue[el_n];
      n_eqv_c = 1;
      mn_co = colCost[c_n];
    } else if (c_pri_fn_v == mx_c_pri_fn_v) {
      //				Matches the best column so far: possibly break tie on cost
      if (mn_co_tie_bk && (colCost[c_n] < mn_co)) {
	cz_c_n = c_n;
	pv_v = CrshARvalue[el_n];
	n_eqv_c = 1;
	mn_co = colCost[c_n];
      }
      n_eqv_c += 1;
    }
  }
}

void HCrash::ltssf_u_da(HModel *ptr_model) {
  if ((cz_r_n != no_ix) && (cz_c_n != no_ix)) {
    ltssf_u_da_af_bs_cg(ptr_model);
  } else {
    ltssf_u_da_af_no_bs_cg();
  }
  //	     If there are no more rows with the current maximum row priority then get the new maximum row priority value
  if (OneD_hdr) {
    if ((crsh_r_pri_hdr[cz_r_pri_v] == no_lk)
	&& (cz_r_pri_v == mx_r_pri_v)) {
      mx_r_pri_v = -HSOL_CONST_I_INF;
      for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++)
	if (crsh_r_pri_hdr[pri_v] != no_lk)
	  mx_r_pri_v = pri_v;
    }
  }
  if (TwoD_hdr) {
    //		TODO Surely this is not necessary with 2-d headers
    if ((crsh_r_pri_mn_r_k[cz_r_pri_v] > numCol)
	&& (cz_r_pri_v == mx_r_pri_v)) {
      mx_r_pri_v = -HSOL_CONST_I_INF;
      for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++)
	if (crsh_r_pri_mn_r_k[pri_v] <= numCol)
	  mx_r_pri_v = pri_v;
    }
  }
}

void HCrash::ltssf_u_da_af_bs_cg(HModel *ptr_model) {
  model = ptr_model;
  const int *Astart = &model->Astart[0];
  const int *Aindex = &model->Aindex[0];
  const double *Avalue = &model->Avalue[0];;
  //	ltssf_rp_r_k();
  for (int r_el_n = CrshARstart[cz_r_n]; r_el_n < CrshARstart[cz_r_n + 1];
       r_el_n++) {
    int c_n = CrshARindex[r_el_n];
    if (crsh_act_c[c_n] == crsh_vr_st_no_act)
      continue;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      if (crsh_act_r[r_n] == crsh_vr_st_no_act)
	continue;
      //     Remove the row from the linked list with this number of active entries
      int prev_r_n;
      int r_k = crsh_r_k[r_n];
      if (OneD_hdr) {
	int nx_r_n = crsh_r_k_lkf[r_n];
	if (r_n == crsh_r_k_hdr[r_k]) {
	  prev_r_n = no_lk;
	  crsh_r_k_hdr[r_k] = nx_r_n;
	} else {
	  prev_r_n = crsh_r_k_lkb[r_n];
	  crsh_r_k_lkf[prev_r_n] = nx_r_n;
	}
	if (nx_r_n != no_lk)
	  crsh_r_k_lkb[nx_r_n] = prev_r_n;
	//     Reduce the number of active entries in this row by one and...
	r_k -= 1;
	crsh_r_k[r_n] = r_k;
	if (r_k > 0) {
	  //     ... either add the row as the header of the list with one fewer number of active entries...
	  nx_r_n = crsh_r_k_hdr[r_k];
	  crsh_r_k_hdr[r_k] = r_n;
	  crsh_r_k_lkb[r_n] = no_lk;
	  crsh_r_k_lkf[r_n] = nx_r_n;
	  if (nx_r_n != no_lk)
	    crsh_r_k_lkb[nx_r_n] = r_n;
	} else {
	  //     ...or, if the count is zero, the row leaves the active submatrix and the priority value data structure
	  crsh_act_r[r_n] = crsh_vr_st_no_act;
	  //     Unless this is the chosen row, the logical remains basic.
	  int pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
	  nx_r_n = crsh_r_pri_lkf[r_n];
	  if (r_n == crsh_r_pri_hdr[pri_v]) {
	    prev_r_n = no_lk;
	    crsh_r_pri_hdr[pri_v] = nx_r_n;
	  } else {
	    prev_r_n = crsh_r_pri_lkb[r_n];
	    crsh_r_pri_lkf[prev_r_n] = nx_r_n;
	  }
	  if (nx_r_n != no_lk)
	    crsh_r_pri_lkb[nx_r_n] = prev_r_n;
	}
      }
      if (TwoD_hdr) {
	if (OneD_hdr) {
	  //TODO Remove this r_k += 1; when only using TwoD_hdr
	  r_k += 1;
	}
	int pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
	if (Rp_TwoD_da) {
	  ltssf_rp_pri_k_da();
	  printf(
		 "1: Remove row %d of pri %d from linked list with %d entries\n",
		 r_n, pri_v, r_k);
	}
	int hdr_ix = pri_v * (numCol + 1) + r_k;
	//     Remove the row from the linked list with this number of active entries
	int prev_r_n;
	int nx_r_n = crsh_r_pri_k_lkf[r_n];
	if (r_n == crsh_r_pri_k_hdr[hdr_ix]) {
	  prev_r_n = no_lk;
	  crsh_r_pri_k_hdr[hdr_ix] = nx_r_n;
	} else {
	  prev_r_n = crsh_r_pri_k_lkb[r_n];
	  crsh_r_pri_k_lkf[prev_r_n] = nx_r_n;
	}
	if (nx_r_n != no_lk)
	  crsh_r_pri_k_lkb[nx_r_n] = prev_r_n;
	if ((crsh_r_pri_k_hdr[hdr_ix] == no_lk)
	    && (crsh_r_pri_mn_r_k[pri_v] == r_k)) {
	  //This was the only row of minimum row count so look for the next row count with non-null header
	  //					Set crsh_r_pri_mn_r_k to numCol + 1 in case r_k=numCol so priority is cleared
	  crsh_r_pri_mn_r_k[pri_v] = numCol + 1;
	  for (int qy_k = r_k + 1; qy_k < numCol + 1; qy_k++) {
	    int hdr_ix = pri_v * (numCol + 1) + qy_k;
	    if (crsh_r_pri_k_hdr[hdr_ix] != no_lk) {
	      crsh_r_pri_mn_r_k[pri_v] = qy_k;
	      break;
	    }
	  }
	  //					if (crsh_r_pri_mn_r_k[pri_v] == numCol + 1) printf("Priority %d now has no rows\n", pri_v);
	}
	//     Reduce the number of active entries in this row by one and...
	r_k -= 1;
	crsh_r_k[r_n] = r_k;
	if (r_k > 0) {
	  //     ... either add the row as the header of the list with one fewer number of active entries...
	  if (Rp_TwoD_da) {
	    ltssf_rp_pri_k_da();
	    printf(
		   "Add row %d of pri %d to linked list with %d entries\n",
		   r_n, pri_v, r_k);
	  }
	  int hdr_ix = pri_v * (numCol + 1) + r_k;
	  nx_r_n = crsh_r_pri_k_hdr[hdr_ix];
	  crsh_r_pri_k_hdr[hdr_ix] = r_n;
	  crsh_r_pri_k_lkf[r_n] = nx_r_n;
	  if (nx_r_n != no_lk)
	    crsh_r_pri_k_lkb[nx_r_n] = r_n;
	  if (crsh_r_pri_mn_r_k[pri_v] > r_k) {
	    //There is now a row of smaller count for this priority
	    crsh_r_pri_mn_r_k[pri_v] = r_k;
	  }
	} else {
	  if (Rp_TwoD_da) {
	    ltssf_rp_pri_k_da();
	    printf(
		   "2: Remove row %d of pri %d and count %d from active submatrix\n",
		   r_n, pri_v, r_k);
	  }
	  //     ...or, if the count is zero, the row leaves the active submatrix...
	  crsh_act_r[r_n] = crsh_vr_st_no_act;
	  //     ... and has already left the priority value and count data structure
	}
      }
    }
    //     The column leaves the active submatrix
    crsh_act_c[c_n] = crsh_vr_st_no_act;
    //     Unless this is the chosen column, the structural remains nonbasic.
  }
}

void HCrash::ltssf_u_da_af_no_bs_cg() {
  // The basis has not changed.  The row becomes inactive and the
  // number of active entries in each column in which it has entries
  // is reduced by one.
  for (int r_el_n = CrshARstart[cz_r_n]; r_el_n < CrshARstart[cz_r_n + 1]; r_el_n++) {
    int c_n = CrshARindex[r_el_n];
    //printf("Deactivating chosen row %2d: Entry %4d is column %2d has activity %2d\n", cz_r_n, r_el_n, c_n, crsh_act_c[c_n]);
    if (crsh_act_c[c_n] == crsh_vr_st_no_act)
      continue;
    crsh_c_k[c_n] -= 1;
    // If the number of active entries in the column is zeroed then it becomes inactive.
    if (crsh_c_k[c_n] == 0)
      crsh_act_c[c_n] = crsh_vr_st_no_act;
  }
  int r_n = cz_r_n;
  //	     Remove the row from the linked list with this number of active entries
  if (OneD_hdr) {
    int r_k = crsh_r_k[r_n];
    int prev_r_n;
    int nx_r_n = crsh_r_k_lkf[r_n];
    if (r_n == crsh_r_k_hdr[r_k]) {
      prev_r_n = no_lk;
      crsh_r_k_hdr[r_k] = nx_r_n;
    } else {
      prev_r_n = crsh_r_k_lkb[r_n];
      crsh_r_k_lkf[prev_r_n] = nx_r_n;
    }
    if (nx_r_n != no_lk)
      crsh_r_k_lkb[nx_r_n] = prev_r_n;
    crsh_act_r[r_n] = crsh_vr_st_no_act;
    int pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
    nx_r_n = crsh_r_pri_lkf[r_n];
    if (r_n == crsh_r_pri_hdr[pri_v]) {
      prev_r_n = no_lk;
      crsh_r_pri_hdr[pri_v] = nx_r_n;
    } else {
      prev_r_n = crsh_r_pri_lkb[r_n];
      crsh_r_pri_lkf[prev_r_n] = nx_r_n;
    }
    if (nx_r_n != no_lk)
      crsh_r_pri_lkb[nx_r_n] = prev_r_n;
  }
  if (TwoD_hdr) {
    //	     Remove the row from the linked list with this priority and number of active entries
    crsh_act_r[r_n] = crsh_vr_st_no_act;
    int pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
    int r_k = crsh_r_k[r_n];
    int hdr_ix = pri_v * (numCol + 1) + r_k;
    int prev_r_n;
    int nx_r_n = crsh_r_pri_k_lkf[r_n];
    if (r_n == crsh_r_pri_k_hdr[hdr_ix]) {
      prev_r_n = no_lk;
      crsh_r_pri_k_hdr[hdr_ix] = nx_r_n;
    } else {
      prev_r_n = crsh_r_pri_k_lkb[r_n];
      crsh_r_pri_k_lkf[prev_r_n] = nx_r_n;
    }
    if (nx_r_n != no_lk)
      crsh_r_pri_k_lkb[nx_r_n] = prev_r_n;
    if ((crsh_r_pri_k_hdr[hdr_ix] == no_lk)
	&& (crsh_r_pri_mn_r_k[pri_v] == r_k)) {
      //No more rows of this count for this priority value: look for next highest
      //Set crsh_r_pri_mn_r_k to numCol + 1 in case r_k=numCol so priority is cleared
      crsh_r_pri_mn_r_k[pri_v] = numCol + 1;
      for (int qy_k = r_k + 1; qy_k < numCol + 1; qy_k++) {
	int hdr_ix = pri_v * (numCol + 1) + qy_k;
	if (crsh_r_pri_k_hdr[hdr_ix] != no_lk) {
	  crsh_r_pri_mn_r_k[pri_v] = qy_k;
	  break;
	}
      }
      //if (crsh_r_pri_mn_r_k[pri_v] == numCol + 1)	printf("Priority %d now has no rows\n", pri_v);
    }
  }
}

void HCrash::ltssf_iz_da(HModel *ptr_model) {
  //	bool ImpliedDualLTSSF = false;
  //	ImpliedDualLTSSF = true;
  model = ptr_model;
  const int *Astart = &model->Astart[0];
  const int *Aindex = &model->Aindex[0];
  const double *Avalue = &model->Avalue[0];;
  int numEl = Astart[numCol];
  //	const double *primalColLowerImplied = model->getprimalColLowerImplied();
  //	const double *primalColUpperImplied = model->getprimalColUpperImplied();
  //	const double *primalRowLowerImplied = model->getprimalRowLowerImplied();
  //	const double *primalRowUpperImplied = model->getprimalRowUpperImplied();
  //
  //	const double *dualColLowerImplied = model->getdualColLowerImplied();
  //	const double *dualColUpperImplied = model->getdualColUpperImplied();
  //	const double *dualRowLowerImplied = model->getdualRowLowerImplied();
  //	const double *dualRowUpperImplied = model->getdualRowUpperImplied();
  
  //Allocate the crash variable type arrays
  crsh_r_ty_pri_v.resize(crsh_l_vr_ty);
  crsh_c_ty_pri_v.resize(crsh_l_vr_ty);
  //Basis-preserving crash:
  //    crsh_r_ty_pri_v[crsh_vr_ty_non_bc] = 1;
  //    crsh_r_ty_pri_v[crsh_vr_ty_bc] =     0;
  //    crsh_c_ty_pri_v[crsh_vr_ty_non_bc] = 0;
  //    crsh_c_ty_pri_v[crsh_vr_ty_bc] =     1;
  //Standard crash:
  crsh_r_ty_pri_v[crsh_vr_ty_fx] = 3;
  crsh_r_ty_pri_v[crsh_vr_ty_2_sd] = 2;
  crsh_r_ty_pri_v[crsh_vr_ty_1_sd] = 1;
  crsh_r_ty_pri_v[crsh_vr_ty_fr] = 0;
  crsh_c_ty_pri_v[crsh_vr_ty_fx] = 0;
  crsh_c_ty_pri_v[crsh_vr_ty_2_sd] = 1;
  crsh_c_ty_pri_v[crsh_vr_ty_1_sd] = 2;
  crsh_c_ty_pri_v[crsh_vr_ty_fr] = 3;
  
  //Allocate the arrays required for crash
  crsh_r_k.resize(numRow);
  crsh_c_k.resize(numCol);
  
  if (TwoD_hdr) {
    //    Add one since need [0..crsh_mx_pri_v]
    crsh_r_pri_k_hdr.resize((crsh_mx_pri_v + 1) * (numCol + 1));
    crsh_r_pri_k_lkf.resize(numRow);
    crsh_r_pri_k_lkb.resize(numRow);
    crsh_r_pri_mn_r_k.resize(crsh_mx_pri_v + 1);
  }
  
  if (OneD_hdr) {
    crsh_r_pri_hdr.resize(crsh_mx_pri_v + 1);
    crsh_r_pri_lkb.resize(numRow);
    crsh_r_pri_lkf.resize(numRow);
    
    //    Add one since numbers of active entries can be 0..numCol.
    crsh_r_k_hdr.resize(numCol + 1);
    crsh_r_k_lkb.resize(numRow);
    crsh_r_k_lkf.resize(numRow);
  }
  crsh_mtx_c_mx_abs_v.resize(numCol);
  CrshARvalue.resize(numEl);
  CrshARindex.resize(numEl);
  CrshARstart.resize(numRow + 1);
  crsh_act_r.resize(numRow);
  crsh_act_c.resize(numCol);
  
  crsh_iz_vr_ty(ptr_model);
  
  mx_r_pri = crsh_mn_pri_v;
  for (int r_n = 0; r_n < numRow; r_n++) mx_r_pri = max(mx_r_pri, crsh_r_ty_pri_v[crsh_r_ty[r_n]]);
  mx_c_pri = crsh_mn_pri_v;
  for (int c_n = 0; c_n < numCol; c_n++) mx_c_pri = max(mx_c_pri, crsh_c_ty_pri_v[crsh_c_ty[c_n]]);
  
  if ((!alw_al_bs_cg) && (mx_r_pri + mx_c_pri <= crsh_mx_pri_v))
    return;
  if (OneD_hdr) {
    //Initialise the linked lists of priority sets
    for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++)
      crsh_r_pri_hdr[pri_v] = no_lk;
    for (int c_n = 0; c_n < numCol + 1; c_n++)
      crsh_r_k_hdr[c_n] = no_lk;
  }
  if (TwoD_hdr) {
    for (int c_n = 0; c_n < numCol + 1; c_n++)
      for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++)
	crsh_r_pri_k_hdr[pri_v * (numCol + 1) + c_n] = no_lk;
  }
  mn_abs_pv_v = HSOL_CONST_INF;
  mn_rlv_pv_v = HSOL_CONST_INF;
  n_abs_pv_no_ok = 0;
  n_rlv_pv_no_ok = 0;
  //Determine the status and type of each row
  for (int r_n = 0; r_n < numRow; r_n++) {
    int vr_ty = crsh_r_ty[r_n];
    int pri_v = crsh_r_ty_pri_v[vr_ty];
    if (pri_v == crsh_no_act_pri_v) {
      //	Rows with no priority value are free - and should be removed by presolve
      //	They will be basic
      crsh_act_r[r_n] = crsh_vr_st_no_act;
    } else {
      crsh_act_r[r_n] = crsh_vr_st_act;
    }
    //	Initialise the count of active entries in the row
    crsh_r_k[r_n] = 0;
  }
  
  //	Determine the status and type of each column and
  //	compute the number of active entries in each row and column
  for (int c_n = 0; c_n < numCol; c_n++) {
    int vr_ty = crsh_c_ty[c_n];
    int pri_v = crsh_c_ty_pri_v[vr_ty];
    int c_n_en = Astart[c_n + 1] - Astart[c_n];
    if ((pri_v == crsh_no_act_pri_v) || (c_n_en == 0)) {
      //	Columns with no priority value are fixed or zero - and should be removed by presolve
      //	They will be nonbasic
      crsh_act_c[c_n] = crsh_vr_st_no_act;
    } else {
      crsh_act_c[c_n] = crsh_vr_st_act;
    }
    //		TODO Count the original number of columns of this type
    if (crsh_act_c[c_n] == crsh_vr_st_act) {
      //			The column is active: count the active rows and
      //			update the number of active entries in those rows
      //
      //			Initialise the count of active entries in the row
      crsh_c_k[c_n] = 0;
      for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
	int r_n = Aindex[el_n];
	if (crsh_act_r[r_n] != crsh_vr_st_no_act) {
	  //			The row is active so increase the number of active entries in the column and row
	  crsh_c_k[c_n] += 1;
	  crsh_r_k[r_n] += 1;
	}
      }
    }
    //Error surely, since it ends up doubling hte numbers of columns of each type
    //    crsh_vr_ty_og_n_c[vr_ty] += 1;
  }
  //Now that the row counts are known, make any zero rows non-active.
  //Form linked list of active rows with each priority.
  //	for (int r_n = 0; r_n < numRow; r_n++) {
  //		int vr_ty = crsh_r_ty[r_n];
  //		int pri_v = crsh_r_ty_pri_v[vr_ty];
  //		printf("Row %d [%10g,%10g] has type %d; priority %d; count = %d\n", r_n, rowLower[r_n], rowUpper[r_n], vr_ty, pri_v, crsh_r_k[r_n]);
  //	}
  
  if (TwoD_hdr) {
    for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++)
      crsh_r_pri_mn_r_k[pri_v] = numCol + 1;
  }
  
  for (int r_n = 0; r_n < numRow; r_n++) {
    if (crsh_act_r[r_n] == crsh_vr_st_no_act)
      continue;
    //    	Row appears active but check whether it is empty
    if (crsh_r_k[r_n] <= 0) {
      crsh_act_r[r_n] = crsh_vr_st_no_act;
      continue;
    }
    //		Add as the header of the appropriate priority value list
    int pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
    if (OneD_hdr) {
      int nx_r_n = crsh_r_pri_hdr[pri_v];
      
      //		printf("Row %2d; Ty = %1d; Pri = %1d; Nx = %2d\n", r_n, crsh_r_ty[r_n], pri_v, nx_r_n);
      crsh_r_pri_hdr[pri_v] = r_n;
      crsh_r_pri_lkb[r_n] = no_lk;
      crsh_r_pri_lkf[r_n] = nx_r_n;
      if (nx_r_n != no_lk)
	crsh_r_pri_lkb[nx_r_n] = r_n;
    }
    int r_k = crsh_r_k[r_n];
    if (TwoD_hdr) {
      //Add as the header of rows of a given priority with this count
      int nx_r_n = crsh_r_pri_k_hdr[pri_v * (numCol + 1) + r_k];
      crsh_r_pri_k_hdr[pri_v * (numCol + 1) + r_k] = r_n;
      crsh_r_pri_k_lkb[r_n] = no_lk;
      crsh_r_pri_k_lkf[r_n] = nx_r_n;
      if (nx_r_n != no_lk)
	crsh_r_pri_k_lkb[nx_r_n] = r_n;
      //Update the minimum row count for this priority value
      crsh_r_pri_mn_r_k[pri_v] = min(r_k, crsh_r_pri_mn_r_k[pri_v]);
    }
    if (OneD_hdr) {
      //		Add as the header of the appropriate list of rows with a given count
      
      int nx_r_n = crsh_r_k_hdr[r_k];
      crsh_r_k_hdr[r_k] = r_n;
      crsh_r_k_lkb[r_n] = no_lk;
      crsh_r_k_lkf[r_n] = nx_r_n;
      if (nx_r_n != no_lk)
	crsh_r_k_lkb[nx_r_n] = r_n;
    }
  }
  //	Set up the row-wise matrix.
  //	Although not essential, it is convenient to avoid columns which
  //	are not active since the row counts determined above can be used.
  CrshARstart[0] = 0;
  for (int r_n = 0; r_n < numRow; r_n++) {
    CrshARstart[r_n + 1] = CrshARstart[r_n] + crsh_r_k[r_n];
  }
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_act_c[c_n] == crsh_vr_st_no_act)
      continue;
    //	Column is active so find its largest entry in updating the row counts
    crsh_mtx_c_mx_abs_v[c_n] = 0.0;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      //	    	printf("Col %d; Row %d: Activity %d;  CrshARstart = %d; |Value| = %g Max = %g\n", c_n, r_n, crsh_act_r[r_n],CrshARstart[r_n], abs(Avalue[el_n]), crsh_mtx_c_mx_abs_v[c_n]);
      crsh_mtx_c_mx_abs_v[c_n] = max(abs(Avalue[el_n]),
				     crsh_mtx_c_mx_abs_v[c_n]);
      if (crsh_act_r[r_n] == crsh_vr_st_no_act)
	continue;
      int r_el_n = CrshARstart[r_n];
      CrshARindex[r_el_n] = c_n;
      CrshARvalue[r_el_n] = Avalue[el_n];
      CrshARstart[r_n] += 1;
    }
  }
  for (int r_n = 0; r_n < numRow; r_n++) {
    CrshARstart[r_n] -= crsh_r_k[r_n];
  }
  //	printf("\nCol-wise matrix:\n");
  //	for (int c_n = 0; c_n < numCol; c_n++) {
  //		printf("Col %2d (%3d->-%3d): ", c_n, Astart[c_n], Astart[c_n + 1]-1);
  //		for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
  //			int r_n = Aindex[el_n];
  //			printf("%2d ", r_n);
  //		}
  //		printf("\n");
  //	}
  //	printf("\nRow-wise matrix:\n");
  //	for (int r_n = 0; r_n < numRow; r_n++) {
  //		printf("Row %2d (%3d->-%3d): ", r_n, CrshARstart[r_n], CrshARstart[r_n + 1]-1);
  //		for (int el_n = CrshARstart[r_n]; el_n < CrshARstart[r_n + 1]; el_n++) {
  //			int c_n = CrshARindex[el_n];
  //			printf("%2d ", c_n);
  //		}
  //		printf("\n");
  //	}
  //	printf("ltssf_iz_da complete\n");cout<<flush;
}

void HCrash::ltssf_ck_da() {
  bool er_fd = false;
  if (OneD_hdr) {
    for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++) {
      
      int r_n = crsh_r_pri_hdr[pri_v];
      //		printf("Priority %d: Header %d\n", pri_v, r_n);
      if (r_n == no_lk)
	continue;
      do {
	int r_ty = crsh_r_ty[r_n];
	//			printf("Row %d: Type %d\n", r_n, r_ty);
	int ck_pri_v = crsh_r_ty_pri_v[r_ty];
	if (ck_pri_v != pri_v) {
	  er_fd = true;
	  printf(
		 "ERROR: Row %d: Type %d has priority value error %d; %d\n",
		 r_n, r_ty, ck_pri_v, pri_v);
	}
	int nx_r_n = crsh_r_pri_lkf[r_n];
	if (nx_r_n != no_lk) {
	  if (crsh_r_pri_lkb[nx_r_n] != r_n) {
	    er_fd = true;
	    printf(
		   "ERROR: Row %d: Type %d has back link error %d; %d\n",
		   r_n, r_ty, nx_r_n, crsh_r_pri_lkf[r_n]);
	  }
	}
	r_n = nx_r_n;
      } while (r_n != no_lk);
    }
    for (int k = 0; k < numCol + 1; k++) {
      int r_n = crsh_r_k_hdr[k];
      //		printf("#Entries %d: Header %d\n", k, r_n);
      if (r_n != no_lk)
	do {
	  //			printf("Row %d: k %d\n", r_n, k);
	  int ck_k = 0;
	  for (int el_n = CrshARstart[r_n];
	       el_n < CrshARstart[r_n + 1]; el_n++) {
	    int c_n = CrshARindex[el_n];
	    if (crsh_act_c[c_n] != crsh_vr_st_no_act)
	      ck_k += 1;
	  }
	  if (ck_k != k) {
	    er_fd = true;
	    printf(
		   "ERROR: Row %d has number of entries error: True = %d; Updated = %d\n",
		   r_n, ck_k, k);
	  }
	  int nx_r_n = crsh_r_k_lkf[r_n];
	  if (nx_r_n != no_lk) {
	    if (crsh_r_k_lkb[nx_r_n] != r_n) {
	      er_fd = true;
	      printf(
		     "ERROR: Row %d: Count %d has back link error nx_r_n = %d; crsh_r_k_lkf[r_n] = %d\n",
		     r_n, k, nx_r_n, crsh_r_k_lkf[r_n]);
	    }
	  }
	  r_n = nx_r_n;
	} while (r_n != no_lk);
    }
  }
  if (TwoD_hdr) {
    for (int r_n = 0; r_n < numRow; r_n++) {
      if (crsh_act_r[r_n] == crsh_vr_st_no_act)
	continue;
      int k = crsh_r_k[r_n];
      int ck_k = 0;
      for (int el_n = CrshARstart[r_n]; el_n < CrshARstart[r_n + 1];
	   el_n++) {
	int c_n = CrshARindex[el_n];
	if (crsh_act_c[c_n] != crsh_vr_st_no_act)
	  ck_k += 1;
      }
      if (ck_k != k) {
	er_fd = true;
	printf(
	       "ERROR: Row %d has number of entries error: True = %d; Updated = %d\n",
	       r_n, ck_k, k);
      }
    }
    for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++) {
      int mn_r_k = numCol + 1;
      for (int k = 1; k < numCol + 1; k++) {
	//Check the rows with this priority value and count
	int r_n = crsh_r_pri_k_hdr[pri_v * (numCol + 1) + k];
	if (r_n == no_ix)
	  continue;
	do {
	  int ck_pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
	  if (ck_pri_v != pri_v) {
	    er_fd = true;
	    printf(
		   "ERROR: Row %d has ck_pri_v = %d but pri_v = %d\n",
		   r_n, ck_pri_v, pri_v);
	  }
	  int ck_k = crsh_r_k[r_n];
	  if (ck_k != k) {
	    er_fd = true;
	    printf("ERROR: Row %d has ck_k = %d but k = %d\n", r_n,
		   ck_k, k);
	  }
	  int nx_r_n = crsh_r_pri_k_lkf[r_n];
	  if (nx_r_n != no_lk) {
	    int prev_nx_r_n = crsh_r_pri_k_lkb[nx_r_n];
	    if (prev_nx_r_n != r_n) {
	      er_fd = true;
	      printf(
		     "ERROR: Back link error for nx_r_n = %d: prev_nx_r_n = %d but r_n = %d\n",
		     nx_r_n, prev_nx_r_n, r_n);
	    }
	  }
	  //					Update the true minimum row count
	  mn_r_k = min(k, mn_r_k);
	  r_n = nx_r_n;
	} while (r_n != no_lk);
      }
      if (crsh_r_pri_mn_r_k[pri_v] != mn_r_k)
	printf(
	       "ERROR: Priority %d has crsh_r_pri_mn_r_k = %d < %d = true mn_r_k\n",
	       pri_v, crsh_r_pri_mn_r_k[pri_v], mn_r_k);
      crsh_r_pri_mn_r_k[pri_v] = mn_r_k;
    }
  }
  //	if (!er_fd)	printf("No error in LTSSF data\n");
}

void HCrash::ltssf_rp_r_pri() {
  printf("\nRow       : ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", r_n);
  }
  printf(" \n");
  printf("crsh_r_ty:  ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_ty[r_n]);
  }
  printf(" \n");
  printf("pri_v:      ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_ty_pri_v[crsh_r_ty[r_n]]);
  }
  printf(" \n");
  printf("act_r:      ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_act_r[r_n]);
  }
  if (OneD_hdr) {
    printf(" \n");
    for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++) {
      int r_n = crsh_r_pri_hdr[pri_v];
      printf("Priority %d: Header %d: Rows ", pri_v, r_n);
      if (r_n == no_lk) {
	printf("\n");
	continue;
      }
      do {
	printf(" %2d: ", r_n);
	r_n = crsh_r_pri_lkf[r_n];
      } while (r_n != no_lk);
      printf("\n");
    }
  }
}

void HCrash::ltssf_rp_r_k() {
  printf("\n        : ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", r_n);
  }
  printf(" \n");
  printf("k: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_k[r_n]);
  }
  printf(" \n");
}


void HCrash::ltssf_rp_pri_k_da() {
  for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++) {
    printf("Pri = %1d: ", pri_v);
    for (int k = 1; k < numCol + 1; k++) {
      int hdr_ix = pri_v * (numCol + 1) + k;
      printf(" %2d", crsh_r_pri_k_hdr[hdr_ix]);
    }
    printf("\n");
  }
  printf("r_n: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", r_n);
  }
  printf("\n");
  printf("Act: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_act_r[r_n]);
  }
  printf("\n");
  printf("Cnt: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_k[r_n]);
  }
  printf("\n");
  printf("Lkf: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_pri_k_lkf[r_n]);
  }
  printf("\n");
  printf("Lkb: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_pri_k_lkb[r_n]);
  }
  printf("\n");
  for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++) {
    for (int k = 1; k < numCol + 1; k++) {
      int hdr_ix = pri_v * (numCol + 1) + k;
      int r_n = crsh_r_pri_k_hdr[hdr_ix];
      //			printf("Pri = %1d; k = %1d: hdr = %d\n", pri_v, k, r_n);
      if (r_n == no_lk)
	continue;
      int nx_r_n = crsh_r_pri_k_lkf[r_n];
      printf("Pri = %1d; k = %1d: %3d %3d <-  ", pri_v, k, r_n, nx_r_n);
      r_n = nx_r_n;
      if (r_n != no_lk) {
	do {
	  int prev_r_n = crsh_r_pri_k_lkb[r_n];
	  int nx_r_n = crsh_r_pri_k_lkf[r_n];
	  printf("  <-%3d %3d %3d ->  ", prev_r_n, r_n, nx_r_n);
	  r_n = nx_r_n;
	} while (r_n != no_lk);
      }
      printf("\n");
    }
  }
}

void HCrash::build_maxheap(double *heap_v, int *heap_i, int n) {
  int i;
  for(i = n/2; i >= 1; i--)
    {
      max_heapify(heap_v, heap_i, i, n);
    }
}

void HCrash::heapsort(double *heap_v, int *heap_i, int n) {
  double temp_v;
  int i, temp_i;
  for (i = n; i >= 2; i--)
    {
      temp_v = heap_v[i];
      heap_v[i] = heap_v[1];
      heap_v[1] = temp_v;
      temp_i = heap_i[i];
      heap_i[i] = heap_i[1];
      heap_i[1] = temp_i;
      max_heapify(heap_v, heap_i, 1, i - 1);
    }
}

void HCrash::max_heapify(double *heap_v, int *heap_i, int i, int n) {
  double temp_v;
  int j, temp_i;
  temp_v = heap_v[i];
  temp_i = heap_i[i];
  j = 2*i;
  while (j <= n)
    {
      if (j < n && heap_v[j+1] > heap_v[j])
	j = j+1;
      if (temp_v > heap_v[j])
	break;
      else if (temp_v <= heap_v[j])
        {
	  heap_v[j/2] = heap_v[j];
	  heap_i[j/2] = heap_i[j];
	  j = 2*j;
        }
    }
  heap_v[j/2] = temp_v;
  heap_i[j/2] = temp_i;
  return;
}
