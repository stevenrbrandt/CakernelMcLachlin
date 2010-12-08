/*  File produced by Kranc */

#define KRANC_C

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

void ML_ADMQuantities_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMQuantities::ML_Jadm","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMQuantities::ML_Jadm.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMQuantities::ML_Madm","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMQuantities::ML_Madm.");
  return;
}

void ML_ADMQuantities_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMQuantities_Body");
  }
  
  if (cctk_iteration % ML_ADMQuantities_calc_every != ML_ADMQuantities_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"grid::coordinates","Grid::coordinates","McLachlan::ML_curv","McLachlan::ML_Gamma","McLachlan::ML_lapse","McLachlan::ML_log_confac","McLachlan::ML_metric","McLachlan::ML_shift","McLachlan::ML_trace_curv","ML_ADMQuantities::ML_Jadm","ML_ADMQuantities::ML_Madm"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADMQuantities", 11, groups);
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  CCTK_REAL const dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL const dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL const dz = CCTK_DELTA_SPACE(2);
  int const di = 1;
  int const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  int const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  CCTK_REAL const dxi = 1.0 / dx;
  CCTK_REAL const dyi = 1.0 / dy;
  CCTK_REAL const dzi = 1.0 / dz;
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  CCTK_REAL const p1o12dx = INV(dx)/12.;
  CCTK_REAL const p1o12dy = INV(dy)/12.;
  CCTK_REAL const p1o12dz = INV(dz)/12.;
  CCTK_REAL const p1o144dxdy = (INV(dx)*INV(dy))/144.;
  CCTK_REAL const p1o144dxdz = (INV(dx)*INV(dz))/144.;
  CCTK_REAL const p1o144dydz = (INV(dy)*INV(dz))/144.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o12dx2 = -pow(dx,-2)/12.;
  CCTK_REAL const pm1o12dy2 = -pow(dy,-2)/12.;
  CCTK_REAL const pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_ADMQuantities,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth1trK = INITVALUE;
    // CCTK_REAL PDstandardNth2trK = INITVALUE;
    // CCTK_REAL PDstandardNth3trK = INITVALUE;
    // CCTK_REAL PDstandardNth1Xt1 = INITVALUE;
    // CCTK_REAL PDstandardNth2Xt1 = INITVALUE;
    // CCTK_REAL PDstandardNth3Xt1 = INITVALUE;
    // CCTK_REAL PDstandardNth1Xt2 = INITVALUE;
    // CCTK_REAL PDstandardNth2Xt2 = INITVALUE;
    // CCTK_REAL PDstandardNth3Xt2 = INITVALUE;
    // CCTK_REAL PDstandardNth1Xt3 = INITVALUE;
    // CCTK_REAL PDstandardNth2Xt3 = INITVALUE;
    // CCTK_REAL PDstandardNth3Xt3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL  alphaL = alpha[index];
    CCTK_REAL  At11L = At11[index];
    CCTK_REAL  At12L = At12[index];
    CCTK_REAL  At13L = At13[index];
    CCTK_REAL  At22L = At22[index];
    CCTK_REAL  At23L = At23[index];
    CCTK_REAL  At33L = At33[index];
    CCTK_REAL  beta1L = beta1[index];
    CCTK_REAL  beta2L = beta2[index];
    CCTK_REAL  beta3L = beta3[index];
    CCTK_REAL  eTttL = (*stress_energy_state) ? (eTtt[index]) : 0.0;
    CCTK_REAL  eTtxL = (*stress_energy_state) ? (eTtx[index]) : 0.0;
    CCTK_REAL  eTtyL = (*stress_energy_state) ? (eTty[index]) : 0.0;
    CCTK_REAL  eTtzL = (*stress_energy_state) ? (eTtz[index]) : 0.0;
    CCTK_REAL  eTxxL = (*stress_energy_state) ? (eTxx[index]) : 0.0;
    CCTK_REAL  eTxyL = (*stress_energy_state) ? (eTxy[index]) : 0.0;
    CCTK_REAL  eTxzL = (*stress_energy_state) ? (eTxz[index]) : 0.0;
    CCTK_REAL  eTyyL = (*stress_energy_state) ? (eTyy[index]) : 0.0;
    CCTK_REAL  eTyzL = (*stress_energy_state) ? (eTyz[index]) : 0.0;
    CCTK_REAL  eTzzL = (*stress_energy_state) ? (eTzz[index]) : 0.0;
    CCTK_REAL  gt11L = gt11[index];
    CCTK_REAL  gt12L = gt12[index];
    CCTK_REAL  gt13L = gt13[index];
    CCTK_REAL  gt22L = gt22[index];
    CCTK_REAL  gt23L = gt23[index];
    CCTK_REAL  gt33L = gt33[index];
    CCTK_REAL  phiL = phi[index];
    CCTK_REAL  trKL = trK[index];
    CCTK_REAL  xL = x[index];
    CCTK_REAL  Xt1L = Xt1[index];
    CCTK_REAL  Xt2L = Xt2[index];
    CCTK_REAL  Xt3L = Xt3[index];
    CCTK_REAL  yL = y[index];
    CCTK_REAL  zL = z[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    CCTK_REAL const PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    CCTK_REAL const PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    CCTK_REAL const PDstandardNth11gt11 = PDstandardNth11(gt11, i, j, k);
    CCTK_REAL const PDstandardNth22gt11 = PDstandardNth22(gt11, i, j, k);
    CCTK_REAL const PDstandardNth33gt11 = PDstandardNth33(gt11, i, j, k);
    CCTK_REAL const PDstandardNth12gt11 = PDstandardNth12(gt11, i, j, k);
    CCTK_REAL const PDstandardNth13gt11 = PDstandardNth13(gt11, i, j, k);
    CCTK_REAL const PDstandardNth23gt11 = PDstandardNth23(gt11, i, j, k);
    CCTK_REAL const PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    CCTK_REAL const PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    CCTK_REAL const PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    CCTK_REAL const PDstandardNth11gt12 = PDstandardNth11(gt12, i, j, k);
    CCTK_REAL const PDstandardNth22gt12 = PDstandardNth22(gt12, i, j, k);
    CCTK_REAL const PDstandardNth33gt12 = PDstandardNth33(gt12, i, j, k);
    CCTK_REAL const PDstandardNth12gt12 = PDstandardNth12(gt12, i, j, k);
    CCTK_REAL const PDstandardNth13gt12 = PDstandardNth13(gt12, i, j, k);
    CCTK_REAL const PDstandardNth23gt12 = PDstandardNth23(gt12, i, j, k);
    CCTK_REAL const PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    CCTK_REAL const PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    CCTK_REAL const PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    CCTK_REAL const PDstandardNth11gt13 = PDstandardNth11(gt13, i, j, k);
    CCTK_REAL const PDstandardNth22gt13 = PDstandardNth22(gt13, i, j, k);
    CCTK_REAL const PDstandardNth33gt13 = PDstandardNth33(gt13, i, j, k);
    CCTK_REAL const PDstandardNth12gt13 = PDstandardNth12(gt13, i, j, k);
    CCTK_REAL const PDstandardNth13gt13 = PDstandardNth13(gt13, i, j, k);
    CCTK_REAL const PDstandardNth23gt13 = PDstandardNth23(gt13, i, j, k);
    CCTK_REAL const PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    CCTK_REAL const PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    CCTK_REAL const PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    CCTK_REAL const PDstandardNth11gt22 = PDstandardNth11(gt22, i, j, k);
    CCTK_REAL const PDstandardNth22gt22 = PDstandardNth22(gt22, i, j, k);
    CCTK_REAL const PDstandardNth33gt22 = PDstandardNth33(gt22, i, j, k);
    CCTK_REAL const PDstandardNth12gt22 = PDstandardNth12(gt22, i, j, k);
    CCTK_REAL const PDstandardNth13gt22 = PDstandardNth13(gt22, i, j, k);
    CCTK_REAL const PDstandardNth23gt22 = PDstandardNth23(gt22, i, j, k);
    CCTK_REAL const PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    CCTK_REAL const PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    CCTK_REAL const PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    CCTK_REAL const PDstandardNth11gt23 = PDstandardNth11(gt23, i, j, k);
    CCTK_REAL const PDstandardNth22gt23 = PDstandardNth22(gt23, i, j, k);
    CCTK_REAL const PDstandardNth33gt23 = PDstandardNth33(gt23, i, j, k);
    CCTK_REAL const PDstandardNth12gt23 = PDstandardNth12(gt23, i, j, k);
    CCTK_REAL const PDstandardNth13gt23 = PDstandardNth13(gt23, i, j, k);
    CCTK_REAL const PDstandardNth23gt23 = PDstandardNth23(gt23, i, j, k);
    CCTK_REAL const PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    CCTK_REAL const PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    CCTK_REAL const PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    CCTK_REAL const PDstandardNth11gt33 = PDstandardNth11(gt33, i, j, k);
    CCTK_REAL const PDstandardNth22gt33 = PDstandardNth22(gt33, i, j, k);
    CCTK_REAL const PDstandardNth33gt33 = PDstandardNth33(gt33, i, j, k);
    CCTK_REAL const PDstandardNth12gt33 = PDstandardNth12(gt33, i, j, k);
    CCTK_REAL const PDstandardNth13gt33 = PDstandardNth13(gt33, i, j, k);
    CCTK_REAL const PDstandardNth23gt33 = PDstandardNth23(gt33, i, j, k);
    CCTK_REAL const PDstandardNth1trK = PDstandardNth1(trK, i, j, k);
    CCTK_REAL const PDstandardNth2trK = PDstandardNth2(trK, i, j, k);
    CCTK_REAL const PDstandardNth3trK = PDstandardNth3(trK, i, j, k);
    CCTK_REAL const PDstandardNth1Xt1 = PDstandardNth1(Xt1, i, j, k);
    CCTK_REAL const PDstandardNth2Xt1 = PDstandardNth2(Xt1, i, j, k);
    CCTK_REAL const PDstandardNth3Xt1 = PDstandardNth3(Xt1, i, j, k);
    CCTK_REAL const PDstandardNth1Xt2 = PDstandardNth1(Xt2, i, j, k);
    CCTK_REAL const PDstandardNth2Xt2 = PDstandardNth2(Xt2, i, j, k);
    CCTK_REAL const PDstandardNth3Xt2 = PDstandardNth3(Xt2, i, j, k);
    CCTK_REAL const PDstandardNth1Xt3 = PDstandardNth1(Xt3, i, j, k);
    CCTK_REAL const PDstandardNth2Xt3 = PDstandardNth2(Xt3, i, j, k);
    CCTK_REAL const PDstandardNth3Xt3 = PDstandardNth3(Xt3, i, j, k);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu21 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu31 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu32 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL dgtu111 = -2*(gtu11*gtu21*PDstandardNth1gt12 + 
      gtu11*gtu31*PDstandardNth1gt13 + gtu21*gtu31*PDstandardNth1gt23) - 
      PDstandardNth1gt11*SQR(gtu11) - PDstandardNth1gt22*SQR(gtu21) - 
      PDstandardNth1gt33*SQR(gtu31);
    
    CCTK_REAL dgtu211 = -(gtu11*gtu21*PDstandardNth1gt11) - 
      gtu11*gtu22*PDstandardNth1gt12 - gtu21*gtu31*PDstandardNth1gt13 - 
      gtu11*gtu32*PDstandardNth1gt13 - gtu21*gtu22*PDstandardNth1gt22 - 
      gtu22*gtu31*PDstandardNth1gt23 - gtu21*gtu32*PDstandardNth1gt23 - 
      gtu31*gtu32*PDstandardNth1gt33 - PDstandardNth1gt12*SQR(gtu21);
    
    CCTK_REAL dgtu311 = -(gtu11*gtu31*PDstandardNth1gt11) - 
      gtu21*gtu31*PDstandardNth1gt12 - gtu11*gtu32*PDstandardNth1gt12 - 
      gtu11*gtu33*PDstandardNth1gt13 - gtu21*gtu32*PDstandardNth1gt22 - 
      gtu31*gtu32*PDstandardNth1gt23 - gtu21*gtu33*PDstandardNth1gt23 - 
      gtu31*gtu33*PDstandardNth1gt33 - PDstandardNth1gt13*SQR(gtu31);
    
    CCTK_REAL dgtu221 = -2*(gtu21*gtu22*PDstandardNth1gt12 + 
      gtu21*gtu32*PDstandardNth1gt13 + gtu22*gtu32*PDstandardNth1gt23) - 
      PDstandardNth1gt11*SQR(gtu21) - PDstandardNth1gt22*SQR(gtu22) - 
      PDstandardNth1gt33*SQR(gtu32);
    
    CCTK_REAL dgtu321 = -(gtu21*gtu31*PDstandardNth1gt11) - 
      gtu22*gtu31*PDstandardNth1gt12 - gtu21*gtu32*PDstandardNth1gt12 - 
      gtu31*gtu32*PDstandardNth1gt13 - gtu21*gtu33*PDstandardNth1gt13 - 
      gtu22*gtu32*PDstandardNth1gt22 - gtu22*gtu33*PDstandardNth1gt23 - 
      gtu32*gtu33*PDstandardNth1gt33 - PDstandardNth1gt23*SQR(gtu32);
    
    CCTK_REAL dgtu331 = -2*(gtu31*gtu32*PDstandardNth1gt12 + 
      gtu31*gtu33*PDstandardNth1gt13 + gtu32*gtu33*PDstandardNth1gt23) - 
      PDstandardNth1gt11*SQR(gtu31) - PDstandardNth1gt22*SQR(gtu32) - 
      PDstandardNth1gt33*SQR(gtu33);
    
    CCTK_REAL dgtu112 = -2*(gtu11*gtu21*PDstandardNth2gt12 + 
      gtu11*gtu31*PDstandardNth2gt13 + gtu21*gtu31*PDstandardNth2gt23) - 
      PDstandardNth2gt11*SQR(gtu11) - PDstandardNth2gt22*SQR(gtu21) - 
      PDstandardNth2gt33*SQR(gtu31);
    
    CCTK_REAL dgtu212 = -(gtu11*gtu21*PDstandardNth2gt11) - 
      gtu11*gtu22*PDstandardNth2gt12 - gtu21*gtu31*PDstandardNth2gt13 - 
      gtu11*gtu32*PDstandardNth2gt13 - gtu21*gtu22*PDstandardNth2gt22 - 
      gtu22*gtu31*PDstandardNth2gt23 - gtu21*gtu32*PDstandardNth2gt23 - 
      gtu31*gtu32*PDstandardNth2gt33 - PDstandardNth2gt12*SQR(gtu21);
    
    CCTK_REAL dgtu312 = -(gtu11*gtu31*PDstandardNth2gt11) - 
      gtu21*gtu31*PDstandardNth2gt12 - gtu11*gtu32*PDstandardNth2gt12 - 
      gtu11*gtu33*PDstandardNth2gt13 - gtu21*gtu32*PDstandardNth2gt22 - 
      gtu31*gtu32*PDstandardNth2gt23 - gtu21*gtu33*PDstandardNth2gt23 - 
      gtu31*gtu33*PDstandardNth2gt33 - PDstandardNth2gt13*SQR(gtu31);
    
    CCTK_REAL dgtu222 = -2*(gtu21*gtu22*PDstandardNth2gt12 + 
      gtu21*gtu32*PDstandardNth2gt13 + gtu22*gtu32*PDstandardNth2gt23) - 
      PDstandardNth2gt11*SQR(gtu21) - PDstandardNth2gt22*SQR(gtu22) - 
      PDstandardNth2gt33*SQR(gtu32);
    
    CCTK_REAL dgtu322 = -(gtu21*gtu31*PDstandardNth2gt11) - 
      gtu22*gtu31*PDstandardNth2gt12 - gtu21*gtu32*PDstandardNth2gt12 - 
      gtu31*gtu32*PDstandardNth2gt13 - gtu21*gtu33*PDstandardNth2gt13 - 
      gtu22*gtu32*PDstandardNth2gt22 - gtu22*gtu33*PDstandardNth2gt23 - 
      gtu32*gtu33*PDstandardNth2gt33 - PDstandardNth2gt23*SQR(gtu32);
    
    CCTK_REAL dgtu332 = -2*(gtu31*gtu32*PDstandardNth2gt12 + 
      gtu31*gtu33*PDstandardNth2gt13 + gtu32*gtu33*PDstandardNth2gt23) - 
      PDstandardNth2gt11*SQR(gtu31) - PDstandardNth2gt22*SQR(gtu32) - 
      PDstandardNth2gt33*SQR(gtu33);
    
    CCTK_REAL dgtu113 = -2*(gtu11*gtu21*PDstandardNth3gt12 + 
      gtu11*gtu31*PDstandardNth3gt13 + gtu21*gtu31*PDstandardNth3gt23) - 
      PDstandardNth3gt11*SQR(gtu11) - PDstandardNth3gt22*SQR(gtu21) - 
      PDstandardNth3gt33*SQR(gtu31);
    
    CCTK_REAL dgtu213 = -(gtu11*gtu21*PDstandardNth3gt11) - 
      gtu11*gtu22*PDstandardNth3gt12 - gtu21*gtu31*PDstandardNth3gt13 - 
      gtu11*gtu32*PDstandardNth3gt13 - gtu21*gtu22*PDstandardNth3gt22 - 
      gtu22*gtu31*PDstandardNth3gt23 - gtu21*gtu32*PDstandardNth3gt23 - 
      gtu31*gtu32*PDstandardNth3gt33 - PDstandardNth3gt12*SQR(gtu21);
    
    CCTK_REAL dgtu313 = -(gtu11*gtu31*PDstandardNth3gt11) - 
      gtu21*gtu31*PDstandardNth3gt12 - gtu11*gtu32*PDstandardNth3gt12 - 
      gtu11*gtu33*PDstandardNth3gt13 - gtu21*gtu32*PDstandardNth3gt22 - 
      gtu31*gtu32*PDstandardNth3gt23 - gtu21*gtu33*PDstandardNth3gt23 - 
      gtu31*gtu33*PDstandardNth3gt33 - PDstandardNth3gt13*SQR(gtu31);
    
    CCTK_REAL dgtu223 = -2*(gtu21*gtu22*PDstandardNth3gt12 + 
      gtu21*gtu32*PDstandardNth3gt13 + gtu22*gtu32*PDstandardNth3gt23) - 
      PDstandardNth3gt11*SQR(gtu21) - PDstandardNth3gt22*SQR(gtu22) - 
      PDstandardNth3gt33*SQR(gtu32);
    
    CCTK_REAL dgtu323 = -(gtu21*gtu31*PDstandardNth3gt11) - 
      gtu22*gtu31*PDstandardNth3gt12 - gtu21*gtu32*PDstandardNth3gt12 - 
      gtu31*gtu32*PDstandardNth3gt13 - gtu21*gtu33*PDstandardNth3gt13 - 
      gtu22*gtu32*PDstandardNth3gt22 - gtu22*gtu33*PDstandardNth3gt23 - 
      gtu32*gtu33*PDstandardNth3gt33 - PDstandardNth3gt23*SQR(gtu32);
    
    CCTK_REAL dgtu333 = -2*(gtu31*gtu32*PDstandardNth3gt12 + 
      gtu31*gtu33*PDstandardNth3gt13 + gtu32*gtu33*PDstandardNth3gt23) - 
      PDstandardNth3gt11*SQR(gtu31) - PDstandardNth3gt22*SQR(gtu32) - 
      PDstandardNth3gt33*SQR(gtu33);
    
    CCTK_REAL Gtl111 = khalf*PDstandardNth1gt11;
    
    CCTK_REAL Gtl112 = khalf*PDstandardNth2gt11;
    
    CCTK_REAL Gtl113 = khalf*PDstandardNth3gt11;
    
    CCTK_REAL Gtl122 = -(khalf*PDstandardNth1gt22) + PDstandardNth2gt12;
    
    CCTK_REAL Gtl123 = khalf*(-PDstandardNth1gt23 + PDstandardNth2gt13 + 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl133 = -(khalf*PDstandardNth1gt33) + PDstandardNth3gt13;
    
    CCTK_REAL Gtl211 = PDstandardNth1gt12 - khalf*PDstandardNth2gt11;
    
    CCTK_REAL Gtl212 = khalf*PDstandardNth1gt22;
    
    CCTK_REAL Gtl213 = khalf*(PDstandardNth1gt23 - PDstandardNth2gt13 + 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl222 = khalf*PDstandardNth2gt22;
    
    CCTK_REAL Gtl223 = khalf*PDstandardNth3gt22;
    
    CCTK_REAL Gtl233 = -(khalf*PDstandardNth2gt33) + PDstandardNth3gt23;
    
    CCTK_REAL Gtl311 = PDstandardNth1gt13 - khalf*PDstandardNth3gt11;
    
    CCTK_REAL Gtl312 = khalf*(PDstandardNth1gt23 + PDstandardNth2gt13 - 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl313 = khalf*PDstandardNth1gt33;
    
    CCTK_REAL Gtl322 = PDstandardNth2gt23 - khalf*PDstandardNth3gt22;
    
    CCTK_REAL Gtl323 = khalf*PDstandardNth2gt33;
    
    CCTK_REAL Gtl333 = khalf*PDstandardNth3gt33;
    
    CCTK_REAL Gtlu111 = Gtl111*gtu11 + Gtl112*gtu21 + Gtl113*gtu31;
    
    CCTK_REAL Gtlu112 = Gtl111*gtu21 + Gtl112*gtu22 + Gtl113*gtu32;
    
    CCTK_REAL Gtlu113 = Gtl111*gtu31 + Gtl112*gtu32 + Gtl113*gtu33;
    
    CCTK_REAL Gtlu121 = Gtl112*gtu11 + Gtl122*gtu21 + Gtl123*gtu31;
    
    CCTK_REAL Gtlu122 = Gtl112*gtu21 + Gtl122*gtu22 + Gtl123*gtu32;
    
    CCTK_REAL Gtlu123 = Gtl112*gtu31 + Gtl122*gtu32 + Gtl123*gtu33;
    
    CCTK_REAL Gtlu131 = Gtl113*gtu11 + Gtl123*gtu21 + Gtl133*gtu31;
    
    CCTK_REAL Gtlu132 = Gtl113*gtu21 + Gtl123*gtu22 + Gtl133*gtu32;
    
    CCTK_REAL Gtlu133 = Gtl113*gtu31 + Gtl123*gtu32 + Gtl133*gtu33;
    
    CCTK_REAL Gtlu211 = Gtl211*gtu11 + Gtl212*gtu21 + Gtl213*gtu31;
    
    CCTK_REAL Gtlu212 = Gtl211*gtu21 + Gtl212*gtu22 + Gtl213*gtu32;
    
    CCTK_REAL Gtlu213 = Gtl211*gtu31 + Gtl212*gtu32 + Gtl213*gtu33;
    
    CCTK_REAL Gtlu221 = Gtl212*gtu11 + Gtl222*gtu21 + Gtl223*gtu31;
    
    CCTK_REAL Gtlu222 = Gtl212*gtu21 + Gtl222*gtu22 + Gtl223*gtu32;
    
    CCTK_REAL Gtlu223 = Gtl212*gtu31 + Gtl222*gtu32 + Gtl223*gtu33;
    
    CCTK_REAL Gtlu231 = Gtl213*gtu11 + Gtl223*gtu21 + Gtl233*gtu31;
    
    CCTK_REAL Gtlu232 = Gtl213*gtu21 + Gtl223*gtu22 + Gtl233*gtu32;
    
    CCTK_REAL Gtlu233 = Gtl213*gtu31 + Gtl223*gtu32 + Gtl233*gtu33;
    
    CCTK_REAL Gtlu311 = Gtl311*gtu11 + Gtl312*gtu21 + Gtl313*gtu31;
    
    CCTK_REAL Gtlu312 = Gtl311*gtu21 + Gtl312*gtu22 + Gtl313*gtu32;
    
    CCTK_REAL Gtlu313 = Gtl311*gtu31 + Gtl312*gtu32 + Gtl313*gtu33;
    
    CCTK_REAL Gtlu321 = Gtl312*gtu11 + Gtl322*gtu21 + Gtl323*gtu31;
    
    CCTK_REAL Gtlu322 = Gtl312*gtu21 + Gtl322*gtu22 + Gtl323*gtu32;
    
    CCTK_REAL Gtlu323 = Gtl312*gtu31 + Gtl322*gtu32 + Gtl323*gtu33;
    
    CCTK_REAL Gtlu331 = Gtl313*gtu11 + Gtl323*gtu21 + Gtl333*gtu31;
    
    CCTK_REAL Gtlu332 = Gtl313*gtu21 + Gtl323*gtu22 + Gtl333*gtu32;
    
    CCTK_REAL Gtlu333 = Gtl313*gtu31 + Gtl323*gtu32 + Gtl333*gtu33;
    
    CCTK_REAL Gt111 = Gtl111*gtu11 + Gtl211*gtu21 + Gtl311*gtu31;
    
    CCTK_REAL Gt211 = Gtl111*gtu21 + Gtl211*gtu22 + Gtl311*gtu32;
    
    CCTK_REAL Gt311 = Gtl111*gtu31 + Gtl211*gtu32 + Gtl311*gtu33;
    
    CCTK_REAL Gt112 = Gtl112*gtu11 + Gtl212*gtu21 + Gtl312*gtu31;
    
    CCTK_REAL Gt212 = Gtl112*gtu21 + Gtl212*gtu22 + Gtl312*gtu32;
    
    CCTK_REAL Gt312 = Gtl112*gtu31 + Gtl212*gtu32 + Gtl312*gtu33;
    
    CCTK_REAL Gt113 = Gtl113*gtu11 + Gtl213*gtu21 + Gtl313*gtu31;
    
    CCTK_REAL Gt213 = Gtl113*gtu21 + Gtl213*gtu22 + Gtl313*gtu32;
    
    CCTK_REAL Gt313 = Gtl113*gtu31 + Gtl213*gtu32 + Gtl313*gtu33;
    
    CCTK_REAL Gt122 = Gtl122*gtu11 + Gtl222*gtu21 + Gtl322*gtu31;
    
    CCTK_REAL Gt222 = Gtl122*gtu21 + Gtl222*gtu22 + Gtl322*gtu32;
    
    CCTK_REAL Gt322 = Gtl122*gtu31 + Gtl222*gtu32 + Gtl322*gtu33;
    
    CCTK_REAL Gt123 = Gtl123*gtu11 + Gtl223*gtu21 + Gtl323*gtu31;
    
    CCTK_REAL Gt223 = Gtl123*gtu21 + Gtl223*gtu22 + Gtl323*gtu32;
    
    CCTK_REAL Gt323 = Gtl123*gtu31 + Gtl223*gtu32 + Gtl323*gtu33;
    
    CCTK_REAL Gt133 = Gtl133*gtu11 + Gtl233*gtu21 + Gtl333*gtu31;
    
    CCTK_REAL Gt233 = Gtl133*gtu21 + Gtl233*gtu22 + Gtl333*gtu32;
    
    CCTK_REAL Gt333 = Gtl133*gtu31 + Gtl233*gtu32 + Gtl333*gtu33;
    
    CCTK_REAL Xtn1 = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + 
      Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    CCTK_REAL Xtn2 = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + 
      Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    CCTK_REAL Xtn3 = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + 
      Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    CCTK_REAL Rt11 = 3*(Gt111*Gtlu111 + Gt112*Gtlu112 + Gt113*Gtlu113) + 
      2*(Gt211*Gtlu121 + Gt212*Gtlu122 + Gt213*Gtlu123 + Gt311*Gtlu131 + 
      Gt312*Gtlu132 + Gt313*Gtlu133) + Gt211*Gtlu211 + Gt212*Gtlu212 + 
      Gt213*Gtlu213 + Gt311*Gtlu311 + Gt312*Gtlu312 + Gt313*Gtlu313 + 
      gt11L*PDstandardNth1Xt1 + gt12L*PDstandardNth1Xt2 + 
      gt13L*PDstandardNth1Xt3 + khalf*(-(gtu11*PDstandardNth11gt11) - 
      2*gtu21*PDstandardNth12gt11 - 2*gtu31*PDstandardNth13gt11 - 
      gtu22*PDstandardNth22gt11 - 2*gtu32*PDstandardNth23gt11 - 
      gtu33*PDstandardNth33gt11) + Gtl111*Xtn1 + Gtl112*Xtn2 + Gtl113*Xtn3;
    
    CCTK_REAL Rt12 = khalf*(4*(Gt211*Gtlu221 + Gt212*Gtlu222 + 
      Gt213*Gtlu223) + 2*(Gt112*Gtlu111 + Gt122*Gtlu112 + Gt123*Gtlu113 + 
      Gt111*Gtlu121 + Gt212*Gtlu121 + Gt112*Gtlu122 + Gt222*Gtlu122 + 
      Gt113*Gtlu123 + Gt223*Gtlu123 + Gt312*Gtlu131 + Gt322*Gtlu132 + 
      Gt323*Gtlu133 + Gt111*Gtlu211 + Gt112*Gtlu212 + Gt113*Gtlu213 + 
      Gt311*Gtlu231 + Gt312*Gtlu232 + Gt313*Gtlu233 + Gt311*Gtlu321 + 
      Gt312*Gtlu322 + Gt313*Gtlu323) - gtu11*PDstandardNth11gt12 - 
      2*gtu21*PDstandardNth12gt12 - 2*gtu31*PDstandardNth13gt12 + 
      gt12L*PDstandardNth1Xt1 + gt22L*PDstandardNth1Xt2 + 
      gt23L*PDstandardNth1Xt3 - gtu22*PDstandardNth22gt12 - 
      2*gtu32*PDstandardNth23gt12 + gt11L*PDstandardNth2Xt1 + 
      gt12L*PDstandardNth2Xt2 + gt13L*PDstandardNth2Xt3 - 
      gtu33*PDstandardNth33gt12 + Gtl112*Xtn1 + Gtl211*Xtn1 + Gtl122*Xtn2 + 
      Gtl212*Xtn2 + Gtl123*Xtn3 + Gtl213*Xtn3);
    
    CCTK_REAL Rt13 = khalf*(2*(Gt113*Gtlu111 + Gt123*Gtlu112 + 
      Gt133*Gtlu113 + Gt213*Gtlu121 + Gt223*Gtlu122 + Gt233*Gtlu123 + 
      Gt111*Gtlu131 + Gt313*Gtlu131 + Gt112*Gtlu132 + Gt323*Gtlu132 + 
      Gt113*Gtlu133 + Gt333*Gtlu133 + Gt211*Gtlu231 + Gt212*Gtlu232 + 
      Gt213*Gtlu233 + Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*Gtlu313 + 
      Gt211*Gtlu321 + Gt212*Gtlu322 + Gt213*Gtlu323) + 4*(Gt311*Gtlu331 + 
      Gt312*Gtlu332 + Gt313*Gtlu333) - gtu11*PDstandardNth11gt13 - 
      2*gtu21*PDstandardNth12gt13 - 2*gtu31*PDstandardNth13gt13 + 
      gt13L*PDstandardNth1Xt1 + gt23L*PDstandardNth1Xt2 + 
      gt33L*PDstandardNth1Xt3 - gtu22*PDstandardNth22gt13 - 
      2*gtu32*PDstandardNth23gt13 - gtu33*PDstandardNth33gt13 + 
      gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + 
      gt13L*PDstandardNth3Xt3 + Gtl113*Xtn1 + Gtl311*Xtn1 + Gtl123*Xtn2 + 
      Gtl312*Xtn2 + Gtl133*Xtn3 + Gtl313*Xtn3);
    
    CCTK_REAL Rt22 = Gt112*(Gtlu121 + 2*Gtlu211) + Gt122*(Gtlu122 + 
      2*Gtlu212) + Gt123*(Gtlu123 + 2*Gtlu213) + 3*(Gt212*Gtlu221 + 
      Gt222*Gtlu222 + Gt223*Gtlu223) + 2*(Gt312*Gtlu231 + Gt322*Gtlu232 + 
      Gt323*Gtlu233) + Gt312*Gtlu321 + Gt322*Gtlu322 + Gt323*Gtlu323 + 
      gt12L*PDstandardNth2Xt1 + gt22L*PDstandardNth2Xt2 + 
      gt23L*PDstandardNth2Xt3 + khalf*(-(gtu11*PDstandardNth11gt22) - 
      2*gtu21*PDstandardNth12gt22 - 2*gtu31*PDstandardNth13gt22 - 
      gtu22*PDstandardNth22gt22 - 2*gtu32*PDstandardNth23gt22 - 
      gtu33*PDstandardNth33gt22) + Gtl212*Xtn1 + Gtl222*Xtn2 + Gtl223*Xtn3;
    
    CCTK_REAL Rt23 = khalf*(2*(Gt112*Gtlu131 + Gt122*Gtlu132 + 
      Gt123*Gtlu133 + Gt113*Gtlu211 + Gt123*Gtlu212 + Gt133*Gtlu213 + 
      Gt213*Gtlu221 + Gt223*Gtlu222 + Gt233*Gtlu223 + Gt212*Gtlu231 + 
      Gt313*Gtlu231 + Gt222*Gtlu232 + Gt323*Gtlu232 + Gt223*Gtlu233 + 
      Gt333*Gtlu233 + Gt112*Gtlu311 + Gt122*Gtlu312 + Gt123*Gtlu313 + 
      Gt212*Gtlu321 + Gt222*Gtlu322 + Gt223*Gtlu323) + 4*(Gt312*Gtlu331 + 
      Gt322*Gtlu332 + Gt323*Gtlu333) - gtu11*PDstandardNth11gt23 - 
      2*gtu21*PDstandardNth12gt23 - 2*gtu31*PDstandardNth13gt23 - 
      gtu22*PDstandardNth22gt23 - 2*gtu32*PDstandardNth23gt23 + 
      gt13L*PDstandardNth2Xt1 + gt23L*PDstandardNth2Xt2 + 
      gt33L*PDstandardNth2Xt3 - gtu33*PDstandardNth33gt23 + 
      gt12L*PDstandardNth3Xt1 + gt22L*PDstandardNth3Xt2 + 
      gt23L*PDstandardNth3Xt3 + Gtl213*Xtn1 + Gtl312*Xtn1 + Gtl223*Xtn2 + 
      Gtl322*Xtn2 + Gtl233*Xtn3 + Gtl323*Xtn3);
    
    CCTK_REAL Rt33 = Gt113*(Gtlu131 + 2*Gtlu311) + Gt123*(Gtlu132 + 
      2*Gtlu312) + Gt133*(Gtlu133 + 2*Gtlu313) + Gt213*(Gtlu231 + 2*Gtlu321) 
      + Gt223*(Gtlu232 + 2*Gtlu322) + Gt233*(Gtlu233 + 2*Gtlu323) + 
      3*(Gt313*Gtlu331 + Gt323*Gtlu332 + Gt333*Gtlu333) + 
      khalf*(-(gtu11*PDstandardNth11gt33) - 2*gtu21*PDstandardNth12gt33 - 
      2*gtu31*PDstandardNth13gt33 - gtu22*PDstandardNth22gt33 - 
      2*gtu32*PDstandardNth23gt33 - gtu33*PDstandardNth33gt33) + 
      gt13L*PDstandardNth3Xt1 + gt23L*PDstandardNth3Xt2 + 
      gt33L*PDstandardNth3Xt3 + Gtl313*Xtn1 + Gtl323*Xtn2 + Gtl333*Xtn3;
    
    CCTK_REAL trRt = gtu11*Rt11 + gtu22*Rt22 + 2*(gtu21*Rt12 + gtu31*Rt13 
      + gtu32*Rt23) + gtu33*Rt33;
    
    CCTK_REAL ephi = IfThen(conformalMethod,pow(phiL,-khalf),exp(phiL));
    
    CCTK_REAL Atm11 = At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    CCTK_REAL Atm21 = At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    CCTK_REAL Atm31 = At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    CCTK_REAL Atm12 = At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    CCTK_REAL Atm22 = At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    CCTK_REAL Atm32 = At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    CCTK_REAL Atm13 = At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    CCTK_REAL Atm23 = At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    CCTK_REAL Atm33 = At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    CCTK_REAL rho = pow(alphaL,-2)*(eTttL - 2*(beta2L*eTtyL + 
      beta3L*eTtzL) + 2*(beta1L*(-eTtxL + beta2L*eTxyL + beta3L*eTxzL) + 
      beta2L*beta3L*eTyzL) + eTxxL*SQR(beta1L) + eTyyL*SQR(beta2L) + 
      eTzzL*SQR(beta3L));
    
    CCTK_REAL S1 = (-eTtxL + beta1L*eTxxL + beta2L*eTxyL + 
      beta3L*eTxzL)*INV(alphaL);
    
    CCTK_REAL S2 = (-eTtyL + beta1L*eTxyL + beta2L*eTyyL + 
      beta3L*eTyzL)*INV(alphaL);
    
    CCTK_REAL S3 = (-eTtzL + beta1L*eTxzL + beta2L*eTyzL + 
      beta3L*eTzzL)*INV(alphaL);
    
    CCTK_REAL MadmL = 
      0.01989436788648691697111047042156429525431*(-((Gt111*Gtlu111 + 
      Gt112*Gtlu112 + Gt113*Gtlu113 + Gt211*Gtlu121 + Gt212*Gtlu122 + 
      Gt213*Gtlu123 + Gt311*Gtlu131 + Gt312*Gtlu132 + Gt313*Gtlu133)*gtu11) - 
      (Gt112*Gtlu111 + Gt122*Gtlu112 + Gt123*Gtlu113 + Gt212*Gtlu121 + 
      Gt222*Gtlu122 + Gt223*Gtlu123 + Gt312*Gtlu131 + Gt322*Gtlu132 + 
      Gt323*Gtlu133)*gtu21 - (Gt111*Gtlu211 + Gt112*Gtlu212 + Gt113*Gtlu213 + 
      Gt211*Gtlu221 + Gt212*Gtlu222 + Gt213*Gtlu223 + Gt311*Gtlu231 + 
      Gt312*Gtlu232 + Gt313*Gtlu233)*gtu21 - (Gt112*Gtlu211 + Gt122*Gtlu212 + 
      Gt123*Gtlu213 + Gt212*Gtlu221 + Gt222*Gtlu222 + Gt223*Gtlu223 + 
      Gt312*Gtlu231 + Gt322*Gtlu232 + Gt323*Gtlu233)*gtu22 - (Gt113*Gtlu111 + 
      Gt123*Gtlu112 + Gt133*Gtlu113 + Gt213*Gtlu121 + Gt223*Gtlu122 + 
      Gt233*Gtlu123 + Gt313*Gtlu131 + Gt323*Gtlu132 + Gt333*Gtlu133)*gtu31 - 
      (Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*Gtlu313 + Gt211*Gtlu321 + 
      Gt212*Gtlu322 + Gt213*Gtlu323 + Gt311*Gtlu331 + Gt312*Gtlu332 + 
      Gt313*Gtlu333)*gtu31 - (Gt113*Gtlu211 + Gt123*Gtlu212 + Gt133*Gtlu213 + 
      Gt213*Gtlu221 + Gt223*Gtlu222 + Gt233*Gtlu223 + Gt313*Gtlu231 + 
      Gt323*Gtlu232 + Gt333*Gtlu233)*gtu32 - (Gt112*Gtlu311 + Gt122*Gtlu312 + 
      Gt123*Gtlu313 + Gt212*Gtlu321 + Gt222*Gtlu322 + Gt223*Gtlu323 + 
      Gt312*Gtlu331 + Gt322*Gtlu332 + Gt323*Gtlu333)*gtu32 - (Gt113*Gtlu311 + 
      Gt123*Gtlu312 + Gt133*Gtlu313 + Gt213*Gtlu321 + Gt223*Gtlu322 + 
      Gt233*Gtlu323 + Gt313*Gtlu331 + Gt323*Gtlu332 + Gt333*Gtlu333)*gtu33 + 
      trRt - ephi*trRt + pow(ephi,5)*(2*Atm12*Atm21 + 2.*Atm13*Atm31 + 
      2.*Atm23*Atm32 + 50.26548245743669181540229413247204614715*rho + 
      SQR(Atm11) + SQR(Atm22) + SQR(Atm33) - 
      0.6666666666666666666666666666666666666667*SQR(trKL)));
    
    CCTK_REAL Jadm1L = 
      0.01989436788648691697111047042156429525431*(2*Atm23 - 2*Atm32 + 
      (-(At11L*dgtu113) - 2*At12L*dgtu213 - At22L*dgtu223 - 2*At13L*dgtu313 - 
      2*At23L*dgtu323 - At33L*dgtu333 + kfourthird*PDstandardNth3trK + 
      50.26548245743669181540229413247204614715*S3)*yL + (At11L*dgtu112 + 
      At22L*dgtu222 + 2*(At12L*dgtu212 + At13L*dgtu312 + At23L*dgtu322) + 
      At33L*dgtu332 - kfourthird*PDstandardNth2trK - 
      50.26548245743669181540229413247204614715*S2)*zL)*pow(ephi,6);
    
    CCTK_REAL Jadm2L = 
      0.01989436788648691697111047042156429525431*(-2*Atm13 + 2*Atm31 + 
      (At11L*dgtu113 + At22L*dgtu223 + 2*(At12L*dgtu213 + At13L*dgtu313 + 
      At23L*dgtu323) + At33L*dgtu333 - kfourthird*PDstandardNth3trK - 
      50.26548245743669181540229413247204614715*S3)*xL + (-(At11L*dgtu111) - 
      2*At12L*dgtu211 - At22L*dgtu221 - 2*At13L*dgtu311 - 2*At23L*dgtu321 - 
      At33L*dgtu331 + kfourthird*PDstandardNth1trK + 
      50.26548245743669181540229413247204614715*S1)*zL)*pow(ephi,6);
    
    CCTK_REAL Jadm3L = 
      0.01989436788648691697111047042156429525431*(2*Atm12 - 2*Atm21 + 
      (-(At11L*dgtu112) - 2*At12L*dgtu212 - At22L*dgtu222 - 2*At13L*dgtu312 - 
      2*At23L*dgtu322 - At33L*dgtu332 + kfourthird*PDstandardNth2trK + 
      50.26548245743669181540229413247204614715*S2)*xL + (At11L*dgtu111 + 
      At22L*dgtu221 + 2*(At12L*dgtu211 + At13L*dgtu311 + At23L*dgtu321) + 
      At33L*dgtu331 - kfourthird*PDstandardNth1trK - 
      50.26548245743669181540229413247204614715*S1)*yL)*pow(ephi,6);
    
    
    /* Copy local copies back to grid functions */
    Jadm1[index] = Jadm1L;
    Jadm2[index] = Jadm2L;
    Jadm3[index] = Jadm3L;
    Madm[index] = MadmL;
  }
  LC_ENDLOOP3 (ML_ADMQuantities);
}

void ML_ADMQuantities(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADMQuantities_Body);
}
