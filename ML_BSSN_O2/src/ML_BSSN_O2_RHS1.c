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

void ML_BSSN_O2_RHS1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O2_RHS1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O2_RHS1_calc_every != ML_BSSN_O2_RHS1_calc_offset)
  {
    return;
  }
  
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
  CCTK_REAL const p1o2dx = khalf*INV(dx);
  CCTK_REAL const p1o2dy = khalf*INV(dy);
  CCTK_REAL const p1o2dz = khalf*INV(dz);
  CCTK_REAL const p1o4dxdy = (INV(dx)*INV(dy))/4.;
  CCTK_REAL const p1o4dxdz = (INV(dx)*INV(dz))/4.;
  CCTK_REAL const p1o4dydz = (INV(dy)*INV(dz))/4.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1odx2 = pow(dx,-2);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1ody2 = pow(dy,-2);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const p1odz2 = pow(dz,-2);
  CCTK_REAL const pm1o2dx = -(khalf*INV(dx));
  CCTK_REAL const pm1o2dy = -(khalf*INV(dy));
  CCTK_REAL const pm1o2dz = -(khalf*INV(dz));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_O2_RHS1,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth1alpha = INITVALUE;
    // CCTK_REAL PDstandardNth2alpha = INITVALUE;
    // CCTK_REAL PDstandardNth3alpha = INITVALUE;
    // CCTK_REAL PDstandardNth11alpha = INITVALUE;
    // CCTK_REAL PDstandardNth22alpha = INITVALUE;
    // CCTK_REAL PDstandardNth33alpha = INITVALUE;
    // CCTK_REAL PDstandardNth12alpha = INITVALUE;
    // CCTK_REAL PDstandardNth13alpha = INITVALUE;
    // CCTK_REAL PDstandardNth23alpha = INITVALUE;
    // CCTK_REAL PDstandardNth1beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth11beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth22beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth33beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth12beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth13beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth23beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth1beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth11beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth22beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth33beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth12beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth13beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth23beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth1beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth11beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth22beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth33beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth12beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth13beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth23beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth1phi = INITVALUE;
    // CCTK_REAL PDstandardNth2phi = INITVALUE;
    // CCTK_REAL PDstandardNth3phi = INITVALUE;
    // CCTK_REAL PDstandardNth1trK = INITVALUE;
    // CCTK_REAL PDstandardNth2trK = INITVALUE;
    // CCTK_REAL PDstandardNth3trK = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL  AL = A[index];
    CCTK_REAL  alphaL = alpha[index];
    CCTK_REAL  At11L = At11[index];
    CCTK_REAL  At12L = At12[index];
    CCTK_REAL  At13L = At13[index];
    CCTK_REAL  At22L = At22[index];
    CCTK_REAL  At23L = At23[index];
    CCTK_REAL  At33L = At33[index];
    CCTK_REAL  B1L = B1[index];
    CCTK_REAL  B2L = B2[index];
    CCTK_REAL  B3L = B3[index];
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
    CCTK_REAL  rL = r[index];
    CCTK_REAL  trKL = trK[index];
    CCTK_REAL  trKrhsL = trKrhs[index];
    CCTK_REAL  Xt1L = Xt1[index];
    CCTK_REAL  Xt1rhsL = Xt1rhs[index];
    CCTK_REAL  Xt2L = Xt2[index];
    CCTK_REAL  Xt2rhsL = Xt2rhs[index];
    CCTK_REAL  Xt3L = Xt3[index];
    CCTK_REAL  Xt3rhsL = Xt3rhs[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandardNth1alpha = PDstandardNth1(alpha, i, j, k);
    CCTK_REAL const PDstandardNth2alpha = PDstandardNth2(alpha, i, j, k);
    CCTK_REAL const PDstandardNth3alpha = PDstandardNth3(alpha, i, j, k);
    CCTK_REAL const PDstandardNth11alpha = PDstandardNth11(alpha, i, j, k);
    CCTK_REAL const PDstandardNth22alpha = PDstandardNth22(alpha, i, j, k);
    CCTK_REAL const PDstandardNth33alpha = PDstandardNth33(alpha, i, j, k);
    CCTK_REAL const PDstandardNth12alpha = PDstandardNth12(alpha, i, j, k);
    CCTK_REAL const PDstandardNth13alpha = PDstandardNth13(alpha, i, j, k);
    CCTK_REAL const PDstandardNth23alpha = PDstandardNth23(alpha, i, j, k);
    CCTK_REAL const PDstandardNth1beta1 = PDstandardNth1(beta1, i, j, k);
    CCTK_REAL const PDstandardNth2beta1 = PDstandardNth2(beta1, i, j, k);
    CCTK_REAL const PDstandardNth3beta1 = PDstandardNth3(beta1, i, j, k);
    CCTK_REAL const PDstandardNth11beta1 = PDstandardNth11(beta1, i, j, k);
    CCTK_REAL const PDstandardNth22beta1 = PDstandardNth22(beta1, i, j, k);
    CCTK_REAL const PDstandardNth33beta1 = PDstandardNth33(beta1, i, j, k);
    CCTK_REAL const PDstandardNth12beta1 = PDstandardNth12(beta1, i, j, k);
    CCTK_REAL const PDstandardNth13beta1 = PDstandardNth13(beta1, i, j, k);
    CCTK_REAL const PDstandardNth23beta1 = PDstandardNth23(beta1, i, j, k);
    CCTK_REAL const PDstandardNth1beta2 = PDstandardNth1(beta2, i, j, k);
    CCTK_REAL const PDstandardNth2beta2 = PDstandardNth2(beta2, i, j, k);
    CCTK_REAL const PDstandardNth3beta2 = PDstandardNth3(beta2, i, j, k);
    CCTK_REAL const PDstandardNth11beta2 = PDstandardNth11(beta2, i, j, k);
    CCTK_REAL const PDstandardNth22beta2 = PDstandardNth22(beta2, i, j, k);
    CCTK_REAL const PDstandardNth33beta2 = PDstandardNth33(beta2, i, j, k);
    CCTK_REAL const PDstandardNth12beta2 = PDstandardNth12(beta2, i, j, k);
    CCTK_REAL const PDstandardNth13beta2 = PDstandardNth13(beta2, i, j, k);
    CCTK_REAL const PDstandardNth23beta2 = PDstandardNth23(beta2, i, j, k);
    CCTK_REAL const PDstandardNth1beta3 = PDstandardNth1(beta3, i, j, k);
    CCTK_REAL const PDstandardNth2beta3 = PDstandardNth2(beta3, i, j, k);
    CCTK_REAL const PDstandardNth3beta3 = PDstandardNth3(beta3, i, j, k);
    CCTK_REAL const PDstandardNth11beta3 = PDstandardNth11(beta3, i, j, k);
    CCTK_REAL const PDstandardNth22beta3 = PDstandardNth22(beta3, i, j, k);
    CCTK_REAL const PDstandardNth33beta3 = PDstandardNth33(beta3, i, j, k);
    CCTK_REAL const PDstandardNth12beta3 = PDstandardNth12(beta3, i, j, k);
    CCTK_REAL const PDstandardNth13beta3 = PDstandardNth13(beta3, i, j, k);
    CCTK_REAL const PDstandardNth23beta3 = PDstandardNth23(beta3, i, j, k);
    CCTK_REAL const PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    CCTK_REAL const PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    CCTK_REAL const PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    CCTK_REAL const PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    CCTK_REAL const PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    CCTK_REAL const PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    CCTK_REAL const PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    CCTK_REAL const PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    CCTK_REAL const PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    CCTK_REAL const PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    CCTK_REAL const PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    CCTK_REAL const PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    CCTK_REAL const PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    CCTK_REAL const PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    CCTK_REAL const PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    CCTK_REAL const PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    CCTK_REAL const PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    CCTK_REAL const PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    CCTK_REAL const PDstandardNth1phi = PDstandardNth1(phi, i, j, k);
    CCTK_REAL const PDstandardNth2phi = PDstandardNth2(phi, i, j, k);
    CCTK_REAL const PDstandardNth3phi = PDstandardNth3(phi, i, j, k);
    CCTK_REAL const PDstandardNth1trK = PDstandardNth1(trK, i, j, k);
    CCTK_REAL const PDstandardNth2trK = PDstandardNth2(trK, i, j, k);
    CCTK_REAL const PDstandardNth3trK = PDstandardNth3(trK, i, j, k);
    
    /* Calculate temporaries and grid functions */
    int dir1 = Sign(beta1L);
    
    int dir2 = Sign(beta2L);
    
    int dir3 = Sign(beta3L);
    
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu21 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu31 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu32 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL Gt111 = khalf*(gtu11*PDstandardNth1gt11 + 
      2*(gtu21*PDstandardNth1gt12 + gtu31*PDstandardNth1gt13) - 
      gtu21*PDstandardNth2gt11 - gtu31*PDstandardNth3gt11);
    
    CCTK_REAL Gt211 = khalf*(gtu21*PDstandardNth1gt11 + 
      2*(gtu22*PDstandardNth1gt12 + gtu32*PDstandardNth1gt13) - 
      gtu22*PDstandardNth2gt11 - gtu32*PDstandardNth3gt11);
    
    CCTK_REAL Gt311 = khalf*(gtu31*PDstandardNth1gt11 + 
      2*(gtu32*PDstandardNth1gt12 + gtu33*PDstandardNth1gt13) - 
      gtu32*PDstandardNth2gt11 - gtu33*PDstandardNth3gt11);
    
    CCTK_REAL Gt112 = khalf*(gtu21*PDstandardNth1gt22 + 
      gtu11*PDstandardNth2gt11 + gtu31*(PDstandardNth1gt23 + 
      PDstandardNth2gt13 - PDstandardNth3gt12));
    
    CCTK_REAL Gt212 = khalf*(gtu22*PDstandardNth1gt22 + 
      gtu21*PDstandardNth2gt11 + gtu32*(PDstandardNth1gt23 + 
      PDstandardNth2gt13 - PDstandardNth3gt12));
    
    CCTK_REAL Gt312 = khalf*(gtu32*PDstandardNth1gt22 + 
      gtu31*PDstandardNth2gt11 + gtu33*(PDstandardNth1gt23 + 
      PDstandardNth2gt13 - PDstandardNth3gt12));
    
    CCTK_REAL Gt113 = khalf*(gtu31*PDstandardNth1gt33 + 
      gtu11*PDstandardNth3gt11 + gtu21*(PDstandardNth1gt23 - 
      PDstandardNth2gt13 + PDstandardNth3gt12));
    
    CCTK_REAL Gt213 = khalf*(gtu32*PDstandardNth1gt33 + 
      gtu21*PDstandardNth3gt11 + gtu22*(PDstandardNth1gt23 - 
      PDstandardNth2gt13 + PDstandardNth3gt12));
    
    CCTK_REAL Gt313 = khalf*(gtu33*PDstandardNth1gt33 + 
      gtu31*PDstandardNth3gt11 + gtu32*(PDstandardNth1gt23 - 
      PDstandardNth2gt13 + PDstandardNth3gt12));
    
    CCTK_REAL Gt122 = khalf*(gtu11*(-PDstandardNth1gt22 + 
      2*PDstandardNth2gt12) + gtu21*PDstandardNth2gt22 + 
      gtu31*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    CCTK_REAL Gt222 = khalf*(gtu21*(-PDstandardNth1gt22 + 
      2*PDstandardNth2gt12) + gtu22*PDstandardNth2gt22 + 
      gtu32*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    CCTK_REAL Gt322 = khalf*(gtu31*(-PDstandardNth1gt22 + 
      2*PDstandardNth2gt12) + gtu32*PDstandardNth2gt22 + 
      gtu33*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    CCTK_REAL Gt123 = khalf*(gtu31*PDstandardNth2gt33 + 
      gtu11*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) + 
      gtu21*PDstandardNth3gt22);
    
    CCTK_REAL Gt223 = khalf*(gtu32*PDstandardNth2gt33 + 
      gtu21*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) + 
      gtu22*PDstandardNth3gt22);
    
    CCTK_REAL Gt323 = khalf*(gtu33*PDstandardNth2gt33 + 
      gtu31*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) + 
      gtu32*PDstandardNth3gt22);
    
    CCTK_REAL Gt133 = khalf*(-(gtu11*PDstandardNth1gt33) - 
      gtu21*PDstandardNth2gt33 + 2*gtu11*PDstandardNth3gt13 + 
      2*gtu21*PDstandardNth3gt23 + gtu31*PDstandardNth3gt33);
    
    CCTK_REAL Gt233 = khalf*(-(gtu21*PDstandardNth1gt33) - 
      gtu22*PDstandardNth2gt33 + 2*gtu21*PDstandardNth3gt13 + 
      2*gtu22*PDstandardNth3gt23 + gtu32*PDstandardNth3gt33);
    
    CCTK_REAL Gt333 = khalf*(-(gtu31*PDstandardNth1gt33) - 
      gtu32*PDstandardNth2gt33 + 2*gtu31*PDstandardNth3gt13 + 
      2*gtu32*PDstandardNth3gt23 + gtu33*PDstandardNth3gt33);
    
    CCTK_REAL Xtn1 = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + 
      Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    CCTK_REAL Xtn2 = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + 
      Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    CCTK_REAL Xtn3 = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + 
      Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    CCTK_REAL fac1 = IfThen(conformalMethod,-(khalf*INV(phiL)),1);
    
    CCTK_REAL cdphi1 = fac1*PDstandardNth1phi;
    
    CCTK_REAL cdphi2 = fac1*PDstandardNth2phi;
    
    CCTK_REAL cdphi3 = fac1*PDstandardNth3phi;
    
    CCTK_REAL Atm11 = At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    CCTK_REAL Atm21 = At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    CCTK_REAL Atm31 = At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    CCTK_REAL Atm12 = At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    CCTK_REAL Atm22 = At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    CCTK_REAL Atm32 = At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    CCTK_REAL Atm13 = At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    CCTK_REAL Atm23 = At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    CCTK_REAL Atm33 = At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    CCTK_REAL Atu11 = Atm11*gtu11 + Atm12*gtu21 + Atm13*gtu31;
    
    CCTK_REAL Atu21 = Atm11*gtu21 + Atm12*gtu22 + Atm13*gtu32;
    
    CCTK_REAL Atu31 = Atm11*gtu31 + Atm12*gtu32 + Atm13*gtu33;
    
    CCTK_REAL Atu22 = Atm21*gtu21 + Atm22*gtu22 + Atm23*gtu32;
    
    CCTK_REAL Atu32 = Atm21*gtu31 + Atm22*gtu32 + Atm23*gtu33;
    
    CCTK_REAL Atu33 = Atm31*gtu31 + Atm32*gtu32 + Atm33*gtu33;
    
    CCTK_REAL e4phi = IfThen(conformalMethod,pow(phiL,-2),exp(4*phiL));
    
    CCTK_REAL em4phi = INV(e4phi);
    
    CCTK_REAL rho = pow(alphaL,-2)*(eTttL - 2*(beta2L*eTtyL + 
      beta3L*eTtzL) + 2*(beta1L*(-eTtxL + beta2L*eTxyL + beta3L*eTxzL) + 
      beta2L*beta3L*eTyzL) + eTxxL*SQR(beta1L) + eTyyL*SQR(beta2L) + 
      eTzzL*SQR(beta3L));
    
    CCTK_REAL trS = em4phi*(eTxxL*gtu11 + eTyyL*gtu22 + 2*(eTxyL*gtu21 + 
      eTxzL*gtu31 + eTyzL*gtu32) + eTzzL*gtu33);
    
    CCTK_REAL S1 = (-eTtxL + beta1L*eTxxL + beta2L*eTxyL + 
      beta3L*eTxzL)*INV(alphaL);
    
    CCTK_REAL S2 = (-eTtyL + beta1L*eTxyL + beta2L*eTyyL + 
      beta3L*eTyzL)*INV(alphaL);
    
    CCTK_REAL S3 = (-eTtzL + beta1L*eTxzL + beta2L*eTyzL + 
      beta3L*eTzzL)*INV(alphaL);
    
    trKrhsL = PDupwindNth1(trK, i, j, k)*beta1L + PDupwindNth2(trK, i, 
      j, k)*beta2L + PDupwindNth3(trK, i, j, k)*beta3L - 
      em4phi*(gtu11*PDstandardNth11alpha + gtu22*PDstandardNth22alpha + 
      gtu33*(PDstandardNth33alpha + 2*cdphi3*PDstandardNth3alpha) + 
      2*(gtu21*PDstandardNth12alpha + gtu31*(PDstandardNth13alpha + 
      cdphi1*PDstandardNth3alpha) + gtu32*(PDstandardNth23alpha + 
      cdphi2*PDstandardNth3alpha)) + PDstandardNth1alpha*(2*(cdphi1*gtu11 + 
      cdphi2*gtu21 + cdphi3*gtu31) - Xtn1) + 
      PDstandardNth2alpha*(2*(cdphi1*gtu21 + cdphi2*gtu22 + cdphi3*gtu32) - 
      Xtn2) - PDstandardNth3alpha*Xtn3) + alphaL*(2*(Atm12*Atm21 + 
      Atm13*Atm31 + Atm23*Atm32) + 
      12.56637061435917295385057353311801153679*(rho + trS) + SQR(Atm11) + 
      SQR(Atm22) + SQR(Atm33) + kthird*SQR(trKL));
    
    CCTK_REAL phirhsL = PDupwindNth1(phi, i, j, k)*beta1L + 
      PDupwindNth2(phi, i, j, k)*beta2L + PDupwindNth3(phi, i, j, 
      k)*beta3L + (PDstandardNth1beta1 + PDstandardNth2beta2 + 
      PDstandardNth3beta3)*IfThen(conformalMethod,-(kthird*phiL),0.16666666666666666) 
      + alphaL*trKL*IfThen(conformalMethod,kthird*phiL,-0.16666666666666666);
    
    CCTK_REAL gt11rhsL = -2*alphaL*At11L + PDupwindNth1(gt11, i, j, 
      k)*beta1L + PDupwindNth2(gt11, i, j, k)*beta2L + PDupwindNth3(gt11, 
      i, j, k)*beta3L + 2*(gt12L*PDstandardNth1beta2 + 
      gt13L*PDstandardNth1beta3) + gt11L*(kfourthird*PDstandardNth1beta1 - 
      ktwothird*PDstandardNth2beta2 - ktwothird*PDstandardNth3beta3);
    
    CCTK_REAL gt12rhsL = -2*alphaL*At12L + PDupwindNth1(gt12, i, j, 
      k)*beta1L + PDupwindNth2(gt12, i, j, k)*beta2L + PDupwindNth3(gt12, 
      i, j, k)*beta3L + gt22L*PDstandardNth1beta2 + 
      gt23L*PDstandardNth1beta3 + gt11L*PDstandardNth2beta1 + 
      gt13L*PDstandardNth2beta3 + gt12L*(kthird*(PDstandardNth1beta1 + 
      PDstandardNth2beta2) - ktwothird*PDstandardNth3beta3);
    
    CCTK_REAL gt13rhsL = -2*alphaL*At13L + PDupwindNth1(gt13, i, j, 
      k)*beta1L + PDupwindNth2(gt13, i, j, k)*beta2L + PDupwindNth3(gt13, 
      i, j, k)*beta3L + gt23L*PDstandardNth1beta2 + 
      gt33L*PDstandardNth1beta3 + gt11L*PDstandardNth3beta1 + 
      gt12L*PDstandardNth3beta2 + gt13L*(-(ktwothird*PDstandardNth2beta2) + 
      kthird*(PDstandardNth1beta1 + PDstandardNth3beta3));
    
    CCTK_REAL gt22rhsL = -2*alphaL*At22L + PDupwindNth1(gt22, i, j, 
      k)*beta1L + PDupwindNth2(gt22, i, j, k)*beta2L + PDupwindNth3(gt22, 
      i, j, k)*beta3L + 2*(gt12L*PDstandardNth2beta1 + 
      gt23L*PDstandardNth2beta3) + gt22L*(-(ktwothird*PDstandardNth1beta1) + 
      kfourthird*PDstandardNth2beta2 - ktwothird*PDstandardNth3beta3);
    
    CCTK_REAL gt23rhsL = -2*alphaL*At23L + PDupwindNth1(gt23, i, j, 
      k)*beta1L + PDupwindNth2(gt23, i, j, k)*beta2L + PDupwindNth3(gt23, 
      i, j, k)*beta3L + gt13L*PDstandardNth2beta1 + 
      gt33L*PDstandardNth2beta3 + gt12L*PDstandardNth3beta1 + 
      gt22L*PDstandardNth3beta2 + gt23L*(-(ktwothird*PDstandardNth1beta1) + 
      kthird*(PDstandardNth2beta2 + PDstandardNth3beta3));
    
    CCTK_REAL gt33rhsL = -2*alphaL*At33L + PDupwindNth1(gt33, i, j, 
      k)*beta1L + PDupwindNth2(gt33, i, j, k)*beta2L + PDupwindNth3(gt33, 
      i, j, k)*beta3L - gt33L*ktwothird*PDstandardNth1beta1 - 
      gt33L*ktwothird*PDstandardNth2beta2 + 2*gt13L*PDstandardNth3beta1 + 
      2*gt23L*PDstandardNth3beta2 + gt33L*kfourthird*PDstandardNth3beta3;
    
    Xt1rhsL = kthird*(7*(gtu21*PDstandardNth12beta1 + 
      gtu31*PDstandardNth13beta1) + gtu11*(4*PDstandardNth11beta1 + 
      PDstandardNth12beta2 + PDstandardNth13beta3) + 
      gtu21*(PDstandardNth22beta2 + PDstandardNth23beta3) + 
      3*(PDupwindNth1(Xt1, i, j, k)*beta1L + PDupwindNth2(Xt1, i, j, 
      k)*beta2L + PDupwindNth3(Xt1, i, j, k)*beta3L + 
      gtu22*PDstandardNth22beta1 + gtu33*PDstandardNth33beta1) + 
      gtu31*(PDstandardNth23beta2 + PDstandardNth33beta3) - 
      6*(Atu11*PDstandardNth1alpha + Atu21*PDstandardNth2alpha + 
      Atu31*PDstandardNth3alpha) + 6*(gtu32*PDstandardNth23beta1 + 
      alphaL*(6*(Atu11*cdphi1 + Atu21*cdphi2 + Atu31*cdphi3) + Atu11*Gt111 + 
      Atu22*Gt122 + 2*(Atu21*Gt112 + Atu31*Gt113 + Atu32*Gt123) + Atu33*Gt133 
      - ktwothird*(gtu11*PDstandardNth1trK + gtu21*PDstandardNth2trK + 
      gtu31*PDstandardNth3trK))) - 
      150.7964473723100754462068823974161384415*alphaL*(gtu11*S1 + gtu21*S2 + 
      gtu31*S3) + (-3*PDstandardNth1beta1 + 2*(PDstandardNth1beta1 + 
      PDstandardNth2beta2 + PDstandardNth3beta3))*Xtn1 - 
      3*(PDstandardNth2beta1*Xtn2 + PDstandardNth3beta1*Xtn3));
    
    Xt2rhsL = kthird*(gtu21*(PDstandardNth11beta1 + 7*PDstandardNth12beta2 
      + PDstandardNth13beta3) + gtu22*(PDstandardNth12beta1 + 
      4*PDstandardNth22beta2 + PDstandardNth23beta3) + 3*(PDupwindNth1(Xt2, 
      i, j, k)*beta1L + PDupwindNth2(Xt2, i, j, k)*beta2L + 
      PDupwindNth3(Xt2, i, j, k)*beta3L + gtu11*PDstandardNth11beta2 + 
      gtu33*PDstandardNth33beta2) + gtu32*(PDstandardNth13beta1 + 
      7*PDstandardNth23beta2 + PDstandardNth33beta3) - 
      6*(Atu21*PDstandardNth1alpha + Atu22*PDstandardNth2alpha + 
      Atu32*PDstandardNth3alpha) + 6*(gtu31*PDstandardNth13beta2 + 
      alphaL*(6*(Atu21*cdphi1 + Atu22*cdphi2 + Atu32*cdphi3) + Atu11*Gt211 + 
      Atu22*Gt222 + 2*(Atu21*Gt212 + Atu31*Gt213 + Atu32*Gt223) + Atu33*Gt233 
      - ktwothird*(gtu21*PDstandardNth1trK + gtu22*PDstandardNth2trK + 
      gtu32*PDstandardNth3trK))) - 
      150.7964473723100754462068823974161384415*alphaL*(gtu21*S1 + gtu22*S2 + 
      gtu32*S3) + 2*(PDstandardNth1beta1 + PDstandardNth2beta2 + 
      PDstandardNth3beta3)*Xtn2 - 3*(PDstandardNth1beta2*Xtn1 + 
      PDstandardNth2beta2*Xtn2 + PDstandardNth3beta2*Xtn3));
    
    Xt3rhsL = kthird*(gtu31*(PDstandardNth11beta1 + PDstandardNth12beta2 + 
      7*PDstandardNth13beta3) + 3*(PDupwindNth1(Xt3, i, j, k)*beta1L + 
      PDupwindNth2(Xt3, i, j, k)*beta2L + PDupwindNth3(Xt3, i, j, 
      k)*beta3L + gtu11*PDstandardNth11beta3 + gtu22*PDstandardNth22beta3) + 
      gtu32*(PDstandardNth12beta1 + PDstandardNth22beta2 + 
      7*PDstandardNth23beta3) + gtu33*(PDstandardNth13beta1 + 
      PDstandardNth23beta2 + 4*PDstandardNth33beta3) - 
      6*(Atu31*PDstandardNth1alpha + Atu32*PDstandardNth2alpha + 
      Atu33*PDstandardNth3alpha) + 6*(gtu21*PDstandardNth12beta3 + 
      alphaL*(6*(Atu31*cdphi1 + Atu32*cdphi2 + Atu33*cdphi3) + Atu11*Gt311 + 
      Atu22*Gt322 + 2*(Atu21*Gt312 + Atu31*Gt313 + Atu32*Gt323) + Atu33*Gt333 
      - ktwothird*(gtu31*PDstandardNth1trK + gtu32*PDstandardNth2trK + 
      gtu33*PDstandardNth3trK))) - 
      150.7964473723100754462068823974161384415*alphaL*(gtu31*S1 + gtu32*S2 + 
      gtu33*S3) + 2*(PDstandardNth1beta1 + PDstandardNth2beta2 + 
      PDstandardNth3beta3)*Xtn3 - 3*(PDstandardNth1beta3*Xtn1 + 
      PDstandardNth2beta3*Xtn2 + PDstandardNth3beta3*Xtn3));
    
    CCTK_REAL eta = BetaDriver*IfThen(rL > 
      SpatialBetaDriverRadius,SpatialBetaDriverRadius*INV(rL),1);
    
    CCTK_REAL alpharhsL = (PDupwindNth1(alpha, i, j, k)*beta1L + 
      PDupwindNth2(alpha, i, j, k)*beta2L + PDupwindNth3(alpha, i, j, 
      k)*beta3L)*LapseAdvectionCoeff + harmonicF*(AL*(-1 + 
      LapseAdvectionCoeff) - LapseAdvectionCoeff*trKL)*pow(alphaL,harmonicN);
    
    CCTK_REAL ArhsL = (-1 + LapseAdvectionCoeff)*(AL*AlphaDriver - 
      trKrhsL);
    
    CCTK_REAL beta1rhsL = (PDupwindNth1(beta1, i, j, k)*beta1L + 
      PDupwindNth2(beta1, i, j, k)*beta2L + PDupwindNth3(beta1, i, j, 
      k)*beta3L)*ShiftAdvectionCoeff + B1L*ShiftGammaCoeff;
    
    CCTK_REAL beta2rhsL = (PDupwindNth1(beta2, i, j, k)*beta1L + 
      PDupwindNth2(beta2, i, j, k)*beta2L + PDupwindNth3(beta2, i, j, 
      k)*beta3L)*ShiftAdvectionCoeff + B2L*ShiftGammaCoeff;
    
    CCTK_REAL beta3rhsL = (PDupwindNth1(beta3, i, j, k)*beta1L + 
      PDupwindNth2(beta3, i, j, k)*beta2L + PDupwindNth3(beta3, i, j, 
      k)*beta3L)*ShiftAdvectionCoeff + B3L*ShiftGammaCoeff;
    
    CCTK_REAL B1rhsL = -(B1L*eta) + ((PDupwindNth1(B1, i, j, k) - 
      PDupwindNth1(Xt1, i, j, k))*beta1L + (PDupwindNth2(B1, i, j, k) - 
      PDupwindNth2(Xt1, i, j, k))*beta2L + (PDupwindNth3(B1, i, j, k) - 
      PDupwindNth3(Xt1, i, j, k))*beta3L)*ShiftAdvectionCoeff + Xt1rhsL;
    
    CCTK_REAL B2rhsL = -(B2L*eta) + ((PDupwindNth1(B2, i, j, k) - 
      PDupwindNth1(Xt2, i, j, k))*beta1L + (PDupwindNth2(B2, i, j, k) - 
      PDupwindNth2(Xt2, i, j, k))*beta2L + (PDupwindNth3(B2, i, j, k) - 
      PDupwindNth3(Xt2, i, j, k))*beta3L)*ShiftAdvectionCoeff + Xt2rhsL;
    
    CCTK_REAL B3rhsL = -(B3L*eta) + ((PDupwindNth1(B3, i, j, k) - 
      PDupwindNth1(Xt3, i, j, k))*beta1L + (PDupwindNth2(B3, i, j, k) - 
      PDupwindNth2(Xt3, i, j, k))*beta2L + (PDupwindNth3(B3, i, j, k) - 
      PDupwindNth3(Xt3, i, j, k))*beta3L)*ShiftAdvectionCoeff + Xt3rhsL;
    
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
    B1rhs[index] = B1rhsL;
    B2rhs[index] = B2rhsL;
    B3rhs[index] = B3rhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    gt11rhs[index] = gt11rhsL;
    gt12rhs[index] = gt12rhsL;
    gt13rhs[index] = gt13rhsL;
    gt22rhs[index] = gt22rhsL;
    gt23rhs[index] = gt23rhsL;
    gt33rhs[index] = gt33rhsL;
    phirhs[index] = phirhsL;
    trKrhs[index] = trKrhsL;
    Xt1rhs[index] = Xt1rhsL;
    Xt2rhs[index] = Xt2rhsL;
    Xt3rhs[index] = Xt3rhsL;
  }
  LC_ENDLOOP3 (ML_BSSN_O2_RHS1);
}

void ML_BSSN_O2_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_O2_RHS1_Body);
}
