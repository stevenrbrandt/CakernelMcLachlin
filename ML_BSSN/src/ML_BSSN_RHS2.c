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

void ML_BSSN_RHS2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  /* Declare predefined quantities */
  // CCTK_REAL p1o12dx = INITVALUE;
  // CCTK_REAL p1o12dy = INITVALUE;
  // CCTK_REAL p1o12dz = INITVALUE;
  // CCTK_REAL p1o144dxdy = INITVALUE;
  // CCTK_REAL p1o144dxdz = INITVALUE;
  // CCTK_REAL p1o144dydz = INITVALUE;
  // CCTK_REAL p1odx = INITVALUE;
  // CCTK_REAL p1ody = INITVALUE;
  // CCTK_REAL p1odz = INITVALUE;
  // CCTK_REAL pm1o12dx2 = INITVALUE;
  // CCTK_REAL pm1o12dy2 = INITVALUE;
  // CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_RHS2_Body");
  }
  
  if (cctk_iteration % ML_BSSN_RHS2_calc_every != ML_BSSN_RHS2_calc_offset)
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
  LC_LOOP3 (ML_BSSN_RHS2,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    // int subblock_index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int const subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    // CCTK_REAL Atm11 = INITVALUE, Atm12 = INITVALUE, Atm13 = INITVALUE, Atm21 = INITVALUE, Atm22 = INITVALUE, Atm23 = INITVALUE;
    // CCTK_REAL Atm31 = INITVALUE, Atm32 = INITVALUE, Atm33 = INITVALUE;
    // CCTK_REAL Ats11 = INITVALUE, Ats12 = INITVALUE, Ats13 = INITVALUE, Ats22 = INITVALUE, Ats23 = INITVALUE, Ats33 = INITVALUE;
    // CCTK_REAL cdphi1 = INITVALUE, cdphi2 = INITVALUE, cdphi211 = INITVALUE, cdphi212 = INITVALUE, cdphi213 = INITVALUE, cdphi222 = INITVALUE;
    // CCTK_REAL cdphi223 = INITVALUE, cdphi233 = INITVALUE, cdphi3 = INITVALUE;
    // CCTK_REAL detgt = INITVALUE;
    // CCTK_REAL dir1 = INITVALUE, dir2 = INITVALUE, dir3 = INITVALUE;
    // CCTK_REAL e4phi = INITVALUE;
    // CCTK_REAL em4phi = INITVALUE;
    // CCTK_REAL fac1 = INITVALUE, fac2 = INITVALUE;
    // CCTK_REAL g11 = INITVALUE;
    // CCTK_REAL G111 = INITVALUE, G112 = INITVALUE, G113 = INITVALUE;
    // CCTK_REAL g12 = INITVALUE;
    // CCTK_REAL G122 = INITVALUE, G123 = INITVALUE;
    // CCTK_REAL g13 = INITVALUE;
    // CCTK_REAL G133 = INITVALUE, G211 = INITVALUE, G212 = INITVALUE, G213 = INITVALUE;
    // CCTK_REAL g22 = INITVALUE;
    // CCTK_REAL G222 = INITVALUE, G223 = INITVALUE;
    // CCTK_REAL g23 = INITVALUE;
    // CCTK_REAL G233 = INITVALUE, G311 = INITVALUE, G312 = INITVALUE, G313 = INITVALUE, G322 = INITVALUE, G323 = INITVALUE;
    // CCTK_REAL g33 = INITVALUE;
    // CCTK_REAL G333 = INITVALUE;
    // CCTK_REAL Gt111 = INITVALUE, Gt112 = INITVALUE, Gt113 = INITVALUE, Gt122 = INITVALUE, Gt123 = INITVALUE, Gt133 = INITVALUE;
    // CCTK_REAL Gt211 = INITVALUE, Gt212 = INITVALUE, Gt213 = INITVALUE, Gt222 = INITVALUE, Gt223 = INITVALUE, Gt233 = INITVALUE;
    // CCTK_REAL Gt311 = INITVALUE, Gt312 = INITVALUE, Gt313 = INITVALUE, Gt322 = INITVALUE, Gt323 = INITVALUE, Gt333 = INITVALUE;
    // CCTK_REAL Gtl111 = INITVALUE, Gtl112 = INITVALUE, Gtl113 = INITVALUE, Gtl122 = INITVALUE, Gtl123 = INITVALUE, Gtl133 = INITVALUE;
    // CCTK_REAL Gtl211 = INITVALUE, Gtl212 = INITVALUE, Gtl213 = INITVALUE, Gtl222 = INITVALUE, Gtl223 = INITVALUE, Gtl233 = INITVALUE;
    // CCTK_REAL Gtl311 = INITVALUE, Gtl312 = INITVALUE, Gtl313 = INITVALUE, Gtl322 = INITVALUE, Gtl323 = INITVALUE, Gtl333 = INITVALUE;
    // CCTK_REAL Gtlu111 = INITVALUE, Gtlu112 = INITVALUE, Gtlu113 = INITVALUE, Gtlu121 = INITVALUE, Gtlu122 = INITVALUE, Gtlu123 = INITVALUE;
    // CCTK_REAL Gtlu131 = INITVALUE, Gtlu132 = INITVALUE, Gtlu133 = INITVALUE, Gtlu211 = INITVALUE, Gtlu212 = INITVALUE, Gtlu213 = INITVALUE;
    // CCTK_REAL Gtlu221 = INITVALUE, Gtlu222 = INITVALUE, Gtlu223 = INITVALUE, Gtlu231 = INITVALUE, Gtlu232 = INITVALUE, Gtlu233 = INITVALUE;
    // CCTK_REAL Gtlu311 = INITVALUE, Gtlu312 = INITVALUE, Gtlu313 = INITVALUE, Gtlu321 = INITVALUE, Gtlu322 = INITVALUE, Gtlu323 = INITVALUE;
    // CCTK_REAL Gtlu331 = INITVALUE, Gtlu332 = INITVALUE, Gtlu333 = INITVALUE;
    // CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    // CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    // CCTK_REAL R11 = INITVALUE, R12 = INITVALUE, R13 = INITVALUE, R22 = INITVALUE, R23 = INITVALUE, R33 = INITVALUE;
    // CCTK_REAL rho = INITVALUE;
    // CCTK_REAL Rphi11 = INITVALUE, Rphi12 = INITVALUE, Rphi13 = INITVALUE, Rphi22 = INITVALUE, Rphi23 = INITVALUE, Rphi33 = INITVALUE;
    // CCTK_REAL Rt11 = INITVALUE, Rt12 = INITVALUE, Rt13 = INITVALUE, Rt22 = INITVALUE, Rt23 = INITVALUE, Rt33 = INITVALUE;
    // CCTK_REAL trAts = INITVALUE;
    // CCTK_REAL trS = INITVALUE;
    // CCTK_REAL Xtn1 = INITVALUE, Xtn2 = INITVALUE, Xtn3 = INITVALUE;
    
    /* Declare local copies of grid functions */
    // CCTK_REAL AL = INITVALUE;
    // CCTK_REAL alphaL = INITVALUE, alpharhsL = INITVALUE;
    // CCTK_REAL ArhsL = INITVALUE;
    // CCTK_REAL At11L = INITVALUE, At11rhsL = INITVALUE, At12L = INITVALUE, At12rhsL = INITVALUE, At13L = INITVALUE, At13rhsL = INITVALUE;
    // CCTK_REAL At22L = INITVALUE, At22rhsL = INITVALUE, At23L = INITVALUE, At23rhsL = INITVALUE, At33L = INITVALUE, At33rhsL = INITVALUE;
    // CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    // CCTK_REAL eTttL = INITVALUE;
    // CCTK_REAL eTtxL = INITVALUE;
    // CCTK_REAL eTtyL = INITVALUE;
    // CCTK_REAL eTtzL = INITVALUE;
    // CCTK_REAL eTxxL = INITVALUE;
    // CCTK_REAL eTxyL = INITVALUE;
    // CCTK_REAL eTxzL = INITVALUE;
    // CCTK_REAL eTyyL = INITVALUE;
    // CCTK_REAL eTyzL = INITVALUE;
    // CCTK_REAL eTzzL = INITVALUE;
    // CCTK_REAL gt11L = INITVALUE, gt12L = INITVALUE, gt13L = INITVALUE, gt22L = INITVALUE, gt23L = INITVALUE, gt33L = INITVALUE;
    // CCTK_REAL phiL = INITVALUE;
    // CCTK_REAL trKL = INITVALUE, trKrhsL = INITVALUE;
    // CCTK_REAL Xt1L = INITVALUE, Xt2L = INITVALUE, Xt3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
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
    // CCTK_REAL PDstandardNth1beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth1beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta3 = INITVALUE;
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
    // CCTK_REAL PDstandardNth1phi = INITVALUE;
    // CCTK_REAL PDstandardNth2phi = INITVALUE;
    // CCTK_REAL PDstandardNth3phi = INITVALUE;
    // CCTK_REAL PDstandardNth11phi = INITVALUE;
    // CCTK_REAL PDstandardNth22phi = INITVALUE;
    // CCTK_REAL PDstandardNth33phi = INITVALUE;
    // CCTK_REAL PDstandardNth12phi = INITVALUE;
    // CCTK_REAL PDstandardNth13phi = INITVALUE;
    // CCTK_REAL PDstandardNth23phi = INITVALUE;
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
    CCTK_REAL const AL = A[index];
    CCTK_REAL const alphaL = alpha[index];
    CCTK_REAL const At11L = At11[index];
    CCTK_REAL const At12L = At12[index];
    CCTK_REAL const At13L = At13[index];
    CCTK_REAL const At22L = At22[index];
    CCTK_REAL const At23L = At23[index];
    CCTK_REAL const At33L = At33[index];
    CCTK_REAL const beta1L = beta1[index];
    CCTK_REAL const beta2L = beta2[index];
    CCTK_REAL const beta3L = beta3[index];
    CCTK_REAL const eTttL = (*stress_energy_state) ? (eTtt[index]) : 0.0;
    CCTK_REAL const eTtxL = (*stress_energy_state) ? (eTtx[index]) : 0.0;
    CCTK_REAL const eTtyL = (*stress_energy_state) ? (eTty[index]) : 0.0;
    CCTK_REAL const eTtzL = (*stress_energy_state) ? (eTtz[index]) : 0.0;
    CCTK_REAL const eTxxL = (*stress_energy_state) ? (eTxx[index]) : 0.0;
    CCTK_REAL const eTxyL = (*stress_energy_state) ? (eTxy[index]) : 0.0;
    CCTK_REAL const eTxzL = (*stress_energy_state) ? (eTxz[index]) : 0.0;
    CCTK_REAL const eTyyL = (*stress_energy_state) ? (eTyy[index]) : 0.0;
    CCTK_REAL const eTyzL = (*stress_energy_state) ? (eTyz[index]) : 0.0;
    CCTK_REAL const eTzzL = (*stress_energy_state) ? (eTzz[index]) : 0.0;
    CCTK_REAL const gt11L = gt11[index];
    CCTK_REAL const gt12L = gt12[index];
    CCTK_REAL const gt13L = gt13[index];
    CCTK_REAL const gt22L = gt22[index];
    CCTK_REAL const gt23L = gt23[index];
    CCTK_REAL const gt33L = gt33[index];
    CCTK_REAL const phiL = phi[index];
    CCTK_REAL const trKL = trK[index];
    CCTK_REAL const Xt1L = Xt1[index];
    CCTK_REAL const Xt2L = Xt2[index];
    CCTK_REAL const Xt3L = Xt3[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
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
    CCTK_REAL const PDstandardNth1beta2 = PDstandardNth1(beta2, i, j, k);
    CCTK_REAL const PDstandardNth2beta2 = PDstandardNth2(beta2, i, j, k);
    CCTK_REAL const PDstandardNth3beta2 = PDstandardNth3(beta2, i, j, k);
    CCTK_REAL const PDstandardNth1beta3 = PDstandardNth1(beta3, i, j, k);
    CCTK_REAL const PDstandardNth2beta3 = PDstandardNth2(beta3, i, j, k);
    CCTK_REAL const PDstandardNth3beta3 = PDstandardNth3(beta3, i, j, k);
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
    CCTK_REAL const PDstandardNth1phi = PDstandardNth1(phi, i, j, k);
    CCTK_REAL const PDstandardNth2phi = PDstandardNth2(phi, i, j, k);
    CCTK_REAL const PDstandardNth3phi = PDstandardNth3(phi, i, j, k);
    CCTK_REAL const PDstandardNth11phi = PDstandardNth11(phi, i, j, k);
    CCTK_REAL const PDstandardNth22phi = PDstandardNth22(phi, i, j, k);
    CCTK_REAL const PDstandardNth33phi = PDstandardNth33(phi, i, j, k);
    CCTK_REAL const PDstandardNth12phi = PDstandardNth12(phi, i, j, k);
    CCTK_REAL const PDstandardNth13phi = PDstandardNth13(phi, i, j, k);
    CCTK_REAL const PDstandardNth23phi = PDstandardNth23(phi, i, j, k);
    CCTK_REAL const PDstandardNth1Xt1 = PDstandardNth1(Xt1, i, j, k);
    CCTK_REAL const PDstandardNth2Xt1 = PDstandardNth2(Xt1, i, j, k);
    CCTK_REAL const PDstandardNth3Xt1 = PDstandardNth3(Xt1, i, j, k);
    CCTK_REAL const PDstandardNth1Xt2 = PDstandardNth1(Xt2, i, j, k);
    CCTK_REAL const PDstandardNth2Xt2 = PDstandardNth2(Xt2, i, j, k);
    CCTK_REAL const PDstandardNth3Xt2 = PDstandardNth3(Xt2, i, j, k);
    CCTK_REAL const PDstandardNth1Xt3 = PDstandardNth1(Xt3, i, j, k);
    CCTK_REAL const PDstandardNth2Xt3 = PDstandardNth2(Xt3, i, j, k);
    CCTK_REAL const PDstandardNth3Xt3 = PDstandardNth3(Xt3, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    int const dir1  =  Sign(beta1L);
    
    int const dir2  =  Sign(beta2L);
    
    int const dir3  =  Sign(beta3L);
    
    CCTK_REAL const detgt  =  1;
    
    CCTK_REAL const gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL const gtu21  =  (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL const gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL const gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL const gtu32  =  (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL const gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL const Gtl111  =  khalf*PDstandardNth1gt11;
    
    CCTK_REAL const Gtl112  =  khalf*PDstandardNth2gt11;
    
    CCTK_REAL const Gtl113  =  khalf*PDstandardNth3gt11;
    
    CCTK_REAL const Gtl122  =  -(khalf*PDstandardNth1gt22) + PDstandardNth2gt12;
    
    CCTK_REAL const Gtl123  =  khalf*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12);
    
    CCTK_REAL const Gtl133  =  -(khalf*PDstandardNth1gt33) + PDstandardNth3gt13;
    
    CCTK_REAL const Gtl211  =  PDstandardNth1gt12 - khalf*PDstandardNth2gt11;
    
    CCTK_REAL const Gtl212  =  khalf*PDstandardNth1gt22;
    
    CCTK_REAL const Gtl213  =  khalf*(PDstandardNth1gt23 - PDstandardNth2gt13 + PDstandardNth3gt12);
    
    CCTK_REAL const Gtl222  =  khalf*PDstandardNth2gt22;
    
    CCTK_REAL const Gtl223  =  khalf*PDstandardNth3gt22;
    
    CCTK_REAL const Gtl233  =  -(khalf*PDstandardNth2gt33) + PDstandardNth3gt23;
    
    CCTK_REAL const Gtl311  =  PDstandardNth1gt13 - khalf*PDstandardNth3gt11;
    
    CCTK_REAL const Gtl312  =  khalf*(PDstandardNth1gt23 + PDstandardNth2gt13 - PDstandardNth3gt12);
    
    CCTK_REAL const Gtl313  =  khalf*PDstandardNth1gt33;
    
    CCTK_REAL const Gtl322  =  PDstandardNth2gt23 - khalf*PDstandardNth3gt22;
    
    CCTK_REAL const Gtl323  =  khalf*PDstandardNth2gt33;
    
    CCTK_REAL const Gtl333  =  khalf*PDstandardNth3gt33;
    
    CCTK_REAL const Gtlu111  =  Gtl111*gtu11 + Gtl112*gtu21 + Gtl113*gtu31;
    
    CCTK_REAL const Gtlu112  =  Gtl111*gtu21 + Gtl112*gtu22 + Gtl113*gtu32;
    
    CCTK_REAL const Gtlu113  =  Gtl111*gtu31 + Gtl112*gtu32 + Gtl113*gtu33;
    
    CCTK_REAL const Gtlu121  =  Gtl112*gtu11 + Gtl122*gtu21 + Gtl123*gtu31;
    
    CCTK_REAL const Gtlu122  =  Gtl112*gtu21 + Gtl122*gtu22 + Gtl123*gtu32;
    
    CCTK_REAL const Gtlu123  =  Gtl112*gtu31 + Gtl122*gtu32 + Gtl123*gtu33;
    
    CCTK_REAL const Gtlu131  =  Gtl113*gtu11 + Gtl123*gtu21 + Gtl133*gtu31;
    
    CCTK_REAL const Gtlu132  =  Gtl113*gtu21 + Gtl123*gtu22 + Gtl133*gtu32;
    
    CCTK_REAL const Gtlu133  =  Gtl113*gtu31 + Gtl123*gtu32 + Gtl133*gtu33;
    
    CCTK_REAL const Gtlu211  =  Gtl211*gtu11 + Gtl212*gtu21 + Gtl213*gtu31;
    
    CCTK_REAL const Gtlu212  =  Gtl211*gtu21 + Gtl212*gtu22 + Gtl213*gtu32;
    
    CCTK_REAL const Gtlu213  =  Gtl211*gtu31 + Gtl212*gtu32 + Gtl213*gtu33;
    
    CCTK_REAL const Gtlu221  =  Gtl212*gtu11 + Gtl222*gtu21 + Gtl223*gtu31;
    
    CCTK_REAL const Gtlu222  =  Gtl212*gtu21 + Gtl222*gtu22 + Gtl223*gtu32;
    
    CCTK_REAL const Gtlu223  =  Gtl212*gtu31 + Gtl222*gtu32 + Gtl223*gtu33;
    
    CCTK_REAL const Gtlu231  =  Gtl213*gtu11 + Gtl223*gtu21 + Gtl233*gtu31;
    
    CCTK_REAL const Gtlu232  =  Gtl213*gtu21 + Gtl223*gtu22 + Gtl233*gtu32;
    
    CCTK_REAL const Gtlu233  =  Gtl213*gtu31 + Gtl223*gtu32 + Gtl233*gtu33;
    
    CCTK_REAL const Gtlu311  =  Gtl311*gtu11 + Gtl312*gtu21 + Gtl313*gtu31;
    
    CCTK_REAL const Gtlu312  =  Gtl311*gtu21 + Gtl312*gtu22 + Gtl313*gtu32;
    
    CCTK_REAL const Gtlu313  =  Gtl311*gtu31 + Gtl312*gtu32 + Gtl313*gtu33;
    
    CCTK_REAL const Gtlu321  =  Gtl312*gtu11 + Gtl322*gtu21 + Gtl323*gtu31;
    
    CCTK_REAL const Gtlu322  =  Gtl312*gtu21 + Gtl322*gtu22 + Gtl323*gtu32;
    
    CCTK_REAL const Gtlu323  =  Gtl312*gtu31 + Gtl322*gtu32 + Gtl323*gtu33;
    
    CCTK_REAL const Gtlu331  =  Gtl313*gtu11 + Gtl323*gtu21 + Gtl333*gtu31;
    
    CCTK_REAL const Gtlu332  =  Gtl313*gtu21 + Gtl323*gtu22 + Gtl333*gtu32;
    
    CCTK_REAL const Gtlu333  =  Gtl313*gtu31 + Gtl323*gtu32 + Gtl333*gtu33;
    
    CCTK_REAL const Gt111  =  Gtl111*gtu11 + Gtl211*gtu21 + Gtl311*gtu31;
    
    CCTK_REAL const Gt211  =  Gtl111*gtu21 + Gtl211*gtu22 + Gtl311*gtu32;
    
    CCTK_REAL const Gt311  =  Gtl111*gtu31 + Gtl211*gtu32 + Gtl311*gtu33;
    
    CCTK_REAL const Gt112  =  Gtl112*gtu11 + Gtl212*gtu21 + Gtl312*gtu31;
    
    CCTK_REAL const Gt212  =  Gtl112*gtu21 + Gtl212*gtu22 + Gtl312*gtu32;
    
    CCTK_REAL const Gt312  =  Gtl112*gtu31 + Gtl212*gtu32 + Gtl312*gtu33;
    
    CCTK_REAL const Gt113  =  Gtl113*gtu11 + Gtl213*gtu21 + Gtl313*gtu31;
    
    CCTK_REAL const Gt213  =  Gtl113*gtu21 + Gtl213*gtu22 + Gtl313*gtu32;
    
    CCTK_REAL const Gt313  =  Gtl113*gtu31 + Gtl213*gtu32 + Gtl313*gtu33;
    
    CCTK_REAL const Gt122  =  Gtl122*gtu11 + Gtl222*gtu21 + Gtl322*gtu31;
    
    CCTK_REAL const Gt222  =  Gtl122*gtu21 + Gtl222*gtu22 + Gtl322*gtu32;
    
    CCTK_REAL const Gt322  =  Gtl122*gtu31 + Gtl222*gtu32 + Gtl322*gtu33;
    
    CCTK_REAL const Gt123  =  Gtl123*gtu11 + Gtl223*gtu21 + Gtl323*gtu31;
    
    CCTK_REAL const Gt223  =  Gtl123*gtu21 + Gtl223*gtu22 + Gtl323*gtu32;
    
    CCTK_REAL const Gt323  =  Gtl123*gtu31 + Gtl223*gtu32 + Gtl323*gtu33;
    
    CCTK_REAL const Gt133  =  Gtl133*gtu11 + Gtl233*gtu21 + Gtl333*gtu31;
    
    CCTK_REAL const Gt233  =  Gtl133*gtu21 + Gtl233*gtu22 + Gtl333*gtu32;
    
    CCTK_REAL const Gt333  =  Gtl133*gtu31 + Gtl233*gtu32 + Gtl333*gtu33;
    
    CCTK_REAL const Xtn1  =  Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    CCTK_REAL const Xtn2  =  Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    CCTK_REAL const Xtn3  =  Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    CCTK_REAL const Rt11  =  3*(Gt111*Gtlu111 + Gt112*Gtlu112 + Gt113*Gtlu113) + 
        2*(Gt211*Gtlu121 + Gt212*Gtlu122 + Gt213*Gtlu123 + Gt311*Gtlu131 + Gt312*Gtlu132 + Gt313*Gtlu133) + Gt211*Gtlu211 + 
        Gt212*Gtlu212 + Gt213*Gtlu213 + Gt311*Gtlu311 + Gt312*Gtlu312 + Gt313*Gtlu313 + gt11L*PDstandardNth1Xt1 + 
        gt12L*PDstandardNth1Xt2 + gt13L*PDstandardNth1Xt3 + 
        khalf*(-(gtu11*PDstandardNth11gt11) - 2*gtu21*PDstandardNth12gt11 - 2*gtu31*PDstandardNth13gt11 - 
           gtu22*PDstandardNth22gt11 - 2*gtu32*PDstandardNth23gt11 - gtu33*PDstandardNth33gt11) + Gtl111*Xtn1 + 
        Gtl112*Xtn2 + Gtl113*Xtn3;
    
    CCTK_REAL const Rt12  =  khalf*(4*(Gt211*Gtlu221 + Gt212*Gtlu222 + Gt213*Gtlu223) + 
          2*(Gt112*Gtlu111 + Gt122*Gtlu112 + Gt123*Gtlu113 + Gt111*Gtlu121 + Gt212*Gtlu121 + Gt112*Gtlu122 + 
             Gt222*Gtlu122 + Gt113*Gtlu123 + Gt223*Gtlu123 + Gt312*Gtlu131 + Gt322*Gtlu132 + Gt323*Gtlu133 + 
             Gt111*Gtlu211 + Gt112*Gtlu212 + Gt113*Gtlu213 + Gt311*Gtlu231 + Gt312*Gtlu232 + Gt313*Gtlu233 + 
             Gt311*Gtlu321 + Gt312*Gtlu322 + Gt313*Gtlu323) - gtu11*PDstandardNth11gt12 - 2*gtu21*PDstandardNth12gt12 - 
          2*gtu31*PDstandardNth13gt12 + gt12L*PDstandardNth1Xt1 + gt22L*PDstandardNth1Xt2 + gt23L*PDstandardNth1Xt3 - 
          gtu22*PDstandardNth22gt12 - 2*gtu32*PDstandardNth23gt12 + gt11L*PDstandardNth2Xt1 + gt12L*PDstandardNth2Xt2 + 
          gt13L*PDstandardNth2Xt3 - gtu33*PDstandardNth33gt12 + Gtl112*Xtn1 + Gtl211*Xtn1 + Gtl122*Xtn2 + Gtl212*Xtn2 + 
          Gtl123*Xtn3 + Gtl213*Xtn3);
    
    CCTK_REAL const Rt13  =  khalf*(2*(Gt113*Gtlu111 + Gt123*Gtlu112 + Gt133*Gtlu113 + Gt213*Gtlu121 + Gt223*Gtlu122 + Gt233*Gtlu123 + 
             Gt111*Gtlu131 + Gt313*Gtlu131 + Gt112*Gtlu132 + Gt323*Gtlu132 + Gt113*Gtlu133 + Gt333*Gtlu133 + 
             Gt211*Gtlu231 + Gt212*Gtlu232 + Gt213*Gtlu233 + Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*Gtlu313 + 
             Gt211*Gtlu321 + Gt212*Gtlu322 + Gt213*Gtlu323) + 4*(Gt311*Gtlu331 + Gt312*Gtlu332 + Gt313*Gtlu333) - 
          gtu11*PDstandardNth11gt13 - 2*gtu21*PDstandardNth12gt13 - 2*gtu31*PDstandardNth13gt13 + gt13L*PDstandardNth1Xt1 + 
          gt23L*PDstandardNth1Xt2 + gt33L*PDstandardNth1Xt3 - gtu22*PDstandardNth22gt13 - 2*gtu32*PDstandardNth23gt13 - 
          gtu33*PDstandardNth33gt13 + gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + gt13L*PDstandardNth3Xt3 + 
          Gtl113*Xtn1 + Gtl311*Xtn1 + Gtl123*Xtn2 + Gtl312*Xtn2 + Gtl133*Xtn3 + Gtl313*Xtn3);
    
    CCTK_REAL const Rt22  =  Gt112*(Gtlu121 + 2*Gtlu211) + Gt122*(Gtlu122 + 2*Gtlu212) + Gt123*(Gtlu123 + 2*Gtlu213) + 
        3*(Gt212*Gtlu221 + Gt222*Gtlu222 + Gt223*Gtlu223) + 2*(Gt312*Gtlu231 + Gt322*Gtlu232 + Gt323*Gtlu233) + 
        Gt312*Gtlu321 + Gt322*Gtlu322 + Gt323*Gtlu323 + gt12L*PDstandardNth2Xt1 + gt22L*PDstandardNth2Xt2 + 
        gt23L*PDstandardNth2Xt3 + khalf*(-(gtu11*PDstandardNth11gt22) - 2*gtu21*PDstandardNth12gt22 - 
           2*gtu31*PDstandardNth13gt22 - gtu22*PDstandardNth22gt22 - 2*gtu32*PDstandardNth23gt22 - gtu33*PDstandardNth33gt22
           ) + Gtl212*Xtn1 + Gtl222*Xtn2 + Gtl223*Xtn3;
    
    CCTK_REAL const Rt23  =  khalf*(2*(Gt112*Gtlu131 + Gt122*Gtlu132 + Gt123*Gtlu133 + Gt113*Gtlu211 + Gt123*Gtlu212 + Gt133*Gtlu213 + 
             Gt213*Gtlu221 + Gt223*Gtlu222 + Gt233*Gtlu223 + Gt212*Gtlu231 + Gt313*Gtlu231 + Gt222*Gtlu232 + 
             Gt323*Gtlu232 + Gt223*Gtlu233 + Gt333*Gtlu233 + Gt112*Gtlu311 + Gt122*Gtlu312 + Gt123*Gtlu313 + 
             Gt212*Gtlu321 + Gt222*Gtlu322 + Gt223*Gtlu323) + 4*(Gt312*Gtlu331 + Gt322*Gtlu332 + Gt323*Gtlu333) - 
          gtu11*PDstandardNth11gt23 - 2*gtu21*PDstandardNth12gt23 - 2*gtu31*PDstandardNth13gt23 - 
          gtu22*PDstandardNth22gt23 - 2*gtu32*PDstandardNth23gt23 + gt13L*PDstandardNth2Xt1 + gt23L*PDstandardNth2Xt2 + 
          gt33L*PDstandardNth2Xt3 - gtu33*PDstandardNth33gt23 + gt12L*PDstandardNth3Xt1 + gt22L*PDstandardNth3Xt2 + 
          gt23L*PDstandardNth3Xt3 + Gtl213*Xtn1 + Gtl312*Xtn1 + Gtl223*Xtn2 + Gtl322*Xtn2 + Gtl233*Xtn3 + Gtl323*Xtn3);
    
    CCTK_REAL const Rt33  =  Gt113*(Gtlu131 + 2*Gtlu311) + Gt123*(Gtlu132 + 2*Gtlu312) + Gt133*(Gtlu133 + 2*Gtlu313) + 
        Gt213*(Gtlu231 + 2*Gtlu321) + Gt223*(Gtlu232 + 2*Gtlu322) + Gt233*(Gtlu233 + 2*Gtlu323) + 
        3*(Gt313*Gtlu331 + Gt323*Gtlu332 + Gt333*Gtlu333) + 
        khalf*(-(gtu11*PDstandardNth11gt33) - 2*gtu21*PDstandardNth12gt33 - 2*gtu31*PDstandardNth13gt33 - 
           gtu22*PDstandardNth22gt33 - 2*gtu32*PDstandardNth23gt33 - gtu33*PDstandardNth33gt33) + gt13L*PDstandardNth3Xt1 + 
        gt23L*PDstandardNth3Xt2 + gt33L*PDstandardNth3Xt3 + Gtl313*Xtn1 + Gtl323*Xtn2 + Gtl333*Xtn3;
    
    CCTK_REAL const fac1  =  IfThen(conformalMethod,-(khalf*INV(phiL)),1);
    
    CCTK_REAL const cdphi1  =  fac1*PDstandardNth1phi;
    
    CCTK_REAL const cdphi2  =  fac1*PDstandardNth2phi;
    
    CCTK_REAL const cdphi3  =  fac1*PDstandardNth3phi;
    
    CCTK_REAL const fac2  =  IfThen(conformalMethod,khalf*pow(phiL,-2),0);
    
    CCTK_REAL const cdphi211  =  -(fac1*(-PDstandardNth11phi + Gt111*PDstandardNth1phi + Gt211*PDstandardNth2phi + 
             Gt311*PDstandardNth3phi)) + fac2*SQR(PDstandardNth1phi);
    
    CCTK_REAL const cdphi212  =  fac2*PDstandardNth1phi*PDstandardNth2phi - 
        fac1*(-PDstandardNth12phi + Gt112*PDstandardNth1phi + Gt212*PDstandardNth2phi + Gt312*PDstandardNth3phi);
    
    CCTK_REAL const cdphi213  =  fac2*PDstandardNth1phi*PDstandardNth3phi - 
        fac1*(-PDstandardNth13phi + Gt113*PDstandardNth1phi + Gt213*PDstandardNth2phi + Gt313*PDstandardNth3phi);
    
    CCTK_REAL const cdphi222  =  -(fac1*(Gt122*PDstandardNth1phi - PDstandardNth22phi + Gt222*PDstandardNth2phi + 
             Gt322*PDstandardNth3phi)) + fac2*SQR(PDstandardNth2phi);
    
    CCTK_REAL const cdphi223  =  fac2*PDstandardNth2phi*PDstandardNth3phi - 
        fac1*(Gt123*PDstandardNth1phi - PDstandardNth23phi + Gt223*PDstandardNth2phi + Gt323*PDstandardNth3phi);
    
    CCTK_REAL const cdphi233  =  -(fac1*(Gt133*PDstandardNth1phi + Gt233*PDstandardNth2phi - PDstandardNth33phi + 
             Gt333*PDstandardNth3phi)) + fac2*SQR(PDstandardNth3phi);
    
    CCTK_REAL const Rphi11  =  -2*(cdphi211 + 2*(-1 + gt11L*gtu11)*SQR(cdphi1) + 
          gt11L*(cdphi211*gtu11 + 4*(cdphi1*(cdphi2*gtu21 + cdphi3*gtu31) + cdphi2*cdphi3*gtu32) + cdphi233*gtu33 + 
             gtu22*(cdphi222 + 2*SQR(cdphi2)) + 2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu33*SQR(cdphi3))));
    
    CCTK_REAL const Rphi12  =  -2*(cdphi212 + cdphi1*(cdphi2*(-2 + 4*gt12L*gtu21) + 4*cdphi3*gt12L*gtu31) + 
          gt12L*(cdphi211*gtu11 + 4*cdphi2*cdphi3*gtu32 + 
             2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu11*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2)) + 
             gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL const Rphi13  =  -2*(cdphi213 + cdphi1*(4*cdphi2*gt13L*gtu21 + cdphi3*(-2 + 4*gt13L*gtu31)) + 
          gt13L*(cdphi211*gtu11 + 4*cdphi2*cdphi3*gtu32 + 
             2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu11*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2)) + 
             gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL const Rphi22  =  -2*(cdphi222 + 2*(-1 + gt22L*gtu22)*SQR(cdphi2) + 
          gt22L*(cdphi222*gtu22 + 4*(cdphi1*cdphi3*gtu31 + cdphi2*(cdphi1*gtu21 + cdphi3*gtu32)) + cdphi233*gtu33 + 
             gtu11*(cdphi211 + 2*SQR(cdphi1)) + 2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu33*SQR(cdphi3))));
    
    CCTK_REAL const Rphi23  =  -2*(cdphi223 + cdphi2*(4*cdphi1*gt23L*gtu21 + cdphi3*(-2 + 4*gt23L*gtu32)) + 
          gt23L*(cdphi222*gtu22 + 4*cdphi1*cdphi3*gtu31 + gtu11*(cdphi211 + 2*SQR(cdphi1)) + 
             2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu22*SQR(cdphi2)) + gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL const Rphi33  =  -2*(cdphi233 + gt33L*((4*cdphi1*cdphi2 + 2*cdphi212)*gtu21 + 4*cdphi3*(cdphi1*gtu31 + cdphi2*gtu32) + 
             2*(cdphi213*gtu31 + cdphi223*gtu32) + cdphi233*gtu33 + gtu11*(cdphi211 + 2*SQR(cdphi1)) + 
             gtu22*(cdphi222 + 2*SQR(cdphi2))) + 2*(-1 + gt33L*gtu33)*SQR(cdphi3));
    
    CCTK_REAL const Atm11  =  At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    CCTK_REAL const Atm21  =  At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    CCTK_REAL const Atm31  =  At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    CCTK_REAL const Atm12  =  At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    CCTK_REAL const Atm22  =  At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    CCTK_REAL const Atm32  =  At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    CCTK_REAL const Atm13  =  At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    CCTK_REAL const Atm23  =  At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    CCTK_REAL const Atm33  =  At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    CCTK_REAL const e4phi  =  IfThen(conformalMethod,pow(phiL,-2),exp(4*phiL));
    
    CCTK_REAL const em4phi  =  INV(e4phi);
    
    CCTK_REAL const g11  =  e4phi*gt11L;
    
    CCTK_REAL const g12  =  e4phi*gt12L;
    
    CCTK_REAL const g13  =  e4phi*gt13L;
    
    CCTK_REAL const g22  =  e4phi*gt22L;
    
    CCTK_REAL const g23  =  e4phi*gt23L;
    
    CCTK_REAL const g33  =  e4phi*gt33L;
    
    CCTK_REAL const gu11  =  em4phi*gtu11;
    
    CCTK_REAL const gu21  =  em4phi*gtu21;
    
    CCTK_REAL const gu31  =  em4phi*gtu31;
    
    CCTK_REAL const gu22  =  em4phi*gtu22;
    
    CCTK_REAL const gu32  =  em4phi*gtu32;
    
    CCTK_REAL const gu33  =  em4phi*gtu33;
    
    CCTK_REAL const G111  =  Gt111 + cdphi1*(4 - 2*gt11L*gtu11) - 2*gt11L*(cdphi2*gtu21 + cdphi3*gtu31);
    
    CCTK_REAL const G211  =  Gt211 - 2*gt11L*(cdphi1*gtu21 + cdphi2*gtu22 + cdphi3*gtu32);
    
    CCTK_REAL const G311  =  Gt311 - 2*gt11L*(cdphi1*gtu31 + cdphi2*gtu32 + cdphi3*gtu33);
    
    CCTK_REAL const G112  =  Gt112 + cdphi2*(2 - 2*gt12L*gtu21) - 2*gt12L*(cdphi1*gtu11 + cdphi3*gtu31);
    
    CCTK_REAL const G212  =  Gt212 + cdphi1*(2 - 2*gt12L*gtu21) - 2*gt12L*(cdphi2*gtu22 + cdphi3*gtu32);
    
    CCTK_REAL const G312  =  Gt312 - 2*gt12L*(cdphi1*gtu31 + cdphi2*gtu32 + cdphi3*gtu33);
    
    CCTK_REAL const G113  =  Gt113 - 2*gt13L*(cdphi1*gtu11 + cdphi2*gtu21) + cdphi3*(2 - 2*gt13L*gtu31);
    
    CCTK_REAL const G213  =  Gt213 - 2*gt13L*(cdphi1*gtu21 + cdphi2*gtu22 + cdphi3*gtu32);
    
    CCTK_REAL const G313  =  Gt313 + cdphi1*(2 - 2*gt13L*gtu31) - 2*gt13L*(cdphi2*gtu32 + cdphi3*gtu33);
    
    CCTK_REAL const G122  =  Gt122 - 2*gt22L*(cdphi1*gtu11 + cdphi2*gtu21 + cdphi3*gtu31);
    
    CCTK_REAL const G222  =  Gt222 + cdphi2*(4 - 2*gt22L*gtu22) - 2*gt22L*(cdphi1*gtu21 + cdphi3*gtu32);
    
    CCTK_REAL const G322  =  Gt322 - 2*gt22L*(cdphi1*gtu31 + cdphi2*gtu32 + cdphi3*gtu33);
    
    CCTK_REAL const G123  =  Gt123 - 2*gt23L*(cdphi1*gtu11 + cdphi2*gtu21 + cdphi3*gtu31);
    
    CCTK_REAL const G223  =  Gt223 - 2*gt23L*(cdphi1*gtu21 + cdphi2*gtu22) + cdphi3*(2 - 2*gt23L*gtu32);
    
    CCTK_REAL const G323  =  Gt323 + cdphi2*(2 - 2*gt23L*gtu32) - 2*gt23L*(cdphi1*gtu31 + cdphi3*gtu33);
    
    CCTK_REAL const G133  =  Gt133 - 2*gt33L*(cdphi1*gtu11 + cdphi2*gtu21 + cdphi3*gtu31);
    
    CCTK_REAL const G233  =  Gt233 - 2*gt33L*(cdphi1*gtu21 + cdphi2*gtu22 + cdphi3*gtu32);
    
    CCTK_REAL const G333  =  Gt333 - 2*gt33L*(cdphi1*gtu31 + cdphi2*gtu32) + cdphi3*(4 - 2*gt33L*gtu33);
    
    CCTK_REAL const R11  =  Rphi11 + Rt11;
    
    CCTK_REAL const R12  =  Rphi12 + Rt12;
    
    CCTK_REAL const R13  =  Rphi13 + Rt13;
    
    CCTK_REAL const R22  =  Rphi22 + Rt22;
    
    CCTK_REAL const R23  =  Rphi23 + Rt23;
    
    CCTK_REAL const R33  =  Rphi33 + Rt33;
    
    CCTK_REAL const rho  =  pow(alphaL,-2)*(eTttL - 2*(beta2L*eTtyL + beta3L*eTtzL) + 
          2*(beta1L*(-eTtxL + beta2L*eTxyL + beta3L*eTxzL) + beta2L*beta3L*eTyzL) + eTxxL*SQR(beta1L) + eTyyL*SQR(beta2L) + 
          eTzzL*SQR(beta3L));
    
    CCTK_REAL const trS  =  eTxxL*gu11 + eTyyL*gu22 + 2*(eTxyL*gu21 + eTxzL*gu31 + eTyzL*gu32) + eTzzL*gu33;
    
    CCTK_REAL const trKrhsL  =  PDupwindNth1(trK, i, j, k)*beta1L + PDupwindNth2(trK, i, j, k)*beta2L + 
        PDupwindNth3(trK, i, j, k)*beta3L + (G111*gu11 + G122*gu22 + 2.*(G112*gu21 + G113*gu31 + G123*gu32) + G133*gu33)*
         PDstandardNth1alpha - 2.*(gu21*PDstandardNth12alpha + gu31*PDstandardNth13alpha + gu32*PDstandardNth23alpha) + 
        (G211*gu11 + G222*gu22 + 2.*(G212*gu21 + G213*gu31 + G223*gu32) + G233*gu33)*PDstandardNth2alpha - 
        1.*(gu11*PDstandardNth11alpha + gu22*PDstandardNth22alpha + gu33*PDstandardNth33alpha) + 
        (G311*gu11 + G322*gu22 + 2.*(G313*gu31 + G323*gu32) + G333*gu33)*PDstandardNth3alpha + 
        2.*(alphaL*(Atm12*Atm21 + Atm13*Atm31 + Atm23*Atm32) + G312*gu21*PDstandardNth3alpha) + 
        alphaL*(12.56637061435917295385057353311801153679*(rho + trS) + SQR(Atm11) + SQR(Atm22) + SQR(Atm33) + 
           0.3333333333333333333333333333333333333333*SQR(trKL));
    
    CCTK_REAL const Ats11  =  -PDstandardNth11alpha + G111*PDstandardNth1alpha + G211*PDstandardNth2alpha + G311*PDstandardNth3alpha + 
        alphaL*R11;
    
    CCTK_REAL const Ats12  =  -PDstandardNth12alpha + G112*PDstandardNth1alpha + G212*PDstandardNth2alpha + G312*PDstandardNth3alpha + 
        alphaL*R12;
    
    CCTK_REAL const Ats13  =  -PDstandardNth13alpha + G113*PDstandardNth1alpha + G213*PDstandardNth2alpha + G313*PDstandardNth3alpha + 
        alphaL*R13;
    
    CCTK_REAL const Ats22  =  G122*PDstandardNth1alpha - PDstandardNth22alpha + G222*PDstandardNth2alpha + G322*PDstandardNth3alpha + 
        alphaL*R22;
    
    CCTK_REAL const Ats23  =  G123*PDstandardNth1alpha - PDstandardNth23alpha + G223*PDstandardNth2alpha + G323*PDstandardNth3alpha + 
        alphaL*R23;
    
    CCTK_REAL const Ats33  =  G133*PDstandardNth1alpha + G233*PDstandardNth2alpha - PDstandardNth33alpha + G333*PDstandardNth3alpha + 
        alphaL*R33;
    
    CCTK_REAL const trAts  =  Ats11*gu11 + Ats22*gu22 + 2*(Ats12*gu21 + Ats13*gu31 + Ats23*gu32) + Ats33*gu33;
    
    CCTK_REAL const At11rhsL  =  -2.*alphaL*(At11L*Atm11 + At12L*Atm21 + At13L*Atm31) + PDupwindNth1(At11, i, j, k)*beta1L + 
        PDupwindNth2(At11, i, j, k)*beta2L + PDupwindNth3(At11, i, j, k)*beta3L + 
        2.*(At12L*PDstandardNth1beta2 + At13L*PDstandardNth1beta3) + 
        At11L*(1.333333333333333333333333333333333333333*PDstandardNth1beta1 - 
           0.6666666666666666666666666666666666666667*(PDstandardNth2beta2 + PDstandardNth3beta3) + alphaL*trKL) + 
        em4phi*(Ats11 - 0.3333333333333333333333333333333333333333*g11*trAts + 
           alphaL*(-25.13274122871834590770114706623602307358*eTxxL + 8.377580409572781969233715688745341024526*g11*trS));
    
    CCTK_REAL const At12rhsL  =  -2.*alphaL*(At11L*Atm12 + At12L*Atm22 + At13L*Atm32) + PDupwindNth1(At12, i, j, k)*beta1L + 
        PDupwindNth2(At12, i, j, k)*beta2L + PDupwindNth3(At12, i, j, k)*beta3L + At22L*PDstandardNth1beta2 + 
        At23L*PDstandardNth1beta3 + At11L*PDstandardNth2beta1 + At13L*PDstandardNth2beta3 + 
        At12L*(0.3333333333333333333333333333333333333333*(PDstandardNth1beta1 + PDstandardNth2beta2) - 
           0.6666666666666666666666666666666666666667*PDstandardNth3beta3 + alphaL*trKL) + 
        em4phi*(Ats12 - 0.3333333333333333333333333333333333333333*g12*trAts + 
           alphaL*(-25.13274122871834590770114706623602307358*eTxyL + 8.377580409572781969233715688745341024526*g12*trS));
    
    CCTK_REAL const At13rhsL  =  -2.*alphaL*(At11L*Atm13 + At12L*Atm23 + At13L*Atm33) + PDupwindNth1(At13, i, j, k)*beta1L + 
        PDupwindNth2(At13, i, j, k)*beta2L + PDupwindNth3(At13, i, j, k)*beta3L + At23L*PDstandardNth1beta2 + 
        At33L*PDstandardNth1beta3 + At11L*PDstandardNth3beta1 + At12L*PDstandardNth3beta2 + 
        At13L*(-0.6666666666666666666666666666666666666667*PDstandardNth2beta2 + 
           0.3333333333333333333333333333333333333333*(PDstandardNth1beta1 + PDstandardNth3beta3) + alphaL*trKL) + 
        em4phi*(Ats13 - 0.3333333333333333333333333333333333333333*g13*trAts + 
           alphaL*(-25.13274122871834590770114706623602307358*eTxzL + 8.377580409572781969233715688745341024526*g13*trS));
    
    CCTK_REAL const At22rhsL  =  -2.*alphaL*(At12L*Atm12 + At22L*Atm22 + At23L*Atm32) + PDupwindNth1(At22, i, j, k)*beta1L + 
        PDupwindNth2(At22, i, j, k)*beta2L + PDupwindNth3(At22, i, j, k)*beta3L + 
        2.*(At12L*PDstandardNth2beta1 + At23L*PDstandardNth2beta3) + 
        At22L*(1.333333333333333333333333333333333333333*PDstandardNth2beta2 - 
           0.6666666666666666666666666666666666666667*(PDstandardNth1beta1 + PDstandardNth3beta3) + alphaL*trKL) + 
        em4phi*(Ats22 - 0.3333333333333333333333333333333333333333*g22*trAts + 
           alphaL*(-25.13274122871834590770114706623602307358*eTyyL + 8.377580409572781969233715688745341024526*g22*trS));
    
    CCTK_REAL const At23rhsL  =  -2.*alphaL*(At12L*Atm13 + At22L*Atm23 + At23L*Atm33) + PDupwindNth1(At23, i, j, k)*beta1L + 
        PDupwindNth2(At23, i, j, k)*beta2L + PDupwindNth3(At23, i, j, k)*beta3L + At13L*PDstandardNth2beta1 + 
        At33L*PDstandardNth2beta3 + At12L*PDstandardNth3beta1 + At22L*PDstandardNth3beta2 + 
        At23L*(-0.6666666666666666666666666666666666666667*PDstandardNth1beta1 + 
           0.3333333333333333333333333333333333333333*(PDstandardNth2beta2 + PDstandardNth3beta3) + alphaL*trKL) + 
        em4phi*(Ats23 - 0.3333333333333333333333333333333333333333*g23*trAts + 
           alphaL*(-25.13274122871834590770114706623602307358*eTyzL + 8.377580409572781969233715688745341024526*g23*trS));
    
    CCTK_REAL const At33rhsL  =  -2.*alphaL*(At13L*Atm13 + At23L*Atm23 + At33L*Atm33) + PDupwindNth1(At33, i, j, k)*beta1L + 
        PDupwindNth2(At33, i, j, k)*beta2L + PDupwindNth3(At33, i, j, k)*beta3L + 
        2.*(At13L*PDstandardNth3beta1 + At23L*PDstandardNth3beta2) + 
        At33L*(-0.6666666666666666666666666666666666666667*(PDstandardNth1beta1 + PDstandardNth2beta2) + 
           1.333333333333333333333333333333333333333*PDstandardNth3beta3 + alphaL*trKL) + 
        em4phi*(Ats33 - 0.3333333333333333333333333333333333333333*g33*trAts + 
           alphaL*(-25.13274122871834590770114706623602307358*eTzzL + 8.377580409572781969233715688745341024526*g33*trS));
    
    CCTK_REAL const alpharhsL  =  (PDupwindNth1(alpha, i, j, k)*beta1L + PDupwindNth2(alpha, i, j, k)*beta2L + 
           PDupwindNth3(alpha, i, j, k)*beta3L)*LapseAdvectionCoeff + 
        harmonicF*(AL*(-1 + LapseAdvectionCoeff) - LapseAdvectionCoeff*trKL)*pow(alphaL,harmonicN);
    
    CCTK_REAL const ArhsL  =  (-1 + LapseAdvectionCoeff)*(AL*AlphaDriver - trKrhsL);
    
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
    trKrhs[index] = trKrhsL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSN_RHS2);
}

void ML_BSSN_RHS2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_RHS2_Body);
}
