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

void ML_BSSN_constraints_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_constraints_Body");
  }
  
  if (cctk_iteration % ML_BSSN_constraints_calc_every != ML_BSSN_constraints_calc_offset)
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
  LC_LOOP3 (ML_BSSN_constraints,
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
    // CCTK_REAL cdphi1 = INITVALUE, cdphi2 = INITVALUE, cdphi211 = INITVALUE, cdphi212 = INITVALUE, cdphi213 = INITVALUE, cdphi222 = INITVALUE;
    // CCTK_REAL cdphi223 = INITVALUE, cdphi233 = INITVALUE, cdphi3 = INITVALUE;
    // CCTK_REAL detgt = INITVALUE;
    // CCTK_REAL e4phi = INITVALUE;
    // CCTK_REAL em4phi = INITVALUE;
    // CCTK_REAL fac1 = INITVALUE, fac2 = INITVALUE;
    // CCTK_REAL Gt111 = INITVALUE, Gt112 = INITVALUE, Gt113 = INITVALUE, Gt122 = INITVALUE, Gt123 = INITVALUE, Gt133 = INITVALUE;
    // CCTK_REAL Gt211 = INITVALUE, Gt212 = INITVALUE, Gt213 = INITVALUE, Gt222 = INITVALUE, Gt223 = INITVALUE, Gt233 = INITVALUE;
    // CCTK_REAL Gt311 = INITVALUE, Gt312 = INITVALUE, Gt313 = INITVALUE, Gt322 = INITVALUE, Gt323 = INITVALUE, Gt333 = INITVALUE;
    // CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    // CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    // CCTK_REAL R11 = INITVALUE, R12 = INITVALUE, R13 = INITVALUE, R22 = INITVALUE, R23 = INITVALUE, R33 = INITVALUE;
    // CCTK_REAL rho = INITVALUE;
    // CCTK_REAL Rphi11 = INITVALUE, Rphi12 = INITVALUE, Rphi13 = INITVALUE, Rphi22 = INITVALUE, Rphi23 = INITVALUE, Rphi33 = INITVALUE;
    // CCTK_REAL Rt11 = INITVALUE, Rt12 = INITVALUE, Rt13 = INITVALUE, Rt22 = INITVALUE, Rt23 = INITVALUE, Rt33 = INITVALUE;
    // CCTK_REAL S1 = INITVALUE, S2 = INITVALUE, S3 = INITVALUE;
    // CCTK_REAL T00 = INITVALUE, T01 = INITVALUE, T02 = INITVALUE, T03 = INITVALUE, T11 = INITVALUE, T12 = INITVALUE;
    // CCTK_REAL T13 = INITVALUE, T22 = INITVALUE, T23 = INITVALUE, T33 = INITVALUE;
    // CCTK_REAL trR = INITVALUE;
    
    /* Declare local copies of grid functions */
    // CCTK_REAL alphaL = INITVALUE;
    // CCTK_REAL At11L = INITVALUE, At12L = INITVALUE, At13L = INITVALUE, At22L = INITVALUE, At23L = INITVALUE, At33L = INITVALUE;
    // CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    // CCTK_REAL cAL = INITVALUE;
    // CCTK_REAL cSL = INITVALUE;
    // CCTK_REAL cXt1L = INITVALUE, cXt2L = INITVALUE, cXt3L = INITVALUE;
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
    // CCTK_REAL HL = INITVALUE;
    // CCTK_REAL M1L = INITVALUE, M2L = INITVALUE, M3L = INITVALUE;
    // CCTK_REAL phiL = INITVALUE;
    // CCTK_REAL trKL = INITVALUE;
    // CCTK_REAL Xt1L = INITVALUE, Xt2L = INITVALUE, Xt3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth1At11 = INITVALUE;
    // CCTK_REAL PDstandardNth2At11 = INITVALUE;
    // CCTK_REAL PDstandardNth3At11 = INITVALUE;
    // CCTK_REAL PDstandardNth1At12 = INITVALUE;
    // CCTK_REAL PDstandardNth2At12 = INITVALUE;
    // CCTK_REAL PDstandardNth3At12 = INITVALUE;
    // CCTK_REAL PDstandardNth1At13 = INITVALUE;
    // CCTK_REAL PDstandardNth2At13 = INITVALUE;
    // CCTK_REAL PDstandardNth3At13 = INITVALUE;
    // CCTK_REAL PDstandardNth1At22 = INITVALUE;
    // CCTK_REAL PDstandardNth2At22 = INITVALUE;
    // CCTK_REAL PDstandardNth3At22 = INITVALUE;
    // CCTK_REAL PDstandardNth1At23 = INITVALUE;
    // CCTK_REAL PDstandardNth2At23 = INITVALUE;
    // CCTK_REAL PDstandardNth3At23 = INITVALUE;
    // CCTK_REAL PDstandardNth1At33 = INITVALUE;
    // CCTK_REAL PDstandardNth2At33 = INITVALUE;
    // CCTK_REAL PDstandardNth3At33 = INITVALUE;
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
    CCTK_REAL const eTttL = (stress_energy_state) ? (eTtt[index]) : 0.0;
    CCTK_REAL const eTtxL = (stress_energy_state) ? (eTtx[index]) : 0.0;
    CCTK_REAL const eTtyL = (stress_energy_state) ? (eTty[index]) : 0.0;
    CCTK_REAL const eTtzL = (stress_energy_state) ? (eTtz[index]) : 0.0;
    CCTK_REAL const eTxxL = (stress_energy_state) ? (eTxx[index]) : 0.0;
    CCTK_REAL const eTxyL = (stress_energy_state) ? (eTxy[index]) : 0.0;
    CCTK_REAL const eTxzL = (stress_energy_state) ? (eTxz[index]) : 0.0;
    CCTK_REAL const eTyyL = (stress_energy_state) ? (eTyy[index]) : 0.0;
    CCTK_REAL const eTyzL = (stress_energy_state) ? (eTyz[index]) : 0.0;
    CCTK_REAL const eTzzL = (stress_energy_state) ? (eTzz[index]) : 0.0;
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
    CCTK_REAL const PDstandardNth1At11 = PDstandardNth1(At11, i, j, k);
    CCTK_REAL const PDstandardNth2At11 = PDstandardNth2(At11, i, j, k);
    CCTK_REAL const PDstandardNth3At11 = PDstandardNth3(At11, i, j, k);
    CCTK_REAL const PDstandardNth1At12 = PDstandardNth1(At12, i, j, k);
    CCTK_REAL const PDstandardNth2At12 = PDstandardNth2(At12, i, j, k);
    CCTK_REAL const PDstandardNth3At12 = PDstandardNth3(At12, i, j, k);
    CCTK_REAL const PDstandardNth1At13 = PDstandardNth1(At13, i, j, k);
    CCTK_REAL const PDstandardNth2At13 = PDstandardNth2(At13, i, j, k);
    CCTK_REAL const PDstandardNth3At13 = PDstandardNth3(At13, i, j, k);
    CCTK_REAL const PDstandardNth1At22 = PDstandardNth1(At22, i, j, k);
    CCTK_REAL const PDstandardNth2At22 = PDstandardNth2(At22, i, j, k);
    CCTK_REAL const PDstandardNth3At22 = PDstandardNth3(At22, i, j, k);
    CCTK_REAL const PDstandardNth1At23 = PDstandardNth1(At23, i, j, k);
    CCTK_REAL const PDstandardNth2At23 = PDstandardNth2(At23, i, j, k);
    CCTK_REAL const PDstandardNth3At23 = PDstandardNth3(At23, i, j, k);
    CCTK_REAL const PDstandardNth1At33 = PDstandardNth1(At33, i, j, k);
    CCTK_REAL const PDstandardNth2At33 = PDstandardNth2(At33, i, j, k);
    CCTK_REAL const PDstandardNth3At33 = PDstandardNth3(At33, i, j, k);
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
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL const detgt  =  1;
    
    CCTK_REAL const gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL const gtu21  =  (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL const gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL const gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL const gtu32  =  (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL const gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL const Gt111  =  khalf*(gtu11*PDstandardNth1gt11 + 2*(gtu21*PDstandardNth1gt12 + gtu31*PDstandardNth1gt13) - 
          gtu21*PDstandardNth2gt11 - gtu31*PDstandardNth3gt11);
    
    CCTK_REAL const Gt211  =  khalf*(gtu21*PDstandardNth1gt11 + 2*(gtu22*PDstandardNth1gt12 + gtu32*PDstandardNth1gt13) - 
          gtu22*PDstandardNth2gt11 - gtu32*PDstandardNth3gt11);
    
    CCTK_REAL const Gt311  =  khalf*(gtu31*PDstandardNth1gt11 + 2*(gtu32*PDstandardNth1gt12 + gtu33*PDstandardNth1gt13) - 
          gtu32*PDstandardNth2gt11 - gtu33*PDstandardNth3gt11);
    
    CCTK_REAL const Gt112  =  khalf*(gtu21*PDstandardNth1gt22 + gtu11*PDstandardNth2gt11 + 
          gtu31*(PDstandardNth1gt23 + PDstandardNth2gt13 - PDstandardNth3gt12));
    
    CCTK_REAL const Gt212  =  khalf*(gtu22*PDstandardNth1gt22 + gtu21*PDstandardNth2gt11 + 
          gtu32*(PDstandardNth1gt23 + PDstandardNth2gt13 - PDstandardNth3gt12));
    
    CCTK_REAL const Gt312  =  khalf*(gtu32*PDstandardNth1gt22 + gtu31*PDstandardNth2gt11 + 
          gtu33*(PDstandardNth1gt23 + PDstandardNth2gt13 - PDstandardNth3gt12));
    
    CCTK_REAL const Gt113  =  khalf*(gtu31*PDstandardNth1gt33 + gtu11*PDstandardNth3gt11 + 
          gtu21*(PDstandardNth1gt23 - PDstandardNth2gt13 + PDstandardNth3gt12));
    
    CCTK_REAL const Gt213  =  khalf*(gtu32*PDstandardNth1gt33 + gtu21*PDstandardNth3gt11 + 
          gtu22*(PDstandardNth1gt23 - PDstandardNth2gt13 + PDstandardNth3gt12));
    
    CCTK_REAL const Gt313  =  khalf*(gtu33*PDstandardNth1gt33 + gtu31*PDstandardNth3gt11 + 
          gtu32*(PDstandardNth1gt23 - PDstandardNth2gt13 + PDstandardNth3gt12));
    
    CCTK_REAL const Gt122  =  khalf*(gtu11*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + gtu21*PDstandardNth2gt22 + 
          gtu31*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    CCTK_REAL const Gt222  =  khalf*(gtu21*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + gtu22*PDstandardNth2gt22 + 
          gtu32*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    CCTK_REAL const Gt322  =  khalf*(gtu31*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + gtu32*PDstandardNth2gt22 + 
          gtu33*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    CCTK_REAL const Gt123  =  khalf*(gtu31*PDstandardNth2gt33 + gtu11*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) + 
          gtu21*PDstandardNth3gt22);
    
    CCTK_REAL const Gt223  =  khalf*(gtu32*PDstandardNth2gt33 + gtu21*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) + 
          gtu22*PDstandardNth3gt22);
    
    CCTK_REAL const Gt323  =  khalf*(gtu33*PDstandardNth2gt33 + gtu31*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) + 
          gtu32*PDstandardNth3gt22);
    
    CCTK_REAL const Gt133  =  khalf*(-(gtu11*PDstandardNth1gt33) - gtu21*PDstandardNth2gt33 + 2*gtu11*PDstandardNth3gt13 + 
          2*gtu21*PDstandardNth3gt23 + gtu31*PDstandardNth3gt33);
    
    CCTK_REAL const Gt233  =  khalf*(-(gtu21*PDstandardNth1gt33) - gtu22*PDstandardNth2gt33 + 2*gtu21*PDstandardNth3gt13 + 
          2*gtu22*PDstandardNth3gt23 + gtu32*PDstandardNth3gt33);
    
    CCTK_REAL const Gt333  =  khalf*(-(gtu31*PDstandardNth1gt33) - gtu32*PDstandardNth2gt33 + 2*gtu31*PDstandardNth3gt13 + 
          2*gtu32*PDstandardNth3gt23 + gtu33*PDstandardNth3gt33);
    
    CCTK_REAL const Rt11  =  -(gtu11*khalf*PDstandardNth11gt11) + gtu21*
         (2*Gt211*Gt212*gt22L + 4*Gt112*gt13L*Gt311 + 2*Gt113*gt11L*Gt312 + 2*gt13L*Gt312*Gt313 + 2*gt13L*Gt211*Gt322 + 
           2*gt13L*Gt311*Gt323 + 2*Gt311*Gt312*gt33L - PDstandardNth12gt11) - gtu31*PDstandardNth13gt11 + 
        gt11L*PDstandardNth1Xt1 + gt12L*(4*Gt111*Gt212*gtu21 + 2*Gt211*Gt222*gtu21 + 2*Gt212*Gt222*gtu22 + 
           4*Gt113*Gt211*gtu31 + 4*Gt113*Gt212*gtu32 + 4*Gt113*Gt213*gtu33 + PDstandardNth1Xt2) + 
        gt13L*(4*Gt111*Gt312*gtu21 + 2*Gt212*Gt312*gtu21 + 4*Gt112*Gt312*gtu22 + 4*Gt113*Gt311*gtu31 + 
           4*Gt113*Gt312*gtu32 + 4*Gt113*Gt313*gtu33 + PDstandardNth1Xt3) - gtu22*khalf*PDstandardNth22gt11 - 
        gtu32*PDstandardNth23gt11 - gtu33*khalf*PDstandardNth33gt11 + 
        Gt111*(6*Gt113*gt11L*gtu31 + 4*gt12L*Gt213*gtu31 + gt11L*Xt1L) + 
        Gt211*(2*Gt112*gt11L*gtu11 + 4*Gt111*gt12L*gtu11 + 2*gt11L*Gt122*gtu21 + 2*gt11L*Gt123*gtu31 + gt12L*Xt1L) + 
        Gt311*(4*Gt111*gt13L*gtu11 + 2*gt12L*Gt213*gtu11 + 2*gt13L*Gt313*gtu11 + 2*gt11L*Gt123*gtu21 + 
           2*gt11L*Gt133*gtu31 + gt13L*Xt1L) + gt12L*Gt212*Xt2L + gt13L*Gt312*Xt2L + 
        Gt112*(6*Gt111*gt11L*gtu21 + 4*gt12L*Gt211*gtu21 + 4*gt12L*Gt212*gtu22 + 2*gt11L*Gt213*gtu31 + 
           6*Gt113*gt11L*gtu32 + gt11L*Xt2L) + Gt113*gt11L*Xt3L + 
        Gt213*(2*gt11L*Gt122*gtu32 + 4*Gt112*gt12L*gtu32 + 2*gt11L*Gt123*gtu33 + gt12L*Xt3L) + 
        Gt313*(4*Gt111*gt13L*gtu31 + 2*gt12L*Gt213*gtu31 + 2*gt11L*Gt123*gtu32 + 4*Gt112*gt13L*gtu32 + 
           2*gt12L*Gt223*gtu32 + 2*gt11L*Gt133*gtu33 + gt13L*Xt3L) + 3*gt11L*gtu11*SQR(Gt111) + 3*gt11L*gtu22*SQR(Gt112) + 
        3*gt11L*gtu33*SQR(Gt113) + gt22L*gtu11*SQR(Gt211) + gt22L*gtu22*SQR(Gt212) + 
        2*(gt12L*Gt211*Gt212*gtu11 + Gt113*gt11L*Gt311*gtu11 + Gt211*gt23L*Gt311*gtu11 + gt13L*Gt211*Gt312*gtu11 + 
           Gt112*gt11L*Gt212*gtu21 + gt12L*Gt223*Gt311*gtu21 + Gt212*gt23L*Gt311*gtu21 + gt12L*Gt213*Gt312*gtu21 + 
           Gt211*gt23L*Gt312*gtu21 + gt11L*Gt122*Gt212*gtu22 + gt11L*Gt123*Gt312*gtu22 + gt12L*Gt223*Gt312*gtu22 + 
           Gt212*gt23L*Gt312*gtu22 + gt13L*Gt212*Gt322*gtu22 + gt13L*Gt312*Gt323*gtu22 + gt12L*Gt212*Gt213*gtu31 + 
           gt12L*Gt211*Gt223*gtu31 + Gt211*Gt213*gt22L*gtu31 + gt12L*Gt233*Gt311*gtu31 + Gt213*gt23L*Gt311*gtu31 + 
           gt13L*Gt213*Gt312*gtu31 + Gt113*gt11L*Gt313*gtu31 + Gt211*gt23L*Gt313*gtu31 + gt13L*Gt211*Gt323*gtu31 + 
           gt13L*Gt311*Gt333*gtu31 + Gt311*Gt313*gt33L*gtu31 + gt11L*Gt123*Gt212*gtu32 + gt12L*Gt213*Gt222*gtu32 + 
           gt12L*Gt212*Gt223*gtu32 + Gt212*Gt213*gt22L*gtu32 + gt11L*Gt133*Gt312*gtu32 + gt12L*Gt233*Gt312*gtu32 + 
           Gt213*gt23L*Gt312*gtu32 + Gt212*gt23L*Gt313*gtu32 + gt13L*Gt213*Gt322*gtu32 + gt13L*Gt212*Gt323*gtu32 + 
           gt13L*Gt313*Gt323*gtu32 + gt13L*Gt312*Gt333*gtu32 + Gt312*Gt313*gt33L*gtu32 + gt12L*Gt213*Gt223*gtu33 + 
           gt12L*Gt233*Gt313*gtu33 + Gt213*gt23L*Gt313*gtu33 + gt13L*Gt213*Gt323*gtu33 + gt13L*Gt313*Gt333*gtu33 + 
           gt12L*gtu21*SQR(Gt212)) + gt22L*gtu33*SQR(Gt213) + gt33L*gtu11*SQR(Gt311) + gt33L*gtu22*SQR(Gt312) + 
        2*gt13L*gtu31*SQR(Gt313) + gt33L*gtu33*SQR(Gt313);
    
    CCTK_REAL const Rt12  =  khalf*(-(gtu11*PDstandardNth11gt12) - 2*gtu21*PDstandardNth12gt12 - 2*gtu31*PDstandardNth13gt12 + 
          gt12L*PDstandardNth1Xt1 + gt22L*PDstandardNth1Xt2 + gt23L*PDstandardNth1Xt3 - gtu22*PDstandardNth22gt12 - 
          2*gtu32*PDstandardNth23gt12 + gt11L*PDstandardNth2Xt1 + gt12L*PDstandardNth2Xt2 + gt13L*PDstandardNth2Xt3 - 
          gtu33*PDstandardNth33gt12 + (Gt111*gt12L + Gt211*gt22L + gt23L*Gt311)*Xt1L + 
          (Gt112*gt11L + gt12L*Gt212 + gt13L*Gt312)*Xt1L + (Gt112*gt12L + Gt212*gt22L + gt23L*Gt312)*Xt2L + 
          (gt11L*Gt122 + gt12L*Gt222 + gt13L*Gt322)*Xt2L + (Gt113*gt12L + Gt213*gt22L + gt23L*Gt313)*Xt3L + 
          (gt11L*Gt123 + gt12L*Gt223 + gt13L*Gt323)*Xt3L + 
          2*gtu21*(Gt112*gt11L*Gt222 + Gt112*Gt211*gt22L + Gt211*Gt222*gt22L + 2*Gt122*gt13L*Gt311 + Gt112*gt23L*Gt311 + 
             Gt222*gt23L*Gt311 + gt13L*Gt222*Gt312 + Gt213*gt22L*Gt312 + Gt212*gt23L*Gt312 + gt23L*Gt312*Gt313 + 
             Gt113*gt11L*Gt322 + Gt211*gt23L*Gt322 + gt13L*Gt313*Gt322 + 
             Gt111*(2*gt11L*Gt122 + Gt112*gt12L + gt12L*Gt222 + gt13L*Gt322) + 
             gt12L*(2*Gt122*Gt211 + Gt112*Gt212 + Gt212*Gt222 + Gt113*Gt312 + Gt213*Gt322) + Gt311*Gt322*gt33L + 
             gt22L*SQR(Gt212)) + 2*((Gt123*gt12L*Gt211 + Gt113*gt12L*Gt212 + 2*Gt112*gt12L*Gt213 + gt12L*Gt212*Gt223 + 
                Gt212*Gt213*gt22L + Gt211*Gt223*gt22L + gt12L*Gt133*Gt311 + gt22L*Gt233*Gt311 + Gt113*gt13L*Gt312 + 
                gt12L*Gt233*Gt312 + Gt213*gt23L*Gt312 + gt11L*(2*Gt112*Gt113 + Gt123*Gt212 + Gt133*Gt312) + 
                2*Gt112*gt13L*Gt313 + Gt212*gt23L*Gt313 + Gt111*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + 
                gt13L*Gt212*Gt323 + Gt211*gt23L*Gt323 + gt23L*Gt311*Gt333 + gt13L*Gt312*Gt333 + Gt312*Gt313*gt33L)*gtu31 + 
             (Gt123*gt12L*Gt212 + 2*Gt122*gt12L*Gt213 + Gt113*gt12L*Gt222 + gt12L*Gt222*Gt223 + Gt213*Gt222*gt22L + 
                Gt212*Gt223*gt22L + gt12L*Gt133*Gt312 + gt22L*Gt233*Gt312 + 2*Gt122*gt13L*Gt313 + Gt222*gt23L*Gt313 + 
                Gt112*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + Gt113*gt13L*Gt322 + gt12L*Gt233*Gt322 + 
                Gt213*gt23L*Gt322 + gt11L*(2*Gt113*Gt122 + Gt123*Gt222 + Gt133*Gt322) + gt13L*Gt222*Gt323 + 
                Gt212*gt23L*Gt323 + gt23L*Gt312*Gt333 + gt13L*Gt322*Gt333 + Gt313*Gt322*gt33L)*gtu32 + 
             gtu11*(3*Gt112*gt12L*Gt211 + 2*Gt211*Gt212*gt22L + Gt113*gt12L*Gt311 + 2*Gt112*gt13L*Gt311 + 
                Gt213*gt22L*Gt311 + Gt212*gt23L*Gt311 + gt13L*Gt212*Gt312 + gt12L*Gt213*Gt312 + 2*Gt211*gt23L*Gt312 + 
                gt11L*(2*Gt111*Gt112 + Gt112*Gt212 + Gt113*Gt312) + 
                Gt111*(gt12L*Gt212 + Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + gt23L*Gt311*Gt313 + gt13L*Gt312*Gt313 + 
                Gt311*Gt312*gt33L + gt12L*SQR(Gt111) + gt12L*SQR(Gt212))) + 
          2*gtu22*(gt11L*Gt122*Gt222 + 2*Gt212*Gt222*gt22L + 2*Gt122*gt13L*Gt312 + Gt223*gt22L*Gt312 + Gt222*gt23L*Gt312 + 
             gt11L*Gt123*Gt322 + gt13L*Gt222*Gt322 + 2*Gt212*gt23L*Gt322 + 
             Gt112*(2*gt11L*Gt122 + gt12L*Gt222 + Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322) + gt23L*Gt312*Gt323 + 
             gt13L*Gt322*Gt323 + Gt312*Gt322*gt33L + gt12L*SQR(Gt112) + 
             gt12L*(3*Gt122*Gt212 + Gt123*Gt312 + Gt223*Gt322 + SQR(Gt222))) + 
          2*gtu33*(gt11L*Gt123*Gt223 + 2*Gt213*Gt223*gt22L + 2*Gt123*gt13L*Gt313 + gt22L*Gt233*Gt313 + Gt223*gt23L*Gt313 + 
             gt11L*Gt133*Gt323 + gt13L*Gt223*Gt323 + 2*Gt213*gt23L*Gt323 + 
             Gt113*(2*gt11L*Gt123 + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323) + gt23L*Gt313*Gt333 + 
             gt13L*Gt323*Gt333 + Gt313*Gt323*gt33L + gt12L*SQR(Gt113) + 
             gt12L*(3*Gt123*Gt213 + Gt133*Gt313 + Gt233*Gt323 + SQR(Gt223))) + 
          2*gtu21*(Gt122*gt12L*Gt211 + 3*Gt112*gt12L*Gt212 + gt12L*Gt212*Gt222 + Gt211*Gt222*gt22L + Gt123*gt12L*Gt311 + 
             Gt223*gt22L*Gt311 + 3*Gt112*gt13L*Gt312 + gt12L*Gt223*Gt312 + 2*Gt212*gt23L*Gt312 + 
             Gt111*(Gt112*gt12L + Gt212*gt22L + gt23L*Gt312) + gt13L*Gt212*Gt322 + Gt211*gt23L*Gt322 + gt23L*Gt311*Gt323 + 
             gt13L*Gt312*Gt323 + gt11L*(Gt122*Gt212 + Gt123*Gt312 + 2*SQR(Gt112)) + gt22L*SQR(Gt212) + gt33L*SQR(Gt312)) + 
          2*gtu31*(Gt112*gt11L*Gt223 + Gt113*Gt211*gt22L + Gt212*Gt213*gt22L + Gt211*Gt223*gt22L + 2*Gt123*gt13L*Gt311 + 
             Gt113*gt23L*Gt311 + Gt223*gt23L*Gt311 + gt13L*Gt223*Gt312 + Gt213*gt23L*Gt312 + Gt213*gt22L*Gt313 + 
             Gt113*gt11L*Gt323 + Gt211*gt23L*Gt323 + gt13L*Gt313*Gt323 + 
             Gt111*(2*gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + gt13L*Gt323) + 
             gt12L*(2*Gt123*Gt211 + Gt112*Gt213 + Gt212*Gt223 + Gt113*Gt313 + Gt213*Gt323) + Gt311*Gt323*gt33L + 
             gt23L*SQR(Gt313)) + 2*gtu32*(gt11L*Gt122*Gt223 + Gt113*Gt212*gt22L + Gt213*Gt222*gt22L + Gt212*Gt223*gt22L + 
             2*Gt123*gt13L*Gt312 + Gt113*gt23L*Gt312 + Gt223*gt23L*Gt312 + Gt223*gt22L*Gt313 + gt13L*Gt223*Gt322 + 
             Gt213*gt23L*Gt322 + gt11L*Gt123*Gt323 + Gt212*gt23L*Gt323 + gt23L*Gt313*Gt323 + 
             Gt112*(2*gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + gt13L*Gt323) + 
             gt12L*(Gt122*Gt213 + Gt123*(2*Gt212 + Gt313) + Gt223*(Gt222 + Gt323)) + Gt312*Gt323*gt33L + gt13L*SQR(Gt323)));
    
    CCTK_REAL const Rt13  =  khalf*(-(gtu11*PDstandardNth11gt13) - 2*gtu21*PDstandardNth12gt13 - 2*gtu31*PDstandardNth13gt13 + 
          gt13L*PDstandardNth1Xt1 + gt23L*PDstandardNth1Xt2 + gt33L*PDstandardNth1Xt3 - gtu22*PDstandardNth22gt13 - 
          2*gtu32*PDstandardNth23gt13 - gtu33*PDstandardNth33gt13 + gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + 
          gt13L*PDstandardNth3Xt3 + (Gt113*gt11L + gt12L*Gt213 + gt13L*Gt313)*Xt1L + 
          (Gt111*gt13L + Gt211*gt23L + Gt311*gt33L)*Xt1L + (gt11L*Gt123 + gt12L*Gt223 + gt13L*Gt323)*Xt2L + 
          (Gt112*gt13L + Gt212*gt23L + Gt312*gt33L)*Xt2L + (gt11L*Gt133 + gt12L*Gt233 + gt13L*Gt333)*Xt3L + 
          (Gt113*gt13L + Gt213*gt23L + Gt313*gt33L)*Xt3L + 
          2*((Gt122*gt13L*Gt211 + 2*Gt113*gt12L*Gt212 + Gt112*gt12L*Gt213 + gt12L*Gt213*Gt222 + Gt212*Gt213*gt22L + 
                Gt211*Gt222*gt23L + Gt123*gt13L*Gt311 + Gt223*gt23L*Gt311 + 2*Gt113*gt13L*Gt312 + Gt213*gt23L*Gt312 + 
                Gt112*gt13L*Gt313 + gt12L*Gt223*Gt313 + Gt212*gt23L*Gt313 + 
                gt11L*(2*Gt112*Gt113 + Gt122*Gt213 + Gt123*Gt313) + gt13L*Gt213*Gt322 + gt13L*Gt313*Gt323 + 
                Gt312*Gt313*gt33L + Gt211*Gt322*gt33L + Gt311*Gt323*gt33L + Gt111*(Gt112*gt13L + Gt212*gt23L + Gt312*gt33L))
               *gtu21 + (Gt122*gt13L*Gt213 + gt11L*Gt122*Gt233 + Gt212*gt22L*Gt233 + Gt113*Gt212*gt23L + 
                Gt213*Gt222*gt23L + 2*Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + Gt123*gt13L*Gt313 + Gt223*gt23L*Gt313 + 
                gt13L*Gt233*Gt322 + gt11L*Gt123*Gt333 + Gt212*gt23L*Gt333 + gt13L*Gt323*Gt333 + 
                Gt112*(2*gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + gt13L*Gt333) + 
                gt12L*(2*Gt133*Gt212 + Gt222*Gt233 + Gt223*Gt333) + Gt113*Gt312*gt33L + Gt213*Gt322*gt33L + 
                Gt313*Gt323*gt33L + Gt312*Gt333*gt33L)*gtu32 + 
             gtu21*(2*Gt123*gt12L*Gt211 + Gt112*gt13L*Gt212 + gt12L*Gt212*Gt223 + Gt211*Gt223*gt22L + Gt112*Gt211*gt23L + 
                2*Gt123*gt13L*Gt311 + Gt223*gt23L*Gt311 + Gt113*gt13L*Gt312 + gt13L*Gt223*Gt312 + Gt213*gt23L*Gt312 + 
                gt12L*Gt213*Gt323 + Gt211*gt23L*Gt323 + gt13L*Gt313*Gt323 + 
                gt11L*(2*Gt111*Gt123 + Gt112*Gt223 + Gt113*Gt323) + Gt111*(Gt112*gt13L + gt12L*Gt223 + gt13L*Gt323) + 
                Gt112*Gt311*gt33L + Gt212*Gt312*gt33L + Gt312*Gt313*gt33L + Gt311*Gt323*gt33L + gt23L*SQR(Gt212))) + 
          2*gtu32*(Gt123*gt13L*Gt212 + 2*Gt123*gt12L*Gt213 + Gt113*gt12L*Gt223 + Gt213*Gt223*gt22L + Gt212*Gt223*gt23L + 
             Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + 2*Gt123*gt13L*Gt313 + Gt223*gt23L*Gt313 + Gt113*gt13L*Gt323 + 
             gt13L*Gt223*Gt323 + gt12L*Gt233*Gt323 + Gt213*gt23L*Gt323 + 
             gt11L*(2*Gt113*Gt123 + Gt123*Gt223 + Gt133*Gt323) + gt13L*Gt323*Gt333 + Gt212*Gt323*gt33L + 
             Gt313*Gt323*gt33L + Gt312*Gt333*gt33L + Gt112*(Gt113*gt13L + Gt213*gt23L + Gt313*gt33L) + gt12L*SQR(Gt223)) + 
          2*gtu11*(2*Gt113*gt12L*Gt211 + Gt112*gt13L*Gt211 + gt12L*Gt212*Gt213 + Gt211*Gt213*gt22L + Gt211*Gt212*gt23L + 
             3*Gt113*gt13L*Gt311 + 2*Gt213*gt23L*Gt311 + gt13L*Gt213*Gt312 + gt12L*Gt213*Gt313 + Gt211*gt23L*Gt313 + 
             gt11L*(2*Gt111*Gt113 + Gt112*Gt213 + Gt113*Gt313) + Gt211*Gt312*gt33L + 2*Gt311*Gt313*gt33L + 
             Gt111*(gt12L*Gt213 + Gt211*gt23L + gt13L*Gt313 + Gt311*gt33L) + gt13L*SQR(Gt111) + gt13L*SQR(Gt313)) + 
          2*gtu31*(Gt112*gt13L*Gt213 + Gt112*gt11L*Gt233 + Gt211*gt22L*Gt233 + Gt113*Gt211*gt23L + Gt212*Gt213*gt23L + 
             2*Gt133*gt13L*Gt311 + Gt233*gt23L*Gt311 + gt13L*Gt233*Gt312 + Gt113*gt13L*Gt313 + Gt213*gt23L*Gt313 + 
             Gt113*gt11L*Gt333 + Gt211*gt23L*Gt333 + gt13L*Gt313*Gt333 + 
             Gt111*(2*gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + gt13L*Gt333) + 
             gt12L*(2*Gt133*Gt211 + Gt212*Gt233 + Gt213*Gt333) + Gt113*Gt311*gt33L + Gt213*Gt312*gt33L + 
             Gt311*Gt333*gt33L + gt33L*SQR(Gt313)) + 
          2*gtu31*(Gt123*gt13L*Gt211 + 3*Gt113*gt12L*Gt213 + gt12L*Gt213*Gt223 + Gt211*Gt223*gt23L + Gt133*gt13L*Gt311 + 
             Gt233*gt23L*Gt311 + 3*Gt113*gt13L*Gt313 + gt12L*Gt233*Gt313 + 2*Gt213*gt23L*Gt313 + gt13L*Gt213*Gt323 + 
             gt13L*Gt313*Gt333 + Gt211*Gt323*gt33L + Gt311*Gt333*gt33L + Gt111*(Gt113*gt13L + Gt213*gt23L + Gt313*gt33L) + 
             gt11L*(Gt123*Gt213 + Gt133*Gt313 + 2*SQR(Gt113)) + gt22L*SQR(Gt213) + gt33L*SQR(Gt313)) + 
          2*gtu22*(2*Gt123*gt12L*Gt212 + Gt122*gt13L*Gt212 + gt12L*Gt222*Gt223 + Gt212*Gt223*gt22L + Gt212*Gt222*gt23L + 
             3*Gt123*gt13L*Gt312 + 2*Gt223*gt23L*Gt312 + gt13L*Gt223*Gt322 + gt12L*Gt223*Gt323 + Gt212*gt23L*Gt323 + 
             gt11L*(2*Gt112*Gt123 + Gt122*Gt223 + Gt123*Gt323) + Gt212*Gt322*gt33L + 2*Gt312*Gt323*gt33L + 
             Gt112*(gt12L*Gt223 + Gt212*gt23L + gt13L*Gt323 + Gt312*gt33L) + gt13L*SQR(Gt112) + gt13L*SQR(Gt323)) + 
          2*gtu33*(2*gt12L*Gt133*Gt213 + Gt123*gt13L*Gt213 + gt11L*Gt123*Gt233 + gt12L*Gt223*Gt233 + Gt213*gt22L*Gt233 + 
             Gt213*Gt223*gt23L + 3*Gt133*gt13L*Gt313 + 2*Gt233*gt23L*Gt313 + gt13L*Gt233*Gt323 + gt11L*Gt133*Gt333 + 
             gt12L*Gt233*Gt333 + Gt213*gt23L*Gt333 + Gt213*Gt323*gt33L + 2*Gt313*Gt333*gt33L + 
             Gt113*(2*gt11L*Gt133 + gt12L*Gt233 + Gt213*gt23L + gt13L*Gt333 + Gt313*gt33L) + gt13L*SQR(Gt113) + 
             gt13L*SQR(Gt333)));
    
    CCTK_REAL const Rt22  =  4*(Gt122*gt12L*Gt212*gtu21 + Gt112*gt12L*Gt222*gtu21 + Gt122*gt12L*Gt222*gtu22 + Gt123*gt12L*Gt212*gtu31 + 
           Gt123*gt12L*Gt222*gtu32 + Gt123*gt12L*Gt223*gtu33) - gtu11*khalf*PDstandardNth11gt22 + 
        gtu21*(6*Gt212*Gt222*gt22L + 2*Gt122*gt23L*Gt311 + 2*Gt122*gt13L*Gt312 + 4*Gt222*gt23L*Gt312 + 
           2*Gt113*gt12L*Gt322 + 2*gt23L*Gt312*Gt323 + 2*Gt312*Gt322*gt33L - PDstandardNth12gt22) + 
        gtu31*(6*Gt212*Gt223*gt22L + 2*Gt123*gt13L*Gt312 + 2*Gt112*gt23L*Gt313 + 2*Gt113*gt12L*Gt323 + 
           2*gt23L*Gt312*Gt333 + 2*Gt312*Gt323*gt33L - PDstandardNth13gt22) - gtu22*khalf*PDstandardNth22gt22 + 
        gtu32*(4*Gt122*gt12L*Gt223 + 2*Gt123*Gt212*gt22L + 2*gt12L*Gt133*Gt322 + 4*Gt223*gt23L*Gt322 + 
           2*Gt123*gt12L*Gt323 + 4*Gt222*gt23L*Gt323 + 2*gt23L*Gt322*Gt333 + 2*Gt322*Gt323*gt33L - PDstandardNth23gt22) + 
        gt12L*(2*Gt111*Gt123*gtu31 + 4*Gt112*Gt223*gtu31 + 2*Gt113*Gt122*gtu32 + 2*Gt113*Gt123*gtu33 + PDstandardNth2Xt1) + 
        gt22L*(2*Gt122*Gt213*gtu32 + 6*Gt222*Gt223*gtu32 + 2*Gt123*Gt213*gtu33 + PDstandardNth2Xt2) + 
        gt23L*(4*Gt212*Gt322*gtu21 + 2*Gt313*Gt322*gtu21 + 4*Gt222*Gt322*gtu22 + 2*Gt123*Gt311*gtu31 + 
           4*Gt212*Gt323*gtu31 + 2*Gt313*Gt323*gtu31 + 2*Gt122*Gt313*gtu32 + 2*Gt123*Gt313*gtu33 + 4*Gt223*Gt323*gtu33 + 
           2*Gt323*Gt333*gtu33 + PDstandardNth2Xt3) - gtu33*khalf*PDstandardNth33gt22 + Gt212*gt22L*Xt1L + 
        Gt112*(2*Gt111*gt12L*gtu11 + 4*gt12L*Gt212*gtu11 + 2*gt11L*Gt122*gtu21 + 2*Gt122*gt12L*gtu22 + 
           2*gt11L*Gt123*gtu31 + 2*Gt123*gt12L*gtu32 + gt12L*Xt1L) + 
        Gt312*(2*Gt213*gt22L*gtu11 + 4*Gt212*gt23L*gtu11 + 2*gt23L*Gt313*gtu11 + 2*Gt123*gt12L*gtu21 + 
           2*Gt122*gt23L*gtu22 + 2*gt12L*Gt133*gtu31 + 2*gt22L*Gt233*gtu31 + 4*Gt223*gt23L*gtu31 + 2*Gt123*gt23L*gtu32 + 
           gt23L*Xt1L) + Gt122*gt12L*Xt2L + Gt222*gt22L*Xt2L + gt23L*Gt322*Xt2L + Gt123*gt12L*Xt3L + Gt223*gt22L*Xt3L + 
        gt23L*Gt323*Xt3L + gt11L*gtu11*SQR(Gt112) + 2*
         (Gt112*Gt211*gt22L*gtu11 + Gt112*gt23L*Gt311*gtu11 + Gt113*gt12L*Gt312*gtu11 + Gt112*gt13L*Gt312*gtu11 + 
           Gt111*Gt122*gt12L*gtu21 + Gt122*Gt211*gt22L*gtu21 + Gt112*Gt212*gt22L*gtu21 + Gt223*gt22L*Gt312*gtu21 + 
           Gt112*gt23L*Gt312*gtu21 + Gt112*gt13L*Gt322*gtu21 + Gt213*gt22L*Gt322*gtu21 + Gt122*Gt212*gt22L*gtu22 + 
           Gt123*gt12L*Gt322*gtu22 + Gt122*gt13L*Gt322*gtu22 + Gt223*gt22L*Gt322*gtu22 + gt23L*Gt322*Gt323*gtu22 + 
           Gt112*Gt113*gt12L*gtu31 + Gt123*Gt211*gt22L*gtu31 + Gt112*Gt213*gt22L*gtu31 + Gt112*gt13L*Gt323*gtu31 + 
           Gt213*gt22L*Gt323*gtu31 + gt11L*Gt122*Gt123*gtu32 + Gt123*gt13L*Gt322*gtu32 + gt22L*Gt233*Gt322*gtu32 + 
           Gt122*gt13L*Gt323*gtu32 + Gt223*gt22L*Gt323*gtu32 + gt12L*Gt133*Gt323*gtu33 + Gt123*gt13L*Gt323*gtu33 + 
           gt22L*Gt233*Gt323*gtu33 + gt12L*gtu21*SQR(Gt112)) + gt11L*gtu22*SQR(Gt122) + gt11L*gtu33*SQR(Gt123) + 
        3*gt22L*gtu11*SQR(Gt212) + 3*gt22L*gtu22*SQR(Gt222) + 3*gt22L*gtu33*SQR(Gt223) + gt33L*gtu11*SQR(Gt312) + 
        gt33L*gtu22*SQR(Gt322) + 2*gt23L*gtu32*SQR(Gt323) + gt33L*gtu33*SQR(Gt323);
    
    CCTK_REAL const Rt23  =  khalf*(-(gtu11*PDstandardNth11gt23) - 2*gtu21*PDstandardNth12gt23 - 2*gtu31*PDstandardNth13gt23 - 
          gtu22*PDstandardNth22gt23 - 2*gtu32*PDstandardNth23gt23 + gt13L*PDstandardNth2Xt1 + gt23L*PDstandardNth2Xt2 + 
          gt33L*PDstandardNth2Xt3 - gtu33*PDstandardNth33gt23 + gt12L*PDstandardNth3Xt1 + gt22L*PDstandardNth3Xt2 + 
          gt23L*PDstandardNth3Xt3 + (Gt113*gt12L + Gt213*gt22L + gt23L*Gt313)*Xt1L + 
          (Gt112*gt13L + Gt212*gt23L + Gt312*gt33L)*Xt1L + (Gt123*gt12L + Gt223*gt22L + gt23L*Gt323)*Xt2L + 
          (Gt122*gt13L + Gt222*gt23L + Gt322*gt33L)*Xt2L + (gt12L*Gt133 + gt22L*Gt233 + gt23L*Gt333)*Xt3L + 
          (Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)*Xt3L + 
          2*((Gt112*gt11L*Gt123 + Gt111*Gt123*gt12L + Gt111*Gt122*gt13L + Gt123*gt12L*Gt212 + Gt112*gt13L*Gt222 + 
                2*Gt112*gt12L*Gt223 + Gt123*Gt211*gt22L + 2*Gt212*Gt223*gt22L + Gt122*Gt211*gt23L + Gt212*Gt222*gt23L + 
                Gt123*gt23L*Gt311 + Gt123*gt13L*Gt312 + 2*Gt223*gt23L*Gt312 + Gt113*gt13L*Gt322 + Gt213*gt23L*Gt322 + 
                Gt113*gt12L*Gt323 + Gt112*gt13L*Gt323 + Gt213*gt22L*Gt323 + Gt212*gt23L*Gt323 + gt23L*Gt313*Gt323 + 
                Gt122*Gt311*gt33L + Gt222*Gt312*gt33L + Gt313*Gt322*gt33L + Gt312*Gt323*gt33L)*gtu21 + 
             (Gt112*gt11L*Gt133 + Gt111*gt12L*Gt133 + Gt111*Gt123*gt13L + gt12L*Gt133*Gt212 + Gt112*gt13L*Gt223 + 
                Gt133*Gt211*gt22L + 2*Gt112*gt12L*Gt233 + 2*Gt212*gt22L*Gt233 + Gt123*Gt211*gt23L + Gt212*Gt223*gt23L + 
                Gt133*gt23L*Gt311 + Gt133*gt13L*Gt312 + 2*Gt233*gt23L*Gt312 + Gt113*gt13L*Gt323 + Gt213*gt23L*Gt323 + 
                Gt113*gt12L*Gt333 + Gt112*gt13L*Gt333 + Gt213*gt22L*Gt333 + Gt212*gt23L*Gt333 + gt23L*Gt313*Gt333 + 
                Gt123*Gt311*gt33L + Gt223*Gt312*gt33L + Gt313*Gt323*gt33L + Gt312*Gt333*gt33L)*gtu31 + 
             gtu21*(Gt113*gt11L*Gt122 + Gt122*gt13L*Gt212 + 2*Gt122*gt12L*Gt213 + Gt113*gt12L*Gt222 + Gt113*Gt212*gt22L + 
                2*Gt213*Gt222*gt22L + Gt212*Gt222*gt23L + Gt123*gt13L*Gt312 + Gt113*gt23L*Gt312 + Gt223*gt23L*Gt312 + 
                Gt123*gt12L*Gt313 + Gt122*gt13L*Gt313 + Gt223*gt22L*Gt313 + Gt222*gt23L*Gt313 + Gt113*gt13L*Gt322 + 
                2*Gt213*gt23L*Gt322 + gt23L*Gt313*Gt323 + Gt212*Gt322*gt33L + Gt313*Gt322*gt33L + Gt312*Gt323*gt33L + 
                Gt112*(Gt113*gt12L + Gt212*gt23L + Gt312*gt33L) + gt13L*SQR(Gt112))) + 
          2*gtu31*(2*Gt213*Gt223*gt22L + Gt112*Gt213*gt23L + Gt212*Gt223*gt23L + Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + 
             gt12L*Gt133*Gt313 + gt22L*Gt233*Gt313 + Gt223*gt23L*Gt313 + Gt123*(2*gt12L*Gt213 + gt13L*(Gt212 + Gt313)) + 
             2*Gt213*gt23L*Gt323 + Gt113*(gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + 
                gt13L*Gt323) + gt23L*Gt313*Gt333 + Gt112*Gt313*gt33L + Gt212*Gt323*gt33L + Gt313*Gt323*gt33L + 
             Gt312*Gt333*gt33L + gt12L*SQR(Gt113)) + 
          2*gtu11*(Gt112*Gt113*gt11L + Gt111*Gt113*gt12L + Gt111*Gt112*gt13L + Gt113*gt12L*Gt212 + Gt112*gt13L*Gt212 + 
             2*Gt112*gt12L*Gt213 + Gt113*Gt211*gt22L + 2*Gt212*Gt213*gt22L + Gt112*Gt211*gt23L + Gt113*gt23L*Gt311 + 
             2*Gt113*gt13L*Gt312 + 3*Gt213*gt23L*Gt312 + Gt113*gt12L*Gt313 + Gt112*gt13L*Gt313 + Gt213*gt22L*Gt313 + 
             Gt212*gt23L*Gt313 + Gt112*Gt311*gt33L + Gt212*Gt312*gt33L + 2*Gt312*Gt313*gt33L + gt23L*SQR(Gt212) + 
             gt23L*SQR(Gt313)) + 2*gtu22*(gt11L*Gt122*Gt123 + Gt112*Gt123*gt12L + Gt112*Gt122*gt13L + Gt123*gt12L*Gt222 + 
             Gt122*gt13L*Gt222 + 2*Gt122*gt12L*Gt223 + Gt123*Gt212*gt22L + 2*Gt222*Gt223*gt22L + Gt122*Gt212*gt23L + 
             Gt123*gt23L*Gt312 + 2*Gt123*gt13L*Gt322 + 3*Gt223*gt23L*Gt322 + Gt123*gt12L*Gt323 + Gt122*gt13L*Gt323 + 
             Gt223*gt22L*Gt323 + Gt222*gt23L*Gt323 + Gt122*Gt312*gt33L + Gt222*Gt322*gt33L + 2*Gt322*Gt323*gt33L + 
             gt23L*SQR(Gt222) + gt23L*SQR(Gt323)) + 2*gtu32*
           (gt11L*Gt122*Gt133 + Gt112*gt12L*Gt133 + Gt112*Gt123*gt13L + gt12L*Gt133*Gt222 + Gt122*gt13L*Gt223 + 
             Gt133*Gt212*gt22L + 2*Gt122*gt12L*Gt233 + 2*Gt222*gt22L*Gt233 + Gt123*Gt212*gt23L + Gt222*Gt223*gt23L + 
             Gt133*gt23L*Gt312 + Gt133*gt13L*Gt322 + 2*Gt233*gt23L*Gt322 + Gt123*gt13L*Gt323 + Gt223*gt23L*Gt323 + 
             Gt123*gt12L*Gt333 + Gt122*gt13L*Gt333 + Gt223*gt22L*Gt333 + Gt222*gt23L*Gt333 + gt23L*Gt323*Gt333 + 
             Gt123*Gt312*gt33L + Gt223*Gt322*gt33L + Gt322*Gt333*gt33L + gt33L*SQR(Gt323)) + 
          2*gtu32*(Gt113*Gt123*gt12L + Gt113*Gt122*gt13L + Gt123*gt13L*Gt222 + 3*Gt123*gt12L*Gt223 + Gt123*Gt213*gt22L + 
             Gt122*Gt213*gt23L + Gt222*Gt223*gt23L + Gt123*gt23L*Gt313 + Gt133*gt13L*Gt322 + Gt233*gt23L*Gt322 + 
             gt12L*Gt133*Gt323 + 2*Gt123*gt13L*Gt323 + gt22L*Gt233*Gt323 + 3*Gt223*gt23L*Gt323 + gt23L*Gt323*Gt333 + 
             Gt122*Gt313*gt33L + Gt222*Gt323*gt33L + Gt322*Gt333*gt33L + gt11L*SQR(Gt123) + 2*gt22L*SQR(Gt223) + 
             gt33L*SQR(Gt323)) + 2*gtu33*(gt11L*Gt123*Gt133 + Gt113*gt12L*Gt133 + Gt113*Gt123*gt13L + gt12L*Gt133*Gt223 + 
             Gt123*gt13L*Gt223 + Gt133*Gt213*gt22L + 2*Gt123*gt12L*Gt233 + 2*Gt223*gt22L*Gt233 + Gt123*Gt213*gt23L + 
             Gt133*gt23L*Gt313 + 2*Gt133*gt13L*Gt323 + 3*Gt233*gt23L*Gt323 + gt12L*Gt133*Gt333 + Gt123*gt13L*Gt333 + 
             gt22L*Gt233*Gt333 + Gt223*gt23L*Gt333 + Gt123*Gt313*gt33L + Gt223*Gt323*gt33L + 2*Gt323*Gt333*gt33L + 
             gt23L*SQR(Gt223) + gt23L*SQR(Gt333)));
    
    CCTK_REAL const Rt33  =  4*(Gt123*gt13L*Gt323*gtu22 + Gt223*gt23L*Gt323*gtu22 + Gt133*gt13L*Gt313*gtu31 + Gt233*gt23L*Gt313*gtu31 + 
           Gt113*gt13L*Gt333*gtu31 + Gt133*gt13L*Gt323*gtu32 + Gt233*gt23L*Gt323*gtu32 + Gt123*gt13L*Gt333*gtu32 + 
           Gt133*gt13L*Gt333*gtu33) + gtu21*(2*Gt212*Gt223*gt23L + 4*Gt123*gt13L*Gt313 + 4*Gt223*gt23L*Gt313 + 
           4*Gt113*gt13L*Gt323 + 4*Gt213*gt23L*Gt323 + 2*Gt123*Gt311*gt33L - PDstandardNth12gt33) + 
        gtu31*(4*Gt213*gt23L*Gt333 + 2*Gt233*Gt312*gt33L + 6*Gt313*Gt333*gt33L - PDstandardNth13gt33) - 
        gtu22*khalf*PDstandardNth22gt33 + gtu32*(4*Gt223*gt23L*Gt333 + 2*Gt123*Gt313*gt33L + 6*Gt323*Gt333*gt33L - 
           PDstandardNth23gt33) - gtu33*khalf*PDstandardNth33gt33 + gt13L*PDstandardNth3Xt1 + gt23L*PDstandardNth3Xt2 + 
        gt33L*(2*Gt213*Gt322*gtu21 + 6*Gt313*Gt323*gtu21 + 2*Gt123*Gt312*gtu22 + 2*Gt133*Gt311*gtu31 + 
           2*Gt133*Gt312*gtu32 + 2*Gt133*Gt313*gtu33 + PDstandardNth3Xt3) + Gt113*gt13L*Xt1L + Gt213*gt23L*Xt1L + 
        Gt313*gt33L*Xt1L + Gt123*gt13L*Xt2L + Gt223*gt23L*Xt2L + Gt323*gt33L*Xt2L + Gt133*gt13L*Xt3L + Gt333*gt33L*Xt3L + 
        Gt233*(4*gt23L*Gt333*gtu33 + 2*Gt323*gt33L*gtu33 + gt23L*Xt3L) + 
        gtu11*(2*Gt212*Gt213*gt23L + 4*Gt113*gt13L*Gt313 + 4*Gt213*gt23L*Gt313 + 2*Gt113*Gt311*gt33L + 
           2*Gt213*Gt312*gt33L - khalf*PDstandardNth11gt33 + gt11L*SQR(Gt113)) + 
        2*(Gt111*Gt113*gt13L*gtu11 + Gt113*gt12L*Gt213*gtu11 + Gt112*gt13L*Gt213*gtu11 + Gt113*Gt211*gt23L*gtu11 + 
           Gt113*gt11L*Gt123*gtu21 + Gt112*Gt113*gt13L*gtu21 + Gt111*Gt123*gt13L*gtu21 + Gt123*gt12L*Gt213*gtu21 + 
           Gt122*gt13L*Gt213*gtu21 + Gt113*gt12L*Gt223*gtu21 + Gt112*gt13L*Gt223*gtu21 + Gt213*Gt223*gt22L*gtu21 + 
           Gt123*Gt211*gt23L*gtu21 + Gt113*Gt212*gt23L*gtu21 + Gt213*Gt222*gt23L*gtu21 + Gt113*Gt312*gt33L*gtu21 + 
           Gt223*Gt312*gt33L*gtu21 + Gt112*Gt123*gt13L*gtu22 + Gt123*gt12L*Gt223*gtu22 + Gt122*gt13L*Gt223*gtu22 + 
           Gt123*Gt212*gt23L*gtu22 + Gt222*Gt223*gt23L*gtu22 + Gt223*Gt322*gt33L*gtu22 + Gt113*gt11L*Gt133*gtu31 + 
           Gt111*Gt133*gt13L*gtu31 + gt12L*Gt133*Gt213*gtu31 + Gt123*gt13L*Gt213*gtu31 + Gt113*gt12L*Gt233*gtu31 + 
           Gt112*gt13L*Gt233*gtu31 + Gt213*gt22L*Gt233*gtu31 + Gt133*Gt211*gt23L*gtu31 + Gt113*Gt213*gt23L*gtu31 + 
           Gt213*Gt223*gt23L*gtu31 + Gt212*Gt233*gt23L*gtu31 + Gt113*Gt313*gt33L*gtu31 + Gt213*Gt323*gt33L*gtu31 + 
           gt11L*Gt123*Gt133*gtu32 + Gt113*Gt123*gt13L*gtu32 + Gt112*Gt133*gt13L*gtu32 + gt12L*Gt133*Gt223*gtu32 + 
           Gt123*gt13L*Gt223*gtu32 + Gt123*gt12L*Gt233*gtu32 + Gt122*gt13L*Gt233*gtu32 + Gt223*gt22L*Gt233*gtu32 + 
           Gt133*Gt212*gt23L*gtu32 + Gt123*Gt213*gt23L*gtu32 + Gt222*Gt233*gt23L*gtu32 + Gt233*Gt322*gt33L*gtu32 + 
           Gt223*Gt323*gt33L*gtu32 + Gt113*Gt133*gt13L*gtu33 + gt12L*Gt133*Gt233*gtu33 + Gt123*gt13L*Gt233*gtu33 + 
           Gt133*Gt213*gt23L*gtu33 + Gt223*Gt233*gt23L*gtu33 + gt13L*gtu31*SQR(Gt113)) + gt11L*gtu22*SQR(Gt123) + 
        gt11L*gtu33*SQR(Gt133) + gt22L*gtu11*SQR(Gt213) + gt22L*gtu22*SQR(Gt223) + 2*gt23L*gtu32*SQR(Gt223) + 
        gt22L*gtu33*SQR(Gt233) + 3*gt33L*gtu11*SQR(Gt313) + 3*gt33L*gtu22*SQR(Gt323) + 3*gt33L*gtu33*SQR(Gt333);
    
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
    
    CCTK_REAL const e4phi  =  IfThen(conformalMethod,pow(phiL,-2),exp(4*phiL));
    
    CCTK_REAL const em4phi  =  INV(e4phi);
    
    CCTK_REAL const gu11  =  em4phi*gtu11;
    
    CCTK_REAL const gu21  =  em4phi*gtu21;
    
    CCTK_REAL const gu31  =  em4phi*gtu31;
    
    CCTK_REAL const gu22  =  em4phi*gtu22;
    
    CCTK_REAL const gu32  =  em4phi*gtu32;
    
    CCTK_REAL const gu33  =  em4phi*gtu33;
    
    CCTK_REAL const R11  =  Rphi11 + Rt11;
    
    CCTK_REAL const R12  =  Rphi12 + Rt12;
    
    CCTK_REAL const R13  =  Rphi13 + Rt13;
    
    CCTK_REAL const R22  =  Rphi22 + Rt22;
    
    CCTK_REAL const R23  =  Rphi23 + Rt23;
    
    CCTK_REAL const R33  =  Rphi33 + Rt33;
    
    CCTK_REAL const trR  =  gu11*R11 + gu22*R22 + 2*(gu21*R12 + gu31*R13 + gu32*R23) + gu33*R33;
    
    CCTK_REAL const Atm11  =  At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    CCTK_REAL const Atm21  =  At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    CCTK_REAL const Atm31  =  At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    CCTK_REAL const Atm12  =  At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    CCTK_REAL const Atm22  =  At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    CCTK_REAL const Atm32  =  At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    CCTK_REAL const Atm13  =  At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    CCTK_REAL const Atm23  =  At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    CCTK_REAL const Atm33  =  At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    CCTK_REAL const T00  =  eTttL;
    
    CCTK_REAL const T01  =  eTtxL;
    
    CCTK_REAL const T02  =  eTtyL;
    
    CCTK_REAL const T03  =  eTtzL;
    
    CCTK_REAL const T11  =  eTxxL;
    
    CCTK_REAL const T12  =  eTxyL;
    
    CCTK_REAL const T13  =  eTxzL;
    
    CCTK_REAL const T22  =  eTyyL;
    
    CCTK_REAL const T23  =  eTyzL;
    
    CCTK_REAL const T33  =  eTzzL;
    
    CCTK_REAL const rho  =  pow(alphaL,-2)*(T00 - 2*(beta2L*T02 + beta3L*T03) + 
          2*(beta1L*(-T01 + beta2L*T12 + beta3L*T13) + beta2L*beta3L*T23) + T11*SQR(beta1L) + T22*SQR(beta2L) + 
          T33*SQR(beta3L));
    
    CCTK_REAL const S1  =  (-T01 + beta1L*T11 + beta2L*T12 + beta3L*T13)*INV(alphaL);
    
    CCTK_REAL const S2  =  (-T02 + beta1L*T12 + beta2L*T22 + beta3L*T23)*INV(alphaL);
    
    CCTK_REAL const S3  =  (-T03 + beta1L*T13 + beta2L*T23 + beta3L*T33)*INV(alphaL);
    
    CCTK_REAL const HL  =  -2.*(Atm12*Atm21 + Atm13*Atm31 + Atm23*Atm32) - 50.26548245743669181540229413247204614715*rho + trR - 
        1.*(SQR(Atm11) + SQR(Atm22) + SQR(Atm33)) + 0.6666666666666666666666666666666666666667*SQR(trKL);
    
    CCTK_REAL const M1L  =  -2.*((At12L*Gt211 + At13L*Gt311)*gtu11 + At11L*Gt123*gtu32) + 
        At12L*((6.*cdphi1 - 1.*Gt111)*gtu21 - 3.*Gt213*gtu31 + 6.*(cdphi2*gtu22 + cdphi3*gtu32) - 1.*Gt233*gtu33) - 
        1.*((At22L*Gt212 + At12L*(Gt112 + Gt222) + At23L*Gt312 + At13L*Gt322)*gtu22 + 
           At13L*(Gt111*gtu31 + Gt112*gtu32 + Gt113*gtu33)) + 
        At11L*((6.*cdphi1 - 2.*Gt111)*gtu11 + 6.*(cdphi2*gtu21 + cdphi3*gtu31) - 1.*(Gt122*gtu22 + Gt133*gtu33)) + 
        gtu11*PDstandardNth1At11 - 0.6666666666666666666666666666666666666667*PDstandardNth1trK + 
        gtu21*(-1.*(At22L*Gt211 + At23L*Gt311) - 3.*(At11L*Gt112 + At12L*Gt212 + At13L*Gt312) + PDstandardNth1At12 + 
           PDstandardNth2At11) + gtu22*PDstandardNth2At12 + 
        gtu31*(-3.*At11L*Gt113 - 1.*(At23L*Gt211 + At33L*Gt311) + At13L*(6.*cdphi1 - 3.*Gt313) + PDstandardNth1At13 + 
           PDstandardNth3At11) + gtu32*(-2.*At12L*Gt223 - 
           1.*(At12L*Gt113 + At22L*Gt213 + At33L*Gt312 + At23L*(Gt212 + Gt313)) + At13L*(6.*cdphi2 - 2.*Gt323) + 
           PDstandardNth2At13 + PDstandardNth3At12) + 
        gtu33*(-1.*(At23L*Gt213 + At33L*Gt313) + At13L*(6.*cdphi3 - 1.*Gt333) + PDstandardNth3At13) - 
        25.13274122871834590770114706623602307358*S1;
    
    CCTK_REAL const M2L  =  At12L*((6.*cdphi1 - 1.*Gt111)*gtu11 - 2.*Gt122*gtu22 + 6.*(cdphi2*gtu21 + cdphi3*gtu31) - 3.*Gt123*gtu32 - 
           1.*Gt133*gtu33) + At22L*((6.*cdphi2 - 2.*Gt222)*gtu22 - 2.*Gt213*gtu31 + 6.*cdphi3*gtu32 - 1.*Gt233*gtu33) + 
        At23L*(-2.*Gt322*gtu22 - 1.*Gt333*gtu33 + 6.*(cdphi1*gtu31 + cdphi2*gtu32 + cdphi3*gtu33)) - 
        1.*((At11L*Gt112 + At22L*Gt211 + At12L*Gt212 + At23L*Gt311 + At13L*Gt312)*gtu11 + 
           Gt122*(At11L*gtu21 + At13L*gtu32) + (At23L*Gt223 + At33L*Gt323)*gtu33 + At13L*(Gt112*gtu31 + Gt123*gtu33)) + 
        gtu11*PDstandardNth1At12 + gtu21*(At22L*(6.*cdphi1 - 3.*Gt212) + At12L*(-3.*Gt112 - 1.*Gt222) - 3.*At23L*Gt312 - 
           1.*At13L*Gt322 + PDstandardNth1At22 + PDstandardNth2At12) + gtu22*PDstandardNth2At22 - 
        0.6666666666666666666666666666666666666667*PDstandardNth2trK + 
        gtu31*(At12L*(-2.*Gt113 - 1.*Gt223) - 2.*At23L*Gt313 - 1.*(At11L*Gt123 + At23L*Gt212 + At33L*Gt312 + At13L*Gt323) + 
           PDstandardNth1At23 + PDstandardNth3At12) + 
        gtu32*(-1.*(At23L*Gt222 + At33L*Gt322) - 3.*(At22L*Gt223 + At23L*Gt323) + PDstandardNth2At23 + 
           PDstandardNth3At22) + gtu33*PDstandardNth3At23 - 25.13274122871834590770114706623602307358*S2;
    
    CCTK_REAL const M3L  =  -1.*((At11L*Gt113 + At23L*Gt211 + At12L*Gt213 + At33L*Gt311)*gtu11 + 
           (At22L*Gt223 + At33L*Gt322 + At23L*Gt323)*gtu22 + At11L*Gt133*gtu31 + 
           At12L*(Gt113*gtu21 + Gt123*gtu22 + Gt133*gtu32)) + 
        At13L*((6.*cdphi1 - 1.*(Gt111 + Gt313))*gtu11 - 1.*Gt122*gtu22 + 6.*(cdphi2*gtu21 + cdphi3*gtu31) - 
           3.*Gt123*gtu32 - 2.*Gt133*gtu33) + At23L*((6.*cdphi2 - 1.*Gt222)*gtu22 - 3.*Gt213*gtu31 + 6.*cdphi3*gtu32 - 
           2.*Gt233*gtu33) + gtu11*PDstandardNth1At13 + 
        gtu21*(-2.*(At13L*Gt112 + At33L*Gt312) + At23L*(6.*cdphi1 - 2.*Gt212 - 1.*Gt313) - 
           1.*(At11L*Gt123 + At22L*Gt213 + At12L*Gt223 + At13L*Gt323) + PDstandardNth1At23 + PDstandardNth2At13) + 
        gtu22*PDstandardNth2At23 + gtu31*(-1.*At12L*Gt233 + At33L*(6.*cdphi1 - 3.*Gt313) + At13L*(-3.*Gt113 - 1.*Gt333) + 
           PDstandardNth1At33 + PDstandardNth3At13) + 
        gtu32*(-1.*At22L*Gt233 + At33L*(6.*cdphi2 - 3.*Gt323) + At23L*(-3.*Gt223 - 1.*Gt333) + PDstandardNth2At33 + 
           PDstandardNth3At23) + gtu33*(At33L*(6.*cdphi3 - 2.*Gt333) + PDstandardNth3At33) - 
        0.6666666666666666666666666666666666666667*PDstandardNth3trK - 25.13274122871834590770114706623602307358*S3;
    
    CCTK_REAL const cSL  =  Log(detgt);
    
    CCTK_REAL const cXt1L  =  Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33 - Xt1L;
    
    CCTK_REAL const cXt2L  =  Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33 - Xt2L;
    
    CCTK_REAL const cXt3L  =  Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33 - Xt3L;
    
    CCTK_REAL const cAL  =  At11L*gtu11 + At22L*gtu22 + 2*(At12L*gtu21 + At13L*gtu31 + At23L*gtu32) + At33L*gtu33;
    
    
    /* Copy local copies back to grid functions */
    cA[index] = cAL;
    cS[index] = cSL;
    cXt1[index] = cXt1L;
    cXt2[index] = cXt2L;
    cXt3[index] = cXt3L;
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSN_constraints);
}

void ML_BSSN_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_constraints_Body);
}
