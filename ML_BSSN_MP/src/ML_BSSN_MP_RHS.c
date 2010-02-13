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

void ML_BSSN_MP_RHS_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_RHS_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_RHS_calc_every != ML_BSSN_MP_RHS_calc_offset)
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
  LC_LOOP3 (ML_BSSN_MP_RHS,
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
    // CCTK_REAL Atu11 = INITVALUE, Atu21 = INITVALUE, Atu22 = INITVALUE, Atu31 = INITVALUE, Atu32 = INITVALUE, Atu33 = INITVALUE;
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
    // CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    // CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    // CCTK_REAL R11 = INITVALUE, R12 = INITVALUE, R13 = INITVALUE, R22 = INITVALUE, R23 = INITVALUE, R33 = INITVALUE;
    // CCTK_REAL rho = INITVALUE;
    // CCTK_REAL Rphi11 = INITVALUE, Rphi12 = INITVALUE, Rphi13 = INITVALUE, Rphi22 = INITVALUE, Rphi23 = INITVALUE, Rphi33 = INITVALUE;
    // CCTK_REAL Rt11 = INITVALUE, Rt12 = INITVALUE, Rt13 = INITVALUE, Rt22 = INITVALUE, Rt23 = INITVALUE, Rt33 = INITVALUE;
    // CCTK_REAL S1 = INITVALUE, S2 = INITVALUE, S3 = INITVALUE;
    // CCTK_REAL T00 = INITVALUE, T01 = INITVALUE, T02 = INITVALUE, T03 = INITVALUE, T11 = INITVALUE, T12 = INITVALUE;
    // CCTK_REAL T13 = INITVALUE, T22 = INITVALUE, T23 = INITVALUE, T33 = INITVALUE;
    // CCTK_REAL trAts = INITVALUE;
    // CCTK_REAL trS = INITVALUE;
    // CCTK_REAL Xtn1 = INITVALUE, Xtn2 = INITVALUE, Xtn3 = INITVALUE;
    
    /* Declare local copies of grid functions */
    // CCTK_REAL AL = INITVALUE;
    // CCTK_REAL alphaL = INITVALUE, alpharhsL = INITVALUE;
    // CCTK_REAL ArhsL = INITVALUE;
    // CCTK_REAL At11L = INITVALUE, At11rhsL = INITVALUE, At12L = INITVALUE, At12rhsL = INITVALUE, At13L = INITVALUE, At13rhsL = INITVALUE;
    // CCTK_REAL At22L = INITVALUE, At22rhsL = INITVALUE, At23L = INITVALUE, At23rhsL = INITVALUE, At33L = INITVALUE, At33rhsL = INITVALUE;
    // CCTK_REAL B1L = INITVALUE, B1rhsL = INITVALUE, B2L = INITVALUE, B2rhsL = INITVALUE, B3L = INITVALUE, B3rhsL = INITVALUE;
    // CCTK_REAL beta1L = INITVALUE, beta1rhsL = INITVALUE, beta2L = INITVALUE, beta2rhsL = INITVALUE, beta3L = INITVALUE, beta3rhsL = INITVALUE;
    // CCTK_REAL dJ111L = INITVALUE, dJ112L = INITVALUE, dJ113L = INITVALUE, dJ122L = INITVALUE, dJ123L = INITVALUE, dJ133L = INITVALUE;
    // CCTK_REAL dJ211L = INITVALUE, dJ212L = INITVALUE, dJ213L = INITVALUE, dJ222L = INITVALUE, dJ223L = INITVALUE, dJ233L = INITVALUE;
    // CCTK_REAL dJ311L = INITVALUE, dJ312L = INITVALUE, dJ313L = INITVALUE, dJ322L = INITVALUE, dJ323L = INITVALUE, dJ333L = INITVALUE;
    // CCTK_REAL etaL = INITVALUE;
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
    // CCTK_REAL gt11L = INITVALUE, gt11rhsL = INITVALUE, gt12L = INITVALUE, gt12rhsL = INITVALUE, gt13L = INITVALUE, gt13rhsL = INITVALUE;
    // CCTK_REAL gt22L = INITVALUE, gt22rhsL = INITVALUE, gt23L = INITVALUE, gt23rhsL = INITVALUE, gt33L = INITVALUE, gt33rhsL = INITVALUE;
    // CCTK_REAL J11L = INITVALUE, J12L = INITVALUE, J13L = INITVALUE, J21L = INITVALUE, J22L = INITVALUE, J23L = INITVALUE;
    // CCTK_REAL J31L = INITVALUE, J32L = INITVALUE, J33L = INITVALUE;
    // CCTK_REAL phiL = INITVALUE, phirhsL = INITVALUE;
    // CCTK_REAL trKL = INITVALUE, trKrhsL = INITVALUE;
    // CCTK_REAL Xt1L = INITVALUE, Xt1rhsL = INITVALUE, Xt2L = INITVALUE, Xt2rhsL = INITVALUE, Xt3L = INITVALUE, Xt3rhsL = INITVALUE;
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
    CCTK_REAL const AL = A[index];
    CCTK_REAL const alphaL = alpha[index];
    CCTK_REAL const At11L = At11[index];
    CCTK_REAL const At12L = At12[index];
    CCTK_REAL const At13L = At13[index];
    CCTK_REAL const At22L = At22[index];
    CCTK_REAL const At23L = At23[index];
    CCTK_REAL const At33L = At33[index];
    CCTK_REAL const B1L = B1[index];
    CCTK_REAL const B2L = B2[index];
    CCTK_REAL const B3L = B3[index];
    CCTK_REAL const beta1L = beta1[index];
    CCTK_REAL const beta2L = beta2[index];
    CCTK_REAL const beta3L = beta3[index];
    CCTK_REAL const dJ111L = dJ111[index];
    CCTK_REAL const dJ112L = dJ112[index];
    CCTK_REAL const dJ113L = dJ113[index];
    CCTK_REAL const dJ122L = dJ122[index];
    CCTK_REAL const dJ123L = dJ123[index];
    CCTK_REAL const dJ133L = dJ133[index];
    CCTK_REAL const dJ211L = dJ211[index];
    CCTK_REAL const dJ212L = dJ212[index];
    CCTK_REAL const dJ213L = dJ213[index];
    CCTK_REAL const dJ222L = dJ222[index];
    CCTK_REAL const dJ223L = dJ223[index];
    CCTK_REAL const dJ233L = dJ233[index];
    CCTK_REAL const dJ311L = dJ311[index];
    CCTK_REAL const dJ312L = dJ312[index];
    CCTK_REAL const dJ313L = dJ313[index];
    CCTK_REAL const dJ322L = dJ322[index];
    CCTK_REAL const dJ323L = dJ323[index];
    CCTK_REAL const dJ333L = dJ333[index];
    CCTK_REAL const etaL = eta[index];
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
    CCTK_REAL const J11L = J11[index];
    CCTK_REAL const J12L = J12[index];
    CCTK_REAL const J13L = J13[index];
    CCTK_REAL const J21L = J21[index];
    CCTK_REAL const J22L = J22[index];
    CCTK_REAL const J23L = J23[index];
    CCTK_REAL const J31L = J31[index];
    CCTK_REAL const J32L = J32[index];
    CCTK_REAL const J33L = J33[index];
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
    
    CCTK_REAL const Gt111  =  khalf*((gtu11*J11L - gtu21*J12L - gtu31*J13L)*PDstandardNth1gt11 + 
          (gtu11*J21L - gtu21*J22L - gtu31*J23L)*PDstandardNth2gt11 + 
          (gtu11*J31L - gtu21*J32L - gtu31*J33L)*PDstandardNth3gt11 + 
          2*(J11L*(gtu21*PDstandardNth1gt12 + gtu31*PDstandardNth1gt13) + 
             J21L*(gtu21*PDstandardNth2gt12 + gtu31*PDstandardNth2gt13) + 
             J31L*(gtu21*PDstandardNth3gt12 + gtu31*PDstandardNth3gt13)));
    
    CCTK_REAL const Gt211  =  khalf*((gtu21*J11L - gtu22*J12L - gtu32*J13L)*PDstandardNth1gt11 + 
          (gtu21*J21L - gtu22*J22L - gtu32*J23L)*PDstandardNth2gt11 + 
          (gtu21*J31L - gtu22*J32L - gtu32*J33L)*PDstandardNth3gt11 + 
          2*(J11L*(gtu22*PDstandardNth1gt12 + gtu32*PDstandardNth1gt13) + 
             J21L*(gtu22*PDstandardNth2gt12 + gtu32*PDstandardNth2gt13) + 
             J31L*(gtu22*PDstandardNth3gt12 + gtu32*PDstandardNth3gt13)));
    
    CCTK_REAL const Gt311  =  khalf*((gtu31*J11L - gtu32*J12L - gtu33*J13L)*PDstandardNth1gt11 + 
          (gtu31*J21L - gtu32*J22L - gtu33*J23L)*PDstandardNth2gt11 + 
          (gtu31*J31L - gtu32*J32L - gtu33*J33L)*PDstandardNth3gt11 + 
          2*(J11L*(gtu32*PDstandardNth1gt12 + gtu33*PDstandardNth1gt13) + 
             J21L*(gtu32*PDstandardNth2gt12 + gtu33*PDstandardNth2gt13) + 
             J31L*(gtu32*PDstandardNth3gt12 + gtu33*PDstandardNth3gt13)));
    
    CCTK_REAL const Gt112  =  khalf*(gtu11*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
          gtu21*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22) + 
          gtu31*(-(J13L*PDstandardNth1gt12) + J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23));
    
    CCTK_REAL const Gt212  =  khalf*(gtu21*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
          gtu22*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22) + 
          gtu32*(-(J13L*PDstandardNth1gt12) + J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23));
    
    CCTK_REAL const Gt312  =  khalf*(gtu31*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
          gtu32*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22) + 
          gtu33*(-(J13L*PDstandardNth1gt12) + J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23));
    
    CCTK_REAL const Gt113  =  khalf*(gtu11*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
          gtu21*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23) + gtu31*(J11L*PDstandardNth1gt33 + J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL const Gt213  =  khalf*(gtu21*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
          gtu22*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23) + gtu32*(J11L*PDstandardNth1gt33 + J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL const Gt313  =  khalf*(gtu31*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
          gtu32*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23) + gtu33*(J11L*PDstandardNth1gt33 + J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL const Gt122  =  khalf*(gtu11*(-(J11L*PDstandardNth1gt22) + 2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
             J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - J31L*PDstandardNth3gt22) + 
          gtu21*(J12L*PDstandardNth1gt22 + J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
          gtu31*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22 - 
             2*(J12L*PDstandardNth1gt23 + J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL const Gt222  =  khalf*(gtu21*(-(J11L*PDstandardNth1gt22) + 2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
             J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - J31L*PDstandardNth3gt22) + 
          gtu22*(J12L*PDstandardNth1gt22 + J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
          gtu32*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22 - 
             2*(J12L*PDstandardNth1gt23 + J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL const Gt322  =  khalf*(gtu31*(-(J11L*PDstandardNth1gt22) + 2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
             J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - J31L*PDstandardNth3gt22) + 
          gtu32*(J12L*PDstandardNth1gt22 + J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
          gtu33*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22 - 
             2*(J12L*PDstandardNth1gt23 + J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL const Gt123  =  khalf*(gtu21*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
          gtu11*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt23) + gtu31*(J12L*PDstandardNth1gt33 + J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL const Gt223  =  khalf*(gtu22*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
          gtu21*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt23) + gtu32*(J12L*PDstandardNth1gt33 + J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL const Gt323  =  khalf*(gtu32*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
          gtu31*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt23) + gtu33*(J12L*PDstandardNth1gt33 + J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL const Gt133  =  khalf*(gtu11*(-(J11L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
             J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - J31L*PDstandardNth3gt33) + 
          gtu21*(-(J12L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
             J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - J32L*PDstandardNth3gt33) + 
          gtu31*(J13L*PDstandardNth1gt33 + J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL const Gt233  =  khalf*(gtu21*(-(J11L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
             J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - J31L*PDstandardNth3gt33) + 
          gtu22*(-(J12L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
             J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - J32L*PDstandardNth3gt33) + 
          gtu32*(J13L*PDstandardNth1gt33 + J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL const Gt333  =  khalf*(gtu31*(-(J11L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
             J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - J31L*PDstandardNth3gt33) + 
          gtu32*(-(J12L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
             J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - J32L*PDstandardNth3gt33) + 
          gtu33*(J13L*PDstandardNth1gt33 + J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL const Xtn1  =  Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    CCTK_REAL const Xtn2  =  Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    CCTK_REAL const Xtn3  =  Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    CCTK_REAL const Rt11  =  (Gt113*(gt13L*Gt312 + 3*(gt12L*Gt212 + gt13L*Gt312)) + 
           gt12L*(Gt213*(4*Gt112 + 2*Gt222) + Gt212*(Gt113 + 2*Gt223) + 2*(Gt233*Gt312 + Gt223*Gt313)) + 
           gt11L*(6*Gt112*Gt113 + 2*(Gt122*Gt213 + Gt133*Gt312 + Gt123*(Gt212 + Gt313))) + 
           gt13L*(2*Gt213*Gt322 + Gt313*(4*Gt112 + 2*Gt323)) + 
           2*(Gt213*(Gt212*gt22L + gt23L*Gt312) + Gt212*(gt23L*Gt313 + gt13L*Gt323) + Gt312*(gt13L*Gt333 + Gt313*gt33L)))*
         gtu32 + J11L*(gt11L*PDstandardNth1Xt1 + gt12L*PDstandardNth1Xt2 + gt13L*PDstandardNth1Xt3) + 
        J21L*(gt11L*PDstandardNth2Xt1 + gt12L*PDstandardNth2Xt2 + gt13L*PDstandardNth2Xt3) + 
        J31L*(gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + gt13L*PDstandardNth3Xt3) + 
        (Gt111*gt11L + gt12L*Gt211 + gt13L*Gt311)*Xtn1 + (Gt112*gt11L + gt12L*Gt212 + gt13L*Gt312)*Xtn2 + 
        (Gt113*gt11L + gt12L*Gt213 + gt13L*Gt313)*Xtn3 + 
        gtu21*(Gt112*(gt13L*Gt311 + 3*(gt12L*Gt211 + gt13L*Gt311)) + 
           gt11L*(Gt112*(6*Gt111 + 2*Gt212) + 2*(Gt122*Gt211 + Gt123*Gt311 + Gt113*Gt312)) + 
           2*(Gt212*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + gt13L*(Gt211*Gt322 + Gt311*Gt323)) + 
           Gt312*(gt13L*(4*Gt111 + 2*Gt313) + 2*(Gt211*gt23L + Gt311*gt33L)) + 
           gt12L*(4*Gt111*Gt212 + Gt211*(Gt112 + 2*Gt222) + 2*(Gt223*Gt311 + Gt213*Gt312 + SQR(Gt212)))) + 
        gtu11*(4*Gt111*(gt12L*Gt211 + gt13L*Gt311) + 
           2*(gt12L*(Gt211*Gt212 + Gt213*Gt311) + Gt211*(gt23L*Gt311 + gt13L*Gt312) + gt13L*Gt311*Gt313) + 
           gt11L*(2*(Gt112*Gt211 + Gt113*Gt311) + 3*SQR(Gt111)) + gt22L*SQR(Gt211) + gt33L*SQR(Gt311)) + 
        gtu22*(4*Gt112*(gt12L*Gt212 + gt13L*Gt312) + 
           2*(gt12L*(Gt212*Gt222 + Gt223*Gt312) + Gt212*(gt23L*Gt312 + gt13L*Gt322) + gt13L*Gt312*Gt323) + 
           gt11L*(2*(Gt122*Gt212 + Gt123*Gt312) + 3*SQR(Gt112)) + gt22L*SQR(Gt212) + gt33L*SQR(Gt312)) + 
        gtu33*(4*Gt113*(gt12L*Gt213 + gt13L*Gt313) + 
           2*(gt12L*(Gt213*Gt223 + Gt233*Gt313) + Gt213*(gt23L*Gt313 + gt13L*Gt323) + gt13L*Gt313*Gt333) + 
           gt11L*(2*(Gt123*Gt213 + Gt133*Gt313) + 3*SQR(Gt113)) + gt22L*SQR(Gt213) + gt33L*SQR(Gt313)) + 
        gtu31*(Gt113*(gt13L*Gt311 + 3*(gt12L*Gt211 + gt13L*Gt311)) + 
           gt11L*(2*(Gt123*Gt211 + Gt112*Gt213 + Gt133*Gt311) + Gt113*(6*Gt111 + 2*Gt313)) + 
           gt12L*(Gt211*(Gt113 + 2*Gt223) + 2*Gt233*Gt311 + Gt213*(4*Gt111 + 2*(Gt212 + Gt313))) + 
           2*(Gt213*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + gt13L*Gt311*Gt333 + Gt313*(Gt211*gt23L + Gt311*gt33L)) + 
           gt13L*(4*Gt111*Gt313 + 2*(Gt211*Gt323 + SQR(Gt313)))) + 
        khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt11 + J12L*J21L*PDstandardNth12gt11 + J11L*J22L*PDstandardNth12gt11 + 
                 J12L*J31L*PDstandardNth13gt11 + J11L*J32L*PDstandardNth13gt11 + dJ112L*PDstandardNth1gt11 + 
                 J21L*J22L*PDstandardNth22gt11 + J22L*J31L*PDstandardNth23gt11 + J21L*J32L*PDstandardNth23gt11 + 
                 dJ212L*PDstandardNth2gt11 + J31L*J32L*PDstandardNth33gt11 + dJ312L*PDstandardNth3gt11) + 
              gtu31*(J11L*J13L*PDstandardNth11gt11 + J13L*J21L*PDstandardNth12gt11 + J11L*J23L*PDstandardNth12gt11 + 
                 J13L*J31L*PDstandardNth13gt11 + J11L*J33L*PDstandardNth13gt11 + dJ113L*PDstandardNth1gt11 + 
                 J21L*J23L*PDstandardNth22gt11 + J23L*J31L*PDstandardNth23gt11 + J21L*J33L*PDstandardNth23gt11 + 
                 dJ213L*PDstandardNth2gt11 + J31L*J33L*PDstandardNth33gt11 + dJ313L*PDstandardNth3gt11) + 
              gtu32*(J12L*J13L*PDstandardNth11gt11 + J13L*J22L*PDstandardNth12gt11 + J12L*J23L*PDstandardNth12gt11 + 
                 J13L*J32L*PDstandardNth13gt11 + J12L*J33L*PDstandardNth13gt11 + dJ123L*PDstandardNth1gt11 + 
                 J22L*J23L*PDstandardNth22gt11 + J23L*J32L*PDstandardNth23gt11 + J22L*J33L*PDstandardNth23gt11 + 
                 dJ223L*PDstandardNth2gt11 + J32L*J33L*PDstandardNth33gt11 + dJ323L*PDstandardNth3gt11)) - 
           gtu11*(2*J11L*J21L*PDstandardNth12gt11 + 2*J11L*J31L*PDstandardNth13gt11 + dJ111L*PDstandardNth1gt11 + 
              2*J21L*J31L*PDstandardNth23gt11 + dJ211L*PDstandardNth2gt11 + dJ311L*PDstandardNth3gt11 + 
              PDstandardNth11gt11*SQR(J11L) + PDstandardNth22gt11*SQR(J21L) + PDstandardNth33gt11*SQR(J31L)) - 
           gtu22*(2*J12L*J22L*PDstandardNth12gt11 + 2*J12L*J32L*PDstandardNth13gt11 + dJ122L*PDstandardNth1gt11 + 
              2*J22L*J32L*PDstandardNth23gt11 + dJ222L*PDstandardNth2gt11 + dJ322L*PDstandardNth3gt11 + 
              PDstandardNth11gt11*SQR(J12L) + PDstandardNth22gt11*SQR(J22L) + PDstandardNth33gt11*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt11 + 2*J13L*J33L*PDstandardNth13gt11 + dJ133L*PDstandardNth1gt11 + 
              2*J23L*J33L*PDstandardNth23gt11 + dJ233L*PDstandardNth2gt11 + dJ333L*PDstandardNth3gt11 + 
              PDstandardNth11gt11*SQR(J13L) + PDstandardNth22gt11*SQR(J23L) + PDstandardNth33gt11*SQR(J33L)));
    
    CCTK_REAL const Rt12  =  khalf*((gt12L*J11L + gt11L*J12L)*PDstandardNth1Xt1 + (gt22L*J11L + gt12L*J12L)*PDstandardNth1Xt2 + 
          (gt23L*J11L + gt13L*J12L)*PDstandardNth1Xt3 + (gt12L*J21L + gt11L*J22L)*PDstandardNth2Xt1 + 
          (gt22L*J21L + gt12L*J22L)*PDstandardNth2Xt2 + (gt23L*J21L + gt13L*J22L)*PDstandardNth2Xt3 - 
          2*(gtu21*(J11L*J12L*PDstandardNth11gt12 + J12L*J21L*PDstandardNth12gt12 + J11L*J22L*PDstandardNth12gt12 + 
                J12L*J31L*PDstandardNth13gt12 + J11L*J32L*PDstandardNth13gt12 + dJ112L*PDstandardNth1gt12 + 
                J21L*J22L*PDstandardNth22gt12 + J22L*J31L*PDstandardNth23gt12 + J21L*J32L*PDstandardNth23gt12 + 
                dJ212L*PDstandardNth2gt12 + J31L*J32L*PDstandardNth33gt12 + dJ312L*PDstandardNth3gt12) + 
             gtu31*(J11L*J13L*PDstandardNth11gt12 + J13L*J21L*PDstandardNth12gt12 + J11L*J23L*PDstandardNth12gt12 + 
                J13L*J31L*PDstandardNth13gt12 + J11L*J33L*PDstandardNth13gt12 + dJ113L*PDstandardNth1gt12 + 
                J21L*J23L*PDstandardNth22gt12 + J23L*J31L*PDstandardNth23gt12 + J21L*J33L*PDstandardNth23gt12 + 
                dJ213L*PDstandardNth2gt12 + J31L*J33L*PDstandardNth33gt12 + dJ313L*PDstandardNth3gt12) + 
             gtu32*(J12L*J13L*PDstandardNth11gt12 + J13L*J22L*PDstandardNth12gt12 + J12L*J23L*PDstandardNth12gt12 + 
                J13L*J32L*PDstandardNth13gt12 + J12L*J33L*PDstandardNth13gt12 + dJ123L*PDstandardNth1gt12 + 
                J22L*J23L*PDstandardNth22gt12 + J23L*J32L*PDstandardNth23gt12 + J22L*J33L*PDstandardNth23gt12 + 
                dJ223L*PDstandardNth2gt12 + J32L*J33L*PDstandardNth33gt12 + dJ323L*PDstandardNth3gt12)) + 
          (gt12L*J31L + gt11L*J32L)*PDstandardNth3Xt1 + (gt22L*J31L + gt12L*J32L)*PDstandardNth3Xt2 + 
          (gt23L*J31L + gt13L*J32L)*PDstandardNth3Xt3 + 
          (Gt112*gt11L + Gt111*gt12L + gt12L*Gt212 + Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)*Xtn1 + 
          (gt11L*Gt122 + Gt112*gt12L + gt12L*Gt222 + Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322)*Xtn2 + 
          (gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323)*Xtn3 + 
          2*((Gt123*gt12L*Gt211 + Gt113*gt12L*Gt212 + 2*Gt112*gt12L*Gt213 + gt12L*Gt212*Gt223 + Gt212*Gt213*gt22L + 
                Gt211*Gt223*gt22L + gt12L*Gt133*Gt311 + gt22L*Gt233*Gt311 + Gt113*gt13L*Gt312 + gt12L*Gt233*Gt312 + 
                Gt213*gt23L*Gt312 + gt11L*(2*Gt112*Gt113 + Gt123*Gt212 + Gt133*Gt312) + 2*Gt112*gt13L*Gt313 + 
                Gt212*gt23L*Gt313 + Gt111*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + gt13L*Gt212*Gt323 + 
                Gt211*gt23L*Gt323 + gt23L*Gt311*Gt333 + gt13L*Gt312*Gt333 + Gt312*Gt313*gt33L)*gtu31 + 
             (Gt123*gt12L*Gt212 + 2*Gt122*gt12L*Gt213 + Gt113*gt12L*Gt222 + gt12L*Gt222*Gt223 + Gt213*Gt222*gt22L + 
                Gt212*Gt223*gt22L + gt12L*Gt133*Gt312 + gt22L*Gt233*Gt312 + 2*Gt122*gt13L*Gt313 + Gt222*gt23L*Gt313 + 
                Gt112*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + Gt113*gt13L*Gt322 + gt12L*Gt233*Gt322 + 
                Gt213*gt23L*Gt322 + gt11L*(2*Gt113*Gt122 + Gt123*Gt222 + Gt133*Gt322) + gt13L*Gt222*Gt323 + 
                Gt212*gt23L*Gt323 + gt23L*Gt312*Gt333 + gt13L*Gt322*Gt333 + Gt313*Gt322*gt33L)*gtu32 + 
             gtu11*(3*Gt112*gt12L*Gt211 + 2*Gt211*Gt212*gt22L + Gt113*gt12L*Gt311 + 2*Gt112*gt13L*Gt311 + 
                Gt213*gt22L*Gt311 + Gt212*gt23L*Gt311 + gt13L*Gt212*Gt312 + gt12L*Gt213*Gt312 + 2*Gt211*gt23L*Gt312 + 
                gt11L*(2*Gt111*Gt112 + Gt112*Gt212 + Gt113*Gt312) + 
                Gt111*(gt12L*Gt212 + Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + gt23L*Gt311*Gt313 + gt13L*Gt312*Gt313 + 
                Gt311*Gt312*gt33L + gt12L*SQR(Gt111) + gt12L*SQR(Gt212)) + 
             gtu21*(Gt112*gt11L*Gt222 + Gt112*Gt211*gt22L + Gt211*Gt222*gt22L + 2*Gt122*gt13L*Gt311 + Gt112*gt23L*Gt311 + 
                Gt222*gt23L*Gt311 + gt13L*Gt222*Gt312 + Gt213*gt22L*Gt312 + Gt212*gt23L*Gt312 + gt23L*Gt312*Gt313 + 
                Gt113*gt11L*Gt322 + Gt211*gt23L*Gt322 + gt13L*Gt313*Gt322 + 
                Gt111*(2*gt11L*Gt122 + Gt112*gt12L + gt12L*Gt222 + gt13L*Gt322) + 
                gt12L*(2*Gt122*Gt211 + Gt112*Gt212 + Gt212*Gt222 + Gt113*Gt312 + Gt213*Gt322) + Gt311*Gt322*gt33L + 
                gt22L*SQR(Gt212)) + gtu22*(gt11L*Gt122*Gt222 + 2*Gt212*Gt222*gt22L + 2*Gt122*gt13L*Gt312 + 
                Gt223*gt22L*Gt312 + Gt222*gt23L*Gt312 + gt11L*Gt123*Gt322 + gt13L*Gt222*Gt322 + 2*Gt212*gt23L*Gt322 + 
                Gt112*(2*gt11L*Gt122 + gt12L*Gt222 + Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322) + gt23L*Gt312*Gt323 + 
                gt13L*Gt322*Gt323 + Gt312*Gt322*gt33L + gt12L*SQR(Gt112) + 
                gt12L*(3*Gt122*Gt212 + Gt123*Gt312 + Gt223*Gt322 + SQR(Gt222))) + 
             gtu33*(gt11L*Gt123*Gt223 + 2*Gt213*Gt223*gt22L + 2*Gt123*gt13L*Gt313 + gt22L*Gt233*Gt313 + Gt223*gt23L*Gt313 + 
                gt11L*Gt133*Gt323 + gt13L*Gt223*Gt323 + 2*Gt213*gt23L*Gt323 + 
                Gt113*(2*gt11L*Gt123 + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323) + gt23L*Gt313*Gt333 + 
                gt13L*Gt323*Gt333 + Gt313*Gt323*gt33L + gt12L*SQR(Gt113) + 
                gt12L*(3*Gt123*Gt213 + Gt133*Gt313 + Gt233*Gt323 + SQR(Gt223))) + 
             gtu21*(Gt122*gt12L*Gt211 + 3*Gt112*gt12L*Gt212 + gt12L*Gt212*Gt222 + Gt211*Gt222*gt22L + Gt123*gt12L*Gt311 + 
                Gt223*gt22L*Gt311 + 3*Gt112*gt13L*Gt312 + gt12L*Gt223*Gt312 + 2*Gt212*gt23L*Gt312 + 
                Gt111*(Gt112*gt12L + Gt212*gt22L + gt23L*Gt312) + gt13L*Gt212*Gt322 + Gt211*gt23L*Gt322 + 
                gt23L*Gt311*Gt323 + gt13L*Gt312*Gt323 + gt11L*(Gt122*Gt212 + Gt123*Gt312 + 2*SQR(Gt112)) + 
                gt22L*SQR(Gt212) + gt33L*SQR(Gt312)) + 
             gtu31*(Gt112*gt11L*Gt223 + Gt113*Gt211*gt22L + Gt212*Gt213*gt22L + Gt211*Gt223*gt22L + 2*Gt123*gt13L*Gt311 + 
                Gt113*gt23L*Gt311 + Gt223*gt23L*Gt311 + gt13L*Gt223*Gt312 + Gt213*gt23L*Gt312 + Gt213*gt22L*Gt313 + 
                Gt113*gt11L*Gt323 + Gt211*gt23L*Gt323 + gt13L*Gt313*Gt323 + 
                Gt111*(2*gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + gt13L*Gt323) + 
                gt12L*(2*Gt123*Gt211 + Gt112*Gt213 + Gt212*Gt223 + Gt113*Gt313 + Gt213*Gt323) + Gt311*Gt323*gt33L + 
                gt23L*SQR(Gt313)) + gtu32*(gt11L*Gt122*Gt223 + Gt113*Gt212*gt22L + Gt213*Gt222*gt22L + Gt212*Gt223*gt22L + 
                2*Gt123*gt13L*Gt312 + Gt113*gt23L*Gt312 + Gt223*gt23L*Gt312 + Gt223*gt22L*Gt313 + gt13L*Gt223*Gt322 + 
                Gt213*gt23L*Gt322 + gt11L*Gt123*Gt323 + Gt212*gt23L*Gt323 + gt23L*Gt313*Gt323 + 
                Gt112*(2*gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + gt13L*Gt323) + 
                gt12L*(Gt122*Gt213 + Gt123*(2*Gt212 + Gt313) + Gt223*(Gt222 + Gt323)) + Gt312*Gt323*gt33L + gt13L*SQR(Gt323)
                )) - gtu11*(2*J11L*J21L*PDstandardNth12gt12 + 2*J11L*J31L*PDstandardNth13gt12 + dJ111L*PDstandardNth1gt12 + 
             2*J21L*J31L*PDstandardNth23gt12 + dJ211L*PDstandardNth2gt12 + dJ311L*PDstandardNth3gt12 + 
             PDstandardNth11gt12*SQR(J11L) + PDstandardNth22gt12*SQR(J21L) + PDstandardNth33gt12*SQR(J31L)) - 
          gtu22*(2*J12L*J22L*PDstandardNth12gt12 + 2*J12L*J32L*PDstandardNth13gt12 + dJ122L*PDstandardNth1gt12 + 
             2*J22L*J32L*PDstandardNth23gt12 + dJ222L*PDstandardNth2gt12 + dJ322L*PDstandardNth3gt12 + 
             PDstandardNth11gt12*SQR(J12L) + PDstandardNth22gt12*SQR(J22L) + PDstandardNth33gt12*SQR(J32L)) - 
          gtu33*(2*J13L*J23L*PDstandardNth12gt12 + 2*J13L*J33L*PDstandardNth13gt12 + dJ133L*PDstandardNth1gt12 + 
             2*J23L*J33L*PDstandardNth23gt12 + dJ233L*PDstandardNth2gt12 + dJ333L*PDstandardNth3gt12 + 
             PDstandardNth11gt12*SQR(J13L) + PDstandardNth22gt12*SQR(J23L) + PDstandardNth33gt12*SQR(J33L)));
    
    CCTK_REAL const Rt13  =  khalf*((gt13L*J11L + gt11L*J13L)*PDstandardNth1Xt1 + (gt23L*J11L + gt12L*J13L)*PDstandardNth1Xt2 + 
          (gt33L*J11L + gt13L*J13L)*PDstandardNth1Xt3 + (gt13L*J21L + gt11L*J23L)*PDstandardNth2Xt1 + 
          (gt23L*J21L + gt12L*J23L)*PDstandardNth2Xt2 + (gt33L*J21L + gt13L*J23L)*PDstandardNth2Xt3 - 
          2*(gtu21*(J11L*J12L*PDstandardNth11gt13 + J12L*J21L*PDstandardNth12gt13 + J11L*J22L*PDstandardNth12gt13 + 
                J12L*J31L*PDstandardNth13gt13 + J11L*J32L*PDstandardNth13gt13 + dJ112L*PDstandardNth1gt13 + 
                J21L*J22L*PDstandardNth22gt13 + J22L*J31L*PDstandardNth23gt13 + J21L*J32L*PDstandardNth23gt13 + 
                dJ212L*PDstandardNth2gt13 + J31L*J32L*PDstandardNth33gt13 + dJ312L*PDstandardNth3gt13) + 
             gtu31*(J11L*J13L*PDstandardNth11gt13 + J13L*J21L*PDstandardNth12gt13 + J11L*J23L*PDstandardNth12gt13 + 
                J13L*J31L*PDstandardNth13gt13 + J11L*J33L*PDstandardNth13gt13 + dJ113L*PDstandardNth1gt13 + 
                J21L*J23L*PDstandardNth22gt13 + J23L*J31L*PDstandardNth23gt13 + J21L*J33L*PDstandardNth23gt13 + 
                dJ213L*PDstandardNth2gt13 + J31L*J33L*PDstandardNth33gt13 + dJ313L*PDstandardNth3gt13) + 
             gtu32*(J12L*J13L*PDstandardNth11gt13 + J13L*J22L*PDstandardNth12gt13 + J12L*J23L*PDstandardNth12gt13 + 
                J13L*J32L*PDstandardNth13gt13 + J12L*J33L*PDstandardNth13gt13 + dJ123L*PDstandardNth1gt13 + 
                J22L*J23L*PDstandardNth22gt13 + J23L*J32L*PDstandardNth23gt13 + J22L*J33L*PDstandardNth23gt13 + 
                dJ223L*PDstandardNth2gt13 + J32L*J33L*PDstandardNth33gt13 + dJ323L*PDstandardNth3gt13)) + 
          (gt13L*J31L + gt11L*J33L)*PDstandardNth3Xt1 + (gt23L*J31L + gt12L*J33L)*PDstandardNth3Xt2 + 
          (gt33L*J31L + gt13L*J33L)*PDstandardNth3Xt3 + 
          (Gt113*gt11L + Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L + gt13L*Gt313 + Gt311*gt33L)*Xtn1 + 
          (gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + gt13L*Gt323 + Gt312*gt33L)*Xtn2 + 
          (gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + Gt213*gt23L + gt13L*Gt333 + Gt313*gt33L)*Xtn3 + 
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
                Gt112*Gt311*gt33L + Gt212*Gt312*gt33L + Gt312*Gt313*gt33L + Gt311*Gt323*gt33L + gt23L*SQR(Gt212)) + 
             gtu32*(Gt123*gt13L*Gt212 + 2*Gt123*gt12L*Gt213 + Gt113*gt12L*Gt223 + Gt213*Gt223*gt22L + Gt212*Gt223*gt23L + 
                Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + 2*Gt123*gt13L*Gt313 + Gt223*gt23L*Gt313 + Gt113*gt13L*Gt323 + 
                gt13L*Gt223*Gt323 + gt12L*Gt233*Gt323 + Gt213*gt23L*Gt323 + 
                gt11L*(2*Gt113*Gt123 + Gt123*Gt223 + Gt133*Gt323) + gt13L*Gt323*Gt333 + Gt212*Gt323*gt33L + 
                Gt313*Gt323*gt33L + Gt312*Gt333*gt33L + Gt112*(Gt113*gt13L + Gt213*gt23L + Gt313*gt33L) + gt12L*SQR(Gt223))\
              + gtu11*(2*Gt113*gt12L*Gt211 + Gt112*gt13L*Gt211 + gt12L*Gt212*Gt213 + Gt211*Gt213*gt22L + 
                Gt211*Gt212*gt23L + 3*Gt113*gt13L*Gt311 + 2*Gt213*gt23L*Gt311 + gt13L*Gt213*Gt312 + gt12L*Gt213*Gt313 + 
                Gt211*gt23L*Gt313 + gt11L*(2*Gt111*Gt113 + Gt112*Gt213 + Gt113*Gt313) + Gt211*Gt312*gt33L + 
                2*Gt311*Gt313*gt33L + Gt111*(gt12L*Gt213 + Gt211*gt23L + gt13L*Gt313 + Gt311*gt33L) + gt13L*SQR(Gt111) + 
                gt13L*SQR(Gt313)) + gtu31*(Gt112*gt13L*Gt213 + Gt112*gt11L*Gt233 + Gt211*gt22L*Gt233 + Gt113*Gt211*gt23L + 
                Gt212*Gt213*gt23L + 2*Gt133*gt13L*Gt311 + Gt233*gt23L*Gt311 + gt13L*Gt233*Gt312 + Gt113*gt13L*Gt313 + 
                Gt213*gt23L*Gt313 + Gt113*gt11L*Gt333 + Gt211*gt23L*Gt333 + gt13L*Gt313*Gt333 + 
                Gt111*(2*gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + gt13L*Gt333) + 
                gt12L*(2*Gt133*Gt211 + Gt212*Gt233 + Gt213*Gt333) + Gt113*Gt311*gt33L + Gt213*Gt312*gt33L + 
                Gt311*Gt333*gt33L + gt33L*SQR(Gt313)) + 
             gtu31*(Gt123*gt13L*Gt211 + 3*Gt113*gt12L*Gt213 + gt12L*Gt213*Gt223 + Gt211*Gt223*gt23L + Gt133*gt13L*Gt311 + 
                Gt233*gt23L*Gt311 + 3*Gt113*gt13L*Gt313 + gt12L*Gt233*Gt313 + 2*Gt213*gt23L*Gt313 + gt13L*Gt213*Gt323 + 
                gt13L*Gt313*Gt333 + Gt211*Gt323*gt33L + Gt311*Gt333*gt33L + 
                Gt111*(Gt113*gt13L + Gt213*gt23L + Gt313*gt33L) + gt11L*(Gt123*Gt213 + Gt133*Gt313 + 2*SQR(Gt113)) + 
                gt22L*SQR(Gt213) + gt33L*SQR(Gt313)) + 
             gtu22*(2*Gt123*gt12L*Gt212 + Gt122*gt13L*Gt212 + gt12L*Gt222*Gt223 + Gt212*Gt223*gt22L + Gt212*Gt222*gt23L + 
                3*Gt123*gt13L*Gt312 + 2*Gt223*gt23L*Gt312 + gt13L*Gt223*Gt322 + gt12L*Gt223*Gt323 + Gt212*gt23L*Gt323 + 
                gt11L*(2*Gt112*Gt123 + Gt122*Gt223 + Gt123*Gt323) + Gt212*Gt322*gt33L + 2*Gt312*Gt323*gt33L + 
                Gt112*(gt12L*Gt223 + Gt212*gt23L + gt13L*Gt323 + Gt312*gt33L) + gt13L*SQR(Gt112) + gt13L*SQR(Gt323)) + 
             gtu33*(2*gt12L*Gt133*Gt213 + Gt123*gt13L*Gt213 + gt11L*Gt123*Gt233 + gt12L*Gt223*Gt233 + Gt213*gt22L*Gt233 + 
                Gt213*Gt223*gt23L + 3*Gt133*gt13L*Gt313 + 2*Gt233*gt23L*Gt313 + gt13L*Gt233*Gt323 + gt11L*Gt133*Gt333 + 
                gt12L*Gt233*Gt333 + Gt213*gt23L*Gt333 + Gt213*Gt323*gt33L + 2*Gt313*Gt333*gt33L + 
                Gt113*(2*gt11L*Gt133 + gt12L*Gt233 + Gt213*gt23L + gt13L*Gt333 + Gt313*gt33L) + gt13L*SQR(Gt113) + 
                gt13L*SQR(Gt333))) - gtu11*(2*J11L*J21L*PDstandardNth12gt13 + 2*J11L*J31L*PDstandardNth13gt13 + 
             dJ111L*PDstandardNth1gt13 + 2*J21L*J31L*PDstandardNth23gt13 + dJ211L*PDstandardNth2gt13 + 
             dJ311L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J11L) + PDstandardNth22gt13*SQR(J21L) + 
             PDstandardNth33gt13*SQR(J31L)) - gtu22*(2*J12L*J22L*PDstandardNth12gt13 + 2*J12L*J32L*PDstandardNth13gt13 + 
             dJ122L*PDstandardNth1gt13 + 2*J22L*J32L*PDstandardNth23gt13 + dJ222L*PDstandardNth2gt13 + 
             dJ322L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J12L) + PDstandardNth22gt13*SQR(J22L) + 
             PDstandardNth33gt13*SQR(J32L)) - gtu33*(2*J13L*J23L*PDstandardNth12gt13 + 2*J13L*J33L*PDstandardNth13gt13 + 
             dJ133L*PDstandardNth1gt13 + 2*J23L*J33L*PDstandardNth23gt13 + dJ233L*PDstandardNth2gt13 + 
             dJ333L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J13L) + PDstandardNth22gt13*SQR(J23L) + 
             PDstandardNth33gt13*SQR(J33L)));
    
    CCTK_REAL const Rt22  =  (Gt223*(3*Gt112*gt12L + 6*Gt212*gt22L + 4*gt23L*Gt312) + 
           Gt123*(Gt112*gt11L + gt12L*(2*Gt111 + 4*Gt212) + 2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + 
           Gt112*(gt11L*Gt123 + gt12L*(2*Gt113 + Gt223) + 2*(Gt213*gt22L + gt23L*Gt313) + gt13L*Gt323) + 
           2*(Gt113*gt12L*Gt323 + Gt312*(gt12L*Gt133 + gt22L*Gt233 + gt23L*Gt333)) + 
           Gt323*(Gt112*gt13L + 4*Gt212*gt23L + 2*(Gt213*gt22L + gt23L*Gt313 + Gt312*gt33L)))*gtu31 + 
        J12L*(gt12L*PDstandardNth1Xt1 + gt22L*PDstandardNth1Xt2 + gt23L*PDstandardNth1Xt3) + 
        J22L*(gt12L*PDstandardNth2Xt1 + gt22L*PDstandardNth2Xt2 + gt23L*PDstandardNth2Xt3) + 
        J32L*(gt12L*PDstandardNth3Xt1 + gt22L*PDstandardNth3Xt2 + gt23L*PDstandardNth3Xt3) + 
        (Gt112*gt12L + Gt212*gt22L + gt23L*Gt312)*Xtn1 + (Gt122*gt12L + Gt222*gt22L + gt23L*Gt322)*Xtn2 + 
        (Gt123*gt12L + Gt223*gt22L + gt23L*Gt323)*Xtn3 + 
        gtu21*(Gt222*(3*Gt112*gt12L + 6*Gt212*gt22L + 4*gt23L*Gt312) + 
           Gt122*(Gt112*gt11L + gt12L*(2*Gt111 + 4*Gt212) + 2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + 
           Gt112*(gt11L*Gt122 + gt12L*Gt222 + 2*(Gt212*gt22L + gt23L*Gt312) + gt13L*Gt322) + 
           Gt322*(Gt112*gt13L + 4*Gt212*gt23L + 2*(Gt213*gt22L + gt23L*Gt313 + Gt312*gt33L)) + 
           2*(Gt312*(Gt123*gt12L + Gt223*gt22L + gt23L*Gt323) + gt12L*(Gt113*Gt322 + SQR(Gt112)))) + 
        gtu11*(Gt112*(gt12L*(2*Gt111 + 4*Gt212) + 2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + 
           Gt312*(2*(Gt113*gt12L + Gt213*gt22L) + gt23L*(4*Gt212 + 2*Gt313)) + gt11L*SQR(Gt112) + 3*gt22L*SQR(Gt212) + 
           gt33L*SQR(Gt312)) + gtu22*(Gt122*(gt12L*(2*Gt112 + 4*Gt222) + 2*(Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322)) + 
           Gt322*(2*(Gt123*gt12L + Gt223*gt22L) + gt23L*(4*Gt222 + 2*Gt323)) + gt11L*SQR(Gt122) + 3*gt22L*SQR(Gt222) + 
           gt33L*SQR(Gt322)) + gtu33*(Gt123*(gt12L*(2*Gt113 + 4*Gt223) + 2*(Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323)) + 
           Gt323*(2*(gt12L*Gt133 + gt22L*Gt233) + gt23L*(4*Gt223 + 2*Gt333)) + gt11L*SQR(Gt123) + 3*gt22L*SQR(Gt223) + 
           gt33L*SQR(Gt323)) + gtu32*(gt22L*(2*(Gt122*Gt213 + Gt233*Gt322) + Gt223*(6*Gt222 + 2*Gt323)) + 
           4*(gt12L*(Gt123*Gt222 + Gt122*Gt223) + gt23L*(Gt223*Gt322 + Gt222*Gt323)) + 
           2*(Gt123*(Gt112*gt12L + Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322) + gt12L*(Gt133*Gt322 + Gt123*Gt323) + 
              Gt122*(gt11L*Gt123 + Gt113*gt12L + gt23L*Gt313 + gt13L*Gt323) + Gt322*(gt23L*Gt333 + Gt323*gt33L) + 
              gt23L*SQR(Gt323))) + khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt22 + J12L*J21L*PDstandardNth12gt22 + 
                 J11L*J22L*PDstandardNth12gt22 + J12L*J31L*PDstandardNth13gt22 + J11L*J32L*PDstandardNth13gt22 + 
                 dJ112L*PDstandardNth1gt22 + J21L*J22L*PDstandardNth22gt22 + J22L*J31L*PDstandardNth23gt22 + 
                 J21L*J32L*PDstandardNth23gt22 + dJ212L*PDstandardNth2gt22 + J31L*J32L*PDstandardNth33gt22 + 
                 dJ312L*PDstandardNth3gt22) + gtu31*(J11L*J13L*PDstandardNth11gt22 + J13L*J21L*PDstandardNth12gt22 + 
                 J11L*J23L*PDstandardNth12gt22 + J13L*J31L*PDstandardNth13gt22 + J11L*J33L*PDstandardNth13gt22 + 
                 dJ113L*PDstandardNth1gt22 + J21L*J23L*PDstandardNth22gt22 + J23L*J31L*PDstandardNth23gt22 + 
                 J21L*J33L*PDstandardNth23gt22 + dJ213L*PDstandardNth2gt22 + J31L*J33L*PDstandardNth33gt22 + 
                 dJ313L*PDstandardNth3gt22) + gtu32*(J12L*J13L*PDstandardNth11gt22 + J13L*J22L*PDstandardNth12gt22 + 
                 J12L*J23L*PDstandardNth12gt22 + J13L*J32L*PDstandardNth13gt22 + J12L*J33L*PDstandardNth13gt22 + 
                 dJ123L*PDstandardNth1gt22 + J22L*J23L*PDstandardNth22gt22 + J23L*J32L*PDstandardNth23gt22 + 
                 J22L*J33L*PDstandardNth23gt22 + dJ223L*PDstandardNth2gt22 + J32L*J33L*PDstandardNth33gt22 + 
                 dJ323L*PDstandardNth3gt22)) - gtu11*
            (2*J11L*J21L*PDstandardNth12gt22 + 2*J11L*J31L*PDstandardNth13gt22 + dJ111L*PDstandardNth1gt22 + 
              2*J21L*J31L*PDstandardNth23gt22 + dJ211L*PDstandardNth2gt22 + dJ311L*PDstandardNth3gt22 + 
              PDstandardNth11gt22*SQR(J11L) + PDstandardNth22gt22*SQR(J21L) + PDstandardNth33gt22*SQR(J31L)) - 
           gtu22*(2*J12L*J22L*PDstandardNth12gt22 + 2*J12L*J32L*PDstandardNth13gt22 + dJ122L*PDstandardNth1gt22 + 
              2*J22L*J32L*PDstandardNth23gt22 + dJ222L*PDstandardNth2gt22 + dJ322L*PDstandardNth3gt22 + 
              PDstandardNth11gt22*SQR(J12L) + PDstandardNth22gt22*SQR(J22L) + PDstandardNth33gt22*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt22 + 2*J13L*J33L*PDstandardNth13gt22 + dJ133L*PDstandardNth1gt22 + 
              2*J23L*J33L*PDstandardNth23gt22 + dJ233L*PDstandardNth2gt22 + dJ333L*PDstandardNth3gt22 + 
              PDstandardNth11gt22*SQR(J13L) + PDstandardNth22gt22*SQR(J23L) + PDstandardNth33gt22*SQR(J33L)));
    
    CCTK_REAL const Rt23  =  khalf*((gt13L*J12L + gt12L*J13L)*PDstandardNth1Xt1 + (gt23L*J12L + gt22L*J13L)*PDstandardNth1Xt2 + 
          (gt33L*J12L + gt23L*J13L)*PDstandardNth1Xt3 + (gt13L*J22L + gt12L*J23L)*PDstandardNth2Xt1 + 
          (gt23L*J22L + gt22L*J23L)*PDstandardNth2Xt2 + (gt33L*J22L + gt23L*J23L)*PDstandardNth2Xt3 - 
          2*(gtu21*(J11L*J12L*PDstandardNth11gt23 + J12L*J21L*PDstandardNth12gt23 + J11L*J22L*PDstandardNth12gt23 + 
                J12L*J31L*PDstandardNth13gt23 + J11L*J32L*PDstandardNth13gt23 + dJ112L*PDstandardNth1gt23 + 
                J21L*J22L*PDstandardNth22gt23 + J22L*J31L*PDstandardNth23gt23 + J21L*J32L*PDstandardNth23gt23 + 
                dJ212L*PDstandardNth2gt23 + J31L*J32L*PDstandardNth33gt23 + dJ312L*PDstandardNth3gt23) + 
             gtu31*(J11L*J13L*PDstandardNth11gt23 + J13L*J21L*PDstandardNth12gt23 + J11L*J23L*PDstandardNth12gt23 + 
                J13L*J31L*PDstandardNth13gt23 + J11L*J33L*PDstandardNth13gt23 + dJ113L*PDstandardNth1gt23 + 
                J21L*J23L*PDstandardNth22gt23 + J23L*J31L*PDstandardNth23gt23 + J21L*J33L*PDstandardNth23gt23 + 
                dJ213L*PDstandardNth2gt23 + J31L*J33L*PDstandardNth33gt23 + dJ313L*PDstandardNth3gt23) + 
             gtu32*(J12L*J13L*PDstandardNth11gt23 + J13L*J22L*PDstandardNth12gt23 + J12L*J23L*PDstandardNth12gt23 + 
                J13L*J32L*PDstandardNth13gt23 + J12L*J33L*PDstandardNth13gt23 + dJ123L*PDstandardNth1gt23 + 
                J22L*J23L*PDstandardNth22gt23 + J23L*J32L*PDstandardNth23gt23 + J22L*J33L*PDstandardNth23gt23 + 
                dJ223L*PDstandardNth2gt23 + J32L*J33L*PDstandardNth33gt23 + dJ323L*PDstandardNth3gt23)) + 
          (gt13L*J32L + gt12L*J33L)*PDstandardNth3Xt1 + (gt23L*J32L + gt22L*J33L)*PDstandardNth3Xt2 + 
          (gt33L*J32L + gt23L*J33L)*PDstandardNth3Xt3 + 
          (Gt113*gt12L + Gt112*gt13L + Gt213*gt22L + Gt212*gt23L + gt23L*Gt313 + Gt312*gt33L)*Xtn1 + 
          (Gt123*gt12L + Gt122*gt13L + Gt223*gt22L + Gt222*gt23L + gt23L*Gt323 + Gt322*gt33L)*Xtn2 + 
          (gt12L*Gt133 + Gt123*gt13L + gt22L*Gt233 + Gt223*gt23L + gt23L*Gt333 + Gt323*gt33L)*Xtn3 + 
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
                Gt112*(Gt113*gt12L + Gt212*gt23L + Gt312*gt33L) + gt13L*SQR(Gt112)) + 
             gtu31*(2*Gt213*Gt223*gt22L + Gt112*Gt213*gt23L + Gt212*Gt223*gt23L + Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + 
                gt12L*Gt133*Gt313 + gt22L*Gt233*Gt313 + Gt223*gt23L*Gt313 + Gt123*(2*gt12L*Gt213 + gt13L*(Gt212 + Gt313)) + 
                2*Gt213*gt23L*Gt323 + Gt113*(gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + 
                   gt13L*Gt323) + gt23L*Gt313*Gt333 + Gt112*Gt313*gt33L + Gt212*Gt323*gt33L + Gt313*Gt323*gt33L + 
                Gt312*Gt333*gt33L + gt12L*SQR(Gt113)) + 
             gtu11*(Gt112*Gt113*gt11L + Gt111*Gt113*gt12L + Gt111*Gt112*gt13L + Gt113*gt12L*Gt212 + Gt112*gt13L*Gt212 + 
                2*Gt112*gt12L*Gt213 + Gt113*Gt211*gt22L + 2*Gt212*Gt213*gt22L + Gt112*Gt211*gt23L + Gt113*gt23L*Gt311 + 
                2*Gt113*gt13L*Gt312 + 3*Gt213*gt23L*Gt312 + Gt113*gt12L*Gt313 + Gt112*gt13L*Gt313 + Gt213*gt22L*Gt313 + 
                Gt212*gt23L*Gt313 + Gt112*Gt311*gt33L + Gt212*Gt312*gt33L + 2*Gt312*Gt313*gt33L + gt23L*SQR(Gt212) + 
                gt23L*SQR(Gt313)) + gtu22*(gt11L*Gt122*Gt123 + Gt112*Gt123*gt12L + Gt112*Gt122*gt13L + Gt123*gt12L*Gt222 + 
                Gt122*gt13L*Gt222 + 2*Gt122*gt12L*Gt223 + Gt123*Gt212*gt22L + 2*Gt222*Gt223*gt22L + Gt122*Gt212*gt23L + 
                Gt123*gt23L*Gt312 + 2*Gt123*gt13L*Gt322 + 3*Gt223*gt23L*Gt322 + Gt123*gt12L*Gt323 + Gt122*gt13L*Gt323 + 
                Gt223*gt22L*Gt323 + Gt222*gt23L*Gt323 + Gt122*Gt312*gt33L + Gt222*Gt322*gt33L + 2*Gt322*Gt323*gt33L + 
                gt23L*SQR(Gt222) + gt23L*SQR(Gt323)) + 
             gtu32*(gt11L*Gt122*Gt133 + Gt112*gt12L*Gt133 + Gt112*Gt123*gt13L + gt12L*Gt133*Gt222 + Gt122*gt13L*Gt223 + 
                Gt133*Gt212*gt22L + 2*Gt122*gt12L*Gt233 + 2*Gt222*gt22L*Gt233 + Gt123*Gt212*gt23L + Gt222*Gt223*gt23L + 
                Gt133*gt23L*Gt312 + Gt133*gt13L*Gt322 + 2*Gt233*gt23L*Gt322 + Gt123*gt13L*Gt323 + Gt223*gt23L*Gt323 + 
                Gt123*gt12L*Gt333 + Gt122*gt13L*Gt333 + Gt223*gt22L*Gt333 + Gt222*gt23L*Gt333 + gt23L*Gt323*Gt333 + 
                Gt123*Gt312*gt33L + Gt223*Gt322*gt33L + Gt322*Gt333*gt33L + gt33L*SQR(Gt323)) + 
             gtu32*(Gt113*Gt123*gt12L + Gt113*Gt122*gt13L + Gt123*gt13L*Gt222 + 3*Gt123*gt12L*Gt223 + Gt123*Gt213*gt22L + 
                Gt122*Gt213*gt23L + Gt222*Gt223*gt23L + Gt123*gt23L*Gt313 + Gt133*gt13L*Gt322 + Gt233*gt23L*Gt322 + 
                gt12L*Gt133*Gt323 + 2*Gt123*gt13L*Gt323 + gt22L*Gt233*Gt323 + 3*Gt223*gt23L*Gt323 + gt23L*Gt323*Gt333 + 
                Gt122*Gt313*gt33L + Gt222*Gt323*gt33L + Gt322*Gt333*gt33L + gt11L*SQR(Gt123) + 2*gt22L*SQR(Gt223) + 
                gt33L*SQR(Gt323)) + gtu33*(gt11L*Gt123*Gt133 + Gt113*gt12L*Gt133 + Gt113*Gt123*gt13L + gt12L*Gt133*Gt223 + 
                Gt123*gt13L*Gt223 + Gt133*Gt213*gt22L + 2*Gt123*gt12L*Gt233 + 2*Gt223*gt22L*Gt233 + Gt123*Gt213*gt23L + 
                Gt133*gt23L*Gt313 + 2*Gt133*gt13L*Gt323 + 3*Gt233*gt23L*Gt323 + gt12L*Gt133*Gt333 + Gt123*gt13L*Gt333 + 
                gt22L*Gt233*Gt333 + Gt223*gt23L*Gt333 + Gt123*Gt313*gt33L + Gt223*Gt323*gt33L + 2*Gt323*Gt333*gt33L + 
                gt23L*SQR(Gt223) + gt23L*SQR(Gt333))) - 
          gtu11*(2*J11L*J21L*PDstandardNth12gt23 + 2*J11L*J31L*PDstandardNth13gt23 + dJ111L*PDstandardNth1gt23 + 
             2*J21L*J31L*PDstandardNth23gt23 + dJ211L*PDstandardNth2gt23 + dJ311L*PDstandardNth3gt23 + 
             PDstandardNth11gt23*SQR(J11L) + PDstandardNth22gt23*SQR(J21L) + PDstandardNth33gt23*SQR(J31L)) - 
          gtu22*(2*J12L*J22L*PDstandardNth12gt23 + 2*J12L*J32L*PDstandardNth13gt23 + dJ122L*PDstandardNth1gt23 + 
             2*J22L*J32L*PDstandardNth23gt23 + dJ222L*PDstandardNth2gt23 + dJ322L*PDstandardNth3gt23 + 
             PDstandardNth11gt23*SQR(J12L) + PDstandardNth22gt23*SQR(J22L) + PDstandardNth33gt23*SQR(J32L)) - 
          gtu33*(2*J13L*J23L*PDstandardNth12gt23 + 2*J13L*J33L*PDstandardNth13gt23 + dJ133L*PDstandardNth1gt23 + 
             2*J23L*J33L*PDstandardNth23gt23 + dJ233L*PDstandardNth2gt23 + dJ333L*PDstandardNth3gt23 + 
             PDstandardNth11gt23*SQR(J13L) + PDstandardNth22gt23*SQR(J23L) + PDstandardNth33gt23*SQR(J33L)));
    
    CCTK_REAL const Rt33  =  (4*((Gt123*gt13L + Gt223*gt23L)*Gt313 + (Gt113*gt13L + Gt213*gt23L)*Gt323) + 
           (2*Gt213*Gt322 + 6*Gt313*Gt323)*gt33L + 2*
            (gt13L*(Gt122*Gt213 + Gt112*Gt223) + Gt213*(Gt223*gt22L + Gt222*gt23L) + 
              Gt123*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L + Gt311*gt33L) + Gt223*(Gt212*gt23L + Gt312*gt33L) + 
              Gt113*(gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L)))*gtu21 + 
        J13L*(gt13L*PDstandardNth1Xt1 + gt23L*PDstandardNth1Xt2 + gt33L*PDstandardNth1Xt3) + 
        J23L*(gt13L*PDstandardNth2Xt1 + gt23L*PDstandardNth2Xt2 + gt33L*PDstandardNth2Xt3) + 
        J33L*(gt13L*PDstandardNth3Xt1 + gt23L*PDstandardNth3Xt2 + gt33L*PDstandardNth3Xt3) + 
        (Gt113*gt13L + Gt213*gt23L + Gt313*gt33L)*Xtn1 + (Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)*Xtn2 + 
        (Gt133*gt13L + Gt233*gt23L + Gt333*gt33L)*Xtn3 + 
        gtu31*(Gt133*(Gt113*gt11L + 2*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L) + 4*gt13L*Gt313) + 
           Gt333*(3*Gt113*gt13L + 4*Gt213*gt23L + 6*Gt313*gt33L) + 
           Gt233*(Gt113*gt12L + 4*gt23L*Gt313 + 2*(Gt213*gt22L + Gt212*gt23L + Gt312*gt33L)) + 
           Gt113*(gt11L*Gt133 + gt12L*Gt233 + gt13L*Gt333 + 2*(Gt213*gt23L + Gt313*gt33L)) + 
           2*(Gt133*Gt311*gt33L + Gt213*(Gt223*gt23L + Gt323*gt33L) + gt13L*(Gt123*Gt213 + Gt112*Gt233 + SQR(Gt113)))) + 
        gtu32*(4*((Gt133*gt13L + Gt233*gt23L)*Gt323 + (Gt123*gt13L + Gt223*gt23L)*Gt333) + 
           Gt323*(2*Gt223 + 6*Gt333)*gt33L + 2*(Gt133*(Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L) + 
              Gt123*(gt11L*Gt133 + gt13L*(Gt113 + Gt223) + gt12L*Gt233 + Gt213*gt23L + Gt313*gt33L) + 
              Gt233*(Gt122*gt13L + Gt223*gt22L + Gt222*gt23L + Gt322*gt33L) + gt23L*SQR(Gt223))) + 
        gtu11*(4*(Gt113*gt13L + Gt213*gt23L)*Gt313 + 
           2*(Gt113*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L + Gt311*gt33L) + 
              Gt213*(Gt112*gt13L + Gt212*gt23L + Gt312*gt33L)) + gt11L*SQR(Gt113) + gt22L*SQR(Gt213) + 3*gt33L*SQR(Gt313))\
         + gtu22*(4*(Gt123*gt13L + Gt223*gt23L)*Gt323 + 
           2*(Gt123*(Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L) + 
              Gt223*(Gt122*gt13L + Gt222*gt23L + Gt322*gt33L)) + gt11L*SQR(Gt123) + gt22L*SQR(Gt223) + 3*gt33L*SQR(Gt323))\
         + gtu33*(4*(Gt133*gt13L + Gt233*gt23L)*Gt333 + 
           2*(Gt133*(Gt113*gt13L + gt12L*Gt233 + Gt213*gt23L + Gt313*gt33L) + 
              Gt233*(Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)) + gt11L*SQR(Gt133) + gt22L*SQR(Gt233) + 3*gt33L*SQR(Gt333))\
         + khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt33 + J12L*J21L*PDstandardNth12gt33 + 
                 J11L*J22L*PDstandardNth12gt33 + J12L*J31L*PDstandardNth13gt33 + J11L*J32L*PDstandardNth13gt33 + 
                 dJ112L*PDstandardNth1gt33 + J21L*J22L*PDstandardNth22gt33 + J22L*J31L*PDstandardNth23gt33 + 
                 J21L*J32L*PDstandardNth23gt33 + dJ212L*PDstandardNth2gt33 + J31L*J32L*PDstandardNth33gt33 + 
                 dJ312L*PDstandardNth3gt33) + gtu31*(J11L*J13L*PDstandardNth11gt33 + J13L*J21L*PDstandardNth12gt33 + 
                 J11L*J23L*PDstandardNth12gt33 + J13L*J31L*PDstandardNth13gt33 + J11L*J33L*PDstandardNth13gt33 + 
                 dJ113L*PDstandardNth1gt33 + J21L*J23L*PDstandardNth22gt33 + J23L*J31L*PDstandardNth23gt33 + 
                 J21L*J33L*PDstandardNth23gt33 + dJ213L*PDstandardNth2gt33 + J31L*J33L*PDstandardNth33gt33 + 
                 dJ313L*PDstandardNth3gt33) + gtu32*(J12L*J13L*PDstandardNth11gt33 + J13L*J22L*PDstandardNth12gt33 + 
                 J12L*J23L*PDstandardNth12gt33 + J13L*J32L*PDstandardNth13gt33 + J12L*J33L*PDstandardNth13gt33 + 
                 dJ123L*PDstandardNth1gt33 + J22L*J23L*PDstandardNth22gt33 + J23L*J32L*PDstandardNth23gt33 + 
                 J22L*J33L*PDstandardNth23gt33 + dJ223L*PDstandardNth2gt33 + J32L*J33L*PDstandardNth33gt33 + 
                 dJ323L*PDstandardNth3gt33)) - gtu11*
            (2*J11L*J21L*PDstandardNth12gt33 + 2*J11L*J31L*PDstandardNth13gt33 + dJ111L*PDstandardNth1gt33 + 
              2*J21L*J31L*PDstandardNth23gt33 + dJ211L*PDstandardNth2gt33 + dJ311L*PDstandardNth3gt33 + 
              PDstandardNth11gt33*SQR(J11L) + PDstandardNth22gt33*SQR(J21L) + PDstandardNth33gt33*SQR(J31L)) - 
           gtu22*(2*J12L*J22L*PDstandardNth12gt33 + 2*J12L*J32L*PDstandardNth13gt33 + dJ122L*PDstandardNth1gt33 + 
              2*J22L*J32L*PDstandardNth23gt33 + dJ222L*PDstandardNth2gt33 + dJ322L*PDstandardNth3gt33 + 
              PDstandardNth11gt33*SQR(J12L) + PDstandardNth22gt33*SQR(J22L) + PDstandardNth33gt33*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt33 + 2*J13L*J33L*PDstandardNth13gt33 + dJ133L*PDstandardNth1gt33 + 
              2*J23L*J33L*PDstandardNth23gt33 + dJ233L*PDstandardNth2gt33 + dJ333L*PDstandardNth3gt33 + 
              PDstandardNth11gt33*SQR(J13L) + PDstandardNth22gt33*SQR(J23L) + PDstandardNth33gt33*SQR(J33L)));
    
    CCTK_REAL const fac1  =  IfThen(conformalMethod,-(khalf*INV(phiL)),1);
    
    CCTK_REAL const cdphi1  =  fac1*(J11L*PDstandardNth1phi + J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    CCTK_REAL const cdphi2  =  fac1*(J12L*PDstandardNth1phi + J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    CCTK_REAL const cdphi3  =  fac1*(J13L*PDstandardNth1phi + J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
    CCTK_REAL const fac2  =  IfThen(conformalMethod,khalf*pow(phiL,-2),0);
    
    CCTK_REAL const cdphi211  =  fac1*((dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L)*PDstandardNth1phi + 
           2*(J11L*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + J21L*J31L*PDstandardNth23phi) + 
           (dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L)*PDstandardNth2phi + 
           (dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L)*PDstandardNth3phi + PDstandardNth11phi*SQR(J11L) + 
           PDstandardNth22phi*SQR(J21L) + PDstandardNth33phi*SQR(J31L)) + 
        fac2*SQR(J11L*PDstandardNth1phi + J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    CCTK_REAL const cdphi212  =  fac2*(J11L*PDstandardNth1phi + J21L*PDstandardNth2phi + J31L*PDstandardNth3phi)*
         (J12L*PDstandardNth1phi + J22L*PDstandardNth2phi + J32L*PDstandardNth3phi) + 
        fac1*(J12L*(J11L*PDstandardNth11phi + J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + 
           J11L*(J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
           (dJ112L - Gt112*J11L - Gt212*J12L - Gt312*J13L)*PDstandardNth1phi + 
           J22L*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi) + 
           (dJ212L - Gt112*J21L - Gt212*J22L - Gt312*J23L)*PDstandardNth2phi + 
           J32L*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi) + 
           (dJ312L - Gt112*J31L - Gt212*J32L - Gt312*J33L)*PDstandardNth3phi);
    
    CCTK_REAL const cdphi213  =  fac2*(J11L*PDstandardNth1phi + J21L*PDstandardNth2phi + J31L*PDstandardNth3phi)*
         (J13L*PDstandardNth1phi + J23L*PDstandardNth2phi + J33L*PDstandardNth3phi) + 
        fac1*(J13L*(J11L*PDstandardNth11phi + J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + 
           J11L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
           (dJ113L - Gt113*J11L - Gt213*J12L - Gt313*J13L)*PDstandardNth1phi + 
           J23L*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi) + 
           (dJ213L - Gt113*J21L - Gt213*J22L - Gt313*J23L)*PDstandardNth2phi + 
           J33L*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi) + 
           (dJ313L - Gt113*J31L - Gt213*J32L - Gt313*J33L)*PDstandardNth3phi);
    
    CCTK_REAL const cdphi222  =  fac1*((dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L)*PDstandardNth1phi + 
           2*(J12L*(J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + J22L*J32L*PDstandardNth23phi) + 
           (dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L)*PDstandardNth2phi + 
           (dJ322L - Gt122*J31L - Gt222*J32L - Gt322*J33L)*PDstandardNth3phi + PDstandardNth11phi*SQR(J12L) + 
           PDstandardNth22phi*SQR(J22L) + PDstandardNth33phi*SQR(J32L)) + 
        fac2*SQR(J12L*PDstandardNth1phi + J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    CCTK_REAL const cdphi223  =  fac2*(J12L*PDstandardNth1phi + J22L*PDstandardNth2phi + J32L*PDstandardNth3phi)*
         (J13L*PDstandardNth1phi + J23L*PDstandardNth2phi + J33L*PDstandardNth3phi) + 
        fac1*(J13L*(J12L*PDstandardNth11phi + J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
           J12L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
           (dJ123L - Gt123*J11L - Gt223*J12L - Gt323*J13L)*PDstandardNth1phi + 
           J23L*(J22L*PDstandardNth22phi + J32L*PDstandardNth23phi) + 
           (dJ223L - Gt123*J21L - Gt223*J22L - Gt323*J23L)*PDstandardNth2phi + 
           J33L*(J22L*PDstandardNth23phi + J32L*PDstandardNth33phi) + 
           (dJ323L - Gt123*J31L - Gt223*J32L - Gt323*J33L)*PDstandardNth3phi);
    
    CCTK_REAL const cdphi233  =  fac1*((dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L)*PDstandardNth1phi + 
           2*(J13L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + J23L*J33L*PDstandardNth23phi) + 
           (dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L)*PDstandardNth2phi + 
           (dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L)*PDstandardNth3phi + PDstandardNth11phi*SQR(J13L) + 
           PDstandardNth22phi*SQR(J23L) + PDstandardNth33phi*SQR(J33L)) + 
        fac2*SQR(J13L*PDstandardNth1phi + J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
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
    
    CCTK_REAL const Atu11  =  Atm11*gtu11 + Atm12*gtu21 + Atm13*gtu31;
    
    CCTK_REAL const Atu21  =  Atm11*gtu21 + Atm12*gtu22 + Atm13*gtu32;
    
    CCTK_REAL const Atu31  =  Atm11*gtu31 + Atm12*gtu32 + Atm13*gtu33;
    
    CCTK_REAL const Atu22  =  Atm21*gtu21 + Atm22*gtu22 + Atm23*gtu32;
    
    CCTK_REAL const Atu32  =  Atm21*gtu31 + Atm22*gtu32 + Atm23*gtu33;
    
    CCTK_REAL const Atu33  =  Atm31*gtu31 + Atm32*gtu32 + Atm33*gtu33;
    
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
    
    CCTK_REAL const trS  =  gu11*T11 + gu22*T22 + 2*(gu21*T12 + gu31*T13 + gu32*T23) + gu33*T33;
    
    CCTK_REAL const phirhsL  =  PDupwindNth1(phi, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(phi, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(phi, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) + 
        (J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
           J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
           J33L*PDstandardNth3beta3)*IfThen(conformalMethod,-(kthird*phiL),0.16666666666666666) + 
        alphaL*trKL*IfThen(conformalMethod,kthird*phiL,-0.16666666666666666);
    
    CCTK_REAL const gt11rhsL  =  -2*alphaL*At11L + PDupwindNth1(gt11, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(gt11, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(gt11, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) - 
        gt11L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + 
        2*(J11L*(gt11L*PDstandardNth1beta1 + gt12L*PDstandardNth1beta2 + gt13L*PDstandardNth1beta3) + 
           J21L*(gt11L*PDstandardNth2beta1 + gt12L*PDstandardNth2beta2 + gt13L*PDstandardNth2beta3) + 
           J31L*(gt11L*PDstandardNth3beta1 + gt12L*PDstandardNth3beta2 + gt13L*PDstandardNth3beta3));
    
    CCTK_REAL const gt12rhsL  =  -2*alphaL*At12L + PDupwindNth1(gt12, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(gt12, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(gt12, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) + 
        (gt12L*J11L + gt11L*J12L)*PDstandardNth1beta1 + (gt22L*J11L + gt12L*J12L)*PDstandardNth1beta2 + 
        (gt23L*J11L + gt13L*J12L)*PDstandardNth1beta3 + (gt12L*J21L + gt11L*J22L)*PDstandardNth2beta1 + 
        (gt22L*J21L + gt12L*J22L)*PDstandardNth2beta2 + (gt23L*J21L + gt13L*J22L)*PDstandardNth2beta3 + 
        (gt12L*J31L + gt11L*J32L)*PDstandardNth3beta1 + (gt22L*J31L + gt12L*J32L)*PDstandardNth3beta2 + 
        (gt23L*J31L + gt13L*J32L)*PDstandardNth3beta3 - 
        gt12L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3);
    
    CCTK_REAL const gt13rhsL  =  -2*alphaL*At13L + PDupwindNth1(gt13, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(gt13, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(gt13, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) + 
        (gt13L*J11L + gt11L*J13L)*PDstandardNth1beta1 + (gt23L*J11L + gt12L*J13L)*PDstandardNth1beta2 + 
        (gt33L*J11L + gt13L*J13L)*PDstandardNth1beta3 + (gt13L*J21L + gt11L*J23L)*PDstandardNth2beta1 + 
        (gt23L*J21L + gt12L*J23L)*PDstandardNth2beta2 + (gt33L*J21L + gt13L*J23L)*PDstandardNth2beta3 + 
        (gt13L*J31L + gt11L*J33L)*PDstandardNth3beta1 + (gt23L*J31L + gt12L*J33L)*PDstandardNth3beta2 + 
        (gt33L*J31L + gt13L*J33L)*PDstandardNth3beta3 - 
        gt13L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3);
    
    CCTK_REAL const gt22rhsL  =  -2*alphaL*At22L + PDupwindNth1(gt22, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(gt22, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(gt22, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) - 
        gt22L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + 
        2*(J12L*(gt12L*PDstandardNth1beta1 + gt22L*PDstandardNth1beta2 + gt23L*PDstandardNth1beta3) + 
           J22L*(gt12L*PDstandardNth2beta1 + gt22L*PDstandardNth2beta2 + gt23L*PDstandardNth2beta3) + 
           J32L*(gt12L*PDstandardNth3beta1 + gt22L*PDstandardNth3beta2 + gt23L*PDstandardNth3beta3));
    
    CCTK_REAL const gt23rhsL  =  -2*alphaL*At23L + PDupwindNth1(gt23, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(gt23, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(gt23, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) + 
        (gt13L*J12L + gt12L*J13L - gt23L*J11L*ktwothird)*PDstandardNth1beta1 + 
        (gt22L*J13L + gt23L*J12L*kthird)*PDstandardNth1beta2 + (gt33L*J12L + gt23L*J13L*kthird)*PDstandardNth1beta3 + 
        (gt13L*J22L + gt12L*J23L - gt23L*J21L*ktwothird)*PDstandardNth2beta1 + 
        (gt22L*J23L + gt23L*J22L*kthird)*PDstandardNth2beta2 + (gt33L*J22L + gt23L*J23L*kthird)*PDstandardNth2beta3 + 
        (gt13L*J32L + gt12L*J33L - gt23L*J31L*ktwothird)*PDstandardNth3beta1 + 
        (gt22L*J33L + gt23L*J32L*kthird)*PDstandardNth3beta2 + (gt33L*J32L + gt23L*J33L*kthird)*PDstandardNth3beta3;
    
    CCTK_REAL const gt33rhsL  =  -2*alphaL*At33L + PDupwindNth1(gt33, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(gt33, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(gt33, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) - 
        gt33L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + 
        2*(J13L*(gt13L*PDstandardNth1beta1 + gt23L*PDstandardNth1beta2 + gt33L*PDstandardNth1beta3) + 
           J23L*(gt13L*PDstandardNth2beta1 + gt23L*PDstandardNth2beta2 + gt33L*PDstandardNth2beta3) + 
           J33L*(gt13L*PDstandardNth3beta1 + gt23L*PDstandardNth3beta2 + gt33L*PDstandardNth3beta3));
    
    CCTK_REAL const Xt1rhsL  =  PDupwindNth1(Xt1, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(Xt1, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(Xt1, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) - 
        2*((Atu11*J11L + Atu21*J12L + Atu31*J13L)*PDstandardNth1alpha + 
           (Atu11*J21L + Atu21*J22L + Atu31*J23L)*PDstandardNth2alpha + 
           (Atu11*J31L + Atu21*J32L + Atu31*J33L)*PDstandardNth3alpha) + 
        2*(gtu21*(J11L*J12L*PDstandardNth11beta1 + J12L*J21L*PDstandardNth12beta1 + J11L*J22L*PDstandardNth12beta1 + 
              J12L*J31L*PDstandardNth13beta1 + J11L*J32L*PDstandardNth13beta1 + dJ112L*PDstandardNth1beta1 + 
              J21L*J22L*PDstandardNth22beta1 + J22L*J31L*PDstandardNth23beta1 + J21L*J32L*PDstandardNth23beta1 + 
              dJ212L*PDstandardNth2beta1 + J31L*J32L*PDstandardNth33beta1 + dJ312L*PDstandardNth3beta1) + 
           gtu31*(J11L*J13L*PDstandardNth11beta1 + J13L*J21L*PDstandardNth12beta1 + J11L*J23L*PDstandardNth12beta1 + 
              J13L*J31L*PDstandardNth13beta1 + J11L*J33L*PDstandardNth13beta1 + dJ113L*PDstandardNth1beta1 + 
              J21L*J23L*PDstandardNth22beta1 + J23L*J31L*PDstandardNth23beta1 + J21L*J33L*PDstandardNth23beta1 + 
              dJ213L*PDstandardNth2beta1 + J31L*J33L*PDstandardNth33beta1 + dJ313L*PDstandardNth3beta1) + 
           gtu32*(J12L*J13L*PDstandardNth11beta1 + J13L*J22L*PDstandardNth12beta1 + J12L*J23L*PDstandardNth12beta1 + 
              J13L*J32L*PDstandardNth13beta1 + J12L*J33L*PDstandardNth13beta1 + dJ123L*PDstandardNth1beta1 + 
              J22L*J23L*PDstandardNth22beta1 + J23L*J32L*PDstandardNth23beta1 + J22L*J33L*PDstandardNth23beta1 + 
              dJ223L*PDstandardNth2beta1 + J32L*J33L*PDstandardNth33beta1 + dJ323L*PDstandardNth3beta1) + 
           alphaL*(6*(Atu11*cdphi1 + Atu21*cdphi2 + Atu31*cdphi3) + Atu11*Gt111 + 2*Atu21*Gt112 + 2*Atu31*Gt113 + 
              Atu22*Gt122 + 2*Atu32*Gt123 + Atu33*Gt133 - 
              ktwothird*((gtu11*J11L + gtu21*J12L + gtu31*J13L)*PDstandardNth1trK + 
                 (gtu11*J21L + gtu21*J22L + gtu31*J23L)*PDstandardNth2trK + 
                 (gtu11*J31L + gtu21*J32L + gtu31*J33L)*PDstandardNth3trK))) - 
        50.26548245743669181540229413247204614715*alphaL*(gtu11*S1 + gtu21*S2 + gtu31*S3) + 
        ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3)*Xtn1 - 
        PDstandardNth1beta1*(J11L*Xtn1 + J12L*Xtn2 + J13L*Xtn3) - PDstandardNth2beta1*(J21L*Xtn1 + J22L*Xtn2 + J23L*Xtn3) - 
        PDstandardNth3beta1*(J31L*Xtn1 + J32L*Xtn2 + J33L*Xtn3) + 
        gtu11*(2*J11L*J21L*PDstandardNth12beta1 + 2*J11L*J31L*PDstandardNth13beta1 + dJ111L*PDstandardNth1beta1 + 
           2*J21L*J31L*PDstandardNth23beta1 + dJ211L*PDstandardNth2beta1 + dJ311L*PDstandardNth3beta1 + 
           PDstandardNth11beta1*SQR(J11L) + PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L)) + 
        gtu22*(2*J12L*J22L*PDstandardNth12beta1 + 2*J12L*J32L*PDstandardNth13beta1 + dJ122L*PDstandardNth1beta1 + 
           2*J22L*J32L*PDstandardNth23beta1 + dJ222L*PDstandardNth2beta1 + dJ322L*PDstandardNth3beta1 + 
           PDstandardNth11beta1*SQR(J12L) + PDstandardNth22beta1*SQR(J22L) + PDstandardNth33beta1*SQR(J32L)) + 
        gtu33*(2*J13L*J23L*PDstandardNth12beta1 + 2*J13L*J33L*PDstandardNth13beta1 + dJ133L*PDstandardNth1beta1 + 
           2*J23L*J33L*PDstandardNth23beta1 + dJ233L*PDstandardNth2beta1 + dJ333L*PDstandardNth3beta1 + 
           PDstandardNth11beta1*SQR(J13L) + PDstandardNth22beta1*SQR(J23L) + PDstandardNth33beta1*SQR(J33L)) + 
        kthird*(gtu11*(J11L*J12L*PDstandardNth11beta2 + J11L*J13L*PDstandardNth11beta3 + 2*J11L*J21L*PDstandardNth12beta1 + 
              J12L*J21L*PDstandardNth12beta2 + J11L*J22L*PDstandardNth12beta2 + J13L*J21L*PDstandardNth12beta3 + 
              J11L*J23L*PDstandardNth12beta3 + 2*J11L*J31L*PDstandardNth13beta1 + J12L*J31L*PDstandardNth13beta2 + 
              J11L*J32L*PDstandardNth13beta2 + J13L*J31L*PDstandardNth13beta3 + J11L*J33L*PDstandardNth13beta3 + 
              dJ111L*PDstandardNth1beta1 + dJ112L*PDstandardNth1beta2 + dJ113L*PDstandardNth1beta3 + 
              J21L*J22L*PDstandardNth22beta2 + J21L*J23L*PDstandardNth22beta3 + 2*J21L*J31L*PDstandardNth23beta1 + 
              J22L*J31L*PDstandardNth23beta2 + J21L*J32L*PDstandardNth23beta2 + J23L*J31L*PDstandardNth23beta3 + 
              J21L*J33L*PDstandardNth23beta3 + dJ211L*PDstandardNth2beta1 + dJ212L*PDstandardNth2beta2 + 
              dJ213L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta2 + J31L*J33L*PDstandardNth33beta3 + 
              dJ311L*PDstandardNth3beta1 + dJ312L*PDstandardNth3beta2 + dJ313L*PDstandardNth3beta3 + 
              PDstandardNth11beta1*SQR(J11L) + PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L)) + 
           gtu21*(J11L*J12L*PDstandardNth11beta1 + J12L*J13L*PDstandardNth11beta3 + J12L*J21L*PDstandardNth12beta1 + 
              J11L*J22L*PDstandardNth12beta1 + 2*J12L*J22L*PDstandardNth12beta2 + J13L*J22L*PDstandardNth12beta3 + 
              J12L*J23L*PDstandardNth12beta3 + J12L*J31L*PDstandardNth13beta1 + J11L*J32L*PDstandardNth13beta1 + 
              2*J12L*J32L*PDstandardNth13beta2 + J13L*J32L*PDstandardNth13beta3 + J12L*J33L*PDstandardNth13beta3 + 
              dJ112L*PDstandardNth1beta1 + dJ122L*PDstandardNth1beta2 + dJ123L*PDstandardNth1beta3 + 
              J21L*J22L*PDstandardNth22beta1 + J22L*J23L*PDstandardNth22beta3 + J22L*J31L*PDstandardNth23beta1 + 
              J21L*J32L*PDstandardNth23beta1 + 2*J22L*J32L*PDstandardNth23beta2 + J23L*J32L*PDstandardNth23beta3 + 
              J22L*J33L*PDstandardNth23beta3 + dJ212L*PDstandardNth2beta1 + dJ222L*PDstandardNth2beta2 + 
              dJ223L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta1 + J32L*J33L*PDstandardNth33beta3 + 
              dJ312L*PDstandardNth3beta1 + dJ322L*PDstandardNth3beta2 + dJ323L*PDstandardNth3beta3 + 
              PDstandardNth11beta2*SQR(J12L) + PDstandardNth22beta2*SQR(J22L) + PDstandardNth33beta2*SQR(J32L)) + 
           gtu31*(J11L*J13L*PDstandardNth11beta1 + J12L*J13L*PDstandardNth11beta2 + J13L*J21L*PDstandardNth12beta1 + 
              J11L*J23L*PDstandardNth12beta1 + J13L*J22L*PDstandardNth12beta2 + J12L*J23L*PDstandardNth12beta2 + 
              2*J13L*J23L*PDstandardNth12beta3 + J13L*J31L*PDstandardNth13beta1 + J11L*J33L*PDstandardNth13beta1 + 
              J13L*J32L*PDstandardNth13beta2 + J12L*J33L*PDstandardNth13beta2 + 2*J13L*J33L*PDstandardNth13beta3 + 
              dJ113L*PDstandardNth1beta1 + dJ123L*PDstandardNth1beta2 + dJ133L*PDstandardNth1beta3 + 
              J21L*J23L*PDstandardNth22beta1 + J22L*J23L*PDstandardNth22beta2 + J23L*J31L*PDstandardNth23beta1 + 
              J21L*J33L*PDstandardNth23beta1 + J23L*J32L*PDstandardNth23beta2 + J22L*J33L*PDstandardNth23beta2 + 
              2*J23L*J33L*PDstandardNth23beta3 + dJ213L*PDstandardNth2beta1 + dJ223L*PDstandardNth2beta2 + 
              dJ233L*PDstandardNth2beta3 + J31L*J33L*PDstandardNth33beta1 + J32L*J33L*PDstandardNth33beta2 + 
              dJ313L*PDstandardNth3beta1 + dJ323L*PDstandardNth3beta2 + dJ333L*PDstandardNth3beta3 + 
              PDstandardNth11beta3*SQR(J13L) + PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L)));
    
    CCTK_REAL const Xt2rhsL  =  PDupwindNth1(Xt2, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(Xt2, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(Xt2, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) - 
        2*((Atu21*J11L + Atu22*J12L + Atu32*J13L)*PDstandardNth1alpha + 
           (Atu21*J21L + Atu22*J22L + Atu32*J23L)*PDstandardNth2alpha + 
           (Atu21*J31L + Atu22*J32L + Atu32*J33L)*PDstandardNth3alpha) + 
        2*(gtu21*(J11L*J12L*PDstandardNth11beta2 + J12L*J21L*PDstandardNth12beta2 + J11L*J22L*PDstandardNth12beta2 + 
              J12L*J31L*PDstandardNth13beta2 + J11L*J32L*PDstandardNth13beta2 + dJ112L*PDstandardNth1beta2 + 
              J21L*J22L*PDstandardNth22beta2 + J22L*J31L*PDstandardNth23beta2 + J21L*J32L*PDstandardNth23beta2 + 
              dJ212L*PDstandardNth2beta2 + J31L*J32L*PDstandardNth33beta2 + dJ312L*PDstandardNth3beta2) + 
           gtu31*(J11L*J13L*PDstandardNth11beta2 + J13L*J21L*PDstandardNth12beta2 + J11L*J23L*PDstandardNth12beta2 + 
              J13L*J31L*PDstandardNth13beta2 + J11L*J33L*PDstandardNth13beta2 + dJ113L*PDstandardNth1beta2 + 
              J21L*J23L*PDstandardNth22beta2 + J23L*J31L*PDstandardNth23beta2 + J21L*J33L*PDstandardNth23beta2 + 
              dJ213L*PDstandardNth2beta2 + J31L*J33L*PDstandardNth33beta2 + dJ313L*PDstandardNth3beta2) + 
           gtu32*(J12L*J13L*PDstandardNth11beta2 + J13L*J22L*PDstandardNth12beta2 + J12L*J23L*PDstandardNth12beta2 + 
              J13L*J32L*PDstandardNth13beta2 + J12L*J33L*PDstandardNth13beta2 + dJ123L*PDstandardNth1beta2 + 
              J22L*J23L*PDstandardNth22beta2 + J23L*J32L*PDstandardNth23beta2 + J22L*J33L*PDstandardNth23beta2 + 
              dJ223L*PDstandardNth2beta2 + J32L*J33L*PDstandardNth33beta2 + dJ323L*PDstandardNth3beta2) + 
           alphaL*(6*(Atu21*cdphi1 + Atu22*cdphi2 + Atu32*cdphi3) + Atu11*Gt211 + 2*Atu21*Gt212 + 2*Atu31*Gt213 + 
              Atu22*Gt222 + 2*Atu32*Gt223 + Atu33*Gt233 - 
              ktwothird*((gtu21*J11L + gtu22*J12L + gtu32*J13L)*PDstandardNth1trK + 
                 (gtu21*J21L + gtu22*J22L + gtu32*J23L)*PDstandardNth2trK + 
                 (gtu21*J31L + gtu22*J32L + gtu32*J33L)*PDstandardNth3trK))) - 
        50.26548245743669181540229413247204614715*alphaL*(gtu21*S1 + gtu22*S2 + gtu32*S3) + 
        ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3)*Xtn2 - 
        PDstandardNth1beta2*(J11L*Xtn1 + J12L*Xtn2 + J13L*Xtn3) - PDstandardNth2beta2*(J21L*Xtn1 + J22L*Xtn2 + J23L*Xtn3) - 
        PDstandardNth3beta2*(J31L*Xtn1 + J32L*Xtn2 + J33L*Xtn3) + 
        gtu11*(2*J11L*J21L*PDstandardNth12beta2 + 2*J11L*J31L*PDstandardNth13beta2 + dJ111L*PDstandardNth1beta2 + 
           2*J21L*J31L*PDstandardNth23beta2 + dJ211L*PDstandardNth2beta2 + dJ311L*PDstandardNth3beta2 + 
           PDstandardNth11beta2*SQR(J11L) + PDstandardNth22beta2*SQR(J21L) + PDstandardNth33beta2*SQR(J31L)) + 
        gtu22*(2*J12L*J22L*PDstandardNth12beta2 + 2*J12L*J32L*PDstandardNth13beta2 + dJ122L*PDstandardNth1beta2 + 
           2*J22L*J32L*PDstandardNth23beta2 + dJ222L*PDstandardNth2beta2 + dJ322L*PDstandardNth3beta2 + 
           PDstandardNth11beta2*SQR(J12L) + PDstandardNth22beta2*SQR(J22L) + PDstandardNth33beta2*SQR(J32L)) + 
        gtu33*(2*J13L*J23L*PDstandardNth12beta2 + 2*J13L*J33L*PDstandardNth13beta2 + dJ133L*PDstandardNth1beta2 + 
           2*J23L*J33L*PDstandardNth23beta2 + dJ233L*PDstandardNth2beta2 + dJ333L*PDstandardNth3beta2 + 
           PDstandardNth11beta2*SQR(J13L) + PDstandardNth22beta2*SQR(J23L) + PDstandardNth33beta2*SQR(J33L)) + 
        kthird*(gtu21*(J11L*J12L*PDstandardNth11beta2 + J11L*J13L*PDstandardNth11beta3 + 2*J11L*J21L*PDstandardNth12beta1 + 
              J12L*J21L*PDstandardNth12beta2 + J11L*J22L*PDstandardNth12beta2 + J13L*J21L*PDstandardNth12beta3 + 
              J11L*J23L*PDstandardNth12beta3 + 2*J11L*J31L*PDstandardNth13beta1 + J12L*J31L*PDstandardNth13beta2 + 
              J11L*J32L*PDstandardNth13beta2 + J13L*J31L*PDstandardNth13beta3 + J11L*J33L*PDstandardNth13beta3 + 
              dJ111L*PDstandardNth1beta1 + dJ112L*PDstandardNth1beta2 + dJ113L*PDstandardNth1beta3 + 
              J21L*J22L*PDstandardNth22beta2 + J21L*J23L*PDstandardNth22beta3 + 2*J21L*J31L*PDstandardNth23beta1 + 
              J22L*J31L*PDstandardNth23beta2 + J21L*J32L*PDstandardNth23beta2 + J23L*J31L*PDstandardNth23beta3 + 
              J21L*J33L*PDstandardNth23beta3 + dJ211L*PDstandardNth2beta1 + dJ212L*PDstandardNth2beta2 + 
              dJ213L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta2 + J31L*J33L*PDstandardNth33beta3 + 
              dJ311L*PDstandardNth3beta1 + dJ312L*PDstandardNth3beta2 + dJ313L*PDstandardNth3beta3 + 
              PDstandardNth11beta1*SQR(J11L) + PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L)) + 
           gtu22*(J11L*J12L*PDstandardNth11beta1 + J12L*J13L*PDstandardNth11beta3 + J12L*J21L*PDstandardNth12beta1 + 
              J11L*J22L*PDstandardNth12beta1 + 2*J12L*J22L*PDstandardNth12beta2 + J13L*J22L*PDstandardNth12beta3 + 
              J12L*J23L*PDstandardNth12beta3 + J12L*J31L*PDstandardNth13beta1 + J11L*J32L*PDstandardNth13beta1 + 
              2*J12L*J32L*PDstandardNth13beta2 + J13L*J32L*PDstandardNth13beta3 + J12L*J33L*PDstandardNth13beta3 + 
              dJ112L*PDstandardNth1beta1 + dJ122L*PDstandardNth1beta2 + dJ123L*PDstandardNth1beta3 + 
              J21L*J22L*PDstandardNth22beta1 + J22L*J23L*PDstandardNth22beta3 + J22L*J31L*PDstandardNth23beta1 + 
              J21L*J32L*PDstandardNth23beta1 + 2*J22L*J32L*PDstandardNth23beta2 + J23L*J32L*PDstandardNth23beta3 + 
              J22L*J33L*PDstandardNth23beta3 + dJ212L*PDstandardNth2beta1 + dJ222L*PDstandardNth2beta2 + 
              dJ223L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta1 + J32L*J33L*PDstandardNth33beta3 + 
              dJ312L*PDstandardNth3beta1 + dJ322L*PDstandardNth3beta2 + dJ323L*PDstandardNth3beta3 + 
              PDstandardNth11beta2*SQR(J12L) + PDstandardNth22beta2*SQR(J22L) + PDstandardNth33beta2*SQR(J32L)) + 
           gtu32*(J11L*J13L*PDstandardNth11beta1 + J12L*J13L*PDstandardNth11beta2 + J13L*J21L*PDstandardNth12beta1 + 
              J11L*J23L*PDstandardNth12beta1 + J13L*J22L*PDstandardNth12beta2 + J12L*J23L*PDstandardNth12beta2 + 
              2*J13L*J23L*PDstandardNth12beta3 + J13L*J31L*PDstandardNth13beta1 + J11L*J33L*PDstandardNth13beta1 + 
              J13L*J32L*PDstandardNth13beta2 + J12L*J33L*PDstandardNth13beta2 + 2*J13L*J33L*PDstandardNth13beta3 + 
              dJ113L*PDstandardNth1beta1 + dJ123L*PDstandardNth1beta2 + dJ133L*PDstandardNth1beta3 + 
              J21L*J23L*PDstandardNth22beta1 + J22L*J23L*PDstandardNth22beta2 + J23L*J31L*PDstandardNth23beta1 + 
              J21L*J33L*PDstandardNth23beta1 + J23L*J32L*PDstandardNth23beta2 + J22L*J33L*PDstandardNth23beta2 + 
              2*J23L*J33L*PDstandardNth23beta3 + dJ213L*PDstandardNth2beta1 + dJ223L*PDstandardNth2beta2 + 
              dJ233L*PDstandardNth2beta3 + J31L*J33L*PDstandardNth33beta1 + J32L*J33L*PDstandardNth33beta2 + 
              dJ313L*PDstandardNth3beta1 + dJ323L*PDstandardNth3beta2 + dJ333L*PDstandardNth3beta3 + 
              PDstandardNth11beta3*SQR(J13L) + PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L)));
    
    CCTK_REAL const Xt3rhsL  =  PDupwindNth1(Xt3, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(Xt3, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(Xt3, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) - 
        2*((Atu31*J11L + Atu32*J12L + Atu33*J13L)*PDstandardNth1alpha + 
           (Atu31*J21L + Atu32*J22L + Atu33*J23L)*PDstandardNth2alpha + 
           (Atu31*J31L + Atu32*J32L + Atu33*J33L)*PDstandardNth3alpha) + 
        2*(gtu21*(J11L*J12L*PDstandardNth11beta3 + J12L*J21L*PDstandardNth12beta3 + J11L*J22L*PDstandardNth12beta3 + 
              J12L*J31L*PDstandardNth13beta3 + J11L*J32L*PDstandardNth13beta3 + dJ112L*PDstandardNth1beta3 + 
              J21L*J22L*PDstandardNth22beta3 + J22L*J31L*PDstandardNth23beta3 + J21L*J32L*PDstandardNth23beta3 + 
              dJ212L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta3 + dJ312L*PDstandardNth3beta3) + 
           gtu31*(J11L*J13L*PDstandardNth11beta3 + J13L*J21L*PDstandardNth12beta3 + J11L*J23L*PDstandardNth12beta3 + 
              J13L*J31L*PDstandardNth13beta3 + J11L*J33L*PDstandardNth13beta3 + dJ113L*PDstandardNth1beta3 + 
              J21L*J23L*PDstandardNth22beta3 + J23L*J31L*PDstandardNth23beta3 + J21L*J33L*PDstandardNth23beta3 + 
              dJ213L*PDstandardNth2beta3 + J31L*J33L*PDstandardNth33beta3 + dJ313L*PDstandardNth3beta3) + 
           gtu32*(J12L*J13L*PDstandardNth11beta3 + J13L*J22L*PDstandardNth12beta3 + J12L*J23L*PDstandardNth12beta3 + 
              J13L*J32L*PDstandardNth13beta3 + J12L*J33L*PDstandardNth13beta3 + dJ123L*PDstandardNth1beta3 + 
              J22L*J23L*PDstandardNth22beta3 + J23L*J32L*PDstandardNth23beta3 + J22L*J33L*PDstandardNth23beta3 + 
              dJ223L*PDstandardNth2beta3 + J32L*J33L*PDstandardNth33beta3 + dJ323L*PDstandardNth3beta3) + 
           alphaL*(6*(Atu31*cdphi1 + Atu32*cdphi2 + Atu33*cdphi3) + Atu11*Gt311 + 2*Atu21*Gt312 + 2*Atu31*Gt313 + 
              Atu22*Gt322 + 2*Atu32*Gt323 + Atu33*Gt333 - 
              ktwothird*((gtu31*J11L + gtu32*J12L + gtu33*J13L)*PDstandardNth1trK + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*PDstandardNth2trK + 
                 (gtu31*J31L + gtu32*J32L + gtu33*J33L)*PDstandardNth3trK))) - 
        50.26548245743669181540229413247204614715*alphaL*(gtu31*S1 + gtu32*S2 + gtu33*S3) + 
        ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3)*Xtn3 - 
        PDstandardNth1beta3*(J11L*Xtn1 + J12L*Xtn2 + J13L*Xtn3) - PDstandardNth2beta3*(J21L*Xtn1 + J22L*Xtn2 + J23L*Xtn3) - 
        PDstandardNth3beta3*(J31L*Xtn1 + J32L*Xtn2 + J33L*Xtn3) + 
        gtu11*(2*J11L*J21L*PDstandardNth12beta3 + 2*J11L*J31L*PDstandardNth13beta3 + dJ111L*PDstandardNth1beta3 + 
           2*J21L*J31L*PDstandardNth23beta3 + dJ211L*PDstandardNth2beta3 + dJ311L*PDstandardNth3beta3 + 
           PDstandardNth11beta3*SQR(J11L) + PDstandardNth22beta3*SQR(J21L) + PDstandardNth33beta3*SQR(J31L)) + 
        gtu22*(2*J12L*J22L*PDstandardNth12beta3 + 2*J12L*J32L*PDstandardNth13beta3 + dJ122L*PDstandardNth1beta3 + 
           2*J22L*J32L*PDstandardNth23beta3 + dJ222L*PDstandardNth2beta3 + dJ322L*PDstandardNth3beta3 + 
           PDstandardNth11beta3*SQR(J12L) + PDstandardNth22beta3*SQR(J22L) + PDstandardNth33beta3*SQR(J32L)) + 
        gtu33*(2*J13L*J23L*PDstandardNth12beta3 + 2*J13L*J33L*PDstandardNth13beta3 + dJ133L*PDstandardNth1beta3 + 
           2*J23L*J33L*PDstandardNth23beta3 + dJ233L*PDstandardNth2beta3 + dJ333L*PDstandardNth3beta3 + 
           PDstandardNth11beta3*SQR(J13L) + PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L)) + 
        kthird*(gtu31*(J11L*J12L*PDstandardNth11beta2 + J11L*J13L*PDstandardNth11beta3 + 2*J11L*J21L*PDstandardNth12beta1 + 
              J12L*J21L*PDstandardNth12beta2 + J11L*J22L*PDstandardNth12beta2 + J13L*J21L*PDstandardNth12beta3 + 
              J11L*J23L*PDstandardNth12beta3 + 2*J11L*J31L*PDstandardNth13beta1 + J12L*J31L*PDstandardNth13beta2 + 
              J11L*J32L*PDstandardNth13beta2 + J13L*J31L*PDstandardNth13beta3 + J11L*J33L*PDstandardNth13beta3 + 
              dJ111L*PDstandardNth1beta1 + dJ112L*PDstandardNth1beta2 + dJ113L*PDstandardNth1beta3 + 
              J21L*J22L*PDstandardNth22beta2 + J21L*J23L*PDstandardNth22beta3 + 2*J21L*J31L*PDstandardNth23beta1 + 
              J22L*J31L*PDstandardNth23beta2 + J21L*J32L*PDstandardNth23beta2 + J23L*J31L*PDstandardNth23beta3 + 
              J21L*J33L*PDstandardNth23beta3 + dJ211L*PDstandardNth2beta1 + dJ212L*PDstandardNth2beta2 + 
              dJ213L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta2 + J31L*J33L*PDstandardNth33beta3 + 
              dJ311L*PDstandardNth3beta1 + dJ312L*PDstandardNth3beta2 + dJ313L*PDstandardNth3beta3 + 
              PDstandardNth11beta1*SQR(J11L) + PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L)) + 
           gtu32*(J11L*J12L*PDstandardNth11beta1 + J12L*J13L*PDstandardNth11beta3 + J12L*J21L*PDstandardNth12beta1 + 
              J11L*J22L*PDstandardNth12beta1 + 2*J12L*J22L*PDstandardNth12beta2 + J13L*J22L*PDstandardNth12beta3 + 
              J12L*J23L*PDstandardNth12beta3 + J12L*J31L*PDstandardNth13beta1 + J11L*J32L*PDstandardNth13beta1 + 
              2*J12L*J32L*PDstandardNth13beta2 + J13L*J32L*PDstandardNth13beta3 + J12L*J33L*PDstandardNth13beta3 + 
              dJ112L*PDstandardNth1beta1 + dJ122L*PDstandardNth1beta2 + dJ123L*PDstandardNth1beta3 + 
              J21L*J22L*PDstandardNth22beta1 + J22L*J23L*PDstandardNth22beta3 + J22L*J31L*PDstandardNth23beta1 + 
              J21L*J32L*PDstandardNth23beta1 + 2*J22L*J32L*PDstandardNth23beta2 + J23L*J32L*PDstandardNth23beta3 + 
              J22L*J33L*PDstandardNth23beta3 + dJ212L*PDstandardNth2beta1 + dJ222L*PDstandardNth2beta2 + 
              dJ223L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta1 + J32L*J33L*PDstandardNth33beta3 + 
              dJ312L*PDstandardNth3beta1 + dJ322L*PDstandardNth3beta2 + dJ323L*PDstandardNth3beta3 + 
              PDstandardNth11beta2*SQR(J12L) + PDstandardNth22beta2*SQR(J22L) + PDstandardNth33beta2*SQR(J32L)) + 
           gtu33*(J11L*J13L*PDstandardNth11beta1 + J12L*J13L*PDstandardNth11beta2 + J13L*J21L*PDstandardNth12beta1 + 
              J11L*J23L*PDstandardNth12beta1 + J13L*J22L*PDstandardNth12beta2 + J12L*J23L*PDstandardNth12beta2 + 
              2*J13L*J23L*PDstandardNth12beta3 + J13L*J31L*PDstandardNth13beta1 + J11L*J33L*PDstandardNth13beta1 + 
              J13L*J32L*PDstandardNth13beta2 + J12L*J33L*PDstandardNth13beta2 + 2*J13L*J33L*PDstandardNth13beta3 + 
              dJ113L*PDstandardNth1beta1 + dJ123L*PDstandardNth1beta2 + dJ133L*PDstandardNth1beta3 + 
              J21L*J23L*PDstandardNth22beta1 + J22L*J23L*PDstandardNth22beta2 + J23L*J31L*PDstandardNth23beta1 + 
              J21L*J33L*PDstandardNth23beta1 + J23L*J32L*PDstandardNth23beta2 + J22L*J33L*PDstandardNth23beta2 + 
              2*J23L*J33L*PDstandardNth23beta3 + dJ213L*PDstandardNth2beta1 + dJ223L*PDstandardNth2beta2 + 
              dJ233L*PDstandardNth2beta3 + J31L*J33L*PDstandardNth33beta1 + J32L*J33L*PDstandardNth33beta2 + 
              dJ313L*PDstandardNth3beta1 + dJ323L*PDstandardNth3beta2 + dJ333L*PDstandardNth3beta3 + 
              PDstandardNth11beta3*SQR(J13L) + PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L)));
    
    CCTK_REAL const trKrhsL  =  PDupwindNth1(trK, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(trK, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(trK, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) + 
        (-(dJ111L*gu11) - 2*dJ112L*gu21 - dJ122L*gu22 - 2*dJ113L*gu31 - 2*dJ123L*gu32 - dJ133L*gu33 + G111*gu11*J11L + 
           2*G112*gu21*J11L + G122*gu22*J11L + 2*G113*gu31*J11L + 2*G123*gu32*J11L + G133*gu33*J11L + G211*gu11*J12L + 
           2*G212*gu21*J12L + G222*gu22*J12L + 2*G213*gu31*J12L + 2*G223*gu32*J12L + G233*gu33*J12L + G311*gu11*J13L + 
           2*G312*gu21*J13L + G322*gu22*J13L + 2*G313*gu31*J13L + 2*G323*gu32*J13L + G333*gu33*J13L)*PDstandardNth1alpha + 
        (-(dJ211L*gu11) - 2*dJ212L*gu21 - dJ222L*gu22 - 2*dJ213L*gu31 - 2*dJ223L*gu32 - dJ233L*gu33 + G111*gu11*J21L + 
           2*G112*gu21*J21L + G122*gu22*J21L + 2*G113*gu31*J21L + 2*G123*gu32*J21L + G133*gu33*J21L + G211*gu11*J22L + 
           2*G212*gu21*J22L + G222*gu22*J22L + 2*G213*gu31*J22L + 2*G223*gu32*J22L + G233*gu33*J22L + G311*gu11*J23L + 
           2*G312*gu21*J23L + G322*gu22*J23L + 2*G313*gu31*J23L + 2*G323*gu32*J23L + G333*gu33*J23L)*PDstandardNth2alpha - 
        dJ311L*gu11*PDstandardNth3alpha - dJ322L*gu22*PDstandardNth3alpha - 2*dJ313L*gu31*PDstandardNth3alpha - 
        2*dJ323L*gu32*PDstandardNth3alpha - dJ333L*gu33*PDstandardNth3alpha + G111*gu11*J31L*PDstandardNth3alpha + 
        G122*gu22*J31L*PDstandardNth3alpha + 2*G113*gu31*J31L*PDstandardNth3alpha + 2*G123*gu32*J31L*PDstandardNth3alpha + 
        G133*gu33*J31L*PDstandardNth3alpha + G211*gu11*J32L*PDstandardNth3alpha + 2*G212*gu21*J32L*PDstandardNth3alpha + 
        G222*gu22*J32L*PDstandardNth3alpha + 2*G213*gu31*J32L*PDstandardNth3alpha + 2*G223*gu32*J32L*PDstandardNth3alpha + 
        G233*gu33*J32L*PDstandardNth3alpha + G311*gu11*J33L*PDstandardNth3alpha + 2*G312*gu21*J33L*PDstandardNth3alpha + 
        G322*gu22*J33L*PDstandardNth3alpha + 2*G313*gu31*J33L*PDstandardNth3alpha + 2*G323*gu32*J33L*PDstandardNth3alpha + 
        G333*gu33*J33L*PDstandardNth3alpha - 2*(gu21*J11L*J12L*PDstandardNth11alpha + gu31*J11L*J13L*PDstandardNth11alpha + 
           gu32*J12L*J13L*PDstandardNth11alpha + gu11*J11L*J21L*PDstandardNth12alpha + 
           gu21*J12L*J21L*PDstandardNth12alpha + gu31*J13L*J21L*PDstandardNth12alpha + 
           gu21*J11L*J22L*PDstandardNth12alpha + gu22*J12L*J22L*PDstandardNth12alpha + 
           gu32*J13L*J22L*PDstandardNth12alpha + gu31*J11L*J23L*PDstandardNth12alpha + 
           gu32*J12L*J23L*PDstandardNth12alpha + gu33*J13L*J23L*PDstandardNth12alpha + 
           gu11*J11L*J31L*PDstandardNth13alpha + gu21*J12L*J31L*PDstandardNth13alpha + 
           gu31*J13L*J31L*PDstandardNth13alpha + gu21*J11L*J32L*PDstandardNth13alpha + 
           gu22*J12L*J32L*PDstandardNth13alpha + gu32*J13L*J32L*PDstandardNth13alpha + 
           gu31*J11L*J33L*PDstandardNth13alpha + gu32*J12L*J33L*PDstandardNth13alpha + 
           gu33*J13L*J33L*PDstandardNth13alpha + gu21*J21L*J22L*PDstandardNth22alpha + 
           gu31*J21L*J23L*PDstandardNth22alpha + gu32*J22L*J23L*PDstandardNth22alpha + 
           gu11*J21L*J31L*PDstandardNth23alpha + gu21*J22L*J31L*PDstandardNth23alpha + 
           gu31*J23L*J31L*PDstandardNth23alpha + gu21*J21L*J32L*PDstandardNth23alpha + 
           gu22*J22L*J32L*PDstandardNth23alpha + gu32*J23L*J32L*PDstandardNth23alpha + 
           gu31*J21L*J33L*PDstandardNth23alpha + gu32*J22L*J33L*PDstandardNth23alpha + 
           gu33*J23L*J33L*PDstandardNth23alpha + gu21*J31L*J32L*PDstandardNth33alpha + 
           gu31*J31L*J33L*PDstandardNth33alpha + gu32*J32L*J33L*PDstandardNth33alpha + dJ312L*gu21*PDstandardNth3alpha) + 
        2*(alphaL*Atm12*Atm21 + alphaL*Atm13*Atm31 + alphaL*Atm23*Atm32 + G112*gu21*J31L*PDstandardNth3alpha) + 
        12.56637061435917295385057353311801153679*alphaL*rho + 12.56637061435917295385057353311801153679*alphaL*trS + 
        alphaL*SQR(Atm11) + alphaL*SQR(Atm22) + alphaL*SQR(Atm33) - gu11*PDstandardNth11alpha*SQR(J11L) - 
        gu22*PDstandardNth11alpha*SQR(J12L) - gu33*PDstandardNth11alpha*SQR(J13L) - gu11*PDstandardNth22alpha*SQR(J21L) - 
        gu22*PDstandardNth22alpha*SQR(J22L) - gu33*PDstandardNth22alpha*SQR(J23L) - gu11*PDstandardNth33alpha*SQR(J31L) - 
        gu22*PDstandardNth33alpha*SQR(J32L) - gu33*PDstandardNth33alpha*SQR(J33L) + alphaL*kthird*SQR(trKL);
    
    CCTK_REAL const Ats11  =  (-dJ111L + G111*J11L + G211*J12L + G311*J13L)*PDstandardNth1alpha - 
        2*(J11L*J21L*PDstandardNth12alpha + J11L*J31L*PDstandardNth13alpha + J21L*J31L*PDstandardNth23alpha) + 
        (-dJ211L + G111*J21L + G211*J22L + G311*J23L)*PDstandardNth2alpha + 
        (-dJ311L + G111*J31L + G211*J32L + G311*J33L)*PDstandardNth3alpha + alphaL*R11 - PDstandardNth11alpha*SQR(J11L) - 
        PDstandardNth22alpha*SQR(J21L) - PDstandardNth33alpha*SQR(J31L);
    
    CCTK_REAL const Ats12  =  -(J11L*J12L*PDstandardNth11alpha) - J12L*J21L*PDstandardNth12alpha - J11L*J22L*PDstandardNth12alpha - 
        J12L*J31L*PDstandardNth13alpha - J11L*J32L*PDstandardNth13alpha + 
        (-dJ112L + G112*J11L + G212*J12L + G312*J13L)*PDstandardNth1alpha - J21L*J22L*PDstandardNth22alpha - 
        J22L*J31L*PDstandardNth23alpha - J21L*J32L*PDstandardNth23alpha + 
        (-dJ212L + G112*J21L + G212*J22L + G312*J23L)*PDstandardNth2alpha - J31L*J32L*PDstandardNth33alpha - 
        dJ312L*PDstandardNth3alpha + G112*J31L*PDstandardNth3alpha + G212*J32L*PDstandardNth3alpha + 
        G312*J33L*PDstandardNth3alpha + alphaL*R12;
    
    CCTK_REAL const Ats13  =  -(J11L*J13L*PDstandardNth11alpha) - J13L*J21L*PDstandardNth12alpha - J11L*J23L*PDstandardNth12alpha - 
        J13L*J31L*PDstandardNth13alpha - J11L*J33L*PDstandardNth13alpha + 
        (-dJ113L + G113*J11L + G213*J12L + G313*J13L)*PDstandardNth1alpha - J21L*J23L*PDstandardNth22alpha - 
        J23L*J31L*PDstandardNth23alpha - J21L*J33L*PDstandardNth23alpha + 
        (-dJ213L + G113*J21L + G213*J22L + G313*J23L)*PDstandardNth2alpha - J31L*J33L*PDstandardNth33alpha - 
        dJ313L*PDstandardNth3alpha + G113*J31L*PDstandardNth3alpha + G213*J32L*PDstandardNth3alpha + 
        G313*J33L*PDstandardNth3alpha + alphaL*R13;
    
    CCTK_REAL const Ats22  =  (-dJ122L + G122*J11L + G222*J12L + G322*J13L)*PDstandardNth1alpha - 
        2*(J12L*J22L*PDstandardNth12alpha + J12L*J32L*PDstandardNth13alpha + J22L*J32L*PDstandardNth23alpha) + 
        (-dJ222L + G122*J21L + G222*J22L + G322*J23L)*PDstandardNth2alpha + 
        (-dJ322L + G122*J31L + G222*J32L + G322*J33L)*PDstandardNth3alpha + alphaL*R22 - PDstandardNth11alpha*SQR(J12L) - 
        PDstandardNth22alpha*SQR(J22L) - PDstandardNth33alpha*SQR(J32L);
    
    CCTK_REAL const Ats23  =  -(J12L*J13L*PDstandardNth11alpha) - J13L*J22L*PDstandardNth12alpha - J12L*J23L*PDstandardNth12alpha - 
        J13L*J32L*PDstandardNth13alpha - J12L*J33L*PDstandardNth13alpha + 
        (-dJ123L + G123*J11L + G223*J12L + G323*J13L)*PDstandardNth1alpha - J22L*J23L*PDstandardNth22alpha - 
        J23L*J32L*PDstandardNth23alpha - J22L*J33L*PDstandardNth23alpha + 
        (-dJ223L + G123*J21L + G223*J22L + G323*J23L)*PDstandardNth2alpha - J32L*J33L*PDstandardNth33alpha - 
        dJ323L*PDstandardNth3alpha + G123*J31L*PDstandardNth3alpha + G223*J32L*PDstandardNth3alpha + 
        G323*J33L*PDstandardNth3alpha + alphaL*R23;
    
    CCTK_REAL const Ats33  =  (-dJ133L + G133*J11L + G233*J12L + G333*J13L)*PDstandardNth1alpha - 
        2*(J13L*J23L*PDstandardNth12alpha + J13L*J33L*PDstandardNth13alpha + J23L*J33L*PDstandardNth23alpha) + 
        (-dJ233L + G133*J21L + G233*J22L + G333*J23L)*PDstandardNth2alpha + 
        (-dJ333L + G133*J31L + G233*J32L + G333*J33L)*PDstandardNth3alpha + alphaL*R33 - PDstandardNth11alpha*SQR(J13L) - 
        PDstandardNth22alpha*SQR(J23L) - PDstandardNth33alpha*SQR(J33L);
    
    CCTK_REAL const trAts  =  Ats11*gu11 + Ats22*gu22 + 2*(Ats12*gu21 + Ats13*gu31 + Ats23*gu32) + Ats33*gu33;
    
    CCTK_REAL const At11rhsL  =  PDupwindNth1(At11, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(At11, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(At11, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) - 
        At11L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + 
        2*(J11L*(At11L*PDstandardNth1beta1 + At12L*PDstandardNth1beta2 + At13L*PDstandardNth1beta3) + 
           J21L*(At11L*PDstandardNth2beta1 + At12L*PDstandardNth2beta2 + At13L*PDstandardNth2beta3) + 
           J31L*(At11L*PDstandardNth3beta1 + At12L*PDstandardNth3beta2 + At13L*PDstandardNth3beta3)) + 
        alphaL*(-2*(At11L*Atm11 + At12L*Atm21 + At13L*Atm31) + At11L*trKL) + 
        em4phi*(Ats11 - g11*kthird*trAts - 25.13274122871834590770114706623602307358*alphaL*(T11 - g11*kthird*trS));
    
    CCTK_REAL const At12rhsL  =  PDupwindNth1(At12, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(At12, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(At12, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) + 
        (0.3333333333333333333333333333333333333333*At12L*J11L + At11L*J12L)*PDstandardNth1beta1 + 
        (At22L*J11L + 0.3333333333333333333333333333333333333333*At12L*J12L)*PDstandardNth1beta2 + 
        (At23L*J11L + At13L*J12L - 0.6666666666666666666666666666666666666667*At12L*J13L)*PDstandardNth1beta3 + 
        (0.3333333333333333333333333333333333333333*At12L*J21L + At11L*J22L)*PDstandardNth2beta1 + 
        (At22L*J21L + 0.3333333333333333333333333333333333333333*At12L*J22L)*PDstandardNth2beta2 + 
        (At23L*J21L + At13L*J22L - 0.6666666666666666666666666666666666666667*At12L*J23L)*PDstandardNth2beta3 + 
        (0.3333333333333333333333333333333333333333*At12L*J31L + At11L*J32L)*PDstandardNth3beta1 + 
        (At22L*J31L + 0.3333333333333333333333333333333333333333*At12L*J32L)*PDstandardNth3beta2 + 
        (At23L*J31L + At13L*J32L - 0.6666666666666666666666666666666666666667*At12L*J33L)*PDstandardNth3beta3 + 
        alphaL*(-2.*(At11L*Atm12 + At12L*Atm22 + At13L*Atm32) + At12L*trKL) + 
        em4phi*(Ats12 - 0.3333333333333333333333333333333333333333*g12*trAts + 
           alphaL*(-25.13274122871834590770114706623602307358*T12 + 8.377580409572781969233715688745341024526*g12*trS));
    
    CCTK_REAL const At13rhsL  =  PDupwindNth1(At13, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(At13, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(At13, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) + 
        (0.3333333333333333333333333333333333333333*At13L*J11L + At11L*J13L)*PDstandardNth1beta1 + 
        (At23L*J11L - 0.6666666666666666666666666666666666666667*At13L*J12L + At12L*J13L)*PDstandardNth1beta2 + 
        (At33L*J11L + 0.3333333333333333333333333333333333333333*At13L*J13L)*PDstandardNth1beta3 + 
        (0.3333333333333333333333333333333333333333*At13L*J21L + At11L*J23L)*PDstandardNth2beta1 + 
        (At23L*J21L - 0.6666666666666666666666666666666666666667*At13L*J22L + At12L*J23L)*PDstandardNth2beta2 + 
        (At33L*J21L + 0.3333333333333333333333333333333333333333*At13L*J23L)*PDstandardNth2beta3 + 
        (0.3333333333333333333333333333333333333333*At13L*J31L + At11L*J33L)*PDstandardNth3beta1 + 
        (At23L*J31L - 0.6666666666666666666666666666666666666667*At13L*J32L + At12L*J33L)*PDstandardNth3beta2 + 
        (At33L*J31L + 0.3333333333333333333333333333333333333333*At13L*J33L)*PDstandardNth3beta3 + 
        alphaL*(-2.*(At11L*Atm13 + At12L*Atm23 + At13L*Atm33) + At13L*trKL) + 
        em4phi*(Ats13 - 0.3333333333333333333333333333333333333333*g13*trAts + 
           alphaL*(-25.13274122871834590770114706623602307358*T13 + 8.377580409572781969233715688745341024526*g13*trS));
    
    CCTK_REAL const At22rhsL  =  PDupwindNth1(At22, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(At22, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(At22, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) - 
        At22L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + 
        2*(J12L*(At12L*PDstandardNth1beta1 + At22L*PDstandardNth1beta2 + At23L*PDstandardNth1beta3) + 
           J22L*(At12L*PDstandardNth2beta1 + At22L*PDstandardNth2beta2 + At23L*PDstandardNth2beta3) + 
           J32L*(At12L*PDstandardNth3beta1 + At22L*PDstandardNth3beta2 + At23L*PDstandardNth3beta3)) + 
        alphaL*(-2*(At12L*Atm12 + At22L*Atm22 + At23L*Atm32) + At22L*trKL) + 
        em4phi*(Ats22 - g22*kthird*trAts - 25.13274122871834590770114706623602307358*alphaL*(T22 - g22*kthird*trS));
    
    CCTK_REAL const At23rhsL  =  PDupwindNth1(At23, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(At23, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(At23, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) + 
        (-0.6666666666666666666666666666666666666667*At23L*J11L + At13L*J12L + At12L*J13L)*PDstandardNth1beta1 + 
        (0.3333333333333333333333333333333333333333*At23L*J12L + At22L*J13L)*PDstandardNth1beta2 + 
        (At33L*J12L + 0.3333333333333333333333333333333333333333*At23L*J13L)*PDstandardNth1beta3 + 
        (-0.6666666666666666666666666666666666666667*At23L*J21L + At13L*J22L + At12L*J23L)*PDstandardNth2beta1 + 
        (0.3333333333333333333333333333333333333333*At23L*J22L + At22L*J23L)*PDstandardNth2beta2 + 
        (At33L*J22L + 0.3333333333333333333333333333333333333333*At23L*J23L)*PDstandardNth2beta3 + 
        (-0.6666666666666666666666666666666666666667*At23L*J31L + At13L*J32L + At12L*J33L)*PDstandardNth3beta1 + 
        (0.3333333333333333333333333333333333333333*At23L*J32L + At22L*J33L)*PDstandardNth3beta2 + 
        (At33L*J32L + 0.3333333333333333333333333333333333333333*At23L*J33L)*PDstandardNth3beta3 + 
        alphaL*(-2.*(At12L*Atm13 + At22L*Atm23 + At23L*Atm33) + At23L*trKL) + 
        em4phi*(Ats23 - 0.3333333333333333333333333333333333333333*g23*trAts + 
           alphaL*(-25.13274122871834590770114706623602307358*T23 + 8.377580409572781969233715688745341024526*g23*trS));
    
    CCTK_REAL const At33rhsL  =  PDupwindNth1(At33, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
        PDupwindNth2(At33, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
        PDupwindNth3(At33, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L) - 
        At33L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
           J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
           J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + 
        2*(J13L*(At13L*PDstandardNth1beta1 + At23L*PDstandardNth1beta2 + At33L*PDstandardNth1beta3) + 
           J23L*(At13L*PDstandardNth2beta1 + At23L*PDstandardNth2beta2 + At33L*PDstandardNth2beta3) + 
           J33L*(At13L*PDstandardNth3beta1 + At23L*PDstandardNth3beta2 + At33L*PDstandardNth3beta3)) + 
        alphaL*(-2*(At13L*Atm13 + At23L*Atm23 + At33L*Atm33) + At33L*trKL) + 
        em4phi*(Ats33 - g33*kthird*trAts - 25.13274122871834590770114706623602307358*alphaL*(T33 - g33*kthird*trS));
    
    CCTK_REAL const alpharhsL  =  (PDupwindNth1(alpha, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
           PDupwindNth2(alpha, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
           PDupwindNth3(alpha, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L))*LapseAdvectionCoeff + 
        harmonicF*(AL*(-1 + LapseAdvectionCoeff) - LapseAdvectionCoeff*trKL)*pow(alphaL,harmonicN);
    
    CCTK_REAL const ArhsL  =  (-1 + LapseAdvectionCoeff)*(AL*AlphaDriver - trKrhsL);
    
    CCTK_REAL const beta1rhsL  =  (PDupwindNth1(beta1, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
           PDupwindNth2(beta1, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
           PDupwindNth3(beta1, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L))*ShiftAdvectionCoeff + 
        B1L*ShiftGammaCoeff;
    
    CCTK_REAL const beta2rhsL  =  (PDupwindNth1(beta2, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
           PDupwindNth2(beta2, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
           PDupwindNth3(beta2, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L))*ShiftAdvectionCoeff + 
        B2L*ShiftGammaCoeff;
    
    CCTK_REAL const beta3rhsL  =  (PDupwindNth1(beta3, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
           PDupwindNth2(beta3, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
           PDupwindNth3(beta3, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L))*ShiftAdvectionCoeff + 
        B3L*ShiftGammaCoeff;
    
    CCTK_REAL const B1rhsL  =  -(B1L*etaL) + (beta1L*((PDupwindNth1(B1, i, j, k) - PDupwindNth1(Xt1, i, j, k))*J11L + 
              (PDupwindNth2(B1, i, j, k) - PDupwindNth2(Xt1, i, j, k))*J21L + 
              (PDupwindNth3(B1, i, j, k) - PDupwindNth3(Xt1, i, j, k))*J31L) + 
           beta2L*((PDupwindNth1(B1, i, j, k) - PDupwindNth1(Xt1, i, j, k))*J12L + 
              (PDupwindNth2(B1, i, j, k) - PDupwindNth2(Xt1, i, j, k))*J22L + 
              (PDupwindNth3(B1, i, j, k) - PDupwindNth3(Xt1, i, j, k))*J32L) + 
           beta3L*((PDupwindNth1(B1, i, j, k) - PDupwindNth1(Xt1, i, j, k))*J13L + 
              (PDupwindNth2(B1, i, j, k) - PDupwindNth2(Xt1, i, j, k))*J23L + 
              (PDupwindNth3(B1, i, j, k) - PDupwindNth3(Xt1, i, j, k))*J33L))*ShiftAdvectionCoeff + Xt1rhsL;
    
    CCTK_REAL const B2rhsL  =  -(B2L*etaL) + (beta1L*((PDupwindNth1(B2, i, j, k) - PDupwindNth1(Xt2, i, j, k))*J11L + 
              (PDupwindNth2(B2, i, j, k) - PDupwindNth2(Xt2, i, j, k))*J21L + 
              (PDupwindNth3(B2, i, j, k) - PDupwindNth3(Xt2, i, j, k))*J31L) + 
           beta2L*((PDupwindNth1(B2, i, j, k) - PDupwindNth1(Xt2, i, j, k))*J12L + 
              (PDupwindNth2(B2, i, j, k) - PDupwindNth2(Xt2, i, j, k))*J22L + 
              (PDupwindNth3(B2, i, j, k) - PDupwindNth3(Xt2, i, j, k))*J32L) + 
           beta3L*((PDupwindNth1(B2, i, j, k) - PDupwindNth1(Xt2, i, j, k))*J13L + 
              (PDupwindNth2(B2, i, j, k) - PDupwindNth2(Xt2, i, j, k))*J23L + 
              (PDupwindNth3(B2, i, j, k) - PDupwindNth3(Xt2, i, j, k))*J33L))*ShiftAdvectionCoeff + Xt2rhsL;
    
    CCTK_REAL const B3rhsL  =  -(B3L*etaL) + (beta1L*((PDupwindNth1(B3, i, j, k) - PDupwindNth1(Xt3, i, j, k))*J11L + 
              (PDupwindNth2(B3, i, j, k) - PDupwindNth2(Xt3, i, j, k))*J21L + 
              (PDupwindNth3(B3, i, j, k) - PDupwindNth3(Xt3, i, j, k))*J31L) + 
           beta2L*((PDupwindNth1(B3, i, j, k) - PDupwindNth1(Xt3, i, j, k))*J12L + 
              (PDupwindNth2(B3, i, j, k) - PDupwindNth2(Xt3, i, j, k))*J22L + 
              (PDupwindNth3(B3, i, j, k) - PDupwindNth3(Xt3, i, j, k))*J32L) + 
           beta3L*((PDupwindNth1(B3, i, j, k) - PDupwindNth1(Xt3, i, j, k))*J13L + 
              (PDupwindNth2(B3, i, j, k) - PDupwindNth2(Xt3, i, j, k))*J23L + 
              (PDupwindNth3(B3, i, j, k) - PDupwindNth3(Xt3, i, j, k))*J33L))*ShiftAdvectionCoeff + Xt3rhsL;
    
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
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
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSN_MP_RHS);
}

void ML_BSSN_MP_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_RHS_Body);
}
