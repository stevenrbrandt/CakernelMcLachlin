/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (32-bit) (April 20, 2007) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

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

void ML_BSSN_MP_RHS_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  
  /* Declare finite differencing variables */
  CCTK_REAL dx = INITVALUE, dy = INITVALUE, dz = INITVALUE;
  CCTK_REAL dxi = INITVALUE, dyi = INITVALUE, dzi = INITVALUE;
  CCTK_REAL khalf = INITVALUE, kthird = INITVALUE, ktwothird = INITVALUE, kfourthird = INITVALUE, keightthird = INITVALUE;
  CCTK_REAL hdxi = INITVALUE, hdyi = INITVALUE, hdzi = INITVALUE;
  
  
  /* Declare predefined quantities */
  CCTK_REAL p1o12dx = INITVALUE;
  CCTK_REAL p1o12dy = INITVALUE;
  CCTK_REAL p1o12dz = INITVALUE;
  CCTK_REAL p1o144dxdy = INITVALUE;
  CCTK_REAL p1o144dxdz = INITVALUE;
  CCTK_REAL p1o144dydz = INITVALUE;
  CCTK_REAL p1o2dx = INITVALUE;
  CCTK_REAL p1o2dy = INITVALUE;
  CCTK_REAL p1o2dz = INITVALUE;
  CCTK_REAL pm1o12dx2 = INITVALUE;
  CCTK_REAL pm1o12dy2 = INITVALUE;
  CCTK_REAL pm1o12dz2 = INITVALUE;
  CCTK_REAL pm1o2dx = INITVALUE;
  CCTK_REAL pm1o2dy = INITVALUE;
  CCTK_REAL pm1o2dz = INITVALUE;
  
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
  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);
  dxi = 1.0 / dx;
  dyi = 1.0 / dy;
  dzi = 1.0 / dz;
  khalf = 0.5;
  kthird = 1/3.0;
  ktwothird = 2.0/3.0;
  kfourthird = 4.0/3.0;
  keightthird = 8.0/3.0;
  hdxi = 0.5 * dxi;
  hdyi = 0.5 * dyi;
  hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  p1o12dx = INV(dx)/12.;
  p1o12dy = INV(dy)/12.;
  p1o12dz = INV(dz)/12.;
  p1o144dxdy = (INV(dx)*INV(dy))/144.;
  p1o144dxdz = (INV(dx)*INV(dz))/144.;
  p1o144dydz = (INV(dy)*INV(dz))/144.;
  p1o2dx = khalf*INV(dx);
  p1o2dy = khalf*INV(dy);
  p1o2dz = khalf*INV(dz);
  pm1o12dx2 = -pow(dx,-2)/12.;
  pm1o12dy2 = -pow(dy,-2)/12.;
  pm1o12dz2 = -pow(dz,-2)/12.;
  pm1o2dx = -(khalf*INV(dx));
  pm1o2dy = -(khalf*INV(dy));
  pm1o2dz = -(khalf*INV(dz));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_MP_RHS,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lssh[CCTK_LSSH_IDX(0,0)],cctk_lssh[CCTK_LSSH_IDX(0,1)],cctk_lssh[CCTK_LSSH_IDX(0,2)])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL Atm11 = INITVALUE, Atm12 = INITVALUE, Atm13 = INITVALUE, Atm21 = INITVALUE, Atm22 = INITVALUE, Atm23 = INITVALUE;
    CCTK_REAL Atm31 = INITVALUE, Atm32 = INITVALUE, Atm33 = INITVALUE;
    CCTK_REAL Ats11 = INITVALUE, Ats12 = INITVALUE, Ats13 = INITVALUE, Ats22 = INITVALUE, Ats23 = INITVALUE, Ats33 = INITVALUE;
    CCTK_REAL Atu11 = INITVALUE, Atu21 = INITVALUE, Atu22 = INITVALUE, Atu31 = INITVALUE, Atu32 = INITVALUE, Atu33 = INITVALUE;
    CCTK_REAL detgt = INITVALUE;
    CCTK_REAL e4phi = INITVALUE;
    CCTK_REAL em4phi = INITVALUE;
    CCTK_REAL g11 = INITVALUE;
    CCTK_REAL G111 = INITVALUE, G112 = INITVALUE, G113 = INITVALUE;
    CCTK_REAL g12 = INITVALUE;
    CCTK_REAL G122 = INITVALUE, G123 = INITVALUE;
    CCTK_REAL g13 = INITVALUE;
    CCTK_REAL G133 = INITVALUE, G211 = INITVALUE, G212 = INITVALUE, G213 = INITVALUE;
    CCTK_REAL g22 = INITVALUE;
    CCTK_REAL G222 = INITVALUE, G223 = INITVALUE;
    CCTK_REAL g23 = INITVALUE;
    CCTK_REAL G233 = INITVALUE, G311 = INITVALUE, G312 = INITVALUE, G313 = INITVALUE, G322 = INITVALUE, G323 = INITVALUE;
    CCTK_REAL g33 = INITVALUE;
    CCTK_REAL G333 = INITVALUE;
    CCTK_REAL Gt111 = INITVALUE, Gt112 = INITVALUE, Gt113 = INITVALUE, Gt122 = INITVALUE, Gt123 = INITVALUE, Gt133 = INITVALUE;
    CCTK_REAL Gt211 = INITVALUE, Gt212 = INITVALUE, Gt213 = INITVALUE, Gt222 = INITVALUE, Gt223 = INITVALUE, Gt233 = INITVALUE;
    CCTK_REAL Gt311 = INITVALUE, Gt312 = INITVALUE, Gt313 = INITVALUE, Gt322 = INITVALUE, Gt323 = INITVALUE, Gt333 = INITVALUE;
    CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL R11 = INITVALUE, R12 = INITVALUE, R13 = INITVALUE, R22 = INITVALUE, R23 = INITVALUE, R33 = INITVALUE;
    CCTK_REAL Rphi11 = INITVALUE, Rphi12 = INITVALUE, Rphi13 = INITVALUE, Rphi22 = INITVALUE, Rphi23 = INITVALUE, Rphi33 = INITVALUE;
    CCTK_REAL Rt11 = INITVALUE, Rt12 = INITVALUE, Rt13 = INITVALUE, Rt22 = INITVALUE, Rt23 = INITVALUE, Rt33 = INITVALUE;
    CCTK_REAL trAts = INITVALUE;
    CCTK_REAL Xtn1 = INITVALUE, Xtn2 = INITVALUE, Xtn3 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL AL = INITVALUE;
    CCTK_REAL alphaL = INITVALUE, alpharhsL = INITVALUE;
    CCTK_REAL ArhsL = INITVALUE;
    CCTK_REAL At11L = INITVALUE, At11rhsL = INITVALUE, At12L = INITVALUE, At12rhsL = INITVALUE, At13L = INITVALUE, At13rhsL = INITVALUE;
    CCTK_REAL At22L = INITVALUE, At22rhsL = INITVALUE, At23L = INITVALUE, At23rhsL = INITVALUE, At33L = INITVALUE, At33rhsL = INITVALUE;
    CCTK_REAL B1L = INITVALUE, B1rhsL = INITVALUE, B2L = INITVALUE, B2rhsL = INITVALUE, B3L = INITVALUE, B3rhsL = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta1rhsL = INITVALUE, beta2L = INITVALUE, beta2rhsL = INITVALUE, beta3L = INITVALUE, beta3rhsL = INITVALUE;
    CCTK_REAL dJ111L = INITVALUE, dJ112L = INITVALUE, dJ113L = INITVALUE, dJ122L = INITVALUE, dJ123L = INITVALUE, dJ133L = INITVALUE;
    CCTK_REAL dJ211L = INITVALUE, dJ212L = INITVALUE, dJ213L = INITVALUE, dJ222L = INITVALUE, dJ223L = INITVALUE, dJ233L = INITVALUE;
    CCTK_REAL dJ311L = INITVALUE, dJ312L = INITVALUE, dJ313L = INITVALUE, dJ322L = INITVALUE, dJ323L = INITVALUE, dJ333L = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt11rhsL = INITVALUE, gt12L = INITVALUE, gt12rhsL = INITVALUE, gt13L = INITVALUE, gt13rhsL = INITVALUE;
    CCTK_REAL gt22L = INITVALUE, gt22rhsL = INITVALUE, gt23L = INITVALUE, gt23rhsL = INITVALUE, gt33L = INITVALUE, gt33rhsL = INITVALUE;
    CCTK_REAL J11L = INITVALUE, J12L = INITVALUE, J13L = INITVALUE, J21L = INITVALUE, J22L = INITVALUE, J23L = INITVALUE;
    CCTK_REAL J31L = INITVALUE, J32L = INITVALUE, J33L = INITVALUE;
    CCTK_REAL phiL = INITVALUE, phirhsL = INITVALUE;
    CCTK_REAL trKL = INITVALUE, trKrhsL = INITVALUE;
    CCTK_REAL Xt1L = INITVALUE, Xt1rhsL = INITVALUE, Xt2L = INITVALUE, Xt2rhsL = INITVALUE, Xt3L = INITVALUE, Xt3rhsL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandardNth1alpha = INITVALUE;
    CCTK_REAL PDstandardNth2alpha = INITVALUE;
    CCTK_REAL PDstandardNth3alpha = INITVALUE;
    CCTK_REAL PDstandardNth11alpha = INITVALUE;
    CCTK_REAL PDstandardNth22alpha = INITVALUE;
    CCTK_REAL PDstandardNth33alpha = INITVALUE;
    CCTK_REAL PDstandardNth12alpha = INITVALUE;
    CCTK_REAL PDstandardNth13alpha = INITVALUE;
    CCTK_REAL PDstandardNth23alpha = INITVALUE;
    CCTK_REAL PDstandardNth1At11 = INITVALUE;
    CCTK_REAL PDstandardNth2At11 = INITVALUE;
    CCTK_REAL PDstandardNth3At11 = INITVALUE;
    CCTK_REAL PDstandardNth1At12 = INITVALUE;
    CCTK_REAL PDstandardNth2At12 = INITVALUE;
    CCTK_REAL PDstandardNth3At12 = INITVALUE;
    CCTK_REAL PDstandardNth1At13 = INITVALUE;
    CCTK_REAL PDstandardNth2At13 = INITVALUE;
    CCTK_REAL PDstandardNth3At13 = INITVALUE;
    CCTK_REAL PDstandardNth1At22 = INITVALUE;
    CCTK_REAL PDstandardNth2At22 = INITVALUE;
    CCTK_REAL PDstandardNth3At22 = INITVALUE;
    CCTK_REAL PDstandardNth1At23 = INITVALUE;
    CCTK_REAL PDstandardNth2At23 = INITVALUE;
    CCTK_REAL PDstandardNth3At23 = INITVALUE;
    CCTK_REAL PDstandardNth1At33 = INITVALUE;
    CCTK_REAL PDstandardNth2At33 = INITVALUE;
    CCTK_REAL PDstandardNth3At33 = INITVALUE;
    CCTK_REAL PDstandardNth1B1 = INITVALUE;
    CCTK_REAL PDstandardNth2B1 = INITVALUE;
    CCTK_REAL PDstandardNth3B1 = INITVALUE;
    CCTK_REAL PDstandardNth1B2 = INITVALUE;
    CCTK_REAL PDstandardNth2B2 = INITVALUE;
    CCTK_REAL PDstandardNth3B2 = INITVALUE;
    CCTK_REAL PDstandardNth1B3 = INITVALUE;
    CCTK_REAL PDstandardNth2B3 = INITVALUE;
    CCTK_REAL PDstandardNth3B3 = INITVALUE;
    CCTK_REAL PDstandardNth1beta1 = INITVALUE;
    CCTK_REAL PDstandardNth2beta1 = INITVALUE;
    CCTK_REAL PDstandardNth3beta1 = INITVALUE;
    CCTK_REAL PDstandardNth11beta1 = INITVALUE;
    CCTK_REAL PDstandardNth22beta1 = INITVALUE;
    CCTK_REAL PDstandardNth33beta1 = INITVALUE;
    CCTK_REAL PDstandardNth12beta1 = INITVALUE;
    CCTK_REAL PDstandardNth13beta1 = INITVALUE;
    CCTK_REAL PDstandardNth23beta1 = INITVALUE;
    CCTK_REAL PDstandardNth1beta2 = INITVALUE;
    CCTK_REAL PDstandardNth2beta2 = INITVALUE;
    CCTK_REAL PDstandardNth3beta2 = INITVALUE;
    CCTK_REAL PDstandardNth11beta2 = INITVALUE;
    CCTK_REAL PDstandardNth22beta2 = INITVALUE;
    CCTK_REAL PDstandardNth33beta2 = INITVALUE;
    CCTK_REAL PDstandardNth12beta2 = INITVALUE;
    CCTK_REAL PDstandardNth13beta2 = INITVALUE;
    CCTK_REAL PDstandardNth23beta2 = INITVALUE;
    CCTK_REAL PDstandardNth1beta3 = INITVALUE;
    CCTK_REAL PDstandardNth2beta3 = INITVALUE;
    CCTK_REAL PDstandardNth3beta3 = INITVALUE;
    CCTK_REAL PDstandardNth11beta3 = INITVALUE;
    CCTK_REAL PDstandardNth22beta3 = INITVALUE;
    CCTK_REAL PDstandardNth33beta3 = INITVALUE;
    CCTK_REAL PDstandardNth12beta3 = INITVALUE;
    CCTK_REAL PDstandardNth13beta3 = INITVALUE;
    CCTK_REAL PDstandardNth23beta3 = INITVALUE;
    CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    CCTK_REAL PDstandardNth11gt11 = INITVALUE;
    CCTK_REAL PDstandardNth22gt11 = INITVALUE;
    CCTK_REAL PDstandardNth33gt11 = INITVALUE;
    CCTK_REAL PDstandardNth12gt11 = INITVALUE;
    CCTK_REAL PDstandardNth13gt11 = INITVALUE;
    CCTK_REAL PDstandardNth23gt11 = INITVALUE;
    CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    CCTK_REAL PDstandardNth11gt12 = INITVALUE;
    CCTK_REAL PDstandardNth22gt12 = INITVALUE;
    CCTK_REAL PDstandardNth33gt12 = INITVALUE;
    CCTK_REAL PDstandardNth12gt12 = INITVALUE;
    CCTK_REAL PDstandardNth13gt12 = INITVALUE;
    CCTK_REAL PDstandardNth23gt12 = INITVALUE;
    CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    CCTK_REAL PDstandardNth11gt13 = INITVALUE;
    CCTK_REAL PDstandardNth22gt13 = INITVALUE;
    CCTK_REAL PDstandardNth33gt13 = INITVALUE;
    CCTK_REAL PDstandardNth12gt13 = INITVALUE;
    CCTK_REAL PDstandardNth13gt13 = INITVALUE;
    CCTK_REAL PDstandardNth23gt13 = INITVALUE;
    CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    CCTK_REAL PDstandardNth11gt22 = INITVALUE;
    CCTK_REAL PDstandardNth22gt22 = INITVALUE;
    CCTK_REAL PDstandardNth33gt22 = INITVALUE;
    CCTK_REAL PDstandardNth12gt22 = INITVALUE;
    CCTK_REAL PDstandardNth13gt22 = INITVALUE;
    CCTK_REAL PDstandardNth23gt22 = INITVALUE;
    CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    CCTK_REAL PDstandardNth11gt23 = INITVALUE;
    CCTK_REAL PDstandardNth22gt23 = INITVALUE;
    CCTK_REAL PDstandardNth33gt23 = INITVALUE;
    CCTK_REAL PDstandardNth12gt23 = INITVALUE;
    CCTK_REAL PDstandardNth13gt23 = INITVALUE;
    CCTK_REAL PDstandardNth23gt23 = INITVALUE;
    CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    CCTK_REAL PDstandardNth11gt33 = INITVALUE;
    CCTK_REAL PDstandardNth22gt33 = INITVALUE;
    CCTK_REAL PDstandardNth33gt33 = INITVALUE;
    CCTK_REAL PDstandardNth12gt33 = INITVALUE;
    CCTK_REAL PDstandardNth13gt33 = INITVALUE;
    CCTK_REAL PDstandardNth23gt33 = INITVALUE;
    CCTK_REAL PDstandardNth1phi = INITVALUE;
    CCTK_REAL PDstandardNth2phi = INITVALUE;
    CCTK_REAL PDstandardNth3phi = INITVALUE;
    CCTK_REAL PDstandardNth11phi = INITVALUE;
    CCTK_REAL PDstandardNth22phi = INITVALUE;
    CCTK_REAL PDstandardNth33phi = INITVALUE;
    CCTK_REAL PDstandardNth12phi = INITVALUE;
    CCTK_REAL PDstandardNth13phi = INITVALUE;
    CCTK_REAL PDstandardNth23phi = INITVALUE;
    CCTK_REAL PDstandardNth1trK = INITVALUE;
    CCTK_REAL PDstandardNth2trK = INITVALUE;
    CCTK_REAL PDstandardNth3trK = INITVALUE;
    CCTK_REAL PDstandardNth1Xt1 = INITVALUE;
    CCTK_REAL PDstandardNth2Xt1 = INITVALUE;
    CCTK_REAL PDstandardNth3Xt1 = INITVALUE;
    CCTK_REAL PDstandardNth1Xt2 = INITVALUE;
    CCTK_REAL PDstandardNth2Xt2 = INITVALUE;
    CCTK_REAL PDstandardNth3Xt2 = INITVALUE;
    CCTK_REAL PDstandardNth1Xt3 = INITVALUE;
    CCTK_REAL PDstandardNth2Xt3 = INITVALUE;
    CCTK_REAL PDstandardNth3Xt3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    AL = A[index];
    alphaL = alpha[index];
    At11L = At11[index];
    At12L = At12[index];
    At13L = At13[index];
    At22L = At22[index];
    At23L = At23[index];
    At33L = At33[index];
    B1L = B1[index];
    B2L = B2[index];
    B3L = B3[index];
    beta1L = beta1[index];
    beta2L = beta2[index];
    beta3L = beta3[index];
    dJ111L = dJ111[index];
    dJ112L = dJ112[index];
    dJ113L = dJ113[index];
    dJ122L = dJ122[index];
    dJ123L = dJ123[index];
    dJ133L = dJ133[index];
    dJ211L = dJ211[index];
    dJ212L = dJ212[index];
    dJ213L = dJ213[index];
    dJ222L = dJ222[index];
    dJ223L = dJ223[index];
    dJ233L = dJ233[index];
    dJ311L = dJ311[index];
    dJ312L = dJ312[index];
    dJ313L = dJ313[index];
    dJ322L = dJ322[index];
    dJ323L = dJ323[index];
    dJ333L = dJ333[index];
    gt11L = gt11[index];
    gt12L = gt12[index];
    gt13L = gt13[index];
    gt22L = gt22[index];
    gt23L = gt23[index];
    gt33L = gt33[index];
    J11L = J11[index];
    J12L = J12[index];
    J13L = J13[index];
    J21L = J21[index];
    J22L = J22[index];
    J23L = J23[index];
    J31L = J31[index];
    J32L = J32[index];
    J33L = J33[index];
    phiL = phi[index];
    trKL = trK[index];
    trKrhsL = trKrhs[index];
    Xt1L = Xt1[index];
    Xt1rhsL = Xt1rhs[index];
    Xt2L = Xt2[index];
    Xt2rhsL = Xt2rhs[index];
    Xt3L = Xt3[index];
    Xt3rhsL = Xt3rhs[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandardNth1alpha = PDstandardNth1(alpha, i, j, k);
    PDstandardNth2alpha = PDstandardNth2(alpha, i, j, k);
    PDstandardNth3alpha = PDstandardNth3(alpha, i, j, k);
    PDstandardNth11alpha = PDstandardNth11(alpha, i, j, k);
    PDstandardNth22alpha = PDstandardNth22(alpha, i, j, k);
    PDstandardNth33alpha = PDstandardNth33(alpha, i, j, k);
    PDstandardNth12alpha = PDstandardNth12(alpha, i, j, k);
    PDstandardNth13alpha = PDstandardNth13(alpha, i, j, k);
    PDstandardNth23alpha = PDstandardNth23(alpha, i, j, k);
    PDstandardNth1At11 = PDstandardNth1(At11, i, j, k);
    PDstandardNth2At11 = PDstandardNth2(At11, i, j, k);
    PDstandardNth3At11 = PDstandardNth3(At11, i, j, k);
    PDstandardNth1At12 = PDstandardNth1(At12, i, j, k);
    PDstandardNth2At12 = PDstandardNth2(At12, i, j, k);
    PDstandardNth3At12 = PDstandardNth3(At12, i, j, k);
    PDstandardNth1At13 = PDstandardNth1(At13, i, j, k);
    PDstandardNth2At13 = PDstandardNth2(At13, i, j, k);
    PDstandardNth3At13 = PDstandardNth3(At13, i, j, k);
    PDstandardNth1At22 = PDstandardNth1(At22, i, j, k);
    PDstandardNth2At22 = PDstandardNth2(At22, i, j, k);
    PDstandardNth3At22 = PDstandardNth3(At22, i, j, k);
    PDstandardNth1At23 = PDstandardNth1(At23, i, j, k);
    PDstandardNth2At23 = PDstandardNth2(At23, i, j, k);
    PDstandardNth3At23 = PDstandardNth3(At23, i, j, k);
    PDstandardNth1At33 = PDstandardNth1(At33, i, j, k);
    PDstandardNth2At33 = PDstandardNth2(At33, i, j, k);
    PDstandardNth3At33 = PDstandardNth3(At33, i, j, k);
    PDstandardNth1B1 = PDstandardNth1(B1, i, j, k);
    PDstandardNth2B1 = PDstandardNth2(B1, i, j, k);
    PDstandardNth3B1 = PDstandardNth3(B1, i, j, k);
    PDstandardNth1B2 = PDstandardNth1(B2, i, j, k);
    PDstandardNth2B2 = PDstandardNth2(B2, i, j, k);
    PDstandardNth3B2 = PDstandardNth3(B2, i, j, k);
    PDstandardNth1B3 = PDstandardNth1(B3, i, j, k);
    PDstandardNth2B3 = PDstandardNth2(B3, i, j, k);
    PDstandardNth3B3 = PDstandardNth3(B3, i, j, k);
    PDstandardNth1beta1 = PDstandardNth1(beta1, i, j, k);
    PDstandardNth2beta1 = PDstandardNth2(beta1, i, j, k);
    PDstandardNth3beta1 = PDstandardNth3(beta1, i, j, k);
    PDstandardNth11beta1 = PDstandardNth11(beta1, i, j, k);
    PDstandardNth22beta1 = PDstandardNth22(beta1, i, j, k);
    PDstandardNth33beta1 = PDstandardNth33(beta1, i, j, k);
    PDstandardNth12beta1 = PDstandardNth12(beta1, i, j, k);
    PDstandardNth13beta1 = PDstandardNth13(beta1, i, j, k);
    PDstandardNth23beta1 = PDstandardNth23(beta1, i, j, k);
    PDstandardNth1beta2 = PDstandardNth1(beta2, i, j, k);
    PDstandardNth2beta2 = PDstandardNth2(beta2, i, j, k);
    PDstandardNth3beta2 = PDstandardNth3(beta2, i, j, k);
    PDstandardNth11beta2 = PDstandardNth11(beta2, i, j, k);
    PDstandardNth22beta2 = PDstandardNth22(beta2, i, j, k);
    PDstandardNth33beta2 = PDstandardNth33(beta2, i, j, k);
    PDstandardNth12beta2 = PDstandardNth12(beta2, i, j, k);
    PDstandardNth13beta2 = PDstandardNth13(beta2, i, j, k);
    PDstandardNth23beta2 = PDstandardNth23(beta2, i, j, k);
    PDstandardNth1beta3 = PDstandardNth1(beta3, i, j, k);
    PDstandardNth2beta3 = PDstandardNth2(beta3, i, j, k);
    PDstandardNth3beta3 = PDstandardNth3(beta3, i, j, k);
    PDstandardNth11beta3 = PDstandardNth11(beta3, i, j, k);
    PDstandardNth22beta3 = PDstandardNth22(beta3, i, j, k);
    PDstandardNth33beta3 = PDstandardNth33(beta3, i, j, k);
    PDstandardNth12beta3 = PDstandardNth12(beta3, i, j, k);
    PDstandardNth13beta3 = PDstandardNth13(beta3, i, j, k);
    PDstandardNth23beta3 = PDstandardNth23(beta3, i, j, k);
    PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    PDstandardNth11gt11 = PDstandardNth11(gt11, i, j, k);
    PDstandardNth22gt11 = PDstandardNth22(gt11, i, j, k);
    PDstandardNth33gt11 = PDstandardNth33(gt11, i, j, k);
    PDstandardNth12gt11 = PDstandardNth12(gt11, i, j, k);
    PDstandardNth13gt11 = PDstandardNth13(gt11, i, j, k);
    PDstandardNth23gt11 = PDstandardNth23(gt11, i, j, k);
    PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    PDstandardNth11gt12 = PDstandardNth11(gt12, i, j, k);
    PDstandardNth22gt12 = PDstandardNth22(gt12, i, j, k);
    PDstandardNth33gt12 = PDstandardNth33(gt12, i, j, k);
    PDstandardNth12gt12 = PDstandardNth12(gt12, i, j, k);
    PDstandardNth13gt12 = PDstandardNth13(gt12, i, j, k);
    PDstandardNth23gt12 = PDstandardNth23(gt12, i, j, k);
    PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    PDstandardNth11gt13 = PDstandardNth11(gt13, i, j, k);
    PDstandardNth22gt13 = PDstandardNth22(gt13, i, j, k);
    PDstandardNth33gt13 = PDstandardNth33(gt13, i, j, k);
    PDstandardNth12gt13 = PDstandardNth12(gt13, i, j, k);
    PDstandardNth13gt13 = PDstandardNth13(gt13, i, j, k);
    PDstandardNth23gt13 = PDstandardNth23(gt13, i, j, k);
    PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    PDstandardNth11gt22 = PDstandardNth11(gt22, i, j, k);
    PDstandardNth22gt22 = PDstandardNth22(gt22, i, j, k);
    PDstandardNth33gt22 = PDstandardNth33(gt22, i, j, k);
    PDstandardNth12gt22 = PDstandardNth12(gt22, i, j, k);
    PDstandardNth13gt22 = PDstandardNth13(gt22, i, j, k);
    PDstandardNth23gt22 = PDstandardNth23(gt22, i, j, k);
    PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    PDstandardNth11gt23 = PDstandardNth11(gt23, i, j, k);
    PDstandardNth22gt23 = PDstandardNth22(gt23, i, j, k);
    PDstandardNth33gt23 = PDstandardNth33(gt23, i, j, k);
    PDstandardNth12gt23 = PDstandardNth12(gt23, i, j, k);
    PDstandardNth13gt23 = PDstandardNth13(gt23, i, j, k);
    PDstandardNth23gt23 = PDstandardNth23(gt23, i, j, k);
    PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    PDstandardNth11gt33 = PDstandardNth11(gt33, i, j, k);
    PDstandardNth22gt33 = PDstandardNth22(gt33, i, j, k);
    PDstandardNth33gt33 = PDstandardNth33(gt33, i, j, k);
    PDstandardNth12gt33 = PDstandardNth12(gt33, i, j, k);
    PDstandardNth13gt33 = PDstandardNth13(gt33, i, j, k);
    PDstandardNth23gt33 = PDstandardNth23(gt33, i, j, k);
    PDstandardNth1phi = PDstandardNth1(phi, i, j, k);
    PDstandardNth2phi = PDstandardNth2(phi, i, j, k);
    PDstandardNth3phi = PDstandardNth3(phi, i, j, k);
    PDstandardNth11phi = PDstandardNth11(phi, i, j, k);
    PDstandardNth22phi = PDstandardNth22(phi, i, j, k);
    PDstandardNth33phi = PDstandardNth33(phi, i, j, k);
    PDstandardNth12phi = PDstandardNth12(phi, i, j, k);
    PDstandardNth13phi = PDstandardNth13(phi, i, j, k);
    PDstandardNth23phi = PDstandardNth23(phi, i, j, k);
    PDstandardNth1trK = PDstandardNth1(trK, i, j, k);
    PDstandardNth2trK = PDstandardNth2(trK, i, j, k);
    PDstandardNth3trK = PDstandardNth3(trK, i, j, k);
    PDstandardNth1Xt1 = PDstandardNth1(Xt1, i, j, k);
    PDstandardNth2Xt1 = PDstandardNth2(Xt1, i, j, k);
    PDstandardNth3Xt1 = PDstandardNth3(Xt1, i, j, k);
    PDstandardNth1Xt2 = PDstandardNth1(Xt2, i, j, k);
    PDstandardNth2Xt2 = PDstandardNth2(Xt2, i, j, k);
    PDstandardNth3Xt2 = PDstandardNth3(Xt2, i, j, k);
    PDstandardNth1Xt3 = PDstandardNth1(Xt3, i, j, k);
    PDstandardNth2Xt3 = PDstandardNth2(Xt3, i, j, k);
    PDstandardNth3Xt3 = PDstandardNth3(Xt3, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    detgt  =  1;
    
    gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    gtu21  =  (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    gtu32  =  (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    Gt111  =  khalf*((gtu11*J11L - gtu21*J12L - gtu31*J13L)*PDstandardNth1gt11 + 
          (gtu11*J21L - gtu21*J22L - gtu31*J23L)*PDstandardNth2gt11 + 
          (gtu11*J31L - gtu21*J32L - gtu31*J33L)*PDstandardNth3gt11 + 
          2*(J11L*(gtu21*PDstandardNth1gt12 + gtu31*PDstandardNth1gt13) + 
             J21L*(gtu21*PDstandardNth2gt12 + gtu31*PDstandardNth2gt13) + 
             J31L*(gtu21*PDstandardNth3gt12 + gtu31*PDstandardNth3gt13)));
    
    Gt211  =  khalf*((gtu21*J11L - gtu22*J12L - gtu32*J13L)*PDstandardNth1gt11 + 
          (gtu21*J21L - gtu22*J22L - gtu32*J23L)*PDstandardNth2gt11 + 
          (gtu21*J31L - gtu22*J32L - gtu32*J33L)*PDstandardNth3gt11 + 
          2*(J11L*(gtu22*PDstandardNth1gt12 + gtu32*PDstandardNth1gt13) + 
             J21L*(gtu22*PDstandardNth2gt12 + gtu32*PDstandardNth2gt13) + 
             J31L*(gtu22*PDstandardNth3gt12 + gtu32*PDstandardNth3gt13)));
    
    Gt311  =  khalf*((gtu31*J11L - gtu32*J12L - gtu33*J13L)*PDstandardNth1gt11 + 
          (gtu31*J21L - gtu32*J22L - gtu33*J23L)*PDstandardNth2gt11 + 
          (gtu31*J31L - gtu32*J32L - gtu33*J33L)*PDstandardNth3gt11 + 
          2*(J11L*(gtu32*PDstandardNth1gt12 + gtu33*PDstandardNth1gt13) + 
             J21L*(gtu32*PDstandardNth2gt12 + gtu33*PDstandardNth2gt13) + 
             J31L*(gtu32*PDstandardNth3gt12 + gtu33*PDstandardNth3gt13)));
    
    Gt112  =  khalf*(gtu11*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + 
             J32L*PDstandardNth3gt11) + 
          gtu21*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
             J31L*PDstandardNth3gt22) + 
          gtu31*(-(J13L*PDstandardNth1gt12) + J12L*PDstandardNth1gt13 + 
             J11L*PDstandardNth1gt23 - J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 - 
             J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23));
    
    Gt212  =  khalf*(gtu21*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + 
             J32L*PDstandardNth3gt11) + 
          gtu22*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
             J31L*PDstandardNth3gt22) + 
          gtu32*(-(J13L*PDstandardNth1gt12) + J12L*PDstandardNth1gt13 + 
             J11L*PDstandardNth1gt23 - J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 - 
             J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23));
    
    Gt312  =  khalf*(gtu31*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + 
             J32L*PDstandardNth3gt11) + 
          gtu32*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
             J31L*PDstandardNth3gt22) + 
          gtu33*(-(J13L*PDstandardNth1gt12) + J12L*PDstandardNth1gt13 + 
             J11L*PDstandardNth1gt23 - J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 - 
             J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23));
    
    Gt113  =  khalf*(gtu11*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + 
             J33L*PDstandardNth3gt11) + 
          gtu21*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
             J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
             J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23) + 
          gtu31*(J11L*PDstandardNth1gt33 + J21L*PDstandardNth2gt33 + 
             J31L*PDstandardNth3gt33));
    
    Gt213  =  khalf*(gtu21*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + 
             J33L*PDstandardNth3gt11) + 
          gtu22*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
             J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
             J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23) + 
          gtu32*(J11L*PDstandardNth1gt33 + J21L*PDstandardNth2gt33 + 
             J31L*PDstandardNth3gt33));
    
    Gt313  =  khalf*(gtu31*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + 
             J33L*PDstandardNth3gt11) + 
          gtu32*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
             J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
             J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23) + 
          gtu33*(J11L*PDstandardNth1gt33 + J21L*PDstandardNth2gt33 + 
             J31L*PDstandardNth3gt33));
    
    Gt122  =  khalf*(gtu11*(-(J11L*PDstandardNth1gt22) + 
             2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
             J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
             J31L*PDstandardNth3gt22) + 
          gtu21*(J12L*PDstandardNth1gt22 + J22L*PDstandardNth2gt22 + 
             J32L*PDstandardNth3gt22) - 
          gtu31*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
             J33L*PDstandardNth3gt22 - 
             2*(J12L*PDstandardNth1gt23 + J22L*PDstandardNth2gt23 + 
                J32L*PDstandardNth3gt23)));
    
    Gt222  =  khalf*(gtu21*(-(J11L*PDstandardNth1gt22) + 
             2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
             J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
             J31L*PDstandardNth3gt22) + 
          gtu22*(J12L*PDstandardNth1gt22 + J22L*PDstandardNth2gt22 + 
             J32L*PDstandardNth3gt22) - 
          gtu32*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
             J33L*PDstandardNth3gt22 - 
             2*(J12L*PDstandardNth1gt23 + J22L*PDstandardNth2gt23 + 
                J32L*PDstandardNth3gt23)));
    
    Gt322  =  khalf*(gtu31*(-(J11L*PDstandardNth1gt22) + 
             2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
             J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
             J31L*PDstandardNth3gt22) + 
          gtu32*(J12L*PDstandardNth1gt22 + J22L*PDstandardNth2gt22 + 
             J32L*PDstandardNth3gt22) - 
          gtu33*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
             J33L*PDstandardNth3gt22 - 
             2*(J12L*PDstandardNth1gt23 + J22L*PDstandardNth2gt23 + 
                J32L*PDstandardNth3gt23)));
    
    Gt123  =  khalf*(gtu21*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
             J33L*PDstandardNth3gt22) + 
          gtu11*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
             J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
             J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt23) + 
          gtu31*(J12L*PDstandardNth1gt33 + J22L*PDstandardNth2gt33 + 
             J32L*PDstandardNth3gt33));
    
    Gt223  =  khalf*(gtu22*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
             J33L*PDstandardNth3gt22) + 
          gtu21*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
             J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
             J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt23) + 
          gtu32*(J12L*PDstandardNth1gt33 + J22L*PDstandardNth2gt33 + 
             J32L*PDstandardNth3gt33));
    
    Gt323  =  khalf*(gtu32*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
             J33L*PDstandardNth3gt22) + 
          gtu31*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
             J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
             J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt23) + 
          gtu33*(J12L*PDstandardNth1gt33 + J22L*PDstandardNth2gt33 + 
             J32L*PDstandardNth3gt33));
    
    Gt133  =  khalf*(gtu11*(-(J11L*PDstandardNth1gt33) + 
             2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
             J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt33) + 
          gtu21*(-(J12L*PDstandardNth1gt33) + 
             2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
             J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
             J32L*PDstandardNth3gt33) + 
          gtu31*(J13L*PDstandardNth1gt33 + J23L*PDstandardNth2gt33 + 
             J33L*PDstandardNth3gt33));
    
    Gt233  =  khalf*(gtu21*(-(J11L*PDstandardNth1gt33) + 
             2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
             J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt33) + 
          gtu22*(-(J12L*PDstandardNth1gt33) + 
             2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
             J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
             J32L*PDstandardNth3gt33) + 
          gtu32*(J13L*PDstandardNth1gt33 + J23L*PDstandardNth2gt33 + 
             J33L*PDstandardNth3gt33));
    
    Gt333  =  khalf*(gtu31*(-(J11L*PDstandardNth1gt33) + 
             2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
             J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt33) + 
          gtu32*(-(J12L*PDstandardNth1gt33) + 
             2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
             J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
             J32L*PDstandardNth3gt33) + 
          gtu33*(J13L*PDstandardNth1gt33 + J23L*PDstandardNth2gt33 + 
             J33L*PDstandardNth3gt33));
    
    Xtn1  =  Gt111*gtu11 + Gt122*gtu22 + 
        2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    Xtn2  =  Gt211*gtu11 + Gt222*gtu22 + 
        2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    Xtn3  =  Gt311*gtu11 + Gt322*gtu22 + 
        2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    Rt11  =  (Gt113*(gt13L*Gt312 + 3*(gt12L*Gt212 + gt13L*Gt312)) + 
           gt12L*(Gt213*(4*Gt112 + 2*Gt222) + Gt212*(Gt113 + 2*Gt223) + 
              2*(Gt233*Gt312 + Gt223*Gt313)) + 
           gt11L*(6*Gt112*Gt113 + 2*(Gt122*Gt213 + Gt133*Gt312 + 
                 Gt123*(Gt212 + Gt313))) + 
           gt13L*(2*Gt213*Gt322 + Gt313*(4*Gt112 + 2*Gt323)) + 
           2*(Gt213*(Gt212*gt22L + gt23L*Gt312) + 
              Gt212*(gt23L*Gt313 + gt13L*Gt323) + Gt312*(gt13L*Gt333 + Gt313*gt33L))
           )*gtu32 + J11L*(gt11L*PDstandardNth1Xt1 + gt12L*PDstandardNth1Xt2 + 
           gt13L*PDstandardNth1Xt3) + 
        J21L*(gt11L*PDstandardNth2Xt1 + gt12L*PDstandardNth2Xt2 + 
           gt13L*PDstandardNth2Xt3) + 
        J31L*(gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + 
           gt13L*PDstandardNth3Xt3) + 
        (Gt111*gt11L + gt12L*Gt211 + gt13L*Gt311)*Xtn1 + 
        (Gt112*gt11L + gt12L*Gt212 + gt13L*Gt312)*Xtn2 + 
        (Gt113*gt11L + gt12L*Gt213 + gt13L*Gt313)*Xtn3 + 
        gtu21*(Gt112*(gt13L*Gt311 + 3*(gt12L*Gt211 + gt13L*Gt311)) + 
           gt11L*(Gt112*(6*Gt111 + 2*Gt212) + 
              2*(Gt122*Gt211 + Gt123*Gt311 + Gt113*Gt312)) + 
           2*(Gt212*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + 
              gt13L*(Gt211*Gt322 + Gt311*Gt323)) + 
           Gt312*(gt13L*(4*Gt111 + 2*Gt313) + 2*(Gt211*gt23L + Gt311*gt33L)) + 
           gt12L*(4*Gt111*Gt212 + Gt211*(Gt112 + 2*Gt222) + 
              2*(Gt223*Gt311 + Gt213*Gt312 + SQR(Gt212)))) + 
        gtu11*(4*Gt111*(gt12L*Gt211 + gt13L*Gt311) + 
           2*(gt12L*(Gt211*Gt212 + Gt213*Gt311) + 
              Gt211*(gt23L*Gt311 + gt13L*Gt312) + gt13L*Gt311*Gt313) + 
           gt11L*(2*(Gt112*Gt211 + Gt113*Gt311) + 3*SQR(Gt111)) + 
           gt22L*SQR(Gt211) + gt33L*SQR(Gt311)) + 
        gtu22*(4*Gt112*(gt12L*Gt212 + gt13L*Gt312) + 
           2*(gt12L*(Gt212*Gt222 + Gt223*Gt312) + 
              Gt212*(gt23L*Gt312 + gt13L*Gt322) + gt13L*Gt312*Gt323) + 
           gt11L*(2*(Gt122*Gt212 + Gt123*Gt312) + 3*SQR(Gt112)) + 
           gt22L*SQR(Gt212) + gt33L*SQR(Gt312)) + 
        gtu33*(4*Gt113*(gt12L*Gt213 + gt13L*Gt313) + 
           2*(gt12L*(Gt213*Gt223 + Gt233*Gt313) + 
              Gt213*(gt23L*Gt313 + gt13L*Gt323) + gt13L*Gt313*Gt333) + 
           gt11L*(2*(Gt123*Gt213 + Gt133*Gt313) + 3*SQR(Gt113)) + 
           gt22L*SQR(Gt213) + gt33L*SQR(Gt313)) + 
        gtu31*(Gt113*(gt13L*Gt311 + 3*(gt12L*Gt211 + gt13L*Gt311)) + 
           gt11L*(2*(Gt123*Gt211 + Gt112*Gt213 + Gt133*Gt311) + 
              Gt113*(6*Gt111 + 2*Gt313)) + 
           gt12L*(Gt211*(Gt113 + 2*Gt223) + 2*Gt233*Gt311 + 
              Gt213*(4*Gt111 + 2*(Gt212 + Gt313))) + 
           2*(Gt213*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + gt13L*Gt311*Gt333 + 
              Gt313*(Gt211*gt23L + Gt311*gt33L)) + 
           gt13L*(4*Gt111*Gt313 + 2*(Gt211*Gt323 + SQR(Gt313)))) + 
        khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt11 + 
                 J12L*J21L*PDstandardNth12gt11 + J11L*J22L*PDstandardNth12gt11 + 
                 J12L*J31L*PDstandardNth13gt11 + J11L*J32L*PDstandardNth13gt11 + 
                 dJ112L*PDstandardNth1gt11 + J21L*J22L*PDstandardNth22gt11 + 
                 J22L*J31L*PDstandardNth23gt11 + J21L*J32L*PDstandardNth23gt11 + 
                 dJ212L*PDstandardNth2gt11 + J31L*J32L*PDstandardNth33gt11 + 
                 dJ312L*PDstandardNth3gt11) + 
              gtu31*(J11L*J13L*PDstandardNth11gt11 + 
                 J13L*J21L*PDstandardNth12gt11 + J11L*J23L*PDstandardNth12gt11 + 
                 J13L*J31L*PDstandardNth13gt11 + J11L*J33L*PDstandardNth13gt11 + 
                 dJ113L*PDstandardNth1gt11 + J21L*J23L*PDstandardNth22gt11 + 
                 J23L*J31L*PDstandardNth23gt11 + J21L*J33L*PDstandardNth23gt11 + 
                 dJ213L*PDstandardNth2gt11 + J31L*J33L*PDstandardNth33gt11 + 
                 dJ313L*PDstandardNth3gt11) + 
              gtu32*(J12L*J13L*PDstandardNth11gt11 + 
                 J13L*J22L*PDstandardNth12gt11 + J12L*J23L*PDstandardNth12gt11 + 
                 J13L*J32L*PDstandardNth13gt11 + J12L*J33L*PDstandardNth13gt11 + 
                 dJ123L*PDstandardNth1gt11 + J22L*J23L*PDstandardNth22gt11 + 
                 J23L*J32L*PDstandardNth23gt11 + J22L*J33L*PDstandardNth23gt11 + 
                 dJ223L*PDstandardNth2gt11 + J32L*J33L*PDstandardNth33gt11 + 
                 dJ323L*PDstandardNth3gt11)) - 
           gtu11*(2*J11L*J21L*PDstandardNth12gt11 + 
              2*J11L*J31L*PDstandardNth13gt11 + dJ111L*PDstandardNth1gt11 + 
              2*J21L*J31L*PDstandardNth23gt11 + dJ211L*PDstandardNth2gt11 + 
              dJ311L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J11L) + 
              PDstandardNth22gt11*SQR(J21L) + PDstandardNth33gt11*SQR(J31L)) - 
           gtu22*(2*J12L*J22L*PDstandardNth12gt11 + 
              2*J12L*J32L*PDstandardNth13gt11 + dJ122L*PDstandardNth1gt11 + 
              2*J22L*J32L*PDstandardNth23gt11 + dJ222L*PDstandardNth2gt11 + 
              dJ322L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J12L) + 
              PDstandardNth22gt11*SQR(J22L) + PDstandardNth33gt11*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt11 + 
              2*J13L*J33L*PDstandardNth13gt11 + dJ133L*PDstandardNth1gt11 + 
              2*J23L*J33L*PDstandardNth23gt11 + dJ233L*PDstandardNth2gt11 + 
              dJ333L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J13L) + 
              PDstandardNth22gt11*SQR(J23L) + PDstandardNth33gt11*SQR(J33L)));
    
    Rt12  =  khalf*((gt12L*J11L + gt11L*J12L)*PDstandardNth1Xt1 + 
          (gt22L*J11L + gt12L*J12L)*PDstandardNth1Xt2 + 
          (gt23L*J11L + gt13L*J12L)*PDstandardNth1Xt3 + 
          (gt12L*J21L + gt11L*J22L)*PDstandardNth2Xt1 + 
          (gt22L*J21L + gt12L*J22L)*PDstandardNth2Xt2 + 
          (gt23L*J21L + gt13L*J22L)*PDstandardNth2Xt3 - 
          2*(gtu21*(J11L*J12L*PDstandardNth11gt12 + J12L*J21L*PDstandardNth12gt12 + 
                J11L*J22L*PDstandardNth12gt12 + J12L*J31L*PDstandardNth13gt12 + 
                J11L*J32L*PDstandardNth13gt12 + dJ112L*PDstandardNth1gt12 + 
                J21L*J22L*PDstandardNth22gt12 + J22L*J31L*PDstandardNth23gt12 + 
                J21L*J32L*PDstandardNth23gt12 + dJ212L*PDstandardNth2gt12 + 
                J31L*J32L*PDstandardNth33gt12 + dJ312L*PDstandardNth3gt12) + 
             gtu31*(J11L*J13L*PDstandardNth11gt12 + J13L*J21L*PDstandardNth12gt12 + 
                J11L*J23L*PDstandardNth12gt12 + J13L*J31L*PDstandardNth13gt12 + 
                J11L*J33L*PDstandardNth13gt12 + dJ113L*PDstandardNth1gt12 + 
                J21L*J23L*PDstandardNth22gt12 + J23L*J31L*PDstandardNth23gt12 + 
                J21L*J33L*PDstandardNth23gt12 + dJ213L*PDstandardNth2gt12 + 
                J31L*J33L*PDstandardNth33gt12 + dJ313L*PDstandardNth3gt12) + 
             gtu32*(J12L*J13L*PDstandardNth11gt12 + J13L*J22L*PDstandardNth12gt12 + 
                J12L*J23L*PDstandardNth12gt12 + J13L*J32L*PDstandardNth13gt12 + 
                J12L*J33L*PDstandardNth13gt12 + dJ123L*PDstandardNth1gt12 + 
                J22L*J23L*PDstandardNth22gt12 + J23L*J32L*PDstandardNth23gt12 + 
                J22L*J33L*PDstandardNth23gt12 + dJ223L*PDstandardNth2gt12 + 
                J32L*J33L*PDstandardNth33gt12 + dJ323L*PDstandardNth3gt12)) + 
          (gt12L*J31L + gt11L*J32L)*PDstandardNth3Xt1 + 
          (gt22L*J31L + gt12L*J32L)*PDstandardNth3Xt2 + 
          (gt23L*J31L + gt13L*J32L)*PDstandardNth3Xt3 + 
          (Gt112*gt11L + Gt111*gt12L + gt12L*Gt212 + Gt211*gt22L + gt23L*Gt311 + 
             gt13L*Gt312)*Xtn1 + (gt11L*Gt122 + Gt112*gt12L + gt12L*Gt222 + 
             Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322)*Xtn2 + 
          (gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + 
             gt13L*Gt323)*Xtn3 + 2*((Gt123*gt12L*Gt211 + Gt113*gt12L*Gt212 + 
                2*Gt112*gt12L*Gt213 + gt12L*Gt212*Gt223 + Gt212*Gt213*gt22L + 
                Gt211*Gt223*gt22L + gt12L*Gt133*Gt311 + gt22L*Gt233*Gt311 + 
                Gt113*gt13L*Gt312 + gt12L*Gt233*Gt312 + Gt213*gt23L*Gt312 + 
                gt11L*(2*Gt112*Gt113 + Gt123*Gt212 + Gt133*Gt312) + 
                2*Gt112*gt13L*Gt313 + Gt212*gt23L*Gt313 + 
                Gt111*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + 
                gt13L*Gt212*Gt323 + Gt211*gt23L*Gt323 + gt23L*Gt311*Gt333 + 
                gt13L*Gt312*Gt333 + Gt312*Gt313*gt33L)*gtu31 + 
             (Gt123*gt12L*Gt212 + 2*Gt122*gt12L*Gt213 + Gt113*gt12L*Gt222 + 
                gt12L*Gt222*Gt223 + Gt213*Gt222*gt22L + Gt212*Gt223*gt22L + 
                gt12L*Gt133*Gt312 + gt22L*Gt233*Gt312 + 2*Gt122*gt13L*Gt313 + 
                Gt222*gt23L*Gt313 + Gt112*
                 (Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + Gt113*gt13L*Gt322 + 
                gt12L*Gt233*Gt322 + Gt213*gt23L*Gt322 + 
                gt11L*(2*Gt113*Gt122 + Gt123*Gt222 + Gt133*Gt322) + 
                gt13L*Gt222*Gt323 + Gt212*gt23L*Gt323 + gt23L*Gt312*Gt333 + 
                gt13L*Gt322*Gt333 + Gt313*Gt322*gt33L)*gtu32 + 
             gtu11*(3*Gt112*gt12L*Gt211 + 2*Gt211*Gt212*gt22L + Gt113*gt12L*Gt311 + 
                2*Gt112*gt13L*Gt311 + Gt213*gt22L*Gt311 + Gt212*gt23L*Gt311 + 
                gt13L*Gt212*Gt312 + gt12L*Gt213*Gt312 + 2*Gt211*gt23L*Gt312 + 
                gt11L*(2*Gt111*Gt112 + Gt112*Gt212 + Gt113*Gt312) + 
                Gt111*(gt12L*Gt212 + Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + 
                gt23L*Gt311*Gt313 + gt13L*Gt312*Gt313 + Gt311*Gt312*gt33L + 
                gt12L*SQR(Gt111) + gt12L*SQR(Gt212)) + 
             gtu21*(Gt112*gt11L*Gt222 + Gt112*Gt211*gt22L + Gt211*Gt222*gt22L + 
                2*Gt122*gt13L*Gt311 + Gt112*gt23L*Gt311 + Gt222*gt23L*Gt311 + 
                gt13L*Gt222*Gt312 + Gt213*gt22L*Gt312 + Gt212*gt23L*Gt312 + 
                gt23L*Gt312*Gt313 + Gt113*gt11L*Gt322 + Gt211*gt23L*Gt322 + 
                gt13L*Gt313*Gt322 + Gt111*
                 (2*gt11L*Gt122 + Gt112*gt12L + gt12L*Gt222 + gt13L*Gt322) + 
                gt12L*(2*Gt122*Gt211 + Gt112*Gt212 + Gt212*Gt222 + Gt113*Gt312 + 
                   Gt213*Gt322) + Gt311*Gt322*gt33L + gt22L*SQR(Gt212)) + 
             gtu22*(gt11L*Gt122*Gt222 + 2*Gt212*Gt222*gt22L + 2*Gt122*gt13L*Gt312 + 
                Gt223*gt22L*Gt312 + Gt222*gt23L*Gt312 + gt11L*Gt123*Gt322 + 
                gt13L*Gt222*Gt322 + 2*Gt212*gt23L*Gt322 + 
                Gt112*(2*gt11L*Gt122 + gt12L*Gt222 + Gt212*gt22L + gt23L*Gt312 + 
                   gt13L*Gt322) + gt23L*Gt312*Gt323 + gt13L*Gt322*Gt323 + 
                Gt312*Gt322*gt33L + gt12L*SQR(Gt112) + 
                gt12L*(3*Gt122*Gt212 + Gt123*Gt312 + Gt223*Gt322 + SQR(Gt222))) + 
             gtu33*(gt11L*Gt123*Gt223 + 2*Gt213*Gt223*gt22L + 2*Gt123*gt13L*Gt313 + 
                gt22L*Gt233*Gt313 + Gt223*gt23L*Gt313 + gt11L*Gt133*Gt323 + 
                gt13L*Gt223*Gt323 + 2*Gt213*gt23L*Gt323 + 
                Gt113*(2*gt11L*Gt123 + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + 
                   gt13L*Gt323) + gt23L*Gt313*Gt333 + gt13L*Gt323*Gt333 + 
                Gt313*Gt323*gt33L + gt12L*SQR(Gt113) + 
                gt12L*(3*Gt123*Gt213 + Gt133*Gt313 + Gt233*Gt323 + SQR(Gt223))) + 
             gtu21*(Gt122*gt12L*Gt211 + 3*Gt112*gt12L*Gt212 + gt12L*Gt212*Gt222 + 
                Gt211*Gt222*gt22L + Gt123*gt12L*Gt311 + Gt223*gt22L*Gt311 + 
                3*Gt112*gt13L*Gt312 + gt12L*Gt223*Gt312 + 2*Gt212*gt23L*Gt312 + 
                Gt111*(Gt112*gt12L + Gt212*gt22L + gt23L*Gt312) + 
                gt13L*Gt212*Gt322 + Gt211*gt23L*Gt322 + gt23L*Gt311*Gt323 + 
                gt13L*Gt312*Gt323 + gt11L*
                 (Gt122*Gt212 + Gt123*Gt312 + 2*SQR(Gt112)) + gt22L*SQR(Gt212) + 
                gt33L*SQR(Gt312)) + gtu31*
              (Gt112*gt11L*Gt223 + Gt113*Gt211*gt22L + Gt212*Gt213*gt22L + 
                Gt211*Gt223*gt22L + 2*Gt123*gt13L*Gt311 + Gt113*gt23L*Gt311 + 
                Gt223*gt23L*Gt311 + gt13L*Gt223*Gt312 + Gt213*gt23L*Gt312 + 
                Gt213*gt22L*Gt313 + Gt113*gt11L*Gt323 + Gt211*gt23L*Gt323 + 
                gt13L*Gt313*Gt323 + Gt111*
                 (2*gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + gt13L*Gt323) + 
                gt12L*(2*Gt123*Gt211 + Gt112*Gt213 + Gt212*Gt223 + Gt113*Gt313 + 
                   Gt213*Gt323) + Gt311*Gt323*gt33L + gt23L*SQR(Gt313)) + 
             gtu32*(gt11L*Gt122*Gt223 + Gt113*Gt212*gt22L + Gt213*Gt222*gt22L + 
                Gt212*Gt223*gt22L + 2*Gt123*gt13L*Gt312 + Gt113*gt23L*Gt312 + 
                Gt223*gt23L*Gt312 + Gt223*gt22L*Gt313 + gt13L*Gt223*Gt322 + 
                Gt213*gt23L*Gt322 + gt11L*Gt123*Gt323 + Gt212*gt23L*Gt323 + 
                gt23L*Gt313*Gt323 + Gt112*
                 (2*gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + gt13L*Gt323) + 
                gt12L*(Gt122*Gt213 + Gt123*(2*Gt212 + Gt313) + 
                   Gt223*(Gt222 + Gt323)) + Gt312*Gt323*gt33L + gt13L*SQR(Gt323)))\
           - gtu11*(2*J11L*J21L*PDstandardNth12gt12 + 
             2*J11L*J31L*PDstandardNth13gt12 + dJ111L*PDstandardNth1gt12 + 
             2*J21L*J31L*PDstandardNth23gt12 + dJ211L*PDstandardNth2gt12 + 
             dJ311L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J11L) + 
             PDstandardNth22gt12*SQR(J21L) + PDstandardNth33gt12*SQR(J31L)) - 
          gtu22*(2*J12L*J22L*PDstandardNth12gt12 + 
             2*J12L*J32L*PDstandardNth13gt12 + dJ122L*PDstandardNth1gt12 + 
             2*J22L*J32L*PDstandardNth23gt12 + dJ222L*PDstandardNth2gt12 + 
             dJ322L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J12L) + 
             PDstandardNth22gt12*SQR(J22L) + PDstandardNth33gt12*SQR(J32L)) - 
          gtu33*(2*J13L*J23L*PDstandardNth12gt12 + 
             2*J13L*J33L*PDstandardNth13gt12 + dJ133L*PDstandardNth1gt12 + 
             2*J23L*J33L*PDstandardNth23gt12 + dJ233L*PDstandardNth2gt12 + 
             dJ333L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J13L) + 
             PDstandardNth22gt12*SQR(J23L) + PDstandardNth33gt12*SQR(J33L)));
    
    Rt13  =  khalf*((gt13L*J11L + gt11L*J13L)*PDstandardNth1Xt1 + 
          (gt23L*J11L + gt12L*J13L)*PDstandardNth1Xt2 + 
          (gt33L*J11L + gt13L*J13L)*PDstandardNth1Xt3 + 
          (gt13L*J21L + gt11L*J23L)*PDstandardNth2Xt1 + 
          (gt23L*J21L + gt12L*J23L)*PDstandardNth2Xt2 + 
          (gt33L*J21L + gt13L*J23L)*PDstandardNth2Xt3 - 
          2*(gtu21*(J11L*J12L*PDstandardNth11gt13 + J12L*J21L*PDstandardNth12gt13 + 
                J11L*J22L*PDstandardNth12gt13 + J12L*J31L*PDstandardNth13gt13 + 
                J11L*J32L*PDstandardNth13gt13 + dJ112L*PDstandardNth1gt13 + 
                J21L*J22L*PDstandardNth22gt13 + J22L*J31L*PDstandardNth23gt13 + 
                J21L*J32L*PDstandardNth23gt13 + dJ212L*PDstandardNth2gt13 + 
                J31L*J32L*PDstandardNth33gt13 + dJ312L*PDstandardNth3gt13) + 
             gtu31*(J11L*J13L*PDstandardNth11gt13 + J13L*J21L*PDstandardNth12gt13 + 
                J11L*J23L*PDstandardNth12gt13 + J13L*J31L*PDstandardNth13gt13 + 
                J11L*J33L*PDstandardNth13gt13 + dJ113L*PDstandardNth1gt13 + 
                J21L*J23L*PDstandardNth22gt13 + J23L*J31L*PDstandardNth23gt13 + 
                J21L*J33L*PDstandardNth23gt13 + dJ213L*PDstandardNth2gt13 + 
                J31L*J33L*PDstandardNth33gt13 + dJ313L*PDstandardNth3gt13) + 
             gtu32*(J12L*J13L*PDstandardNth11gt13 + J13L*J22L*PDstandardNth12gt13 + 
                J12L*J23L*PDstandardNth12gt13 + J13L*J32L*PDstandardNth13gt13 + 
                J12L*J33L*PDstandardNth13gt13 + dJ123L*PDstandardNth1gt13 + 
                J22L*J23L*PDstandardNth22gt13 + J23L*J32L*PDstandardNth23gt13 + 
                J22L*J33L*PDstandardNth23gt13 + dJ223L*PDstandardNth2gt13 + 
                J32L*J33L*PDstandardNth33gt13 + dJ323L*PDstandardNth3gt13)) + 
          (gt13L*J31L + gt11L*J33L)*PDstandardNth3Xt1 + 
          (gt23L*J31L + gt12L*J33L)*PDstandardNth3Xt2 + 
          (gt33L*J31L + gt13L*J33L)*PDstandardNth3Xt3 + 
          (Gt113*gt11L + Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L + gt13L*Gt313 + 
             Gt311*gt33L)*Xtn1 + (gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + 
             Gt212*gt23L + gt13L*Gt323 + Gt312*gt33L)*Xtn2 + 
          (gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + Gt213*gt23L + gt13L*Gt333 + 
             Gt313*gt33L)*Xtn3 + 2*((Gt122*gt13L*Gt211 + 2*Gt113*gt12L*Gt212 + 
                Gt112*gt12L*Gt213 + gt12L*Gt213*Gt222 + Gt212*Gt213*gt22L + 
                Gt211*Gt222*gt23L + Gt123*gt13L*Gt311 + Gt223*gt23L*Gt311 + 
                2*Gt113*gt13L*Gt312 + Gt213*gt23L*Gt312 + Gt112*gt13L*Gt313 + 
                gt12L*Gt223*Gt313 + Gt212*gt23L*Gt313 + 
                gt11L*(2*Gt112*Gt113 + Gt122*Gt213 + Gt123*Gt313) + 
                gt13L*Gt213*Gt322 + gt13L*Gt313*Gt323 + Gt312*Gt313*gt33L + 
                Gt211*Gt322*gt33L + Gt311*Gt323*gt33L + 
                Gt111*(Gt112*gt13L + Gt212*gt23L + Gt312*gt33L))*gtu21 + 
             (Gt122*gt13L*Gt213 + gt11L*Gt122*Gt233 + Gt212*gt22L*Gt233 + 
                Gt113*Gt212*gt23L + Gt213*Gt222*gt23L + 2*Gt133*gt13L*Gt312 + 
                Gt233*gt23L*Gt312 + Gt123*gt13L*Gt313 + Gt223*gt23L*Gt313 + 
                gt13L*Gt233*Gt322 + gt11L*Gt123*Gt333 + Gt212*gt23L*Gt333 + 
                gt13L*Gt323*Gt333 + Gt112*
                 (2*gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + gt13L*Gt333) + 
                gt12L*(2*Gt133*Gt212 + Gt222*Gt233 + Gt223*Gt333) + 
                Gt113*Gt312*gt33L + Gt213*Gt322*gt33L + Gt313*Gt323*gt33L + 
                Gt312*Gt333*gt33L)*gtu32 + 
             gtu21*(2*Gt123*gt12L*Gt211 + Gt112*gt13L*Gt212 + gt12L*Gt212*Gt223 + 
                Gt211*Gt223*gt22L + Gt112*Gt211*gt23L + 2*Gt123*gt13L*Gt311 + 
                Gt223*gt23L*Gt311 + Gt113*gt13L*Gt312 + gt13L*Gt223*Gt312 + 
                Gt213*gt23L*Gt312 + gt12L*Gt213*Gt323 + Gt211*gt23L*Gt323 + 
                gt13L*Gt313*Gt323 + gt11L*
                 (2*Gt111*Gt123 + Gt112*Gt223 + Gt113*Gt323) + 
                Gt111*(Gt112*gt13L + gt12L*Gt223 + gt13L*Gt323) + 
                Gt112*Gt311*gt33L + Gt212*Gt312*gt33L + Gt312*Gt313*gt33L + 
                Gt311*Gt323*gt33L + gt23L*SQR(Gt212)) + 
             gtu32*(Gt123*gt13L*Gt212 + 2*Gt123*gt12L*Gt213 + Gt113*gt12L*Gt223 + 
                Gt213*Gt223*gt22L + Gt212*Gt223*gt23L + Gt133*gt13L*Gt312 + 
                Gt233*gt23L*Gt312 + 2*Gt123*gt13L*Gt313 + Gt223*gt23L*Gt313 + 
                Gt113*gt13L*Gt323 + gt13L*Gt223*Gt323 + gt12L*Gt233*Gt323 + 
                Gt213*gt23L*Gt323 + gt11L*
                 (2*Gt113*Gt123 + Gt123*Gt223 + Gt133*Gt323) + gt13L*Gt323*Gt333 + 
                Gt212*Gt323*gt33L + Gt313*Gt323*gt33L + Gt312*Gt333*gt33L + 
                Gt112*(Gt113*gt13L + Gt213*gt23L + Gt313*gt33L) + gt12L*SQR(Gt223))\
              + gtu11*(2*Gt113*gt12L*Gt211 + Gt112*gt13L*Gt211 + 
                gt12L*Gt212*Gt213 + Gt211*Gt213*gt22L + Gt211*Gt212*gt23L + 
                3*Gt113*gt13L*Gt311 + 2*Gt213*gt23L*Gt311 + gt13L*Gt213*Gt312 + 
                gt12L*Gt213*Gt313 + Gt211*gt23L*Gt313 + 
                gt11L*(2*Gt111*Gt113 + Gt112*Gt213 + Gt113*Gt313) + 
                Gt211*Gt312*gt33L + 2*Gt311*Gt313*gt33L + 
                Gt111*(gt12L*Gt213 + Gt211*gt23L + gt13L*Gt313 + Gt311*gt33L) + 
                gt13L*SQR(Gt111) + gt13L*SQR(Gt313)) + 
             gtu31*(Gt112*gt13L*Gt213 + Gt112*gt11L*Gt233 + Gt211*gt22L*Gt233 + 
                Gt113*Gt211*gt23L + Gt212*Gt213*gt23L + 2*Gt133*gt13L*Gt311 + 
                Gt233*gt23L*Gt311 + gt13L*Gt233*Gt312 + Gt113*gt13L*Gt313 + 
                Gt213*gt23L*Gt313 + Gt113*gt11L*Gt333 + Gt211*gt23L*Gt333 + 
                gt13L*Gt313*Gt333 + Gt111*
                 (2*gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + gt13L*Gt333) + 
                gt12L*(2*Gt133*Gt211 + Gt212*Gt233 + Gt213*Gt333) + 
                Gt113*Gt311*gt33L + Gt213*Gt312*gt33L + Gt311*Gt333*gt33L + 
                gt33L*SQR(Gt313)) + gtu31*
              (Gt123*gt13L*Gt211 + 3*Gt113*gt12L*Gt213 + gt12L*Gt213*Gt223 + 
                Gt211*Gt223*gt23L + Gt133*gt13L*Gt311 + Gt233*gt23L*Gt311 + 
                3*Gt113*gt13L*Gt313 + gt12L*Gt233*Gt313 + 2*Gt213*gt23L*Gt313 + 
                gt13L*Gt213*Gt323 + gt13L*Gt313*Gt333 + Gt211*Gt323*gt33L + 
                Gt311*Gt333*gt33L + Gt111*
                 (Gt113*gt13L + Gt213*gt23L + Gt313*gt33L) + 
                gt11L*(Gt123*Gt213 + Gt133*Gt313 + 2*SQR(Gt113)) + 
                gt22L*SQR(Gt213) + gt33L*SQR(Gt313)) + 
             gtu22*(2*Gt123*gt12L*Gt212 + Gt122*gt13L*Gt212 + gt12L*Gt222*Gt223 + 
                Gt212*Gt223*gt22L + Gt212*Gt222*gt23L + 3*Gt123*gt13L*Gt312 + 
                2*Gt223*gt23L*Gt312 + gt13L*Gt223*Gt322 + gt12L*Gt223*Gt323 + 
                Gt212*gt23L*Gt323 + gt11L*
                 (2*Gt112*Gt123 + Gt122*Gt223 + Gt123*Gt323) + Gt212*Gt322*gt33L + 
                2*Gt312*Gt323*gt33L + 
                Gt112*(gt12L*Gt223 + Gt212*gt23L + gt13L*Gt323 + Gt312*gt33L) + 
                gt13L*SQR(Gt112) + gt13L*SQR(Gt323)) + 
             gtu33*(2*gt12L*Gt133*Gt213 + Gt123*gt13L*Gt213 + gt11L*Gt123*Gt233 + 
                gt12L*Gt223*Gt233 + Gt213*gt22L*Gt233 + Gt213*Gt223*gt23L + 
                3*Gt133*gt13L*Gt313 + 2*Gt233*gt23L*Gt313 + gt13L*Gt233*Gt323 + 
                gt11L*Gt133*Gt333 + gt12L*Gt233*Gt333 + Gt213*gt23L*Gt333 + 
                Gt213*Gt323*gt33L + 2*Gt313*Gt333*gt33L + 
                Gt113*(2*gt11L*Gt133 + gt12L*Gt233 + Gt213*gt23L + gt13L*Gt333 + 
                   Gt313*gt33L) + gt13L*SQR(Gt113) + gt13L*SQR(Gt333))) - 
          gtu11*(2*J11L*J21L*PDstandardNth12gt13 + 
             2*J11L*J31L*PDstandardNth13gt13 + dJ111L*PDstandardNth1gt13 + 
             2*J21L*J31L*PDstandardNth23gt13 + dJ211L*PDstandardNth2gt13 + 
             dJ311L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J11L) + 
             PDstandardNth22gt13*SQR(J21L) + PDstandardNth33gt13*SQR(J31L)) - 
          gtu22*(2*J12L*J22L*PDstandardNth12gt13 + 
             2*J12L*J32L*PDstandardNth13gt13 + dJ122L*PDstandardNth1gt13 + 
             2*J22L*J32L*PDstandardNth23gt13 + dJ222L*PDstandardNth2gt13 + 
             dJ322L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J12L) + 
             PDstandardNth22gt13*SQR(J22L) + PDstandardNth33gt13*SQR(J32L)) - 
          gtu33*(2*J13L*J23L*PDstandardNth12gt13 + 
             2*J13L*J33L*PDstandardNth13gt13 + dJ133L*PDstandardNth1gt13 + 
             2*J23L*J33L*PDstandardNth23gt13 + dJ233L*PDstandardNth2gt13 + 
             dJ333L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J13L) + 
             PDstandardNth22gt13*SQR(J23L) + PDstandardNth33gt13*SQR(J33L)));
    
    Rt22  =  (Gt223*(3*Gt112*gt12L + 6*Gt212*gt22L + 4*gt23L*Gt312) + 
           Gt123*(Gt112*gt11L + gt12L*(2*Gt111 + 4*Gt212) + 
              2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + 
           Gt112*(gt11L*Gt123 + gt12L*(2*Gt113 + Gt223) + 
              2*(Gt213*gt22L + gt23L*Gt313) + gt13L*Gt323) + 
           2*(Gt113*gt12L*Gt323 + Gt312*
               (gt12L*Gt133 + gt22L*Gt233 + gt23L*Gt333)) + 
           Gt323*(Gt112*gt13L + 4*Gt212*gt23L + 
              2*(Gt213*gt22L + gt23L*Gt313 + Gt312*gt33L)))*gtu31 + 
        J12L*(gt12L*PDstandardNth1Xt1 + gt22L*PDstandardNth1Xt2 + 
           gt23L*PDstandardNth1Xt3) + 
        J22L*(gt12L*PDstandardNth2Xt1 + gt22L*PDstandardNth2Xt2 + 
           gt23L*PDstandardNth2Xt3) + 
        J32L*(gt12L*PDstandardNth3Xt1 + gt22L*PDstandardNth3Xt2 + 
           gt23L*PDstandardNth3Xt3) + 
        (Gt112*gt12L + Gt212*gt22L + gt23L*Gt312)*Xtn1 + 
        (Gt122*gt12L + Gt222*gt22L + gt23L*Gt322)*Xtn2 + 
        (Gt123*gt12L + Gt223*gt22L + gt23L*Gt323)*Xtn3 + 
        gtu21*(Gt222*(3*Gt112*gt12L + 6*Gt212*gt22L + 4*gt23L*Gt312) + 
           Gt122*(Gt112*gt11L + gt12L*(2*Gt111 + 4*Gt212) + 
              2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + 
           Gt112*(gt11L*Gt122 + gt12L*Gt222 + 2*(Gt212*gt22L + gt23L*Gt312) + 
              gt13L*Gt322) + Gt322*(Gt112*gt13L + 4*Gt212*gt23L + 
              2*(Gt213*gt22L + gt23L*Gt313 + Gt312*gt33L)) + 
           2*(Gt312*(Gt123*gt12L + Gt223*gt22L + gt23L*Gt323) + 
              gt12L*(Gt113*Gt322 + SQR(Gt112)))) + 
        gtu11*(Gt112*(gt12L*(2*Gt111 + 4*Gt212) + 
              2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + 
           Gt312*(2*(Gt113*gt12L + Gt213*gt22L) + gt23L*(4*Gt212 + 2*Gt313)) + 
           gt11L*SQR(Gt112) + 3*gt22L*SQR(Gt212) + gt33L*SQR(Gt312)) + 
        gtu22*(Gt122*(gt12L*(2*Gt112 + 4*Gt222) + 
              2*(Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322)) + 
           Gt322*(2*(Gt123*gt12L + Gt223*gt22L) + gt23L*(4*Gt222 + 2*Gt323)) + 
           gt11L*SQR(Gt122) + 3*gt22L*SQR(Gt222) + gt33L*SQR(Gt322)) + 
        gtu33*(Gt123*(gt12L*(2*Gt113 + 4*Gt223) + 
              2*(Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323)) + 
           Gt323*(2*(gt12L*Gt133 + gt22L*Gt233) + gt23L*(4*Gt223 + 2*Gt333)) + 
           gt11L*SQR(Gt123) + 3*gt22L*SQR(Gt223) + gt33L*SQR(Gt323)) + 
        gtu32*(gt22L*(2*(Gt122*Gt213 + Gt233*Gt322) + Gt223*(6*Gt222 + 2*Gt323)) + 
           4*(gt12L*(Gt123*Gt222 + Gt122*Gt223) + 
              gt23L*(Gt223*Gt322 + Gt222*Gt323)) + 
           2*(Gt123*(Gt112*gt12L + Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322) + 
              gt12L*(Gt133*Gt322 + Gt123*Gt323) + 
              Gt122*(gt11L*Gt123 + Gt113*gt12L + gt23L*Gt313 + gt13L*Gt323) + 
              Gt322*(gt23L*Gt333 + Gt323*gt33L) + gt23L*SQR(Gt323))) + 
        khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt22 + 
                 J12L*J21L*PDstandardNth12gt22 + J11L*J22L*PDstandardNth12gt22 + 
                 J12L*J31L*PDstandardNth13gt22 + J11L*J32L*PDstandardNth13gt22 + 
                 dJ112L*PDstandardNth1gt22 + J21L*J22L*PDstandardNth22gt22 + 
                 J22L*J31L*PDstandardNth23gt22 + J21L*J32L*PDstandardNth23gt22 + 
                 dJ212L*PDstandardNth2gt22 + J31L*J32L*PDstandardNth33gt22 + 
                 dJ312L*PDstandardNth3gt22) + 
              gtu31*(J11L*J13L*PDstandardNth11gt22 + 
                 J13L*J21L*PDstandardNth12gt22 + J11L*J23L*PDstandardNth12gt22 + 
                 J13L*J31L*PDstandardNth13gt22 + J11L*J33L*PDstandardNth13gt22 + 
                 dJ113L*PDstandardNth1gt22 + J21L*J23L*PDstandardNth22gt22 + 
                 J23L*J31L*PDstandardNth23gt22 + J21L*J33L*PDstandardNth23gt22 + 
                 dJ213L*PDstandardNth2gt22 + J31L*J33L*PDstandardNth33gt22 + 
                 dJ313L*PDstandardNth3gt22) + 
              gtu32*(J12L*J13L*PDstandardNth11gt22 + 
                 J13L*J22L*PDstandardNth12gt22 + J12L*J23L*PDstandardNth12gt22 + 
                 J13L*J32L*PDstandardNth13gt22 + J12L*J33L*PDstandardNth13gt22 + 
                 dJ123L*PDstandardNth1gt22 + J22L*J23L*PDstandardNth22gt22 + 
                 J23L*J32L*PDstandardNth23gt22 + J22L*J33L*PDstandardNth23gt22 + 
                 dJ223L*PDstandardNth2gt22 + J32L*J33L*PDstandardNth33gt22 + 
                 dJ323L*PDstandardNth3gt22)) - 
           gtu11*(2*J11L*J21L*PDstandardNth12gt22 + 
              2*J11L*J31L*PDstandardNth13gt22 + dJ111L*PDstandardNth1gt22 + 
              2*J21L*J31L*PDstandardNth23gt22 + dJ211L*PDstandardNth2gt22 + 
              dJ311L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J11L) + 
              PDstandardNth22gt22*SQR(J21L) + PDstandardNth33gt22*SQR(J31L)) - 
           gtu22*(2*J12L*J22L*PDstandardNth12gt22 + 
              2*J12L*J32L*PDstandardNth13gt22 + dJ122L*PDstandardNth1gt22 + 
              2*J22L*J32L*PDstandardNth23gt22 + dJ222L*PDstandardNth2gt22 + 
              dJ322L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J12L) + 
              PDstandardNth22gt22*SQR(J22L) + PDstandardNth33gt22*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt22 + 
              2*J13L*J33L*PDstandardNth13gt22 + dJ133L*PDstandardNth1gt22 + 
              2*J23L*J33L*PDstandardNth23gt22 + dJ233L*PDstandardNth2gt22 + 
              dJ333L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J13L) + 
              PDstandardNth22gt22*SQR(J23L) + PDstandardNth33gt22*SQR(J33L)));
    
    Rt23  =  khalf*((gt13L*J12L + gt12L*J13L)*PDstandardNth1Xt1 + 
          (gt23L*J12L + gt22L*J13L)*PDstandardNth1Xt2 + 
          (gt33L*J12L + gt23L*J13L)*PDstandardNth1Xt3 + 
          (gt13L*J22L + gt12L*J23L)*PDstandardNth2Xt1 + 
          (gt23L*J22L + gt22L*J23L)*PDstandardNth2Xt2 + 
          (gt33L*J22L + gt23L*J23L)*PDstandardNth2Xt3 - 
          2*(gtu21*(J11L*J12L*PDstandardNth11gt23 + J12L*J21L*PDstandardNth12gt23 + 
                J11L*J22L*PDstandardNth12gt23 + J12L*J31L*PDstandardNth13gt23 + 
                J11L*J32L*PDstandardNth13gt23 + dJ112L*PDstandardNth1gt23 + 
                J21L*J22L*PDstandardNth22gt23 + J22L*J31L*PDstandardNth23gt23 + 
                J21L*J32L*PDstandardNth23gt23 + dJ212L*PDstandardNth2gt23 + 
                J31L*J32L*PDstandardNth33gt23 + dJ312L*PDstandardNth3gt23) + 
             gtu31*(J11L*J13L*PDstandardNth11gt23 + J13L*J21L*PDstandardNth12gt23 + 
                J11L*J23L*PDstandardNth12gt23 + J13L*J31L*PDstandardNth13gt23 + 
                J11L*J33L*PDstandardNth13gt23 + dJ113L*PDstandardNth1gt23 + 
                J21L*J23L*PDstandardNth22gt23 + J23L*J31L*PDstandardNth23gt23 + 
                J21L*J33L*PDstandardNth23gt23 + dJ213L*PDstandardNth2gt23 + 
                J31L*J33L*PDstandardNth33gt23 + dJ313L*PDstandardNth3gt23) + 
             gtu32*(J12L*J13L*PDstandardNth11gt23 + J13L*J22L*PDstandardNth12gt23 + 
                J12L*J23L*PDstandardNth12gt23 + J13L*J32L*PDstandardNth13gt23 + 
                J12L*J33L*PDstandardNth13gt23 + dJ123L*PDstandardNth1gt23 + 
                J22L*J23L*PDstandardNth22gt23 + J23L*J32L*PDstandardNth23gt23 + 
                J22L*J33L*PDstandardNth23gt23 + dJ223L*PDstandardNth2gt23 + 
                J32L*J33L*PDstandardNth33gt23 + dJ323L*PDstandardNth3gt23)) + 
          (gt13L*J32L + gt12L*J33L)*PDstandardNth3Xt1 + 
          (gt23L*J32L + gt22L*J33L)*PDstandardNth3Xt2 + 
          (gt33L*J32L + gt23L*J33L)*PDstandardNth3Xt3 + 
          (Gt113*gt12L + Gt112*gt13L + Gt213*gt22L + Gt212*gt23L + gt23L*Gt313 + 
             Gt312*gt33L)*Xtn1 + (Gt123*gt12L + Gt122*gt13L + Gt223*gt22L + 
             Gt222*gt23L + gt23L*Gt323 + Gt322*gt33L)*Xtn2 + 
          (gt12L*Gt133 + Gt123*gt13L + gt22L*Gt233 + Gt223*gt23L + gt23L*Gt333 + 
             Gt323*gt33L)*Xtn3 + 2*((Gt112*gt11L*Gt123 + Gt111*Gt123*gt12L + 
                Gt111*Gt122*gt13L + Gt123*gt12L*Gt212 + Gt112*gt13L*Gt222 + 
                2*Gt112*gt12L*Gt223 + Gt123*Gt211*gt22L + 2*Gt212*Gt223*gt22L + 
                Gt122*Gt211*gt23L + Gt212*Gt222*gt23L + Gt123*gt23L*Gt311 + 
                Gt123*gt13L*Gt312 + 2*Gt223*gt23L*Gt312 + Gt113*gt13L*Gt322 + 
                Gt213*gt23L*Gt322 + Gt113*gt12L*Gt323 + Gt112*gt13L*Gt323 + 
                Gt213*gt22L*Gt323 + Gt212*gt23L*Gt323 + gt23L*Gt313*Gt323 + 
                Gt122*Gt311*gt33L + Gt222*Gt312*gt33L + Gt313*Gt322*gt33L + 
                Gt312*Gt323*gt33L)*gtu21 + 
             (Gt112*gt11L*Gt133 + Gt111*gt12L*Gt133 + Gt111*Gt123*gt13L + 
                gt12L*Gt133*Gt212 + Gt112*gt13L*Gt223 + Gt133*Gt211*gt22L + 
                2*Gt112*gt12L*Gt233 + 2*Gt212*gt22L*Gt233 + Gt123*Gt211*gt23L + 
                Gt212*Gt223*gt23L + Gt133*gt23L*Gt311 + Gt133*gt13L*Gt312 + 
                2*Gt233*gt23L*Gt312 + Gt113*gt13L*Gt323 + Gt213*gt23L*Gt323 + 
                Gt113*gt12L*Gt333 + Gt112*gt13L*Gt333 + Gt213*gt22L*Gt333 + 
                Gt212*gt23L*Gt333 + gt23L*Gt313*Gt333 + Gt123*Gt311*gt33L + 
                Gt223*Gt312*gt33L + Gt313*Gt323*gt33L + Gt312*Gt333*gt33L)*gtu31 + 
             gtu21*(Gt113*gt11L*Gt122 + Gt122*gt13L*Gt212 + 2*Gt122*gt12L*Gt213 + 
                Gt113*gt12L*Gt222 + Gt113*Gt212*gt22L + 2*Gt213*Gt222*gt22L + 
                Gt212*Gt222*gt23L + Gt123*gt13L*Gt312 + Gt113*gt23L*Gt312 + 
                Gt223*gt23L*Gt312 + Gt123*gt12L*Gt313 + Gt122*gt13L*Gt313 + 
                Gt223*gt22L*Gt313 + Gt222*gt23L*Gt313 + Gt113*gt13L*Gt322 + 
                2*Gt213*gt23L*Gt322 + gt23L*Gt313*Gt323 + Gt212*Gt322*gt33L + 
                Gt313*Gt322*gt33L + Gt312*Gt323*gt33L + 
                Gt112*(Gt113*gt12L + Gt212*gt23L + Gt312*gt33L) + gt13L*SQR(Gt112))\
              + gtu31*(2*Gt213*Gt223*gt22L + Gt112*Gt213*gt23L + 
                Gt212*Gt223*gt23L + Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + 
                gt12L*Gt133*Gt313 + gt22L*Gt233*Gt313 + Gt223*gt23L*Gt313 + 
                Gt123*(2*gt12L*Gt213 + gt13L*(Gt212 + Gt313)) + 
                2*Gt213*gt23L*Gt323 + 
                Gt113*(gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt213*gt22L + 
                   gt23L*Gt313 + gt13L*Gt323) + gt23L*Gt313*Gt333 + 
                Gt112*Gt313*gt33L + Gt212*Gt323*gt33L + Gt313*Gt323*gt33L + 
                Gt312*Gt333*gt33L + gt12L*SQR(Gt113)) + 
             gtu11*(Gt112*Gt113*gt11L + Gt111*Gt113*gt12L + Gt111*Gt112*gt13L + 
                Gt113*gt12L*Gt212 + Gt112*gt13L*Gt212 + 2*Gt112*gt12L*Gt213 + 
                Gt113*Gt211*gt22L + 2*Gt212*Gt213*gt22L + Gt112*Gt211*gt23L + 
                Gt113*gt23L*Gt311 + 2*Gt113*gt13L*Gt312 + 3*Gt213*gt23L*Gt312 + 
                Gt113*gt12L*Gt313 + Gt112*gt13L*Gt313 + Gt213*gt22L*Gt313 + 
                Gt212*gt23L*Gt313 + Gt112*Gt311*gt33L + Gt212*Gt312*gt33L + 
                2*Gt312*Gt313*gt33L + gt23L*SQR(Gt212) + gt23L*SQR(Gt313)) + 
             gtu22*(gt11L*Gt122*Gt123 + Gt112*Gt123*gt12L + Gt112*Gt122*gt13L + 
                Gt123*gt12L*Gt222 + Gt122*gt13L*Gt222 + 2*Gt122*gt12L*Gt223 + 
                Gt123*Gt212*gt22L + 2*Gt222*Gt223*gt22L + Gt122*Gt212*gt23L + 
                Gt123*gt23L*Gt312 + 2*Gt123*gt13L*Gt322 + 3*Gt223*gt23L*Gt322 + 
                Gt123*gt12L*Gt323 + Gt122*gt13L*Gt323 + Gt223*gt22L*Gt323 + 
                Gt222*gt23L*Gt323 + Gt122*Gt312*gt33L + Gt222*Gt322*gt33L + 
                2*Gt322*Gt323*gt33L + gt23L*SQR(Gt222) + gt23L*SQR(Gt323)) + 
             gtu32*(gt11L*Gt122*Gt133 + Gt112*gt12L*Gt133 + Gt112*Gt123*gt13L + 
                gt12L*Gt133*Gt222 + Gt122*gt13L*Gt223 + Gt133*Gt212*gt22L + 
                2*Gt122*gt12L*Gt233 + 2*Gt222*gt22L*Gt233 + Gt123*Gt212*gt23L + 
                Gt222*Gt223*gt23L + Gt133*gt23L*Gt312 + Gt133*gt13L*Gt322 + 
                2*Gt233*gt23L*Gt322 + Gt123*gt13L*Gt323 + Gt223*gt23L*Gt323 + 
                Gt123*gt12L*Gt333 + Gt122*gt13L*Gt333 + Gt223*gt22L*Gt333 + 
                Gt222*gt23L*Gt333 + gt23L*Gt323*Gt333 + Gt123*Gt312*gt33L + 
                Gt223*Gt322*gt33L + Gt322*Gt333*gt33L + gt33L*SQR(Gt323)) + 
             gtu32*(Gt113*Gt123*gt12L + Gt113*Gt122*gt13L + Gt123*gt13L*Gt222 + 
                3*Gt123*gt12L*Gt223 + Gt123*Gt213*gt22L + Gt122*Gt213*gt23L + 
                Gt222*Gt223*gt23L + Gt123*gt23L*Gt313 + Gt133*gt13L*Gt322 + 
                Gt233*gt23L*Gt322 + gt12L*Gt133*Gt323 + 2*Gt123*gt13L*Gt323 + 
                gt22L*Gt233*Gt323 + 3*Gt223*gt23L*Gt323 + gt23L*Gt323*Gt333 + 
                Gt122*Gt313*gt33L + Gt222*Gt323*gt33L + Gt322*Gt333*gt33L + 
                gt11L*SQR(Gt123) + 2*gt22L*SQR(Gt223) + gt33L*SQR(Gt323)) + 
             gtu33*(gt11L*Gt123*Gt133 + Gt113*gt12L*Gt133 + Gt113*Gt123*gt13L + 
                gt12L*Gt133*Gt223 + Gt123*gt13L*Gt223 + Gt133*Gt213*gt22L + 
                2*Gt123*gt12L*Gt233 + 2*Gt223*gt22L*Gt233 + Gt123*Gt213*gt23L + 
                Gt133*gt23L*Gt313 + 2*Gt133*gt13L*Gt323 + 3*Gt233*gt23L*Gt323 + 
                gt12L*Gt133*Gt333 + Gt123*gt13L*Gt333 + gt22L*Gt233*Gt333 + 
                Gt223*gt23L*Gt333 + Gt123*Gt313*gt33L + Gt223*Gt323*gt33L + 
                2*Gt323*Gt333*gt33L + gt23L*SQR(Gt223) + gt23L*SQR(Gt333))) - 
          gtu11*(2*J11L*J21L*PDstandardNth12gt23 + 
             2*J11L*J31L*PDstandardNth13gt23 + dJ111L*PDstandardNth1gt23 + 
             2*J21L*J31L*PDstandardNth23gt23 + dJ211L*PDstandardNth2gt23 + 
             dJ311L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J11L) + 
             PDstandardNth22gt23*SQR(J21L) + PDstandardNth33gt23*SQR(J31L)) - 
          gtu22*(2*J12L*J22L*PDstandardNth12gt23 + 
             2*J12L*J32L*PDstandardNth13gt23 + dJ122L*PDstandardNth1gt23 + 
             2*J22L*J32L*PDstandardNth23gt23 + dJ222L*PDstandardNth2gt23 + 
             dJ322L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J12L) + 
             PDstandardNth22gt23*SQR(J22L) + PDstandardNth33gt23*SQR(J32L)) - 
          gtu33*(2*J13L*J23L*PDstandardNth12gt23 + 
             2*J13L*J33L*PDstandardNth13gt23 + dJ133L*PDstandardNth1gt23 + 
             2*J23L*J33L*PDstandardNth23gt23 + dJ233L*PDstandardNth2gt23 + 
             dJ333L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J13L) + 
             PDstandardNth22gt23*SQR(J23L) + PDstandardNth33gt23*SQR(J33L)));
    
    Rt33  =  (4*((Gt123*gt13L + Gt223*gt23L)*Gt313 + 
              (Gt113*gt13L + Gt213*gt23L)*Gt323) + 
           (2*Gt213*Gt322 + 6*Gt313*Gt323)*gt33L + 
           2*(gt13L*(Gt122*Gt213 + Gt112*Gt223) + 
              Gt213*(Gt223*gt22L + Gt222*gt23L) + 
              Gt123*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L + Gt311*gt33L) + 
              Gt223*(Gt212*gt23L + Gt312*gt33L) + 
              Gt113*(gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + 
                 Gt312*gt33L)))*gtu21 + 
        J13L*(gt13L*PDstandardNth1Xt1 + gt23L*PDstandardNth1Xt2 + 
           gt33L*PDstandardNth1Xt3) + 
        J23L*(gt13L*PDstandardNth2Xt1 + gt23L*PDstandardNth2Xt2 + 
           gt33L*PDstandardNth2Xt3) + 
        J33L*(gt13L*PDstandardNth3Xt1 + gt23L*PDstandardNth3Xt2 + 
           gt33L*PDstandardNth3Xt3) + 
        (Gt113*gt13L + Gt213*gt23L + Gt313*gt33L)*Xtn1 + 
        (Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)*Xtn2 + 
        (Gt133*gt13L + Gt233*gt23L + Gt333*gt33L)*Xtn3 + 
        gtu31*(Gt133*(Gt113*gt11L + 2*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L) + 
              4*gt13L*Gt313) + Gt333*
            (3*Gt113*gt13L + 4*Gt213*gt23L + 6*Gt313*gt33L) + 
           Gt233*(Gt113*gt12L + 4*gt23L*Gt313 + 
              2*(Gt213*gt22L + Gt212*gt23L + Gt312*gt33L)) + 
           Gt113*(gt11L*Gt133 + gt12L*Gt233 + gt13L*Gt333 + 
              2*(Gt213*gt23L + Gt313*gt33L)) + 
           2*(Gt133*Gt311*gt33L + Gt213*(Gt223*gt23L + Gt323*gt33L) + 
              gt13L*(Gt123*Gt213 + Gt112*Gt233 + SQR(Gt113)))) + 
        gtu32*(4*((Gt133*gt13L + Gt233*gt23L)*Gt323 + 
              (Gt123*gt13L + Gt223*gt23L)*Gt333) + 
           Gt323*(2*Gt223 + 6*Gt333)*gt33L + 
           2*(Gt133*(Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L) + 
              Gt123*(gt11L*Gt133 + gt13L*(Gt113 + Gt223) + gt12L*Gt233 + 
                 Gt213*gt23L + Gt313*gt33L) + 
              Gt233*(Gt122*gt13L + Gt223*gt22L + Gt222*gt23L + Gt322*gt33L) + 
              gt23L*SQR(Gt223))) + gtu11*
         (4*(Gt113*gt13L + Gt213*gt23L)*Gt313 + 
           2*(Gt113*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L + Gt311*gt33L) + 
              Gt213*(Gt112*gt13L + Gt212*gt23L + Gt312*gt33L)) + gt11L*SQR(Gt113) + 
           gt22L*SQR(Gt213) + 3*gt33L*SQR(Gt313)) + 
        gtu22*(4*(Gt123*gt13L + Gt223*gt23L)*Gt323 + 
           2*(Gt123*(Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L) + 
              Gt223*(Gt122*gt13L + Gt222*gt23L + Gt322*gt33L)) + gt11L*SQR(Gt123) + 
           gt22L*SQR(Gt223) + 3*gt33L*SQR(Gt323)) + 
        gtu33*(4*(Gt133*gt13L + Gt233*gt23L)*Gt333 + 
           2*(Gt133*(Gt113*gt13L + gt12L*Gt233 + Gt213*gt23L + Gt313*gt33L) + 
              Gt233*(Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)) + gt11L*SQR(Gt133) + 
           gt22L*SQR(Gt233) + 3*gt33L*SQR(Gt333)) + 
        khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt33 + 
                 J12L*J21L*PDstandardNth12gt33 + J11L*J22L*PDstandardNth12gt33 + 
                 J12L*J31L*PDstandardNth13gt33 + J11L*J32L*PDstandardNth13gt33 + 
                 dJ112L*PDstandardNth1gt33 + J21L*J22L*PDstandardNth22gt33 + 
                 J22L*J31L*PDstandardNth23gt33 + J21L*J32L*PDstandardNth23gt33 + 
                 dJ212L*PDstandardNth2gt33 + J31L*J32L*PDstandardNth33gt33 + 
                 dJ312L*PDstandardNth3gt33) + 
              gtu31*(J11L*J13L*PDstandardNth11gt33 + 
                 J13L*J21L*PDstandardNth12gt33 + J11L*J23L*PDstandardNth12gt33 + 
                 J13L*J31L*PDstandardNth13gt33 + J11L*J33L*PDstandardNth13gt33 + 
                 dJ113L*PDstandardNth1gt33 + J21L*J23L*PDstandardNth22gt33 + 
                 J23L*J31L*PDstandardNth23gt33 + J21L*J33L*PDstandardNth23gt33 + 
                 dJ213L*PDstandardNth2gt33 + J31L*J33L*PDstandardNth33gt33 + 
                 dJ313L*PDstandardNth3gt33) + 
              gtu32*(J12L*J13L*PDstandardNth11gt33 + 
                 J13L*J22L*PDstandardNth12gt33 + J12L*J23L*PDstandardNth12gt33 + 
                 J13L*J32L*PDstandardNth13gt33 + J12L*J33L*PDstandardNth13gt33 + 
                 dJ123L*PDstandardNth1gt33 + J22L*J23L*PDstandardNth22gt33 + 
                 J23L*J32L*PDstandardNth23gt33 + J22L*J33L*PDstandardNth23gt33 + 
                 dJ223L*PDstandardNth2gt33 + J32L*J33L*PDstandardNth33gt33 + 
                 dJ323L*PDstandardNth3gt33)) - 
           gtu11*(2*J11L*J21L*PDstandardNth12gt33 + 
              2*J11L*J31L*PDstandardNth13gt33 + dJ111L*PDstandardNth1gt33 + 
              2*J21L*J31L*PDstandardNth23gt33 + dJ211L*PDstandardNth2gt33 + 
              dJ311L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J11L) + 
              PDstandardNth22gt33*SQR(J21L) + PDstandardNth33gt33*SQR(J31L)) - 
           gtu22*(2*J12L*J22L*PDstandardNth12gt33 + 
              2*J12L*J32L*PDstandardNth13gt33 + dJ122L*PDstandardNth1gt33 + 
              2*J22L*J32L*PDstandardNth23gt33 + dJ222L*PDstandardNth2gt33 + 
              dJ322L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J12L) + 
              PDstandardNth22gt33*SQR(J22L) + PDstandardNth33gt33*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt33 + 
              2*J13L*J33L*PDstandardNth13gt33 + dJ133L*PDstandardNth1gt33 + 
              2*J23L*J33L*PDstandardNth23gt33 + dJ233L*PDstandardNth2gt33 + 
              dJ333L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J13L) + 
              PDstandardNth22gt33*SQR(J23L) + PDstandardNth33gt33*SQR(J33L)));
    
    Rphi11  =  -2*((dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L)*
            PDstandardNth1phi + 2*(J11L*
               (J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + 
              J21L*J31L*PDstandardNth23phi) + 
           (dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L)*PDstandardNth2phi + 
           (dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L)*PDstandardNth3phi + 
           PDstandardNth11phi*SQR(J11L) + PDstandardNth22phi*SQR(J21L) + 
           PDstandardNth33phi*SQR(J31L) + 
           gt11L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + 
                 gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*
               PDstandardNth1phi + (2*
                  (dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + 
                 gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*
               PDstandardNth2phi + (dJ322L*gtu22 + 
                 2*(dJ313L*gtu31 + dJ323L*gtu32) + 
                 (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + 
                 (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + 
                    (Gt213*gtu31 + Gt223*gtu32)*J32L + 
                    (Gt313*gtu31 + Gt323*gtu32)*J33L))*PDstandardNth3phi + 
              2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*
                  PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*
                  PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + 
                    gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + 
                    (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + 
                    (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + 
                    gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + 
                 gtu33*SQR(J13L)) + PDstandardNth22phi*
               (gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + 
                 gtu33*SQR(J33L)))) - 
        4*gt11L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + 
                 (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*
               PDstandardNth3phi + PDstandardNth1phi*
               (((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*PDstandardNth2phi\
                  + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*PDstandardNth3phi))
             + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + 
              gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L))*
            SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + 
              gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + 
           (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + 
              gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L))*
            SQR(PDstandardNth3phi)) + 
        4*SQR(J11L*PDstandardNth1phi + J21L*PDstandardNth2phi + 
           J31L*PDstandardNth3phi);
    
    Rphi12  =  4*(J11L*PDstandardNth1phi + J21L*PDstandardNth2phi + 
           J31L*PDstandardNth3phi)*(J12L*PDstandardNth1phi + 
           J22L*PDstandardNth2phi + J32L*PDstandardNth3phi) - 
        2*(J12L*(J11L*PDstandardNth11phi + J21L*PDstandardNth12phi + 
              J31L*PDstandardNth13phi) + 
           J11L*(J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
           (dJ112L - Gt112*J11L - Gt212*J12L - Gt312*J13L)*PDstandardNth1phi + 
           J22L*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi) + 
           (dJ212L - Gt112*J21L - Gt212*J22L - Gt312*J23L)*PDstandardNth2phi + 
           J32L*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi) + 
           (dJ312L - Gt112*J31L - Gt212*J32L - Gt312*J33L)*PDstandardNth3phi + 
           gt12L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + 
                 gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*
               PDstandardNth1phi + (2*
                  (dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + 
                 gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*
               PDstandardNth2phi + (dJ322L*gtu22 + 
                 2*(dJ313L*gtu31 + dJ323L*gtu32) + 
                 (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + 
                 (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + 
                    (Gt213*gtu31 + Gt223*gtu32)*J32L + 
                    (Gt313*gtu31 + Gt323*gtu32)*J33L))*PDstandardNth3phi + 
              2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*
                  PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*
                  PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + 
                    gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + 
                    (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + 
                    (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + 
                    gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + 
                 gtu33*SQR(J13L)) + PDstandardNth22phi*
               (gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + 
                 gtu33*SQR(J33L)))) - 
        4*gt12L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + 
                 (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*
               PDstandardNth3phi + PDstandardNth1phi*
               (((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*PDstandardNth2phi\
                  + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*PDstandardNth3phi))
             + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + 
              gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L))*
            SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + 
              gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + 
           (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + 
              gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L))*
            SQR(PDstandardNth3phi));
    
    Rphi13  =  4*(J11L*PDstandardNth1phi + J21L*PDstandardNth2phi + 
           J31L*PDstandardNth3phi)*(J13L*PDstandardNth1phi + 
           J23L*PDstandardNth2phi + J33L*PDstandardNth3phi) - 
        2*(J13L*(J11L*PDstandardNth11phi + J21L*PDstandardNth12phi + 
              J31L*PDstandardNth13phi) + 
           J11L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
           (dJ113L - Gt113*J11L - Gt213*J12L - Gt313*J13L)*PDstandardNth1phi + 
           J23L*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi) + 
           (dJ213L - Gt113*J21L - Gt213*J22L - Gt313*J23L)*PDstandardNth2phi + 
           J33L*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi) + 
           (dJ313L - Gt113*J31L - Gt213*J32L - Gt313*J33L)*PDstandardNth3phi + 
           gt13L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + 
                 gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*
               PDstandardNth1phi + (2*
                  (dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + 
                 gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*
               PDstandardNth2phi + (dJ322L*gtu22 + 
                 2*(dJ313L*gtu31 + dJ323L*gtu32) + 
                 (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + 
                 (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + 
                    (Gt213*gtu31 + Gt223*gtu32)*J32L + 
                    (Gt313*gtu31 + Gt323*gtu32)*J33L))*PDstandardNth3phi + 
              2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*
                  PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*
                  PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + 
                    gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + 
                    (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + 
                    (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + 
                    gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + 
                 gtu33*SQR(J13L)) + PDstandardNth22phi*
               (gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + 
                 gtu33*SQR(J33L)))) - 
        4*gt13L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + 
                 (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*
               PDstandardNth3phi + PDstandardNth1phi*
               (((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*PDstandardNth2phi\
                  + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*PDstandardNth3phi))
             + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + 
              gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L))*
            SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + 
              gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + 
           (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + 
              gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L))*
            SQR(PDstandardNth3phi));
    
    Rphi22  =  -2*((dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L)*
            PDstandardNth1phi + 2*(J12L*
               (J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
              J22L*J32L*PDstandardNth23phi) + 
           (dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L)*PDstandardNth2phi + 
           (dJ322L - Gt122*J31L - Gt222*J32L - Gt322*J33L)*PDstandardNth3phi + 
           PDstandardNth11phi*SQR(J12L) + PDstandardNth22phi*SQR(J22L) + 
           PDstandardNth33phi*SQR(J32L) + 
           gt22L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + 
                 gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*
               PDstandardNth1phi + (2*
                  (dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + 
                 gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*
               PDstandardNth2phi + (dJ322L*gtu22 + 
                 2*(dJ313L*gtu31 + dJ323L*gtu32) + 
                 (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + 
                 (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + 
                    (Gt213*gtu31 + Gt223*gtu32)*J32L + 
                    (Gt313*gtu31 + Gt323*gtu32)*J33L))*PDstandardNth3phi + 
              2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*
                  PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*
                  PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + 
                    gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + 
                    (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + 
                    (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + 
                    gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + 
                 gtu33*SQR(J13L)) + PDstandardNth22phi*
               (gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + 
                 gtu33*SQR(J33L)))) - 
        4*gt22L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + 
                 (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*
               PDstandardNth3phi + PDstandardNth1phi*
               (((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*PDstandardNth2phi\
                  + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*PDstandardNth3phi))
             + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + 
              gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L))*
            SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + 
              gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + 
           (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + 
              gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L))*
            SQR(PDstandardNth3phi)) + 
        4*SQR(J12L*PDstandardNth1phi + J22L*PDstandardNth2phi + 
           J32L*PDstandardNth3phi);
    
    Rphi23  =  4*(J12L*PDstandardNth1phi + J22L*PDstandardNth2phi + 
           J32L*PDstandardNth3phi)*(J13L*PDstandardNth1phi + 
           J23L*PDstandardNth2phi + J33L*PDstandardNth3phi) - 
        2*(J13L*(J12L*PDstandardNth11phi + J22L*PDstandardNth12phi + 
              J32L*PDstandardNth13phi) + 
           J12L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
           (dJ123L - Gt123*J11L - Gt223*J12L - Gt323*J13L)*PDstandardNth1phi + 
           J23L*(J22L*PDstandardNth22phi + J32L*PDstandardNth23phi) + 
           (dJ223L - Gt123*J21L - Gt223*J22L - Gt323*J23L)*PDstandardNth2phi + 
           J33L*(J22L*PDstandardNth23phi + J32L*PDstandardNth33phi) + 
           (dJ323L - Gt123*J31L - Gt223*J32L - Gt323*J33L)*PDstandardNth3phi + 
           gt23L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + 
                 gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*
               PDstandardNth1phi + (2*
                  (dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + 
                 gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*
               PDstandardNth2phi + (dJ322L*gtu22 + 
                 2*(dJ313L*gtu31 + dJ323L*gtu32) + 
                 (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + 
                 (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + 
                    (Gt213*gtu31 + Gt223*gtu32)*J32L + 
                    (Gt313*gtu31 + Gt323*gtu32)*J33L))*PDstandardNth3phi + 
              2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*
                  PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*
                  PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + 
                    gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + 
                    (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + 
                    (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + 
                    gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + 
                 gtu33*SQR(J13L)) + PDstandardNth22phi*
               (gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + 
                 gtu33*SQR(J33L)))) - 
        4*gt23L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + 
                 (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*
               PDstandardNth3phi + PDstandardNth1phi*
               (((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*PDstandardNth2phi\
                  + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*PDstandardNth3phi))
             + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + 
              gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L))*
            SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + 
              gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + 
           (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + 
              gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L))*
            SQR(PDstandardNth3phi));
    
    Rphi33  =  -2*((dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L)*
            PDstandardNth1phi + 2*(J13L*
               (J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
              J23L*J33L*PDstandardNth23phi) + 
           (dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L)*PDstandardNth2phi + 
           (dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L)*PDstandardNth3phi + 
           PDstandardNth11phi*SQR(J13L) + PDstandardNth22phi*SQR(J23L) + 
           PDstandardNth33phi*SQR(J33L) + 
           gt33L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + 
                 gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*
               PDstandardNth1phi + (2*
                  (dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + 
                 gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + 
                    (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*
               PDstandardNth2phi + (dJ322L*gtu22 + 
                 2*(dJ313L*gtu31 + dJ323L*gtu32) + 
                 (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + 
                 (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + 
                    (Gt213*gtu31 + Gt223*gtu32)*J32L + 
                    (Gt313*gtu31 + Gt323*gtu32)*J33L))*PDstandardNth3phi + 
              2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*
                  PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*
                  PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + 
                    gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + 
                    (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + 
                    (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + 
                    gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + 
                 gtu33*SQR(J13L)) + PDstandardNth22phi*
               (gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + 
                 gtu33*SQR(J33L)))) - 
        4*gt33L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + 
                 (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*
               PDstandardNth3phi + PDstandardNth1phi*
               (((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*PDstandardNth2phi\
                  + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + 
                    (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*PDstandardNth3phi))
             + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + 
              gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L))*
            SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + 
              gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + 
           (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + 
              gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L))*
            SQR(PDstandardNth3phi)) + 
        4*SQR(J13L*PDstandardNth1phi + J23L*PDstandardNth2phi + 
           J33L*PDstandardNth3phi);
    
    Atm11  =  At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    Atm21  =  At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    Atm31  =  At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    Atm12  =  At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    Atm22  =  At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    Atm32  =  At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    Atm13  =  At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    Atm23  =  At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    Atm33  =  At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    Atu11  =  Atm11*gtu11 + Atm12*gtu21 + Atm13*gtu31;
    
    Atu21  =  Atm11*gtu21 + Atm12*gtu22 + Atm13*gtu32;
    
    Atu31  =  Atm11*gtu31 + Atm12*gtu32 + Atm13*gtu33;
    
    Atu22  =  Atm21*gtu21 + Atm22*gtu22 + Atm23*gtu32;
    
    Atu32  =  Atm21*gtu31 + Atm22*gtu32 + Atm23*gtu33;
    
    Atu33  =  Atm31*gtu31 + Atm32*gtu32 + Atm33*gtu33;
    
    e4phi  =  exp(4*phiL);
    
    em4phi  =  INV(e4phi);
    
    g11  =  e4phi*gt11L;
    
    g12  =  e4phi*gt12L;
    
    g13  =  e4phi*gt13L;
    
    g22  =  e4phi*gt22L;
    
    g23  =  e4phi*gt23L;
    
    g33  =  e4phi*gt33L;
    
    gu11  =  em4phi*gtu11;
    
    gu21  =  em4phi*gtu21;
    
    gu31  =  em4phi*gtu31;
    
    gu22  =  em4phi*gtu22;
    
    gu32  =  em4phi*gtu32;
    
    gu33  =  em4phi*gtu33;
    
    G111  =  Gt111 + 2*(-((-2*J11L + gt11L*gtu11*J11L + gt11L*gtu21*J12L + 
                gt11L*gtu31*J13L)*PDstandardNth1phi) - 
           (-2*J21L + gt11L*gtu11*J21L + gt11L*gtu21*J22L + gt11L*gtu31*J23L)*
            PDstandardNth2phi - (-2*J31L + gt11L*gtu11*J31L + gt11L*gtu21*J32L + 
              gt11L*gtu31*J33L)*PDstandardNth3phi);
    
    G211  =  Gt211 + 2*gt11L*(-((gtu21*J11L + gtu22*J12L + gtu32*J13L)*
              PDstandardNth1phi) - (gtu21*J21L + gtu22*J22L + gtu32*J23L)*
            PDstandardNth2phi - (gtu21*J31L + gtu22*J32L + gtu32*J33L)*
            PDstandardNth3phi);
    
    G311  =  Gt311 + 2*gt11L*(-((gtu31*J11L + gtu32*J12L + gtu33*J13L)*
              PDstandardNth1phi) - (gtu31*J21L + gtu32*J22L + gtu33*J23L)*
            PDstandardNth2phi - (gtu31*J31L + gtu32*J32L + gtu33*J33L)*
            PDstandardNth3phi);
    
    G112  =  Gt112 + 2*((J12L - gt12L*(gtu11*J11L + gtu21*J12L + gtu31*J13L))*
            PDstandardNth1phi + (J22L - 
              gt12L*(gtu11*J21L + gtu21*J22L + gtu31*J23L))*PDstandardNth2phi + 
           (J32L - gt12L*(gtu11*J31L + gtu21*J32L + gtu31*J33L))*PDstandardNth3phi);
    
    G212  =  Gt212 + 2*((J11L - gt12L*(gtu21*J11L + gtu22*J12L + gtu32*J13L))*
            PDstandardNth1phi + (J21L - 
              gt12L*(gtu21*J21L + gtu22*J22L + gtu32*J23L))*PDstandardNth2phi + 
           (J31L - gt12L*(gtu21*J31L + gtu22*J32L + gtu32*J33L))*PDstandardNth3phi);
    
    G312  =  Gt312 + 2*gt12L*(-((gtu31*J11L + gtu32*J12L + gtu33*J13L)*
              PDstandardNth1phi) - (gtu31*J21L + gtu32*J22L + gtu33*J23L)*
            PDstandardNth2phi - (gtu31*J31L + gtu32*J32L + gtu33*J33L)*
            PDstandardNth3phi);
    
    G113  =  Gt113 + 2*((J13L - gt13L*(gtu11*J11L + gtu21*J12L + gtu31*J13L))*
            PDstandardNth1phi + (J23L - 
              gt13L*(gtu11*J21L + gtu21*J22L + gtu31*J23L))*PDstandardNth2phi + 
           (J33L - gt13L*(gtu11*J31L + gtu21*J32L + gtu31*J33L))*PDstandardNth3phi);
    
    G213  =  Gt213 + 2*gt13L*(-((gtu21*J11L + gtu22*J12L + gtu32*J13L)*
              PDstandardNth1phi) - (gtu21*J21L + gtu22*J22L + gtu32*J23L)*
            PDstandardNth2phi - (gtu21*J31L + gtu22*J32L + gtu32*J33L)*
            PDstandardNth3phi);
    
    G313  =  Gt313 + 2*((J11L - gt13L*(gtu31*J11L + gtu32*J12L + gtu33*J13L))*
            PDstandardNth1phi + (J21L - 
              gt13L*(gtu31*J21L + gtu32*J22L + gtu33*J23L))*PDstandardNth2phi + 
           (J31L - gt13L*(gtu31*J31L + gtu32*J32L + gtu33*J33L))*PDstandardNth3phi);
    
    G122  =  Gt122 + 2*gt22L*(-((gtu11*J11L + gtu21*J12L + gtu31*J13L)*
              PDstandardNth1phi) - (gtu11*J21L + gtu21*J22L + gtu31*J23L)*
            PDstandardNth2phi - (gtu11*J31L + gtu21*J32L + gtu31*J33L)*
            PDstandardNth3phi);
    
    G222  =  Gt222 + 2*(-((-2*J12L + gt22L*(gtu21*J11L + gtu22*J12L + gtu32*J13L))*
              PDstandardNth1phi) - (gt22L*gtu21*J21L - 2*J22L + gt22L*gtu22*J22L + 
              gt22L*gtu32*J23L)*PDstandardNth2phi - 
           (gt22L*gtu21*J31L - 2*J32L + gt22L*gtu22*J32L + gt22L*gtu32*J33L)*
            PDstandardNth3phi);
    
    G322  =  Gt322 + 2*gt22L*(-((gtu31*J11L + gtu32*J12L + gtu33*J13L)*
              PDstandardNth1phi) - (gtu31*J21L + gtu32*J22L + gtu33*J23L)*
            PDstandardNth2phi - (gtu31*J31L + gtu32*J32L + gtu33*J33L)*
            PDstandardNth3phi);
    
    G123  =  Gt123 + 2*gt23L*(-((gtu11*J11L + gtu21*J12L + gtu31*J13L)*
              PDstandardNth1phi) - (gtu11*J21L + gtu21*J22L + gtu31*J23L)*
            PDstandardNth2phi - (gtu11*J31L + gtu21*J32L + gtu31*J33L)*
            PDstandardNth3phi);
    
    G223  =  Gt223 + 2*((J13L - gt23L*(gtu21*J11L + gtu22*J12L + gtu32*J13L))*
            PDstandardNth1phi + (J23L - 
              gt23L*(gtu21*J21L + gtu22*J22L + gtu32*J23L))*PDstandardNth2phi + 
           (J33L - gt23L*(gtu21*J31L + gtu22*J32L + gtu32*J33L))*PDstandardNth3phi);
    
    G323  =  Gt323 + 2*((J12L - gt23L*(gtu31*J11L + gtu32*J12L + gtu33*J13L))*
            PDstandardNth1phi + (J22L - 
              gt23L*(gtu31*J21L + gtu32*J22L + gtu33*J23L))*PDstandardNth2phi + 
           (J32L - gt23L*(gtu31*J31L + gtu32*J32L + gtu33*J33L))*PDstandardNth3phi);
    
    G133  =  Gt133 + 2*gt33L*(-((gtu11*J11L + gtu21*J12L + gtu31*J13L)*
              PDstandardNth1phi) - (gtu11*J21L + gtu21*J22L + gtu31*J23L)*
            PDstandardNth2phi - (gtu11*J31L + gtu21*J32L + gtu31*J33L)*
            PDstandardNth3phi);
    
    G233  =  Gt233 + 2*gt33L*(-((gtu21*J11L + gtu22*J12L + gtu32*J13L)*
              PDstandardNth1phi) - (gtu21*J21L + gtu22*J22L + gtu32*J23L)*
            PDstandardNth2phi - (gtu21*J31L + gtu22*J32L + gtu32*J33L)*
            PDstandardNth3phi);
    
    G333  =  Gt333 + 2*(-((-2*J13L + gt33L*(gtu31*J11L + gtu32*J12L + gtu33*J13L))*
              PDstandardNth1phi) - (gt33L*gtu31*J21L + gt33L*gtu32*J22L - 2*J23L + 
              gt33L*gtu33*J23L)*PDstandardNth2phi - 
           (gt33L*gtu31*J31L + gt33L*gtu32*J32L - 2*J33L + gt33L*gtu33*J33L)*
            PDstandardNth3phi);
    
    R11  =  Rphi11 + Rt11;
    
    R12  =  Rphi12 + Rt12;
    
    R13  =  Rphi13 + Rt13;
    
    R22  =  Rphi22 + Rt22;
    
    R23  =  Rphi23 + Rt23;
    
    R33  =  Rphi33 + Rt33;
    
    phirhsL  =  (beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1phi + 
        (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2phi + 
        (J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
           J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
           J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
           J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
           J33L*PDstandardNth3beta3)/6. + 
        (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3phi - 
        (alphaL*trKL)/6.;
    
    gt11rhsL  =  -2*alphaL*At11L + (beta1L*J11L + beta2L*J12L + beta3L*J13L)*
         PDstandardNth1gt11 + (beta1L*J21L + beta2L*J22L + beta3L*J23L)*
         PDstandardNth2gt11 - gt11L*ktwothird*
         (J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
           J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
           J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
           J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
           J33L*PDstandardNth3beta3) + 
        2*(J11L*(gt11L*PDstandardNth1beta1 + gt12L*PDstandardNth1beta2 + 
              gt13L*PDstandardNth1beta3) + 
           J21L*(gt11L*PDstandardNth2beta1 + gt12L*PDstandardNth2beta2 + 
              gt13L*PDstandardNth2beta3) + 
           J31L*(gt11L*PDstandardNth3beta1 + gt12L*PDstandardNth3beta2 + 
              gt13L*PDstandardNth3beta3)) + 
        (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3gt11;
    
    gt12rhsL  =  -2*alphaL*At12L + (gt12L*J11L + gt11L*J12L)*PDstandardNth1beta1 + 
        (gt22L*J11L + gt12L*J12L)*PDstandardNth1beta2 + 
        (gt23L*J11L + gt13L*J12L)*PDstandardNth1beta3 + 
        (beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1gt12 + 
        (gt12L*J21L + gt11L*J22L)*PDstandardNth2beta1 + 
        (gt22L*J21L + gt12L*J22L)*PDstandardNth2beta2 + 
        (gt23L*J21L + gt13L*J22L)*PDstandardNth2beta3 + 
        (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2gt12 + 
        (gt12L*J31L + gt11L*J32L)*PDstandardNth3beta1 + 
        (gt22L*J31L + gt12L*J32L)*PDstandardNth3beta2 + 
        (gt23L*J31L + gt13L*J32L)*PDstandardNth3beta3 - 
        gt12L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
           J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
           J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
           J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
           J33L*PDstandardNth3beta3) + 
        (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3gt12;
    
    gt13rhsL  =  -2*alphaL*At13L + (gt13L*J11L + gt11L*J13L)*PDstandardNth1beta1 + 
        (gt23L*J11L + gt12L*J13L)*PDstandardNth1beta2 + 
        (gt33L*J11L + gt13L*J13L)*PDstandardNth1beta3 + 
        (beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1gt13 + 
        (gt13L*J21L + gt11L*J23L)*PDstandardNth2beta1 + 
        (gt23L*J21L + gt12L*J23L)*PDstandardNth2beta2 + 
        (gt33L*J21L + gt13L*J23L)*PDstandardNth2beta3 + 
        (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2gt13 + 
        (gt13L*J31L + gt11L*J33L)*PDstandardNth3beta1 + 
        (gt23L*J31L + gt12L*J33L)*PDstandardNth3beta2 + 
        (gt33L*J31L + gt13L*J33L)*PDstandardNth3beta3 - 
        gt13L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
           J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
           J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
           J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
           J33L*PDstandardNth3beta3) + 
        (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3gt13;
    
    gt22rhsL  =  -2*alphaL*At22L + (beta1L*J11L + beta2L*J12L + beta3L*J13L)*
         PDstandardNth1gt22 + (beta1L*J21L + beta2L*J22L + beta3L*J23L)*
         PDstandardNth2gt22 - gt22L*ktwothird*
         (J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
           J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
           J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
           J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
           J33L*PDstandardNth3beta3) + 
        2*(J12L*(gt12L*PDstandardNth1beta1 + gt22L*PDstandardNth1beta2 + 
              gt23L*PDstandardNth1beta3) + 
           J22L*(gt12L*PDstandardNth2beta1 + gt22L*PDstandardNth2beta2 + 
              gt23L*PDstandardNth2beta3) + 
           J32L*(gt12L*PDstandardNth3beta1 + gt22L*PDstandardNth3beta2 + 
              gt23L*PDstandardNth3beta3)) + 
        (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3gt22;
    
    gt23rhsL  =  -2*alphaL*At23L + (gt13L*J12L + gt12L*J13L - gt23L*J11L*ktwothird)*
         PDstandardNth1beta1 + (gt22L*J13L + gt23L*J12L*kthird)*
         PDstandardNth1beta2 + (gt33L*J12L + gt23L*J13L*kthird)*
         PDstandardNth1beta3 + (beta1L*J11L + beta2L*J12L + beta3L*J13L)*
         PDstandardNth1gt23 + (gt13L*J22L + gt12L*J23L - gt23L*J21L*ktwothird)*
         PDstandardNth2beta1 + (gt22L*J23L + gt23L*J22L*kthird)*
         PDstandardNth2beta2 + (gt33L*J22L + gt23L*J23L*kthird)*
         PDstandardNth2beta3 + (beta1L*J21L + beta2L*J22L + beta3L*J23L)*
         PDstandardNth2gt23 + (gt13L*J32L + gt12L*J33L - gt23L*J31L*ktwothird)*
         PDstandardNth3beta1 + (gt22L*J33L + gt23L*J32L*kthird)*
         PDstandardNth3beta2 + (gt33L*J32L + gt23L*J33L*kthird)*
         PDstandardNth3beta3 + (beta1L*J31L + beta2L*J32L + beta3L*J33L)*
         PDstandardNth3gt23;
    
    gt33rhsL  =  -2*alphaL*At33L + (beta1L*J11L + beta2L*J12L + beta3L*J13L)*
         PDstandardNth1gt33 + (beta1L*J21L + beta2L*J22L + beta3L*J23L)*
         PDstandardNth2gt33 - gt33L*ktwothird*
         (J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
           J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
           J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
           J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
           J33L*PDstandardNth3beta3) + 
        2*(J13L*(gt13L*PDstandardNth1beta1 + gt23L*PDstandardNth1beta2 + 
              gt33L*PDstandardNth1beta3) + 
           J23L*(gt13L*PDstandardNth2beta1 + gt23L*PDstandardNth2beta2 + 
              gt33L*PDstandardNth2beta3) + 
           J33L*(gt13L*PDstandardNth3beta1 + gt23L*PDstandardNth3beta2 + 
              gt33L*PDstandardNth3beta3)) + 
        (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3gt33;
    
    Xt1rhsL  =  kthird*((6*gtu32*J13L*J22L + 7*gtu31*J11L*J23L)*
           PDstandardNth12beta1 + (7*gtu31*J11L*J33L + 
             6*J13L*(gtu32*J32L + gtu33*J33L))*PDstandardNth13beta1 - 
          6*((Atu11*J11L + Atu21*J12L + Atu31*J13L)*PDstandardNth1alpha + 
             (Atu11*J21L + Atu21*J22L + Atu31*J23L)*PDstandardNth2alpha + 
             Atu11*J31L*PDstandardNth3alpha) + 
          7*((gtu21*J12L + gtu31*J13L)*J21L*PDstandardNth12beta1 + 
             J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11beta1 + 
                gtu21*(J22L*PDstandardNth12beta1 + J32L*PDstandardNth13beta1)) + 
             J31L*((gtu21*J12L + gtu31*J13L)*PDstandardNth13beta1 + 
                gtu21*J22L*PDstandardNth23beta1) + 
             gtu21*(dJ112L*PDstandardNth1beta1 + dJ212L*PDstandardNth2beta1 + 
                dJ312L*PDstandardNth3beta1)) + 
          6*(gtu22*J32L*(J12L*PDstandardNth13beta1 + J22L*PDstandardNth23beta1) + 
             J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12beta1 + 
                gtu32*J22L*PDstandardNth22beta1 + gtu33*J33L*PDstandardNth23beta1)\
              + gtu32*(dJ123L*PDstandardNth1beta1 + dJ223L*PDstandardNth2beta1 + 
                J32L*(J23L*PDstandardNth23beta1 + J33L*PDstandardNth33beta1) + 
                dJ323L*PDstandardNth3beta1)) + 
          alphaL*(12*(Atu21*Gt112 + Atu31*Gt113 + Atu32*Gt123) + 
             6*(Atu11*Gt111 + Atu22*Gt122 + Atu33*Gt133) - 
             4*gtu31*J13L*PDstandardNth1trK + 
             36*((Atu11*J11L + Atu31*J13L)*PDstandardNth1phi + 
                Atu11*(J21L*PDstandardNth2phi + J31L*PDstandardNth3phi))) + 
          PDstandardNth1beta2*(dJ122L*gtu21 + dJ123L*gtu31 + 2*J12L*Xtn1) + 
          PDstandardNth1beta3*(dJ123L*gtu21 + dJ133L*gtu31 + 2*J13L*Xtn1) + 
          PDstandardNth2beta1*(7*dJ213L*gtu31 + 3*dJ233L*gtu33 - J21L*Xtn1) + 
          PDstandardNth2beta2*(dJ222L*gtu21 + dJ223L*gtu31 + 2*J22L*Xtn1) + 
          PDstandardNth2beta3*(dJ223L*gtu21 + dJ233L*gtu31 + 2*J23L*Xtn1) + 
          PDstandardNth3beta1*(7*dJ313L*gtu31 + 3*dJ333L*gtu33 - J31L*Xtn1) + 
          PDstandardNth3beta2*(dJ322L*gtu21 + dJ323L*gtu31 + 2*J32L*Xtn1) + 
          PDstandardNth3beta3*(dJ323L*gtu21 + dJ333L*gtu31 + 2*J33L*Xtn1) + 
          J12L*(gtu11*J11L*PDstandardNth11beta2 + 
             J13L*(6*gtu32*PDstandardNth11beta1 + gtu31*PDstandardNth11beta2 + 
                gtu21*PDstandardNth11beta3) + 
             J22L*(6*gtu22*PDstandardNth12beta1 + 2*gtu21*PDstandardNth12beta2) + 
             J23L*(gtu31*PDstandardNth12beta2 + gtu21*PDstandardNth12beta3) + 
             J33L*(6*gtu32*PDstandardNth13beta1 + gtu31*PDstandardNth13beta2) + 
             36*alphaL*Atu21*PDstandardNth1phi + 
             gtu21*(2*J32L*PDstandardNth13beta2 + J33L*PDstandardNth13beta3 - 
                4*alphaL*PDstandardNth1trK) + 3*beta2L*PDstandardNth1Xt1 - 
             3*PDstandardNth1beta1*Xtn2) + 
          J22L*((gtu11*J11L + gtu31*J13L)*PDstandardNth12beta2 + 
             gtu31*J23L*PDstandardNth22beta2 + 
             J33L*(6*gtu32*PDstandardNth23beta1 + gtu21*PDstandardNth23beta3) + 
             36*alphaL*Atu21*PDstandardNth2phi + 
             gtu21*(J13L*PDstandardNth12beta3 + 7*J21L*PDstandardNth22beta1 + 
                J23L*PDstandardNth22beta3 + 2*J32L*PDstandardNth23beta2 - 
                4*alphaL*PDstandardNth2trK) + 3*beta2L*PDstandardNth2Xt1 - 
             3*PDstandardNth2beta1*Xtn2) + 
          J32L*((gtu11*J11L + gtu31*J13L)*PDstandardNth13beta2 + 
             gtu11*J21L*PDstandardNth23beta2 + gtu31*J33L*PDstandardNth33beta2 + 
             Atu21*(-6*PDstandardNth3alpha + 36*alphaL*PDstandardNth3phi) + 
             gtu21*(J13L*PDstandardNth13beta3 + J23L*PDstandardNth23beta3 + 
                7*(J21L*PDstandardNth23beta1 + J31L*PDstandardNth33beta1) + 
                J33L*PDstandardNth33beta3 - 4*alphaL*PDstandardNth3trK) + 
             3*beta2L*PDstandardNth3Xt1 - 3*PDstandardNth3beta1*Xtn2) + 
          PDstandardNth1beta1*(7*dJ113L*gtu31 + 3*(dJ122L*gtu22 + dJ133L*gtu33) - 
             J11L*Xtn1 - 3*J13L*Xtn3) + 
          J23L*((gtu11*J11L + 2*gtu31*J13L)*PDstandardNth12beta3 + 
             J21L*(7*gtu31*PDstandardNth22beta1 + gtu11*PDstandardNth22beta3) + 
             36*alphaL*Atu31*PDstandardNth2phi + 
             gtu31*(7*J31L*PDstandardNth23beta1 + J32L*PDstandardNth23beta2 + 
                2*J33L*PDstandardNth23beta3 - 4*alphaL*PDstandardNth2trK) + 
             3*beta3L*PDstandardNth2Xt1 - 3*PDstandardNth2beta1*Xtn3) + 
          J33L*((gtu11*J11L + 2*gtu31*J13L)*PDstandardNth13beta3 + 
             gtu11*(J21L*PDstandardNth23beta3 + J31L*PDstandardNth33beta3) + 
             Atu31*(-6*PDstandardNth3alpha + 36*alphaL*PDstandardNth3phi) + 
             gtu31*(J22L*PDstandardNth23beta2 + 
                7*(J21L*PDstandardNth23beta1 + J31L*PDstandardNth33beta1) - 
                4*alphaL*PDstandardNth3trK) + 3*beta3L*PDstandardNth3Xt1 - 
             3*PDstandardNth3beta1*Xtn3) + 
          gtu11*(dJ112L*PDstandardNth1beta2 + dJ113L*PDstandardNth1beta3 + 
             J11L*(J13L*PDstandardNth11beta3 + 
                8*(J21L*PDstandardNth12beta1 + J31L*PDstandardNth13beta1) - 
                4*alphaL*PDstandardNth1trK) + dJ212L*PDstandardNth2beta2 + 
             dJ213L*PDstandardNth2beta3 + 
             J21L*(J12L*PDstandardNth12beta2 + J13L*PDstandardNth12beta3 + 
                J22L*PDstandardNth22beta2 - 4*alphaL*PDstandardNth2trK) + 
             dJ312L*PDstandardNth3beta2 + dJ313L*PDstandardNth3beta3 + 
             J31L*(J12L*PDstandardNth13beta2 + J13L*PDstandardNth13beta3 + 
                8*J21L*PDstandardNth23beta1 + J22L*PDstandardNth23beta2 + 
                J23L*PDstandardNth23beta3 + J32L*PDstandardNth33beta2 - 
                4*alphaL*PDstandardNth3trK) + 
             4*(dJ111L*PDstandardNth1beta1 + dJ211L*PDstandardNth2beta1 + 
                dJ311L*PDstandardNth3beta1 + PDstandardNth11beta1*SQR(J11L) + 
                PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L))) + 
          gtu21*(PDstandardNth11beta2*SQR(J12L) + PDstandardNth22beta2*SQR(J22L) + 
             PDstandardNth33beta2*SQR(J32L)) + 
          gtu31*(PDstandardNth11beta3*SQR(J13L) + PDstandardNth22beta3*SQR(J23L) + 
             PDstandardNth33beta3*SQR(J33L)) + 
          3*((beta1L*J11L + beta3L*J13L)*PDstandardNth1Xt1 + 
             beta1L*(J21L*PDstandardNth2Xt1 + J31L*PDstandardNth3Xt1) + 
             gtu22*(dJ222L*PDstandardNth2beta1 + dJ322L*PDstandardNth3beta1 + 
                PDstandardNth11beta1*SQR(J12L) + PDstandardNth22beta1*SQR(J22L) + 
                PDstandardNth33beta1*SQR(J32L)) + 
             gtu33*(PDstandardNth11beta1*SQR(J13L) + 
                PDstandardNth22beta1*SQR(J23L) + PDstandardNth33beta1*SQR(J33L))));
    
    Xt2rhsL  =  kthird*((7*gtu32*J13L*J22L + 6*gtu31*J11L*J23L)*
           PDstandardNth12beta2 + J11L*
           ((gtu22*J12L + gtu32*J13L)*PDstandardNth11beta1 + 
             gtu22*J32L*PDstandardNth13beta1 + 
             gtu32*(J23L*PDstandardNth12beta1 + J33L*PDstandardNth13beta1) + 
             6*gtu11*(J21L*PDstandardNth12beta2 + J31L*PDstandardNth13beta2) + 
             gtu21*(7*J12L*PDstandardNth11beta2 + J13L*PDstandardNth11beta3 + 
                2*(J21L*PDstandardNth12beta1 + J31L*PDstandardNth13beta1) - 
                4*alphaL*PDstandardNth1trK) + 3*beta1L*PDstandardNth1Xt2) + 
          J21L*((2*gtu21*J31L + gtu22*J32L)*PDstandardNth23beta1 + 
             6*gtu11*J31L*PDstandardNth23beta2 - 4*alphaL*gtu21*PDstandardNth2trK)\
           - 6*((Atu21*J11L + Atu22*J12L + Atu32*J13L)*PDstandardNth1alpha + 
             (Atu21*J21L + Atu22*J22L + Atu32*J23L)*PDstandardNth2alpha + 
             (Atu21*J31L + Atu22*J32L)*PDstandardNth3alpha) + 
          7*(J32L*((gtu21*J11L + gtu32*J13L)*PDstandardNth13beta2 + 
                (gtu21*J21L + gtu32*J23L)*PDstandardNth23beta2 + 
                gtu21*J31L*PDstandardNth33beta2) + 
             gtu21*((J12L*J21L + J11L*J22L)*PDstandardNth12beta2 + 
                dJ112L*PDstandardNth1beta2 + J21L*J22L*PDstandardNth22beta2 + 
                J31L*(J12L*PDstandardNth13beta2 + J22L*PDstandardNth23beta2) + 
                dJ212L*PDstandardNth2beta2 + dJ312L*PDstandardNth3beta2)) + 
          alphaL*(12*(Atu21*Gt212 + Atu31*Gt213 + Atu32*Gt223) + 
             6*(Atu11*Gt211 + Atu22*Gt222 + Atu33*Gt233) + 
             36*((Atu21*J11L + Atu22*J12L)*PDstandardNth1phi + 
                (Atu21*J21L + Atu22*J22L)*PDstandardNth2phi + 
                (Atu21*J31L + Atu22*J32L)*PDstandardNth3phi) - 
             4*gtu21*J31L*PDstandardNth3trK) + 
          PDstandardNth1beta1*(dJ111L*gtu21 + dJ113L*gtu32 + 2*J11L*Xtn2) + 
          PDstandardNth1beta2*(6*dJ113L*gtu31 + 7*dJ123L*gtu32 + 
             3*(dJ111L*gtu11 + dJ133L*gtu33) - 3*J11L*Xtn1 - J12L*Xtn2) + 
          PDstandardNth1beta3*(dJ113L*gtu21 + dJ133L*gtu32 + 2*J13L*Xtn2) + 
          PDstandardNth2beta1*(dJ211L*gtu21 + dJ213L*gtu32 + 2*J21L*Xtn2) + 
          PDstandardNth2beta2*(6*dJ213L*gtu31 + 7*dJ223L*gtu32 + 3*dJ233L*gtu33 - 
             3*J21L*Xtn1 - J22L*Xtn2) + 
          PDstandardNth2beta3*(dJ213L*gtu21 + dJ233L*gtu32 + 2*J23L*Xtn2) + 
          PDstandardNth3beta1*(dJ311L*gtu21 + dJ313L*gtu32 + 2*J31L*Xtn2) + 
          PDstandardNth3beta2*(6*dJ313L*gtu31 + 7*dJ323L*gtu32 + 3*dJ333L*gtu33 - 
             3*J31L*Xtn1 - J32L*Xtn2) + 
          PDstandardNth3beta3*(dJ313L*gtu21 + dJ333L*gtu32 + 2*J33L*Xtn2) + 
          J13L*((6*gtu31*J11L + 7*gtu32*J12L)*PDstandardNth11beta2 + 
             gtu22*J12L*PDstandardNth11beta3 + 
             J21L*(gtu32*PDstandardNth12beta1 + 6*gtu31*PDstandardNth12beta2 + 
                gtu21*PDstandardNth12beta3) + 
             J31L*(gtu32*PDstandardNth13beta1 + 6*gtu31*PDstandardNth13beta2 + 
                gtu21*PDstandardNth13beta3) + 
             alphaL*(36*Atu32*PDstandardNth1phi - 4*gtu32*PDstandardNth1trK) + 
             3*beta3L*PDstandardNth1Xt2 - 3*PDstandardNth1beta2*Xtn3) + 
          J23L*((7*gtu32*J12L + 6*gtu33*J13L)*PDstandardNth12beta2 + 
             (gtu21*J11L + gtu22*J12L + 2*gtu32*J13L)*PDstandardNth12beta3 + 
             J21L*(gtu32*PDstandardNth22beta1 + 6*gtu31*PDstandardNth22beta2 + 
                gtu21*PDstandardNth22beta3) + 
             J31L*(6*gtu31*PDstandardNth23beta2 + gtu21*PDstandardNth23beta3) + 
             36*alphaL*Atu32*PDstandardNth2phi + 
             gtu32*(7*J22L*PDstandardNth22beta2 + J31L*PDstandardNth23beta1 - 
                4*alphaL*PDstandardNth2trK) + 3*beta3L*PDstandardNth2Xt2 - 
             3*PDstandardNth2beta2*Xtn3) + 
          J33L*((7*gtu32*J12L + 6*(gtu31*J11L + gtu33*J13L))*PDstandardNth13beta2 + 
             (gtu21*J11L + gtu22*J12L + 2*gtu32*J13L)*PDstandardNth13beta3 + 
             (7*gtu32*J22L + 6*gtu33*J23L)*PDstandardNth23beta2 + 
             (gtu22*J22L + 2*gtu32*J23L)*PDstandardNth23beta3 + 
             J21L*(gtu32*PDstandardNth23beta1 + 6*gtu31*PDstandardNth23beta2 + 
                gtu21*PDstandardNth23beta3) + 
             J31L*(gtu32*PDstandardNth33beta1 + 6*gtu31*PDstandardNth33beta2 + 
                gtu21*PDstandardNth33beta3) + 
             Atu32*(-6*PDstandardNth3alpha + 36*alphaL*PDstandardNth3phi) + 
             gtu32*(7*J32L*PDstandardNth33beta2 - 4*alphaL*PDstandardNth3trK) + 
             3*beta3L*PDstandardNth3Xt2 - 3*PDstandardNth3beta2*Xtn3) + 
          gtu21*(PDstandardNth11beta1*SQR(J11L) + PDstandardNth22beta1*SQR(J21L) + 
             PDstandardNth33beta1*SQR(J31L)) + 
          gtu22*((J12L*J21L + J11L*J22L)*PDstandardNth12beta1 + 
             dJ112L*PDstandardNth1beta1 + dJ123L*PDstandardNth1beta3 + 
             J12L*(J31L*PDstandardNth13beta1 + 8*J32L*PDstandardNth13beta2 - 
                4*alphaL*PDstandardNth1trK) + dJ212L*PDstandardNth2beta1 + 
             dJ223L*PDstandardNth2beta3 + 
             J22L*(8*J12L*PDstandardNth12beta2 + J13L*PDstandardNth12beta3 + 
                J21L*PDstandardNth22beta1 + J23L*PDstandardNth22beta3 + 
                J31L*PDstandardNth23beta1 - 4*alphaL*PDstandardNth2trK) + 
             dJ312L*PDstandardNth3beta1 + dJ323L*PDstandardNth3beta3 + 
             J32L*(J13L*PDstandardNth13beta3 + 8*J22L*PDstandardNth23beta2 + 
                J23L*PDstandardNth23beta3 + J31L*PDstandardNth33beta1 + 
                J33L*PDstandardNth33beta3 - 4*alphaL*PDstandardNth3trK) + 
             4*(dJ122L*PDstandardNth1beta2 + dJ222L*PDstandardNth2beta2 + 
                dJ322L*PDstandardNth3beta2 + PDstandardNth11beta2*SQR(J12L) + 
                PDstandardNth22beta2*SQR(J22L) + PDstandardNth33beta2*SQR(J32L))) + 
          gtu32*(PDstandardNth11beta3*SQR(J13L) + PDstandardNth22beta3*SQR(J23L) + 
             PDstandardNth33beta3*SQR(J33L)) + 
          3*(beta1L*(J21L*PDstandardNth2Xt2 + J31L*PDstandardNth3Xt2) + 
             beta2L*(J12L*PDstandardNth1Xt2 + J22L*PDstandardNth2Xt2 + 
                J32L*PDstandardNth3Xt2) + 
             gtu11*(dJ211L*PDstandardNth2beta2 + dJ311L*PDstandardNth3beta2 + 
                PDstandardNth11beta2*SQR(J11L) + PDstandardNth22beta2*SQR(J21L) + 
                PDstandardNth33beta2*SQR(J31L)) + 
             gtu33*(PDstandardNth11beta2*SQR(J13L) + 
                PDstandardNth22beta2*SQR(J23L) + PDstandardNth33beta2*SQR(J33L))));
    
    Xt3rhsL  =  kthird*((gtu32*J11L*J22L + gtu33*(J13L*J21L + J11L*J23L))*
           PDstandardNth12beta1 + (6*gtu22*J12L*J22L + 
             7*(J13L*(gtu31*J21L + gtu32*J22L) + gtu31*J11L*J23L))*
           PDstandardNth12beta3 + (gtu32*J11L*J32L + gtu33*(J13L*J31L + J11L*J33L))*
           PDstandardNth13beta1 + (6*gtu22*J12L*J32L + 
             7*(J13L*(gtu31*J31L + gtu32*J32L) + gtu31*J11L*J33L))*
           PDstandardNth13beta3 + J11L*
           ((gtu32*J12L + gtu33*J13L)*PDstandardNth11beta1 + 
             6*gtu11*(J21L*PDstandardNth12beta3 + J31L*PDstandardNth13beta3) + 
             gtu31*(J12L*PDstandardNth11beta2 + 7*J13L*PDstandardNth11beta3 + 
                J22L*PDstandardNth12beta2 + 
                2*(J21L*PDstandardNth12beta1 + J31L*PDstandardNth13beta1) + 
                J32L*PDstandardNth13beta2 - 4*alphaL*PDstandardNth1trK) + 
             3*beta1L*PDstandardNth1Xt3) - 
          6*((Atu31*J11L + Atu32*J12L + Atu33*J13L)*PDstandardNth1alpha + 
             (Atu31*J21L + Atu32*J22L + Atu33*J23L)*PDstandardNth2alpha + 
             Atu31*J31L*PDstandardNth3alpha) + 
          6*((gtu11*J21L*J31L + gtu22*J22L*J32L)*PDstandardNth23beta3 + 
             gtu21*(dJ212L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta3 + 
                dJ312L*PDstandardNth3beta3)) + 
          alphaL*(12*(Atu21*Gt312 + Atu31*Gt313 + Atu32*Gt323) + 
             6*(Atu11*Gt311 + Atu22*Gt322 + Atu33*Gt333) - 
             4*gtu33*J13L*PDstandardNth1trK + 
             36*((Atu31*J11L + Atu33*J13L)*PDstandardNth1phi + 
                Atu31*(J21L*PDstandardNth2phi + J31L*PDstandardNth3phi))) + 
          PDstandardNth2beta3*(3*dJ222L*gtu22 + 7*dJ223L*gtu32 + 4*dJ233L*gtu33 - 
             3*J21L*Xtn1) + PDstandardNth3beta3*
           (3*dJ322L*gtu22 + 7*dJ323L*gtu32 + 4*dJ333L*gtu33 - 3*J31L*Xtn1) + 
          J12L*(J13L*(gtu33*PDstandardNth11beta2 + 7*gtu32*PDstandardNth11beta3) + 
             J21L*(gtu32*PDstandardNth12beta1 + gtu31*PDstandardNth12beta2) + 
             J31L*(gtu32*PDstandardNth13beta1 + gtu31*PDstandardNth13beta2) + 
             gtu33*(J23L*PDstandardNth12beta2 + J33L*PDstandardNth13beta2) + 
             6*gtu21*(J11L*PDstandardNth11beta3 + J21L*PDstandardNth12beta3 + 
                J31L*PDstandardNth13beta3) + 
             alphaL*(36*Atu32*PDstandardNth1phi - 4*gtu32*PDstandardNth1trK) + 
             3*beta2L*PDstandardNth1Xt3 - 3*PDstandardNth1beta3*Xtn2) + 
          J22L*((2*gtu32*J12L + gtu33*J13L)*PDstandardNth12beta2 + 
             J21L*(gtu32*PDstandardNth22beta1 + gtu31*PDstandardNth22beta2) + 
             6*gtu21*(J11L*PDstandardNth12beta3 + J21L*PDstandardNth22beta3 + 
                J31L*PDstandardNth23beta3) + 36*alphaL*Atu32*PDstandardNth2phi + 
             gtu32*(7*J23L*PDstandardNth22beta3 + J31L*PDstandardNth23beta1 + 
                2*J32L*PDstandardNth23beta2 - 4*alphaL*PDstandardNth2trK) + 
             3*beta2L*PDstandardNth2Xt3 - 3*PDstandardNth2beta3*Xtn2) + 
          J32L*((2*gtu32*J12L + gtu33*J13L)*PDstandardNth13beta2 + 
             J21L*(gtu32*PDstandardNth23beta1 + gtu31*PDstandardNth23beta2) + 
             6*gtu21*(J11L*PDstandardNth13beta3 + J21L*PDstandardNth23beta3) + 
             gtu33*J33L*PDstandardNth33beta2 + 
             Atu32*(-6*PDstandardNth3alpha + 36*alphaL*PDstandardNth3phi) + 
             gtu32*(J31L*PDstandardNth33beta1 + 7*J33L*PDstandardNth33beta3 - 
                4*alphaL*PDstandardNth3trK) + 3*beta2L*PDstandardNth3Xt3 - 
             3*PDstandardNth3beta3*Xtn2) + 
          PDstandardNth1beta1*(dJ111L*gtu31 + dJ112L*gtu32 + dJ113L*gtu33 + 
             2*J11L*Xtn3) + PDstandardNth1beta2*
           (dJ112L*gtu31 + dJ122L*gtu32 + dJ123L*gtu33 + 2*J12L*Xtn3) + 
          PDstandardNth1beta3*(6*dJ112L*gtu21 + 3*(dJ111L*gtu11 + dJ122L*gtu22) + 
             7*(dJ113L*gtu31 + dJ123L*gtu32) + 4*dJ133L*gtu33 - 3*J11L*Xtn1 - 
             J13L*Xtn3) + PDstandardNth2beta1*
           (dJ212L*gtu32 + dJ213L*gtu33 + 2*J21L*Xtn3) + 
          PDstandardNth2beta2*(dJ222L*gtu32 + dJ223L*gtu33 + 2*J22L*Xtn3) + 
          PDstandardNth3beta1*(dJ312L*gtu32 + dJ313L*gtu33 + 2*J31L*Xtn3) + 
          PDstandardNth3beta2*(dJ322L*gtu32 + dJ323L*gtu33 + 2*J32L*Xtn3) + 
          J23L*((7*gtu32*J12L + 8*gtu33*J13L)*PDstandardNth12beta3 + 
             7*(gtu31*J21L*PDstandardNth22beta3 + 
                gtu32*J32L*PDstandardNth23beta3) + 
             36*alphaL*Atu33*PDstandardNth2phi + 
             gtu33*(J21L*PDstandardNth22beta1 + J22L*PDstandardNth22beta2 + 
                J31L*PDstandardNth23beta1 + J32L*PDstandardNth23beta2 - 
                4*alphaL*PDstandardNth2trK) + 3*beta3L*PDstandardNth2Xt3 - 
             PDstandardNth2beta3*Xtn3) + 
          J33L*((7*gtu32*J12L + 8*gtu33*J13L)*PDstandardNth13beta3 + 
             7*((gtu31*J21L + gtu32*J22L)*PDstandardNth23beta3 + 
                gtu31*J31L*PDstandardNth33beta3) + 
             Atu33*(-6*PDstandardNth3alpha + 36*alphaL*PDstandardNth3phi) + 
             gtu33*(J21L*PDstandardNth23beta1 + J22L*PDstandardNth23beta2 + 
                8*J23L*PDstandardNth23beta3 + J31L*PDstandardNth33beta1 - 
                4*alphaL*PDstandardNth3trK) + 3*beta3L*PDstandardNth3Xt3 - 
             PDstandardNth3beta3*Xtn3) + 
          gtu31*(dJ211L*PDstandardNth2beta1 + dJ212L*PDstandardNth2beta2 + 
             J31L*(2*J21L*PDstandardNth23beta1 + J22L*PDstandardNth23beta2 + 
                7*J23L*PDstandardNth23beta3 + J32L*PDstandardNth33beta2) + 
             dJ311L*PDstandardNth3beta1 + dJ312L*PDstandardNth3beta2 + 
             7*(dJ213L*PDstandardNth2beta3 + dJ313L*PDstandardNth3beta3) - 
             4*alphaL*(J21L*PDstandardNth2trK + J31L*PDstandardNth3trK) + 
             PDstandardNth11beta1*SQR(J11L) + PDstandardNth22beta1*SQR(J21L) + 
             PDstandardNth33beta1*SQR(J31L)) + 
          gtu32*(PDstandardNth11beta2*SQR(J12L) + PDstandardNth22beta2*SQR(J22L) + 
             PDstandardNth33beta2*SQR(J32L)) + 
          3*(beta3L*J13L*PDstandardNth1Xt3 + 
             beta1L*(J21L*PDstandardNth2Xt3 + J31L*PDstandardNth3Xt3) + 
             gtu11*(dJ211L*PDstandardNth2beta3 + dJ311L*PDstandardNth3beta3 + 
                PDstandardNth11beta3*SQR(J11L) + PDstandardNth22beta3*SQR(J21L) + 
                PDstandardNth33beta3*SQR(J31L)) + 
             gtu22*(PDstandardNth11beta3*SQR(J12L) + 
                PDstandardNth22beta3*SQR(J22L) + PDstandardNth33beta3*SQR(J32L))) + 
          4*gtu33*(PDstandardNth11beta3*SQR(J13L) + 
             PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L)));
    
    trKrhsL  =  (-(dJ111L*gu11) - 2*dJ112L*gu21 - dJ122L*gu22 - 2*dJ113L*gu31 - 
           2*dJ123L*gu32 - dJ133L*gu33 + G111*gu11*J11L + 2*G112*gu21*J11L + 
           G122*gu22*J11L + 2*G113*gu31*J11L + 2*G123*gu32*J11L + G133*gu33*J11L + 
           G211*gu11*J12L + 2*G212*gu21*J12L + G222*gu22*J12L + 2*G213*gu31*J12L + 
           2*G223*gu32*J12L + G233*gu33*J12L + G311*gu11*J13L + 2*G312*gu21*J13L + 
           G322*gu22*J13L + 2*G313*gu31*J13L + 2*G323*gu32*J13L + G333*gu33*J13L)*
         PDstandardNth1alpha + (beta1L*J11L + beta2L*J12L + beta3L*J13L)*
         PDstandardNth1trK + (-(dJ211L*gu11) - 2*dJ212L*gu21 - dJ222L*gu22 - 
           2*dJ213L*gu31 - 2*dJ223L*gu32 - dJ233L*gu33 + G111*gu11*J21L + 
           2*G112*gu21*J21L + G122*gu22*J21L + 2*G113*gu31*J21L + 
           2*G123*gu32*J21L + G133*gu33*J21L + G211*gu11*J22L + 2*G212*gu21*J22L + 
           G222*gu22*J22L + 2*G213*gu31*J22L + 2*G223*gu32*J22L + G233*gu33*J22L + 
           G311*gu11*J23L + 2*G312*gu21*J23L + G322*gu22*J23L + 2*G313*gu31*J23L + 
           2*G323*gu32*J23L + G333*gu33*J23L)*PDstandardNth2alpha + 
        (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2trK - 
        dJ311L*gu11*PDstandardNth3alpha - dJ322L*gu22*PDstandardNth3alpha - 
        2*dJ313L*gu31*PDstandardNth3alpha - 2*dJ323L*gu32*PDstandardNth3alpha - 
        dJ333L*gu33*PDstandardNth3alpha + G111*gu11*J31L*PDstandardNth3alpha + 
        G122*gu22*J31L*PDstandardNth3alpha + 2*G113*gu31*J31L*PDstandardNth3alpha + 
        2*G123*gu32*J31L*PDstandardNth3alpha + G133*gu33*J31L*PDstandardNth3alpha + 
        G211*gu11*J32L*PDstandardNth3alpha + 2*G212*gu21*J32L*PDstandardNth3alpha + 
        G222*gu22*J32L*PDstandardNth3alpha + 2*G213*gu31*J32L*PDstandardNth3alpha + 
        2*G223*gu32*J32L*PDstandardNth3alpha + G233*gu33*J32L*PDstandardNth3alpha + 
        G311*gu11*J33L*PDstandardNth3alpha + 2*G312*gu21*J33L*PDstandardNth3alpha + 
        G322*gu22*J33L*PDstandardNth3alpha + 2*G313*gu31*J33L*PDstandardNth3alpha + 
        2*G323*gu32*J33L*PDstandardNth3alpha + G333*gu33*J33L*PDstandardNth3alpha - 
        2*(gu21*J11L*J12L*PDstandardNth11alpha + 
           gu31*J11L*J13L*PDstandardNth11alpha + 
           gu32*J12L*J13L*PDstandardNth11alpha + 
           gu11*J11L*J21L*PDstandardNth12alpha + 
           gu21*J12L*J21L*PDstandardNth12alpha + 
           gu31*J13L*J21L*PDstandardNth12alpha + 
           gu21*J11L*J22L*PDstandardNth12alpha + 
           gu22*J12L*J22L*PDstandardNth12alpha + 
           gu32*J13L*J22L*PDstandardNth12alpha + 
           gu31*J11L*J23L*PDstandardNth12alpha + 
           gu32*J12L*J23L*PDstandardNth12alpha + 
           gu33*J13L*J23L*PDstandardNth12alpha + 
           gu11*J11L*J31L*PDstandardNth13alpha + 
           gu21*J12L*J31L*PDstandardNth13alpha + 
           gu31*J13L*J31L*PDstandardNth13alpha + 
           gu21*J11L*J32L*PDstandardNth13alpha + 
           gu22*J12L*J32L*PDstandardNth13alpha + 
           gu32*J13L*J32L*PDstandardNth13alpha + 
           gu31*J11L*J33L*PDstandardNth13alpha + 
           gu32*J12L*J33L*PDstandardNth13alpha + 
           gu33*J13L*J33L*PDstandardNth13alpha + 
           gu21*J21L*J22L*PDstandardNth22alpha + 
           gu31*J21L*J23L*PDstandardNth22alpha + 
           gu32*J22L*J23L*PDstandardNth22alpha + 
           gu11*J21L*J31L*PDstandardNth23alpha + 
           gu21*J22L*J31L*PDstandardNth23alpha + 
           gu31*J23L*J31L*PDstandardNth23alpha + 
           gu21*J21L*J32L*PDstandardNth23alpha + 
           gu22*J22L*J32L*PDstandardNth23alpha + 
           gu32*J23L*J32L*PDstandardNth23alpha + 
           gu31*J21L*J33L*PDstandardNth23alpha + 
           gu32*J22L*J33L*PDstandardNth23alpha + 
           gu33*J23L*J33L*PDstandardNth23alpha + 
           gu21*J31L*J32L*PDstandardNth33alpha + 
           gu31*J31L*J33L*PDstandardNth33alpha + 
           gu32*J32L*J33L*PDstandardNth33alpha + dJ312L*gu21*PDstandardNth3alpha) + 
        2*(alphaL*Atm12*Atm21 + alphaL*Atm13*Atm31 + alphaL*Atm23*Atm32 + 
           G112*gu21*J31L*PDstandardNth3alpha) + beta1L*J31L*PDstandardNth3trK + 
        beta2L*J32L*PDstandardNth3trK + beta3L*J33L*PDstandardNth3trK + 
        alphaL*SQR(Atm11) + alphaL*SQR(Atm22) + alphaL*SQR(Atm33) - 
        gu11*PDstandardNth11alpha*SQR(J11L) - gu22*PDstandardNth11alpha*SQR(J12L) - 
        gu33*PDstandardNth11alpha*SQR(J13L) - gu11*PDstandardNth22alpha*SQR(J21L) - 
        gu22*PDstandardNth22alpha*SQR(J22L) - gu33*PDstandardNth22alpha*SQR(J23L) - 
        gu11*PDstandardNth33alpha*SQR(J31L) - gu22*PDstandardNth33alpha*SQR(J32L) - 
        gu33*PDstandardNth33alpha*SQR(J33L) + alphaL*kthird*SQR(trKL);
    
    Ats11  =  (-dJ111L + G111*J11L + G211*J12L + G311*J13L)*PDstandardNth1alpha - 
        2*(J11L*J21L*PDstandardNth12alpha + J11L*J31L*PDstandardNth13alpha + 
           J21L*J31L*PDstandardNth23alpha) + 
        (-dJ211L + G111*J21L + G211*J22L + G311*J23L)*PDstandardNth2alpha + 
        (-dJ311L + G111*J31L + G211*J32L + G311*J33L)*PDstandardNth3alpha + 
        alphaL*R11 - PDstandardNth11alpha*SQR(J11L) - 
        PDstandardNth22alpha*SQR(J21L) - PDstandardNth33alpha*SQR(J31L);
    
    Ats12  =  -(J11L*J12L*PDstandardNth11alpha) - J12L*J21L*PDstandardNth12alpha - 
        J11L*J22L*PDstandardNth12alpha - J12L*J31L*PDstandardNth13alpha - 
        J11L*J32L*PDstandardNth13alpha + 
        (-dJ112L + G112*J11L + G212*J12L + G312*J13L)*PDstandardNth1alpha - 
        J21L*J22L*PDstandardNth22alpha - J22L*J31L*PDstandardNth23alpha - 
        J21L*J32L*PDstandardNth23alpha + 
        (-dJ212L + G112*J21L + G212*J22L + G312*J23L)*PDstandardNth2alpha - 
        J31L*J32L*PDstandardNth33alpha - dJ312L*PDstandardNth3alpha + 
        G112*J31L*PDstandardNth3alpha + G212*J32L*PDstandardNth3alpha + 
        G312*J33L*PDstandardNth3alpha + alphaL*R12;
    
    Ats13  =  -(J11L*J13L*PDstandardNth11alpha) - J13L*J21L*PDstandardNth12alpha - 
        J11L*J23L*PDstandardNth12alpha - J13L*J31L*PDstandardNth13alpha - 
        J11L*J33L*PDstandardNth13alpha + 
        (-dJ113L + G113*J11L + G213*J12L + G313*J13L)*PDstandardNth1alpha - 
        J21L*J23L*PDstandardNth22alpha - J23L*J31L*PDstandardNth23alpha - 
        J21L*J33L*PDstandardNth23alpha + 
        (-dJ213L + G113*J21L + G213*J22L + G313*J23L)*PDstandardNth2alpha - 
        J31L*J33L*PDstandardNth33alpha - dJ313L*PDstandardNth3alpha + 
        G113*J31L*PDstandardNth3alpha + G213*J32L*PDstandardNth3alpha + 
        G313*J33L*PDstandardNth3alpha + alphaL*R13;
    
    Ats22  =  (-dJ122L + G122*J11L + G222*J12L + G322*J13L)*PDstandardNth1alpha - 
        2*(J12L*J22L*PDstandardNth12alpha + J12L*J32L*PDstandardNth13alpha + 
           J22L*J32L*PDstandardNth23alpha) + 
        (-dJ222L + G122*J21L + G222*J22L + G322*J23L)*PDstandardNth2alpha + 
        (-dJ322L + G122*J31L + G222*J32L + G322*J33L)*PDstandardNth3alpha + 
        alphaL*R22 - PDstandardNth11alpha*SQR(J12L) - 
        PDstandardNth22alpha*SQR(J22L) - PDstandardNth33alpha*SQR(J32L);
    
    Ats23  =  -(J12L*J13L*PDstandardNth11alpha) - J13L*J22L*PDstandardNth12alpha - 
        J12L*J23L*PDstandardNth12alpha - J13L*J32L*PDstandardNth13alpha - 
        J12L*J33L*PDstandardNth13alpha + 
        (-dJ123L + G123*J11L + G223*J12L + G323*J13L)*PDstandardNth1alpha - 
        J22L*J23L*PDstandardNth22alpha - J23L*J32L*PDstandardNth23alpha - 
        J22L*J33L*PDstandardNth23alpha + 
        (-dJ223L + G123*J21L + G223*J22L + G323*J23L)*PDstandardNth2alpha - 
        J32L*J33L*PDstandardNth33alpha - dJ323L*PDstandardNth3alpha + 
        G123*J31L*PDstandardNth3alpha + G223*J32L*PDstandardNth3alpha + 
        G323*J33L*PDstandardNth3alpha + alphaL*R23;
    
    Ats33  =  (-dJ133L + G133*J11L + G233*J12L + G333*J13L)*PDstandardNth1alpha - 
        2*(J13L*J23L*PDstandardNth12alpha + J13L*J33L*PDstandardNth13alpha + 
           J23L*J33L*PDstandardNth23alpha) + 
        (-dJ233L + G133*J21L + G233*J22L + G333*J23L)*PDstandardNth2alpha + 
        (-dJ333L + G133*J31L + G233*J32L + G333*J33L)*PDstandardNth3alpha + 
        alphaL*R33 - PDstandardNth11alpha*SQR(J13L) - 
        PDstandardNth22alpha*SQR(J23L) - PDstandardNth33alpha*SQR(J33L);
    
    trAts  =  Ats11*gu11 + Ats22*gu22 + 2*(Ats12*gu21 + Ats13*gu31 + Ats23*gu32) + 
        Ats33*gu33;
    
    At11rhsL  =  (beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1At11 + 
        (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2At11 + 
        (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3At11 - 
        At11L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
           J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
           J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
           J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
           J33L*PDstandardNth3beta3) + 
        2*(J11L*(At11L*PDstandardNth1beta1 + At12L*PDstandardNth1beta2 + 
              At13L*PDstandardNth1beta3) + 
           J21L*(At11L*PDstandardNth2beta1 + At12L*PDstandardNth2beta2 + 
              At13L*PDstandardNth2beta3) + 
           J31L*(At11L*PDstandardNth3beta1 + At12L*PDstandardNth3beta2 + 
              At13L*PDstandardNth3beta3)) + em4phi*(Ats11 - g11*kthird*trAts) + 
        alphaL*(-2*(At11L*Atm11 + At12L*Atm21 + At13L*Atm31) + At11L*trKL);
    
    At12rhsL  =  (beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1At12 + 
        kthird*(-6*alphaL*(At11L*Atm12 + At12L*Atm22 + At13L*Atm32) + 
           J11L*(At12L*PDstandardNth1beta1 + 3*At22L*PDstandardNth1beta2) + 
           J12L*(At12L*PDstandardNth1beta2 + 3*At13L*PDstandardNth1beta3) + 
           At12L*(J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + 
              J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 - 
              2*(J13L*PDstandardNth1beta3 + J23L*PDstandardNth2beta3 + 
                 J33L*PDstandardNth3beta3)) - em4phi*g12*trAts + 
           3*(Ats12*em4phi + (beta1L*J21L + beta2L*J22L + beta3L*J23L)*
               PDstandardNth2At12 + (beta1L*J31L + beta2L*J32L + beta3L*J33L)*
               PDstandardNth3At12 + At11L*
               (J12L*PDstandardNth1beta1 + J22L*PDstandardNth2beta1 + 
                 J32L*PDstandardNth3beta1) + 
              At22L*(J21L*PDstandardNth2beta2 + J31L*PDstandardNth3beta2) + 
              At23L*(J11L*PDstandardNth1beta3 + J21L*PDstandardNth2beta3 + 
                 J31L*PDstandardNth3beta3) + 
              At13L*(J22L*PDstandardNth2beta3 + J32L*PDstandardNth3beta3) + 
              alphaL*At12L*trKL));
    
    At13rhsL  =  (beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1At13 + 
        kthird*(-6*alphaL*(At11L*Atm13 + At12L*Atm23 + At13L*Atm33) + 
           J11L*(At13L*PDstandardNth1beta1 + 
              3*(At23L*PDstandardNth1beta2 + At33L*PDstandardNth1beta3)) + 
           At13L*(J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
              J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 - 
              2*(J12L*PDstandardNth1beta2 + J22L*PDstandardNth2beta2 + 
                 J32L*PDstandardNth3beta2) + J33L*PDstandardNth3beta3) - 
           em4phi*g13*trAts + 3*(Ats13*em4phi + 
              J13L*(At11L*PDstandardNth1beta1 + At12L*PDstandardNth1beta2) + 
              (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2At13 + 
              J23L*(At11L*PDstandardNth2beta1 + At12L*PDstandardNth2beta2) + 
              J21L*(At23L*PDstandardNth2beta2 + At33L*PDstandardNth2beta3) + 
              (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3At13 + 
              J33L*(At11L*PDstandardNth3beta1 + At12L*PDstandardNth3beta2) + 
              J31L*(At23L*PDstandardNth3beta2 + At33L*PDstandardNth3beta3) + 
              alphaL*At13L*trKL));
    
    At22rhsL  =  (beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1At22 + 
        (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2At22 + 
        (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3At22 - 
        At22L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
           J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
           J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
           J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
           J33L*PDstandardNth3beta3) + 
        2*(J12L*(At12L*PDstandardNth1beta1 + At22L*PDstandardNth1beta2 + 
              At23L*PDstandardNth1beta3) + 
           J22L*(At12L*PDstandardNth2beta1 + At22L*PDstandardNth2beta2 + 
              At23L*PDstandardNth2beta3) + 
           J32L*(At12L*PDstandardNth3beta1 + At22L*PDstandardNth3beta2 + 
              At23L*PDstandardNth3beta3)) + em4phi*(Ats22 - g22*kthird*trAts) + 
        alphaL*(-2*(At12L*Atm12 + At22L*Atm22 + At23L*Atm32) + At22L*trKL);
    
    At23rhsL  =  (beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1At23 + 
        kthird*(-6*alphaL*(At12L*Atm13 + At22L*Atm23 + At23L*Atm33) + 
           (-2*At23L*J11L + 3*At12L*J13L)*PDstandardNth1beta1 + 
           J12L*(At23L*PDstandardNth1beta2 + 3*At33L*PDstandardNth1beta3) + 
           At23L*(J13L*PDstandardNth1beta3 + J22L*PDstandardNth2beta2 + 
              J23L*PDstandardNth2beta3 - 
              2*(J21L*PDstandardNth2beta1 + J31L*PDstandardNth3beta1) + 
              J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) - 
           em4phi*g23*trAts + 3*(Ats23*em4phi + 
              (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2At23 + 
              (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3At23 + 
              At13L*(J12L*PDstandardNth1beta1 + J22L*PDstandardNth2beta1 + 
                 J32L*PDstandardNth3beta1) + 
              At12L*(J23L*PDstandardNth2beta1 + J33L*PDstandardNth3beta1) + 
              At22L*(J13L*PDstandardNth1beta2 + J23L*PDstandardNth2beta2 + 
                 J33L*PDstandardNth3beta2) + 
              At33L*(J22L*PDstandardNth2beta3 + J32L*PDstandardNth3beta3) + 
              alphaL*At23L*trKL));
    
    At33rhsL  =  (beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1At33 + 
        (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2At33 + 
        (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3At33 - 
        At33L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
           J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
           J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
           J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
           J33L*PDstandardNth3beta3) + 
        2*(J13L*(At13L*PDstandardNth1beta1 + At23L*PDstandardNth1beta2 + 
              At33L*PDstandardNth1beta3) + 
           J23L*(At13L*PDstandardNth2beta1 + At23L*PDstandardNth2beta2 + 
              At33L*PDstandardNth2beta3) + 
           J33L*(At13L*PDstandardNth3beta1 + At23L*PDstandardNth3beta2 + 
              At33L*PDstandardNth3beta3)) + em4phi*(Ats33 - g33*kthird*trAts) + 
        alphaL*(-2*(At13L*Atm13 + At23L*Atm23 + At33L*Atm33) + At33L*trKL);
    
    alpharhsL  =  LapseAdvectionCoeff*
         ((beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1alpha + 
           (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2alpha + 
           (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3alpha) + 
        harmonicF*(AL*(-1 + LapseAdvectionCoeff) - LapseAdvectionCoeff*trKL)*
         pow(alphaL,harmonicN);
    
    ArhsL  =  (-1 + LapseAdvectionCoeff)*(AL*AlphaDriver - trKrhsL);
    
    beta1rhsL  =  ((beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1beta1 + 
           (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2beta1 + 
           (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3beta1)*
         ShiftAdvectionCoeff + B1L*ShiftGammaCoeff;
    
    beta2rhsL  =  ((beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1beta2 + 
           (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2beta2 + 
           (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3beta2)*
         ShiftAdvectionCoeff + B2L*ShiftGammaCoeff;
    
    beta3rhsL  =  ((beta1L*J11L + beta2L*J12L + beta3L*J13L)*PDstandardNth1beta3 + 
           (beta1L*J21L + beta2L*J22L + beta3L*J23L)*PDstandardNth2beta3 + 
           (beta1L*J31L + beta2L*J32L + beta3L*J33L)*PDstandardNth3beta3)*
         ShiftAdvectionCoeff + B3L*ShiftGammaCoeff;
    
    B1rhsL  =  -(B1L*BetaDriver) + (beta1L*
            (J11L*(PDstandardNth1B1 - PDstandardNth1Xt1) + 
              J21L*(PDstandardNth2B1 - PDstandardNth2Xt1) + 
              J31L*(PDstandardNth3B1 - PDstandardNth3Xt1)) + 
           beta2L*(J12L*(PDstandardNth1B1 - PDstandardNth1Xt1) + 
              J22L*(PDstandardNth2B1 - PDstandardNth2Xt1) + 
              J32L*(PDstandardNth3B1 - PDstandardNth3Xt1)) + 
           beta3L*(J13L*(PDstandardNth1B1 - PDstandardNth1Xt1) + 
              J23L*(PDstandardNth2B1 - PDstandardNth2Xt1) + 
              J33L*(PDstandardNth3B1 - PDstandardNth3Xt1)))*ShiftAdvectionCoeff + 
        Xt1rhsL;
    
    B2rhsL  =  -(B2L*BetaDriver) + (beta1L*
            (J11L*(PDstandardNth1B2 - PDstandardNth1Xt2) + 
              J21L*(PDstandardNth2B2 - PDstandardNth2Xt2) + 
              J31L*(PDstandardNth3B2 - PDstandardNth3Xt2)) + 
           beta2L*(J12L*(PDstandardNth1B2 - PDstandardNth1Xt2) + 
              J22L*(PDstandardNth2B2 - PDstandardNth2Xt2) + 
              J32L*(PDstandardNth3B2 - PDstandardNth3Xt2)) + 
           beta3L*(J13L*(PDstandardNth1B2 - PDstandardNth1Xt2) + 
              J23L*(PDstandardNth2B2 - PDstandardNth2Xt2) + 
              J33L*(PDstandardNth3B2 - PDstandardNth3Xt2)))*ShiftAdvectionCoeff + 
        Xt2rhsL;
    
    B3rhsL  =  -(B3L*BetaDriver) + (beta1L*
            (J11L*(PDstandardNth1B3 - PDstandardNth1Xt3) + 
              J21L*(PDstandardNth2B3 - PDstandardNth2Xt3) + 
              J31L*(PDstandardNth3B3 - PDstandardNth3Xt3)) + 
           beta2L*(J12L*(PDstandardNth1B3 - PDstandardNth1Xt3) + 
              J22L*(PDstandardNth2B3 - PDstandardNth2Xt3) + 
              J32L*(PDstandardNth3B3 - PDstandardNth3Xt3)) + 
           beta3L*(J13L*(PDstandardNth1B3 - PDstandardNth1Xt3) + 
              J23L*(PDstandardNth2B3 - PDstandardNth2Xt3) + 
              J33L*(PDstandardNth3B3 - PDstandardNth3Xt3)))*ShiftAdvectionCoeff + 
        Xt3rhsL;
    
    
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
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_RHS_Body);
}
