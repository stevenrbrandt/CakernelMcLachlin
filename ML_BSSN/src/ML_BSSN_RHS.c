/*  File produced by user eschnett */
/*  Produced with Mathematica Version 7.0 for Mac OS X x86 (64-bit) (November 11, 2008) */

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

void ML_BSSN_RHS_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
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
  CCTK_REAL p1odx = INITVALUE;
  CCTK_REAL p1ody = INITVALUE;
  CCTK_REAL p1odz = INITVALUE;
  CCTK_REAL pm1o12dx2 = INITVALUE;
  CCTK_REAL pm1o12dy2 = INITVALUE;
  CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_RHS_Body");
  }
  
  if (cctk_iteration % ML_BSSN_RHS_calc_every != ML_BSSN_RHS_calc_offset)
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
  p1odx = INV(dx);
  p1ody = INV(dy);
  p1odz = INV(dz);
  pm1o12dx2 = -pow(dx,-2)/12.;
  pm1o12dy2 = -pow(dy,-2)/12.;
  pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_RHS,
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
    CCTK_REAL dir1 = INITVALUE, dir2 = INITVALUE, dir3 = INITVALUE;
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
    CCTK_REAL gt11L = INITVALUE, gt11rhsL = INITVALUE, gt12L = INITVALUE, gt12rhsL = INITVALUE, gt13L = INITVALUE, gt13rhsL = INITVALUE;
    CCTK_REAL gt22L = INITVALUE, gt22rhsL = INITVALUE, gt23L = INITVALUE, gt23rhsL = INITVALUE, gt33L = INITVALUE, gt33rhsL = INITVALUE;
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
    gt11L = gt11[index];
    gt12L = gt12[index];
    gt13L = gt13[index];
    gt22L = gt22[index];
    gt23L = gt23[index];
    gt33L = gt33[index];
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
    dir1  =  Sign(beta1L);
    
    dir2  =  Sign(beta2L);
    
    dir3  =  Sign(beta3L);
    
    detgt  =  1;
    
    gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    gtu21  =  (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    gtu32  =  (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    Gt111  =  khalf*(gtu11*PDstandardNth1gt11 + 2*(gtu21*PDstandardNth1gt12 + gtu31*PDstandardNth1gt13) - 
          gtu21*PDstandardNth2gt11 - gtu31*PDstandardNth3gt11);
    
    Gt211  =  khalf*(gtu21*PDstandardNth1gt11 + 2*(gtu22*PDstandardNth1gt12 + gtu32*PDstandardNth1gt13) - 
          gtu22*PDstandardNth2gt11 - gtu32*PDstandardNth3gt11);
    
    Gt311  =  khalf*(gtu31*PDstandardNth1gt11 + 2*(gtu32*PDstandardNth1gt12 + gtu33*PDstandardNth1gt13) - 
          gtu32*PDstandardNth2gt11 - gtu33*PDstandardNth3gt11);
    
    Gt112  =  khalf*(gtu21*PDstandardNth1gt22 + gtu11*PDstandardNth2gt11 + 
          gtu31*(PDstandardNth1gt23 + PDstandardNth2gt13 - PDstandardNth3gt12));
    
    Gt212  =  khalf*(gtu22*PDstandardNth1gt22 + gtu21*PDstandardNth2gt11 + 
          gtu32*(PDstandardNth1gt23 + PDstandardNth2gt13 - PDstandardNth3gt12));
    
    Gt312  =  khalf*(gtu32*PDstandardNth1gt22 + gtu31*PDstandardNth2gt11 + 
          gtu33*(PDstandardNth1gt23 + PDstandardNth2gt13 - PDstandardNth3gt12));
    
    Gt113  =  khalf*(gtu31*PDstandardNth1gt33 + gtu11*PDstandardNth3gt11 + 
          gtu21*(PDstandardNth1gt23 - PDstandardNth2gt13 + PDstandardNth3gt12));
    
    Gt213  =  khalf*(gtu32*PDstandardNth1gt33 + gtu21*PDstandardNth3gt11 + 
          gtu22*(PDstandardNth1gt23 - PDstandardNth2gt13 + PDstandardNth3gt12));
    
    Gt313  =  khalf*(gtu33*PDstandardNth1gt33 + gtu31*PDstandardNth3gt11 + 
          gtu32*(PDstandardNth1gt23 - PDstandardNth2gt13 + PDstandardNth3gt12));
    
    Gt122  =  khalf*(gtu11*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + gtu21*PDstandardNth2gt22 + 
          gtu31*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    Gt222  =  khalf*(gtu21*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + gtu22*PDstandardNth2gt22 + 
          gtu32*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    Gt322  =  khalf*(gtu31*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + gtu32*PDstandardNth2gt22 + 
          gtu33*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    Gt123  =  khalf*(gtu31*PDstandardNth2gt33 + gtu11*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) + 
          gtu21*PDstandardNth3gt22);
    
    Gt223  =  khalf*(gtu32*PDstandardNth2gt33 + gtu21*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) + 
          gtu22*PDstandardNth3gt22);
    
    Gt323  =  khalf*(gtu33*PDstandardNth2gt33 + gtu31*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) + 
          gtu32*PDstandardNth3gt22);
    
    Gt133  =  khalf*(-(gtu11*PDstandardNth1gt33) - gtu21*PDstandardNth2gt33 + 2*gtu11*PDstandardNth3gt13 + 
          2*gtu21*PDstandardNth3gt23 + gtu31*PDstandardNth3gt33);
    
    Gt233  =  khalf*(-(gtu21*PDstandardNth1gt33) - gtu22*PDstandardNth2gt33 + 2*gtu21*PDstandardNth3gt13 + 
          2*gtu22*PDstandardNth3gt23 + gtu32*PDstandardNth3gt33);
    
    Gt333  =  khalf*(-(gtu31*PDstandardNth1gt33) - gtu32*PDstandardNth2gt33 + 2*gtu31*PDstandardNth3gt13 + 
          2*gtu32*PDstandardNth3gt23 + gtu33*PDstandardNth3gt33);
    
    Xtn1  =  Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    Xtn2  =  Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    Xtn3  =  Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    Rt11  =  -(gtu11*khalf*PDstandardNth11gt11) + gtu21*
         (2*Gt211*Gt212*gt22L + 4*Gt112*gt13L*Gt311 + 2*Gt113*gt11L*Gt312 + 2*gt13L*Gt312*Gt313 + 2*gt13L*Gt211*Gt322 + 
           2*gt13L*Gt311*Gt323 + 2*Gt311*Gt312*gt33L - PDstandardNth12gt11) - gtu31*PDstandardNth13gt11 + 
        gt11L*PDstandardNth1Xt1 + gt12L*(4*Gt111*Gt212*gtu21 + 2*Gt211*Gt222*gtu21 + 2*Gt212*Gt222*gtu22 + 
           4*Gt113*Gt211*gtu31 + 4*Gt113*Gt212*gtu32 + 4*Gt113*Gt213*gtu33 + PDstandardNth1Xt2) + 
        gt13L*(4*Gt111*Gt312*gtu21 + 2*Gt212*Gt312*gtu21 + 4*Gt112*Gt312*gtu22 + 4*Gt113*Gt311*gtu31 + 
           4*Gt113*Gt312*gtu32 + 4*Gt113*Gt313*gtu33 + PDstandardNth1Xt3) - gtu22*khalf*PDstandardNth22gt11 - 
        gtu32*PDstandardNth23gt11 - gtu33*khalf*PDstandardNth33gt11 + 
        Gt111*(6*Gt113*gt11L*gtu31 + 4*gt12L*Gt213*gtu31 + gt11L*Xtn1) + 
        Gt211*(2*Gt112*gt11L*gtu11 + 4*Gt111*gt12L*gtu11 + 2*gt11L*Gt122*gtu21 + 2*gt11L*Gt123*gtu31 + gt12L*Xtn1) + 
        Gt311*(4*Gt111*gt13L*gtu11 + 2*gt12L*Gt213*gtu11 + 2*gt13L*Gt313*gtu11 + 2*gt11L*Gt123*gtu21 + 
           2*gt11L*Gt133*gtu31 + gt13L*Xtn1) + gt12L*Gt212*Xtn2 + gt13L*Gt312*Xtn2 + 
        Gt112*(6*Gt111*gt11L*gtu21 + 4*gt12L*Gt211*gtu21 + 4*gt12L*Gt212*gtu22 + 2*gt11L*Gt213*gtu31 + 
           6*Gt113*gt11L*gtu32 + gt11L*Xtn2) + Gt113*gt11L*Xtn3 + 
        Gt213*(2*gt11L*Gt122*gtu32 + 4*Gt112*gt12L*gtu32 + 2*gt11L*Gt123*gtu33 + gt12L*Xtn3) + 
        Gt313*(4*Gt111*gt13L*gtu31 + 2*gt12L*Gt213*gtu31 + 2*gt11L*Gt123*gtu32 + 4*Gt112*gt13L*gtu32 + 
           2*gt12L*Gt223*gtu32 + 2*gt11L*Gt133*gtu33 + gt13L*Xtn3) + 3*gt11L*gtu11*SQR(Gt111) + 3*gt11L*gtu22*SQR(Gt112) + 
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
    
    Rt12  =  khalf*(-(gtu11*PDstandardNth11gt12) - 2*gtu21*PDstandardNth12gt12 - 2*gtu31*PDstandardNth13gt12 + 
          gt12L*PDstandardNth1Xt1 + gt22L*PDstandardNth1Xt2 + gt23L*PDstandardNth1Xt3 - gtu22*PDstandardNth22gt12 - 
          2*gtu32*PDstandardNth23gt12 + gt11L*PDstandardNth2Xt1 + gt12L*PDstandardNth2Xt2 + gt13L*PDstandardNth2Xt3 - 
          gtu33*PDstandardNth33gt12 + (Gt111*gt12L + Gt211*gt22L + gt23L*Gt311)*Xtn1 + 
          (Gt112*gt11L + gt12L*Gt212 + gt13L*Gt312)*Xtn1 + (Gt112*gt12L + Gt212*gt22L + gt23L*Gt312)*Xtn2 + 
          (gt11L*Gt122 + gt12L*Gt222 + gt13L*Gt322)*Xtn2 + (Gt113*gt12L + Gt213*gt22L + gt23L*Gt313)*Xtn3 + 
          (gt11L*Gt123 + gt12L*Gt223 + gt13L*Gt323)*Xtn3 + 
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
    
    Rt13  =  khalf*(-(gtu11*PDstandardNth11gt13) - 2*gtu21*PDstandardNth12gt13 - 2*gtu31*PDstandardNth13gt13 + 
          gt13L*PDstandardNth1Xt1 + gt23L*PDstandardNth1Xt2 + gt33L*PDstandardNth1Xt3 - gtu22*PDstandardNth22gt13 - 
          2*gtu32*PDstandardNth23gt13 - gtu33*PDstandardNth33gt13 + gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + 
          gt13L*PDstandardNth3Xt3 + (Gt113*gt11L + gt12L*Gt213 + gt13L*Gt313)*Xtn1 + 
          (Gt111*gt13L + Gt211*gt23L + Gt311*gt33L)*Xtn1 + (gt11L*Gt123 + gt12L*Gt223 + gt13L*Gt323)*Xtn2 + 
          (Gt112*gt13L + Gt212*gt23L + Gt312*gt33L)*Xtn2 + (gt11L*Gt133 + gt12L*Gt233 + gt13L*Gt333)*Xtn3 + 
          (Gt113*gt13L + Gt213*gt23L + Gt313*gt33L)*Xtn3 + 
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
    
    Rt22  =  4*(Gt122*gt12L*Gt212*gtu21 + Gt112*gt12L*Gt222*gtu21 + Gt122*gt12L*Gt222*gtu22 + Gt123*gt12L*Gt212*gtu31 + 
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
           2*Gt323*Gt333*gtu33 + PDstandardNth2Xt3) - gtu33*khalf*PDstandardNth33gt22 + Gt212*gt22L*Xtn1 + 
        Gt112*(2*Gt111*gt12L*gtu11 + 4*gt12L*Gt212*gtu11 + 2*gt11L*Gt122*gtu21 + 2*Gt122*gt12L*gtu22 + 
           2*gt11L*Gt123*gtu31 + 2*Gt123*gt12L*gtu32 + gt12L*Xtn1) + 
        Gt312*(2*Gt213*gt22L*gtu11 + 4*Gt212*gt23L*gtu11 + 2*gt23L*Gt313*gtu11 + 2*Gt123*gt12L*gtu21 + 
           2*Gt122*gt23L*gtu22 + 2*gt12L*Gt133*gtu31 + 2*gt22L*Gt233*gtu31 + 4*Gt223*gt23L*gtu31 + 2*Gt123*gt23L*gtu32 + 
           gt23L*Xtn1) + Gt122*gt12L*Xtn2 + Gt222*gt22L*Xtn2 + gt23L*Gt322*Xtn2 + Gt123*gt12L*Xtn3 + Gt223*gt22L*Xtn3 + 
        gt23L*Gt323*Xtn3 + gt11L*gtu11*SQR(Gt112) + 2*
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
    
    Rt23  =  khalf*(-(gtu11*PDstandardNth11gt23) - 2*gtu21*PDstandardNth12gt23 - 2*gtu31*PDstandardNth13gt23 - 
          gtu22*PDstandardNth22gt23 - 2*gtu32*PDstandardNth23gt23 + gt13L*PDstandardNth2Xt1 + gt23L*PDstandardNth2Xt2 + 
          gt33L*PDstandardNth2Xt3 - gtu33*PDstandardNth33gt23 + gt12L*PDstandardNth3Xt1 + gt22L*PDstandardNth3Xt2 + 
          gt23L*PDstandardNth3Xt3 + (Gt113*gt12L + Gt213*gt22L + gt23L*Gt313)*Xtn1 + 
          (Gt112*gt13L + Gt212*gt23L + Gt312*gt33L)*Xtn1 + (Gt123*gt12L + Gt223*gt22L + gt23L*Gt323)*Xtn2 + 
          (Gt122*gt13L + Gt222*gt23L + Gt322*gt33L)*Xtn2 + (gt12L*Gt133 + gt22L*Gt233 + gt23L*Gt333)*Xtn3 + 
          (Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)*Xtn3 + 
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
    
    Rt33  =  4*(Gt123*gt13L*Gt323*gtu22 + Gt223*gt23L*Gt323*gtu22 + Gt133*gt13L*Gt313*gtu31 + Gt233*gt23L*Gt313*gtu31 + 
           Gt113*gt13L*Gt333*gtu31 + Gt133*gt13L*Gt323*gtu32 + Gt233*gt23L*Gt323*gtu32 + Gt123*gt13L*Gt333*gtu32 + 
           Gt133*gt13L*Gt333*gtu33) + gtu21*(2*Gt212*Gt223*gt23L + 4*Gt123*gt13L*Gt313 + 4*Gt223*gt23L*Gt313 + 
           4*Gt113*gt13L*Gt323 + 4*Gt213*gt23L*Gt323 + 2*Gt123*Gt311*gt33L - PDstandardNth12gt33) + 
        gtu31*(4*Gt213*gt23L*Gt333 + 2*Gt233*Gt312*gt33L + 6*Gt313*Gt333*gt33L - PDstandardNth13gt33) - 
        gtu22*khalf*PDstandardNth22gt33 + gtu32*(4*Gt223*gt23L*Gt333 + 2*Gt123*Gt313*gt33L + 6*Gt323*Gt333*gt33L - 
           PDstandardNth23gt33) - gtu33*khalf*PDstandardNth33gt33 + gt13L*PDstandardNth3Xt1 + gt23L*PDstandardNth3Xt2 + 
        gt33L*(2*Gt213*Gt322*gtu21 + 6*Gt313*Gt323*gtu21 + 2*Gt123*Gt312*gtu22 + 2*Gt133*Gt311*gtu31 + 
           2*Gt133*Gt312*gtu32 + 2*Gt133*Gt313*gtu33 + PDstandardNth3Xt3) + Gt113*gt13L*Xtn1 + Gt213*gt23L*Xtn1 + 
        Gt313*gt33L*Xtn1 + Gt123*gt13L*Xtn2 + Gt223*gt23L*Xtn2 + Gt323*gt33L*Xtn2 + Gt133*gt13L*Xtn3 + Gt333*gt33L*Xtn3 + 
        Gt233*(4*gt23L*Gt333*gtu33 + 2*Gt323*gt33L*gtu33 + gt23L*Xtn3) + 
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
    
    Rphi11  =  2*(-PDstandardNth11phi - gt11L*gtu11*PDstandardNth11phi - 2*gt11L*gtu21*PDstandardNth12phi - 
          2*gt11L*gtu31*PDstandardNth13phi - gt11L*gtu22*PDstandardNth22phi - 2*gt11L*gtu32*PDstandardNth23phi - 
          gt11L*gtu33*PDstandardNth33phi + Gt311*PDstandardNth3phi + gt11L*Gt311*gtu11*PDstandardNth3phi + 
          2*gt11L*Gt312*gtu21*PDstandardNth3phi + gt11L*Gt322*gtu22*PDstandardNth3phi + 
          2*gt11L*Gt313*gtu31*PDstandardNth3phi + 2*gt11L*Gt323*gtu32*PDstandardNth3phi + 
          gt11L*Gt333*gtu33*PDstandardNth3phi + PDstandardNth1phi*
           (Gt111 + Gt111*gt11L*gtu11 + 2*Gt112*gt11L*gtu21 + gt11L*Gt122*gtu22 + 2*Gt113*gt11L*gtu31 + 
             2*gt11L*Gt123*gtu32 + gt11L*Gt133*gtu33 - 4*gt11L*gtu21*PDstandardNth2phi - 4*gt11L*gtu31*PDstandardNth3phi) + 
          PDstandardNth2phi*(Gt211 + gt11L*Gt211*gtu11 + 
             gt11L*(2*Gt212*gtu21 + Gt222*gtu22 + 2*Gt213*gtu31 + 2*Gt223*gtu32 + Gt233*gtu33) - 
             4*gt11L*gtu32*PDstandardNth3phi) + (2 - 2*gt11L*gtu11)*SQR(PDstandardNth1phi) - 
          2*gt11L*gtu22*SQR(PDstandardNth2phi) - 2*gt11L*gtu33*SQR(PDstandardNth3phi));
    
    Rphi12  =  2*(-(gt12L*gtu11*PDstandardNth11phi) - PDstandardNth12phi - 2*gt12L*gtu21*PDstandardNth12phi - 
          2*gt12L*gtu31*PDstandardNth13phi - gt12L*gtu22*PDstandardNth22phi - 2*gt12L*gtu32*PDstandardNth23phi - 
          gt12L*gtu33*PDstandardNth33phi + Gt312*PDstandardNth3phi + gt12L*Gt311*gtu11*PDstandardNth3phi + 
          2*gt12L*Gt312*gtu21*PDstandardNth3phi + gt12L*Gt322*gtu22*PDstandardNth3phi + 
          2*gt12L*Gt313*gtu31*PDstandardNth3phi + 2*gt12L*Gt323*gtu32*PDstandardNth3phi + 
          gt12L*Gt333*gtu33*PDstandardNth3phi + PDstandardNth1phi*
           (Gt112 + Gt111*gt12L*gtu11 + 2*Gt112*gt12L*gtu21 + Gt122*gt12L*gtu22 + 2*Gt113*gt12L*gtu31 + 
             2*Gt123*gt12L*gtu32 + gt12L*Gt133*gtu33 + (2 - 4*gt12L*gtu21)*PDstandardNth2phi - 
             4*gt12L*gtu31*PDstandardNth3phi) + PDstandardNth2phi*
           (Gt212 + 2*gt12L*Gt212*gtu21 + gt12L*(Gt211*gtu11 + Gt222*gtu22 + 2*Gt213*gtu31 + 2*Gt223*gtu32 + Gt233*gtu33) - 
             4*gt12L*gtu32*PDstandardNth3phi) - 2*gt12L*gtu11*SQR(PDstandardNth1phi) - 
          2*gt12L*gtu22*SQR(PDstandardNth2phi) - 2*gt12L*gtu33*SQR(PDstandardNth3phi));
    
    Rphi13  =  2*(-PDstandardNth13phi + gt13L*(-(gtu11*PDstandardNth11phi) - 2*gtu21*PDstandardNth12phi - 
             2*gtu31*PDstandardNth13phi) - gt13L*gtu22*PDstandardNth22phi - 2*gt13L*gtu32*PDstandardNth23phi - 
          gt13L*gtu33*PDstandardNth33phi + Gt313*PDstandardNth3phi + gt13L*Gt311*gtu11*PDstandardNth3phi + 
          2*gt13L*Gt312*gtu21*PDstandardNth3phi + gt13L*Gt322*gtu22*PDstandardNth3phi + 
          2*gt13L*Gt313*gtu31*PDstandardNth3phi + 2*gt13L*Gt323*gtu32*PDstandardNth3phi + 
          gt13L*Gt333*gtu33*PDstandardNth3phi + PDstandardNth1phi*
           (Gt113 + Gt111*gt13L*gtu11 + 2*Gt112*gt13L*gtu21 + Gt122*gt13L*gtu22 + 2*Gt113*gt13L*gtu31 + 
             2*Gt123*gt13L*gtu32 + Gt133*gt13L*gtu33 - 4*gt13L*gtu21*PDstandardNth2phi + 
             (2 - 4*gt13L*gtu31)*PDstandardNth3phi) + 
          PDstandardNth2phi*(Gt213 + 2*gt13L*Gt213*gtu31 + 
             gt13L*(Gt211*gtu11 + 2*Gt212*gtu21 + Gt222*gtu22 + 2*Gt223*gtu32 + Gt233*gtu33) - 
             4*gt13L*gtu32*PDstandardNth3phi) - 2*gt13L*gtu11*SQR(PDstandardNth1phi) - 
          2*gt13L*gtu22*SQR(PDstandardNth2phi) - 2*gt13L*gtu33*SQR(PDstandardNth3phi));
    
    Rphi22  =  2*(-PDstandardNth22phi + gt22L*(-(gtu11*PDstandardNth11phi) - 2*gtu21*PDstandardNth12phi - 
             2*gtu31*PDstandardNth13phi - gtu22*PDstandardNth22phi) - 2*gt22L*gtu32*PDstandardNth23phi - 
          gt22L*gtu33*PDstandardNth33phi + Gt322*PDstandardNth3phi + gt22L*Gt311*gtu11*PDstandardNth3phi + 
          2*gt22L*Gt312*gtu21*PDstandardNth3phi + gt22L*Gt322*gtu22*PDstandardNth3phi + 
          2*gt22L*Gt313*gtu31*PDstandardNth3phi + 2*gt22L*Gt323*gtu32*PDstandardNth3phi + 
          gt22L*Gt333*gtu33*PDstandardNth3phi + PDstandardNth1phi*
           (Gt122 + Gt111*gt22L*gtu11 + 2*Gt112*gt22L*gtu21 + Gt122*gt22L*gtu22 + 2*Gt113*gt22L*gtu31 + 
             2*Gt123*gt22L*gtu32 + Gt133*gt22L*gtu33 - 4*gt22L*gtu21*PDstandardNth2phi - 4*gt22L*gtu31*PDstandardNth3phi) + 
          PDstandardNth2phi*(Gt222 + Gt222*gt22L*gtu22 + 
             gt22L*(Gt211*gtu11 + 2*Gt212*gtu21 + 2*Gt213*gtu31 + 2*Gt223*gtu32 + Gt233*gtu33) - 
             4*gt22L*gtu32*PDstandardNth3phi) - 2*gt22L*gtu11*SQR(PDstandardNth1phi) + 
          (2 - 2*gt22L*gtu22)*SQR(PDstandardNth2phi) - 2*gt22L*gtu33*SQR(PDstandardNth3phi));
    
    Rphi23  =  2*(-PDstandardNth23phi + gt23L*(-(gtu11*PDstandardNth11phi) - 2*gtu21*PDstandardNth12phi - 
             2*gtu31*PDstandardNth13phi - gtu22*PDstandardNth22phi - 2*gtu32*PDstandardNth23phi) - 
          gt23L*gtu33*PDstandardNth33phi + Gt323*PDstandardNth3phi + gt23L*Gt311*gtu11*PDstandardNth3phi + 
          2*gt23L*Gt312*gtu21*PDstandardNth3phi + gt23L*Gt322*gtu22*PDstandardNth3phi + 
          2*gt23L*Gt313*gtu31*PDstandardNth3phi + 2*gt23L*Gt323*gtu32*PDstandardNth3phi + 
          gt23L*Gt333*gtu33*PDstandardNth3phi + PDstandardNth1phi*
           (Gt123 + Gt111*gt23L*gtu11 + 2*Gt112*gt23L*gtu21 + Gt122*gt23L*gtu22 + 2*Gt113*gt23L*gtu31 + 
             2*Gt123*gt23L*gtu32 + Gt133*gt23L*gtu33 - 4*gt23L*gtu21*PDstandardNth2phi - 4*gt23L*gtu31*PDstandardNth3phi) + 
          PDstandardNth2phi*(Gt223 + 2*Gt223*gt23L*gtu32 + 
             gt23L*(Gt211*gtu11 + 2*Gt212*gtu21 + Gt222*gtu22 + 2*Gt213*gtu31 + Gt233*gtu33) + 
             (2 - 4*gt23L*gtu32)*PDstandardNth3phi) - 2*gt23L*gtu11*SQR(PDstandardNth1phi) - 
          2*gt23L*gtu22*SQR(PDstandardNth2phi) - 2*gt23L*gtu33*SQR(PDstandardNth3phi));
    
    Rphi33  =  2*(-PDstandardNth33phi + (Gt333 + gt33L*
              (Gt322*gtu22 + 2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33))*PDstandardNth3phi + 
          PDstandardNth2phi*(Gt233 + gt33L*(Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + 
                Gt233*gtu33 - 4*gtu32*PDstandardNth3phi)) + 
          PDstandardNth1phi*(Gt133 + gt33L*(Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + 
                Gt133*gtu33 - 4*(gtu21*PDstandardNth2phi + gtu31*PDstandardNth3phi))) + 2*SQR(PDstandardNth3phi) + 
          gt33L*(-(gtu11*PDstandardNth11phi) - 2*gtu21*PDstandardNth12phi - 2*gtu31*PDstandardNth13phi - 
             gtu22*PDstandardNth22phi - 2*gtu32*PDstandardNth23phi - gtu33*PDstandardNth33phi + 
             Gt311*gtu11*PDstandardNth3phi - 2*gtu11*SQR(PDstandardNth1phi) - 2*gtu22*SQR(PDstandardNth2phi) - 
             2*gtu33*SQR(PDstandardNth3phi)));
    
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
    
    G111  =  Gt111 + (4 - 2*gt11L*gtu11)*PDstandardNth1phi - 2*gt11L*(gtu21*PDstandardNth2phi + gtu31*PDstandardNth3phi);
    
    G211  =  Gt211 - 2*gt11L*(gtu21*PDstandardNth1phi + gtu22*PDstandardNth2phi + gtu32*PDstandardNth3phi);
    
    G311  =  Gt311 - 2*gt11L*(gtu31*PDstandardNth1phi + gtu32*PDstandardNth2phi + gtu33*PDstandardNth3phi);
    
    G112  =  Gt112 + (2 - 2*gt12L*gtu21)*PDstandardNth2phi - 2*gt12L*(gtu11*PDstandardNth1phi + gtu31*PDstandardNth3phi);
    
    G212  =  Gt212 + (2 - 2*gt12L*gtu21)*PDstandardNth1phi - 2*gt12L*(gtu22*PDstandardNth2phi + gtu32*PDstandardNth3phi);
    
    G312  =  Gt312 - 2*gt12L*(gtu31*PDstandardNth1phi + gtu32*PDstandardNth2phi + gtu33*PDstandardNth3phi);
    
    G113  =  Gt113 + 2*PDstandardNth3phi - 2*gt13L*(gtu11*PDstandardNth1phi + gtu21*PDstandardNth2phi + 
           gtu31*PDstandardNth3phi);
    
    G213  =  Gt213 - 2*gt13L*(gtu21*PDstandardNth1phi + gtu22*PDstandardNth2phi + gtu32*PDstandardNth3phi);
    
    G313  =  Gt313 + (2 - 2*gt13L*gtu31)*PDstandardNth1phi - 2*gt13L*(gtu32*PDstandardNth2phi + gtu33*PDstandardNth3phi);
    
    G122  =  Gt122 - 2*gt22L*(gtu11*PDstandardNth1phi + gtu21*PDstandardNth2phi + gtu31*PDstandardNth3phi);
    
    G222  =  Gt222 + (4 - 2*gt22L*gtu22)*PDstandardNth2phi - 2*gt22L*(gtu21*PDstandardNth1phi + gtu32*PDstandardNth3phi);
    
    G322  =  Gt322 - 2*gt22L*(gtu31*PDstandardNth1phi + gtu32*PDstandardNth2phi + gtu33*PDstandardNth3phi);
    
    G123  =  Gt123 - 2*gt23L*(gtu11*PDstandardNth1phi + gtu21*PDstandardNth2phi + gtu31*PDstandardNth3phi);
    
    G223  =  Gt223 + 2*PDstandardNth3phi - 2*gt23L*(gtu21*PDstandardNth1phi + gtu22*PDstandardNth2phi + 
           gtu32*PDstandardNth3phi);
    
    G323  =  Gt323 + (2 - 2*gt23L*gtu32)*PDstandardNth2phi - 2*gt23L*(gtu31*PDstandardNth1phi + gtu33*PDstandardNth3phi);
    
    G133  =  Gt133 - 2*gt33L*(gtu11*PDstandardNth1phi + gtu21*PDstandardNth2phi + gtu31*PDstandardNth3phi);
    
    G233  =  Gt233 - 2*gt33L*(gtu21*PDstandardNth1phi + gtu22*PDstandardNth2phi + gtu32*PDstandardNth3phi);
    
    G333  =  Gt333 + 4*PDstandardNth3phi - 2*gt33L*(gtu31*PDstandardNth1phi + gtu32*PDstandardNth2phi + 
           gtu33*PDstandardNth3phi);
    
    R11  =  Rphi11 + Rt11;
    
    R12  =  Rphi12 + Rt12;
    
    R13  =  Rphi13 + Rt13;
    
    R22  =  Rphi22 + Rt22;
    
    R23  =  Rphi23 + Rt23;
    
    R33  =  Rphi33 + Rt33;
    
    phirhsL  =  (6*(PDupwindNth1(phi, i, j, k)*beta1L + PDupwindNth2(phi, i, j, k)*beta2L + 
             PDupwindNth3(phi, i, j, k)*beta3L) + PDstandardNth1beta1 + PDstandardNth2beta2 + PDstandardNth3beta3 - 
          alphaL*trKL)/6.;
    
    gt11rhsL  =  -2*alphaL*At11L + PDupwindNth1(gt11, i, j, k)*beta1L + PDupwindNth2(gt11, i, j, k)*beta2L + 
        PDupwindNth3(gt11, i, j, k)*beta3L + 2*(gt12L*PDstandardNth1beta2 + gt13L*PDstandardNth1beta3) + 
        gt11L*(kfourthird*PDstandardNth1beta1 - ktwothird*PDstandardNth2beta2 - ktwothird*PDstandardNth3beta3);
    
    gt12rhsL  =  -2*alphaL*At12L + PDupwindNth1(gt12, i, j, k)*beta1L + PDupwindNth2(gt12, i, j, k)*beta2L + 
        PDupwindNth3(gt12, i, j, k)*beta3L + gt22L*PDstandardNth1beta2 + gt23L*PDstandardNth1beta3 + 
        gt11L*PDstandardNth2beta1 + gt13L*PDstandardNth2beta3 + 
        gt12L*(kthird*(PDstandardNth1beta1 + PDstandardNth2beta2) - ktwothird*PDstandardNth3beta3);
    
    gt13rhsL  =  -2*alphaL*At13L + PDupwindNth1(gt13, i, j, k)*beta1L + PDupwindNth2(gt13, i, j, k)*beta2L + 
        PDupwindNth3(gt13, i, j, k)*beta3L + gt23L*PDstandardNth1beta2 + gt33L*PDstandardNth1beta3 + 
        gt11L*PDstandardNth3beta1 + gt12L*PDstandardNth3beta2 + 
        gt13L*(-(ktwothird*PDstandardNth2beta2) + kthird*(PDstandardNth1beta1 + PDstandardNth3beta3));
    
    gt22rhsL  =  -2*alphaL*At22L + PDupwindNth1(gt22, i, j, k)*beta1L + PDupwindNth2(gt22, i, j, k)*beta2L + 
        PDupwindNth3(gt22, i, j, k)*beta3L + 2*(gt12L*PDstandardNth2beta1 + gt23L*PDstandardNth2beta3) + 
        gt22L*(-(ktwothird*PDstandardNth1beta1) + kfourthird*PDstandardNth2beta2 - ktwothird*PDstandardNth3beta3);
    
    gt23rhsL  =  -2*alphaL*At23L + PDupwindNth1(gt23, i, j, k)*beta1L + PDupwindNth2(gt23, i, j, k)*beta2L + 
        PDupwindNth3(gt23, i, j, k)*beta3L + gt13L*PDstandardNth2beta1 + gt33L*PDstandardNth2beta3 + 
        gt12L*PDstandardNth3beta1 + gt22L*PDstandardNth3beta2 + 
        gt23L*(-(ktwothird*PDstandardNth1beta1) + kthird*(PDstandardNth2beta2 + PDstandardNth3beta3));
    
    gt33rhsL  =  -2*alphaL*At33L + PDupwindNth1(gt33, i, j, k)*beta1L + PDupwindNth2(gt33, i, j, k)*beta2L + 
        PDupwindNth3(gt33, i, j, k)*beta3L - gt33L*ktwothird*PDstandardNth1beta1 - gt33L*ktwothird*PDstandardNth2beta2 + 
        2*gt13L*PDstandardNth3beta1 + 2*gt23L*PDstandardNth3beta2 + gt33L*kfourthird*PDstandardNth3beta3;
    
    Xt1rhsL  =  kthird*(7*(gtu21*PDstandardNth12beta1 + gtu31*PDstandardNth13beta1) + 
          gtu11*(4*PDstandardNth11beta1 + PDstandardNth12beta2 + PDstandardNth13beta3) + 
          gtu21*(PDstandardNth22beta2 + PDstandardNth23beta3) + 
          3*(PDupwindNth1(Xt1, i, j, k)*beta1L + PDupwindNth2(Xt1, i, j, k)*beta2L + 
             PDupwindNth3(Xt1, i, j, k)*beta3L + gtu22*PDstandardNth22beta1 + gtu33*PDstandardNth33beta1) + 
          gtu31*(PDstandardNth23beta2 + PDstandardNth33beta3) - 
          6*(Atu11*PDstandardNth1alpha + Atu21*PDstandardNth2alpha + Atu31*PDstandardNth3alpha) + 
          6*(gtu32*PDstandardNth23beta1 + alphaL*(Atu11*Gt111 + Atu22*Gt122 + 2*(Atu21*Gt112 + Atu31*Gt113 + Atu32*Gt123) + 
                Atu33*Gt133 + 6*(Atu11*PDstandardNth1phi + Atu21*PDstandardNth2phi + Atu31*PDstandardNth3phi) - 
                ktwothird*(gtu11*PDstandardNth1trK + gtu21*PDstandardNth2trK + gtu31*PDstandardNth3trK))) + 
          (-3*PDstandardNth1beta1 + 2*(PDstandardNth1beta1 + PDstandardNth2beta2 + PDstandardNth3beta3))*Xtn1 - 
          3*(PDstandardNth2beta1*Xtn2 + PDstandardNth3beta1*Xtn3));
    
    Xt2rhsL  =  kthird*(gtu21*(PDstandardNth11beta1 + 7*PDstandardNth12beta2 + PDstandardNth13beta3) + 
          gtu22*(PDstandardNth12beta1 + 4*PDstandardNth22beta2 + PDstandardNth23beta3) + 
          3*(PDupwindNth1(Xt2, i, j, k)*beta1L + PDupwindNth2(Xt2, i, j, k)*beta2L + 
             PDupwindNth3(Xt2, i, j, k)*beta3L + gtu11*PDstandardNth11beta2 + gtu33*PDstandardNth33beta2) + 
          gtu32*(PDstandardNth13beta1 + 7*PDstandardNth23beta2 + PDstandardNth33beta3) - 
          6*(Atu21*PDstandardNth1alpha + Atu22*PDstandardNth2alpha + Atu32*PDstandardNth3alpha) + 
          6*(gtu31*PDstandardNth13beta2 + alphaL*(Atu11*Gt211 + Atu22*Gt222 + 2*(Atu21*Gt212 + Atu31*Gt213 + Atu32*Gt223) + 
                Atu33*Gt233 + 6*(Atu21*PDstandardNth1phi + Atu22*PDstandardNth2phi + Atu32*PDstandardNth3phi) - 
                ktwothird*(gtu21*PDstandardNth1trK + gtu22*PDstandardNth2trK + gtu32*PDstandardNth3trK))) + 
          2*(PDstandardNth1beta1 + PDstandardNth2beta2 + PDstandardNth3beta3)*Xtn2 - 
          3*(PDstandardNth1beta2*Xtn1 + PDstandardNth2beta2*Xtn2 + PDstandardNth3beta2*Xtn3));
    
    Xt3rhsL  =  kthird*(gtu31*(PDstandardNth11beta1 + PDstandardNth12beta2 + 7*PDstandardNth13beta3) + 
          3*(PDupwindNth1(Xt3, i, j, k)*beta1L + PDupwindNth2(Xt3, i, j, k)*beta2L + 
             PDupwindNth3(Xt3, i, j, k)*beta3L + gtu11*PDstandardNth11beta3 + gtu22*PDstandardNth22beta3) + 
          gtu32*(PDstandardNth12beta1 + PDstandardNth22beta2 + 7*PDstandardNth23beta3) + 
          gtu33*(PDstandardNth13beta1 + PDstandardNth23beta2 + 4*PDstandardNth33beta3) - 
          6*(Atu31*PDstandardNth1alpha + Atu32*PDstandardNth2alpha + Atu33*PDstandardNth3alpha) + 
          6*(gtu21*PDstandardNth12beta3 + alphaL*(Atu11*Gt311 + Atu22*Gt322 + 2*(Atu21*Gt312 + Atu31*Gt313 + Atu32*Gt323) + 
                Atu33*Gt333 + 6*(Atu31*PDstandardNth1phi + Atu32*PDstandardNth2phi + Atu33*PDstandardNth3phi) - 
                ktwothird*(gtu31*PDstandardNth1trK + gtu32*PDstandardNth2trK + gtu33*PDstandardNth3trK))) + 
          2*(PDstandardNth1beta1 + PDstandardNth2beta2 + PDstandardNth3beta3)*Xtn3 - 
          3*(PDstandardNth1beta3*Xtn1 + PDstandardNth2beta3*Xtn2 + PDstandardNth3beta3*Xtn3));
    
    trKrhsL  =  PDupwindNth1(trK, i, j, k)*beta1L + PDupwindNth2(trK, i, j, k)*beta2L + 
        PDupwindNth3(trK, i, j, k)*beta3L - gu11*PDstandardNth11alpha - 2*gu21*PDstandardNth12alpha - 
        2*gu31*PDstandardNth13alpha + (G111*gu11 + 2*G112*gu21 + G122*gu22 + 2*G113*gu31 + 2*G123*gu32 + G133*gu33)*
         PDstandardNth1alpha - gu22*PDstandardNth22alpha - 2*gu32*PDstandardNth23alpha + 
        (G211*gu11 + 2*G212*gu21 + G222*gu22 + 2*G213*gu31 + 2*G223*gu32 + G233*gu33)*PDstandardNth2alpha - 
        gu33*PDstandardNth33alpha + G311*gu11*PDstandardNth3alpha + G322*gu22*PDstandardNth3alpha + 
        2*G313*gu31*PDstandardNth3alpha + 2*G323*gu32*PDstandardNth3alpha + G333*gu33*PDstandardNth3alpha + 
        2*(alphaL*Atm12*Atm21 + alphaL*Atm13*Atm31 + alphaL*Atm23*Atm32 + G312*gu21*PDstandardNth3alpha) + 
        alphaL*SQR(Atm11) + alphaL*SQR(Atm22) + alphaL*SQR(Atm33) + alphaL*kthird*SQR(trKL);
    
    Ats11  =  -PDstandardNth11alpha + G111*PDstandardNth1alpha + G211*PDstandardNth2alpha + G311*PDstandardNth3alpha + 
        alphaL*R11;
    
    Ats12  =  -PDstandardNth12alpha + G112*PDstandardNth1alpha + G212*PDstandardNth2alpha + G312*PDstandardNth3alpha + 
        alphaL*R12;
    
    Ats13  =  -PDstandardNth13alpha + G113*PDstandardNth1alpha + G213*PDstandardNth2alpha + G313*PDstandardNth3alpha + 
        alphaL*R13;
    
    Ats22  =  G122*PDstandardNth1alpha - PDstandardNth22alpha + G222*PDstandardNth2alpha + G322*PDstandardNth3alpha + 
        alphaL*R22;
    
    Ats23  =  G123*PDstandardNth1alpha - PDstandardNth23alpha + G223*PDstandardNth2alpha + G323*PDstandardNth3alpha + 
        alphaL*R23;
    
    Ats33  =  G133*PDstandardNth1alpha + G233*PDstandardNth2alpha - PDstandardNth33alpha + G333*PDstandardNth3alpha + 
        alphaL*R33;
    
    trAts  =  Ats11*gu11 + Ats22*gu22 + 2*(Ats12*gu21 + Ats13*gu31 + Ats23*gu32) + Ats33*gu33;
    
    At11rhsL  =  PDupwindNth1(At11, i, j, k)*beta1L + PDupwindNth2(At11, i, j, k)*beta2L + 
        PDupwindNth3(At11, i, j, k)*beta3L + 2*(At11L*PDstandardNth1beta1 + At12L*PDstandardNth1beta2 + 
           At13L*PDstandardNth1beta3) - At11L*ktwothird*(PDstandardNth1beta1 + PDstandardNth2beta2 + PDstandardNth3beta3) + 
        em4phi*(Ats11 - g11*kthird*trAts) + alphaL*(-2*(At11L*Atm11 + At12L*Atm21 + At13L*Atm31) + At11L*trKL);
    
    At12rhsL  =  -2*alphaL*(At11L*Atm12 + At12L*Atm22 + At13L*Atm32) + PDupwindNth1(At12, i, j, k)*beta1L + 
        PDupwindNth2(At12, i, j, k)*beta2L + PDupwindNth3(At12, i, j, k)*beta3L + At22L*PDstandardNth1beta2 + 
        At23L*PDstandardNth1beta3 + At11L*PDstandardNth2beta1 + At13L*PDstandardNth2beta3 + 
        em4phi*(Ats12 - g12*kthird*trAts) + At12L*(kthird*(PDstandardNth1beta1 + PDstandardNth2beta2) - 
           ktwothird*PDstandardNth3beta3 + alphaL*trKL);
    
    At13rhsL  =  -2*alphaL*(At11L*Atm13 + At12L*Atm23 + At13L*Atm33) + PDupwindNth1(At13, i, j, k)*beta1L + 
        PDupwindNth2(At13, i, j, k)*beta2L + PDupwindNth3(At13, i, j, k)*beta3L + At23L*PDstandardNth1beta2 + 
        At33L*PDstandardNth1beta3 + At11L*PDstandardNth3beta1 + At12L*PDstandardNth3beta2 + 
        em4phi*(Ats13 - g13*kthird*trAts) + At13L*(-(ktwothird*PDstandardNth2beta2) + 
           kthird*(PDstandardNth1beta1 + PDstandardNth3beta3) + alphaL*trKL);
    
    At22rhsL  =  PDupwindNth1(At22, i, j, k)*beta1L + PDupwindNth2(At22, i, j, k)*beta2L + 
        PDupwindNth3(At22, i, j, k)*beta3L + 2*(At12L*PDstandardNth2beta1 + At22L*PDstandardNth2beta2 + 
           At23L*PDstandardNth2beta3) - At22L*ktwothird*(PDstandardNth1beta1 + PDstandardNth2beta2 + PDstandardNth3beta3) + 
        em4phi*(Ats22 - g22*kthird*trAts) + alphaL*(-2*(At12L*Atm12 + At22L*Atm22 + At23L*Atm32) + At22L*trKL);
    
    At23rhsL  =  -2*alphaL*(At12L*Atm13 + At22L*Atm23 + At23L*Atm33) + PDupwindNth1(At23, i, j, k)*beta1L + 
        PDupwindNth2(At23, i, j, k)*beta2L + PDupwindNth3(At23, i, j, k)*beta3L + At13L*PDstandardNth2beta1 + 
        At33L*PDstandardNth2beta3 + At12L*PDstandardNth3beta1 + At22L*PDstandardNth3beta2 + 
        em4phi*(Ats23 - g23*kthird*trAts) + At23L*(-(ktwothird*PDstandardNth1beta1) + 
           kthird*(PDstandardNth2beta2 + PDstandardNth3beta3) + alphaL*trKL);
    
    At33rhsL  =  PDupwindNth1(At33, i, j, k)*beta1L + PDupwindNth2(At33, i, j, k)*beta2L + 
        PDupwindNth3(At33, i, j, k)*beta3L - At33L*ktwothird*
         (PDstandardNth1beta1 + PDstandardNth2beta2 + PDstandardNth3beta3) + 
        2*(At13L*PDstandardNth3beta1 + At23L*PDstandardNth3beta2 + At33L*PDstandardNth3beta3) + 
        em4phi*(Ats33 - g33*kthird*trAts) + alphaL*(-2*(At13L*Atm13 + At23L*Atm23 + At33L*Atm33) + At33L*trKL);
    
    alpharhsL  =  (PDupwindNth1(alpha, i, j, k)*beta1L + PDupwindNth2(alpha, i, j, k)*beta2L + 
           PDupwindNth3(alpha, i, j, k)*beta3L)*LapseAdvectionCoeff + 
        harmonicF*(AL*(-1 + LapseAdvectionCoeff) - LapseAdvectionCoeff*trKL)*pow(alphaL,harmonicN);
    
    ArhsL  =  (-1 + LapseAdvectionCoeff)*(AL*AlphaDriver - trKrhsL);
    
    beta1rhsL  =  (PDupwindNth1(beta1, i, j, k)*beta1L + PDupwindNth2(beta1, i, j, k)*beta2L + 
           PDupwindNth3(beta1, i, j, k)*beta3L)*ShiftAdvectionCoeff + B1L*ShiftGammaCoeff;
    
    beta2rhsL  =  (PDupwindNth1(beta2, i, j, k)*beta1L + PDupwindNth2(beta2, i, j, k)*beta2L + 
           PDupwindNth3(beta2, i, j, k)*beta3L)*ShiftAdvectionCoeff + B2L*ShiftGammaCoeff;
    
    beta3rhsL  =  (PDupwindNth1(beta3, i, j, k)*beta1L + PDupwindNth2(beta3, i, j, k)*beta2L + 
           PDupwindNth3(beta3, i, j, k)*beta3L)*ShiftAdvectionCoeff + B3L*ShiftGammaCoeff;
    
    B1rhsL  =  -(B1L*BetaDriver) + ((PDupwindNth1(B1, i, j, k) - PDupwindNth1(Xt1, i, j, k))*beta1L + 
           (PDupwindNth2(B1, i, j, k) - PDupwindNth2(Xt1, i, j, k))*beta2L + 
           (PDupwindNth3(B1, i, j, k) - PDupwindNth3(Xt1, i, j, k))*beta3L)*ShiftAdvectionCoeff + Xt1rhsL;
    
    B2rhsL  =  -(B2L*BetaDriver) + ((PDupwindNth1(B2, i, j, k) - PDupwindNth1(Xt2, i, j, k))*beta1L + 
           (PDupwindNth2(B2, i, j, k) - PDupwindNth2(Xt2, i, j, k))*beta2L + 
           (PDupwindNth3(B2, i, j, k) - PDupwindNth3(Xt2, i, j, k))*beta3L)*ShiftAdvectionCoeff + Xt2rhsL;
    
    B3rhsL  =  -(B3L*BetaDriver) + ((PDupwindNth1(B3, i, j, k) - PDupwindNth1(Xt3, i, j, k))*beta1L + 
           (PDupwindNth2(B3, i, j, k) - PDupwindNth2(Xt3, i, j, k))*beta2L + 
           (PDupwindNth3(B3, i, j, k) - PDupwindNth3(Xt3, i, j, k))*beta3L)*ShiftAdvectionCoeff + Xt3rhsL;
    
    
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
  LC_ENDLOOP3 (ML_BSSN_RHS);
}

void ML_BSSN_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_RHS_Body);
}
