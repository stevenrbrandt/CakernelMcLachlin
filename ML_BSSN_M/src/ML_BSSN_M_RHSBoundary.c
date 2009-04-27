/*  File produced by user diener */
/*  Produced with Mathematica Version 6.0 for Linux x86 (32-bit) (April 20, 2007) */

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

void ML_BSSN_M_RHSBoundary_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_M_RHSBoundary_Body");
  }
  
  if (cctk_iteration % ML_BSSN_M_RHSBoundary_calc_every != ML_BSSN_M_RHSBoundary_calc_offset)
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
  LC_LOOP3 (ML_BSSN_M_RHSBoundary,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL detgt = INITVALUE;
    CCTK_REAL em4phi = INITVALUE;
    CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL nlen = INITVALUE, nlen2 = INITVALUE;
    CCTK_REAL nn1 = INITVALUE, nn2 = INITVALUE, nn3 = INITVALUE;
    CCTK_REAL nu1 = INITVALUE, nu2 = INITVALUE, nu3 = INITVALUE;
    CCTK_REAL su1 = INITVALUE, su2 = INITVALUE, su3 = INITVALUE;
    CCTK_REAL vg = INITVALUE;
    
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
    CCTK_REAL PDMinus1A = INITVALUE;
    CCTK_REAL PDMinus2A = INITVALUE;
    CCTK_REAL PDMinus3A = INITVALUE;
    CCTK_REAL PDMinus1alpha = INITVALUE;
    CCTK_REAL PDMinus2alpha = INITVALUE;
    CCTK_REAL PDMinus3alpha = INITVALUE;
    CCTK_REAL PDMinus1At11 = INITVALUE;
    CCTK_REAL PDMinus2At11 = INITVALUE;
    CCTK_REAL PDMinus3At11 = INITVALUE;
    CCTK_REAL PDMinus1At12 = INITVALUE;
    CCTK_REAL PDMinus2At12 = INITVALUE;
    CCTK_REAL PDMinus3At12 = INITVALUE;
    CCTK_REAL PDMinus1At13 = INITVALUE;
    CCTK_REAL PDMinus2At13 = INITVALUE;
    CCTK_REAL PDMinus3At13 = INITVALUE;
    CCTK_REAL PDMinus1At22 = INITVALUE;
    CCTK_REAL PDMinus2At22 = INITVALUE;
    CCTK_REAL PDMinus3At22 = INITVALUE;
    CCTK_REAL PDMinus1At23 = INITVALUE;
    CCTK_REAL PDMinus2At23 = INITVALUE;
    CCTK_REAL PDMinus3At23 = INITVALUE;
    CCTK_REAL PDMinus1At33 = INITVALUE;
    CCTK_REAL PDMinus2At33 = INITVALUE;
    CCTK_REAL PDMinus3At33 = INITVALUE;
    CCTK_REAL PDMinus1B1 = INITVALUE;
    CCTK_REAL PDMinus2B1 = INITVALUE;
    CCTK_REAL PDMinus3B1 = INITVALUE;
    CCTK_REAL PDMinus1B2 = INITVALUE;
    CCTK_REAL PDMinus2B2 = INITVALUE;
    CCTK_REAL PDMinus3B2 = INITVALUE;
    CCTK_REAL PDMinus1B3 = INITVALUE;
    CCTK_REAL PDMinus2B3 = INITVALUE;
    CCTK_REAL PDMinus3B3 = INITVALUE;
    CCTK_REAL PDMinus1beta1 = INITVALUE;
    CCTK_REAL PDMinus2beta1 = INITVALUE;
    CCTK_REAL PDMinus3beta1 = INITVALUE;
    CCTK_REAL PDMinus1beta2 = INITVALUE;
    CCTK_REAL PDMinus2beta2 = INITVALUE;
    CCTK_REAL PDMinus3beta2 = INITVALUE;
    CCTK_REAL PDMinus1beta3 = INITVALUE;
    CCTK_REAL PDMinus2beta3 = INITVALUE;
    CCTK_REAL PDMinus3beta3 = INITVALUE;
    CCTK_REAL PDMinus1gt11 = INITVALUE;
    CCTK_REAL PDMinus2gt11 = INITVALUE;
    CCTK_REAL PDMinus3gt11 = INITVALUE;
    CCTK_REAL PDMinus1gt12 = INITVALUE;
    CCTK_REAL PDMinus2gt12 = INITVALUE;
    CCTK_REAL PDMinus3gt12 = INITVALUE;
    CCTK_REAL PDMinus1gt13 = INITVALUE;
    CCTK_REAL PDMinus2gt13 = INITVALUE;
    CCTK_REAL PDMinus3gt13 = INITVALUE;
    CCTK_REAL PDMinus1gt22 = INITVALUE;
    CCTK_REAL PDMinus2gt22 = INITVALUE;
    CCTK_REAL PDMinus3gt22 = INITVALUE;
    CCTK_REAL PDMinus1gt23 = INITVALUE;
    CCTK_REAL PDMinus2gt23 = INITVALUE;
    CCTK_REAL PDMinus3gt23 = INITVALUE;
    CCTK_REAL PDMinus1gt33 = INITVALUE;
    CCTK_REAL PDMinus2gt33 = INITVALUE;
    CCTK_REAL PDMinus3gt33 = INITVALUE;
    CCTK_REAL PDMinus1phi = INITVALUE;
    CCTK_REAL PDMinus2phi = INITVALUE;
    CCTK_REAL PDMinus3phi = INITVALUE;
    CCTK_REAL PDMinus1trK = INITVALUE;
    CCTK_REAL PDMinus2trK = INITVALUE;
    CCTK_REAL PDMinus3trK = INITVALUE;
    CCTK_REAL PDMinus1Xt1 = INITVALUE;
    CCTK_REAL PDMinus2Xt1 = INITVALUE;
    CCTK_REAL PDMinus3Xt1 = INITVALUE;
    CCTK_REAL PDMinus1Xt2 = INITVALUE;
    CCTK_REAL PDMinus2Xt2 = INITVALUE;
    CCTK_REAL PDMinus3Xt2 = INITVALUE;
    CCTK_REAL PDMinus1Xt3 = INITVALUE;
    CCTK_REAL PDMinus2Xt3 = INITVALUE;
    CCTK_REAL PDMinus3Xt3 = INITVALUE;
    
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
    Xt1L = Xt1[index];
    Xt2L = Xt2[index];
    Xt3L = Xt3[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDMinus1A = PDMinus1(A, i, j, k);
    PDMinus2A = PDMinus2(A, i, j, k);
    PDMinus3A = PDMinus3(A, i, j, k);
    PDMinus1alpha = PDMinus1(alpha, i, j, k);
    PDMinus2alpha = PDMinus2(alpha, i, j, k);
    PDMinus3alpha = PDMinus3(alpha, i, j, k);
    PDMinus1At11 = PDMinus1(At11, i, j, k);
    PDMinus2At11 = PDMinus2(At11, i, j, k);
    PDMinus3At11 = PDMinus3(At11, i, j, k);
    PDMinus1At12 = PDMinus1(At12, i, j, k);
    PDMinus2At12 = PDMinus2(At12, i, j, k);
    PDMinus3At12 = PDMinus3(At12, i, j, k);
    PDMinus1At13 = PDMinus1(At13, i, j, k);
    PDMinus2At13 = PDMinus2(At13, i, j, k);
    PDMinus3At13 = PDMinus3(At13, i, j, k);
    PDMinus1At22 = PDMinus1(At22, i, j, k);
    PDMinus2At22 = PDMinus2(At22, i, j, k);
    PDMinus3At22 = PDMinus3(At22, i, j, k);
    PDMinus1At23 = PDMinus1(At23, i, j, k);
    PDMinus2At23 = PDMinus2(At23, i, j, k);
    PDMinus3At23 = PDMinus3(At23, i, j, k);
    PDMinus1At33 = PDMinus1(At33, i, j, k);
    PDMinus2At33 = PDMinus2(At33, i, j, k);
    PDMinus3At33 = PDMinus3(At33, i, j, k);
    PDMinus1B1 = PDMinus1(B1, i, j, k);
    PDMinus2B1 = PDMinus2(B1, i, j, k);
    PDMinus3B1 = PDMinus3(B1, i, j, k);
    PDMinus1B2 = PDMinus1(B2, i, j, k);
    PDMinus2B2 = PDMinus2(B2, i, j, k);
    PDMinus3B2 = PDMinus3(B2, i, j, k);
    PDMinus1B3 = PDMinus1(B3, i, j, k);
    PDMinus2B3 = PDMinus2(B3, i, j, k);
    PDMinus3B3 = PDMinus3(B3, i, j, k);
    PDMinus1beta1 = PDMinus1(beta1, i, j, k);
    PDMinus2beta1 = PDMinus2(beta1, i, j, k);
    PDMinus3beta1 = PDMinus3(beta1, i, j, k);
    PDMinus1beta2 = PDMinus1(beta2, i, j, k);
    PDMinus2beta2 = PDMinus2(beta2, i, j, k);
    PDMinus3beta2 = PDMinus3(beta2, i, j, k);
    PDMinus1beta3 = PDMinus1(beta3, i, j, k);
    PDMinus2beta3 = PDMinus2(beta3, i, j, k);
    PDMinus3beta3 = PDMinus3(beta3, i, j, k);
    PDMinus1gt11 = PDMinus1(gt11, i, j, k);
    PDMinus2gt11 = PDMinus2(gt11, i, j, k);
    PDMinus3gt11 = PDMinus3(gt11, i, j, k);
    PDMinus1gt12 = PDMinus1(gt12, i, j, k);
    PDMinus2gt12 = PDMinus2(gt12, i, j, k);
    PDMinus3gt12 = PDMinus3(gt12, i, j, k);
    PDMinus1gt13 = PDMinus1(gt13, i, j, k);
    PDMinus2gt13 = PDMinus2(gt13, i, j, k);
    PDMinus3gt13 = PDMinus3(gt13, i, j, k);
    PDMinus1gt22 = PDMinus1(gt22, i, j, k);
    PDMinus2gt22 = PDMinus2(gt22, i, j, k);
    PDMinus3gt22 = PDMinus3(gt22, i, j, k);
    PDMinus1gt23 = PDMinus1(gt23, i, j, k);
    PDMinus2gt23 = PDMinus2(gt23, i, j, k);
    PDMinus3gt23 = PDMinus3(gt23, i, j, k);
    PDMinus1gt33 = PDMinus1(gt33, i, j, k);
    PDMinus2gt33 = PDMinus2(gt33, i, j, k);
    PDMinus3gt33 = PDMinus3(gt33, i, j, k);
    PDMinus1phi = PDMinus1(phi, i, j, k);
    PDMinus2phi = PDMinus2(phi, i, j, k);
    PDMinus3phi = PDMinus3(phi, i, j, k);
    PDMinus1trK = PDMinus1(trK, i, j, k);
    PDMinus2trK = PDMinus2(trK, i, j, k);
    PDMinus3trK = PDMinus3(trK, i, j, k);
    PDMinus1Xt1 = PDMinus1(Xt1, i, j, k);
    PDMinus2Xt1 = PDMinus2(Xt1, i, j, k);
    PDMinus3Xt1 = PDMinus3(Xt1, i, j, k);
    PDMinus1Xt2 = PDMinus1(Xt2, i, j, k);
    PDMinus2Xt2 = PDMinus2(Xt2, i, j, k);
    PDMinus3Xt2 = PDMinus3(Xt2, i, j, k);
    PDMinus1Xt3 = PDMinus1(Xt3, i, j, k);
    PDMinus2Xt3 = PDMinus2(Xt3, i, j, k);
    PDMinus3Xt3 = PDMinus3(Xt3, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    detgt  =  1;
    
    gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    gtu21  =  (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    gtu32  =  (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    em4phi  =  exp(-4*phiL);
    
    gu11  =  em4phi*gtu11;
    
    gu21  =  em4phi*gtu21;
    
    gu31  =  em4phi*gtu31;
    
    gu22  =  em4phi*gtu22;
    
    gu32  =  em4phi*gtu32;
    
    gu33  =  em4phi*gtu33;
    
    nn1  =  normal[0];
    
    nn2  =  normal[1];
    
    nn3  =  normal[2];
    
    nu1  =  gu11*nn1 + gu21*nn2 + gu31*nn3;
    
    nu2  =  gu21*nn1 + gu22*nn2 + gu32*nn3;
    
    nu3  =  gu31*nn1 + gu32*nn2 + gu33*nn3;
    
    nlen2  =  nn1*nu1 + nn2*nu2 + nn3*nu3;
    
    nlen  =  pow(nlen2,0.5);
    
    su1  =  nu1*INV(nlen);
    
    su2  =  nu2*INV(nlen);
    
    su3  =  nu3*INV(nlen);
    
    vg  =  pow(harmonicF,0.5);
    
    phirhsL  =  -((PDMinus1phi*su1 + PDMinus2phi*su2 + PDMinus3phi*su3)*vg);
    
    gt11rhsL  =  -(PDMinus1gt11*su1) - PDMinus2gt11*su2 - PDMinus3gt11*su3;
    
    gt12rhsL  =  -(PDMinus1gt12*su1) - PDMinus2gt12*su2 - PDMinus3gt12*su3;
    
    gt13rhsL  =  -(PDMinus1gt13*su1) - PDMinus2gt13*su2 - PDMinus3gt13*su3;
    
    gt22rhsL  =  -(PDMinus1gt22*su1) - PDMinus2gt22*su2 - PDMinus3gt22*su3;
    
    gt23rhsL  =  -(PDMinus1gt23*su1) - PDMinus2gt23*su2 - PDMinus3gt23*su3;
    
    gt33rhsL  =  -(PDMinus1gt33*su1) - PDMinus2gt33*su2 - PDMinus3gt33*su3;
    
    trKrhsL  =  -((PDMinus1trK*su1 + PDMinus2trK*su2 + PDMinus3trK*su3)*vg);
    
    At11rhsL  =  -(PDMinus1At11*su1) - PDMinus2At11*su2 - PDMinus3At11*su3;
    
    At12rhsL  =  -(PDMinus1At12*su1) - PDMinus2At12*su2 - PDMinus3At12*su3;
    
    At13rhsL  =  -(PDMinus1At13*su1) - PDMinus2At13*su2 - PDMinus3At13*su3;
    
    At22rhsL  =  -(PDMinus1At22*su1) - PDMinus2At22*su2 - PDMinus3At22*su3;
    
    At23rhsL  =  -(PDMinus1At23*su1) - PDMinus2At23*su2 - PDMinus3At23*su3;
    
    At33rhsL  =  -(PDMinus1At33*su1) - PDMinus2At33*su2 - PDMinus3At33*su3;
    
    Xt1rhsL  =  -(PDMinus1Xt1*su1) - PDMinus2Xt1*su2 - PDMinus3Xt1*su3;
    
    Xt2rhsL  =  -(PDMinus1Xt2*su1) - PDMinus2Xt2*su2 - PDMinus3Xt2*su3;
    
    Xt3rhsL  =  -(PDMinus1Xt3*su1) - PDMinus2Xt3*su2 - PDMinus3Xt3*su3;
    
    alpharhsL  =  -((PDMinus1alpha*su1 + PDMinus2alpha*su2 + PDMinus3alpha*su3)*vg);
    
    ArhsL  =  -((PDMinus1A*su1 + PDMinus2A*su2 + PDMinus3A*su3)*vg);
    
    beta1rhsL  =  -(PDMinus1beta1*su1) - PDMinus2beta1*su2 - PDMinus3beta1*su3;
    
    beta2rhsL  =  -(PDMinus1beta2*su1) - PDMinus2beta2*su2 - PDMinus3beta2*su3;
    
    beta3rhsL  =  -(PDMinus1beta3*su1) - PDMinus2beta3*su2 - PDMinus3beta3*su3;
    
    B1rhsL  =  -(PDMinus1B1*su1) - PDMinus2B1*su2 - PDMinus3B1*su3;
    
    B2rhsL  =  -(PDMinus1B2*su1) - PDMinus2B2*su2 - PDMinus3B2*su3;
    
    B3rhsL  =  -(PDMinus1B3*su1) - PDMinus2B3*su2 - PDMinus3B3*su3;
    
    
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
  LC_ENDLOOP3 (ML_BSSN_M_RHSBoundary);
}

void ML_BSSN_M_RHSBoundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverBoundary(cctkGH, &ML_BSSN_M_RHSBoundary_Body);
}
