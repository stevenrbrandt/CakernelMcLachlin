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

void ML_BSSN_MP_convertFromADMBaseGamma_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_convertFromADMBaseGamma_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_convertFromADMBaseGamma_calc_every != ML_BSSN_MP_convertFromADMBaseGamma_calc_offset)
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
  LC_LOOP3 (ML_BSSN_MP_convertFromADMBaseGamma,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lssh[CCTK_LSSH_IDX(0,0)],cctk_lssh[CCTK_LSSH_IDX(0,1)],cctk_lssh[CCTK_LSSH_IDX(0,2)])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL detgt = INITVALUE;
    CCTK_REAL Gt111 = INITVALUE, Gt112 = INITVALUE, Gt113 = INITVALUE, Gt122 = INITVALUE, Gt123 = INITVALUE, Gt133 = INITVALUE;
    CCTK_REAL Gt211 = INITVALUE, Gt212 = INITVALUE, Gt213 = INITVALUE, Gt222 = INITVALUE, Gt223 = INITVALUE, Gt233 = INITVALUE;
    CCTK_REAL Gt311 = INITVALUE, Gt312 = INITVALUE, Gt313 = INITVALUE, Gt322 = INITVALUE, Gt323 = INITVALUE, Gt333 = INITVALUE;
    CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL AL = INITVALUE;
    CCTK_REAL alphaL = INITVALUE;
    CCTK_REAL B1L = INITVALUE, B2L = INITVALUE, B3L = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    CCTK_REAL dtalpL = INITVALUE;
    CCTK_REAL dtbetaxL = INITVALUE;
    CCTK_REAL dtbetayL = INITVALUE;
    CCTK_REAL dtbetazL = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt12L = INITVALUE, gt13L = INITVALUE, gt22L = INITVALUE, gt23L = INITVALUE, gt33L = INITVALUE;
    CCTK_REAL J11L = INITVALUE, J12L = INITVALUE, J13L = INITVALUE, J21L = INITVALUE, J22L = INITVALUE, J23L = INITVALUE;
    CCTK_REAL J31L = INITVALUE, J32L = INITVALUE, J33L = INITVALUE;
    CCTK_REAL Xt1L = INITVALUE, Xt2L = INITVALUE, Xt3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandardNth1beta1 = INITVALUE;
    CCTK_REAL PDstandardNth2beta1 = INITVALUE;
    CCTK_REAL PDstandardNth3beta1 = INITVALUE;
    CCTK_REAL PDstandardNth1beta2 = INITVALUE;
    CCTK_REAL PDstandardNth2beta2 = INITVALUE;
    CCTK_REAL PDstandardNth3beta2 = INITVALUE;
    CCTK_REAL PDstandardNth1beta3 = INITVALUE;
    CCTK_REAL PDstandardNth2beta3 = INITVALUE;
    CCTK_REAL PDstandardNth3beta3 = INITVALUE;
    CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    
    /* Assign local copies of grid functions */
    alphaL = alpha[index];
    beta1L = beta1[index];
    beta2L = beta2[index];
    beta3L = beta3[index];
    dtalpL = dtalp[index];
    dtbetaxL = dtbetax[index];
    dtbetayL = dtbetay[index];
    dtbetazL = dtbetaz[index];
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
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandardNth1beta1 = PDstandardNth1(beta1, i, j, k);
    PDstandardNth2beta1 = PDstandardNth2(beta1, i, j, k);
    PDstandardNth3beta1 = PDstandardNth3(beta1, i, j, k);
    PDstandardNth1beta2 = PDstandardNth1(beta2, i, j, k);
    PDstandardNth2beta2 = PDstandardNth2(beta2, i, j, k);
    PDstandardNth3beta2 = PDstandardNth3(beta2, i, j, k);
    PDstandardNth1beta3 = PDstandardNth1(beta3, i, j, k);
    PDstandardNth2beta3 = PDstandardNth2(beta3, i, j, k);
    PDstandardNth3beta3 = PDstandardNth3(beta3, i, j, k);
    PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    
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
    
    Xt1L  =  Gt111*gtu11 + Gt122*gtu22 + 
        2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    Xt2L  =  Gt211*gtu11 + Gt222*gtu22 + 
        2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    Xt3L  =  Gt311*gtu11 + Gt322*gtu22 + 
        2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    AL  =  -(dtalpL*(-1 + LapseAdvectionCoeff)*INV(harmonicF)*pow(alphaL,-harmonicN));
    
    B1L  =  (dtbetaxL - ((beta1L*J11L + beta2L*J12L + beta3L*J13L)*
              PDstandardNth1beta1 + (beta1L*J21L + beta2L*J22L + beta3L*J23L)*
              PDstandardNth2beta1 + (beta1L*J31L + beta2L*J32L + beta3L*J33L)*
              PDstandardNth3beta1)*ShiftAdvectionCoeff)*INV(ShiftGammaCoeff);
    
    B2L  =  (dtbetayL - ((beta1L*J11L + beta2L*J12L + beta3L*J13L)*
              PDstandardNth1beta2 + (beta1L*J21L + beta2L*J22L + beta3L*J23L)*
              PDstandardNth2beta2 + (beta1L*J31L + beta2L*J32L + beta3L*J33L)*
              PDstandardNth3beta2)*ShiftAdvectionCoeff)*INV(ShiftGammaCoeff);
    
    B3L  =  (dtbetazL - ((beta1L*J11L + beta2L*J12L + beta3L*J13L)*
              PDstandardNth1beta3 + (beta1L*J21L + beta2L*J22L + beta3L*J23L)*
              PDstandardNth2beta3 + (beta1L*J31L + beta2L*J32L + beta3L*J33L)*
              PDstandardNth3beta3)*ShiftAdvectionCoeff)*INV(ShiftGammaCoeff);
    
    
    /* Copy local copies back to grid functions */
    A[index] = AL;
    B1[index] = B1L;
    B2[index] = B2L;
    B3[index] = B3L;
    Xt1[index] = Xt1L;
    Xt2[index] = Xt2L;
    Xt3[index] = Xt3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSN_MP_convertFromADMBaseGamma);
}

void ML_BSSN_MP_convertFromADMBaseGamma(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_convertFromADMBaseGamma_Body);
}
