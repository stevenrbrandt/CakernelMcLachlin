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

void ML_BSSN_O8_convertFromADMBaseGamma_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  /* Declare predefined quantities */
  // CCTK_REAL p1o5040dx2 = INITVALUE;
  // CCTK_REAL p1o5040dy2 = INITVALUE;
  // CCTK_REAL p1o5040dz2 = INITVALUE;
  // CCTK_REAL p1o705600dxdy = INITVALUE;
  // CCTK_REAL p1o705600dxdz = INITVALUE;
  // CCTK_REAL p1o705600dydz = INITVALUE;
  // CCTK_REAL p1o840dx = INITVALUE;
  // CCTK_REAL p1o840dy = INITVALUE;
  // CCTK_REAL p1o840dz = INITVALUE;
  // CCTK_REAL p1odx = INITVALUE;
  // CCTK_REAL p1ody = INITVALUE;
  // CCTK_REAL p1odz = INITVALUE;
  // CCTK_REAL pm1o840dx = INITVALUE;
  // CCTK_REAL pm1o840dy = INITVALUE;
  // CCTK_REAL pm1o840dz = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O8_convertFromADMBaseGamma_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O8_convertFromADMBaseGamma_calc_every != ML_BSSN_O8_convertFromADMBaseGamma_calc_offset)
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
  CCTK_REAL const p1o5040dx2 = pow(dx,-2)/5040.;
  CCTK_REAL const p1o5040dy2 = pow(dy,-2)/5040.;
  CCTK_REAL const p1o5040dz2 = pow(dz,-2)/5040.;
  CCTK_REAL const p1o705600dxdy = (INV(dx)*INV(dy))/705600.;
  CCTK_REAL const p1o705600dxdz = (INV(dx)*INV(dz))/705600.;
  CCTK_REAL const p1o705600dydz = (INV(dy)*INV(dz))/705600.;
  CCTK_REAL const p1o840dx = INV(dx)/840.;
  CCTK_REAL const p1o840dy = INV(dy)/840.;
  CCTK_REAL const p1o840dz = INV(dz)/840.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o840dx = -INV(dx)/840.;
  CCTK_REAL const pm1o840dy = -INV(dy)/840.;
  CCTK_REAL const pm1o840dz = -INV(dz)/840.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_O8_convertFromADMBaseGamma,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    // int subblock_index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int const subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    // CCTK_REAL detgt = INITVALUE;
    // CCTK_REAL dir1 = INITVALUE, dir2 = INITVALUE, dir3 = INITVALUE;
    // CCTK_REAL Gt111 = INITVALUE, Gt112 = INITVALUE, Gt113 = INITVALUE, Gt122 = INITVALUE, Gt123 = INITVALUE, Gt133 = INITVALUE;
    // CCTK_REAL Gt211 = INITVALUE, Gt212 = INITVALUE, Gt213 = INITVALUE, Gt222 = INITVALUE, Gt223 = INITVALUE, Gt233 = INITVALUE;
    // CCTK_REAL Gt311 = INITVALUE, Gt312 = INITVALUE, Gt313 = INITVALUE, Gt322 = INITVALUE, Gt323 = INITVALUE, Gt333 = INITVALUE;
    // CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    // CCTK_REAL AL = INITVALUE;
    // CCTK_REAL alphaL = INITVALUE;
    // CCTK_REAL B1L = INITVALUE, B2L = INITVALUE, B3L = INITVALUE;
    // CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    // CCTK_REAL betaxL = INITVALUE;
    // CCTK_REAL betayL = INITVALUE;
    // CCTK_REAL betazL = INITVALUE;
    // CCTK_REAL dtalpL = INITVALUE;
    // CCTK_REAL dtbetaxL = INITVALUE;
    // CCTK_REAL dtbetayL = INITVALUE;
    // CCTK_REAL dtbetazL = INITVALUE;
    // CCTK_REAL gt11L = INITVALUE, gt12L = INITVALUE, gt13L = INITVALUE, gt22L = INITVALUE, gt23L = INITVALUE, gt33L = INITVALUE;
    // CCTK_REAL Xt1L = INITVALUE, Xt2L = INITVALUE, Xt3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
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
    
    /* Assign local copies of grid functions */
    CCTK_REAL const alphaL = alpha[index];
    CCTK_REAL const beta1L = beta1[index];
    CCTK_REAL const beta2L = beta2[index];
    CCTK_REAL const beta3L = beta3[index];
    CCTK_REAL const betaxL = betax[index];
    CCTK_REAL const betayL = betay[index];
    CCTK_REAL const betazL = betaz[index];
    CCTK_REAL const dtalpL = dtalp[index];
    CCTK_REAL const dtbetaxL = dtbetax[index];
    CCTK_REAL const dtbetayL = dtbetay[index];
    CCTK_REAL const dtbetazL = dtbetaz[index];
    CCTK_REAL const gt11L = gt11[index];
    CCTK_REAL const gt12L = gt12[index];
    CCTK_REAL const gt13L = gt13[index];
    CCTK_REAL const gt22L = gt22[index];
    CCTK_REAL const gt23L = gt23[index];
    CCTK_REAL const gt33L = gt33[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
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
    
    CCTK_REAL const Xt1L  =  Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    CCTK_REAL const Xt2L  =  Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    CCTK_REAL const Xt3L  =  Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    CCTK_REAL const AL  =  -(dtalpL*(-1 + LapseAdvectionCoeff)*INV(harmonicF)*pow(alphaL,-harmonicN));
    
    CCTK_REAL const B1L  =  6*IfThen(ShiftGammaCoeff != 0,dtbetaxL*INV(ShiftGammaCoeff),0) + 
        IfThen(ShiftGammaCoeff != 0,(dtbetaxL - PDupwindNth1(betax, i, j, k)*beta1L*ShiftAdvectionCoeff)*
          INV(ShiftGammaCoeff),0) + IfThen(ShiftGammaCoeff != 0,
         (dtbetaxL - PDupwindNth2(betax, i, j, k)*beta2L*ShiftAdvectionCoeff)*INV(ShiftGammaCoeff),0) + 
        IfThen(ShiftGammaCoeff != 0,(dtbetaxL - PDupwindNth3(betax, i, j, k)*beta3L*ShiftAdvectionCoeff)*
          INV(ShiftGammaCoeff),0);
    
    CCTK_REAL const B2L  =  6*IfThen(ShiftGammaCoeff != 0,dtbetayL*INV(ShiftGammaCoeff),0) + 
        IfThen(ShiftGammaCoeff != 0,(dtbetayL - PDupwindNth1(betay, i, j, k)*beta1L*ShiftAdvectionCoeff)*
          INV(ShiftGammaCoeff),0) + IfThen(ShiftGammaCoeff != 0,
         (dtbetayL - PDupwindNth2(betay, i, j, k)*beta2L*ShiftAdvectionCoeff)*INV(ShiftGammaCoeff),0) + 
        IfThen(ShiftGammaCoeff != 0,(dtbetayL - PDupwindNth3(betay, i, j, k)*beta3L*ShiftAdvectionCoeff)*
          INV(ShiftGammaCoeff),0);
    
    CCTK_REAL const B3L  =  6*IfThen(ShiftGammaCoeff != 0,dtbetazL*INV(ShiftGammaCoeff),0) + 
        IfThen(ShiftGammaCoeff != 0,(dtbetazL - PDupwindNth1(betaz, i, j, k)*beta1L*ShiftAdvectionCoeff)*
          INV(ShiftGammaCoeff),0) + IfThen(ShiftGammaCoeff != 0,
         (dtbetazL - PDupwindNth2(betaz, i, j, k)*beta2L*ShiftAdvectionCoeff)*INV(ShiftGammaCoeff),0) + 
        IfThen(ShiftGammaCoeff != 0,(dtbetazL - PDupwindNth3(betaz, i, j, k)*beta3L*ShiftAdvectionCoeff)*
          INV(ShiftGammaCoeff),0);
    
    
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
  LC_ENDLOOP3 (ML_BSSN_O8_convertFromADMBaseGamma);
}

void ML_BSSN_O8_convertFromADMBaseGamma(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_O8_convertFromADMBaseGamma_Body);
}
