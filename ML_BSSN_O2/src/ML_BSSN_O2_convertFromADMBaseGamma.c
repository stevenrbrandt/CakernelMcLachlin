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

void ML_BSSN_O2_convertFromADMBaseGamma_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O2_convertFromADMBaseGamma_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O2_convertFromADMBaseGamma_calc_every != ML_BSSN_O2_convertFromADMBaseGamma_calc_offset)
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
  LC_LOOP3 (ML_BSSN_O2_convertFromADMBaseGamma,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
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
    CCTK_REAL  alphaL = alpha[index];
    CCTK_REAL  beta1L = beta1[index];
    CCTK_REAL  beta2L = beta2[index];
    CCTK_REAL  beta3L = beta3[index];
    CCTK_REAL  betaxL = betax[index];
    CCTK_REAL  betayL = betay[index];
    CCTK_REAL  betazL = betaz[index];
    CCTK_REAL  dtalpL = dtalp[index];
    CCTK_REAL  dtbetaxL = dtbetax[index];
    CCTK_REAL  dtbetayL = dtbetay[index];
    CCTK_REAL  dtbetazL = dtbetaz[index];
    CCTK_REAL  gt11L = gt11[index];
    CCTK_REAL  gt12L = gt12[index];
    CCTK_REAL  gt13L = gt13[index];
    CCTK_REAL  gt22L = gt22[index];
    CCTK_REAL  gt23L = gt23[index];
    CCTK_REAL  gt33L = gt33[index];
    CCTK_REAL  rL = r[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
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
    
    CCTK_REAL Xt1L = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + 
      Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    CCTK_REAL Xt2L = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + 
      Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    CCTK_REAL Xt3L = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + 
      Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    CCTK_REAL AL = -(dtalpL*(-1 + 
      LapseAdvectionCoeff)*INV(harmonicF)*pow(alphaL,-harmonicN));
    
    CCTK_REAL theta = IfThen(rL > SpatialShiftGammaCoeffRadius,exp(1 - 
      rL*INV(SpatialShiftGammaCoeffRadius)),1);
    
    CCTK_REAL B1L = 6*IfThen(ShiftGammaCoeff*theta != 
      0,dtbetaxL*INV(ShiftGammaCoeff)*INV(theta),0) + 
      IfThen(ShiftGammaCoeff*theta != 0,(dtbetaxL - PDupwindNth1(betax, i, 
      j, 
      k)*beta1L*ShiftAdvectionCoeff*theta)*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftGammaCoeff*theta != 0,(dtbetaxL - PDupwindNth2(betax, i, 
      j, 
      k)*beta2L*ShiftAdvectionCoeff*theta)*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftGammaCoeff*theta != 0,(dtbetaxL - PDupwindNth3(betax, i, 
      j, 
      k)*beta3L*ShiftAdvectionCoeff*theta)*INV(ShiftGammaCoeff)*INV(theta),0);
    
    CCTK_REAL B2L = 6*IfThen(ShiftGammaCoeff*theta != 
      0,dtbetayL*INV(ShiftGammaCoeff)*INV(theta),0) + 
      IfThen(ShiftGammaCoeff*theta != 0,(dtbetayL - PDupwindNth1(betay, i, 
      j, 
      k)*beta1L*ShiftAdvectionCoeff*theta)*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftGammaCoeff*theta != 0,(dtbetayL - PDupwindNth2(betay, i, 
      j, 
      k)*beta2L*ShiftAdvectionCoeff*theta)*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftGammaCoeff*theta != 0,(dtbetayL - PDupwindNth3(betay, i, 
      j, 
      k)*beta3L*ShiftAdvectionCoeff*theta)*INV(ShiftGammaCoeff)*INV(theta),0);
    
    CCTK_REAL B3L = 6*IfThen(ShiftGammaCoeff*theta != 
      0,dtbetazL*INV(ShiftGammaCoeff)*INV(theta),0) + 
      IfThen(ShiftGammaCoeff*theta != 0,(dtbetazL - PDupwindNth1(betaz, i, 
      j, 
      k)*beta1L*ShiftAdvectionCoeff*theta)*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftGammaCoeff*theta != 0,(dtbetazL - PDupwindNth2(betaz, i, 
      j, 
      k)*beta2L*ShiftAdvectionCoeff*theta)*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftGammaCoeff*theta != 0,(dtbetazL - PDupwindNth3(betaz, i, 
      j, 
      k)*beta3L*ShiftAdvectionCoeff*theta)*INV(ShiftGammaCoeff)*INV(theta),0);
    
    
    /* Copy local copies back to grid functions */
    A[index] = AL;
    B1[index] = B1L;
    B2[index] = B2L;
    B3[index] = B3L;
    Xt1[index] = Xt1L;
    Xt2[index] = Xt2L;
    Xt3[index] = Xt3L;
  }
  LC_ENDLOOP3 (ML_BSSN_O2_convertFromADMBaseGamma);
}

void ML_BSSN_O2_convertFromADMBaseGamma(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_O2_convertFromADMBaseGamma_Body);
}
