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
#include "Vectors.hh"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

static void ML_BSSN_O2_convertToADMBaseDtLapseShift_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O2_convertToADMBaseDtLapseShift_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O2_convertToADMBaseDtLapseShift_calc_every != ML_BSSN_O2_convertToADMBaseDtLapseShift_calc_offset)
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
  CCTK_REAL const p1o16dx = INV(dx)/16.;
  CCTK_REAL const p1o16dy = INV(dy)/16.;
  CCTK_REAL const p1o16dz = INV(dz)/16.;
  CCTK_REAL const p1o2dx = khalf*INV(dx);
  CCTK_REAL const p1o2dy = khalf*INV(dy);
  CCTK_REAL const p1o2dz = khalf*INV(dz);
  CCTK_REAL const p1o4dx = INV(dx)/4.;
  CCTK_REAL const p1o4dxdy = (INV(dx)*INV(dy))/4.;
  CCTK_REAL const p1o4dxdz = (INV(dx)*INV(dz))/4.;
  CCTK_REAL const p1o4dy = INV(dy)/4.;
  CCTK_REAL const p1o4dydz = (INV(dy)*INV(dz))/4.;
  CCTK_REAL const p1o4dz = INV(dz)/4.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1odx2 = pow(dx,-2);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1ody2 = pow(dy,-2);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const p1odz2 = pow(dz,-2);
  CCTK_REAL const pm1o2dx = -(khalf*INV(dx));
  CCTK_REAL const pm1o2dy = -(khalf*INV(dy));
  CCTK_REAL const pm1o2dz = -(khalf*INV(dz));
  CCTK_REAL const pm1o4dx = -INV(dx)/4.;
  CCTK_REAL const pm1o4dy = -INV(dy)/4.;
  CCTK_REAL const pm1o4dz = -INV(dz)/4.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_O2_convertToADMBaseDtLapseShift,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL_VEC PDupwindNthAnti1alpha = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm1alpha = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti2alpha = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm2alpha = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti3alpha = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm3alpha = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti1beta1 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm1beta1 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti2beta1 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm2beta1 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti3beta1 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm3beta1 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti1beta2 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm1beta2 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti2beta2 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm2beta2 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti3beta2 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm3beta2 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti1beta3 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm1beta3 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti2beta3 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm2beta3 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti3beta3 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm3beta3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL_VEC  AL = vec_load(A[index]);
    CCTK_REAL_VEC  alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC  B1L = vec_load(B1[index]);
    CCTK_REAL_VEC  B2L = vec_load(B2[index]);
    CCTK_REAL_VEC  B3L = vec_load(B3[index]);
    CCTK_REAL_VEC  beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC  beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC  beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC  rL = vec_load(r[index]);
    CCTK_REAL_VEC  trKL = vec_load(trK[index]);
    CCTK_REAL_VEC  Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC  Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC  Xt3L = vec_load(Xt3[index]);
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC const PDupwindNthAnti1alpha = PDupwindNthAnti1(alpha, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm1alpha = PDupwindNthSymm1(alpha, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti2alpha = PDupwindNthAnti2(alpha, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm2alpha = PDupwindNthSymm2(alpha, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti3alpha = PDupwindNthAnti3(alpha, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm3alpha = PDupwindNthSymm3(alpha, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti1beta1 = PDupwindNthAnti1(beta1, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm1beta1 = PDupwindNthSymm1(beta1, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti2beta1 = PDupwindNthAnti2(beta1, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm2beta1 = PDupwindNthSymm2(beta1, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti3beta1 = PDupwindNthAnti3(beta1, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm3beta1 = PDupwindNthSymm3(beta1, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti1beta2 = PDupwindNthAnti1(beta2, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm1beta2 = PDupwindNthSymm1(beta2, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti2beta2 = PDupwindNthAnti2(beta2, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm2beta2 = PDupwindNthSymm2(beta2, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti3beta2 = PDupwindNthAnti3(beta2, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm3beta2 = PDupwindNthSymm3(beta2, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti1beta3 = PDupwindNthAnti1(beta3, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm1beta3 = PDupwindNthSymm1(beta3, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti2beta3 = PDupwindNthAnti2(beta3, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm2beta3 = PDupwindNthSymm2(beta3, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti3beta3 = PDupwindNthAnti3(beta3, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm3beta3 = PDupwindNthSymm3(beta3, i, j, k);
    
    /* Calculate temporaries and grid functions */
    int dir1 = Sign(beta1L);
    
    int dir2 = Sign(beta2L);
    
    int dir3 = Sign(beta3L);
    
    CCTK_REAL_VEC eta = fmin(1,SpatialBetaDriverRadius*INV(rL));
    
    CCTK_REAL_VEC theta = fmin(1,exp(1 - 
      rL*INV(SpatialShiftGammaCoeffRadius)));
    
    CCTK_REAL_VEC dtalpL = 
      LapseAdvectionCoeff*(beta1L*PDupwindNthAnti1alpha + 
      beta2L*PDupwindNthAnti2alpha + beta3L*PDupwindNthAnti3alpha + 
      PDupwindNthSymm1alpha*Abs(beta1L) + PDupwindNthSymm2alpha*Abs(beta2L) + 
      PDupwindNthSymm3alpha*Abs(beta3L)) - harmonicF*(LapseACoeff*(AL - trKL) 
      + trKL)*pow(alphaL,harmonicN);
    
    CCTK_REAL_VEC dtbetaxL = 
      ShiftGammaCoeff*theta*(beta1L*BetaDriver*eta*(-1 + ShiftBCoeff) + 
      ShiftBCoeff*(B1L - Xt1L) + Xt1L) + 
      ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti1beta1 + 
      beta2L*PDupwindNthAnti2beta1 + beta3L*PDupwindNthAnti3beta1 + 
      PDupwindNthSymm1beta1*Abs(beta1L) + PDupwindNthSymm2beta1*Abs(beta2L) + 
      PDupwindNthSymm3beta1*Abs(beta3L));
    
    CCTK_REAL_VEC dtbetayL = 
      ShiftGammaCoeff*theta*(beta2L*BetaDriver*eta*(-1 + ShiftBCoeff) + 
      ShiftBCoeff*(B2L - Xt2L) + Xt2L) + 
      ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti1beta2 + 
      beta2L*PDupwindNthAnti2beta2 + beta3L*PDupwindNthAnti3beta2 + 
      PDupwindNthSymm1beta2*Abs(beta1L) + PDupwindNthSymm2beta2*Abs(beta2L) + 
      PDupwindNthSymm3beta2*Abs(beta3L));
    
    CCTK_REAL_VEC dtbetazL = 
      ShiftGammaCoeff*theta*(beta3L*BetaDriver*eta*(-1 + ShiftBCoeff) + 
      ShiftBCoeff*(B3L - Xt3L) + Xt3L) + 
      ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti1beta3 + 
      beta2L*PDupwindNthAnti2beta3 + beta3L*PDupwindNthAnti3beta3 + 
      PDupwindNthSymm1beta3*Abs(beta1L) + PDupwindNthSymm2beta3*Abs(beta2L) + 
      PDupwindNthSymm3beta3*Abs(beta3L));
    
    
    /* Copy local copies back to grid functions */
    vec_store_nta(dtalp[index],dtalpL);
    vec_store_nta(dtbetax[index],dtbetaxL);
    vec_store_nta(dtbetay[index],dtbetayL);
    vec_store_nta(dtbetaz[index],dtbetazL);
  i += CCTK_REAL_VEC_SIZE-1;
  }
  LC_ENDLOOP3 (ML_BSSN_O2_convertToADMBaseDtLapseShift);
}

extern "C" void ML_BSSN_O2_convertToADMBaseDtLapseShift(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_O2_convertToADMBaseDtLapseShift_Body);
}
