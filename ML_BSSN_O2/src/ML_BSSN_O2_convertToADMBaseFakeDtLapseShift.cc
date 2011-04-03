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
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

static void ML_BSSN_O2_convertToADMBaseFakeDtLapseShift_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O2_convertToADMBaseFakeDtLapseShift_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O2_convertToADMBaseFakeDtLapseShift_calc_every != ML_BSSN_O2_convertToADMBaseFakeDtLapseShift_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::dtlapse","ADMBase::dtshift","grid::coordinates","Grid::coordinates","ML_BSSN_O2::ML_dtlapse","ML_BSSN_O2::ML_dtshift","ML_BSSN_O2::ML_Gamma","ML_BSSN_O2::ML_lapse","ML_BSSN_O2::ML_shift","ML_BSSN_O2::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_O2_convertToADMBaseFakeDtLapseShift", 10, groups);
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  CCTK_REAL const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL const dz = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL const dxi = INV(dx);
  CCTK_REAL const dyi = INV(dy);
  CCTK_REAL const dzi = INV(dz);
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  CCTK_REAL const p1o16dx = 0.0625*INV(dx);
  CCTK_REAL const p1o16dy = 0.0625*INV(dy);
  CCTK_REAL const p1o16dz = 0.0625*INV(dz);
  CCTK_REAL const p1o2dx = 0.5*INV(dx);
  CCTK_REAL const p1o2dy = 0.5*INV(dy);
  CCTK_REAL const p1o2dz = 0.5*INV(dz);
  CCTK_REAL const p1o4dx = 0.25*INV(dx);
  CCTK_REAL const p1o4dxdy = 0.25*INV(dx)*INV(dy);
  CCTK_REAL const p1o4dxdz = 0.25*INV(dx)*INV(dz);
  CCTK_REAL const p1o4dy = 0.25*INV(dy);
  CCTK_REAL const p1o4dydz = 0.25*INV(dy)*INV(dz);
  CCTK_REAL const p1o4dz = 0.25*INV(dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1odx2 = INV(SQR(dx));
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1ody2 = INV(SQR(dy));
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const p1odz2 = INV(SQR(dz));
  CCTK_REAL const pm1o2dx = -0.5*INV(dx);
  CCTK_REAL const pm1o2dy = -0.5*INV(dy);
  CCTK_REAL const pm1o2dz = -0.5*INV(dz);
  CCTK_REAL const pm1o4dx = -0.25*INV(dx);
  CCTK_REAL const pm1o4dy = -0.25*INV(dy);
  CCTK_REAL const pm1o4dz = -0.25*INV(dz);
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_O2_convertToADMBaseFakeDtLapseShift,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    CCTK_REAL AL = A[index];
    CCTK_REAL alphaL = alpha[index];
    CCTK_REAL B1L = B1[index];
    CCTK_REAL B2L = B2[index];
    CCTK_REAL B3L = B3[index];
    CCTK_REAL beta1L = beta1[index];
    CCTK_REAL beta2L = beta2[index];
    CCTK_REAL beta3L = beta3[index];
    CCTK_REAL rL = r[index];
    CCTK_REAL trKL = trK[index];
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt3L = Xt3[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL eta = fmin(1,INV(rL)*ToReal(SpatialBetaDriverRadius));
    
    CCTK_REAL theta = fmin(1,exp(1 - 
      rL*INV(ToReal(SpatialShiftGammaCoeffRadius))));
    
    CCTK_REAL dtalpL = 
      -(pow(alphaL,ToReal(harmonicN))*ToReal(harmonicF)*(trKL + (AL - 
      trKL)*ToReal(LapseACoeff)));
    
    CCTK_REAL dtbetaxL = theta*(Xt1L + beta1L*eta*ToReal(BetaDriver)*(-1 + 
      ToReal(ShiftBCoeff)) + (B1L - 
      Xt1L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL dtbetayL = theta*(Xt2L + beta2L*eta*ToReal(BetaDriver)*(-1 + 
      ToReal(ShiftBCoeff)) + (B2L - 
      Xt2L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL dtbetazL = theta*(Xt3L + beta3L*eta*ToReal(BetaDriver)*(-1 + 
      ToReal(ShiftBCoeff)) + (B3L - 
      Xt3L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    
    /* Copy local copies back to grid functions */
    dtalp[index] = dtalpL;
    dtbetax[index] = dtbetaxL;
    dtbetay[index] = dtbetayL;
    dtbetaz[index] = dtbetazL;
  }
  LC_ENDLOOP3 (ML_BSSN_O2_convertToADMBaseFakeDtLapseShift);
}

extern "C" void ML_BSSN_O2_convertToADMBaseFakeDtLapseShift(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_BSSN_O2_convertToADMBaseFakeDtLapseShift_Body);
}
