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

static void ML_ADM_convertToADMBase_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_convertToADMBase_Body");
  }
  
  if (cctk_iteration % ML_ADM_convertToADMBase_calc_every != ML_ADM_convertToADMBase_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::curv","ADMBase::dtlapse","ADMBase::dtshift","ADMBase::lapse","ADMBase::metric","ADMBase::shift","ML_ADM::ML_curv","ML_ADM::ML_lapse","ML_ADM::ML_metric","ML_ADM::ML_shift"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADM_convertToADMBase", 10, groups);
  
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
  CCTK_REAL const p1o12dx = 0.0833333333333333333333333333333*INV(dx);
  CCTK_REAL const p1o12dy = 0.0833333333333333333333333333333*INV(dy);
  CCTK_REAL const p1o12dz = 0.0833333333333333333333333333333*INV(dz);
  CCTK_REAL const p1o144dxdy = 0.00694444444444444444444444444444*INV(dx)*INV(dy);
  CCTK_REAL const p1o144dxdz = 0.00694444444444444444444444444444*INV(dx)*INV(dz);
  CCTK_REAL const p1o144dydz = 0.00694444444444444444444444444444*INV(dy)*INV(dz);
  CCTK_REAL const pm1o12dx2 = -0.0833333333333333333333333333333*INV(SQR(dx));
  CCTK_REAL const pm1o12dy2 = -0.0833333333333333333333333333333*INV(SQR(dy));
  CCTK_REAL const pm1o12dz2 = -0.0833333333333333333333333333333*INV(SQR(dz));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_ADM_convertToADMBase,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    CCTK_REAL alphaL = alpha[index];
    CCTK_REAL beta1L = beta1[index];
    CCTK_REAL beta2L = beta2[index];
    CCTK_REAL beta3L = beta3[index];
    CCTK_REAL g11L = g11[index];
    CCTK_REAL g12L = g12[index];
    CCTK_REAL g13L = g13[index];
    CCTK_REAL g22L = g22[index];
    CCTK_REAL g23L = g23[index];
    CCTK_REAL g33L = g33[index];
    CCTK_REAL K11L = K11[index];
    CCTK_REAL K12L = K12[index];
    CCTK_REAL K13L = K13[index];
    CCTK_REAL K22L = K22[index];
    CCTK_REAL K23L = K23[index];
    CCTK_REAL K33L = K33[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL gxxL = g11L;
    
    CCTK_REAL gxyL = g12L;
    
    CCTK_REAL gxzL = g13L;
    
    CCTK_REAL gyyL = g22L;
    
    CCTK_REAL gyzL = g23L;
    
    CCTK_REAL gzzL = g33L;
    
    CCTK_REAL kxxL = K11L;
    
    CCTK_REAL kxyL = K12L;
    
    CCTK_REAL kxzL = K13L;
    
    CCTK_REAL kyyL = K22L;
    
    CCTK_REAL kyzL = K23L;
    
    CCTK_REAL kzzL = K33L;
    
    CCTK_REAL alpL = alphaL;
    
    CCTK_REAL dtalpL = 0;
    
    CCTK_REAL betaxL = beta1L;
    
    CCTK_REAL betayL = beta2L;
    
    CCTK_REAL betazL = beta3L;
    
    CCTK_REAL dtbetaxL = 0;
    
    CCTK_REAL dtbetayL = 0;
    
    CCTK_REAL dtbetazL = 0;
    
    
    /* Copy local copies back to grid functions */
    alp[index] = alpL;
    betax[index] = betaxL;
    betay[index] = betayL;
    betaz[index] = betazL;
    dtalp[index] = dtalpL;
    dtbetax[index] = dtbetaxL;
    dtbetay[index] = dtbetayL;
    dtbetaz[index] = dtbetazL;
    gxx[index] = gxxL;
    gxy[index] = gxyL;
    gxz[index] = gxzL;
    gyy[index] = gyyL;
    gyz[index] = gyzL;
    gzz[index] = gzzL;
    kxx[index] = kxxL;
    kxy[index] = kxyL;
    kxz[index] = kxzL;
    kyy[index] = kyyL;
    kyz[index] = kyzL;
    kzz[index] = kzzL;
  }
  LC_ENDLOOP3 (ML_ADM_convertToADMBase);
}

extern "C" void ML_ADM_convertToADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_ADM_convertToADMBase_Body);
}
