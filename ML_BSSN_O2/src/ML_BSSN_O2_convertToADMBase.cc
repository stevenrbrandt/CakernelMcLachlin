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

static void ML_BSSN_O2_convertToADMBase_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O2_convertToADMBase_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O2_convertToADMBase_calc_every != ML_BSSN_O2_convertToADMBase_calc_offset)
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
  LC_LOOP3 (ML_BSSN_O2_convertToADMBase,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    CCTK_REAL_VEC  alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC  At11L = vec_load(At11[index]);
    CCTK_REAL_VEC  At12L = vec_load(At12[index]);
    CCTK_REAL_VEC  At13L = vec_load(At13[index]);
    CCTK_REAL_VEC  At22L = vec_load(At22[index]);
    CCTK_REAL_VEC  At23L = vec_load(At23[index]);
    CCTK_REAL_VEC  At33L = vec_load(At33[index]);
    CCTK_REAL_VEC  beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC  beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC  beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC  gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC  gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC  gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC  gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC  gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC  gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC  gxxL = vec_load(gxx[index]);
    CCTK_REAL_VEC  gxyL = vec_load(gxy[index]);
    CCTK_REAL_VEC  gxzL = vec_load(gxz[index]);
    CCTK_REAL_VEC  gyyL = vec_load(gyy[index]);
    CCTK_REAL_VEC  gyzL = vec_load(gyz[index]);
    CCTK_REAL_VEC  gzzL = vec_load(gzz[index]);
    CCTK_REAL_VEC  phiL = vec_load(phi[index]);
    CCTK_REAL_VEC  trKL = vec_load(trK[index]);
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC e4phi = 
      IfThen(conformalMethod,pow(phiL,-2),exp(4*phiL));
    
    gxxL = e4phi*gt11L;
    
    gxyL = e4phi*gt12L;
    
    gxzL = e4phi*gt13L;
    
    gyyL = e4phi*gt22L;
    
    gyzL = e4phi*gt23L;
    
    gzzL = e4phi*gt33L;
    
    CCTK_REAL_VEC kxxL = At11L*e4phi + gxxL*kthird*trKL;
    
    CCTK_REAL_VEC kxyL = At12L*e4phi + gxyL*kthird*trKL;
    
    CCTK_REAL_VEC kxzL = At13L*e4phi + gxzL*kthird*trKL;
    
    CCTK_REAL_VEC kyyL = At22L*e4phi + gyyL*kthird*trKL;
    
    CCTK_REAL_VEC kyzL = At23L*e4phi + gyzL*kthird*trKL;
    
    CCTK_REAL_VEC kzzL = At33L*e4phi + gzzL*kthird*trKL;
    
    CCTK_REAL_VEC alpL = alphaL;
    
    CCTK_REAL_VEC betaxL = beta1L;
    
    CCTK_REAL_VEC betayL = beta2L;
    
    CCTK_REAL_VEC betazL = beta3L;
    
    
    /* Copy local copies back to grid functions */
    vec_store_nta(alp[index],alpL);
    vec_store_nta(betax[index],betaxL);
    vec_store_nta(betay[index],betayL);
    vec_store_nta(betaz[index],betazL);
    vec_store_nta(gxx[index],gxxL);
    vec_store_nta(gxy[index],gxyL);
    vec_store_nta(gxz[index],gxzL);
    vec_store_nta(gyy[index],gyyL);
    vec_store_nta(gyz[index],gyzL);
    vec_store_nta(gzz[index],gzzL);
    vec_store_nta(kxx[index],kxxL);
    vec_store_nta(kxy[index],kxyL);
    vec_store_nta(kxz[index],kxzL);
    vec_store_nta(kyy[index],kyyL);
    vec_store_nta(kyz[index],kyzL);
    vec_store_nta(kzz[index],kzzL);
  i += CCTK_REAL_VEC_SIZE-1;
  }
  LC_ENDLOOP3 (ML_BSSN_O2_convertToADMBase);
}

extern "C" void ML_BSSN_O2_convertToADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_BSSN_O2_convertToADMBase_Body);
}
