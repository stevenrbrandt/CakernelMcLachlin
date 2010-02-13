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

void ML_ADM_convertFromADMBase_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  /* Declare predefined quantities */
  // CCTK_REAL p1o12dx = INITVALUE;
  // CCTK_REAL p1o12dy = INITVALUE;
  // CCTK_REAL p1o12dz = INITVALUE;
  // CCTK_REAL p1o144dxdy = INITVALUE;
  // CCTK_REAL p1o144dxdz = INITVALUE;
  // CCTK_REAL p1o144dydz = INITVALUE;
  // CCTK_REAL pm1o12dx2 = INITVALUE;
  // CCTK_REAL pm1o12dy2 = INITVALUE;
  // CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_convertFromADMBase_Body");
  }
  
  if (cctk_iteration % ML_ADM_convertFromADMBase_calc_every != ML_ADM_convertFromADMBase_calc_offset)
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
  CCTK_REAL const p1o12dx = INV(dx)/12.;
  CCTK_REAL const p1o12dy = INV(dy)/12.;
  CCTK_REAL const p1o12dz = INV(dz)/12.;
  CCTK_REAL const p1o144dxdy = (INV(dx)*INV(dy))/144.;
  CCTK_REAL const p1o144dxdz = (INV(dx)*INV(dz))/144.;
  CCTK_REAL const p1o144dydz = (INV(dy)*INV(dz))/144.;
  CCTK_REAL const pm1o12dx2 = -pow(dx,-2)/12.;
  CCTK_REAL const pm1o12dy2 = -pow(dy,-2)/12.;
  CCTK_REAL const pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_ADM_convertFromADMBase,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    // int subblock_index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int const subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    
    /* Declare local copies of grid functions */
    // CCTK_REAL alpL = INITVALUE;
    // CCTK_REAL alphaL = INITVALUE;
    // CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    // CCTK_REAL betaxL = INITVALUE;
    // CCTK_REAL betayL = INITVALUE;
    // CCTK_REAL betazL = INITVALUE;
    // CCTK_REAL g11L = INITVALUE, g12L = INITVALUE, g13L = INITVALUE, g22L = INITVALUE, g23L = INITVALUE, g33L = INITVALUE;
    // CCTK_REAL gxxL = INITVALUE;
    // CCTK_REAL gxyL = INITVALUE;
    // CCTK_REAL gxzL = INITVALUE;
    // CCTK_REAL gyyL = INITVALUE;
    // CCTK_REAL gyzL = INITVALUE;
    // CCTK_REAL gzzL = INITVALUE;
    // CCTK_REAL K11L = INITVALUE, K12L = INITVALUE, K13L = INITVALUE, K22L = INITVALUE, K23L = INITVALUE, K33L = INITVALUE;
    // CCTK_REAL kxxL = INITVALUE;
    // CCTK_REAL kxyL = INITVALUE;
    // CCTK_REAL kxzL = INITVALUE;
    // CCTK_REAL kyyL = INITVALUE;
    // CCTK_REAL kyzL = INITVALUE;
    // CCTK_REAL kzzL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    CCTK_REAL const alpL = alp[index];
    CCTK_REAL const betaxL = betax[index];
    CCTK_REAL const betayL = betay[index];
    CCTK_REAL const betazL = betaz[index];
    CCTK_REAL const gxxL = gxx[index];
    CCTK_REAL const gxyL = gxy[index];
    CCTK_REAL const gxzL = gxz[index];
    CCTK_REAL const gyyL = gyy[index];
    CCTK_REAL const gyzL = gyz[index];
    CCTK_REAL const gzzL = gzz[index];
    CCTK_REAL const kxxL = kxx[index];
    CCTK_REAL const kxyL = kxy[index];
    CCTK_REAL const kxzL = kxz[index];
    CCTK_REAL const kyyL = kyy[index];
    CCTK_REAL const kyzL = kyz[index];
    CCTK_REAL const kzzL = kzz[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL const g11L  =  gxxL;
    
    CCTK_REAL const g12L  =  gxyL;
    
    CCTK_REAL const g13L  =  gxzL;
    
    CCTK_REAL const g22L  =  gyyL;
    
    CCTK_REAL const g23L  =  gyzL;
    
    CCTK_REAL const g33L  =  gzzL;
    
    CCTK_REAL const K11L  =  kxxL;
    
    CCTK_REAL const K12L  =  kxyL;
    
    CCTK_REAL const K13L  =  kxzL;
    
    CCTK_REAL const K22L  =  kyyL;
    
    CCTK_REAL const K23L  =  kyzL;
    
    CCTK_REAL const K33L  =  kzzL;
    
    CCTK_REAL const alphaL  =  alpL;
    
    CCTK_REAL const beta1L  =  betaxL;
    
    CCTK_REAL const beta2L  =  betayL;
    
    CCTK_REAL const beta3L  =  betazL;
    
    
    /* Copy local copies back to grid functions */
    alpha[index] = alphaL;
    beta1[index] = beta1L;
    beta2[index] = beta2L;
    beta3[index] = beta3L;
    g11[index] = g11L;
    g12[index] = g12L;
    g13[index] = g13L;
    g22[index] = g22L;
    g23[index] = g23L;
    g33[index] = g33L;
    K11[index] = K11L;
    K12[index] = K12L;
    K13[index] = K13L;
    K22[index] = K22L;
    K23[index] = K23L;
    K33[index] = K33L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_ADM_convertFromADMBase);
}

void ML_ADM_convertFromADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_ADM_convertFromADMBase_Body);
}
