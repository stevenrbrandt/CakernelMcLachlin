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

void ML_BSSN_O2_convertFromADMBase_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O2_convertFromADMBase_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O2_convertFromADMBase_calc_every != ML_BSSN_O2_convertFromADMBase_calc_offset)
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
  LC_LOOP3 (ML_BSSN_O2_convertFromADMBase,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    CCTK_REAL  alpL = alp[index];
    CCTK_REAL  betaxL = betax[index];
    CCTK_REAL  betayL = betay[index];
    CCTK_REAL  betazL = betaz[index];
    CCTK_REAL  gxxL = gxx[index];
    CCTK_REAL  gxyL = gxy[index];
    CCTK_REAL  gxzL = gxz[index];
    CCTK_REAL  gyyL = gyy[index];
    CCTK_REAL  gyzL = gyz[index];
    CCTK_REAL  gzzL = gzz[index];
    CCTK_REAL  kxxL = kxx[index];
    CCTK_REAL  kxyL = kxy[index];
    CCTK_REAL  kxzL = kxz[index];
    CCTK_REAL  kyyL = kyy[index];
    CCTK_REAL  kyzL = kyz[index];
    CCTK_REAL  kzzL = kzz[index];
    CCTK_REAL  phiL = phi[index];
    CCTK_REAL  trKL = trK[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL g11 = gxxL;
    
    CCTK_REAL g12 = gxyL;
    
    CCTK_REAL g13 = gxzL;
    
    CCTK_REAL g22 = gyyL;
    
    CCTK_REAL g23 = gyzL;
    
    CCTK_REAL g33 = gzzL;
    
    CCTK_REAL detg = 2*g12*g13*g23 + g33*(g11*g22 - SQR(g12)) - 
      g22*SQR(g13) - g11*SQR(g23);
    
    CCTK_REAL gu11 = INV(detg)*(g22*g33 - SQR(g23));
    
    CCTK_REAL gu12 = (g13*g23 - g12*g33)*INV(detg);
    
    CCTK_REAL gu13 = (-(g13*g22) + g12*g23)*INV(detg);
    
    CCTK_REAL gu22 = INV(detg)*(g11*g33 - SQR(g13));
    
    CCTK_REAL gu23 = (g12*g13 - g11*g23)*INV(detg);
    
    CCTK_REAL gu33 = INV(detg)*(g11*g22 - SQR(g12));
    
    phiL = 
      IfThen(conformalMethod,pow(detg,-0.16666666666666666),Log(detg)/12.);
    
    CCTK_REAL em4phi = IfThen(conformalMethod,SQR(phiL),exp(-4*phiL));
    
    CCTK_REAL gt11L = em4phi*g11;
    
    CCTK_REAL gt12L = em4phi*g12;
    
    CCTK_REAL gt13L = em4phi*g13;
    
    CCTK_REAL gt22L = em4phi*g22;
    
    CCTK_REAL gt23L = em4phi*g23;
    
    CCTK_REAL gt33L = em4phi*g33;
    
    trKL = gu11*kxxL + gu22*kyyL + 2*(gu12*kxyL + gu13*kxzL + gu23*kyzL) + 
      gu33*kzzL;
    
    CCTK_REAL At11L = em4phi*(kxxL - g11*kthird*trKL);
    
    CCTK_REAL At12L = em4phi*(kxyL - g12*kthird*trKL);
    
    CCTK_REAL At13L = em4phi*(kxzL - g13*kthird*trKL);
    
    CCTK_REAL At22L = em4phi*(kyyL - g22*kthird*trKL);
    
    CCTK_REAL At23L = em4phi*(kyzL - g23*kthird*trKL);
    
    CCTK_REAL At33L = em4phi*(kzzL - g33*kthird*trKL);
    
    CCTK_REAL alphaL = alpL;
    
    CCTK_REAL beta1L = betaxL;
    
    CCTK_REAL beta2L = betayL;
    
    CCTK_REAL beta3L = betazL;
    
    
    /* Copy local copies back to grid functions */
    alpha[index] = alphaL;
    At11[index] = At11L;
    At12[index] = At12L;
    At13[index] = At13L;
    At22[index] = At22L;
    At23[index] = At23L;
    At33[index] = At33L;
    beta1[index] = beta1L;
    beta2[index] = beta2L;
    beta3[index] = beta3L;
    gt11[index] = gt11L;
    gt12[index] = gt12L;
    gt13[index] = gt13L;
    gt22[index] = gt22L;
    gt23[index] = gt23L;
    gt33[index] = gt33L;
    phi[index] = phiL;
    trK[index] = trKL;
  }
  LC_ENDLOOP3 (ML_BSSN_O2_convertFromADMBase);
}

void ML_BSSN_O2_convertFromADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_BSSN_O2_convertFromADMBase_Body);
}
