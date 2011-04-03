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

static void ML_BSSN_MP_O8_convertFromADMBase_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_O8_convertFromADMBase_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_O8_convertFromADMBase_calc_every != ML_BSSN_MP_O8_convertFromADMBase_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::curv","ADMBase::lapse","ADMBase::metric","ADMBase::shift","ML_BSSN_MP_O8::ML_curv","ML_BSSN_MP_O8::ML_lapse","ML_BSSN_MP_O8::ML_log_confac","ML_BSSN_MP_O8::ML_metric","ML_BSSN_MP_O8::ML_shift","ML_BSSN_MP_O8::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_O8_convertFromADMBase", 10, groups);
  
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
  CCTK_REAL const p1o1024dx = 0.0009765625*INV(dx);
  CCTK_REAL const p1o1024dy = 0.0009765625*INV(dy);
  CCTK_REAL const p1o1024dz = 0.0009765625*INV(dz);
  CCTK_REAL const p1o1680dx = 0.000595238095238095238095238095238*INV(dx);
  CCTK_REAL const p1o1680dy = 0.000595238095238095238095238095238*INV(dy);
  CCTK_REAL const p1o1680dz = 0.000595238095238095238095238095238*INV(dz);
  CCTK_REAL const p1o5040dx2 = 0.000198412698412698412698412698413*INV(SQR(dx));
  CCTK_REAL const p1o5040dy2 = 0.000198412698412698412698412698413*INV(SQR(dy));
  CCTK_REAL const p1o5040dz2 = 0.000198412698412698412698412698413*INV(SQR(dz));
  CCTK_REAL const p1o560dx = 0.00178571428571428571428571428571*INV(dx);
  CCTK_REAL const p1o560dy = 0.00178571428571428571428571428571*INV(dy);
  CCTK_REAL const p1o560dz = 0.00178571428571428571428571428571*INV(dz);
  CCTK_REAL const p1o705600dxdy = 1.41723356009070294784580498866e-6*INV(dx)*INV(dy);
  CCTK_REAL const p1o705600dxdz = 1.41723356009070294784580498866e-6*INV(dx)*INV(dz);
  CCTK_REAL const p1o705600dydz = 1.41723356009070294784580498866e-6*INV(dy)*INV(dz);
  CCTK_REAL const p1o840dx = 0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const p1o840dy = 0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const p1o840dz = 0.00119047619047619047619047619048*INV(dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o840dx = -0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const pm1o840dy = -0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const pm1o840dz = -0.00119047619047619047619047619048*INV(dz);
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_MP_O8_convertFromADMBase,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    CCTK_REAL alpL = alp[index];
    CCTK_REAL betaxL = betax[index];
    CCTK_REAL betayL = betay[index];
    CCTK_REAL betazL = betaz[index];
    CCTK_REAL gxxL = gxx[index];
    CCTK_REAL gxyL = gxy[index];
    CCTK_REAL gxzL = gxz[index];
    CCTK_REAL gyyL = gyy[index];
    CCTK_REAL gyzL = gyz[index];
    CCTK_REAL gzzL = gzz[index];
    CCTK_REAL kxxL = kxx[index];
    CCTK_REAL kxyL = kxy[index];
    CCTK_REAL kxzL = kxz[index];
    CCTK_REAL kyyL = kyy[index];
    CCTK_REAL kyzL = kyz[index];
    CCTK_REAL kzzL = kzz[index];
    CCTK_REAL phiL = phi[index];
    CCTK_REAL trKL = trK[index];
    
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
      IfThen(ToReal(conformalMethod),pow(detg,-0.166666666666666666666666666667),0.0833333333333333333333333333333*Log(detg));
    
    CCTK_REAL em4phi = 
      IfThen(ToReal(conformalMethod),SQR(phiL),exp(-4*phiL));
    
    CCTK_REAL gt11L = em4phi*g11;
    
    CCTK_REAL gt12L = em4phi*g12;
    
    CCTK_REAL gt13L = em4phi*g13;
    
    CCTK_REAL gt22L = em4phi*g22;
    
    CCTK_REAL gt23L = em4phi*g23;
    
    CCTK_REAL gt33L = em4phi*g33;
    
    trKL = gu11*kxxL + gu22*kyyL + 2*(gu12*kxyL + gu13*kxzL + gu23*kyzL) + 
      gu33*kzzL;
    
    CCTK_REAL At11L = em4phi*(kxxL - 
      0.333333333333333333333333333333*g11*trKL);
    
    CCTK_REAL At12L = em4phi*(kxyL - 
      0.333333333333333333333333333333*g12*trKL);
    
    CCTK_REAL At13L = em4phi*(kxzL - 
      0.333333333333333333333333333333*g13*trKL);
    
    CCTK_REAL At22L = em4phi*(kyyL - 
      0.333333333333333333333333333333*g22*trKL);
    
    CCTK_REAL At23L = em4phi*(kyzL - 
      0.333333333333333333333333333333*g23*trKL);
    
    CCTK_REAL At33L = em4phi*(kzzL - 
      0.333333333333333333333333333333*g33*trKL);
    
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
  LC_ENDLOOP3 (ML_BSSN_MP_O8_convertFromADMBase);
}

extern "C" void ML_BSSN_MP_O8_convertFromADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_BSSN_MP_O8_convertFromADMBase_Body);
}
