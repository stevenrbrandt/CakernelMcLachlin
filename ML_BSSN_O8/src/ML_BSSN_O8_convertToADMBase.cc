/*  File produced by Kranc */

#define KRANC_C

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"
#include "loopcontrol.h"
#include "vectors.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))

static void ML_BSSN_O8_convertToADMBase_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O8_convertToADMBase_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O8_convertToADMBase_calc_every != ML_BSSN_O8_convertToADMBase_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::curv","ADMBase::lapse","ADMBase::metric","ADMBase::shift","ML_BSSN_O8::ML_curv","ML_BSSN_O8::ML_lapse","ML_BSSN_O8::ML_log_confac","ML_BSSN_O8::ML_metric","ML_BSSN_O8::ML_shift","ML_BSSN_O8::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_O8_convertToADMBase", 10, groups);
  
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const cdi = sizeof(CCTK_REAL) * di;
  ptrdiff_t const cdj = sizeof(CCTK_REAL) * dj;
  ptrdiff_t const cdk = sizeof(CCTK_REAL) * dk;
  CCTK_REAL_VEC const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL_VEC const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL_VEC const dz = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL_VEC const dt = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL_VEC const dxi = INV(dx);
  CCTK_REAL_VEC const dyi = INV(dy);
  CCTK_REAL_VEC const dzi = INV(dz);
  CCTK_REAL_VEC const khalf = ToReal(0.5);
  CCTK_REAL_VEC const kthird = ToReal(1.0/3.0);
  CCTK_REAL_VEC const ktwothird = ToReal(2.0/3.0);
  CCTK_REAL_VEC const kfourthird = ToReal(4.0/3.0);
  CCTK_REAL_VEC const keightthird = ToReal(8.0/3.0);
  CCTK_REAL_VEC const hdxi = kmul(ToReal(0.5), dxi);
  CCTK_REAL_VEC const hdyi = kmul(ToReal(0.5), dyi);
  CCTK_REAL_VEC const hdzi = kmul(ToReal(0.5), dzi);
  
  /* Initialize predefined quantities */
  CCTK_REAL_VEC const p1o1024dx = kmul(INV(dx),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1024dy = kmul(INV(dy),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1024dz = kmul(INV(dz),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1680dx = kmul(INV(dx),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o1680dy = kmul(INV(dy),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o1680dz = kmul(INV(dz),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o5040dx2 = kmul(INV(SQR(dx)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dy2 = kmul(INV(SQR(dy)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dz2 = kmul(INV(SQR(dz)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o560dx = kmul(INV(dx),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o560dy = kmul(INV(dy),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o560dz = kmul(INV(dz),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o705600dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o705600dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o705600dydz = kmul(INV(dy),kmul(INV(dz),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o840dx = kmul(INV(dx),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dy = kmul(INV(dy),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dz = kmul(INV(dz),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1odx = INV(dx);
  CCTK_REAL_VEC const p1ody = INV(dy);
  CCTK_REAL_VEC const p1odz = INV(dz);
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC (ML_BSSN_O8_convertToADMBase,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L = vec_load(At11[index]);
    CCTK_REAL_VEC At12L = vec_load(At12[index]);
    CCTK_REAL_VEC At13L = vec_load(At13[index]);
    CCTK_REAL_VEC At22L = vec_load(At22[index]);
    CCTK_REAL_VEC At23L = vec_load(At23[index]);
    CCTK_REAL_VEC At33L = vec_load(At33[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC gxxL = vec_load(gxx[index]);
    CCTK_REAL_VEC gxyL = vec_load(gxy[index]);
    CCTK_REAL_VEC gxzL = vec_load(gxz[index]);
    CCTK_REAL_VEC gyyL = vec_load(gyy[index]);
    CCTK_REAL_VEC gyzL = vec_load(gyz[index]);
    CCTK_REAL_VEC gzzL = vec_load(gzz[index]);
    CCTK_REAL_VEC phiL = vec_load(phi[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC e4phi = 
      IfThen(conformalMethod,INV(SQR(phiL)),kexp(kmul(phiL,ToReal(4))));
    
    gxxL = kmul(e4phi,gt11L);
    
    gxyL = kmul(e4phi,gt12L);
    
    gxzL = kmul(e4phi,gt13L);
    
    gyyL = kmul(e4phi,gt22L);
    
    gyzL = kmul(e4phi,gt23L);
    
    gzzL = kmul(e4phi,gt33L);
    
    CCTK_REAL_VEC kxxL = 
      kmadd(At11L,e4phi,kmul(gxxL,kmul(trKL,ToReal(0.333333333333333333333333333333))));
    
    CCTK_REAL_VEC kxyL = 
      kmadd(At12L,e4phi,kmul(gxyL,kmul(trKL,ToReal(0.333333333333333333333333333333))));
    
    CCTK_REAL_VEC kxzL = 
      kmadd(At13L,e4phi,kmul(gxzL,kmul(trKL,ToReal(0.333333333333333333333333333333))));
    
    CCTK_REAL_VEC kyyL = 
      kmadd(At22L,e4phi,kmul(gyyL,kmul(trKL,ToReal(0.333333333333333333333333333333))));
    
    CCTK_REAL_VEC kyzL = 
      kmadd(At23L,e4phi,kmul(gyzL,kmul(trKL,ToReal(0.333333333333333333333333333333))));
    
    CCTK_REAL_VEC kzzL = 
      kmadd(At33L,e4phi,kmul(gzzL,kmul(trKL,ToReal(0.333333333333333333333333333333))));
    
    CCTK_REAL_VEC alpL = alphaL;
    
    CCTK_REAL_VEC betaxL = beta1L;
    
    CCTK_REAL_VEC betayL = beta2L;
    
    CCTK_REAL_VEC betazL = beta3L;
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(alp[index],alpL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(betax[index],betaxL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(betay[index],betayL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(betaz[index],betazL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gxx[index],gxxL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gxy[index],gxyL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gxz[index],gxzL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gyy[index],gyyL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gyz[index],gyzL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gzz[index],gzzL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(kxx[index],kxxL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(kxy[index],kxyL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(kxz[index],kxzL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(kyy[index],kyyL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(kyz[index],kyzL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(kzz[index],kzzL,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(alp[index],alpL,elt_count);
      vec_store_nta_partial_hi(betax[index],betaxL,elt_count);
      vec_store_nta_partial_hi(betay[index],betayL,elt_count);
      vec_store_nta_partial_hi(betaz[index],betazL,elt_count);
      vec_store_nta_partial_hi(gxx[index],gxxL,elt_count);
      vec_store_nta_partial_hi(gxy[index],gxyL,elt_count);
      vec_store_nta_partial_hi(gxz[index],gxzL,elt_count);
      vec_store_nta_partial_hi(gyy[index],gyyL,elt_count);
      vec_store_nta_partial_hi(gyz[index],gyzL,elt_count);
      vec_store_nta_partial_hi(gzz[index],gzzL,elt_count);
      vec_store_nta_partial_hi(kxx[index],kxxL,elt_count);
      vec_store_nta_partial_hi(kxy[index],kxyL,elt_count);
      vec_store_nta_partial_hi(kxz[index],kxzL,elt_count);
      vec_store_nta_partial_hi(kyy[index],kyyL,elt_count);
      vec_store_nta_partial_hi(kyz[index],kyzL,elt_count);
      vec_store_nta_partial_hi(kzz[index],kzzL,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(alp[index],alpL,elt_count);
      vec_store_nta_partial_lo(betax[index],betaxL,elt_count);
      vec_store_nta_partial_lo(betay[index],betayL,elt_count);
      vec_store_nta_partial_lo(betaz[index],betazL,elt_count);
      vec_store_nta_partial_lo(gxx[index],gxxL,elt_count);
      vec_store_nta_partial_lo(gxy[index],gxyL,elt_count);
      vec_store_nta_partial_lo(gxz[index],gxzL,elt_count);
      vec_store_nta_partial_lo(gyy[index],gyyL,elt_count);
      vec_store_nta_partial_lo(gyz[index],gyzL,elt_count);
      vec_store_nta_partial_lo(gzz[index],gzzL,elt_count);
      vec_store_nta_partial_lo(kxx[index],kxxL,elt_count);
      vec_store_nta_partial_lo(kxy[index],kxyL,elt_count);
      vec_store_nta_partial_lo(kxz[index],kxzL,elt_count);
      vec_store_nta_partial_lo(kyy[index],kyyL,elt_count);
      vec_store_nta_partial_lo(kyz[index],kyzL,elt_count);
      vec_store_nta_partial_lo(kzz[index],kzzL,elt_count);
      break;
    }
    
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
  }
  LC_ENDLOOP3VEC (ML_BSSN_O8_convertToADMBase);
}

extern "C" void ML_BSSN_O8_convertToADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_BSSN_O8_convertToADMBase_Body);
}
