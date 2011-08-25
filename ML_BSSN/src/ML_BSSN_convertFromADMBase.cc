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

static void ML_BSSN_convertFromADMBase_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_convertFromADMBase_Body");
  }
  
  if (cctk_iteration % ML_BSSN_convertFromADMBase_calc_every != ML_BSSN_convertFromADMBase_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::curv","ADMBase::lapse","ADMBase::metric","ADMBase::shift","ML_BSSN::ML_curv","ML_BSSN::ML_lapse","ML_BSSN::ML_log_confac","ML_BSSN::ML_metric","ML_BSSN::ML_shift","ML_BSSN::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_convertFromADMBase", 10, groups);
  
  
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
  CCTK_REAL_VEC const p1o12dx = kmul(INV(dx),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o12dy = kmul(INV(dy),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o12dz = kmul(INV(dz),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o144dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o144dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o144dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o24dx = kmul(INV(dx),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o24dy = kmul(INV(dy),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o24dz = kmul(INV(dz),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o64dx = kmul(INV(dx),ToReal(0.015625));
  CCTK_REAL_VEC const p1o64dy = kmul(INV(dy),ToReal(0.015625));
  CCTK_REAL_VEC const p1o64dz = kmul(INV(dz),ToReal(0.015625));
  CCTK_REAL_VEC const p1odx = INV(dx);
  CCTK_REAL_VEC const p1ody = INV(dy);
  CCTK_REAL_VEC const p1odz = INV(dz);
  CCTK_REAL_VEC const pm1o12dx2 = kmul(INV(SQR(dx)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dy2 = kmul(INV(SQR(dy)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dz2 = kmul(INV(SQR(dz)),ToReal(-0.0833333333333333333333333333333));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC (ML_BSSN_convertFromADMBase,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alpL = vec_load(alp[index]);
    CCTK_REAL_VEC betaxL = vec_load(betax[index]);
    CCTK_REAL_VEC betayL = vec_load(betay[index]);
    CCTK_REAL_VEC betazL = vec_load(betaz[index]);
    CCTK_REAL_VEC gxxL = vec_load(gxx[index]);
    CCTK_REAL_VEC gxyL = vec_load(gxy[index]);
    CCTK_REAL_VEC gxzL = vec_load(gxz[index]);
    CCTK_REAL_VEC gyyL = vec_load(gyy[index]);
    CCTK_REAL_VEC gyzL = vec_load(gyz[index]);
    CCTK_REAL_VEC gzzL = vec_load(gzz[index]);
    CCTK_REAL_VEC kxxL = vec_load(kxx[index]);
    CCTK_REAL_VEC kxyL = vec_load(kxy[index]);
    CCTK_REAL_VEC kxzL = vec_load(kxz[index]);
    CCTK_REAL_VEC kyyL = vec_load(kyy[index]);
    CCTK_REAL_VEC kyzL = vec_load(kyz[index]);
    CCTK_REAL_VEC kzzL = vec_load(kzz[index]);
    CCTK_REAL_VEC phiL = vec_load(phi[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC g11 = gxxL;
    
    CCTK_REAL_VEC g12 = gxyL;
    
    CCTK_REAL_VEC g13 = gxzL;
    
    CCTK_REAL_VEC g22 = gyyL;
    
    CCTK_REAL_VEC g23 = gyzL;
    
    CCTK_REAL_VEC g33 = gzzL;
    
    CCTK_REAL_VEC detg = 
      knmsub(g22,SQR(g13),knmsub(g11,SQR(g23),kmadd(g33,kmsub(g11,g22,SQR(g12)),kmul(g12,kmul(g13,kmul(g23,ToReal(2)))))));
    
    CCTK_REAL_VEC gu11 = kmul(INV(detg),kmsub(g22,g33,SQR(g23)));
    
    CCTK_REAL_VEC gu12 = kmul(INV(detg),kmsub(g13,g23,kmul(g12,g33)));
    
    CCTK_REAL_VEC gu13 = kmul(INV(detg),kmsub(g12,g23,kmul(g13,g22)));
    
    CCTK_REAL_VEC gu22 = kmul(INV(detg),kmsub(g11,g33,SQR(g13)));
    
    CCTK_REAL_VEC gu23 = kmul(INV(detg),kmsub(g12,g13,kmul(g11,g23)));
    
    CCTK_REAL_VEC gu33 = kmul(INV(detg),kmsub(g11,g22,SQR(g12)));
    
    CCTK_REAL_VEC em4phi;
    
    if (conformalMethod)
    {
      phiL = kpow(detg,-0.166666666666666666666666666667);
      
      em4phi = SQR(phiL);
    }
    else
    {
      phiL = kmul(klog(detg),ToReal(0.0833333333333333333333333333333));
      
      em4phi = kexp(kmul(phiL,ToReal(-4)));
    }
    
    CCTK_REAL_VEC gt11L = kmul(em4phi,g11);
    
    CCTK_REAL_VEC gt12L = kmul(em4phi,g12);
    
    CCTK_REAL_VEC gt13L = kmul(em4phi,g13);
    
    CCTK_REAL_VEC gt22L = kmul(em4phi,g22);
    
    CCTK_REAL_VEC gt23L = kmul(em4phi,g23);
    
    CCTK_REAL_VEC gt33L = kmul(em4phi,g33);
    
    trKL = 
      kmadd(gu11,kxxL,kmadd(gu22,kyyL,kmadd(gu33,kzzL,kmul(kmadd(gu12,kxyL,kmadd(gu13,kxzL,kmul(gu23,kyzL))),ToReal(2)))));
    
    CCTK_REAL_VEC At11L = 
      kmul(em4phi,kmadd(g11,kmul(trKL,ToReal(-0.333333333333333333333333333333)),kxxL));
    
    CCTK_REAL_VEC At12L = 
      kmul(em4phi,kmadd(g12,kmul(trKL,ToReal(-0.333333333333333333333333333333)),kxyL));
    
    CCTK_REAL_VEC At13L = 
      kmul(em4phi,kmadd(g13,kmul(trKL,ToReal(-0.333333333333333333333333333333)),kxzL));
    
    CCTK_REAL_VEC At22L = 
      kmul(em4phi,kmadd(g22,kmul(trKL,ToReal(-0.333333333333333333333333333333)),kyyL));
    
    CCTK_REAL_VEC At23L = 
      kmul(em4phi,kmadd(g23,kmul(trKL,ToReal(-0.333333333333333333333333333333)),kyzL));
    
    CCTK_REAL_VEC At33L = 
      kmul(em4phi,kmadd(g33,kmul(trKL,ToReal(-0.333333333333333333333333333333)),kzzL));
    
    CCTK_REAL_VEC alphaL = alpL;
    
    CCTK_REAL_VEC beta1L = betaxL;
    
    CCTK_REAL_VEC beta2L = betayL;
    
    CCTK_REAL_VEC beta3L = betazL;
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(alpha[index],alphaL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At11[index],At11L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At12[index],At12L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At13[index],At13L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At22[index],At22L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At23[index],At23L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At33[index],At33L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta1[index],beta1L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta2[index],beta2L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta3[index],beta3L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt11[index],gt11L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt12[index],gt12L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt13[index],gt13L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt22[index],gt22L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt23[index],gt23L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt33[index],gt33L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(phi[index],phiL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(trK[index],trKL,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(alpha[index],alphaL,elt_count);
      vec_store_nta_partial_hi(At11[index],At11L,elt_count);
      vec_store_nta_partial_hi(At12[index],At12L,elt_count);
      vec_store_nta_partial_hi(At13[index],At13L,elt_count);
      vec_store_nta_partial_hi(At22[index],At22L,elt_count);
      vec_store_nta_partial_hi(At23[index],At23L,elt_count);
      vec_store_nta_partial_hi(At33[index],At33L,elt_count);
      vec_store_nta_partial_hi(beta1[index],beta1L,elt_count);
      vec_store_nta_partial_hi(beta2[index],beta2L,elt_count);
      vec_store_nta_partial_hi(beta3[index],beta3L,elt_count);
      vec_store_nta_partial_hi(gt11[index],gt11L,elt_count);
      vec_store_nta_partial_hi(gt12[index],gt12L,elt_count);
      vec_store_nta_partial_hi(gt13[index],gt13L,elt_count);
      vec_store_nta_partial_hi(gt22[index],gt22L,elt_count);
      vec_store_nta_partial_hi(gt23[index],gt23L,elt_count);
      vec_store_nta_partial_hi(gt33[index],gt33L,elt_count);
      vec_store_nta_partial_hi(phi[index],phiL,elt_count);
      vec_store_nta_partial_hi(trK[index],trKL,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(alpha[index],alphaL,elt_count);
      vec_store_nta_partial_lo(At11[index],At11L,elt_count);
      vec_store_nta_partial_lo(At12[index],At12L,elt_count);
      vec_store_nta_partial_lo(At13[index],At13L,elt_count);
      vec_store_nta_partial_lo(At22[index],At22L,elt_count);
      vec_store_nta_partial_lo(At23[index],At23L,elt_count);
      vec_store_nta_partial_lo(At33[index],At33L,elt_count);
      vec_store_nta_partial_lo(beta1[index],beta1L,elt_count);
      vec_store_nta_partial_lo(beta2[index],beta2L,elt_count);
      vec_store_nta_partial_lo(beta3[index],beta3L,elt_count);
      vec_store_nta_partial_lo(gt11[index],gt11L,elt_count);
      vec_store_nta_partial_lo(gt12[index],gt12L,elt_count);
      vec_store_nta_partial_lo(gt13[index],gt13L,elt_count);
      vec_store_nta_partial_lo(gt22[index],gt22L,elt_count);
      vec_store_nta_partial_lo(gt23[index],gt23L,elt_count);
      vec_store_nta_partial_lo(gt33[index],gt33L,elt_count);
      vec_store_nta_partial_lo(phi[index],phiL,elt_count);
      vec_store_nta_partial_lo(trK[index],trKL,elt_count);
      break;
    }
    
    /* Copy local copies back to grid functions */
    vec_store_nta(alpha[index],alphaL);
    vec_store_nta(At11[index],At11L);
    vec_store_nta(At12[index],At12L);
    vec_store_nta(At13[index],At13L);
    vec_store_nta(At22[index],At22L);
    vec_store_nta(At23[index],At23L);
    vec_store_nta(At33[index],At33L);
    vec_store_nta(beta1[index],beta1L);
    vec_store_nta(beta2[index],beta2L);
    vec_store_nta(beta3[index],beta3L);
    vec_store_nta(gt11[index],gt11L);
    vec_store_nta(gt12[index],gt12L);
    vec_store_nta(gt13[index],gt13L);
    vec_store_nta(gt22[index],gt22L);
    vec_store_nta(gt23[index],gt23L);
    vec_store_nta(gt33[index],gt33L);
    vec_store_nta(phi[index],phiL);
    vec_store_nta(trK[index],trKL);
  }
  LC_ENDLOOP3VEC (ML_BSSN_convertFromADMBase);
}

extern "C" void ML_BSSN_convertFromADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_BSSN_convertFromADMBase_Body);
}
