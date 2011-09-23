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
  
  switch(fdOrder)
  {
    case 2:
      break;
    
    case 4:
      break;
    
    case 6:
      break;
    
    case 8:
      break;
  }
  
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
  CCTK_REAL_VEC const p1o180dx2 = kmul(INV(SQR(dx)),ToReal(0.00555555555555555555555555555556));
  CCTK_REAL_VEC const p1o180dy2 = kmul(INV(SQR(dy)),ToReal(0.00555555555555555555555555555556));
  CCTK_REAL_VEC const p1o180dz2 = kmul(INV(SQR(dz)),ToReal(0.00555555555555555555555555555556));
  CCTK_REAL_VEC const p1o2dx = kmul(INV(dx),ToReal(0.5));
  CCTK_REAL_VEC const p1o2dy = kmul(INV(dy),ToReal(0.5));
  CCTK_REAL_VEC const p1o2dz = kmul(INV(dz),ToReal(0.5));
  CCTK_REAL_VEC const p1o3600dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.000277777777777777777777777777778)));
  CCTK_REAL_VEC const p1o3600dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.000277777777777777777777777777778)));
  CCTK_REAL_VEC const p1o3600dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.000277777777777777777777777777778)));
  CCTK_REAL_VEC const p1o4dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.25)));
  CCTK_REAL_VEC const p1o4dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.25)));
  CCTK_REAL_VEC const p1o4dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.25)));
  CCTK_REAL_VEC const p1o5040dx2 = kmul(INV(SQR(dx)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dy2 = kmul(INV(SQR(dy)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dz2 = kmul(INV(SQR(dz)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o60dx = kmul(INV(dx),ToReal(0.0166666666666666666666666666667));
  CCTK_REAL_VEC const p1o60dy = kmul(INV(dy),ToReal(0.0166666666666666666666666666667));
  CCTK_REAL_VEC const p1o60dz = kmul(INV(dz),ToReal(0.0166666666666666666666666666667));
  CCTK_REAL_VEC const p1o705600dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o705600dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o705600dydz = kmul(INV(dy),kmul(INV(dz),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o840dx = kmul(INV(dx),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dy = kmul(INV(dy),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dz = kmul(INV(dz),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1odx2 = INV(SQR(dx));
  CCTK_REAL_VEC const p1ody2 = INV(SQR(dy));
  CCTK_REAL_VEC const p1odz2 = INV(SQR(dz));
  CCTK_REAL_VEC const pm1o12dx2 = kmul(INV(SQR(dx)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dy2 = kmul(INV(SQR(dy)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dz2 = kmul(INV(SQR(dz)),ToReal(-0.0833333333333333333333333333333));
  
  /* Jacobian variable pointers */
  bool const use_jacobian = (!CCTK_IsFunctionAliased("MultiPatch_GetMap") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)
                       && strlen(jacobian_group) > 0;
  if (use_jacobian && strlen(jacobian_derivative_group) == 0)
  {
    CCTK_WARN (1, "GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names");
  }
  
  CCTK_REAL const *restrict jacobian_ptrs[9];
  if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_group,
                                                9, jacobian_ptrs);
  
  CCTK_REAL const *restrict const J11 = use_jacobian ? jacobian_ptrs[0] : 0;
  CCTK_REAL const *restrict const J12 = use_jacobian ? jacobian_ptrs[1] : 0;
  CCTK_REAL const *restrict const J13 = use_jacobian ? jacobian_ptrs[2] : 0;
  CCTK_REAL const *restrict const J21 = use_jacobian ? jacobian_ptrs[3] : 0;
  CCTK_REAL const *restrict const J22 = use_jacobian ? jacobian_ptrs[4] : 0;
  CCTK_REAL const *restrict const J23 = use_jacobian ? jacobian_ptrs[5] : 0;
  CCTK_REAL const *restrict const J31 = use_jacobian ? jacobian_ptrs[6] : 0;
  CCTK_REAL const *restrict const J32 = use_jacobian ? jacobian_ptrs[7] : 0;
  CCTK_REAL const *restrict const J33 = use_jacobian ? jacobian_ptrs[8] : 0;
  
  CCTK_REAL const *restrict jacobian_derivative_ptrs[18];
  if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_derivative_group,
                                                18, jacobian_derivative_ptrs);
  
  CCTK_REAL const *restrict const dJ111 = use_jacobian ? jacobian_derivative_ptrs[0] : 0;
  CCTK_REAL const *restrict const dJ112 = use_jacobian ? jacobian_derivative_ptrs[1] : 0;
  CCTK_REAL const *restrict const dJ113 = use_jacobian ? jacobian_derivative_ptrs[2] : 0;
  CCTK_REAL const *restrict const dJ122 = use_jacobian ? jacobian_derivative_ptrs[3] : 0;
  CCTK_REAL const *restrict const dJ123 = use_jacobian ? jacobian_derivative_ptrs[4] : 0;
  CCTK_REAL const *restrict const dJ133 = use_jacobian ? jacobian_derivative_ptrs[5] : 0;
  CCTK_REAL const *restrict const dJ211 = use_jacobian ? jacobian_derivative_ptrs[6] : 0;
  CCTK_REAL const *restrict const dJ212 = use_jacobian ? jacobian_derivative_ptrs[7] : 0;
  CCTK_REAL const *restrict const dJ213 = use_jacobian ? jacobian_derivative_ptrs[8] : 0;
  CCTK_REAL const *restrict const dJ222 = use_jacobian ? jacobian_derivative_ptrs[9] : 0;
  CCTK_REAL const *restrict const dJ223 = use_jacobian ? jacobian_derivative_ptrs[10] : 0;
  CCTK_REAL const *restrict const dJ233 = use_jacobian ? jacobian_derivative_ptrs[11] : 0;
  CCTK_REAL const *restrict const dJ311 = use_jacobian ? jacobian_derivative_ptrs[12] : 0;
  CCTK_REAL const *restrict const dJ312 = use_jacobian ? jacobian_derivative_ptrs[13] : 0;
  CCTK_REAL const *restrict const dJ313 = use_jacobian ? jacobian_derivative_ptrs[14] : 0;
  CCTK_REAL const *restrict const dJ322 = use_jacobian ? jacobian_derivative_ptrs[15] : 0;
  CCTK_REAL const *restrict const dJ323 = use_jacobian ? jacobian_derivative_ptrs[16] : 0;
  CCTK_REAL const *restrict const dJ333 = use_jacobian ? jacobian_derivative_ptrs[17] : 0;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC (ML_ADM_convertToADMBase,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC g11L = vec_load(g11[index]);
    CCTK_REAL_VEC g12L = vec_load(g12[index]);
    CCTK_REAL_VEC g13L = vec_load(g13[index]);
    CCTK_REAL_VEC g22L = vec_load(g22[index]);
    CCTK_REAL_VEC g23L = vec_load(g23[index]);
    CCTK_REAL_VEC g33L = vec_load(g33[index]);
    CCTK_REAL_VEC K11L = vec_load(K11[index]);
    CCTK_REAL_VEC K12L = vec_load(K12[index]);
    CCTK_REAL_VEC K13L = vec_load(K13[index]);
    CCTK_REAL_VEC K22L = vec_load(K22[index]);
    CCTK_REAL_VEC K23L = vec_load(K23[index]);
    CCTK_REAL_VEC K33L = vec_load(K33[index]);
    
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    switch(fdOrder)
    {
      case 2:
        break;
      
      case 4:
        break;
      
      case 6:
        break;
      
      case 8:
        break;
    }
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC gxxL = g11L;
    
    CCTK_REAL_VEC gxyL = g12L;
    
    CCTK_REAL_VEC gxzL = g13L;
    
    CCTK_REAL_VEC gyyL = g22L;
    
    CCTK_REAL_VEC gyzL = g23L;
    
    CCTK_REAL_VEC gzzL = g33L;
    
    CCTK_REAL_VEC kxxL = K11L;
    
    CCTK_REAL_VEC kxyL = K12L;
    
    CCTK_REAL_VEC kxzL = K13L;
    
    CCTK_REAL_VEC kyyL = K22L;
    
    CCTK_REAL_VEC kyzL = K23L;
    
    CCTK_REAL_VEC kzzL = K33L;
    
    CCTK_REAL_VEC alpL = alphaL;
    
    CCTK_REAL_VEC dtalpL = ToReal(0);
    
    CCTK_REAL_VEC betaxL = beta1L;
    
    CCTK_REAL_VEC betayL = beta2L;
    
    CCTK_REAL_VEC betazL = beta3L;
    
    CCTK_REAL_VEC dtbetaxL = ToReal(0);
    
    CCTK_REAL_VEC dtbetayL = ToReal(0);
    
    CCTK_REAL_VEC dtbetazL = ToReal(0);
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(alp[index],alpL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(betax[index],betaxL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(betay[index],betayL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(betaz[index],betazL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(dtalp[index],dtalpL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(dtbetax[index],dtbetaxL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(dtbetay[index],dtbetayL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(dtbetaz[index],dtbetazL,elt_count_lo,elt_count_hi);
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
      vec_store_nta_partial_hi(dtalp[index],dtalpL,elt_count);
      vec_store_nta_partial_hi(dtbetax[index],dtbetaxL,elt_count);
      vec_store_nta_partial_hi(dtbetay[index],dtbetayL,elt_count);
      vec_store_nta_partial_hi(dtbetaz[index],dtbetazL,elt_count);
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
      vec_store_nta_partial_lo(dtalp[index],dtalpL,elt_count);
      vec_store_nta_partial_lo(dtbetax[index],dtbetaxL,elt_count);
      vec_store_nta_partial_lo(dtbetay[index],dtbetayL,elt_count);
      vec_store_nta_partial_lo(dtbetaz[index],dtbetazL,elt_count);
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
    vec_store_nta(dtalp[index],dtalpL);
    vec_store_nta(dtbetax[index],dtbetaxL);
    vec_store_nta(dtbetay[index],dtbetayL);
    vec_store_nta(dtbetaz[index],dtbetazL);
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
  LC_ENDLOOP3VEC (ML_ADM_convertToADMBase);
}

extern "C" void ML_ADM_convertToADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_ADM_convertToADMBase_Body);
}
