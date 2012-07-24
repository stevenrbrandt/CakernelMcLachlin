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
#include "cctk_Loop.h"
#include "loopcontrol.h"
#include "vectors.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define ScalarINV(x) ((CCTK_REAL)1.0 / (x))
#define ScalarSQR(x) ((x) * (x))
#define ScalarCUB(x) ((x) * ScalarSQR(x))
#define ScalarQAD(x) (ScalarSQR(ScalarSQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))
#define QAD(x) (SQR(SQR(x)))

static void ML_ADM_convertFromADMBase_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
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
  CCTK_REAL_VEC const t = ToReal(cctk_time);
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
  CCTK_REAL_VEC const p1o12dx = kdiv(ToReal(0.0833333333333333333333333333333),dx);
  CCTK_REAL_VEC const p1o12dy = kdiv(ToReal(0.0833333333333333333333333333333),dy);
  CCTK_REAL_VEC const p1o12dz = kdiv(ToReal(0.0833333333333333333333333333333),dz);
  CCTK_REAL_VEC const p1o144dxdy = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dy,dx));
  CCTK_REAL_VEC const p1o144dxdz = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dz,dx));
  CCTK_REAL_VEC const p1o144dydz = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dz,dy));
  CCTK_REAL_VEC const p1o180dx2 = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dx,dx));
  CCTK_REAL_VEC const p1o180dy2 = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dy,dy));
  CCTK_REAL_VEC const p1o180dz2 = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dz,dz));
  CCTK_REAL_VEC const p1o2dx = kdiv(ToReal(0.5),dx);
  CCTK_REAL_VEC const p1o2dy = kdiv(ToReal(0.5),dy);
  CCTK_REAL_VEC const p1o2dz = kdiv(ToReal(0.5),dz);
  CCTK_REAL_VEC const p1o3600dxdy = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dy,dx));
  CCTK_REAL_VEC const p1o3600dxdz = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dz,dx));
  CCTK_REAL_VEC const p1o3600dydz = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dz,dy));
  CCTK_REAL_VEC const p1o4dxdy = kdiv(ToReal(0.25),kmul(dy,dx));
  CCTK_REAL_VEC const p1o4dxdz = kdiv(ToReal(0.25),kmul(dz,dx));
  CCTK_REAL_VEC const p1o4dydz = kdiv(ToReal(0.25),kmul(dz,dy));
  CCTK_REAL_VEC const p1o5040dx2 = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dx,dx));
  CCTK_REAL_VEC const p1o5040dy2 = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dy,dy));
  CCTK_REAL_VEC const p1o5040dz2 = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dz,dz));
  CCTK_REAL_VEC const p1o60dx = kdiv(ToReal(0.0166666666666666666666666666667),dx);
  CCTK_REAL_VEC const p1o60dy = kdiv(ToReal(0.0166666666666666666666666666667),dy);
  CCTK_REAL_VEC const p1o60dz = kdiv(ToReal(0.0166666666666666666666666666667),dz);
  CCTK_REAL_VEC const p1o705600dxdy = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dx));
  CCTK_REAL_VEC const p1o705600dxdz = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dz,dx));
  CCTK_REAL_VEC const p1o705600dydz = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dz,dy));
  CCTK_REAL_VEC const p1o840dx = kdiv(ToReal(0.00119047619047619047619047619048),dx);
  CCTK_REAL_VEC const p1o840dy = kdiv(ToReal(0.00119047619047619047619047619048),dy);
  CCTK_REAL_VEC const p1o840dz = kdiv(ToReal(0.00119047619047619047619047619048),dz);
  CCTK_REAL_VEC const p1odx2 = kdiv(ToReal(1),kmul(dx,dx));
  CCTK_REAL_VEC const p1ody2 = kdiv(ToReal(1),kmul(dy,dy));
  CCTK_REAL_VEC const p1odz2 = kdiv(ToReal(1),kmul(dz,dz));
  CCTK_REAL_VEC const pm1o12dx2 = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dx,dx));
  CCTK_REAL_VEC const pm1o12dy2 = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dy,dy));
  CCTK_REAL_VEC const pm1o12dz2 = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dz,dz));
  
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
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC(ML_ADM_convertFromADMBase,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
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
    CCTK_REAL_VEC g11L = gxxL;
    
    CCTK_REAL_VEC g12L = gxyL;
    
    CCTK_REAL_VEC g13L = gxzL;
    
    CCTK_REAL_VEC g22L = gyyL;
    
    CCTK_REAL_VEC g23L = gyzL;
    
    CCTK_REAL_VEC g33L = gzzL;
    
    CCTK_REAL_VEC K11L = kxxL;
    
    CCTK_REAL_VEC K12L = kxyL;
    
    CCTK_REAL_VEC K13L = kxzL;
    
    CCTK_REAL_VEC K22L = kyyL;
    
    CCTK_REAL_VEC K23L = kyzL;
    
    CCTK_REAL_VEC K33L = kzzL;
    
    CCTK_REAL_VEC alphaL = alpL;
    
    CCTK_REAL_VEC beta1L = betaxL;
    
    CCTK_REAL_VEC beta2L = betayL;
    
    CCTK_REAL_VEC beta3L = betazL;
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,lc_imin,lc_imax);
    vec_store_nta_partial(alpha[index],alphaL);
    vec_store_nta_partial(beta1[index],beta1L);
    vec_store_nta_partial(beta2[index],beta2L);
    vec_store_nta_partial(beta3[index],beta3L);
    vec_store_nta_partial(g11[index],g11L);
    vec_store_nta_partial(g12[index],g12L);
    vec_store_nta_partial(g13[index],g13L);
    vec_store_nta_partial(g22[index],g22L);
    vec_store_nta_partial(g23[index],g23L);
    vec_store_nta_partial(g33[index],g33L);
    vec_store_nta_partial(K11[index],K11L);
    vec_store_nta_partial(K12[index],K12L);
    vec_store_nta_partial(K13[index],K13L);
    vec_store_nta_partial(K22[index],K22L);
    vec_store_nta_partial(K23[index],K23L);
    vec_store_nta_partial(K33[index],K33L);
  }
  LC_ENDLOOP3VEC(ML_ADM_convertFromADMBase);
}

extern "C" void ML_ADM_convertFromADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_convertFromADMBase_Body");
  }
  
  if (cctk_iteration % ML_ADM_convertFromADMBase_calc_every != ML_ADM_convertFromADMBase_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "ADMBase::curv",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "ML_ADM::ML_curv",
    "ML_ADM::ML_lapse",
    "ML_ADM::ML_metric",
    "ML_ADM::ML_shift"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADM_convertFromADMBase", 8, groups);
  
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
  
  GenericFD_LoopOverEverything(cctkGH, ML_ADM_convertFromADMBase_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_ADM_convertFromADMBase_Body");
  }
}
