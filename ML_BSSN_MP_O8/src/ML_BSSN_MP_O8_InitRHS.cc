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

static void ML_BSSN_MP_O8_InitRHS_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_O8_InitRHS_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_O8_InitRHS_calc_every != ML_BSSN_MP_O8_InitRHS_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN_MP_O8::ML_curvrhs","ML_BSSN_MP_O8::ML_dtlapserhs","ML_BSSN_MP_O8::ML_dtshiftrhs","ML_BSSN_MP_O8::ML_Gammarhs","ML_BSSN_MP_O8::ML_lapserhs","ML_BSSN_MP_O8::ML_log_confacrhs","ML_BSSN_MP_O8::ML_metricrhs","ML_BSSN_MP_O8::ML_shiftrhs","ML_BSSN_MP_O8::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_O8_InitRHS", 9, groups);
  
  
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
  LC_LOOP3VEC (ML_BSSN_MP_O8_InitRHS,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC phirhsL = ToReal(0);
    
    CCTK_REAL_VEC gt11rhsL = ToReal(0);
    
    CCTK_REAL_VEC gt12rhsL = ToReal(0);
    
    CCTK_REAL_VEC gt13rhsL = ToReal(0);
    
    CCTK_REAL_VEC gt22rhsL = ToReal(0);
    
    CCTK_REAL_VEC gt23rhsL = ToReal(0);
    
    CCTK_REAL_VEC gt33rhsL = ToReal(0);
    
    CCTK_REAL_VEC trKrhsL = ToReal(0);
    
    CCTK_REAL_VEC At11rhsL = ToReal(0);
    
    CCTK_REAL_VEC At12rhsL = ToReal(0);
    
    CCTK_REAL_VEC At13rhsL = ToReal(0);
    
    CCTK_REAL_VEC At22rhsL = ToReal(0);
    
    CCTK_REAL_VEC At23rhsL = ToReal(0);
    
    CCTK_REAL_VEC At33rhsL = ToReal(0);
    
    CCTK_REAL_VEC Xt1rhsL = ToReal(0);
    
    CCTK_REAL_VEC Xt2rhsL = ToReal(0);
    
    CCTK_REAL_VEC Xt3rhsL = ToReal(0);
    
    CCTK_REAL_VEC alpharhsL = ToReal(0);
    
    CCTK_REAL_VEC ArhsL = ToReal(0);
    
    CCTK_REAL_VEC beta1rhsL = ToReal(0);
    
    CCTK_REAL_VEC beta2rhsL = ToReal(0);
    
    CCTK_REAL_VEC beta3rhsL = ToReal(0);
    
    CCTK_REAL_VEC B1rhsL = ToReal(0);
    
    CCTK_REAL_VEC B2rhsL = ToReal(0);
    
    CCTK_REAL_VEC B3rhsL = ToReal(0);
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(alpharhs[index],alpharhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Arhs[index],ArhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At11rhs[index],At11rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At12rhs[index],At12rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At13rhs[index],At13rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At22rhs[index],At22rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At23rhs[index],At23rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At33rhs[index],At33rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B1rhs[index],B1rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B2rhs[index],B2rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B3rhs[index],B3rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta1rhs[index],beta1rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta2rhs[index],beta2rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta3rhs[index],beta3rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt11rhs[index],gt11rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt12rhs[index],gt12rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt13rhs[index],gt13rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt22rhs[index],gt22rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt23rhs[index],gt23rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt33rhs[index],gt33rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(phirhs[index],phirhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(trKrhs[index],trKrhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt1rhs[index],Xt1rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt2rhs[index],Xt2rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt3rhs[index],Xt3rhsL,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(alpharhs[index],alpharhsL,elt_count);
      vec_store_nta_partial_hi(Arhs[index],ArhsL,elt_count);
      vec_store_nta_partial_hi(At11rhs[index],At11rhsL,elt_count);
      vec_store_nta_partial_hi(At12rhs[index],At12rhsL,elt_count);
      vec_store_nta_partial_hi(At13rhs[index],At13rhsL,elt_count);
      vec_store_nta_partial_hi(At22rhs[index],At22rhsL,elt_count);
      vec_store_nta_partial_hi(At23rhs[index],At23rhsL,elt_count);
      vec_store_nta_partial_hi(At33rhs[index],At33rhsL,elt_count);
      vec_store_nta_partial_hi(B1rhs[index],B1rhsL,elt_count);
      vec_store_nta_partial_hi(B2rhs[index],B2rhsL,elt_count);
      vec_store_nta_partial_hi(B3rhs[index],B3rhsL,elt_count);
      vec_store_nta_partial_hi(beta1rhs[index],beta1rhsL,elt_count);
      vec_store_nta_partial_hi(beta2rhs[index],beta2rhsL,elt_count);
      vec_store_nta_partial_hi(beta3rhs[index],beta3rhsL,elt_count);
      vec_store_nta_partial_hi(gt11rhs[index],gt11rhsL,elt_count);
      vec_store_nta_partial_hi(gt12rhs[index],gt12rhsL,elt_count);
      vec_store_nta_partial_hi(gt13rhs[index],gt13rhsL,elt_count);
      vec_store_nta_partial_hi(gt22rhs[index],gt22rhsL,elt_count);
      vec_store_nta_partial_hi(gt23rhs[index],gt23rhsL,elt_count);
      vec_store_nta_partial_hi(gt33rhs[index],gt33rhsL,elt_count);
      vec_store_nta_partial_hi(phirhs[index],phirhsL,elt_count);
      vec_store_nta_partial_hi(trKrhs[index],trKrhsL,elt_count);
      vec_store_nta_partial_hi(Xt1rhs[index],Xt1rhsL,elt_count);
      vec_store_nta_partial_hi(Xt2rhs[index],Xt2rhsL,elt_count);
      vec_store_nta_partial_hi(Xt3rhs[index],Xt3rhsL,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(alpharhs[index],alpharhsL,elt_count);
      vec_store_nta_partial_lo(Arhs[index],ArhsL,elt_count);
      vec_store_nta_partial_lo(At11rhs[index],At11rhsL,elt_count);
      vec_store_nta_partial_lo(At12rhs[index],At12rhsL,elt_count);
      vec_store_nta_partial_lo(At13rhs[index],At13rhsL,elt_count);
      vec_store_nta_partial_lo(At22rhs[index],At22rhsL,elt_count);
      vec_store_nta_partial_lo(At23rhs[index],At23rhsL,elt_count);
      vec_store_nta_partial_lo(At33rhs[index],At33rhsL,elt_count);
      vec_store_nta_partial_lo(B1rhs[index],B1rhsL,elt_count);
      vec_store_nta_partial_lo(B2rhs[index],B2rhsL,elt_count);
      vec_store_nta_partial_lo(B3rhs[index],B3rhsL,elt_count);
      vec_store_nta_partial_lo(beta1rhs[index],beta1rhsL,elt_count);
      vec_store_nta_partial_lo(beta2rhs[index],beta2rhsL,elt_count);
      vec_store_nta_partial_lo(beta3rhs[index],beta3rhsL,elt_count);
      vec_store_nta_partial_lo(gt11rhs[index],gt11rhsL,elt_count);
      vec_store_nta_partial_lo(gt12rhs[index],gt12rhsL,elt_count);
      vec_store_nta_partial_lo(gt13rhs[index],gt13rhsL,elt_count);
      vec_store_nta_partial_lo(gt22rhs[index],gt22rhsL,elt_count);
      vec_store_nta_partial_lo(gt23rhs[index],gt23rhsL,elt_count);
      vec_store_nta_partial_lo(gt33rhs[index],gt33rhsL,elt_count);
      vec_store_nta_partial_lo(phirhs[index],phirhsL,elt_count);
      vec_store_nta_partial_lo(trKrhs[index],trKrhsL,elt_count);
      vec_store_nta_partial_lo(Xt1rhs[index],Xt1rhsL,elt_count);
      vec_store_nta_partial_lo(Xt2rhs[index],Xt2rhsL,elt_count);
      vec_store_nta_partial_lo(Xt3rhs[index],Xt3rhsL,elt_count);
      break;
    }
    
    /* Copy local copies back to grid functions */
    vec_store_nta(alpharhs[index],alpharhsL);
    vec_store_nta(Arhs[index],ArhsL);
    vec_store_nta(At11rhs[index],At11rhsL);
    vec_store_nta(At12rhs[index],At12rhsL);
    vec_store_nta(At13rhs[index],At13rhsL);
    vec_store_nta(At22rhs[index],At22rhsL);
    vec_store_nta(At23rhs[index],At23rhsL);
    vec_store_nta(At33rhs[index],At33rhsL);
    vec_store_nta(B1rhs[index],B1rhsL);
    vec_store_nta(B2rhs[index],B2rhsL);
    vec_store_nta(B3rhs[index],B3rhsL);
    vec_store_nta(beta1rhs[index],beta1rhsL);
    vec_store_nta(beta2rhs[index],beta2rhsL);
    vec_store_nta(beta3rhs[index],beta3rhsL);
    vec_store_nta(gt11rhs[index],gt11rhsL);
    vec_store_nta(gt12rhs[index],gt12rhsL);
    vec_store_nta(gt13rhs[index],gt13rhsL);
    vec_store_nta(gt22rhs[index],gt22rhsL);
    vec_store_nta(gt23rhs[index],gt23rhsL);
    vec_store_nta(gt33rhs[index],gt33rhsL);
    vec_store_nta(phirhs[index],phirhsL);
    vec_store_nta(trKrhs[index],trKrhsL);
    vec_store_nta(Xt1rhs[index],Xt1rhsL);
    vec_store_nta(Xt2rhs[index],Xt2rhsL);
    vec_store_nta(Xt3rhs[index],Xt3rhsL);
  }
  LC_ENDLOOP3VEC (ML_BSSN_MP_O8_InitRHS);
}

extern "C" void ML_BSSN_MP_O8_InitRHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_BSSN_MP_O8_InitRHS_Body);
}
